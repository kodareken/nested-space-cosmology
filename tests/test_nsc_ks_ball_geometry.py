"""Containment controls for KS ball jets; samples are not bound proofs."""
from fractions import Fraction as Q

import mpmath as mp
import numpy as np
import pytest
pytest.importorskip('flint')
from flint import acb, arb, arb_series, ctx

from recursive_horizons.nsc_ks_ball_geometry import (
    MAX_JET_ORDER, SubdivisionNeeded, background_series, inv_a_integrand,
    local_axial_series, plateau_series, restricted_axial_lower_bound,
    shift_series,
)
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction


def _mp_arb(value):
    return arb(str(value))


def _contains_jet(series, values):
    for coeff, value in zip(series.coeffs(), values):
        point = _mp_arb(value)
        assert coeff.contains(point) or coeff.overlaps(point)


def _mp_axial(rho):
    angle = mp.pi / 2 - mp.atan(rho)
    return mp.sqrt(3 * ((1 + rho * rho) * angle - rho) - 1)


def test_restricted_slab_axial_lower_bound_uses_rational_pi_and_atan():
    angle_lo = Q(314, 100) / 4 - Q(1, 65)
    rho = Q(33, 32)
    a2_lo = 3 * ((1 + rho * rho) * angle_lo - rho) - 1
    assert Q(314, 100) < Q(22, 7)
    assert Q(1, 65) > 0
    assert rho * (Q(22, 7) / 4) == Q(363, 448) < 1
    assert a2_lo >= Q(16, 25)
    assert restricted_axial_lower_bound() == Q(4, 5)


def test_background_and_shift_jets_contain_high_precision_point_derivatives():
    order = 6
    samples = (mp.mpf(1), mp.mpf('1.01'), mp.mpf(33) / 32)
    with mp.workdps(60):
        for rho in samples:
            bg = background_series(arb(str(rho)), order)
            shift = shift_series(arb(str(rho)), order)
            r_jet = [mp.diff(lambda t: mp.sqrt(1 + t * t), rho, n) / mp.factorial(n)
                     for n in range(order + 1)]
            a_jet = [mp.diff(_mp_axial, rho, n) / mp.factorial(n)
                     for n in range(order + 1)]
            s_jet = [(-mp.quad(lambda t: 1 / _mp_axial(t), [1, rho])) if n == 0
                     else mp.diff(lambda t: -1 / _mp_axial(t), rho, n - 1) / mp.factorial(n)
                     for n in range(order + 1)]
            _contains_jet(bg.r, r_jet)
            _contains_jet(bg.a, a_jet)
            _contains_jet(bg.inv_a, [1 / v for v in a_jet[:1]] + [
                mp.diff(lambda t: 1 / _mp_axial(t), rho, n) / mp.factorial(n)
                for n in range(1, order + 1)])
            _contains_jet(shift, s_jet)
            assert bg.a[0] >= arb(4) / 5
    ball = arb('1.01 +/- 1e-4')
    bg_i = background_series(ball, 4)
    shift_i = shift_series(ball, 4)
    with mp.workdps(50):
        for rho in (mp.mpf('1.0099'), mp.mpf('1.01'), mp.mpf('1.0101')):
            for n, coeff in enumerate(bg_i.a.coeffs()):
                value = mp.diff(_mp_axial, rho, n) / mp.factorial(n)
                assert coeff.contains(_mp_arb(value)) or coeff.overlaps(_mp_arb(value))
            value = -mp.quad(lambda t: 1 / _mp_axial(t), [1, rho])
            assert shift_i[0].contains(_mp_arb(value)) or shift_i[0].overlaps(_mp_arb(value))


def _mp_eta(distance, inner=.03, outer=.06):
    inner, outer = mp.mpf(inner), mp.mpf(outer)
    if distance <= inner:
        return mp.mpf(1)
    if distance >= outer:
        return mp.mpf(0)
    u = (distance - inner) / (outer - inner)
    return 1 / (1 + mp.exp(1 / (1 - u) - 1 / u))


def _plateau_identity(distance, order):
    return plateau_series(arb_series([distance, arb(1)], prec=order + 1), .03, .06, order)


def test_plateau_interior_and_both_flat_boundary_cells_contain_point_jets():
    order = 6
    interior = arb('0.045')
    left = arb('0.03 +/- 4e-5')
    right = arb('0.06 +/- 4e-5')
    with mp.workdps(50):
        jet = _plateau_identity(interior, order)
        values = [mp.diff(lambda t: _mp_eta(t), mp.mpf('0.045'), n) / mp.factorial(n)
                  for n in range(order + 1)]
        _contains_jet(jet, values)
        for ball, points in (
            (left, (mp.mpf('0.02997'), mp.mpf('0.03'), mp.mpf('0.03003'))),
            (right, (mp.mpf('0.05997'), mp.mpf('0.06'), mp.mpf('0.06003'))),
        ):
            series = _plateau_identity(ball, order)
            for point in points:
                for n, coeff in enumerate(series.coeffs()):
                    value = mp.diff(lambda t, p=point: _mp_eta(t), point, n) / mp.factorial(n)
                    assert coeff.contains(_mp_arb(value)) or coeff.overlaps(_mp_arb(value))
    ones = _plateau_identity(arb('0.02'), 4)
    zeros = _plateau_identity(arb('0.07'), 4)
    assert ones[0] == 1 and all(c.is_zero() for c in ones.coeffs()[1:])
    assert all(c.is_zero() for c in zeros.coeffs())


def test_local_axial_series_matches_binary_mapparms_and_plateau():
    profile = LocalAxialFunction((.001, -.002, .0001), 1.2068780671724793)
    order = 5
    off, scale = map(float, profile._derivatives[0].mapparms())

    def mp_fun(z):
        x = mp.mpf(off) + mp.mpf(scale) * z
        b1, b2 = mp.mpf(0), mp.mpf(0)
        for ak in reversed([mp.mpf(float(c)) for c in profile.coefficients[1:]]):
            b1, b2 = ak + 2 * x * b1 - b2, b1
        poly = mp.mpf(float(profile.coefficients[0])) + x * b1 - b2
        return _mp_eta(abs(z - mp.mpf(float(profile.center))),
                       float(profile.inner), float(profile.outer)) * poly

    with mp.workdps(45):
        for z, ball in (
            (mp.mpf(float(profile.center)), arb(float(profile.center))),
            (mp.mpf(float(profile.center) + .045), arb(float(profile.center) + .045)),
            (mp.mpf(float(profile.center) + .03),
             arb(float(profile.center) + .03) + arb(0, 4e-5)),
            (mp.mpf(float(profile.center) + .06),
             arb(float(profile.center) + .06) + arb(0, 4e-5)),
        ):
            series = local_axial_series(profile, ball, order)
            for n, coeff in enumerate(series.coeffs()):
                value = mp.diff(mp_fun, z, n) / mp.factorial(n)
                assert coeff.contains(_mp_arb(value)) or coeff.overlaps(_mp_arb(value))
    # Binary evaluation roundoff is separate from the analytic jet enclosure.
    sample = profile(float(profile.center + .045), 0)
    analytic = local_axial_series(profile, arb(float(profile.center) + .045), 0)[0]
    assert abs(float(analytic.mid())-sample) < 3e-18
    center_cell = local_axial_series(profile, arb(float(profile.center), .01), 2)
    assert all(c.is_finite() for c in center_cell.coeffs())


def test_context_precision_is_restored_after_success_and_failure():
    ctx.default()
    ctx.dps, ctx.cap = 22, 7
    state = (ctx.prec, ctx.dps, ctx.cap)
    try:
        background_series(1, 4)
        shift_series(arb('1.01'), 3)
        plateau_series(arb_series([arb('0.045'), 1]), .03, .06, 4)
        assert (ctx.prec, ctx.dps, ctx.cap) == state
        for call in (
            lambda: background_series(.5, 2),
            lambda: plateau_series(arb_series([arb('0.045 +/- 0.02'), 1]), .03, .06, 3),
            lambda: local_axial_series(
                LocalAxialFunction((.1,), 0.), arb(0, .02), 2),
        ):
            try:
                call()
            except (ValueError, SubdivisionNeeded):
                pass
            assert (ctx.prec, ctx.dps, ctx.cap) == state
    finally:
        ctx.default()


def test_subdivision_and_analytic_sqrt_flag():
    with np.testing.assert_raises(SubdivisionNeeded):
        plateau_series(arb_series([arb('0.045 +/- 0.02'), 1]), .03, .06, 3)
    # Decimal-radius parsing can enlarge the ball outside the admitted slab.
    with np.testing.assert_raises(ValueError):
        background_series(arb('1.015625 +/- 0.015625'), 2)
    z = acb(-2)
    assert not inv_a_integrand(z, True).is_finite()
    assert inv_a_integrand(acb(1), True).is_finite()
    assert inv_a_integrand(acb(1), False).is_finite()
    assert MAX_JET_ORDER == 16
