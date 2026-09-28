"""Taylor-remainder controls for f_s,z; samples are not bound proofs."""
import json
from pathlib import Path

import mpmath as mp
import pytest
pytest.importorskip('flint')
from flint import arb, arb_series, ctx

from recursive_horizons import nsc_ks_ball_geometry as G
from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_dirac_source_phase_bound import (
    ANGULAR_SQUARE_WEIGHT, CPU_CAP, MAX_CELLS, MAX_DEPTH, SOURCE_SIGNS,
    aggregate_phase_error_uppers, directed_source_phase_z_enclosure,
    distance_series, dyadic_rho_cover, integrand_rho_series,
    integrate_dyadic_cell, pack_real_ball, unpack_real_ball,
    unit_angular_target,
)
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates


CONTROL = Path('results/development/nsc-ks-cutoff-bridge-control.json')


def _family(amplitude=.001, w=(0., 1.), u=(0.,)):
    center = chart_coordinates(1.)[1] + .15
    return CompatibleIncomingMetric(
        (float(amplitude),),
        (CompatibleRadiusDirection(
            LocalAxialFunction(w, center),
            LocalAxialFunction(u, center), .007, .03),),
    )


def _control():
    record = json.loads(CONTROL.read_text())
    return record, _family(), record['z'][0], record['rho_up']


def test_zero_and_constant_profiles_enclose_vanishing_f_z():
    with ctx.workprec(120):
        zero = directed_source_phase_z_enclosure(
            _family(0.0), chart_coordinates(1.)[1] + .15,
            rho_up=1.015625, cpu_limit=5., max_cells=64)
        assert zero.obstruction is None
        for ball in zero.balls:
            assert ball.contains(0)
            assert ball.rad() <= unit_angular_target()
        constant = directed_source_phase_z_enclosure(
            _family(w=(1.,), u=(0.,)), chart_coordinates(1.)[1] + .15,
            rho_up=1.015625, cpu_limit=5., max_cells=64)
        assert constant.obstruction is None
        for ball in constant.balls:
            assert ball.contains(0)
        assert zero.source_tail_error_bound is None
        assert zero.physical_error_bound is None
        assert zero.remainder_error_bound is None
        assert zero.assembly_correction is False


def test_distance_series_contains_independent_integral_and_restores_cap():
    ctx.default()
    ctx.dps, ctx.cap = 22, 7
    state = (ctx.prec, ctx.dps, ctx.cap)
    try:
        with mp.workdps(50):
            def axial(rho):
                angle = mp.pi / 2 - mp.atan(rho)
                return mp.sqrt(3 * ((1 + rho * rho) * angle - rho) - 1)

            rho = 1 + mp.mpf(1) / 128
            value = mp.quad(lambda t: 1 / axial(t)**2, [1, rho])
        point = arb(1) + arb(1) / 128
        with ctx.workprec(120):
            series = distance_series(point, 4)
            sample = arb(str(value)) + arb(0, arb(1) / arb(10)**16)
            assert series[0].overlaps(sample)
        assert (ctx.prec, ctx.dps, ctx.cap) == state
        try:
            G.plateau_series(arb_series([arb('0.045 +/- 0.02'), 1]), .03, .06, 3)
        except G.SubdivisionNeeded:
            pass
        assert (ctx.prec, ctx.dps, ctx.cap) == state
    finally:
        ctx.default()


def test_dyadic_cover_and_rho_one_balls_stay_in_the_slab():
    record = json.loads(CONTROL.read_text())
    scale, cells = dyadic_rho_cover(1.0, record['rho_up'])
    assert scale >= 50
    assert cells[0][0] * arb(2)**(-scale) == 1
    last = arb((cells[-1][1], -scale))
    assert last.contains(arb(record['rho_up'])) or last == arb(record['rho_up'])
    with ctx.workprec(120):
        start, stop = cells[0]
        ball = arb((start + stop, -(scale + 1)), (stop - start, -(scale + 1)))
        assert ball >= 1
        assert ball <= arb(33) / 32
        assert ball.contains(1)


def _mp_eta(distance, inner, outer):
    if distance <= inner:
        return mp.mpf(1)
    if distance >= outer:
        return mp.mpf(0)
    u = (distance - inner) / (outer - inner)
    return 1 / (1 + mp.exp(1 / (1 - u) - 1 / u))


def test_taylor_cell_contains_independent_high_precision_quadrature():
    family = _family()
    model = AnalyticRadiusFamily(family)
    z = float(family.directions[0].w.center)
    scale, _ = dyadic_rho_cover(1.0, 1.015625)
    # [1+1/512, 1+1/256] sits inside the inner normal window.
    start = (1 << scale) + (1 << (scale - 9))
    stop = (1 << scale) + (1 << (scale - 8))
    assert scale >= 9
    with ctx.workprec(120):
        enclosed, extra = integrate_dyadic_cell(model, start, stop, scale, z, 1)
        series = integrand_rho_series(
            model, arb((start + stop, -(scale + 1))), z, 1, 0)
    inner = float(family.directions[0].inner_radius)
    outer = float(family.directions[0].outer_radius)
    center = float(family.directions[0].w.center)
    off, map_scale = map(float, family.directions[0].w._derivatives[0].mapparms())
    amplitude = float(family.amplitudes[0])

    def mp_axial(rho):
        angle = mp.pi / 2 - mp.atan(rho)
        return mp.sqrt(3 * ((1 + rho * rho) * angle - rho) - 1)

    def mp_integrand(rho, sign=1):
        rho = mp.mpf(rho)
        s = -mp.quad(lambda t: 1 / mp_axial(t), [1, rho])
        distance = mp.quad(lambda t: 1 / mp_axial(t)**2, [1, rho])
        chi = _mp_eta(abs(s), inner, outer)
        z_char = z - sign * distance
        x = mp.mpf(off) + mp.mpf(map_scale) * z_char
        window = _mp_eta(abs(z_char - center), .03, .06)
        w = amplitude * window * x
        w_z = amplitude * window * mp.mpf(map_scale)
        radius = mp.sqrt(1 + rho * rho) + chi * (s * w)
        radius_z = chi * s * w_z
        return -radius_z / radius**3

    with mp.workdps(45):
        left = mp.mpf(start) / mp.mpf(2**scale)
        right = mp.mpf(stop) / mp.mpf(2**scale)
        sample = mp.quad(lambda t: mp_integrand(t, 1), [left, right])
        mid = (left + right) / 2
        point = mp_integrand(mid, 1)
    control = arb(0, arb(1) / arb(10)**18)
    assert enclosed.overlaps(arb(str(sample)) + control)
    assert extra >= 0
    assert series[0].overlaps(arb(str(point)) + control)


def test_control_point_enclosure_has_both_signs_and_named_open_limit():
    record, family, z, rho_up = _control()
    assert profile_identity(family) == record['history_identity']
    assert record['signed_angular_square_weights'] == [
        ANGULAR_SQUARE_WEIGHT, ANGULAR_SQUARE_WEIGHT]
    opened = directed_source_phase_z_enclosure(
        family, z, rho_up=rho_up, max_depth=0, max_cells=1, cpu_limit=5.)
    assert opened.status.startswith('OPEN:')
    assert opened.obstruction in {'max_depth', 'max_cells', 'terminal_dyadic_cell'}
    assert opened.coefficient_error_bound is None
    assert opened.diagnostics['gauss_nodes_are_bounds'] is False
    with ctx.workprec(120):
        result = directed_source_phase_z_enclosure(
            family, z, rho_up=rho_up, cpu_limit=CPU_CAP,
            max_depth=MAX_DEPTH, max_cells=MAX_CELLS)
        assert result.obstruction is None
        assert result.signs == SOURCE_SIGNS == (1, -1)
        assert result.profile_identity == record['history_identity']
        assert result.z == z
        assert result.rho_up == rho_up
        target = unit_angular_target()
        for ball in result.balls:
            assert ball.rad() <= target
        newton, beta = aggregate_phase_error_uppers(*result.balls)
        assert newton < arb(1) / arb(10)**11
        assert beta < arb(1) / arb(10)**11
        packed = pack_real_ball(result.balls[0])
        assert unpack_real_ball(packed).overlaps(result.balls[0])
        assert result.cpu_seconds < CPU_CAP
        assert result.source_tail_error_bound is result.physical_error_bound is None
        assert result.physical_local_gate.startswith('OPEN')
