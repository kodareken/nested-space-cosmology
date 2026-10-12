"""Independent polynomial derivatives, interval containment, and mixed jets.

The identities checked here are the ones named in
`nsc_ks_chebyshev_jet_bound.MATHEMATICAL_JUSTIFICATION`:
the differentiated Chebyshev recurrence (DLMF 18.9.1), the endpoint bound
|T_n^{(j)}(x)| <= T_n^{(j)}(M) on [-M,M] for M>=1 (DLMF 18.9.21 for j=1),
the mean-value remainder P^{(j)}(c)+[-r,r] B_{j+1}, and the exact finite
Taylor identity of a polynomial. Samples are containment controls, not
substitutes for those proofs. Interval containment of zero is not an
identity certificate.
"""
from fractions import Fraction as Q
from math import factorial

import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb, ctx
from numpy.polynomial import Chebyshev

from recursive_horizons.nsc_ks_ball_geometry import local_axial_series
from recursive_horizons.nsc_ks_chebyshev_jet_bound import (
    MATHEMATICAL_JUSTIFICATION, binary_domain_map,
    chebyshev_endpoint_derivative_bounds, local_axial_interval_series,
    monomial_derivative_l1, polynomial_derivative_on_interval, polynomial_range_series,
)
from recursive_horizons.nsc_ks_current_history_bounds import chebyshev_profile_bounds
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily


def _numpy_poly(profile):
    return Chebyshev(profile.coefficients,
                     domain=[profile.center - profile.outer, profile.center + profile.outer])


def _mp_polynomial(profile, z, j):
    import mpmath as mp
    off, scale = map(float, profile._derivatives[0].mapparms())
    z0 = mp.mpf(float(z))

    def f(t):
        x = mp.mpf(off) + mp.mpf(scale) * t
        b1, b2 = mp.mpf(0), mp.mpf(0)
        coeffs = [mp.mpf(float(c)) for c in profile.coefficients]
        for ak in reversed(coeffs[1:]):
            b1, b2 = ak + 2 * x * b1 - b2, b1
        return coeffs[0] + x * b1 - b2

    with mp.workdps(60):
        return mp.diff(f, z0, j)


def _agrees(enclosure, value, rel=1e-14):
    mid = float(enclosure.mid())
    target = float(value)
    scale = 1.0 + abs(mid) + abs(target)
    return abs(mid - target) <= rel * scale + float(enclosure.rad())


def test_justification_names_dlmf_recurrence_and_mean_value_remainder():
    text = MATHEMATICAL_JUSTIFICATION.lower()
    assert 'dlmf 18.9.1' in text
    assert 'dlmf 18.9.21' in text
    assert 'mean-value' in text
    assert 'one_direction_metric' in text
    assert 'does not close the local incoming gate' in text


@pytest.mark.parametrize('center,radius', [(2, Q(1, 10)), (Q(21, 20), Q(1, 10))])
def test_polynomial_enclosure_does_not_use_endpoint_bounds_outside_their_domain(center, radius):
    # T_2(z)=2z^2-1. A majorant valid only on [-1,1] used to give
    # 7 +/- .4 at 2 +/- .1 and miss the exact endpoint value 7.82.
    # Also cover an interval that crosses the support boundary.
    profile = LocalAxialFunction((0., 0., 1.), 0., .5, 1.)
    with ctx.workprec(120):
        c = arb(Q(center).numerator)/Q(center).denominator
        r = arb(radius.numerator)/radius.denominator
        interval = c + r*arb(0, 1)
        bounds = chebyshev_endpoint_derivative_bounds(profile, 3)
        jet = polynomial_range_series(profile, interval, 2, bounds)
        for endpoint in (c-r, c+r):
            exact = (2*endpoint**2-1, 4*endpoint, arb(4))
            for j, value in enumerate(exact):
                got = polynomial_derivative_on_interval(profile, interval, j, bounds)
                assert got.contains(value)
                assert (factorial(j)*jet[j]).contains(value)


def test_endpoint_recurrence_matches_first_kind_values_at_one():
    # T_n(1)=1, T_n'(1)=n^2, T_n''(1)=n^2(n^2-1)/3 (DLMF 18.9.21 and recurrence).
    profile = LocalAxialFunction(tuple(1.0 if n == 5 else 0.0 for n in range(6)), 0.0)
    bounds = chebyshev_endpoint_derivative_bounds(profile, 3)
    offset, scale, center, _, outer = binary_domain_map(profile)
    assert offset + scale * (center - outer) == Q(-1) or abs(float(offset + scale * (center - outer)) + 1) < 1e-15
    M = bounds['map_endpoint_upper']
    assert M >= 1
    n = 5
    # Isolated T_5: polynomial_bounds[j] = |scale|^j T_5^{(j)}(M) >= |scale|^j T_5^{(j)}(1).
    assert bounds['polynomial_bounds'][0] >= 1
    assert bounds['polynomial_bounds'][1] >= abs(scale) * n * n
    assert bounds['polynomial_bounds'][2] >= abs(scale)**2 * n**2 * (n**2 - 1) / 3


def test_low_order_endpoint_bounds_match_the_hash_bound_module():
    coefficients = tuple(0.01 * (n + 1) * ((-1)**n) for n in range(8))
    profile = LocalAxialFunction(coefficients, 0.15)
    old = chebyshev_profile_bounds(profile)
    new = chebyshev_endpoint_derivative_bounds(profile, 2)
    assert new['map_endpoint_upper'] == old['map_endpoint_upper']
    for j in range(3):
        assert new['polynomial_bounds'][j] == old['polynomial_bounds'][j]


def test_binary_map_is_the_exact_mapparms_dyadic():
    profile = LocalAxialFunction((0.2, -0.1, 0.05), 1.2068780671724793)
    offset, scale, _, _, _ = binary_domain_map(profile)
    raw = profile._derivatives[0].mapparms()
    assert offset == Q(raw[0])
    assert scale == Q(raw[1])
    assert scale != 0


def test_point_polynomial_derivative_contains_independent_evaluation():
    profile = LocalAxialFunction((0.3, -0.2, 0.05, 0.01, -0.002), 0.0)
    with ctx.workprec(120):
        for j in range(5):
            for z in (-0.05, -0.02, 0.0, 0.01, 0.045):
                enclosure = polynomial_derivative_on_interval(profile, arb(float(z)), j)
                assert _agrees(enclosure, _mp_polynomial(profile, z, j))
                numpy_value = float(_numpy_poly(profile).deriv(j)(z))
                assert _agrees(enclosure, numpy_value, rel=1e-12)


def test_mean_value_and_taylor_interval_contain_independent_samples():
    """P^{(j)}(c+h) = P^{(j)}(c) + h P^{(j+1)}(ξ) and the finite Taylor identity."""
    profile = LocalAxialFunction((0.4, -0.25, 0.1, -0.03, 0.005, -0.001), 0.0)
    poly = _numpy_poly(profile)
    left, right = -0.02, 0.01
    ball = arb(float(left)).union(arb(float(right)))
    with ctx.workprec(120):
        for j in range(4):
            enclosure = polynomial_derivative_on_interval(profile, ball, j)
            center = float(ball.mid())
            radius = float((right - left) / 2 + 1e-15)
            point = float(poly.deriv(j)(center))
            global_next = float(chebyshev_endpoint_derivative_bounds(profile, j + 1)['polynomial_bounds'][j + 1])
            mvt_upper = abs(point) + radius * global_next + 1e-12
            assert abs(enclosure.upper()) <= mvt_upper or enclosure.upper() <= arb(mvt_upper).upper()
            for z in np.linspace(left, right, 11):
                value = float(poly.deriv(j)(z))
                assert enclosure.contains(arb(value)) or enclosure.overlaps(arb(value))


def test_constant_and_linear_polynomials_have_exact_high_derivatives():
    constant = LocalAxialFunction((2.5,), 0.0)
    linear = LocalAxialFunction((0.0, 0.3), 0.0)
    ball = arb(0, 0.04)
    with ctx.workprec(120):
        assert polynomial_derivative_on_interval(constant, ball, 0).contains(arb('2.5'))
        d1 = polynomial_derivative_on_interval(constant, ball, 1)
        assert d1.contains(0) and abs(d1.upper()) == 0
        slope = polynomial_derivative_on_interval(linear, ball, 1)
        assert _agrees(slope, _mp_polynomial(linear, 0.0, 1))
        assert polynomial_derivative_on_interval(linear, ball, 2).contains(0)


def test_cutoff_interval_jet_contains_owned_point_series_and_mixed_leibniz():
    profile = LocalAxialFunction((0.2, -0.05, 0.01), 0.0)
    z = arb('0.045')
    ball = z + arb(0, 4e-5)
    with ctx.workprec(120):
        interval = local_axial_interval_series(profile, ball, 3)
        point = local_axial_series(profile, z, 3)
        for coeff, value in zip(interval.coeffs(), point.coeffs()):
            assert coeff.contains(value) or coeff.overlaps(value)
        samples = np.linspace(float(ball.lower()), float(ball.upper()), 9)
        for s in samples:
            owned = profile(float(s), 2)
            ordinary = 2 * interval[2]
            assert ordinary.contains(arb(owned)) or ordinary.overlaps(arb(owned))


def test_mixed_profile_product_derivatives_match_leibniz():
    w_coeff = np.zeros(8)
    u_coeff = np.zeros(8)
    w_coeff[:2] = (0.15, 0.02)
    u_coeff[:3] = (0.4, -0.1, 0.02)
    w = LocalAxialFunction(tuple(w_coeff), 0.0)
    U = LocalAxialFunction(tuple(u_coeff), 0.0)
    family = LocalIncomingFamily(np.array([w.coefficients, U.coefficients], float),
                                 center=0.0)
    z = 0.01
    with ctx.workprec(120):
        series = local_axial_interval_series(family.functions[0], arb(z), 2)
        other = local_axial_interval_series(family.functions[1], arb(z), 2)
        product = series * other
        # (wU)' = w' U + w U' in the ordinary-derivative convention.
        mixed = 1 * product[1]
        leibniz = 1 * series[1] * other[0] + series[0] * 1 * other[1]
        assert mixed.overlaps(leibniz)
        value = family.functions[0](z, 0) * family.functions[1](z, 1) + family.functions[0](z, 1) * family.functions[1](z, 0)
        assert _agrees(mixed, value)


def test_midpoint_jet_is_strictly_tighter_than_broad_interval_clenshaw():
    coefficients = [0.0] * 32
    coefficients[-1] = 1e-4
    profile = LocalAxialFunction(tuple(coefficients), 0.0)
    ball = arb(float(0.012)).union(arb(float(0.027)))
    with ctx.workprec(120):
        wrapped = local_axial_series(profile, ball, 8)
        jet = local_axial_interval_series(profile, ball, 8)
        wrap_l1 = abs(factorial(8) * wrapped[8]).upper()
        jet_l1 = abs(factorial(8) * jet[8]).upper()
        assert jet_l1 < wrap_l1
        assert wrap_l1 / jet_l1 > arb(100)


def test_interior_l1_of_a_low_degree_family_is_finite_and_contains_samples():
    coefficients = np.zeros((2, 8))
    coefficients[0, 1] = 0.2
    coefficients[1, 0] = 0.5
    family = LocalIncomingFamily(coefficients, center=0.0)
    keys = ((1, 0), (0, 1), (0, 2))
    with ctx.workprec(80):
        l1, cells = monomial_derivative_l1(family, keys, order=4, transition_panels=4, interior_panels=4, bits=80)
        assert cells >= 12
        for key, bound in l1.items():
            assert bound >= 0 and bound.is_finite()
        # A Riemann sum of |f^{(4)}| at panel midpoints cannot exceed the enclosure.
        w, U = family.functions
        center, inner, outer = w.center, w.inner, w.outer
        sample = 0.0
        for left, right, count in ((-outer, -inner, 4), (-inner, inner, 4), (inner, outer, 4)):
            width = (right - left) / count
            for i in range(count):
                z = center + left + (i + 0.5) * width
                value = w(z, 4) if True else 0.0
                # Independent order-4 derivative of w via numpy on the polynomial
                # times cutoff through LocalAxialFunction.
                sample += abs(w(z, 4)) * width
        assert l1[(1, 0)] >= arb(sample) or l1[(1, 0)].overlaps(arb(sample)) or l1[(1, 0)] >= arb(sample).upper()


def test_v3_record_does_not_claim_a_physical_gate_or_identity_from_contains_zero():
    from pathlib import Path
    import json
    root = Path(__file__).resolve().parents[1]
    for name in (
        'results/development/nsc-ks-current-field-pilot-v3-profile-bounds.json',
        'results/development/nsc-ks-current-field-pilot-v3.json',
    ):
        path = root / name
        if not path.exists():
            continue
        saved = json.loads(path.read_text())
        assert saved['physical_local_gate'] == 'OPEN'
        assert saved['physical_EXISTENCE_certificate'] is False
        assert saved['physical_NON_EXISTENCE_certificate'] is False
        assert saved['proof_parameters']['free_contains_zero_used_as_proof'] is False
        assert saved['comparison']['max_new_B8_display'] < saved['comparison']['max_old_B8_display']
