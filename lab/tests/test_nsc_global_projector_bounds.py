"""Uniform-in-k operator-norm majorants of auxiliary P; not a source or state run."""
from math import factorial
import resource

import numpy as np
import pytest

pytest.importorskip('flint')
from flint import arb, ctx

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_global_projector_bounds import (
    ANALYTIC_FACTORS, FROBENIUS_OVER_OPERATOR, NUCLEAR_OVER_OPERATOR,
    global_projector_bounds, rest_gap_lower, star_norm,
)
from recursive_horizons.nsc_ks_spacetime_geometry_jets import spacetime_geometry_jets
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
from recursive_horizons.nsc_scaled_reference_projector import (
    ScalarJet, positive_inv, projector_derivative, scaled_reference_projector,
    symbol_metric_from_geometry,
)
from recursive_horizons.nsc_spatial_reference_symbol import reference_projector


def _cpu():
    return resource.getrusage(resource.RUSAGE_SELF).ru_utime


def _as2(vec):
    array = np.asarray(vec, complex)
    if array.shape == (2, 2):
        return array
    return np.array([[complex(vec[0]), complex(vec[1])],
                     [complex(vec[2]), complex(vec[3])]], complex)


def _opnorm(vec):
    return float(np.linalg.norm(_as2(vec), 2))


def _le(actual, bound, name=''):
    upper = bound.abs_upper() if hasattr(bound, 'abs_upper') else arb(bound)
    assert float(actual) <= float(upper) + 1e-10 * max(1.0, float(upper)), (
        name, actual, float(upper))


def _family(alpha=.001):
    w = LocalAxialFunction((.2, .3, -.1, .07), 1.2)
    U = LocalAxialFunction((.1, -.03, .05), 1.2)
    return CompatibleIncomingMetric((alpha,), (CompatibleRadiusDirection(w, U, .007, .03),))


def _point_geometry(order=4, bits=120):
    return spacetime_geometry_jets(_family(), 1.012, 1.245, order=order, bits=bits)


def test_zero_gap_rejected_without_a_floor():
    with pytest.raises(ValueError, match='rest gap|no denominator floor'):
        global_projector_bounds(mass=0, angular=0, axial=arb('11/10'), radius=arb('7/5'),
                                physical=2, momentum=2, formal_order=2, bits=80)
    with pytest.raises(ValueError, match='rest gap|no denominator floor'):
        global_projector_bounds(mass=0, angular=arb(0, .1), axial=1, radius=2,
                                physical=2, momentum=2, formal_order=2, bits=80)
    with pytest.raises(ValueError, match='rest gap'):
        rest_gap_lower(0, 0, arb('7/5'))


def test_massless_nonzero_ell_is_accepted():
    bounds = global_projector_bounds(
        mass=0, angular=2, axial=arb('11/10'), radius=arb('7/5'),
        physical=3, momentum=3, formal_order=2, bits=80)
    expected = (arb(2) / arb('7/5')).lower()
    assert bounds.g_min > 0
    assert abs(float(bounds.g_min) - float(expected)) < 1e-12
    assert float(bounds.P[0].value) >= 1 - 1e-12
    assert bounds.analytic_factors is ANALYTIC_FACTORS
    assert 'c_min g_min' in ANALYTIC_FACTORS['X_korder1']
    assert '2/' not in ANALYTIC_FACTORS['X_korder1']


def test_directed_denominators_use_lowers_never_uppers_before_inverse():
    a = arb('11/10') + arb(0, '1/100')
    r = arb('7/5') + arb(0, '1/50')
    bounds = global_projector_bounds(
        mass=1, angular=1, axial=a, radius=r, physical=2, momentum=2,
        formal_order=2, bits=120)
    c_box = positive_inv(a)
    assert bounds.c_min <= c_box.lower() or bounds.c_min.overlaps(c_box.lower())
    # 1/a.lower() is c_max, not a conservative c_min.
    c_max = float(positive_inv(a.lower()))
    assert float(bounds.c_min) < c_max - 1e-6
    # g_min uses r_max in the denominator, not r_min.
    r_max = float(r.upper())
    r_min = float(r.lower())
    g_from_rmax = np.sqrt(1.0 + (1.0 / r_max) ** 2)
    g_from_rmin = np.sqrt(1.0 + (1.0 / r_min) ** 2)
    assert float(bounds.g_min) <= g_from_rmax + 1e-10
    assert float(bounds.g_min) < g_from_rmin - 1e-6
    assert FROBENIUS_OVER_OPERATOR > 1
    assert NUCLEAR_OVER_OPERATOR == 2


def test_constant_geometry_higher_P_vanish_and_P0_k_jet_is_nontrivial():
    kwargs = dict(mass=arb(1), angular=arb(2), axial=arb('11/10'), radius=arb('7/5'),
                  physical=4, momentum=4, bits=80)
    bounds = global_projector_bounds(**kwargs)
    assert bounds.P[0].value == 1 or abs(float(bounds.P[0].value) - 1) < 1e-12
    assert bounds.A.value == 1
    assert bounds.X.value == 0
    for j in range(1, 5):
        for t, z, k in bounds.P[j].indices:
            coeff = bounds.P[j].coefficient(t, z, k)
            assert coeff.contains(0) or float(coeff) <= 1e-18
    # P0 still depends on k uniformly.
    assert float(bounds.P[0].coefficient(0, 0, 1)) > 0
    assert float(bounds.P[0].coefficient(0, 0, 2)) > 0
    a, r, m, ell = arb('11/10'), arb('7/5'), arb(1), arb(2)
    c = 1 / a
    v = (m * m + (ell / r) * (ell / r)).sqrt()
    x001 = bounds.X.coefficient(0, 0, 1)
    expected_x = (c * c) / (c * v)
    assert abs(float(x001) - float(expected_x)) < 1e-10
    # AM-GM: ||d_k P0|| bound is 0.75 c/v, not the looser c/v from an extra 2.
    dk = bounds.derivative_bound(0, k=1)
    cv = float(c / v)
    assert 0.5 * cv - 1e-8 <= float(dk) <= 0.75 * cv + 1e-8
    assert float(dk) < cv - 1e-6
    for t, z, k in bounds.X.indices:
        assert bounds.X.coefficient(t, z, k) >= 0


def test_taylor_factorials_match_derivative_bounds():
    bounds = global_projector_bounds(
        mass=1, angular=1, axial=arb('11/10'), radius=arb('6/5'),
        physical=4, momentum=4, bits=80)
    for j, t, z, k in ((0, 0, 0, 2), (0, 1, 0, 1), (1, 0, 1, 0), (2, 0, 0, 0)):
        taylor = bounds.taylor_bound(j, t, z, k)
        deriv = bounds.derivative_bound(j, t, z, k)
        factor = factorial(t) * factorial(z) * factorial(k)
        assert deriv.overlaps(taylor * factor) or abs(float(deriv) - float(taylor) * factor) < 1e-12
    with pytest.raises(ValueError, match='exceeds|invent|owned'):
        bounds.P[4].coefficient(1, 0, 0)
    with pytest.raises(ValueError, match='exceeds|owned|must lie'):
        bounds.derivative_bound(3, z=2)


def test_rejects_callable_metric_and_non_ks_geometry():
    with pytest.raises(TypeError, match='callables|ball'):
        global_projector_bounds(mass=1, angular=1, axial=lambda z: 1, radius=2)
    with pytest.raises(TypeError, match='KSGeometryJets'):
        global_projector_bounds(mass=1, angular=1, geometry={'a': 1})
    with pytest.raises(ValueError, match='cover the requested'):
        global_projector_bounds(mass=1, angular=1, axial=1, radius=2,
                                physical=3, momentum=5, formal_order=4, bits=80)


def _compare_scaled(bounds, mass, angular, geometry, momenta, physical=4, momentum=4, bits=120):
    with ctx.workprec(bits):
        for k0 in momenta:
            sign = 1 if k0 > 0 else -1
            mu = arb(1 / abs(k0))
            jets = scaled_reference_projector(
                mu=mu, sign=sign, mass=mass, angular=angular, geometry=geometry,
                physical=physical, momentum=momentum, bits=bits)
            checks = ((0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 2), (0, 0, 1, 0),
                      (0, 1, 0, 0), (1, 0, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1),
                      (2, 0, 0, 0), (3, 0, 0, 0), (4, 0, 0, 0))
            for j, t, z, k in checks:
                if t + z > physical - j or k > momentum - j:
                    continue
                actual = projector_derivative(jets, j, t=t, z=z, k=k)
                _le(_opnorm(actual), bounds.derivative_bound(j, t, z, k),
                    f'k={k0} sign={sign} P{j} t{t}z{z}k{k}')


def test_all_k_bounds_versus_scaled_projector_both_signs_and_near_zero():
    geometry = _point_geometry(4)
    mass, angular = arb.pi() / 2, arb(5).sqrt()
    bounds = global_projector_bounds(
        mass=mass, angular=angular, geometry=geometry,
        physical=4, momentum=4, bits=120)
    momenta = (-100.0, -1.0, -0.05, 0.05, 1.0, 100.0)
    _compare_scaled(bounds, mass, angular, geometry, momenta)


def test_unscaled_pst_at_finite_k_and_taylor_factorials():
    geometry = _point_geometry(4)
    mass, angular = float(np.pi / 2), float(np.sqrt(5))
    momenta = [-2.0, -0.5, 0.4, 1.7]
    fields = symbol_metric_from_geometry(geometry)
    old = reference_projector(*fields, momenta, [mass] * 4, [angular] * 4)
    bounds = global_projector_bounds(
        mass=mass, angular=angular, geometry=geometry,
        physical=4, momentum=4, bits=120)
    for index, k0 in enumerate(momenta):
        for j in range(5):
            _le(_opnorm(old['orders'][j, index]), bounds.P[j].value, f'PST P{j} at k={k0}')
        for h in (1, 2):
            actual = old['jets'][0].derivative(k=h).value[index]
            _le(_opnorm(actual), bounds.derivative_bound(0, k=h), f'PST d_k^{h} P0 at k={k0}')
            taylor = _opnorm(actual) / factorial(h)
            _le(taylor, bounds.taylor_bound(0, k=h), f'PST Taylor k^{h} P0 at k={k0}')
        actual_z = old['jets'][1].derivative(z=1).value[index]
        _le(_opnorm(actual_z), bounds.derivative_bound(1, z=1), f'PST d_z P1 at k={k0}')
        actual_t = old['jets'][1].derivative(t=1).value[index]
        _le(_opnorm(actual_t), bounds.derivative_bound(1, t=1), f'PST d_T P1 at k={k0}')


def test_geometry_box_contains_point_samples_as_sanity_not_uniformity_proof():
    family = _family()
    with ctx.workprec(120):
        rho = arb(1.012, 1e-6)
        z = arb(1.245, 1e-6)
        box = spacetime_geometry_jets(family, rho, z, order=4, bits=120)
        assert box.coefficients['a'][(0, 0)] > 0
        assert box.coefficients['r'][(0, 0)] > 0
        mass, angular = arb.pi() / 2, arb(5).sqrt()
        bounds = global_projector_bounds(
            mass=mass, angular=angular, geometry=box, physical=4, momentum=4, bits=120)
        for rho_pt, z_pt in ((1.0119995, 1.2449995), (1.012, 1.245), (1.0120005, 1.2450005)):
            point = spacetime_geometry_jets(family, rho_pt, z_pt, order=4, bits=120)
            _compare_scaled(bounds, mass, angular, point, (-1.0, 0.05, 1.0), bits=120)


def test_star_norm_matches_absolute_weyl_coefficient():
    k = ScalarJet.variable(arb(0), 2, 3, 3)
    z = ScalarJet.variable(arb(0), 1, 3, 3)
    product = star_norm(z, k, 1)
    # |i/2| * (z_z k_k) with C(1,0)=1, both first derivatives 1.
    assert abs(float(product.value) - 0.5) < 1e-12
    assert float(star_norm(k, k, 1).value) <= 1e-18


def test_actual_history_pilot_retains_nine_by_five(record_property):
    geometry = spacetime_geometry_jets(_family(), 1.012, 1.245, order=9, bits=160)
    mass, angular = arb.pi() / 2, arb(5).sqrt()
    start = _cpu()
    bounds = global_projector_bounds(
        mass=mass, angular=angular, geometry=geometry,
        physical=9, momentum=5, bits=160)
    cpu = _cpu() - start
    record_property('pilot_cpu_seconds', cpu)
    assert cpu < 90
    assert bounds.retained == ((9, 5), (8, 4), (7, 3), (6, 2), (5, 1))
    for j, jet in enumerate(bounds.P):
        assert (jet.physical, jet.momentum) == bounds.retained[j]
        for t, z, k in jet.indices:
            coeff = jet.coefficient(t, z, k)
            assert coeff.is_finite() and coeff >= 0
    # P4 z^4 is retained; z^5 is not silently dropped to zero, it is refused.
    assert bounds.P[4].coefficient(0, 4, 0).is_finite()
    assert bounds.derivative_bound(4, z=4).is_finite()
    with pytest.raises(ValueError, match='exceeds|owned|must lie'):
        bounds.P[4].coefficient(0, 6, 0)
    with pytest.raises(ValueError, match='exceeds|owned|must lie'):
        bounds.derivative_bound(4, k=2)
    assert bounds.auxiliary and 'P_g + D_g' in bounds.state_law
    assert bounds.history_identity == geometry.history_identity
    assert bounds.k_axis.startswith('ordinary k')
    assert float(bounds.P[0].value) >= 1 - 1e-12
    assert float(bounds.P[1].value) > 0
    jets = scaled_reference_projector(
        mu=arb('1/4'), sign=1, mass=mass, angular=angular,
        geometry=geometry, physical=4, momentum=4, bits=120)
    _le(_opnorm(projector_derivative(jets, 4)), bounds.P[4].value, 'pilot P4 value')
    _le(_opnorm(projector_derivative(jets, 0, z=4)), bounds.derivative_bound(0, z=4),
        'pilot d_z^4 P0')


def test_saved_control_box_massless_and_massive_sanity():
    from recursive_horizons.nsc_ks_source_envelope import axial_center
    center = axial_center()
    w = LocalAxialFunction((0., 1.), center)
    zero = LocalAxialFunction((0.,), center)
    family = CompatibleIncomingMetric((.001,), (CompatibleRadiusDirection(w, zero, .007, .03),))
    with ctx.workprec(120):
        geometry = spacetime_geometry_jets(
            family, arb(1.012, 1e-7), arb(center + .045, 1e-7), order=4, bits=120)
        massive = global_projector_bounds(
            mass=float(np.pi / 2), angular=float(np.sqrt(5)), geometry=geometry,
            physical=4, momentum=4, bits=120)
        massless = global_projector_bounds(
            mass=0, angular=float(np.sqrt(5)), geometry=geometry,
            physical=4, momentum=4, bits=120)
        assert massive.g_min > 0 and massless.g_min > 0
        point = spacetime_geometry_jets(family, 1.012, center + .045, order=4, bits=120)
        _compare_scaled(massive, arb.pi() / 2, arb(5).sqrt(), point, (-0.05, 1.0))
        _compare_scaled(massless, arb(0), arb(5).sqrt(), point, (0.05, -1.0))


def test_intermediate_star_norm_never_pads_uncomputed_derivatives_with_zero():
    from recursive_horizons.nsc_global_projector_bounds import star_norm
    from recursive_horizons.nsc_scaled_reference_projector import ScalarJet
    x=ScalarJet.constant(arb(1),4,4)
    result=star_norm(x,x,3)
    assert (result.physical,result.momentum)==(1,1)
    with pytest.raises(ValueError,match='exceeds'):
        result.derivative(t=2)


def test_fractional_derivative_orders_cannot_be_silently_truncated():
    bounds=global_projector_bounds(mass=1,angular=2,axial=1,radius=2,physical=4,momentum=4)
    with pytest.raises(ValueError,match='integer'):
        bounds.derivative_bound(0,t=1.5)
