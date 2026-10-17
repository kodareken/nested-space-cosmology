"""A one-stage radial diagnosis authenticates its exact local state."""
import hashlib
import math
from fractions import Fraction as Q
from types import SimpleNamespace
import numpy as np
import pytest

from diagnose_nsc_atlas_endpoint_radial import (
    read_seed, panel_step_options, inward_panel_point, physical_field_comparison,
    geometric_mass_sources, centre_lapse_datum, integrate_independent_mass,
    compare_independent_mass, FULL_MASS_ABSOLUTE_TOLERANCE, LAPSE_ABSOLUTE_TOLERANCE,
    MASS_RATE_ABSOLUTE_TOLERANCE,
)
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


def fixture(root):
    path=root/'.build/results/stage.npz'
    path.parent.mkdir(parents=True)
    np.savez_compressed(path,x=np.array(.43),level=np.array(0),
        source_commit=np.array('a'*40),**{'cap.chi':np.array([.1,.2]),'wall.r':np.array(.5)})
    raw=path.read_bytes()
    return path,dict(level=0,source_commit='a'*40,failed_stage_seed=dict(
        path=path.relative_to(root).as_posix(),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),
        source_x=.43))


def test_authenticated_seed_recovers_exact_arrays(tmp_path):
    _,receipt=fixture(tmp_path)
    state,x=read_seed(receipt,tmp_path)
    assert x==.43 and float(state['wall']['r'])==.5
    np.testing.assert_array_equal(state['cap']['chi'],[.1,.2])


def test_changed_seed_or_metadata_is_rejected(tmp_path):
    path,receipt=fixture(tmp_path)
    receipt['level']=1
    with pytest.raises(ValueError,match='metadata'):read_seed(receipt,tmp_path)
    receipt['level']=0
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='bytes changed'):read_seed(receipt,tmp_path)
    receipt['failed_stage_seed']['path']='outside.npz'
    with pytest.raises(ValueError,match='ignored result directory'):read_seed(receipt,tmp_path)


def test_panel_step_check_preserves_tolerance_and_existing_smaller_step():
    options=dict(rtol=2.5e-14,atol=1e-15,first_step=.4)
    assert panel_step_options((1.,1.5),options,None)==options
    assert panel_step_options((1.,1.5),options,.5)==dict(
        rtol=2.5e-14,atol=1e-15,first_step=.25,max_step=.25)
    assert options['first_step']==.4
    assert panel_step_options((1.,1.5),dict(max_step=.1),.5)['max_step']==.1
    for fraction in (0,-1,1.1,np.inf,np.nan):
        with pytest.raises(ValueError,match='fraction'):
            panel_step_options((1.,1.5),options,fraction)
    with pytest.raises(ValueError,match='positive radial panel'):
        panel_step_options((1.,1.),options,.5)


def test_one_sided_diagnostic_moves_only_resolved_endpoint_samples():
    span=(1.,1.25)
    assert span[0]<inward_panel_point(span[0],span)<1.000000000001
    assert 1.249999999999<inward_panel_point(span[1],span)<span[1]
    assert inward_panel_point(1.125,span)==1.125
    with pytest.raises(ValueError,match='not resolved'):
        inward_panel_point(1.,(1.,np.nextafter(1.,2.)))


def _stage(chi, density, Z=Q(2, 5)):
    return dict(rates=dict(cap=dict(chi=np.asarray(chi), canonical_density=np.asarray(density))),
                connection=dict(cap=dict(wall_coordinate=float(Z), wall=dict(
                    f=1., p_trace=.4, h_trace=.3, q_trace=.2, velocity=.1, A=.9, B=.95)),
                    residual=dict(conformal_wall_flux=.2)))


def _exact_squared_field_integral(basis, coefficients, Z):
    # Exact rational polynomial products, independent of numerical Gram factors.
    knots = (Q(0), Q(1, 4), Q(1, 2))
    total = Q(0)
    for cell in range(2):
        h = knots[cell+1]-knots[cell]
        limit = min(Q(1), (Z-knots[cell])/h)
        poly = [sum(row[j]*Q(float(coefficients[i]))
                    for j, i in enumerate(basis.indices(cell)))
                for row in basis.exact_shapes(cell)]
        total += h*sum(a*b*limit**(j+k+1)/Q(j+k+1)
                       for j, a in enumerate(poly) for k, b in enumerate(poly))
    return float(total)


def test_reported_field_norm_matches_independent_polynomial_integral():
    budget = Budget(maximum_evaluations=1_000_000, cpu_seconds=15.)
    basis = FixedCapBasis([Q(0), Q(1, 4), Q(1, 2)], budget=budget)
    c = np.linspace(.1, .7, basis.dimension)
    other = c+np.linspace(-2e-6, 3e-6, basis.dimension)
    left, right = _stage(c, 2*c), _stage(other, 2*other)
    right['connection']['cap']['wall']['h_trace'] += 1e-8
    right['connection']['residual']['conformal_wall_flux'] += 2e-10
    record = physical_field_comparison(basis, left, right, budget)
    for key, multiplier in (('chi', 1.), ('canonical_density', 2.)):
        row = record['cap_rate_fields'][key]
        expected = np.sqrt(_exact_squared_field_integral(basis, multiplier*(c-other), Q(2, 5)))
        assert row['difference_galerkin_L2'] == pytest.approx(expected, rel=2e-13)
        expected_ref = np.sqrt(_exact_squared_field_integral(basis, multiplier*c, Q(2, 5)))
        assert row['reference_galerkin_L2'] == pytest.approx(expected_ref, rel=2e-13)
    assert record['wall_traces']['h_trace']['absolute_difference'] == pytest.approx(1e-8, abs=1e-16)
    assert record['conformal_wall_flux_absolute_difference'] == pytest.approx(2e-10, abs=1e-16)


def test_field_norm_is_invariant_under_basis_rescaling():
    budget = Budget(maximum_evaluations=1_000_000, cpu_seconds=15.)
    basis = FixedCapBasis([Q(0), Q(1, 4), Q(1, 2)], budget=budget)
    scale = np.geomspace(1e-3, 1e3, basis.dimension)
    c = np.linspace(.1, .7, basis.dimension)
    delta = np.linspace(-1e-7, 2e-7, basis.dimension)
    original = physical_field_comparison(basis, _stage(c, 2*c), _stage(c+delta, 2*(c+delta)), budget)

    class RescaledBasis:
        def gram(self, Z, *, budget):
            gram = basis.gram(Z, budget=budget)
            return SimpleNamespace(condition=gram.condition,
                squared_norm=lambda value, *, budget: gram.squared_norm(scale*value, budget=budget))

    rescaled = physical_field_comparison(RescaledBasis(), _stage(c/scale, 2*c/scale),
                                        _stage((c+delta)/scale, 2*(c+delta)/scale), budget)
    for key in ('chi', 'canonical_density'):
        a, b = original['cap_rate_fields'][key], rescaled['cap_rate_fields'][key]
        assert a['reference_galerkin_L2'] == pytest.approx(b['reference_galerkin_L2'], rel=1e-13)
        assert a['difference_galerkin_L2'] == pytest.approx(b['difference_galerkin_L2'], rel=1e-8)
        assert abs(a['coordinate_infinity']-b['coordinate_infinity']) > 1e-8


def test_zero_reference_field_has_no_fabricated_relative_norm():
    budget = Budget(maximum_evaluations=1_000_000, cpu_seconds=15.)
    basis = FixedCapBasis([Q(0), Q(1, 4), Q(1, 2)], budget=budget)
    zero = np.zeros(basis.dimension)
    record = physical_field_comparison(basis, _stage(zero, zero), _stage(zero, zero), budget)
    for field in record['cap_rate_fields'].values():
        assert field['difference_galerkin_L2'] == 0.
        assert field['difference_over_reference'] is None


def _polynomial_sources(c, a0, c1, a1, radius_of, derivative_of):
    def coefficients(panel, coordinate):
        radius = float(radius_of(coordinate))
        return dict(a=a0*radius, s=3*c*radius**2+a0*c*radius**4, da=a1*radius,
                    ds=3*c1*radius**2+(a0*c1+a1*c)*radius**4,
                    jacobian=float(derivative_of(coordinate)))
    return coefficients


def _boundary():
    return dict(panels=2, centre_lapse=0., centre_field=0., wall_p=.1, wall_h=-.2,
                reference_p_trace=.1, reference_h_trace=-.2)


def test_independent_integrator_recovers_split_panel_mass_rate_and_lapse():
    c, a0, c1, a1 = .2, .5, -.15, .3
    left, right, shift = .25, .8, .17
    spans = ((left, .41), (.41, .63), (.63, right))
    sources = _polynomial_sources(c, a0, c1, a1, lambda coordinate: coordinate, lambda coordinate: 1.)
    entering = dict(entering_mass=c*left**3, entering_rate=c1*left**3,
                    entering_lapse=a0*left**2/2+shift)
    previous = None
    for order in (32, 64):
        result = integrate_independent_mass(spans, sources, order, **entering)
        assert result['samples'] == order*len(spans)
        assert result['mass'] == pytest.approx(c*right**3, abs=1e-12)
        assert result['rate'] == pytest.approx(c1*right**3, abs=1e-12)
        assert result['lapse'] == pytest.approx(a0*right**2/2+shift, abs=1e-12)
        if previous is not None:
            assert result['mass']-previous['mass'] == pytest.approx(0., abs=1e-12)
            assert result['rate']-previous['rate'] == pytest.approx(0., abs=1e-12)
        previous = result


def test_independent_integrator_zero_source_limit_carries_entering_data():
    c, c1, shift = .4, -.25, .15
    left, right = .2, .9
    spans = ((left, .47), (.47, right))
    sources = _polynomial_sources(c, 0., c1, 0., lambda coordinate: coordinate, lambda coordinate: 1.)
    result = integrate_independent_mass(
        spans, sources, 32, entering_mass=c*left**3, entering_rate=c1*left**3, entering_lapse=shift)
    assert result['mass'] == pytest.approx(c*right**3, abs=1e-12)
    assert result['rate'] == pytest.approx(c1*right**3, abs=1e-12)
    assert result['lapse'] == pytest.approx(shift, abs=1e-14)


def test_independent_integrator_includes_coordinate_jacobian():
    c, a0, c1, a1 = .15, .35, .08, -.2
    # r=u^2, so the panel coordinate is not the areal radius.
    left_u, right_u = .3, .85
    left_r, right_r = left_u**2, right_u**2
    spans = ((left_u, .52), (.52, right_u))
    sources = _polynomial_sources(c, a0, c1, a1, lambda coordinate: coordinate**2,
                                 lambda coordinate: 2.*coordinate)
    result = integrate_independent_mass(
        spans, sources, 64, entering_mass=c*left_r**3, entering_rate=c1*left_r**3,
        entering_lapse=a0*left_r**2/2)
    assert result['mass'] == pytest.approx(c*right_r**3, abs=1e-12)
    assert result['rate'] == pytest.approx(c1*right_r**3, abs=1e-12)
    assert result['lapse'] == pytest.approx(a0*right_r**2/2, abs=1e-12)
    budget = Budget(maximum_evaluations=1000, cpu_seconds=15.)
    integrate_independent_mass(spans, sources, 32, budget=budget,
                              entering_mass=c*left_r**3, entering_rate=c1*left_r**3,
                              entering_lapse=a0*left_r**2/2)
    assert budget.evaluations == 32*len(spans)


def test_geometric_sources_cover_zero_and_finite_alpha():
    radius, field, momentum, normal = .37, .4, -.2, .15
    df, dp, dh = .3, -.1, .05
    zero = geometric_mass_sources(radius, field, momentum, normal, df, dp, dh, 0.)
    assert (zero['a'], zero['da'], zero['ds']) == (0., 0., 0.)
    assert zero['s'] == pytest.approx(1.5*radius**2)
    alpha = .25
    base = geometric_mass_sources(radius, field, momentum, normal, df, dp, dh, alpha)
    square = momentum**2+normal**2
    argument = math.sqrt(alpha/6)*field
    assert base['a'] == pytest.approx(alpha*square/(2*radius))
    assert base['s'] == pytest.approx(1.5*radius**2*math.cosh(argument)**4+alpha*square/4)
    # A smaller step loses the source difference to cancellation near s=0.2.
    step = 1e-5
    plus = geometric_mass_sources(radius, field+step*df, momentum+step*dp, normal+step*dh, 0., 0., 0., alpha)
    minus = geometric_mass_sources(radius, field-step*df, momentum-step*dp, normal-step*dh, 0., 0., 0., alpha)
    assert base['da'] == pytest.approx((plus['a']-minus['a'])/(2*step), rel=1e-8, abs=1e-10)
    assert base['ds'] == pytest.approx((plus['s']-minus['s'])/(2*step), rel=1e-8, abs=1e-10)
    assert centre_lapse_datum(field, alpha) == pytest.approx(-math.log(math.cosh(argument)))
    assert centre_lapse_datum(field, 0.) == 0.


def _measurement(order, **values):
    base = dict(full_mass=1., lapse=.2, fixed_radius_mass_rate=.3,
                moving_wall_mass_rate=.4, temporal_mass_residual=.5, samples=order)
    base.update(values)
    base['order'] = order
    return base


def test_mass_lapse_comparison_keeps_failures_and_does_not_reuse_value_tolerance_on_rates():
    reference = _measurement(0)
    close_rate = _measurement(32, fixed_radius_mass_rate=reference['fixed_radius_mass_rate']+5e-9)
    other = _measurement(64, fixed_radius_mass_rate=reference['fixed_radius_mass_rate']+5e-9)
    passed = compare_independent_mass(reference, (close_rate, other), _boundary())
    assert passed['outcome'] == 'passed'
    assert passed['failed_channels'] == []
    assert passed['order_differences']['fixed_radius_mass_rate'] == pytest.approx(0., abs=1e-15)
    assert passed['tolerances']['full_mass_absolute'] == FULL_MASS_ABSOLUTE_TOLERANCE == 1e-10
    assert passed['tolerances']['lapse_absolute'] == LAPSE_ABSOLUTE_TOLERANCE == 1e-10
    assert passed['tolerances']['physical_directional_and_temporal_mass_rate_absolute'] == (
        MASS_RATE_ABSOLUTE_TOLERANCE) == 1e-8

    failed_rate = _measurement(32, temporal_mass_residual=reference['temporal_mass_residual']+2e-8)
    failed_mass = _measurement(64, full_mass=reference['full_mass']+2e-10)
    failed = compare_independent_mass(reference, (failed_rate, failed_mass), _boundary())
    assert failed['outcome'] == 'failed'
    assert failed['failed_channels'] == ['order_32.temporal_mass_residual', 'order_64.full_mass']
    assert failed['orders'][0]['differences']['temporal_mass_residual'] == pytest.approx(2e-8)
    assert failed['orders'][1]['differences']['full_mass'] == pytest.approx(2e-10)
    assert failed['order_differences']['full_mass'] == pytest.approx(2e-10)
    assert failed['tolerances']['full_mass_absolute'] == 1e-10
    assert failed['tolerances']['physical_directional_and_temporal_mass_rate_absolute'] == 1e-8
