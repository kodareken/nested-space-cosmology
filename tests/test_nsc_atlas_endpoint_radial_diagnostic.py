"""A one-stage radial diagnosis authenticates its exact local state."""
import hashlib
from fractions import Fraction as Q
from types import SimpleNamespace
import numpy as np
import pytest

from diagnose_nsc_atlas_endpoint_radial import (
    read_seed, panel_step_options, inward_panel_point, physical_field_comparison,
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
