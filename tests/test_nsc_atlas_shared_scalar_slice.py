import math

import numpy as np
import pytest
from scipy.integrate import quad

from recursive_horizons.nsc_atlas_shared_scalar_slice import SharedScalarSlice,alternative_floor_test


@pytest.fixture(scope='module')
def atlas():
    return SharedScalarSlice()


def test_frozen_state_is_the_same_action_not_a_changed_bare_term(atlas):
    assert atlas.scalar.bare_density==0
    assert atlas.densities==pytest.approx([1/32,7/32,21/32])
    for q,r in zip(atlas.fluxes,atlas.regions):
        assert atlas.reference_law.vacuum(q+atlas.scalar.mu*atlas.field)==pytest.approx(r.Lambda)


def test_full_moving_wall_jets_and_shared_mass(atlas):
    assert atlas.reference_walls[0].outside is atlas.reference_walls[1].inside
    for w in atlas.walls:
        for name in ('angular_residual','temporal_residual','force_residual','reference_energy_residual'):
            assert abs(w[name])<2e-14
        for side in ('inside','outside'):
            s=w[side]
            assert abs(s['timelike_residual'])<1e-14
            assert abs(s['temporal_derivative_identity_residual'])<1e-13
            assert s['independent_areal_acceleration']==pytest.approx(w['acceleration'],abs=1e-13)
            assert s['mass_offset']==pytest.approx(atlas.regions[s['region']].mass,abs=1e-14)
    assert atlas.walls[0]['velocity']>1  # areal speed, while normalization is timelike
    assert atlas.walls[0]['outside']['mass_offset']==pytest.approx(.12)
    assert atlas.walls[1]['inside']['mass_offset']==pytest.approx(.12)
    w=atlas.walls[0];s=w['inside'];R=w['radius']
    resting_mass=R*(1-s['beta']**2)/2
    wrong=-resting_mass/R**2+atlas.scalar.kappa*R*atlas.densities[0]/2+s['A']*s['beta']
    assert abs(wrong-w['acceleration'])>.5


def test_bulk_surface_constraints_and_independent_connector(atlas):
    c=atlas.checks()
    assert c['middle_Hamilton_max']<1e-12
    assert c['middle_momentum_max']<1e-12
    assert max(map(abs,c['cap_Hamilton_residuals']))<1e-12
    assert all(abs(v)<1e-12 for row in c['seams'] for v in row.values())
    assert c['middle_f_lower_bound']>.6
    assert abs(c['independent_radius_endpoint_gap'])<1e-10
    assert abs(c['independent_volume_gap'])<1e-9
    assert c['Hamilton_surface_coefficient']==c['required_surface_coefficient']
    # A complex-step derivative catches an incorrect v' that could otherwise
    # cancel out of a constraint computed twice with the same wrong derivative.
    R=np.linspace(1.02,1.88,17)
    vprime=atlas.middle_velocity(R+1e-30j)[0].imag/1e-30
    np.testing.assert_allclose(vprime,atlas.middle_velocity(R)[1],atol=1e-13,rtol=1e-13)


def test_regular_poles_and_compact_parent_cap(atlas):
    c=atlas.cap_geometry;l=atlas.lengths_and_volumes
    assert c['inner_curvature']<0<c['parent_curvature']
    root=math.sqrt(c['parent_curvature'])
    assert math.pi/2<root*c['parent_length']<math.pi
    assert math.sin(root*c['parent_length'])/root==pytest.approx(atlas.outer_radius)
    assert math.cos(root*c['parent_length'])==pytest.approx(-atlas.walls[1]['outside']['beta'])
    assert l['total_length']==pytest.approx(20.39342299980659,abs=1e-10)
    assert all(x>0 for x in l['volumes'])


def test_positive_lapse_realizes_both_temporal_junctions(atlas):
    lo,hi=atlas.checks()['lapse_global_bounds']
    ell=np.linspace(0,atlas.lengths_and_volumes['total_length'],2001)
    assert np.min(atlas.lapse(ell))>=lo>0
    assert np.max(atlas.lapse(ell))<=hi
    assert atlas.lapse(ell[[0,-1]]).tolist()==[1.,1.]
    h=1e-6
    for x,w in zip(atlas.lengths_and_volumes['wall_positions'],atlas.walls):
        assert atlas.lapse(x)==1
        left=-math.log(float(atlas.lapse(x-h)))/h
        right=math.log(float(atlas.lapse(x+h)))/h
        assert left==pytest.approx(w['inside']['A'],abs=1e-7)
        assert right==pytest.approx(w['outside']['A'],abs=1e-7)


def test_one_sided_observers_recover_the_common_field(atlas):
    for j in (0,1):
        row=atlas.observer(j)
        assert row['scalar_value']==pytest.approx(atlas.field,abs=1e-14)
        assert row['scalar_time_derivative']==0
        assert abs(row['normal_square'])<1e-13
        assert abs(row['tidal_residual'])<1e-13
        assert row['inferred_vacuum']==pytest.approx(atlas.scalar.vacuum(atlas.fluxes[j]),abs=1e-13)
    assert atlas.walls[1]['pressure']-atlas.walls[0]['pressure']==pytest.approx(atlas.charge**2)


def test_conditional_pulses_keep_pressure_work_and_initial_offsets(atlas):
    a,b=(atlas.tangent_coefficients(j) for j in (0,1))
    assert [a['inward_flux_coefficient'],b['inward_flux_coefficient']]==pytest.approx([18.75,56.25])
    assert a['sides'][1]['initial_mass_offset']==b['sides'][0]['initial_mass_offset']==.12
    assert a['sides'][1]['initial_density_subtracted_mass_coefficient']<0
    assert b['sides'][0]['initial_density_subtracted_mass_coefficient']<0


def test_alternative_floor_excludes_every_canonical_state_on_this_geometry(atlas):
    result=alternative_floor_test()
    required=3*(atlas.cap_geometry['inner_curvature']+atlas.expansions[0]**2)/atlas.scalar.kappa
    assert required==pytest.approx(1/32,abs=1e-14)
    # Nearest integer to 3/5 is uniquely 1; arbitrary canonical gradients and
    # momenta can only raise, never remove, this all-sector energy floor.
    assert (25/32)*(3/5-1)**2==pytest.approx(1/8)
    assert result['required_grandchild_density']==pytest.approx(required)
    assert result['density_gap']==pytest.approx(3/32)


def test_parent_future_bound_by_independent_quadrature(atlas):
    d=atlas.parent_causal_domain();a0=d['initial_scale_factor']
    hmin=math.sqrt(atlas.scalar.kappa*d['ground_density']/3)
    integral=quad(lambda a:1/(a*math.sqrt(hmin*hmin*a*a-1)),a0,np.inf,epsabs=1e-12)[0]
    assert integral==pytest.approx(d['future_conformal_distance_upper_bound'],abs=1e-12)
    assert d['H_positive_lower_bound']>0
    assert .907<d['protected_initial_parent_volume_fraction_lower_bound']<.908


def test_parent_evolution_closes_in_two_clocks(atlas):
    a=atlas.evolve_parent();b=atlas.evolve_parent(clock='efolds')
    for key in ('proper_time','H','scalar_displacement','scalar_velocity','excitation_density','conformal_distance'):
        assert a[key]==pytest.approx(b[key],abs=1e-9)
    for r in (a,b):
        assert r['max_Friedmann_residual']<1e-9
        assert r['max_scalar_energy_balance_residual']<1e-9
        assert r['minimum_sampled_H']>atlas.parent_causal_domain()['H_positive_lower_bound']
        assert r['Omega_vacuum']+r['Omega_scalar_excitation']+r['Omega_spatial_curvature']==pytest.approx(1,abs=1e-9)
        assert 0<r['excitation_density']<r['initial_excitation']
    assert a['proper_time']==pytest.approx(2.8888455131,abs=1e-9)


def test_declared_domain_is_enforced(atlas):
    with pytest.raises(ValueError):atlas.middle(.9)
    with pytest.raises(ValueError):atlas.evolve_parent(efolds=2)
    with pytest.raises(ValueError):atlas.evolve_parent(clock='unrecognized')
