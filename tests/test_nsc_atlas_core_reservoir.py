import math
import numpy as np
import pytest
from recursive_horizons.nsc_atlas_core_reservoir import CoreReservoir,compact_bump,mass_transport


@pytest.fixture(scope='module')
def compensated():return CoreReservoir().solve()


def test_local_mass_is_supplied_without_negative_absolute_energy(compensated):
    s=compensated;r=s.model.radius
    assert s.depletion==pytest.approx(.1859136874,abs=2e-8)
    mass=s.solution.sol(r*np.arange(5))[0]
    assert mass[1]>9.5e-5
    assert mass[2]<0
    assert mass[4]-mass[3]>7e-7
    assert abs(mass[4])<2e-11
    row=s.state(np.linspace(0,s.model.outer,501))
    assert np.min(row['density'])>=s.model.reference.scalar.vacuum(s.model.reference.parent_flux)


def test_independent_geometric_constraints_and_initial_walls(compensated):
    r=compensated.summary()
    assert r['maximum_sampled_Hamilton_residual']<1e-7
    assert r['maximum_sampled_momentum_residual']<1e-12
    assert abs(r['independent_mass_gap'])<2e-11
    assert abs(r['independent_radius_gap'])<2e-10
    assert abs(r['mass_integral_gap'])<2e-11
    assert abs(r['outer_wall_speed_gap'])<2e-12
    assert abs(r['outer_wall_acceleration_gap'])<2e-12
    assert r['proper_energy_excess']>1e-6  # Proper energy is not geometric mass.
    causal=r['protected_reference_annulus']
    assert .85<causal['fraction_of_actual_initial_parent_volume']<.91
    assert causal['permanent_inner_angle']>causal['initial_support_angle']
    assert causal['permanent_outer_angle']>causal['permanent_inner_angle']
    assert r['vacuum_Lambda_reference_radius_squared']==pytest.approx(12.3482523989285,abs=2e-11)


def test_outward_characteristic_current_and_pressure_work(compensated):
    R=np.linspace(3*compensated.model.radius,compensated.model.outer,501)[1:-1]
    r=compensated.state(R)
    assert np.min(r['outward_energy_current'])>=0
    assert np.max(r['outward_energy_current'])>2e-4
    np.testing.assert_allclose(r['Pi']+r['beta']*r['field_derivative'],0,atol=1e-14)
    assert np.max(r['fixed_radius_rate'])<1e-15
    assert np.min(r['fixed_radius_rate'])< -3e-5
    # Keeping the old K=-H h would miss the nonzero pulse momentum.
    missing=compensated.model.reference.scalar.kappa*r['Pi']*r['beta']*r['field_derivative']
    assert np.max(abs(missing))>1e-4


def test_mass_law_recovers_vacuum_work_and_null_outflow():
    a=mass_transport(2,1,.3,.9,0,0,.5)
    assert a['normal_mass_derivative']==pytest.approx(.3)
    assert abs(a['fixed_radius_rate'])<1e-15
    F=.2;v=.3;b=.9
    a=mass_transport(2,1,v,b,-b*F,F,.5)
    assert a['fixed_radius_rate']==pytest.approx(-2*b*F*F*(b+v)**2)


def test_zero_state_recovers_the_original_parent():
    m=CoreReservoir(core_enabled=False,pulse_fraction=0)
    r=m.solve()
    assert r.depletion==0
    assert abs(r.solution.y[0,-1])<1e-15
    R=np.linspace(0,m.outer,501)
    row=r.state(R)
    np.testing.assert_allclose(row['beta']**2,1-m.reference.cap_geometry['parent_curvature']*R**2,atol=1e-14)
    np.testing.assert_allclose(row['radial_K'],-m.H,atol=1e-14)


def test_field_matching_and_finite_admission(compensated):
    c=compensated.model
    assert compact_bump(0,1,2)==(0.,0.)
    assert compact_bump(1.5,1,2)==(1.,0.)
    assert compact_bump(3,1,2)==(0.,0.)
    assert c.fields(c.outer,compensated.depletion)[:2]==(0.,0.)
    assert c.fields(c.outer+1,compensated.depletion)[:2]==(0.,0.)
    with pytest.raises(ValueError):CoreReservoir(pulse_fraction=-1)


def test_symbolic_constraint_and_flux_identities():
    from scripts.check_nsc_atlas_core_reservoir import identities
    assert set(identities().values())=={'0'}
