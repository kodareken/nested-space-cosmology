import math
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_scalar_core import solve_core,evolve_core,CoolingRun,CoreAtlasEmbedding
from recursive_horizons.nsc_atlas_shared_scalar_slice import SharedScalarSlice
from recursive_horizons.nsc_atlas_flux_conversion import ScalarFlux


@pytest.fixture(scope='module')
def profile():return solve_core()


def test_self_consistent_profile_and_independent_integrals(profile):
    r=profile.summary()
    assert r['central_amplitude']==pytest.approx(1.0214930362,abs=2e-8)
    assert r['radial_mass_norm']==pytest.approx(3.5053297033,abs=2e-8)
    assert abs(r['mass_integral_gap'])<2e-8
    assert abs(r['virial_residual'])<2e-8
    assert abs(r['eigenvalue_identity'])<2e-8
    assert r['mass_fraction_radii']['0.99']==pytest.approx(5.6746887148,abs=2e-7)


def test_frequency_is_not_selected_by_the_initial_atlas():
    a=SharedScalarSlice();b=SharedScalarSlice(scalar=ScalarFlux(omega=2000.))
    np.testing.assert_allclose(a.densities,b.densities,rtol=1e-14,atol=1e-14)
    np.testing.assert_allclose(a.expansions,b.expansions,rtol=1e-14,atol=1e-14)
    assert b.field==pytest.approx(a.field/100)
    assert math.sqrt(b.scalar.vacuum(b.parent_flux)/3)/(2000*.1**2)<.02


def test_free_envelope_recovers_the_analytic_dispersing_gaussian():
    c=CoolingRun(points=1023,step=.04,end=12,h=0,gravity=False,absorber=False)
    result,data=evolve_core(c,retain=True)
    r=data['radius'];sigma=c.width
    norm=solve_core().mass_norm
    A=math.sqrt(norm/(math.sqrt(math.pi)*sigma**3/4))
    exact=A*sigma**3/(sigma*sigma+1j*c.end)**1.5*np.exp(-r*r/(2*(sigma*sigma+1j*c.end)))
    assert np.max(abs(data['complex_radial_field']/r-exact))<2e-10
    assert abs(result['mass_balance'])<2e-10


def test_stationary_profile_does_not_cool_through_an_interior_sink():
    result=evolve_core(CoolingRun(shape='stationary',points=1023,step=.04,end=20,h=0,absorber=False))
    assert abs(result['final_mass']-result['initial_mass'])<1e-10
    assert result['removed_mass']==0
    assert result['late_profile_gap']<.003
    assert abs(result['energy_balance'])<2e-5


def test_core_can_be_embedded_without_discarding_ambient_interference():
    r=CoreAtlasEmbedding().solve()
    assert r['geometric_core_mass']>r['leading_isolated_mass']
    assert r['geometric_mass_parts']['inherited_cross_term']>1.8e-5
    assert abs(r['independent_field_mass_gap'])<2e-11
    assert abs(r['independent_mass_gap'])<2e-11
    assert abs(r['independent_radius_gap'])<2e-10
    assert r['maximum_sampled_Hamilton_residual']<1e-7
    assert r['minimum_core_beta_squared']>.98
    assert r['proper_slice_energy_excess']>r['geometric_mass_in_energy_units']
    assert r['fixed_middle_mass']==.12
    assert r['middle_constraint_residual']<1e-12
    assert max(abs(x) for x in r['middle_seams'].values())<1e-12
    for key in ('angular_residual','temporal_residual','force_residual','energy_residual'):
        assert abs(r['outer_wall'][key])<1e-12


def test_core_taper_adds_no_material_shell():
    e=CoreAtlasEmbedding();R=e.cutoff_radius
    assert e.field(R)==(0.,0.)
    assert e.field(2*R)==(0.,0.)
    assert abs(e.field(R*(1-1e-5))[1])<1e-10


def test_zero_core_recovers_the_actual_reference_slice():
    e=CoreAtlasEmbedding(with_core=False);r=e.solve();s=e.reference
    assert r['geometric_core_mass']==0
    assert r['proper_slice_energy_excess']==0
    assert r['outer_wall_speed_change']==pytest.approx(0,abs=1e-14)
    assert r['parent_pole_to_wall_length']==pytest.approx(s.cap_geometry['parent_length'],abs=2e-9)
    assert r['total_pole_to_pole_length']==pytest.approx(s.lengths_and_volumes['total_length'],abs=2e-9)


def test_relative_phase_changes_parent_mass_without_changing_isolated_core_energy():
    a=CoreAtlasEmbedding(weak_amplitude=.0826621713254).solve()
    b=CoreAtlasEmbedding(weak_amplitude=.0826621713254,phase_sign=-1).solve()
    assert a['leading_isolated_mass']==b['leading_isolated_mass']
    assert a['geometric_mass_parts']['inherited_cross_term']==pytest.approx(-b['geometric_mass_parts']['inherited_cross_term'])
    assert a['geometric_core_mass']>1.8*b['geometric_core_mass']>0
    assert b['maximum_sampled_Hamilton_residual']<1e-7
    assert abs(b['outer_wall']['force_residual'])<1e-12
