import numpy as np
import pytest

from recursive_horizons.nsc_atlas_parent_source import (
    ParentScalarSource, uniform_excited_parent, planar_interface_control, planar_child_energy_integral,
    planar_volume_energies)


@pytest.fixture(scope='module')
def source():
    model=ParentScalarSource()
    background,_=model.prepared.solve(tolerance=1e-9)
    solved,_=model.solve(background,tolerance=1e-9)
    return model,background,solved


def test_zero_source_recovers_regular_geometry_and_no_wall_flux(source):
    _,bg,_=source
    model=ParentScalarSource(amplitude=0.)
    sol,_=model.solve(bg)
    np.testing.assert_allclose(sol.y[[0,2]],1.,atol=1e-12)
    row=model.wall(bg,sol)
    assert row['radius']==pytest.approx(bg.y[0,-1],abs=1e-10)
    assert row['inside']['wall_areal_acceleration']==pytest.approx(1/row['radius'],abs=1e-10)
    assert abs(row['inside']['wall_field_acceleration'])<1e-10


def test_source_gravity_satisfies_physical_constraints(source):
    model,bg,sol=source
    d=model.diagnostics(bg,sol)
    assert .95<d['minimum_factor']<1<d['maximum_factor']<1.05
    assert max(d['weighted_Hamilton_residual'])<1e-6
    assert max(d['Hamilton_residual'])<1e-6
    assert np.max(abs(sol.y[[1,3],0]))<1e-12


def test_energy_and_mass_are_checked_in_the_sourced_metric(source):
    model,bg,sol=source
    e=model.energy(bg,sol)
    assert e['positive_quadratic']>1.
    assert e['total_scalar_excess']==pytest.approx(e['positive_quadratic']+e['interference'])
    for row in model.mass_checks(bg,sol):
        assert row['geometric_mass']>0
        assert abs(row['gap'])<1e-7


def test_complete_instantaneous_junction_and_two_sided_response(source):
    model,bg,sol=source
    row=model.wall(bg,sol)
    for key,value in row.items():
        if key.endswith('_gap'):assert abs(value)<1e-10,key
    assert row['radius']>bg.y[0,-1]
    assert row['inside']['normal_field']==pytest.approx(row['outside']['normal_field'])
    assert row['inside']['wall_field_acceleration']>.1
    assert row['inside']['outward_flux_slope']>0  # does not establish inward matter deposition


def test_frozen_layer_can_pass_gravity_but_fail_scalar_matching(source):
    model,bg,_=source
    frozen,_=model.solve(bg,match_scalar_jet=False)
    d=model.diagnostics(bg,frozen)
    row=model.wall(bg,frozen)
    assert max(d['weighted_Hamilton_residual'])<1e-6
    assert abs(row['angular_Israel_gap'])<1e-10
    assert abs(row['field_acceleration_gap'])>.3
    assert row['field_acceleration_gap']==pytest.approx(row['scalar_second_jet_gap'],abs=1e-10)


def test_added_layer_preserves_scalar_value_and_first_jet(source):
    model,bg,sol=source
    L=bg.p[1]
    f,fp=model.packet(L,L)
    g,gp=model.collar(L,L,sol.p[0])
    assert f==fp==g==gp==0
    g,gp=model.collar(L-model.collar_width,L,sol.p[0])
    assert abs(g)+abs(gp)<1e-25
    eps=1e-7
    g,_=model.collar(L-eps,L,sol.p[0])
    assert 2*g/eps**2==pytest.approx(sol.p[0],rel=1e-5)


def test_surface_mass_relation_uses_binding_and_both_normals(source):
    model,bg,sol=source
    row=model.wall(bg,sol)
    mean_beta=(row['inside']['beta']+row['outside']['beta'])/2
    expected=model.prepared.k*row['radius']**2*mean_beta
    actual=row['outside']['Misner_Sharp_mass']-row['inside']['Misner_Sharp_mass']
    assert actual==pytest.approx(expected,abs=1e-11)


def test_uniform_excited_parent_closes_exact_initial_constraints():
    row=uniform_excited_parent()
    assert row['radius']**2==pytest.approx(1800/467,abs=1e-12)
    for key in ('angular_Israel_gap','mean_force_gap','child_constraint_gap','parent_constraint_gap'):
        assert abs(row[key])<1e-12
    assert row['child_density']==pytest.approx(1/8)
    assert row['parent_initial_density']==pytest.approx(3/8)
    assert row['parent_excitation_density']==pytest.approx(3/32)
    # It starts at a turning radius with positive acceleration, not static equilibrium.
    assert row['initial_wall_areal_acceleration']==pytest.approx(1/row['radius'])
    assert row['scalar_initial_C2_gap']<0  # a declared weak initial corner


def test_planar_retarded_volume_energy_matches_independent_boundary_flux():
    S=uniform_excited_parent()['field_force_jump']
    boundary=planar_interface_control(.01,force_jump=S,omega=20.)
    volume=planar_child_energy_integral(.01,force_jump=S,omega=20.)
    assert volume==pytest.approx(boundary['child_energy_per_area'],abs=1e-14)
    assert boundary['inward_flux']>0
    assert boundary['wall_time_derivative']<0
    assert boundary['wall_normal_derivative']<0
    energies=planar_volume_energies(.01,force_jump=S,omega=20.)
    assert energies['parent_excess']<0
    assert abs(energies['net'])<1e-14


def test_massive_planar_control_recovers_tangent_flux_and_pressure_limits():
    S=uniform_excited_parent()['field_force_jump']
    coarse=planar_interface_control(.01,force_jump=S,omega=20.)
    fine=planar_interface_control(.001,force_jump=S,omega=20.)
    limit=S*S/4
    assert abs(fine['inward_flux']/.001**2-limit)<abs(coarse['inward_flux']/.01**2-limit)/90
    assert fine['pressure_drop']/fine['inward_flux']==pytest.approx(1.,rel=7e-5)


def test_zero_scalar_mixing_removes_transfer_without_removing_the_wall():
    from recursive_horizons.nsc_atlas_scalar_preparation import ScalarPreparation
    from recursive_horizons.nsc_atlas_flux_conversion import ScalarFlux
    row=uniform_excited_parent(ScalarPreparation(scalar=ScalarFlux(epsilon=1.)))
    assert row['radius']>0
    assert row['pressure']>0
    assert row['parent_excitation_density']==0
    assert row['conditional_inward_flux_over_time_squared']==0
