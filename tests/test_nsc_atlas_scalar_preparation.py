import math
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_scalar_preparation import ScalarPreparation
from recursive_horizons.nsc_atlas_flux_conversion import ScalarFlux,isolated_bounce
from recursive_horizons.nsc_atlas_flux_transition import FluxLaw


@pytest.fixture(scope="module")
def prepared():
    model=ScalarPreparation()
    solution,_=model.solve(tolerance=1e-9)
    return model,solution


def test_all_wall_conditions_use_the_same_field_and_opposite_normals(prepared):
    model,solution=prepared
    d=model.diagnostics(solution)
    for name in ("force_residual","full_Israel_component_residual","scalar_value_jump",
                 "scalar_normal_jump","max_boundary_residual"):
        assert abs(d[name])<1e-10
    assert d["child_normal_derivative"]==pytest.approx(d["parent_normal_derivative"])
    assert max(d["Hamilton_max"])<1e-9
    assert max(d["weighted_tangential_Einstein_max"])<3e-8
    relaxed=model.parent_vacuum-model.child_vacuum
    assert abs(relaxed-d["wall_pressure"])>.001


def test_regular_poles_and_ordinary_cap_branch(prepared):
    _,s=prepared
    assert np.max(np.abs(s.y[[0,3,4,7],0]))<1e-12
    np.testing.assert_allclose(s.y[[1,5],0],[1,1],atol=1e-12)
    assert s.y[1,-1]>0>s.y[5,-1]
    x=np.linspace(.0001,1,501)
    assert np.all(s.sol(x)[[0,4]]>0)


def test_zero_mixing_recovers_the_imported_pure_flux_bounce():
    model=ScalarPreparation(scalar=ScalarFlux(epsilon=1.))
    solution,_=model.solve(tolerance=1e-8)
    expected=isolated_bounce(FluxLaw(),1.5)
    assert solution.y[0,-1]==pytest.approx(expected["radius"],abs=2e-9)
    assert np.max(np.abs(solution.y[[2,3,6,7]]))==0
    assert model.scalar_reduction(solution)["field_layer_energy_per_wall_area"]==0


def test_physical_Schur_map_and_bulk_action_are_independent_checks(prepared):
    model,s=prepared
    row=model.scalar_reduction(s)
    assert min(row["child_DtN"],row["parent_DtN"])>0
    assert row["wall_displacement_reconstruction_gap"]<1e-13
    assert abs(row["bulk_boundary_action_gap"])<1e-9
    assert row["field_layer_energy_per_wall_area"]>0


def test_log_transfer_resolves_the_small_child_state_without_subtraction(prepared):
    model,s=prepared
    t=model.pole_transfer(s)
    assert t["log_gain"]>50
    assert t["pole_displacement"]<0
    assert abs(t["logarithmic_derivative_gap"])<1e-8
    assert 1e-48<t["scalar_energy_at_pole"]<1e-45
    assert t["scalar_energy_at_pole"]/model.scalar.prepare()["oscillation_energy"]<1e-43


def test_Lorentzian_wall_is_timelike_and_has_no_scalar_rest_flux(prepared):
    model,s=prepared
    R=s.y[0,-1]
    row=model.lorentzian_wall(R,1.)
    assert row["areal_radius"]>R
    assert row["radius_proper_speed"]>1  # areal-radius derivative is not a local speed
    assert row["scalar_rest_frame_flux"]==0
    eta=.3
    # Radial wall motion in its intrinsic dS3 metric: gamma_etaeta=-R^2.
    assert (-R*R)*(1/R)**2==pytest.approx(-1.)
    assert model.lorentzian_wall(R,eta)["radius_proper_acceleration"]==pytest.approx(math.cosh(eta)/R)


def test_open_child_continuation_and_energy_bound(prepared):
    model,s=prepared
    chi=model.pole_transfer(s)["pole_displacement"]
    row=model.child_response(chi)
    assert row["spatial_curvature"]==-1
    assert 0<row["excitation_density"]<row["initial_excitation_density"]
    assert row["max_extra_H_squared_at_fixed_a"]==pytest.approx(row["initial_excitation_density"]/3,rel=1e-14,abs=0)
    assert row["vacuum_scale_factor"]>row["proper_time"]
