"""Checks for the sourced cap, distributional wall and regular vacuum bridge."""

import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_connection_slice import (
    connection_slice_record, schwarzschild_bridge, smooth_step,
    transition_proper_distance, vacuum_transition, wall_slice_data,
)


def test_smooth_step_has_flat_endpoints_and_correct_reflection():
    x = np.array([-.1, 0, .1, .2, .5, .8, .9, 1, 1.1])
    s, derivative = smooth_step(x)
    reverse, _ = smooth_step(1-x)
    np.testing.assert_allclose(s+reverse, 1, atol=2e-15)
    assert np.all(derivative >= 0)
    assert np.all(derivative[[0, 1, -2, -1]] == 0)
    assert float(smooth_step(.5)[0]) == .5


def test_wall_boost_reconstructs_both_imported_signed_normals():
    wall = wall_slice_data()
    gamma, gv = wall["gamma"], wall["gamma_v"]
    assert gamma*gamma-gv*gv == pytest.approx(1, abs=1e-14)
    assert gamma >= 1 and 0 < wall["v"] < 1
    for side in ("minus", "plus"):
        Rprime = wall["Rprime_"+side]
        V = wall["V_"+side]
        assert gamma*V+gv*Rprime == pytest.approx(wall["Rdot"], abs=2e-14)
        assert gv*V+gamma*Rprime == pytest.approx(wall["beta_"+side], abs=2e-14)
    assert wall["Rprime_plus"] < 0


def test_shell_projection_requires_delta_jacobian_and_matches_constraints():
    wall = wall_slice_data()
    sigma = wall["parameters"].sigma
    E, J = wall["energy_surface"], wall["momentum_surface"]
    assert E == pytest.approx(sigma*wall["gamma"])
    assert J == pytest.approx(sigma*wall["gamma_v"])
    assert E >= abs(J)
    assert -4*(wall["Rprime_plus"]-1)/wall["R"] == pytest.approx(16*math.pi*E)
    assert 2*(wall["V_plus"]-wall["V_minus"])/wall["R"] == pytest.approx(8*math.pi*J)
    # Forgetting the restricted delta or using rest energy changes the data.
    assert abs(E-sigma*wall["gamma"]**2) > 1e-4
    assert abs(E-sigma) > 1e-4


def test_wrong_exterior_orientation_fails_shell_energy_constraint():
    wall = wall_slice_data()
    wrong_jump = abs(wall["Rprime_plus"])-wall["Rprime_minus"]
    residual = -4*wrong_jump/wall["R"]-16*math.pi*wall["energy_surface"]
    assert abs(residual) > 1


def test_bulk_satisfies_vacuum_constraints_and_invariant_mass():
    wall = wall_slice_data()
    fields = vacuum_transition(np.linspace(1, 1.5, 101), wall)
    assert np.max(abs(fields["hamiltonian_residual"])) < 1e-12
    assert np.max(abs(fields["momentum_residual"])) < 1e-12
    np.testing.assert_allclose(fields["mass"], .4, atol=3e-15)
    assert np.all(fields["Rprime"] < 0)


def test_induced_schwarzschild_metric_is_spacelike():
    wall = wall_slice_data()
    fields = vacuum_transition(np.linspace(1, 1.5, 53), wall)
    np.testing.assert_allclose(fields["embedding_spatial_norm"], 1, atol=5e-14)
    np.testing.assert_allclose(fields["embedding_future_normal_norm"], -1, atol=5e-14)


def test_cutoff_is_slicing_freedom_with_the_same_mass_and_wall():
    original = wall_slice_data()
    for R0 in (.95, 1., 1.1):
        wall = dict(original, taper_stop=R0)
        radii = np.linspace(R0, wall["R"], 37)
        fields = vacuum_transition(radii, wall)
        np.testing.assert_allclose(fields["mass"], .4, atol=3e-15)
        assert fields["V"][0] == fields["V_R"][0] == 0
        assert fields["V"][-1] == pytest.approx(original["V_plus"])
        assert wall["parameters"] == original["parameters"]
        assert np.max(abs(fields["hamiltonian_residual"])) < 1e-11


def test_transition_matches_wall_and_time_symmetric_bridge():
    wall = wall_slice_data()
    fields = vacuum_transition(np.array([wall["R"], wall["taper_stop"]]), wall)
    assert fields["Rprime"][0] == pytest.approx(wall["Rprime_plus"], abs=3e-15)
    assert fields["V"][0] == pytest.approx(wall["V_plus"])
    assert fields["Kll"][-1] == fields["Ktheta"][-1] == 0
    u = -math.acosh(math.sqrt(wall["taper_stop"]/(2*wall["parameters"].mass)))
    bridge = schwarzschild_bridge(u, wall["parameters"].mass)
    for key in ("R", "Rprime", "Rsecond"):
        assert fields[key][-1] == pytest.approx(bridge[key], abs=3e-15)


def test_bridge_is_regular_through_minimum_and_has_the_same_mass():
    bridge = schwarzschild_bridge(np.array([-1e-4, 0, 1e-4, 2]), .4)
    assert bridge["R"][1] == pytest.approx(.8)
    assert bridge["Rprime"][1] == 0
    assert bridge["Rsecond"][1] == pytest.approx(1/(4*.4))
    assert np.all(np.diff(bridge["l"]) > 0)
    assert np.all(np.isfinite(bridge["scalar_curvature"]))
    np.testing.assert_allclose(bridge["mass"], .4, atol=5e-15)
    np.testing.assert_allclose(bridge["scalar_curvature"], 0, atol=2e-15)


def test_bridge_metric_from_parametric_derivatives():
    u = np.linspace(-.4, 1.5, 17)
    h = 1e-5
    centre = schwarzschild_bridge(u, .4)
    plus, minus = schwarzschild_bridge(u+h, .4), schwarzschild_bridge(u-h, .4)
    dRdl = (plus["R"]-minus["R"])/(plus["l"]-minus["l"])
    np.testing.assert_allclose(dRdl, centre["Rprime"], atol=2e-10)


def test_independent_finite_difference_constraints_agree():
    checks = connection_slice_record()["independent_checks"]
    for check in checks["finite_difference_constraints"]:
        assert check["hamiltonian_max_ell2"] < 2e-9
        assert check["momentum_max_ell2"] < 2e-9
        assert check["Rsecond_max_ell_error"] < 2e-9
    assert checks["surface_energy_projection_ell_error"] < 1e-14
    assert checks["surface_momentum_projection_ell_error"] < 1e-14
    assert checks["shell_hamiltonian_delta_ell_error"] < 1e-13
    assert checks["shell_momentum_delta_ell_error"] < 1e-13
    assert checks["transition_32_64_distance_max_error_over_ell"] < 2e-13


def test_unit_rescaling_preserves_slice_and_shell_coefficients():
    reference = wall_slice_data()
    fields0 = vacuum_transition(np.linspace(1, 1.5, 11), reference)
    for scale in (.01, 7, 100):
        wall = wall_slice_data(scale)
        fields = vacuum_transition(scale*np.linspace(1, 1.5, 11), wall)
        for key in ("gamma", "v", "V_plus", "Rprime_plus"):
            assert wall[key] == pytest.approx(reference[key], abs=3e-14)
        for key in ("energy_surface", "momentum_surface"):
            assert scale*wall[key] == pytest.approx(reference[key], abs=3e-14)
        for key in ("Rprime", "V"):
            np.testing.assert_allclose(fields[key], fields0[key], atol=3e-14)
        for key in ("Kll", "Ktheta", "Rsecond"):
            np.testing.assert_allclose(scale*fields[key], fields0[key], atol=1e-12)


def test_geometry_record_keeps_one_connected_slice_and_explicit_limits():
    record = connection_slice_record()
    cap, transition, bridge = record["cap"], record["transition"]["samples"], record["bridge"]["samples"]
    assert cap["R_over_ell"][0] == 0
    assert cap["regular_centre"] is True
    assert cap["contains_antitrapped_spheres"] is True
    assert cap["l_over_ell"][-1] == transition["l_over_ell"][0]
    assert transition["l_over_ell"][-1] == pytest.approx(bridge["l_over_ell"][0])
    assert transition["R_over_ell"][-1] == pytest.approx(bridge["R_over_ell"][0])
    assert record["parent_end"]["finite_display_is_physical_cutoff"] is False
    assert record["matching"]["shared_two_body_parent_established"] is False
    assert record["matching"]["causal_channel_established"] is False
    assert record["admission"]["priority_claim"] is False


@pytest.mark.parametrize("scale", [0, -1, float("nan"), float("inf")])
def test_invalid_units_are_rejected(scale):
    with pytest.raises(ValueError):
        wall_slice_data(scale)


def test_domain_exit_and_invalid_distance_order_are_rejected():
    wall = wall_slice_data()
    with pytest.raises(ValueError):
        vacuum_transition(.99, wall)
    with pytest.raises(ValueError):
        vacuum_transition(1.6, wall)
    with pytest.raises(ValueError):
        vacuum_transition(.8, dict(wall, taper_stop=.8))
    with pytest.raises(ValueError):
        transition_proper_distance([1, 1.5], wall)
    with pytest.raises(ValueError):
        transition_proper_distance([1.5, 1], wall, order=1000)
