"""Independent conservation, orientation and domain checks for the finite adapter."""

from dataclasses import FrozenInstanceError
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_moving_boundaries import (
    VacuumWallParameters, moving_boundary_record, perfect_fluid_flux,
    proper_time_samples, vacuum_wall_state,
)


def test_radiation_flux_is_direct_lorentz_contraction():
    for sign in (-1, 1):
        for rho in (0.1, 1.0, 8.0):
            stress = np.diag([rho, rho/3, rho/3, rho/3])
            for speed in (-0.8, -0.2, 0, 0.2, 0.8):
                gamma = 1/math.sqrt(1-speed*speed)
                U = gamma*np.array([1, speed, 0, 0])
                n = sign*gamma*np.array([speed, 1, 0, 0])
                assert perfect_fluid_flux(rho, rho/3, speed, sign) == pytest.approx(U@stress@n)


def test_vacuum_and_comoving_flux_vanish_but_moving_radiation_does_not():
    speeds = np.linspace(-0.9, 0.9, 19)
    assert np.all(perfect_fluid_flux(3, -3, speeds) == 0)
    assert perfect_fluid_flux(3, 1, 0) == 0
    assert perfect_fluid_flux(3, 1, 0.4) > 0
    assert perfect_fluid_flux(3, 1, -0.4) < 0


@pytest.mark.parametrize("speed", [-1, 1, 1.1, float("nan")])
def test_flux_rejects_null_spacelike_and_unrecorded_speeds(speed):
    with pytest.raises(ValueError):
        perfect_fluid_flux(1, 1/3, speed)


def test_negative_density_and_nonunit_normal_are_rejected():
    with pytest.raises(ValueError):
        perfect_fluid_flux(-1, -1/3, 0.2)
    with pytest.raises(ValueError):
        perfect_fluid_flux(1, 1/3, 0.2, 2)


def test_signed_junction_cannot_be_replaced_by_unsigned_roots():
    radii = np.linspace(1.2, 1.8, 11)
    state = vacuum_wall_state(radii)
    assert np.all(state["beta_minus"] > 0)
    assert np.all(state["beta_plus"] < 0)
    assert np.max(np.abs(state["angular_residual"])) < 2e-14
    unsigned = abs(state["beta_minus"])-abs(state["beta_plus"])-2*radii
    assert np.min(abs(unsigned)) > 1


def test_pressure_and_angular_stress_match_the_prescribed_law():
    state = vacuum_wall_state(np.linspace(1.2, 1.8, 19))
    sigma = 1/(2*math.pi)
    np.testing.assert_allclose(state["sigma_required"], sigma, atol=3e-16)
    np.testing.assert_allclose(state["pressure_required"], -sigma, atol=3e-16)


def test_independent_metric_acceleration_and_quadrature_checks():
    checks = moving_boundary_record()["independent_checks"]
    assert checks["lorentz_contraction_max_ell2_error"] < 1e-14
    assert checks["vacuum_flux_max_ell2"] == 0
    assert checks["metric_pressure_max_ell_error"] < 2e-10
    assert checks["metric_Ktau_max_ell_error"] < 2e-9
    assert checks["tangent_normalization_max_error"] < 2e-14
    assert checks["expanded_potential_max_error"] < 2e-14
    assert checks["proper_time_32_64_max_error_over_ell"] < 2e-14


def test_theory_is_immutable_and_time_reversal_does_not_retune_it():
    parameters = VacuumWallParameters()
    with pytest.raises(FrozenInstanceError):
        parameters.mass = 0.8
    expanding = vacuum_wall_state(1.4, parameters, expanding=True)
    contracting = vacuum_wall_state(1.4, parameters, expanding=False)
    assert expanding["Rdot"] == -contracting["Rdot"]
    for key in ("beta_minus", "beta_plus", "Rddot", "pressure_required", "sigma_required"):
        assert expanding[key] == contracting[key]


def test_unit_rescaling_preserves_dimensionless_motion_and_source():
    radii = np.linspace(1.2, 1.8, 9)
    reference = vacuum_wall_state(radii)
    reference_tau = proper_time_samples(radii)
    for scale in (0.01, 7.0, 100):
        parameters = VacuumWallParameters(H=1/scale, k=2/scale, mass=.4*scale)
        state = vacuum_wall_state(scale*radii, parameters)
        for key in ("Rdot", "beta_minus", "beta_plus", "f_minus", "f_plus"):
            np.testing.assert_allclose(state[key], reference[key], rtol=3e-14, atol=3e-14)
        for key in ("Rddot", "sigma_required", "pressure_required"):
            np.testing.assert_allclose(scale*state[key], reference[key], rtol=3e-14, atol=3e-14)
        np.testing.assert_allclose(proper_time_samples(scale*radii, parameters)/scale,
                                   reference_tau, rtol=3e-14, atol=3e-14)


def test_nonreal_motion_and_turning_points_are_explicit_exits():
    with pytest.raises(ValueError, match="turning"):
        vacuum_wall_state(0.5, VacuumWallParameters(H=1, k=1, mass=.01))
    turning_mass = 0.25*math.sqrt(.75)
    with pytest.raises(ValueError, match="turning"):
        vacuum_wall_state(.5, VacuumWallParameters(H=1, k=1, mass=turning_mass))


def test_horizons_are_adapter_exits_not_claimed_singularities():
    with pytest.raises(ValueError, match="horizon"):
        vacuum_wall_state(1.0)
    with pytest.raises(ValueError, match="horizon"):
        vacuum_wall_state(0.8)
    record = moving_boundary_record()
    assert record["vacuum_wall_control"]["domain"]["exclusions_are_physical_singularities"] is False


def test_zero_signed_normal_requires_another_adapter():
    radius = (0.8/.75)**(1/3)
    with pytest.raises(ValueError, match="normal"):
        vacuum_wall_state(radius, VacuumWallParameters(H=1, k=.5, mass=.4))


def test_quadrature_rejects_crossing_a_patch_and_invalid_grids():
    with pytest.raises(ValueError, match="crossing"):
        proper_time_samples(np.array([.9, 1.2]))
    for radii in ([1.5, 1.2], [1.2, 1.2], [[1.2, 1.5]], [1.2]):
        with pytest.raises(ValueError):
            proper_time_samples(radii)
    with pytest.raises(ValueError, match="order"):
        proper_time_samples([1.2, 1.8], order=1000)


def test_quadrature_cannot_jump_over_a_narrow_forbidden_interval():
    # At H=k=1 the critical mass is 2/(3 sqrt(3)); just below it q has
    # a very narrow negative minimum near R=sqrt(2/3).
    parameters = VacuumWallParameters(H=1, k=1, mass=2/(3*math.sqrt(3))*(1-1e-7))
    radii = np.array([.80, .84])
    nodes, _ = np.polynomial.legendre.leggauss(4)
    vacuum_wall_state(radii, parameters)
    vacuum_wall_state(.82+.02*nodes, parameters)
    with pytest.raises(ValueError, match="forbidden interval"):
        proper_time_samples(radii, parameters, order=4)


def test_clock_is_monotone_proper_duration_with_explicit_nonclaims():
    record = moving_boundary_record()
    control = record["vacuum_wall_control"]
    samples = control["samples"]
    assert samples["tau_over_ell"][0] == 0
    assert np.all(np.diff(samples["tau_over_ell"]) > 0)
    assert 0.2 < control["domain"]["proper_duration_over_ell"] < 0.6
    assert control["clock"]["cross_universe_age_established"] is False
    for key in ("embedded_in_shared_parent", "child_formation_established",
                "stability_established", "sibling_signalling_established"):
        assert control[key] is False
    assert record["admission"]["priority_claim"] is False
