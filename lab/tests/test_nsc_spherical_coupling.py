"""Provisional finite-source coupling. These checks are not a regeneration claim."""
import inspect

import numpy as np

from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons.nsc_covariant_operator import CovariantStaticMetric
from recursive_horizons.nsc_spherical_feedback_action import velocity_rhs


def test_locked_chart_counts_gaussian_once_and_not_the_vacuum_branch():
    source = inspect.getsource(coupling)
    assert "VacuumMatchedDiracAction" not in source
    assert "0.0012511140479781876" not in source
    coefficients = coupling.locked_coefficients()
    assert coefficients["V_rel"] == 0
    assert coupling.magnetic_radius_square(coefficients) == np.float64(1) or abs(
        coupling.magnetic_radius_square(coefficients) - 1
    ) < 1e-12
    assert coupling.SEA_SUBTRACTED is False
    assert coupling.GAMMA_REST_INCLUDED is False
    assert coupling.VACUUM_MATCHED_PREREQUISITE is False
    assert coupling.GLOBAL_REGENERATION is False
    assert coupling.ABSOLUTE_VACUUM_CLAIM is False


def test_periodic_derivative_is_real_antisymmetric_and_ap_momentum_is_existing():
    derivative = coupling.periodic_derivative(32, 8.0)
    assert np.allclose(derivative, np.real(derivative))
    assert np.allclose(derivative.T, -derivative, atol=1e-12)
    coordinate = np.arange(32) * (8 / 32)
    wave = np.sin(2 * np.pi * coordinate / 8)
    assert np.allclose(derivative @ wave, (2 * np.pi / 8) * np.cos(2 * np.pi * coordinate / 8), atol=1e-12)
    assert np.allclose(derivative @ np.ones(32), 0, atol=1e-12)
    momentum, metric = coupling.antiperiodic_momentum(32, 8.0)
    assert isinstance(metric, CovariantStaticMetric)
    assert np.allclose(momentum, metric.momentum_matrix)
    assert np.allclose(momentum, momentum.conj().T, atol=1e-12)
    indices = np.arange(32) - 16
    assert np.allclose(metric.momenta, 2 * np.pi * (indices + 0.5) / 8)
    assert not np.allclose(metric.momenta, 2 * np.pi * indices / 8)


def test_column_source_matches_dense_api_and_phi_dot_has_no_multiplicity():
    system, state = coupling.build_system(32)
    difference = coupling.dense_force_difference(system, state)
    assert difference["multiplicity"] == 4
    assert difference["vacuum_subtracted"] is False
    assert max(difference["L"], difference["Q"], difference["beta"], difference["energy"]) < 1e-9
    assert coupling.dense_hamiltonian_difference(system, state) < 1e-9
    rate = coupling.rates(system, state)
    image0, image1 = coupling.apply_dirac(
        state.phi0, state.phi1, system.length_density, state.Q, system.shift,
        system.kappa, system.momentum,
    )
    assert np.allclose(rate.phi0, -1j * image0)
    assert np.allclose(rate.phi1, -1j * image1)
    assert np.max(np.abs(rate.phi0 + system.multiplicity * 1j * image0)) > 1
    values = np.linalg.eigvalsh(coupling.covariance_matrix(system, state))
    assert np.allclose(np.sort(values)[-6:], np.sort(system.occupations))
    assert np.max(np.abs(values[:-6])) < 1e-10
    assert system.preparation["complement_occupation"] == 0
    assert system.preparation["vacuum_identification"] is False
    assert np.max(np.abs(rate.force_beta)) < 1e-10
    varied = state.copy()
    varied.r = state.r * 1.8
    other = coupling.source_from_columns(system, varied)
    base = coupling.source_from_columns(system, state)
    assert np.max(np.abs(other["force_L"] - base["force_L"])) < 1e-10


def test_source_free_radius_solves_the_discrete_constraint_without_a_product_rule():
    system, _state = coupling.build_system(32)
    residual = coupling._radius_residual(system, np.ones(system.points), np.zeros(system.points))
    assert np.max(np.abs(residual)) < 1e-9
    direction = np.sin(2 * np.pi * system.xi / system.length)
    base = np.ones(system.points)
    rho = np.zeros(system.points)
    eps = 1e-6
    numeric = (
        coupling._radius_residual(system, base + eps * direction, rho)
        - coupling._radius_residual(system, base - eps * direction, rho)
    ) / (2 * eps)
    analytic = coupling._radius_jacobian(system, base) @ direction
    assert np.max(np.abs(numeric - analytic)) < 1e-6


def test_initial_state_has_zero_normal_velocity_and_a_solved_radius():
    system, state = coupling.build_system(64)
    solved, info = coupling.solve_initial_radius(system, state)
    assert info["converged"] is True
    assert info["imposed_radius"] is False
    assert info["bracket_positive"] is True
    assert np.std(solved.r) > 1e-2
    residuals = coupling.constraint_residuals(system, solved)
    assert residuals["hamilton_max"] < coupling.TOL_INITIAL_CONSTRAINT
    assert residuals["momentum_max"] < coupling.TOL_INITIAL_CONSTRAINT
    normal, coordinate, shift = coupling.normal_radius_velocity(system, solved)
    assert np.max(np.abs(normal)) < 1e-12
    assert np.max(np.abs(coordinate)) > 1e-3
    assert np.max(np.abs(shift)) > 1e-3
    reduced = coupling.reduced_dirac(system, solved)
    assert reduced["frozen_embedded_exactly"] is False
    assert reduced["old_circulation_transferred"] is False
    assert reduced["max_link_difference"] > 1e-3
    gap = coupling.kinematic_module_gap(system, solved)
    assert gap["r"] < 1e-12
    assert gap["chi"] < 1e-12
    checked = velocity_rhs(
        float(system.length_density[0]), float(solved.Q[0]), 0.0, float(solved.r[0]),
        float((system.derivative @ solved.r)[0]), 0.0, float(system.shift[0]),
        float((system.derivative @ system.shift)[0]), 0.0, 0.0, 0.0, system.A, system.C_W,
    )
    assert np.isfinite(float(checked[1]))


def test_discrete_hamiltonian_fd_signs_and_mutation():
    system, state = coupling.build_system(32)
    solved, info = coupling.solve_initial_radius(system, state)
    assert info["converged"]
    rng = np.random.default_rng(11)
    trial = solved.copy()
    trial.p_Q = 0.05 * rng.normal(size=system.points)
    trial.p_r = 0.05 * rng.normal(size=system.points)
    trial.p_chi = 0.05 * rng.normal(size=system.points)
    trial.chi = 0.02 * np.sin(2 * np.pi * system.xi / system.length)
    errors = coupling.hamiltonian_fd_errors(system, trial)
    assert errors["relative_worst"] < coupling.TOL_FD_RELATIVE
    mutation = coupling.mutation_power(system, solved)
    assert mutation["breaks_balance"] is True
    assert mutation["fieldwork_identity_error"] < 1e-6
    assert abs(mutation["intact_power"]) < 1e-6


def test_bounded_control_reports_constraints_refinement_and_does_not_invent_success(tmp_path):
    record = coupling.run_bounded_control(tmp_path / "control.json")
    assert record["schema"] == coupling.SCHEMA
    assert record["provisional"] is True
    assert record["global_regeneration"] is False
    assert record["absolute_vacuum_claim"] is False
    assert record["vacuum_matched_prerequisite"] is False
    assert record["phi_dot"] == "-i H Phi"
    assert record["dense_column_crosscheck_N32"]["L"] < 1e-9
    assert record["initial_hamilton_max"] < coupling.TOL_INITIAL_CONSTRAINT
    assert record["initial_normal_velocity_max"] < 1e-12
    assert record["initial_coordinate_velocity_max"] > 1e-3
    assert record["hamiltonian_fd"]["relative_worst"] < coupling.TOL_FD_RELATIVE
    assert record["mutation"]["breaks_balance"] is True
    assert record["reduced_operator"]["max_link_difference"] > 1e-3
    assert record["initial_solve"]["imposed_radius"] is False
    assert record["verdict"] in {
        "PROVISIONAL_INITIAL_VALIDATION",
        "PROVISIONAL_CONSTRAINT_DRIFT",
        "STOP_POSITIVE_CHART",
        "FAIL_INITIAL_SOLVE",
        "FAIL_DISCRETE_RHS",
    }
    runs = record["runs"]
    assert "N64_T0.005" in runs
    assert "N64_T0.005_dt_half" in runs
    for name in ("N64_T0.005", "N64_T0.005_dt_half"):
        run = runs[name]
        assert run["success_interval_completed"] or record["verdict"] == "STOP_POSITIVE_CHART"
        if run["success_interval_completed"]:
            assert run["energy_drift"] < 1e-4 * max(1.0, abs(run["checkpoints"][0]["energy"]))
            assert run["unitarity_max"] < 1e-8
            assert run["number_drift"] < 1e-8
    if record["verdict"] == "PROVISIONAL_INITIAL_VALIDATION":
        assert runs["N64_T0.005"]["hamilton_max"] <= coupling.TOL_CONSTRAINT_VALIDATION
        assert runs["N64_T0.005"]["momentum_max"] <= coupling.TOL_CONSTRAINT_VALIDATION
        assert "N128_T0.005" in runs
        assert "N64_T0.05" in runs
        assert runs["N64_T0.05"]["success_interval_completed"] is True
    if record["verdict"] == "STOP_POSITIVE_CHART":
        assert any(run.get("stopped") for run in runs.values() if isinstance(run, dict))
    saved = (tmp_path / "control.json").read_text()
    assert "PROVISIONAL" in saved
    assert "global_regeneration" in saved
