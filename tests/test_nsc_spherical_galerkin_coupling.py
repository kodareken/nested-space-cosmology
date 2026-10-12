"""Fourier–Galerkin coupling. Algebra checks are not a renewal claim."""
import numpy as np
import pytest

from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


def _manufactured(fermions=16, quadrature=64, seed=3):
    grid = galerkin.build_grid(fermions, quadrature=quadrature)
    phi0, phi1 = galerkin.manufactured_columns(fermions, seed=seed)
    state = galerkin.blank_state(grid, phi0, phi1)
    return grid, state


def test_tolerances_stay_at_the_v1_lines():
    assert galerkin.TOL_CONSTRAINT_VALIDATION == 1e-3
    assert galerkin.TOL_CONSTRAINT_VALIDATION == coupling.TOL_CONSTRAINT_VALIDATION
    assert galerkin.TOL_INITIAL_CONSTRAINT == 1e-8
    assert galerkin.RENEWAL is False
    assert galerkin.POST_STEP_FILTER is False
    assert galerkin.CONSTRAINT_PROJECTION is False
    assert galerkin.OLD_CIRCULATION_TRANSFERRED is False


def test_even_geometry_is_rejected_and_odd_interpolation_is_an_isometry():
    with pytest.raises(ValueError):
        galerkin.periodic_interpolation(64, 256, galerkin.PERIOD)
    grid = galerkin.build_grid(32, quadrature=128)
    assert grid.ng == 31 and grid.ng % 2 == 1
    assert grid.nf == 32 and grid.nq == 128
    eye = np.eye(grid.ng)
    assert np.allclose(grid.weight * grid.A_g.T @ grid.A_g, eye, atol=1e-10)
    assert np.allclose(grid.U_f.conj().T @ grid.U_f, np.eye(grid.nf), atol=1e-10)
    generator = np.random.default_rng(1)
    geometry = generator.normal(size=grid.ng)
    probe = generator.normal(size=grid.nq)
    assert np.allclose(
        grid.weight * np.dot(grid.A_g @ geometry, probe),
        np.dot(geometry, grid.weight * grid.A_g.T @ probe),
    )
    column = generator.normal(size=grid.nf) + 1j * generator.normal(size=grid.nf)
    image = generator.normal(size=grid.nq) + 1j * generator.normal(size=grid.nq)
    assert abs(np.vdot(grid.U_f @ column, image) - np.vdot(column, grid.U_f.conj().T @ image)) < 1e-9


def test_derivatives_intertwine_and_the_ap_band_keeps_its_symbol():
    grid = galerkin.build_grid(32, quadrature=128)
    mode = 7
    coarse = np.cos(2 * np.pi * mode * grid.xi_g / grid.length)
    prolonged = galerkin.prolong_geometry(grid, coarse)
    analytic = -(2 * np.pi * mode / grid.length) * np.sin(2 * np.pi * mode * grid.xi_q / grid.length)
    assert np.max(np.abs(grid.derivative @ prolonged - analytic)) < 1e-10
    wave = grid.modes_f[grid.nf // 3]
    sample = np.exp(1j * 2 * np.pi * wave * grid.xi_f / grid.length) / np.sqrt(grid.nf)
    lifted = grid.U_f @ sample
    momentum = 2 * np.pi * wave / grid.length
    assert np.max(np.abs(grid.momentum @ lifted - momentum * lifted)) < 1e-10
    assert abs(np.linalg.norm(lifted) - 1) < 1e-10


def test_hamiltonian_derivatives_and_lifted_work_balance():
    grid, state = _manufactured(32, 128)
    state.r = 2 + 0.05 * np.cos(2 * np.pi * grid.xi_g / grid.length)
    state.Q = state.Q + 0.15 * np.cos(2 * np.pi * (grid.ng // 2) * grid.xi_g / grid.length)
    generator = np.random.default_rng(2)
    state.p_Q = 0.01 * generator.normal(size=grid.ng)
    state.p_r = 0.01 * generator.normal(size=grid.ng)
    state.p_chi = 0.01 * generator.normal(size=grid.ng)
    state.chi = 0.01 * np.sin(4 * np.pi * grid.xi_g / grid.length)
    errors = galerkin.hamiltonian_fd_errors(grid, state)
    assert errors["relative_worst"] < galerkin.TOL_FD_RELATIVE
    balance = galerkin.work_balance(grid, state)
    assert balance["uses_lifted_Qdot"] is True
    assert balance["lifted_identity_error"] < 1e-8
    assert balance["unprojected_identity_error"] > 1


def test_columns_match_the_source_api_and_keep_occupation_eigenvalues():
    grid, state = _manufactured(16, 32)
    report = galerkin.source_column_report(grid, state)
    assert report["dense_forces"]["multiplicity"] == 4
    assert report["dense_forces"]["vacuum_subtracted"] is False
    assert report["hamiltonian_image_max"] < 1e-9
    assert max(report["dense_forces"]["L"], report["dense_forces"]["Q"], report["dense_forces"]["beta"]) < 1e-9
    assert report["occupation_max_abs_difference"] < 1e-10
    assert report["occupation_tail"] < 1e-10
    assert report["declared_occupation_gap"] < 1e-10
    rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    direct = grid.U_f.conj().T @ (-1j * bundle["source"]["image0"])
    assert np.allclose(rate.phi0, direct)
    assert np.max(np.abs(rate.phi0 + grid.fine.multiplicity * 1j * direct)) > 1


def test_mode31_alias_is_a_true_high_mode_and_low_band_brackets_agree():
    alias = galerkin.alias_diagnostic()
    assert alias["low_band_combined_max"] < 1e-12
    assert alias["uses_true_high_mode"] is True
    assert abs(alias["mode31_collocation_64"] / np.pi - 8) < 1e-9
    assert abs(alias["mode31_collocation_64_l2_over_2pi"] - 22.62741699796956) < 1e-9
    assert alias["mode31_collocation_128"] < 1e-10
    assert alias["mode31_galerkin_quadrature_256"] < 1e-9
    assert alias["mode31_equals_folded_partner_on_64"] is True
    assert alias["mode31_differs_from_folded_partner_on_256"] is True
    assert alias["folded_partner_mode"] == -33


def test_in_band_source_is_pointwise_and_a_top_mode_product_is_held_out():
    grid, state = _manufactured(32, 128)
    rho = 0.02 * np.cos(2 * np.pi * grid.xi_q / grid.length)
    solved, info = galerkin.solve_initial_radius(grid, state, rho=rho)
    assert info["converged"] is True
    assert info["imposed_radius"] is False
    assert info["filtered_v1_radius"] is False
    assert info["diagnostics"]["projected_hamilton_max"] < galerkin.TOL_INITIAL_CONSTRAINT
    assert info["diagnostics"]["full_hamilton_max"] < galerkin.TOL_INITIAL_CONSTRAINT
    assert np.min(galerkin.prolong_geometry(grid, solved.r)) > 0
    top = 5.0 * np.cos(2 * np.pi * (grid.ng // 2) * grid.xi_q / grid.length)
    _solved_top, top_info = galerkin.solve_initial_radius(grid, state, rho=top)
    assert top_info["converged"] is True
    assert top_info["diagnostics"]["positive_r"] is True
    assert top_info["diagnostics"]["projected_hamilton_max"] < galerkin.TOL_INITIAL_CONSTRAINT
    assert top_info["diagnostics"]["full_hamilton_max"] > galerkin.TOL_CONSTRAINT_VALIDATION
    assert top_info["diagnostics"]["held_out_hamilton_max"] > galerkin.TOL_CONSTRAINT_VALIDATION


def test_physical_initial_slice_does_not_treat_the_projected_residual_as_success():
    loaded = galerkin.load_physical_columns(64)
    assert loaded[3] == []
    _grid, solved, info = galerkin._initial_slice(64, 256, loaded)
    assert info["converged"] is True
    assert info["filtered_v1_radius"] is False
    diagnostics = info["diagnostics"]
    assert diagnostics["positive_r"] and diagnostics["positive_Q"]
    assert diagnostics["projected_hamilton_max"] < galerkin.TOL_INITIAL_CONSTRAINT
    assert diagnostics["full_hamilton_max"] > galerkin.TOL_CONSTRAINT_VALIDATION
    assert diagnostics["held_out_hamilton_max"] > galerkin.TOL_CONSTRAINT_VALIDATION
    assert diagnostics["hamilton_modes"]["unresolved_fraction"] > 0.5
    assert info["velocities"]["chart_proper_max"] == 0
    assert solved.p_Q.shape == (63,)
    assert np.max(np.abs(solved.p_r)) == 0
    links = galerkin.physical_link_report(galerkin.build_grid(64, quadrature=256), solved)
    assert links["B_phys_equals_old_B"] is False
    assert links["old_circulation_transferred"] is False
    assert links["copied_from_old_B"] is False
    assert links["max_difference_from_old_B"] > 1


def test_preparation_contract_follows_the_coupling_owner():
    phi0, phi1, preparation, problems = galerkin.load_physical_columns(64)
    assert problems == []
    assert preparation["phase_on_minus_spinor_column"] is True
    assert preparation["phase_applied_to_odd_lobe"] is False
    xi = np.arange(64) * (galerkin.PERIOD / 64)
    rejected = dict(preparation)
    rejected["phase_applied_to_odd_lobe"] = True
    assert galerkin.preparation_problems(rejected, phi0, phi1, xi)


def test_one_step_uses_the_subspace_rate_and_keeps_the_positive_chart():
    grid, state = _manufactured(16, 64)
    state.r = np.full(grid.ng, 3.5)
    step = galerkin.rk4_step(grid, state, 1e-4)
    prolonged = galerkin.prolong_state(grid, step)
    assert np.min(prolonged.r) > 0 and np.min(prolonged.Q) > 0
    assert step.r.shape == (grid.ng,)
    assert step.phi0.shape == (grid.nf, 6)
