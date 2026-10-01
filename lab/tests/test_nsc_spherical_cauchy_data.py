"""Initial data for a nonzero shift current. No evolution."""
import json
import os

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np

from recursive_horizons.nsc_spherical_cauchy_data import (
    RECORD_PATH,
    build_record,
    shift_momentum,
    shift_residual,
)
from recursive_horizons.nsc_spherical_coupling import TOL_INITIAL_CONSTRAINT
from recursive_horizons.nsc_spherical_galerkin_coupling import build_grid


def test_closed_form_antiderivative_keeps_a_nonzero_mean():
    grid = build_grid(32)
    wavenumber = 2.0 * np.pi / grid.length
    current = 0.3 + np.sin(wavenumber * grid.xi_q)
    snapshot = current.copy()
    momentum, info = shift_momentum(grid, current)
    residual = shift_residual(grid, momentum, current)
    assert np.max(np.abs(current - snapshot)) == 0.0
    assert info["source_mean_subtracted"] is False
    assert info["antiderivative_accepted"] is True
    assert abs(info["p_Q_mean"]) <= 1e-12
    assert abs(float(np.mean(residual)) - 0.3) < 1e-10
    assert np.max(np.abs(residual - np.mean(residual))) < 1e-9
    q0 = info["Q"]
    closed = -np.cos(wavenumber * grid.xi_q) / (wavenumber * q0)
    # The constant in the current does not change the mean-zero antiderivative.
    pure_sine, _sine_info = shift_momentum(grid, np.sin(wavenumber * grid.xi_q))
    assert np.max(np.abs(momentum - pure_sine)) < 1e-10
    assert np.max(np.abs(grid.A_g @ momentum - closed)) < 1e-8


def test_constant_current_has_no_periodic_momentum():
    grid = build_grid(32)
    current = np.ones(grid.nq)
    momentum, info = shift_momentum(grid, current)
    residual = shift_residual(grid, momentum, current)
    assert info["source_mean_subtracted"] is False
    assert abs(info["current_mean"] - 1.0) < 1e-15
    assert info["p_Q_max"] < 1e-10
    assert np.max(np.abs(residual - 1.0)) < 1e-10


def test_record_controls_rank6_packet_and_embedding_source():
    record = build_record(64)
    assert RECORD_PATH.is_file()
    saved = json.loads(RECORD_PATH.read_text())
    assert saved["schema"] == record["schema"]
    assert saved["action_changed"] is False
    assert saved["coefficients_changed"] is False
    controls = record["closed_form_controls"]
    assert controls["mixed_current_unchanged"] is True
    assert controls["source_mean_subtracted"] is False
    assert controls["closed_form_p_Q_gap"] < 1e-9
    assert abs(controls["residual_mean"] - 0.3) < 1e-10
    assert controls["residual_minus_mean_max"] < 1e-9
    assert controls["constant_p_Q_shift_gap"] < 1e-9
    assert controls["constant_p_Q_hamilton_gap_at_zero_p_chi"] < 1e-9
    assert controls["constant_p_Q_hamilton_gap_at_nonzero_p_chi"] > 1e-4
    assert controls["constant_p_Q_chi_rate_gap"] < 1e-8
    assert controls["pure_mean_p_Q_max"] < 1e-10
    assert controls["pure_mean_residual_gap"] < 1e-10

    carrier = record["positive_carrier"]
    assert carrier["caller_state_unchanged"] is True
    assert carrier["columns_unchanged"] is True
    assert carrier["source_mean_subtracted"] is False
    assert carrier["force_beta_plus_multiplicity_Pmom_max"] < 1e-12
    assert carrier["occupations_equal_within_regional_pairs"] is True
    assert carrier["current_minus_expected"] < 1e-10
    assert carrier["current_mean"] < -0.5
    assert carrier["mean_compatible"] is False
    assert carrier["p_Q_max"] < 1e-8
    assert carrier["projected_shift_beyond_mean_max"] <= TOL_INITIAL_CONSTRAINT
    assert abs(carrier["projected_momentum_max"] - abs(carrier["current_mean"])) < 1e-8
    assert carrier["solves_full_shift"] is False
    assert carrier["solves_projected_hamilton"] is True
    assert carrier["positive_r"] is True and carrier["positive_Q"] is True
    assert carrier["chi_max"] == 0.0
    assert carrier["p_r_max"] == 0.0 and carrier["p_chi_max"] == 0.0
    assert carrier["chart_proper_max"] <= TOL_INITIAL_CONSTRAINT
    assert carrier["chi_rate_gap"] <= TOL_INITIAL_CONSTRAINT
    assert any("current mean" in item for item in carrier["restrictions"])
    assert carrier["radius_converged"] is True
    assert carrier["rho_independent_of_r_max"] < 1e-8

    rank6 = record["rank6_packet"]
    assert rank6["preparation_problems"] == []
    assert rank6["mean_compatible"] is True
    assert rank6["current_max"] <= TOL_INITIAL_CONSTRAINT
    assert rank6["p_Q_max"] <= TOL_INITIAL_CONSTRAINT
    assert rank6["solves_projected_hamilton"] is True
    assert rank6["solves_projected_shift_except_mean"] is True
    assert rank6["solves_full_shift"] is True
    assert rank6["chart_proper_max"] <= TOL_INITIAL_CONSTRAINT
    assert rank6["positive_r"] is True and rank6["r_min_quadrature"] > 1.0
    assert rank6["source_mean_subtracted"] is False
    witness = record["saved_v5_rank6_witness"]
    assert witness["nq1024_current_max"] < 1e-10
    assert witness["nq1024_chart_proper_velocity_max"] == 0.0

    classical = record["classical_embedding_integral"]
    assert classical["pairwise_occupations"]["equal_within_regional_pairs"] is True
    assert abs(classical["momentum_integral"]) < 1e-12
    assert abs(classical["mean_of_minus_multiplicity_density"]) < 1e-12

    embedding = record["finite_window_embedding"]
    assert embedding["sampling"]["band_retained_phi0"] > 0.99
    assert embedding["columns_unchanged"] is True
    assert embedding["occupations_equal_within_regional_pairs"] is True
    assert embedding["source_mean_subtracted"] is False
    assert embedding["force_beta_plus_multiplicity_Pmom_max"] < 1e-12
    assert embedding["mean_compatible"] is False
    assert embedding["current_max"] > 1.0
    assert 1e-4 < abs(embedding["current_mean"]) < 1e-2
    assert embedding["p_Q_max"] > 1.0
    assert embedding["projected_shift_beyond_mean_max"] <= TOL_INITIAL_CONSTRAINT
    assert abs(embedding["projected_momentum_max"] - abs(embedding["current_mean"])) < 1e-8
    assert embedding["full_momentum_max"] > 0.1
    assert embedding["solves_full_shift"] is False
    assert embedding["solves_projected_shift_except_mean"] is True
    assert embedding["solves_projected_hamilton"] is True
    assert embedding["hamilton_gap_against_zero_p_Q"] <= TOL_INITIAL_CONSTRAINT
    assert embedding["positive_r"] is True and embedding["positive_Q"] is True
    assert embedding["r_min_quadrature"] > 1.0
    assert embedding["chi_max"] == 0.0
    assert embedding["chart_proper_max"] <= TOL_INITIAL_CONSTRAINT
    assert embedding["chi_rate_gap"] <= TOL_INITIAL_CONSTRAINT
    assert embedding["chi_rate_fine_max"] > 1.0
    assert embedding["rho_independent_of_r_max"] < 1e-8
    assert any("current mean" in item for item in embedding["restrictions"])
    assert any("geometry subspace" in item for item in embedding["restrictions"])

    finer = record["finite_window_embedding_nq_8nf_momentum"]
    assert finer["radius_solved"] is False
    assert finer["source_mean_subtracted"] is False
    assert finer["antiderivative_accepted"] is True
    assert 1e-4 < abs(finer["current_mean"]) < 1e-2
    assert finer["projected_shift_beyond_mean_max"] <= TOL_INITIAL_CONSTRAINT
    assert finer["full_momentum_max"] > 0.1
    assert finer["p_Q_max"] > 1.0
