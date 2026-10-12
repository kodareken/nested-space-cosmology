"""Seam-regular compression onto the same frozen window as the v1 representative."""
import hashlib
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

from pathlib import Path

import numpy as np

from recursive_horizons import nsc_finite_window_embedding as v1
from recursive_horizons.nsc_finite_window_embedding_smooth import (
    COLUMN_PHASE,
    FLAT_POWER,
    MINUS_RIGHT_PHASE,
    PLUS_RIGHT_PHASE,
    POLY_TERMS,
    RECORD_PATH,
    TOLERANCES,
    _LAB_ROOT,
    classical_derivative_agreement,
    column_profile,
    continuum_compression,
    discrete_compression,
    half_density_columns,
    seam_report,
    sine_reference_residual,
)
from recursive_horizons.nsc_spherical_coupling import PERIOD, apply_dirac


def test_v1_sources_are_unchanged():
    payload = json.loads(v1.RECORD_PATH.read_text())
    assert payload["schema"] == "NSC-FINITE-WINDOW-EMBEDDING-v1"
    assert payload["status"] == "CONSTRUCTED_COMPRESSION_NOT_INVARIANT"
    for relative, digest in payload["source_hashes"].items():
        path = Path(v1._LAB_ROOT) / relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


def test_constraint_integral_still_sees_the_v1_sine_lobes():
    reference = sine_reference_residual(128)
    assert reference["residual_norm"] < 1e-12


def test_flat_lobes_match_the_window_under_quadrature_refinement():
    coarse = continuum_compression(128)
    fine = continuum_compression(256)
    gap = np.linalg.norm(coarse["matrix"] - fine["matrix"])
    assert fine["matrix_error"] <= TOLERANCES["continuum_matrix_error"]
    assert fine["gram_defect"] <= TOLERANCES["continuum_gram_defect"]
    assert fine["hermiticity_defect"] <= TOLERANCES["continuum_hermiticity_defect"]
    assert fine["corner_frobenius"] <= TOLERANCES["continuum_block_error"]
    assert max(fine["onsite_frobenius"]) <= TOLERANCES["continuum_block_error"]
    assert max(fine["link_frobenius"]) <= TOLERANCES["continuum_block_error"]
    assert gap <= TOLERANCES["quadrature_matrix_gap"]
    assert coarse["matrix_error"] <= TOLERANCES["continuum_matrix_error"]
    for row in fine["support"]:
        assert row["interval"][1] - row["interval"][0] == 2
        assert row["mass_outside_interval"] == 0.0
        assert abs(row["mass"] - 1.0) < 1e-12
    assert min(fine["leak_fraction"]) > 0.9
    current = fine["sign_control_current"]
    assert abs(current["real"] - 0.25) <= TOLERANCES["continuum_matrix_error"]
    assert abs(current["imag"]) <= TOLERANCES["continuum_matrix_error"]
    assert "not the current of the full Dirac evolution" in current["meaning"]
    curvature = fine["subspace_population_second_derivative"]
    assert curvature["minus_two_leak_norm_squared"] < -10.0
    assert "not a closed evolution" in curvature["meaning"]


def test_seams_are_flat_through_third_derivative():
    assert FLAT_POWER >= 4
    assert POLY_TERMS == 6
    assert MINUS_RIGHT_PHASE == 0.0
    assert abs(np.sin(PLUS_RIGHT_PHASE)) > 0.5
    assert np.isfinite(COLUMN_PHASE)
    report = seam_report()
    assert report["cauchy_jet_max_through_order_3"] <= TOLERANCES["seam_cauchy_jet"]
    assert report["max_exact_seam_value"] <= TOLERANCES["seam_value"]
    assert report["max_exact_seam_derivative"] <= TOLERANCES["seam_derivative"]
    assert report["max_derivative_at_1e-3"] <= TOLERANCES["near_seam_derivative"]
    assert report["min_derivative_shrink"] >= TOLERANCES["derivative_shrink_per_decade"]
    agreement = classical_derivative_agreement()
    assert agreement["max_abs"] <= TOLERANCES["classical_derivative_agreement"]
    outside, outside_derivative = column_profile(np.array([5.0, 6.5, 7.5]), 0, "plus")
    assert np.max(np.abs(outside)) == 0.0
    assert np.max(np.abs(outside_derivative)) == 0.0


def test_half_density_columns_use_the_sigma2_frame():
    dx = PERIOD / 64
    xi = np.arange(64, dtype=float) * dx
    phi0, phi1 = half_density_columns(xi, dx)
    assert phi0.shape == (64, 6)
    assert np.allclose(phi1[:, 0], 1j * phi0[:, 0])
    assert np.allclose(phi1[:, 1], -1j * phi0[:, 1])
    assert apply_dirac.__module__ == "recursive_horizons.nsc_spherical_coupling"


def test_existing_dirac_compression_converges_at_n256_and_n512():
    coarse = discrete_compression(256)
    fine = discrete_compression(512)
    assert coarse["sampled_mass_outside_packet_arc"] == 0.0
    assert fine["sampled_mass_outside_packet_arc"] == 0.0
    assert coarse["matrix_error"] <= TOLERANCES["discrete_matrix_error_n256"]
    assert fine["matrix_error"] <= TOLERANCES["discrete_matrix_error_n512"]
    assert coarse["gram_defect"] <= TOLERANCES["discrete_gram_defect_n256"]
    assert fine["gram_defect"] <= TOLERANCES["discrete_gram_defect_n512"]
    assert coarse["matrix_error"] / fine["matrix_error"] >= TOLERANCES["discrete_error_ratio_256_over_512"]
    assert fine["bridge_image_max"] < coarse["bridge_image_max"]
    assert min(fine["leak_fraction"]) > 0.9
    assert fine["operator"] == "nsc_spherical_coupling.apply_dirac"


def test_saved_record_matches_these_sources():
    payload = json.loads(RECORD_PATH.read_text())
    assert payload["schema"] == "NSC-FINITE-WINDOW-EMBEDDING-SMOOTH-v2"
    assert payload["status"] == "SEAM_REGULAR_COMPRESSION_NOT_INVARIANT"
    assert payload["target_B_substituted"] is False
    assert payload["field_equations_changed"] is False
    assert payload["hand_inserted_matrix_couplings"] is False
    assert payload["all_nsc_claimed"] is False
    assert payload["interpretation"]["full_dirac_evolution_is_this_compression"] is False
    assert payload["interpretation"]["invariant_subspace"] is False
    assert payload["interpretation"]["complement_leak_required_small"] is False
    assert payload["interpretation"]["projected_initial_generator_is_window"] is True
    assert payload["quadrature_refinement"][-1]["matrix_error"] <= TOLERANCES["continuum_matrix_error"]
    assert payload["representative"]["poly_terms"] == 6
    assert payload["representative"]["flat_power"] == 4
    for relative, digest in payload["source_hashes"].items():
        path = Path(_LAB_ROOT) / relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
