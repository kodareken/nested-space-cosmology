"""Compression of the fixed Dirac operator onto the frozen finite window."""
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

from recursive_horizons.nsc_finite_window_embedding import (
    COLUMN_PHASE,
    DEGREE,
    MINUS_RIGHT_PHASE,
    PLUS_RIGHT_PHASE,
    RECORD_PATH,
    _LAB_ROOT,
    continuum_compression,
    discrete_compression,
    frozen_window,
    inherited_law_reuse,
    real_lobe_overlap_bound,
    restrictive_probe,
)


def test_restrictive_packet_matches_onsite_and_not_the_link():
    probe = restrictive_probe(128)
    assert probe["max_onsite"] < 1e-4
    assert probe["min_link"] > 5.0
    assert probe["conjugate"]["conjugate_defect"] < 1e-8
    assert probe["frozen_conjugate"]["conjugate_defect"] > 0.18
    assert min(probe["leak_fraction"]) > 0.9
    assert probe["same_real_envelope"] is True


def test_real_lobes_at_recorded_carrier_cannot_match_plus_link():
    bound = real_lobe_overlap_bound()
    assert bound["obstruction"] is True
    assert bound["gap"] > 0.01
    assert bound["carrier_k"] == 1.0


def test_sine_lobes_with_one_plus_phase_reproduce_the_window():
    built, symbolic = frozen_window()
    assert symbolic.shape == (6, 6)
    assert symbolic == symbolic.H
    report = continuum_compression(96)
    assert report["matrix_error"] < 1e-12
    assert report["gram_defect"] < 1e-12
    assert report["hermiticity_defect"] < 1e-12
    assert report["corner_frobenius"] < 1e-12
    assert max(report["onsite_frobenius"]) < 1e-12
    assert max(report["link_frobenius"]) < 1e-12
    assert DEGREE == 4
    assert MINUS_RIGHT_PHASE == 0.0
    assert abs(np.sin(PLUS_RIGHT_PHASE)) > 0.5
    assert report["conjugate"]["conjugate_defect"] == np.float64(
        abs(built[1, 2] - np.conjugate(built[0, 3]))
    ) or abs(
        report["conjugate"]["conjugate_defect"] - abs(built[1, 2] - np.conjugate(built[0, 3]))
    ) < 1e-9
    for row in report["support"]:
        assert row["interval"][1] - row["interval"][0] == 2
        assert row["mass_outside_interval"] == 0.0
        assert abs(row["mass"] - 1.0) < 1e-12
    assert min(report["leak_fraction"]) > 0.98
    current = report["sign_control_current"]
    assert abs(current["real"] - 0.25) < 1e-12
    assert abs(current["imag"]) < 1e-12
    assert report["subspace_population_second_derivative"]["minus_two_leak_norm_squared"] < -100.0


def test_discrete_dirac_samples_approach_the_continuum_compression():
    coarse = discrete_compression(128)
    fine = discrete_compression(256)
    assert coarse["gram_defect"] < 1e-12
    assert fine["gram_defect"] < 1e-12
    assert fine["matrix_error"] < coarse["matrix_error"] / 3.0
    assert fine["matrix_error"] < 0.08
    assert coarse["sampled_mass_outside_packet_arc"] == 0.0
    assert fine["sampled_mass_outside_packet_arc"] == 0.0


def test_inherited_law_is_the_existing_identity():
    reused = inherited_law_reuse()
    assert reused["recomputed"] is False
    assert reused["normalized_shift_exact"] is True
    assert reused["omega"] == "3/2"
    assert reused["depth"] == 3
    # The column phase is mode data. It is not a new matrix entry.
    assert np.isfinite(COLUMN_PHASE)


def test_saved_record_matches_these_sources():
    payload = json.loads(RECORD_PATH.read_text())
    assert payload["schema"] == "NSC-FINITE-WINDOW-EMBEDDING-v1"
    assert payload["status"] == "CONSTRUCTED_COMPRESSION_NOT_INVARIANT"
    assert payload["target_B_substituted"] is False
    assert payload["field_equations_changed"] is False
    assert payload["hand_inserted_matrix_couplings"] is False
    assert payload["all_nsc_claimed"] is False
    assert payload["interpretation"]["closed_stage1_evolution_retained"] is False
    assert payload["interpretation"]["instantaneous_projected_generator_is_J"] is True
    assert payload["construction"]["matrix_error"] < 1e-12
    for relative, digest in payload["source_hashes"].items():
        path = Path(_LAB_ROOT) / relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
