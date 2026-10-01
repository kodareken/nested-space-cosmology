"""Independent saved-domain controls; sealed v2 bytes are never rewritten.

The temporal solver below uses the Galerkin matrix owner but none of the
consumer's evolution, cross-control, covariance, or projection functions.
"""
from __future__ import annotations

from copy import deepcopy
import json
import time

import numpy as np
import pytest

import derive_nsc_coupled_local_response_audit as audit
from derive_nsc_coupled_local_response_audit import (
    OUTPUT, OUTPUT_NPZ, SEALED, covariance, load_inputs, max_abs, occupation,
    replay, saved_consistency_errors, sha256,
)


@pytest.fixture(scope="module")
def inputs():
    return load_inputs()


@pytest.fixture(scope="module")
def independent(inputs):
    before = {path: sha256(path) for path in SEALED}
    result = replay(inputs, deadline=time.process_time() + 300.0)
    assert before == {path: sha256(path) for path in SEALED}
    return result


def test_same_original_preparation_and_geometry(inputs):
    assert saved_consistency_errors(inputs) == []
    assert inputs["times"].size == 8
    assert inputs["times"][0] == pytest.approx(0.05)
    assert inputs["times"][-1] == pytest.approx(0.085)
    assert not np.array_equal(inputs["observer"], inputs["columns"][:, :2])
    projected = inputs["observer"].conj().T @ inputs["columns"]
    exterior = inputs["columns"] - inputs["observer"] @ projected
    cross = (projected * inputs["weights"]) @ exterior.conj().T
    assert float(np.linalg.norm(cross)) == pytest.approx(0.4422387105894591, abs=1e-12)


def test_probe_and_full_cross_are_distinct_saved_domains(inputs):
    saved = inputs["saved"]
    full = saved["covariance_retained"] - saved["covariance_without_cross"]
    assert full.shape == (8, 2, 2)
    assert saved["probe_cross_occupation_movement"].shape == (3, 2)
    assert saved["probe_times"][-1] == pytest.approx(0.06)
    assert max_abs(occupation(full)) == pytest.approx(0.2043417632796783, abs=1e-15)
    assert max_abs(saved["probe_cross_occupation_movement"]) == pytest.approx(0.08797142267734728, abs=1e-15)
    assert float(np.max(np.linalg.norm(full, axis=(1, 2)))) == pytest.approx(0.2351316538667172, abs=1e-15)


def test_independent_full_split_reproduces_the_reference(inputs, independent):
    gap = max_abs(independent["occupation_full"] - inputs["saved"]["occupation_full"])
    assert gap < 2e-12
    assert independent["dense_propagator_history_bytes"] == 0
    assert independent["column_width"] == 12
    assert independent["free_exterior_width"] == 6
    assert independent["exterior_projection_defect"] < 2e-12


def test_cross_effect_against_an_independent_full_split(inputs, independent):
    saved = inputs["saved"]
    streamed_cross = saved["covariance_retained"] - saved["covariance_without_cross"]
    error = max_abs(occupation(streamed_cross) - occupation(independent["covariance_cross"]))
    effect = max_abs(occupation(independent["covariance_cross"]))
    assert effect == pytest.approx(0.20437297122219308, abs=1e-11)
    assert 1e-4 < error < 1.3e-4
    assert error / effect < 0.001
    # The v2 algebraic superposition residual is not a reduction-error estimate.
    roundoff = inputs["record"]["phases"]["full_controls"]["omissions"]["initial_cross"]["numerical_error"]
    assert error > 1e8 * roundoff


def test_drive_omission_uses_the_same_projected_preparation(inputs, independent):
    saved = inputs["saved"]
    error = max_abs(saved["occupation_outside_drive_off"] - independent["occupation_drive_off"])
    effect = max_abs(independent["occupation_drive_off"] - independent["occupation_full"])
    assert effect == pytest.approx(0.1549746431210934, abs=1e-11)
    assert 2.8e-4 < error < 2.9e-4
    assert error / effect < 0.002
    # Dropping the actual cross keeps both diagonal blocks, unlike dropping the drive.
    exterior_covariance = occupation(independent["covariance_without_cross"]) - independent["occupation_drive_off"]
    assert float(np.min(exterior_covariance)) >= -1e-12
    assert max_abs(exterior_covariance) > 100 * error


def test_memory_omission_has_an_independent_exterior_column_recurrence(inputs, independent):
    saved = inputs["saved"]
    assert max_abs(independent["occupation_memory_off"] - saved["occupation_memory_off"]) < 2e-12
    effect = max_abs(independent["occupation_memory_off"] - independent["occupation_full"])
    assert effect == pytest.approx(0.06952775462900784, abs=1e-11)
    assert effect > 200 * max_abs(saved["occupation_retained"] - independent["occupation_full"])


@pytest.mark.parametrize("mutation,expected", [
    ("observer_reset", "binding_observer"),
    ("initial_state_reset", "binding_initial_state"),
    ("geometry_substitution", "binding_geometry"),
    ("weights_refit", "original_weights"),
    ("full_cross_corruption", "full_cross_occupation_diagonal_movement"),
    ("probe_scalar_substitution", "full_cross_occupation_diagonal_movement"),
])
def test_audit_rejects_meaningful_source_and_cross_mutations(inputs, mutation, expected):
    altered = deepcopy(inputs)
    if mutation == "observer_reset":
        altered["observer"] = np.array(altered["columns"][:, :2], copy=True)
    elif mutation == "initial_state_reset":
        altered["columns"][0, 0] += 1e-6
    elif mutation == "geometry_substitution":
        altered["geometry"] *= 1.001
    elif mutation == "weights_refit":
        altered["weights"][0] = 0.74
    elif mutation == "full_cross_corruption":
        altered["saved"]["covariance_without_cross"][-1, 1, 1] += 0.001
    elif mutation == "probe_scalar_substitution":
        panel = altered["record"]["phases"]
        panel["full_controls"]["omissions"]["initial_cross"]["occupation_diagonal_movement"] = panel["probe_controls"]["omissions"]["initial_cross"]["occupation_diagonal_movement"]
    assert expected in saved_consistency_errors(altered)


def test_supplement_binds_the_saved_series_and_replay(inputs, independent):
    record = json.loads(OUTPUT.read_text())
    assert record["payload_npz_sha256"] == sha256(OUTPUT_NPZ)
    assert record["code_sha256"] == sha256(audit.__file__)
    assert record["source_bindings"] == {str(path.relative_to(audit.LAB)): digest for path, digest in SEALED.items()}
    assert record["cpu_seconds"] < record["cpu_budget_seconds"] <= 300.0
    assert record["payload_bytes"] < 64 * 1024 * 1024
    assert record["payload_bytes"] == OUTPUT.stat().st_size + OUTPUT_NPZ.stat().st_size
    assert record["temporal_consumer_helpers_called"] is False
    with np.load(OUTPUT_NPZ, allow_pickle=False) as stored:
        for name in ("covariance_full", "covariance_without_cross", "occupation_memory_off", "occupation_drive_off"):
            np.testing.assert_allclose(stored[name], independent[name], rtol=0.0, atol=2e-12)
        saved_cross = inputs["saved"]["covariance_retained"] - inputs["saved"]["covariance_without_cross"]
        np.testing.assert_array_equal(stored["saved_full_cross_occupation_movement"], occupation(saved_cross))
    for path, expected in SEALED.items():
        assert sha256(path) == expected
