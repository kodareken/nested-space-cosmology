"""Independent shared-column algebra, finalized-source and saved-control checks."""
from copy import deepcopy
import json
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

import numpy as np
import pytest

from recursive_horizons import nsc_conformal_local_response_successor as successor
from recursive_horizons.nsc_conformal_local_response import ConformalSchedule
from recursive_horizons.nsc_coupled_local_response import evolve_streamed, sha256_file
from recursive_horizons.nsc_evolving_reduction import evolve_retained_region
from recursive_horizons.nsc_spherical_galerkin_coupling import manufactured_columns


def test_shared_source_pieces_reproduce_unsplit_state_drive_and_covariance_controls():
    nf, nq, length = 16, 64, 8.
    phi0, phi1 = manufactured_columns(nf, seed=73)
    original = np.vstack((phi0, phi1))
    observer = original[:, :2]
    rotated = np.array(original, copy=True)
    theta = .4
    rotated[:, 0] = np.cos(theta) * original[:, 0] + np.sin(theta) * original[:, 2]
    rotated[:, 2] = -np.sin(theta) * original[:, 0] + np.cos(theta) * original[:, 2]
    weights = np.array([.75, .75, .5, .5, .25, .25])
    x = np.arange(nf - 1) * length / (nf - 1)
    q = .3 + .02 * np.cos(2 * np.pi * x / length)
    times = np.array([.2, .205, .21])
    schedule = ConformalSchedule({"nf": nf, "nq": nq, "length": length, "times": times,
                                 "Q": np.array([q, q * 1.005, q * 1.01])})
    shared = successor.shared_stream(schedule, observer, times, rotated, weights)
    direct, occupation, _ = evolve_streamed(schedule, observer, times, rotated, weights, propagator_substeps=2)
    _, drive_off, _ = evolve_streamed(schedule, observer, times, rotated, weights, propagator_substeps=2, outside_drive=False)
    np.testing.assert_allclose(shared["amplitudes"], direct["amplitudes"], rtol=0, atol=2e-13)
    np.testing.assert_allclose(shared["occupation"], occupation, rtol=0, atol=2e-13)
    np.testing.assert_allclose(shared["occupation_drive_off"], drive_off, rtol=0, atol=2e-13)
    covariance = (rotated * weights) @ rotated.conj().T
    dropped = evolve_retained_region(schedule, observer, times, covariance=covariance,
        drop_cross_covariance=True, backend="streamed", propagator_substeps=2)
    np.testing.assert_allclose(shared["covariance_without_cross"], dropped["covariance_total"], rtol=0, atol=2e-12)
    assert shared["result"]["history"].shape == (3, 30, 12)
    assert shared["original_column_history"].shape == (3, 30, 6)
    assert shared["result"]["allocation"]["time_indexed_exterior_propagator_bytes"] == 0


def test_successor_binds_actual_source_original_observer_and_exact_handoff():
    inputs = successor.load_inputs()
    assert inputs["times"].size == 21
    assert inputs["times"][0] == pytest.approx(.2)
    assert inputs["times"][-1] == pytest.approx(.3)
    assert not np.array_equal(inputs["observer"], inputs["columns"][:, :2])
    assert successor.initial_blocks(inputs)["C_AE_frobenius"] == pytest.approx(.338636575692358, abs=1e-12)
    physical = successor.physical_scope(inputs)
    assert physical["leader_unchanged"] is True
    assert physical["shell_leader_share_min"] > .65
    assert physical["probability_leader_share_min"] > .655
    assert [item["x"] for item in physical["surfaces"]] == [0, 2, 4]
    assert all(piece["within_one_percent"] for item in physical["surfaces"] for piece in item["indicators"].values())


@pytest.mark.parametrize("mutation,expected", [
    ("partial", "unfinished_episode"), ("reset", "source_or_state_reset"),
    ("transport", "regional_exchange_unresolved"), ("chart", "chart_or_completion:nf512_dt_0.0005"),
])
def test_failed_source_or_physical_assessment_holds_production(mutation, expected):
    episode = json.loads(successor.SOURCE_JSON.read_text())
    transport = json.loads(successor.TRANSPORT_JSON.read_text())
    bindings = successor.source_bindings()
    assert successor.readiness_errors(episode, transport, bindings) == []
    if mutation == "partial":
        episode["all_target_reached"] = False
    elif mutation == "reset":
        episode["midtrajectory_reset"] = True
    elif mutation == "transport":
        transport["resolved_regional_surface_exchange"] = False
    elif mutation == "chart":
        episode["results"][successor.PRIMARY]["positive_chart"] = False
    assert expected in successor.readiness_errors(episode, transport, bindings)


@pytest.fixture(scope="module")
def production():
    if not successor.OUTPUT_JSON.exists():
        pytest.skip("production record not yet written")
    record = json.loads(successor.OUTPUT_JSON.read_text())
    if record.get("status") == "RUNNING":
        pytest.skip("production is still running")
    with np.load(successor.OUTPUT_NPZ, allow_pickle=False) as saved:
        arrays = {name: np.array(saved[name]) for name in saved.files}
    return record, arrays, successor.load_inputs()


def test_completed_production_has_own_error_controls_and_immutable_sources(production):
    record, arrays, inputs = production
    assert successor.consistency_errors(record, arrays, inputs) == []
    assert record["status"] == "MEASURED_SAME_REALIZATION_LOCAL_RESPONSE"
    assert record["all_numerical_movements_within_one_percent"] is True
    assert record["state_error_certified"] is False
    assert record["production_geometry_rerun"] is False
    assert record["source_code_sha256"] == sha256_file(successor.__file__)
    assert record["payload_sha256"] == sha256_file(successor.OUTPUT_NPZ)
    assert record["payload_bytes"] == successor.OUTPUT_JSON.stat().st_size + successor.OUTPUT_NPZ.stat().st_size
    assert record["cpu_seconds"] < record["confirmation_cpu_ceiling"] <= 1800
    assert all(row["within_one_percent"] for row in record["numerical_controls"].values())
    for name in ("memory", "outside_drive", "initial_cross"):
        assert record["controls"][name]["occupation_movement"] > 20 * record["occupation_error"]["max_abs"]
    cross = arrays["covariance_retained"] - arrays["covariance_without_cross"]
    error = float(np.max(abs(np.diagonal(cross - arrays["cross_covariance_full"], axis1=1, axis2=2))))
    assert error == pytest.approx(record["numerical_controls"]["cross_against_independent_full"]["movement"], abs=1e-12)
    memory_error = float(np.max(abs(arrays["occupation_memory_off"] - arrays["occupation_memory_off_refined"][::2])))
    assert memory_error == record["numerical_controls"]["memory_off_temporal_indicator"]["movement"]
    drive_error = float(np.max(abs(arrays["occupation_drive_off"] - arrays["occupation_drive_off_full"])))
    assert drive_error == record["numerical_controls"]["drive_off_against_independent_full"]["movement"]
    # Each numerical comparison uses its own control effect and covariance units.
    for key, effect in (("drive_off_against_independent_full", record["controls"]["outside_drive"]["occupation_movement"]),
                       ("memory_off_temporal_indicator", record["controls"]["memory"]["occupation_movement"]),
                       ("cross_against_independent_full", record["controls"]["initial_cross"]["occupation_movement"])):
        item = record["numerical_controls"][key]
        assert item["effect"] == effect
        assert item["fraction"] == item["movement"] / effect


def test_saved_control_mutations_are_rejected(production):
    record, arrays, inputs = production
    changed = deepcopy(arrays)
    changed["covariance_without_cross"][:, 0, 0] += .001
    assert "cross_movement" in successor.consistency_errors(record, changed, inputs)
    changed = deepcopy(arrays)
    changed["occupation_drive_off"] += .001
    assert "outside_drive_movement" in successor.consistency_errors(record, changed, inputs)
    changed = deepcopy(record)
    changed["observer_id"] = "0" * 64
    assert "observer_id" in successor.consistency_errors(changed, arrays, inputs)
    changed = deepcopy(record)
    changed["state_error_certified"] = True
    assert "domain_claim" in successor.consistency_errors(changed, arrays, inputs)
