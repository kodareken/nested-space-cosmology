"""Authentication, timestamp joins, and bounded existing-owner metric replay."""
import copy
import json

import pytest

import assess_nsc_discovery_crossing as consumer


def test_handoff_rejects_partial_dictionary_and_predecessor_mismatch():
    original = {"case_id": "case", "coordinate_time": 3.0, "ordinal": 2,
                "arrays_sha256": "payload", "array_sha256": {"Q": "q", "normal_clocks": "clock"}}
    handoff = {**original, "ordinal": 0,
               "predecessor_chunk": {"case_id": "case", "ordinal": 2, "arrays_sha256": "payload"}}
    assert consumer.validate_handoff(handoff, original)["no_state_momentum_clock_basis_reset"]
    bad = copy.deepcopy(handoff)
    bad["array_sha256"].pop("normal_clocks")
    with pytest.raises(ValueError, match="hash dictionaries differ"):
        consumer.validate_handoff(bad, original)
    bad = copy.deepcopy(handoff)
    bad["predecessor_chunk"]["arrays_sha256"] = "another"
    with pytest.raises(ValueError, match="predecessor or time mismatch"):
        consumer.validate_handoff(bad, original)


def test_timestamp_join_does_not_zip_across_missing_cadence():
    left = [{"time": 5.8}, {"time": 5.85}, {"time": 7.95}, {"time": 8.0}]
    right = [{"time": 5.799999999999002}, {"time": 8.0}]
    pairs, coverage = consumer.match_observations(left, right)
    assert [(round(a["time"], 7), round(b["time"], 7)) for a, b in pairs] == [(5.8, 5.8), (8.0, 8.0)]
    assert coverage["left_only_times"] == [5.85, 7.95]
    gap = consumer.cadence_gaps(right)
    assert len(gap) == 1 and gap[0]["right"] == 8.0
    with pytest.raises(ValueError, match="duplicate"):
        consumer.match_observations(left, right + [{"time": 8.0}])


@pytest.fixture(scope="module")
def report():
    return consumer.assess()


def test_authenticated_positive_metric_and_existing_curvature_controls(report):
    assert report["cpu_seconds"] < 30.0
    assert report["input_sources_unchanged"] is True
    assert report["bindings"]["crossing"]["producing_commit"] == consumer.CROSSING_COMMIT
    assert "post-run" in report["bindings"]["crossing"]["binding_timing"]
    assert report["bindings"]["original_episode"]["producing_commit"].startswith("1c9e770")
    assert all(row["positive_metric"] for row in report["rows"])
    assert all(row["no_state_momentum_clock_basis_reset"] for row in report["T3_handoffs"].values())
    assert all(not consumer.Path(path).is_absolute() for path in report["source_hashes"] | report["input_hashes"])
    exact = consumer.tidal.exact_suite()
    assert exact["minkowski"]["flat_curvature_max_abs"] == 0.0
    assert exact["bertotti_robinson"]["Ricci2"] > 0.0
    assert exact["bertotti_robinson"]["K"] > 0.0
    assert report["scope"]["R4_Rh_Weyl_resolved"] is False
    assert report["scope"]["ambient_echo_proven"] is False
    half = next(row for row in report["rows"] if row["case_id"] == "nf256_coupled_dt0.0005")
    assert half["worldline_x2"]["R_0101"] == pytest.approx(-9.665013073e22, rel=1e-8)
    assert half["worldline_x2"]["R_0202"] == pytest.approx(3.067404935e24, rel=1e-8)
    assert half["cadence_gaps"][-1]["right"] == 8.0


def test_creation_and_readonly_replay_reject_changed_bound_hash(report, tmp_path, monkeypatch):
    destination = tmp_path / "assessment.json"
    monkeypatch.setattr(consumer, "OUTPUT", destination)
    consumer.write_record(report, destination)
    before = destination.read_bytes()
    assert consumer.check_record(destination)["recomputed_metric_jets"]
    assert destination.read_bytes() == before
    with pytest.raises(FileExistsError):
        consumer.write_record(report, destination)
    altered = json.loads(before)
    path = next(iter(altered["input_hashes"]))
    altered["input_hashes"][path] = "0"*64
    destination.write_text(json.dumps(altered))
    with pytest.raises(ValueError, match="bound hash changed"):
        consumer.check_record(destination)
