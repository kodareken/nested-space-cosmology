"""Tests for the resolution convergence check.

The open question was whether the x1.58 saturation is a genuine asymptote or a
grid artefact, said to need more points than were available. The points were
available. These tests pin the trajectory agreement, the honest limit of two
resolutions, and the correction to the memory inference it rested on.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-resolution-convergence-v3.json"
MEMORY = DEVELOPMENT / "nsc-rn-grid-memory-v1.json"
NOTE = LAB / "docs/nsc-rn-resolution-convergence.md"
MEMORY_NOTE = LAB / "docs/nsc-rn-grid-memory.md"
NARROW_FAIL_NOTE = LAB / "docs/nsc-rn-narrow-failure.md"


@pytest.fixture(scope="module")
def record():
    return json.loads(RECORD.read_text(encoding="utf-8"))


class TestTheTrajectoriesAgree:
    def test_all_three_pairs_agree(self, record):
        assert record["verdict"] == "trajectories_agree_between_resolutions"

    def test_the_worst_difference_is_small(self, record):
        assert record["worst_relative_difference"] < 1e-4

    def test_every_case_is_compared(self, record):
        assert len(record["comparisons"]) == 3

    def test_each_pair_covers_several_stations(self, record):
        for comparison in record["comparisons"]:
            assert comparison["matched_stations"] >= 6

    def test_every_field_is_finite(self, record):
        for comparison in record["comparisons"]:
            for row in comparison["rows"]:
                assert np.isfinite(row["fine"])
                assert np.isfinite(row["coarse"])

    def test_the_runs_actually_differ_in_resolution(self, record):
        for comparison in record["comparisons"]:
            fine = DEVELOPMENT / comparison["fine_run"] / "station-0.1.npz"
            coarse = DEVELOPMENT / comparison["coarse_run"] / "station-0.1.npz"
            assert fine.is_file() and coarse.is_file()
            with np.load(fine) as a, np.load(coarse) as b:
                assert a["radius"].size > 5 * b["radius"].size


class TestTheDifferencesAccumulate:
    def test_they_grow_with_time(self, record):
        """Integration differences, not a per-step noise floor."""
        for comparison in record["comparisons"]:
            first = comparison["first_relative_difference"]
            worst = comparison["worst_relative_difference"]
            assert worst > first

    def test_the_latest_station_is_the_worst(self, record):
        for comparison in record["comparisons"]:
            rows = comparison["rows"]
            assert rows[-1]["relative_difference"] == pytest.approx(
                comparison["worst_relative_difference"])


class TestTheClockTrap:
    def test_truncated_final_stations_are_dropped(self, record):
        for comparison in record["comparisons"]:
            assert comparison["fine_stations_dropped_for_clock_offset"] == 1
            assert comparison["coarse_stations_dropped_for_clock_offset"] == 0

    def test_the_dropped_station_really_was_misnamed(self, record):
        for comparison in record["comparisons"]:
            assert comparison["worst_clock_offset"] > 0.03

    def test_the_series_is_smooth_after_the_drop(self, record):
        """The x28 jump was the trap, not the dynamics."""
        for comparison in record["comparisons"]:
            values = [r["relative_difference"] for r in comparison["rows"]][1:]
            for previous, current in zip(values, values[1:]):
                assert current < 12.0 * previous

    def test_the_driver_keys_on_the_recorded_clock(self):
        source = (LAB / "scripts/derive_nsc_rn_resolution_convergence.py").read_text(
            encoding="utf-8")
        assert "coupled_t" in source
        assert "CLOCK_TOLERANCE" in source

    def test_the_note_records_the_trap(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "nearly produced a false result" in note
        assert "0.043" in note


class TestTheLimitIsStated:
    def test_two_resolutions_bound_nothing(self, record):
        assert "Two points bound nothing" in record["scope_note"]
        assert "no convergence order" in record["scope_note"]

    def test_it_is_not_claimed_as_physics(self, record):
        assert "not evidence of physics" in record["scope_note"]
        assert "does not measure the actual universe" in record["scope_note"]

    def test_it_is_not_a_resolution_artefact_claim_in_reverse(self):
        """Agreement is evidence against an artefact, not proof of convergence."""
        assert "not a demonstration that the quantity is\nconverged" in \
            NOTE.read_text(encoding="utf-8") or \
            "not a demonstration that the quantity is" in NOTE.read_text(encoding="utf-8")


class TestTheMemoryInferenceIsWithdrawn:
    def test_memory_is_flat_in_points(self):
        memory = json.loads(MEMORY.read_text(encoding="utf-8"))
        rows = sorted(memory["rows"], key=lambda r: r["points"])
        assert rows[-1]["commit_delta_bytes"] < rows[0]["commit_total_bytes"]

    def test_the_note_says_the_inference_was_wrong(self):
        note = MEMORY_NOTE.read_text(encoding="utf-8")
        assert "never tested" in note
        assert "It was never tested" in note

    def test_both_errors_are_named_as_the_same_kind(self):
        assert "Both were wrong, and they were wrong the same way" in             NOTE.read_text(encoding="utf-8")

    def test_the_narrow_failure_note_is_corrected(self):
        note = NARROW_FAIL_NOTE.read_text(encoding="utf-8")
        assert "resources question" in note
        corrected = NOTE.read_text(encoding="utf-8")
        assert "nobody measured" in corrected
        assert "was wrong for the same reason" in corrected

class TestWhatRemainsOpen:
    def test_the_subhorizon_gap_is_still_unresolvable(self):
        memory = json.loads(MEMORY.read_text(encoding="utf-8"))
        largest = max(memory["rows"], key=lambda r: r["points"])
        assert largest["gap_in_cells"] < 1.0

    def test_the_note_says_it_was_not_for_the_given_reason(self):
        assert "It was simply not for the reason given" in \
            NOTE.read_text(encoding="utf-8")


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py",
         "scripts/derive_nsc_rn_resolution_convergence.py",
         "--pairs", "nsc-rn-inner-packet-n8000-def1040:nsc-rn-inner-packet-curv"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode == 1
    assert "need at least two" in result.stderr or \
        "refusing to overwrite" in result.stderr