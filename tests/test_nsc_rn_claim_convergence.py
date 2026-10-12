"""Tests for the claim-convergence check.

The compression gain was corrected from x1.58 to x2.05 because a gain read
off dM/dr is not converged when the field it differentiates is. That
correction was applied to the gain alone, and three other headline numbers
had never been rechecked under the same refinement.

These tests pin what the check found, including the part that is not clean:
the retained fraction converges at one edge and drifts by about one percent
at the two deepest. A test that only asserted convergence would have hidden
that, which is the mistake the check exists to prevent.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-claim-convergence-v2.json"
V1 = DEVELOPMENT / "nsc-rn-claim-convergence-v1.json"
GAIN_CONVERGENCE = DEVELOPMENT / "nsc-rn-resolution-convergence-v3.json"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestTheComparisonIsMatchedInTime:
    def test_every_row_is_matched(self, record):
        for comparison in record["comparisons"]:
            for row in comparison["rows"]:
                assert row["matched"]

    def test_it_says_why(self):
        flat = _load(RECORD)["what_this_measures"]
        assert "grows with time" in flat

    def test_three_edges_were_compared(self, record):
        assert len(record["comparisons"]) == 3


class TestWhatItFound:
    def test_it_converges_at_the_shallowest_of_the_three(self, record):
        first = record["comparisons"][0]
        assert first["worst_relative_difference"] < 1e-3
        assert first["converged"]

    def test_it_drifts_at_the_two_deepest(self, record):
        """This is the part a bare verdict would have hidden."""
        for comparison in record["comparisons"][1:]:
            assert not comparison["converged"]
            assert 1e-3 < comparison["worst_relative_difference"] < 1e-1

    def test_the_drift_is_small_and_the_record_says_so(self, record):
        assert record["worst_relative_difference"] < 0.02
        assert "drifting, not wrong by much" in record["reading_the_verdict"]

    def test_it_is_contrasted_with_the_gain(self, record):
        assert "thirty percent" in record["reading_the_verdict"]
        assert "thirty percent the gain moved" in record["reading_the_verdict"]


class TestWhatItDoesNotCover:
    def test_it_says_the_floor_and_relaxation_are_not_addressed(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "late-time floor or the relaxation" in flat
        assert "longer run" in flat

    def test_the_superseded_record_is_retained(self):
        first = _load(V1)
        assert "all_converged" in first
        assert "worst_relative_difference" not in first


class TestConsistencyWithTheEarlierFinding:
    def test_the_trajectory_convergence_is_still_on_record(self):
        assert _load(GAIN_CONVERGENCE)["verdict"] == \
            "trajectories_agree_between_resolutions"

    def test_the_retained_fraction_uses_that_trajectory(self, record):
        """It is a mass-stock quantity, which is why it is close to converged."""
        flat = record["what_this_measures"]
        assert "read off a coarse grid" in flat