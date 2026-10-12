"""Tests for the stability audit of the deeper-edge domains.

The first version of this audit required the mass function to be
non-decreasing everywhere and flagged about 1800 steps per run, including in
the default-edge run that predates the whole exercise. The reversals were all
at r >= 19.7, in the far field the oscillation audit already characterised as
grid-scale ringing. The criterion was wrong, not the march.

These tests pin both the corrected criterion and the distinction that caused
it, so the same error is not made again with a different radius.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
STABILITY = DEVELOPMENT / "nsc-rn-edge-stability-v2.json"
SUPERSEDED = DEVELOPMENT / "nsc-rn-edge-stability-v1.json"
DRIVER = LAB / "scripts/derive_nsc_rn_edge_stability.py"


@pytest.fixture(scope="module")
def audit():
    return json.loads(STABILITY.read_text(encoding="utf-8"))


class TestTheDomainsAreStable:
    def test_every_deeper_domain_is_stable(self, audit):
        assert audit["all_runs_stable"]

    def test_it_covers_the_default_edge_and_the_whole_sweep(self, audit):
        edges = [run["r_inner"] for run in audit["runs"]]
        assert len(edges) == 4
        assert min(edges) <= 0.77

    def test_every_recorded_field_is_finite(self, audit):
        for run in audit["runs"]:
            assert run["all_fields_finite"]
            assert run["non_finite_counts"] == {}

    def test_the_lapse_stays_positive(self, audit):
        for run in audit["runs"]:
            assert run["lapse_always_positive"]

    def test_the_radius_grid_stays_ordered(self, audit):
        for run in audit["runs"]:
            assert run["radius_increasing"]

    def test_every_run_reached_its_target(self, audit):
        assert all(run["stop_reason"] == "target" for run in audit["runs"])


class TestTheCriterionIsRegional:
    def test_no_reversals_below_the_limit_anywhere(self, audit):
        for run in audit["runs"]:
            assert run["mass_decreases_below_limit"] == 0

    def test_the_far_field_reversals_are_present_and_identified(self, audit):
        for run in audit["runs"]:
            assert run["mass_reversals_far_field"] > 0
            assert run["far_field_reversals_are_known_ringing"]

    def test_every_far_field_reversal_is_above_the_limit(self, audit):
        for run in audit["runs"]:
            assert run["nearest_far_field_reversal_radius"] >= run["monotonic_below_radius"]

    def test_the_limit_is_declared(self, audit):
        assert all(run["monotonic_below_radius"] == 16.0 for run in audit["runs"])

    def test_the_default_edge_run_shows_the_same_far_field_behaviour(self, audit):
        """The ringing is a property of the family, not of the deep edges."""
        counts = [run["mass_reversals_far_field"] for run in audit["runs"]]
        assert min(counts) > 1500
        assert max(counts) < 2000


class TestAgainstRepeatingTheMistake:
    def test_the_superseded_record_is_retained(self):
        first = json.loads(SUPERSEDED.read_text(encoding="utf-8"))
        assert first["all_runs_stable"] is False
        assert "mass_decreases" in first["runs"][0]

    def test_the_driver_says_why_the_criterion_is_regional(self):
        source = DRIVER.read_text(encoding="utf-8")
        assert "criterion was" in source
        assert "wrong rather than the march" in source

    def test_the_driver_names_the_already_measured_ringing(self):
        source = DRIVER.read_text(encoding="utf-8")
        assert "oscillation audit" in source


class TestScope:
    def test_stability_is_not_correctness(self, audit):
        assert "Stability is not correctness" in audit["scope_note"]

    def test_it_says_it_is_not_the_universe(self, audit):
        assert "does not measure the actual universe" in audit["scope_note"]


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_edge_stability.py",
         "--runs", "nsc-rn-inner-packet-exc085", "nsc-rn-inner-packet-exc080"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr