"""Tests for the inner-edge sweep.

The sweep asks whether the effect of deepening the grid's inner edge
continues or turns over. It continues for the compression gain, and it does
not for the relaxation, for a reason that is itself a finding: the deeper
runs have not finished accumulating.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
SWEEP = DEVELOPMENT / "nsc-rn-edge-sweep-v1.json"
EXCISION_NOTE = LAB / "docs/nsc-rn-excision-move.md"
EDGE_NOTE = LAB / "docs/nsc-rn-edge-sweep.md"
DRIVER = LAB / "scripts/derive_nsc_rn_edge_sweep.py"


@pytest.fixture(scope="module")
def sweep():
    return json.loads(SWEEP.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def block(sweep):
    return sweep["sweep"]


class TestTheSweepIsAnEdgeSweep:
    def test_it_spans_distinct_edges_in_ascending_order(self, block):
        edges = block["r_inner"]
        assert edges == sorted(edges, reverse=True)
        assert len(set(edges)) == len(edges)

    def test_it_reaches_close_to_the_inner_horizon(self, block):
        assert min(block["gap_to_inner_horizon"]) < 0.02

    def test_everything_else_is_held_fixed(self, sweep):
        for key, values in sweep["held_fixed"].items():
            assert len(values) == 1, (key, values)

    def test_the_default_edge_is_included_and_labelled(self, sweep):
        assert any(r["is_default_edge"] for r in sweep["runs"])
        assert not any(r["is_default_edge"] for r in sweep["runs"]
                       if r["excision_requested"] is not None)

    def test_every_run_reached_its_target(self, sweep):
        assert all(r["stop_reason"] == "target" for r in sweep["runs"])


class TestTheCompressionTrend:
    def test_the_gain_increases_monotonically_with_depth(self, block):
        assert block["gain_increases_with_depth"]

    def test_the_gain_is_still_growing_at_the_deepest_edge(self, block):
        assert block["amplitude_gain"][-1] > block["amplitude_gain"][0]

    def test_the_midpoint_gain_is_a_lower_bound_by_a_wide_margin(self, block):
        """The headline x1.20 was truncated by about a third."""
        assert block["amplitude_gain"][0] < 1.25
        assert block["amplitude_gain"][-1] > 1.5
        assert block["amplitude_gain"][-1] / block["amplitude_gain"][0] > 1.25


class TestTheRelaxationIsNotATrend:
    def test_the_relaxation_ratio_is_not_monotone(self, block):
        assert not block["relaxation_weakens_with_depth"]

    def test_the_deep_runs_exceed_one_and_have_not_relaxed(self, block):
        ratios = block["R4_edge_final_over_initial"]
        assert any(ratio > 1.0 for ratio in ratios)

    def test_the_note_withdraws_the_weakening_claim(self):
        note = EDGE_NOTE.read_text(encoding="utf-8")
        assert "withdrawn" in note
        assert "too short to measure a relaxation" in note

    def test_the_verdict_is_about_the_relaxation_only(self, sweep, block):
        assert block["verdict"] == "effect_does_not_continue_monotonically"
        assert block["gain_increases_with_depth"]
        assert "not a failure of the trend, it is the trend" in EDGE_NOTE.read_text(encoding="utf-8")


class TestTheConfoundIsNamed:
    def test_the_deeper_edge_starts_in_more_curved_geometry(self, sweep):
        firsts = [r["R4_edge_first"] for r in sweep["runs"]]
        assert firsts[-1] > 8.0 * firsts[0]

    def test_the_note_says_the_two_causes_are_not_separated(self):
        note = EDGE_NOTE.read_text(encoding="utf-8")
        assert "does not separate the two" in note
        assert "matched comparison" in note


class TestCurvatureStaysBounded:
    def test_no_edge_diverges(self, sweep):
        for run in sweep["runs"]:
            assert run["R4_edge_max_over_initial"] < 10.0


class TestScope:
    def test_every_edge_is_still_a_boundary(self):
        note = EDGE_NOTE.read_text(encoding="utf-8")
        assert "still a **boundary**" in note
        assert "is still a boundary" in note

    def test_it_says_it_is_not_the_universe(self, sweep):
        assert "does not measure the actual universe" in sweep["scope_note"]

    def test_it_says_the_effect_is_not_converged(self):
        note = EDGE_NOTE.read_text(encoding="utf-8")
        assert "has not converged with respect to domain depth" in note

    def test_the_previous_note_points_at_this_one(self):
        previous = EXCISION_NOTE.read_text(encoding="utf-8")
        assert "`0.80` and `0.77`" in previous


class TestTheDriver:
    def test_it_names_runs_explicitly(self):
        assert '"--runs", nargs="+", required=True' in DRIVER.read_text(encoding="utf-8")

    def test_it_refuses_runs_that_differ_in_something_else(self):
        import subprocess
        import sys
        result = subprocess.run(
            [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_edge_sweep.py",
             "--runs", "nsc-rn-inner-packet-curv", "nsc-rn-inner-packet-c3-fine",
             "--out", "_tmp-edge-guard.json"],
            cwd=LAB.parent, capture_output=True, text=True)
        assert result.returncode == 2
        assert "not an edge sweep" in result.stderr

    def test_it_refuses_to_overwrite(self):
        import subprocess
        import sys
        result = subprocess.run(
            [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_edge_sweep.py",
             "--runs", "nsc-rn-inner-packet-exc085", "nsc-rn-inner-packet-exc080",
             "--out", "nsc-rn-edge-sweep-v1.json"],
            cwd=LAB.parent, capture_output=True, text=True)
        assert result.returncode != 0
        assert "refusing to overwrite" in result.stderr