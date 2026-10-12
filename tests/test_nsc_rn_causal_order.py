"""Tests for the causal ordering.

The question is whether time runs one way, nest to nest, never back. That is
a claim about causal structure and is answered from the characteristic speeds
in the evolved metric, not from a profile.

The delicate part is the reversal test. The first version used an absolute
tolerance near 1e-9 and so flagged a step of 1.6e-06, which is four
thousandths of a grid cell, as the horizon changing direction. These tests
pin the resolution-aware version and the reason it is needed.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-causal-order-v4.json"
V1 = DEVELOPMENT / "nsc-rn-causal-order-v1.json"
V2 = DEVELOPMENT / "nsc-rn-causal-order-v2.json"
TRAPPING = DEVELOPMENT / "nsc-rn-trapping-horizon-v1.json"
NOTE = LAB / "docs/nsc-rn-causal-order.md"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestTheOrderingIsOneWay:
    def test_every_run_is_one_way(self, record):
        assert record["all_one_way"]
        assert record["runs_one_way"] == [r["run"] for r in record["runs"]]

    def test_exactly_one_zero_everywhere(self, record):
        for run in record["runs"]:
            assert run["verdict"]["zero_count"] == 1
            assert run["verdict"]["zero_count_constant"]

    def test_the_inside_is_sealed_everywhere(self, record):
        for run in record["runs"]:
            assert run["verdict"]["sealed_everywhere"]

    def test_the_ingoing_characteristic_never_vanishes(self, record):
        """So the zero stays simple rather than becoming degenerate."""
        for run in record["runs"]:
            assert run["verdict"]["ingoing_never_vanishes"]

    def test_the_recorded_family_covers_the_whole_family(self, record):
        family = [r for r in record["runs"] if "neutral-family" in r["run"]]
        assert len(family) == 6
        assert {r["mu"] for r in family} == {1.002, 1.010, 1.04}

    def test_it_covers_more_than_one_mass_ratio(self, record):
        assert len({r["mu"] for r in record["runs"]}) >= 2


class TestTheReversalTestIsResolutionAware:
    def test_the_motions_are_reported_in_cells(self, record):
        for run in record["runs"]:
            assert "total_excursion_in_cells" in run["verdict"]
            assert "largest_inward_step_in_cells" in run["verdict"]

    def test_the_threshold_is_half_a_cell(self, record):
        for run in record["runs"]:
            verdict = run["verdict"]
            assert verdict["reversal_threshold"] == pytest.approx(
                0.5 * verdict["grid_spacing"])

    def test_a_sub_cell_step_is_not_a_reversal(self, record):
        """The recorded-family horizon is static; the inner-packet one moves
        by a tenth of a cell and neither counts as reversing."""
        for run in record["runs"]:
            verdict = run["verdict"]
            if verdict["largest_inward_step"] > 0.0:
                assert verdict["largest_inward_step_in_cells"] < 0.5
                assert not verdict["resolved_reversal"]

    def test_the_first_records_are_retained_and_distinguishable(self):
        """v2 used an absolute tolerance and reported a reversal."""
        v2 = _load(V2)
        inner = [r for r in v2["runs"] if "inner-packet" in r["run"]]
        assert inner and not inner[0]["verdict"]["horizon_moves_outward_only"]
        assert _load(V1)["runs"][0]["verdict"]["one_way"]

    def test_the_note_explains_the_first_tolerance(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "absolute `1e-9`" in note
        assert "sub-resolution variation" in note


class TestTheHorizonIsNotTheEdge:
    def test_the_edge_is_reported_separately(self, record):
        assert record["runs"][0]["rows"][0]["r_inner"] is not None

    def test_the_note_says_the_edge_is_not_a_horizon(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "inside the trapped region" in note
        assert "cannot escape" in note

    def test_the_criterion_is_never_copied_from_the_reference(self, record):
        trapping = _load(TRAPPING)
        assert trapping["criterion"]["event_horizon_copied_from_trapping"] is False


class TestTheLimits:
    def test_it_does_not_claim_the_nesting_iterates(self, record):
        assert "defines one child interval" in record["what_this_cannot_answer"]
        assert "no claim about one is made" in record["what_this_cannot_answer"]

    def test_it_does_not_claim_a_global_event_horizon(self, record):
        assert "not a proof of a global event horizon" in record["scope_note"]

    def test_it_says_it_is_not_the_universe(self, record):
        assert "does not measure the actual universe" in record["scope_note"]

    def test_the_horizon_motion_is_marked_unresolved(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "0.12" in note
        assert "no claim is made about" in note


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_causal_order.py",
         "--runs", "nsc-rn-inner-packet-def1040-long:1.04"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr