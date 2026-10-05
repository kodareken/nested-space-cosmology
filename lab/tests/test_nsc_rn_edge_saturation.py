"""Tests for the edge saturation.

The previous note claimed the compression gain was still rising and had not
turned over. It had: the gain saturates near x1.58, and the apparent
non-convergence was a sub-cell boundary read as depth. These tests pin the
saturation, the resolution reading that explains it, and the correction.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
SATURATION = DEVELOPMENT / "nsc-rn-edge-sweep-v2.json"
SUPERSEDED = DEVELOPMENT / "nsc-rn-edge-sweep-v1.json"
NOTE = LAB / "docs/nsc-rn-edge-saturation.md"
PREVIOUS_NOTE = LAB / "docs/nsc-rn-deep-long.md"
DRIVER = LAB / "scripts/derive_nsc_rn_edge_sweep.py"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def sweep():
    return _load(SATURATION)


@pytest.fixture(scope="module")
def block(sweep):
    return sweep["sweep"]


class TestTheGainSaturates:
    def test_it_is_monotone_across_the_whole_sweep(self, block):
        assert block["gain_increases_with_depth"]

    def test_it_flattens_at_the_deepest_edges(self, block):
        gains = block["amplitude_gain"]
        assert gains[-1] - gains[-2] < 0.01
        assert gains[-1] - gains[-2] < 0.1 * (gains[1] - gains[0])

    def test_the_saturation_value_is_about_one_and_a_half(self, block):
        assert block["amplitude_gain"][-1] == pytest.approx(1.58, abs=0.02)

    def test_the_midpoint_is_a_lower_bound_by_about_a_third(self, block):
        first, last = block["amplitude_gain"][0], block["amplitude_gain"][-1]
        assert last / first == pytest.approx(1.32, rel=0.05)


class TestTheResolutionReading:
    def test_the_gap_is_reported_in_cells(self, block):
        assert len(block["gap_in_cells"]) == len(block["r_inner"])

    def test_the_deepest_edges_are_sub_cell(self, block):
        assert min(block["gap_in_cells"]) < 0.5

    def test_three_of_six_are_not_resolvable(self, block):
        assert block["edge_resolvable_from_horizon"].count(False) == 3

    def test_the_gain_still_rises_where_it_is_resolvable(self, block):
        assert block["gain_increases_with_depth_where_resolvable"]

    def test_the_record_says_it_was_resolution_limited(self, block):
        assert "resolution-limited before it is depth-limited" in block["resolution_limit"]

    def test_the_driver_carries_the_cells_field(self):
        source = DRIVER.read_text(encoding="utf-8")
        assert "gap_in_cells" in source
        assert "edge_resolvable_from_horizon" in source


class TestTheRetainedFractionIsEdgeIndependent:
    def test_it_is_flat_to_three_figures(self, block):
        retained = block["fraction_below_inner_edge"]
        assert max(retained) - min(retained) < 0.005

    def test_it_sits_near_ninety_four_percent(self, block):
        assert all(0.93 < value < 0.95 for value in block["fraction_below_inner_edge"])


class TestTheRelaxationIsNotATrend:
    def test_it_does_not_order_with_depth(self, block):
        ratios = block["R4_edge_final_over_initial"]
        assert not block["relaxation_weakens_with_depth"]

    def test_it_scatters_over_a_wide_band(self, block):
        ratios = block["R4_edge_final_over_initial"]
        assert max(ratios) - min(ratios) > 0.3

    def test_the_note_declines_to_call_it_a_trend(self):
        assert "not a trend and none is claimed" in NOTE.read_text(encoding="utf-8")


class TestTheCorrection:
    def test_the_previous_claim_is_withdrawn_in_the_new_note(self):
        assert "That was wrong" in NOTE.read_text(encoding="utf-8")

    def test_the_new_note_says_what_the_old_one_got_wrong(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "resolution-limited before it was depth-limited" in note
        assert "different statement from" in note

    def test_the_superseded_record_is_retained(self):
        first = _load(SUPERSEDED)
        assert "gap_in_cells" not in first["sweep"]


class TestTheGuardHeld:
    def test_it_refuses_a_mixed_target(self):
        import subprocess
        import sys
        result = subprocess.run(
            [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_edge_sweep.py",
             "--runs", "nsc-rn-inner-packet-curv",
             "nsc-rn-inner-packet-exc085-long",
             "--out", "_tmp-edge-mixed.json"],
            cwd=LAB.parent, capture_output=True, text=True)
        assert result.returncode == 2
        assert "not an edge sweep" in result.stderr

    def test_the_note_records_that_the_guard_fired(self):
        assert "This is not an edge sweep" in NOTE.read_text(encoding="utf-8")


class TestScope:
    def test_every_edge_is_still_a_boundary(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "still a **boundary**" in note
        assert "does not make them horizons" in note

    def test_it_says_it_is_not_the_universe(self):
        assert "does not measure the\nactual universe" in NOTE.read_text(encoding="utf-8")

    def test_it_says_the_saturation_is_a_grid_property(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "not a property of anything else" in note

    def test_it_names_the_open_question_as_resolution(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "resolution\nquestion, not a depth one" in note or \
            "resolution question, not a depth one" in note


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_edge_sweep.py",
         "--runs", "nsc-rn-inner-packet-exc085-long", "nsc-rn-inner-packet-exc080-long",
         "--out", "nsc-rn-edge-sweep-v2.json"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr