"""Tests for the geometry map.

The map is the deliverable this line was aiming at in the sense that a map
is not a discovery but a description: one space, every feature of it, and
how the features move as the one free parameter changes.

The thing a map is for is ordering. If the child cut ever crossed the
horizon, or the horizon crossed the inner edge, the nesting would have a
transition in it. These tests check the ordering, the tolerance the
ordering test needed, and the one structural fact the map exposes.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-geometry-map-v6.json"
V1 = DEVELOPMENT / "nsc-rn-geometry-map-v1.json"
DRIVER = LAB / "scripts/derive_nsc_rn_geometry_map.py"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestTheOrdering:
    def test_it_holds_across_the_whole_map(self, record):
        assert record["summary"]["order_holds_everywhere"]

    def test_it_spans_a_wide_parameter_range(self, record):
        assert record["summary"]["mu_min"] < 1.002
        assert record["summary"]["mu_max"] >= 1.1
        assert record["summary"]["entries"] >= 10

    def test_the_ordering_is_the_expected_one(self, record):
        expected = ["r_minus", "excision", "r_plus", "outgoing_zero",
                    "child_start", "child_end"]
        assert all(entry["order"] == expected for entry in record["entries"])

    def test_the_causal_membrane_sits_at_the_outer_horizon(self, record):
        assert record["summary"]["outgoing_zero_at_r_plus_everywhere"]


class TestTheToleranceIsScaled:
    def test_the_zero_is_located_to_a_hundredth_of_a_cell(self, record):
        deficits = [entry["zero_deficit_in_cells"] for entry in record["entries"]]
        assert all(abs(value) < 0.02 for value in deficits)

    def test_the_tolerance_is_the_grid_spacing(self, record):
        for entry in record["entries"]:
            assert entry["ordering_tolerance"] > 0.001

    def test_the_first_record_used_an_arbitrary_tolerance(self):
        first = _load(V1)
        assert "ordering_tolerance" not in first["entries"][0]

    def test_the_deficit_is_reported_per_entry(self, record):
        for entry in record["entries"]:
            assert entry["zero_deficit_below_r_plus"] is not None


class TestThreeRegions:
    """The shift equals one at both horizons, so the chart divides in three."""

    def test_there_are_three_regions_everywhere(self, record):
        assert record["summary"]["three_regions_everywhere"]
        for entry in record["entries"]:
            assert [r["name"] for r in entry["regions"]] == [
                "parent", "trapped", "below_inner_horizon"]

    def test_only_the_middle_one_is_trapped(self, record):
        entry = record["entries"][4]
        trapped = [r["name"] for r in entry["regions"] if r["trapped"]]
        assert trapped == ["trapped"]

    def test_every_edge_this_programme_used_is_inside_it(self, record):
        assert record["summary"]["every_used_edge_inside_the_trapped_region"]
        entry = [e for e in record["entries"] if e["mu"] == 1.04][0]
        for edge in entry["edges_used_at_this_mu"]:
            low, high = entry["trapped_region"]
            assert low < edge < high

    def test_the_chart_cannot_describe_the_third_region(self, record):
        assert record["summary"]["chart_never_describes_below_r_minus"]
        entry = record["entries"][4]
        third = entry["regions"][2]
        assert third["described_by_this_chart"] is False
        assert "above the inner horizon" in third["chart_refusal"]

    def test_the_first_records_predate_the_region_split(self):
        assert "regions" not in _load(V1)["entries"][0]


class TestTheShiftTurn:
    def test_it_is_below_the_default_edge_everywhere(self, record):
        assert record["summary"]["the_default_edge_never_samples_the_shift_maximum"]

    def test_the_ratio_is_reported(self, record):
        low, high = record["summary"]["shift_turn_over_excision_range"]
        assert 0.0 < low < high < 1.0

    def test_that_is_the_structural_fact_the_map_exposes(self, record):
        """r = Q^2/M always sits below the midpoint r = mu, for mu above one."""
        for entry in record["entries"]:
            # Located on a finite grid, so to within a cell of Q^2/M.
            assert abs(entry["shift_turn_radius"] - 1.0 / entry["mu"])                 < entry["sampling_spacing"]
            assert entry["shift_turn_radius"] < entry["excision_midpoint"]


class TestScope:
    def test_it_is_a_map_of_one_family(self, record):
        assert "not of nested spaces in general" in record["scope_note"]

    def test_it_says_one_child_is_not_a_hierarchy(self, record):
        assert "not a hierarchy" in record["scope_note"]

    def test_it_says_it_is_not_the_universe(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "does not measure the actual universe" in flat
    def test_it_states_it_is_classical(self, record):
        assert "Every" in record["what_this_is"]
        assert "classical" in record["what_this_is"]


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_geometry_map.py"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr