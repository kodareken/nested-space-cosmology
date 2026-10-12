"""Tests for the compression-bounding re-reading.

The saturation was described as deepening, which says where and not why.
The re-reading puts the gain against the room actually available inside the
trapped region, and finds the ceiling sits where the boundary condition
reverses.

These tests pin the arithmetic and, importantly, pin that this is a
structural limit rather than a resolution one, since that is the claim being
made about the earlier saturation.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-compression-bounding-v2.json"
MAP = DEVELOPMENT / "nsc-rn-geometry-map-v6.json"
EDGE = DEVELOPMENT / "nsc-rn-edge-sweep-v2.json"
FINE = DEVELOPMENT / "nsc-rn-layer-shape-n8000-v1.json"
NOTE = LAB / "docs/nsc-rn-compression-bounding.md"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestTheArithmetic:
    def test_the_percentages_span_the_trapped_width(self, record):
        fractions = [r["fraction_of_trapped_width_used"] for r in record["rows"]]
        assert min(fractions) == pytest.approx(0.50, abs=0.01)
        assert max(fractions) == pytest.approx(0.998, abs=0.005)

    def test_the_gains_come_from_the_fine_records(self, record):
        """Not from the coarse sweep: the two differ by up to thirty percent."""
        coarse = dict(zip(_load(EDGE)["sweep"]["r_inner"],
                          _load(EDGE)["sweep"]["amplitude_gain"]))
        for row in record["rows"]:
            assert row["gain"] > coarse[row["r_inner"]] * 0.99

    def test_the_record_declares_which_records_it_read(self):
        text = (LAB / "scripts/derive_nsc_rn_compression_bounding.py").read_text(
            encoding="utf-8")
        assert "nsc-rn-inner-packet-shell-%s-v1.json" in text
        assert "mix two resolutions" in text

    def test_the_horizons_match_the_map(self, record):
        chart = [e for e in _load(MAP)["entries"] if abs(e["mu"] - 1.04) < 1e-9][0]
        assert record["r_minus"] == pytest.approx(chart["r_minus"])
        assert record["r_plus"] == pytest.approx(chart["r_plus"])

    def test_the_gain_turns_over_at_the_far_end(self, record):
        assert record["turns_over_at_the_far_end"]

    def test_it_does_not_rise_all_the_way(self, record):
        """It rises and then reverses; asserting a clean rise would be wrong."""
        assert record["gain_rises_with_distance_to_horizon"] is False


class TestTheCeilingIsStructural:
    def test_the_deepest_edge_is_at_the_boundary_condition(self, record):
        assert record["deepest_edge_distance_to_inner_horizon"] < 0.002

    def test_the_peak_is_before_the_turn_over(self, record):
        rows = sorted(record["rows"], key=lambda r: r["fraction_of_trapped_width_used"])
        peak = max(rows, key=lambda r: r["gain"])
        assert peak["r_inner"] == 0.77
        assert peak["fraction_of_trapped_width_used"] < 0.98

    def test_it_is_not_attributed_to_resolution(self, record):
        flat = record["reading"].replace(chr(10), " ")
        assert "not because of the grid" in flat
        assert record["why_it_matters"].startswith("The saturation was previously")


class TestWhatItDoesNotClaim:
    def test_it_does_not_claim_maximality_in_general(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "adds no new run" in flat

    def test_the_note_bounds_it_to_this_construction(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "for an excision-based construction on this chart" in flat

    def test_it_names_what_would_not_be_bounded(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not excise at all, is not bounded here" in flat