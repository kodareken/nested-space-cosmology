"""Tests for what the child is.

The nesting question asked whether the sealed region is the child. The
construction declares the answer in prose, and it is not the one the numbers
suggested: the child is an interval above the outer horizon, and the trapped
region is something the construction calls the excision.

These tests read the declarations and the map, so the answer cannot drift
away from what the code actually says.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
import derive_nsc_rn_parent as driver  # noqa: E402

REFERENCE = LAB / "src/recursive_horizons/nsc_rn_reference.py"
MAP = LAB / "results/development/nsc-rn-geometry-map-v6.json"
NOTE = LAB / "docs/nsc-rn-what-the-child-is.md"
MU = 1.04


@pytest.fixture(scope="module")
def declarations():
    text = REFERENCE.read_text(encoding="utf-8")
    out = {}
    for line in text.split("\n"):
        for key in ("excision", "parent", "child"):
            if line.strip().startswith('"%s":' % key):
                out.setdefault(key, line.strip())
    return out


@pytest.fixture(scope="module")
def regions():
    record = json.loads(MAP.read_text(encoding="utf-8"))
    entry = [e for e in record["entries"] if e["mu"] == MU][0]
    return entry


class TestTheChildIsNotTheInterior:
    def test_the_child_sits_above_the_outer_horizon(self, regions):
        assert regions["child_start"] > regions["r_plus"]

    def test_the_trapped_region_is_below_the_horizon(self, regions):
        assert (regions["r_minus"] < regions["excision_midpoint"]
                < regions["r_plus"])
        assert regions["trapped_region"] == [regions["r_minus"], regions["r_plus"]]

    def test_the_child_and_the_trapped_region_are_disjoint(self, regions):
        assert regions["child_start"] > regions["r_plus"]

    def test_the_gap_between_horizon_and_child_is_a_gap_sixteenth(self, regions):
        gap = regions["r_plus"] - regions["r_minus"]
        assert regions["child_start"] - regions["r_plus"] == pytest.approx(
            gap / 16.0, rel=1e-3)

    def test_the_nesting_order_is_parent_child_horizon_trapped(self, regions):
        assert (regions["child_end"] > regions["child_start"] > regions["r_plus"]
                > regions["excision_midpoint"] > regions["r_minus"])


class TestTheConstructionSaysSo:
    def test_the_child_is_declared_near_the_outer_horizon(self, declarations):
        assert "near the outer horizon" in declarations["child"]

    def test_the_parent_is_declared_the_exterior(self, declarations):
        assert "large-r exterior" in declarations["parent"]

    def test_the_excision_is_declared_between_the_horizons(self, declarations):
        assert "between r_minus and r_plus" in declarations["excision"]

    def test_the_note_matches_the_declarations(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "not the interior" in flat
        assert "calls the" in flat and "excision" in flat


class TestTheGapThatRemains:
    def test_the_parent_child_cut_has_no_physical_signature(self):
        """It is a label. That is the load-bearing gap now."""
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "label with no physical signature" in flat

    def test_the_horizon_is_below_the_child_not_at_its_boundary(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "membrane between the child and the trapped region" in flat

    def test_it_does_not_claim_the_child_holds_a_universe(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "not an occupancy" in flat

    def test_it_says_it_is_not_the_universe(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not measure the actual universe" in flat