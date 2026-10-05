"""Tests for the parent/child cut scan.

The construction asserts a parent/child nesting and the scan asks whether
the cut carries any signature. It does not, at any mass ratio or in any
field the chart carries.

The statistic matters as much as the result. The first version compared a
value at the cut against a local median, which is exactly degenerate over a
linear stretch: the median of an odd window of nodes is the centre node, so
a smooth field scores zero and so would a symmetric bump on the cut, which
is the case worth catching. These tests pin the corrected statistic so the
same blindness cannot come back.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-child-cut-v2.json"
V1 = DEVELOPMENT / "nsc-rn-child-cut-v1.json"
NOTE = LAB / "docs/nsc-rn-child-cut.md"
DRIVER = LAB / "scripts/derive_nsc_rn_child_cut.py"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestNoSignatureAtTheCut:
    def test_no_field_has_one_anywhere(self, record):
        assert record["no_feature_everywhere"]
        for entry in record["entries"]:
            assert entry["fields_with_a_feature"] == []

    def test_every_field_is_below_the_chart_maximum(self, record):
        for entry in record["entries"]:
            for block in entry["results"].values():
                assert block["curvature_excess"] < 1.0
                assert not block["detected"]

    def test_all_five_fields_were_tested(self, record):
        for entry in record["entries"]:
            assert len(entry["fields_tested"]) == 5
            assert "shift" in entry["fields_tested"]
            assert "speed_plus" in entry["fields_tested"]

    def test_all_three_mass_ratios_were_tested(self, record):
        assert len(record["entries"]) == 3
        assert all(entry["cut_in_parent"] for entry in record["entries"])


class TestTheStatisticIsNotBlind:
    def test_it_is_the_second_difference(self, record):
        assert "second difference" in record["method"]["statistic"]

    def test_the_reason_the_first_one_failed_is_recorded(self, record):
        assert "why_not_a_value_comparison" in record["method"]
        assert "median of an odd window" in record["method"][
            "why_not_a_value_comparison"]

    def test_the_first_record_is_retained(self):
        first = _load(V1)
        assert "deviation_from_local_median" in first["entries"][0]["results"][
            "speed_plus"]
        assert first["entries"][0]["results"]["speed_plus"][
            "deviation_from_local_median"] == 0.0

    def test_the_driver_no_longer_carries_the_dead_statistic(self):
        text = DRIVER.read_text(encoding="utf-8")
        body = text.split("def _entry")[1]
        assert "deviation_from_local_median" not in body


class TestWhatItMeans:
    def test_the_nesting_is_bookkeeping_in_this_chart(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "nesting is **bookkeeping**" in flat

    def test_it_says_the_condition_was_never_supplied(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "has not been supplied" in flat

    def test_it_does_not_claim_no_such_boundary_can_exist(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not show that a parent/child boundary is unreal" in flat
        assert "nothing here shows that one cannot exist" in flat

    def test_it_says_it_is_not_the_universe(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not measure the actual universe" in flat