"""Tests for the nesting-chain document.

That document answers the nested-space question link by link, and it quotes
numbers from eight records. An assembled summary drifts from its sources, so
these tests read the records and check the document against them.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
DOC = LAB / "docs/nsc-rn-nesting-chain.md"
CAUSAL = DEVELOPMENT / "nsc-rn-causal-order-v4.json"
TRAPPING = DEVELOPMENT / "nsc-rn-trapping-horizon-v1.json"
MAP = DEVELOPMENT / "nsc-rn-geometry-map-v6.json"
EDGE = DEVELOPMENT / "nsc-rn-edge-sweep-v2.json"
LAYER = DEVELOPMENT / "nsc-rn-layer-shape-n8000-v1.json"
DEEP = DEVELOPMENT / "nsc-rn-inner-packet-shell-exc0760-long-v1.json"
LATE = DEVELOPMENT / "nsc-rn-late-time-v1.json"
OUTWARD = DEVELOPMENT / "nsc-rn-outward-scaling-v2.json"
BELOW = DEVELOPMENT / "nsc-rn-below-horizon-v1.json"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def doc():
    return DOC.read_text(encoding="utf-8")


class TestTheNumbersMatchTheirRecords:
    def test_the_causal_verdict(self, doc):
        assert _load(CAUSAL)["all_one_way"] is True
        assert "all nine runs" in doc

    def test_the_horizon_sits_on_r_plus(self, doc):
        record = _load(TRAPPING)
        assert record["all_below_child_start"] is True
        assert record["child_interval_is_the_trapped_region"] is False
        assert "disjoint" in doc

    def test_the_ordering_is_stable(self, doc):
        record = _load(MAP)
        assert record["summary"]["order_holds_everywhere"]
        assert record["summary"]["three_regions_everywhere"]
        assert "no transition in it" in doc

    def test_the_saturation_is_the_corrected_value(self, doc):
        """The coarse value x1.58 is superseded and must not be quoted as live."""
        coarse = _load(EDGE)["sweep"]["amplitude_gain"][-1]
        fine = _load(LAYER)
        assert coarse == pytest.approx(1.5849, abs=1e-3)
        assert "x2.05" in doc
        assert "reported as `x1.58`" in doc

    def test_the_retained_fraction(self, doc):
        value = _load(DEEP)["mass_budget"]["fraction_below_inner_edge"]
        assert 0.93 < value < 0.95
        assert "94%" in doc

    def test_the_relaxation_fraction(self, doc):
        value = _load(DEEP)["accumulation_after_arrival"]["amplitude_after_peak_fraction"]
        assert value < 0.06
        assert "2.5%" in doc

    def test_the_late_time_floor(self, doc):
        record = _load(LATE)
        assert record["turns_over_everywhere"]
        assert record["never_reaches_zero_everywhere"]
        assert "turns over" in doc and "grows again" in doc

    def test_the_shell_sides(self, doc):
        assert _load(LAYER)["shell_sides"] == [2, 2, 2, 2, 2, 3]
        assert "three sides appear" in doc

    def test_the_outward_exponent(self, doc):
        value = _load(OUTWARD)["exponents"]["1.04"]
        assert value == pytest.approx(-0.5, abs=2e-3)
        assert "r^(-1/2)" in doc

    def test_the_second_zero_exists(self, doc):
        assert _load(BELOW)["two_zeroes_everywhere"]
        assert "vanishes again at `r_minus`" in doc


class TestTheLinksAreMarked:
    def test_every_link_is_present(self, doc):
        for link in range(1, 8):
            assert f"Link {link} " in doc, link

    def test_the_open_links_are_marked_unanswerable(self, doc):
        assert "Unanswerable from the states" in doc
        assert "Refuted as posed" in doc

    def test_the_child_question_is_left_open(self, doc):
        assert "a question about labels" in doc
        assert "The data cannot" in doc

    def test_it_says_the_measured_region_is_the_trapped_one(self, doc):
        assert "inside the trapped region" in doc
        assert "it is where the trapped region is" in doc


class TestTheCorrectionsAreCarried:
    def test_the_gain_correction_is_named(self, doc):
        assert "not converged when the field it differentiates is" in doc

    def test_the_sub_cell_depth_mistake_is_named(self, doc):
        assert "under half a" in doc

    def test_the_wrong_clock_count_is_named(self, doc):
        assert "four wrong numbers" in doc

    def test_the_outward_assumption_is_named(self, doc):
        assert "assumed exponential and is" in doc


class TestScope:
    def test_it_is_not_a_scientific_result(self, doc):
        assert "not a scientific" in doc
        assert "does not measure the actual universe" in doc

    def test_it_marks_what_failed(self, doc):
        assert "What was wrong along the way" in doc
