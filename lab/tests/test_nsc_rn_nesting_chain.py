"""Tests for the nesting-chain document.

That document walks the nested-space question link by link and marks each
link measured, unanswerable, or refuted. It quotes numbers from several
records, so it can drift from them. These tests read the records and check
the document against them, which is the only thing that keeps an assembled
summary honest.
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
SETTLE = DEVELOPMENT / "nsc-rn-packet-sweep-settle-v3.json"
EDGE = DEVELOPMENT / "nsc-rn-edge-sweep-v2.json"
LAYER = DEVELOPMENT / "nsc-rn-layer-shape-v4.json"
EXCISION = DEVELOPMENT / "nsc-rn-excision-move.md" if False else \
    LAB / "docs/nsc-rn-edge-saturation.md"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def doc():
    return DOC.read_text(encoding="utf-8")


class TestTheQuotedNumbersMatchTheRecords:
    def test_the_causal_verdict_matches(self, doc):
        record = _load(CAUSAL)
        assert record["all_one_way"] is True
        assert record["runs_total"] == 9
        assert "one-way in all nine runs" in doc

    def test_the_trapping_verdict_matches(self, doc):
        record = _load(TRAPPING)
        assert record["all_below_child_start"] is True
        assert record["child_interval_is_the_trapped_region"] is False
        assert "outside the sealed region" in doc

    def test_the_saturation_number_matches(self, doc):
        """The saturation is the deepest edge of the six-edge sweep."""
        record = _load(EDGE)
        gains = record["sweep"]["amplitude_gain"]
        assert gains[-1] == pytest.approx(1.5849, abs=1e-3)
        assert record["sweep"]["gain_increases_with_depth"]
        assert "1.58" in doc

    def test_the_relaxation_numbers_match(self, doc):
        """The 4% figure is from the deep-edge long runs, not the release sweep."""
        record = _load(SETTLE)
        release_sweep_fractions = record["sweep"]["after_peak_fraction"]
        assert all(value > 0.4 for value in release_sweep_fractions)
        deep = [_load(DEVELOPMENT / "nsc-rn-inner-packet-shell-exc077-long-v1.json"),
                _load(DEVELOPMENT / "nsc-rn-inner-packet-shell-exc0760-long-v1.json")]
        for entry in deep:
            assert entry["accumulation_after_arrival"][
                "amplitude_after_peak_fraction"] < 0.05
        assert "relaxes to `4%`" in doc

    def test_the_retained_fraction_matches(self, doc):
        record = _load(EDGE)
        retained = record["sweep"]["fraction_below_inner_edge"]
        assert min(retained) > 0.93 and max(retained) < 0.95
        assert "94%" in doc

    def test_the_layer_widths_match(self, doc):
        record = _load(LAYER)
        widths = record["half_widths"]
        assert widths[0] / widths[-2] == pytest.approx(2.27, rel=0.02)
        assert "2.27" in doc


class TestTheLinksAreMarkedHonestly:
    def test_it_does_not_claim_the_cut_is_physical(self, doc):
        """The child/parent cut has no signature; that is the finding."""
        assert "no discontinuity" in doc
        assert "a label on a smooth metric" in doc

    def test_it_says_the_child_and_sealed_region_are_disjoint(self, doc):
        assert "disjoint" in doc
        assert "cannot say which" in doc

    def test_it_refutes_nest_to_nest(self, doc):
        assert "Refuted as posed for this construction" in doc
        assert "no second nest" in doc

    def test_it_separates_the_new_content_from_the_known_result(self, doc):
        assert "does not say what falls through" in doc

    def test_it_keeps_the_balance_result_as_consistency_not_discovery(self, doc):
        assert "criterion on a known result" in doc or \
            "balance-at-horizon" in doc


class TestTheLimitsAreStated:
    def test_it_is_not_a_scientific_result(self, doc):
        assert "not** a scientific result" in doc
        assert "does not measure the actual universe" in doc

    def test_it_names_what_survives(self, doc):
        assert "What survives" in doc
        assert "not the nested picture as described" in doc


class TestTheChainIsComplete:
    def test_every_link_is_marked(self, doc):
        for link in range(1, 8):
            assert f"Link {link} " in doc, link

    def test_each_link_is_marked_with_a_status(self, doc):
        assert doc.count("**Measured") >= 3
        assert "**Unanswerable" in doc
        assert "**Refuted" in doc