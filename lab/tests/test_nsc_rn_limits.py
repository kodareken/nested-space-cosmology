"""Tests for the instrument-limits document.

That document bounds everything else in this line, so it is the one whose
errors would matter most. These tests read the records each limit is drawn
from and check the document against them, and they pin the two questions it
declares unanswerable here, so they cannot quietly start being answered.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
DOC = LAB / "docs/nsc-rn-limits-of-the-instrument.md"
CHILD_CUT = DEVELOPMENT / "nsc-rn-child-cut-v2.json"
DOMAIN = DEVELOPMENT / "nsc-rn-causal-order-v4.json"
OUTWARD = DEVELOPMENT / "nsc-rn-outward-scaling-v2.json"
EDGE = DEVELOPMENT / "nsc-rn-edge-sweep-v2.json"
STATION = DEVELOPMENT / "nsc-rn-inner-packet-exc0760-long/station-0.npz"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def flat():
    return DOC.read_text(encoding="utf-8").replace(chr(10), " ")


class TestTheOneDegreeOfFreedomLimit:
    def test_the_recorded_state_is_radial(self):
        with np.load(STATION) as payload:
            assert payload["coupled_phi"].shape[2] == 1
            for key in ("coupled_mass", "radius", "frozen_lapse",
                        "frozen_shift", "R4", "Ricci2"):
                assert payload[key].ndim == 1

    def test_the_document_binds_it_correctly(self, flat):
        assert "every result here is radial structure" in flat
        assert "carries a sphere and supplies a radius" in flat


class TestTheChildCutLimit:
    def test_the_scan_found_nothing(self):
        assert _load(CHILD_CUT)["no_feature_everywhere"] is True

    def test_the_document_matches(self, flat):
        assert "Curvature excess runs `0.29` to `0.83`" in flat
        assert "asserts a nesting and supplies no physics" in flat


class TestTheTrappedRegionLimit:
    def test_every_measurement_was_taken_there(self, flat):
        record = _load(DOMAIN)
        assert record["all_one_way"] is True
        assert "taken there" in flat

    def test_it_declares_the_third_region_unvisited_not_unreachable(self, flat):
        assert "not merely unvisited" in flat
        assert "new construction rather than" in flat


class TestTheNoExpansionLimit:
    def test_the_outward_exponent_is_minus_one_half(self):
        record = _load(OUTWARD)
        assert record["exponents"]["1.04"] == pytest.approx(-0.5, abs=2e-3)

    def test_the_document_does_not_claim_a_cosmological_result(self, flat):
        assert "speaks to dark energy" in flat
        assert "would be the easy error" in flat


class TestWhatIsLeftIsNotAMeasurement:
    def test_it_declares_the_occupancy_question_open(self, flat):
        assert "not an occupancy" in flat

    def test_it_declares_the_intent_question_open(self, flat):
        assert "question about the construction, not about the code" in flat

    def test_the_boundary_of_the_instrument_is_stated(self, flat):
        assert "Bounded:" in flat

    def test_it_says_what_the_instrument_is_good_for(self, flat):
        assert "not what 2005 answers" in flat


class TestScope:
    def test_it_is_not_adding_a_result(self, flat):
        assert "Its purpose is to bound the other results" in flat

    def test_it_says_it_is_not_the_universe(self, flat):
        assert "does not measure the actual universe" in flat

    def test_it_does_not_criticise_the_results(self, flat):
        assert "not a criticism of the results" in flat