"""Tests for the outward scaling.

The claim tested is that the outward direction is far more sensitive than
the inward one, exponentially so. The measurement is the local exponent of
the outgoing speed's approach to its asymptote: an exponential would show
that exponent growing without bound, a power law shows it settling.

These tests pin the negative result and, importantly, pin what it is a
negative result *about*. This chart has no expansion, so it cannot test a
cosmological claim, and the record must not be read as if it could.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-outward-scaling-v2.json"
NOTE = LAB / "docs/nsc-rn-outward-scaling.md"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestTheOutwardRegionIsAlgebraic:
    def test_the_exponent_settles_on_minus_one_half(self, record):
        assert record["all_outward_algebraic"]
        for value in record["exponents"].values():
            assert value == pytest.approx(-0.5, abs=2e-3)

    def test_it_converges_rather_than_growing(self, record):
        """An exponential would diverge; a power law approaches a constant.

        The exponent here moves toward minus one half and stops, so the
        distance to that value shrinks at every step rather than growing.
        """
        for entry in record["entries"]:
            slopes = [p["local_exponent"] for p in entry["outward_local_exponents"]]
            distances = [abs(slope + 0.5) for slope in slopes]
            assert all(b < a for a, b in zip(distances, distances[1:]))
            assert entry["outward_exponent_converged"]

    def test_it_is_not_exponential_anywhere(self, record):
        for entry in record["entries"]:
            assert entry["outward_is_exponential"] is False

    def test_the_record_says_what_exponential_would_have_looked_like(self, record):
        flat = record["entries"][0]["exponent_test"].replace(chr(10), " ")
        assert "the local exponent growing without bound" in flat

    def test_the_measure_is_made_inward_too(self, record):
        for entry in record["entries"]:
            assert len(entry["inward_local_exponents"]) >= 3


class TestTheResultIsBoundedProperly:
    def test_it_does_not_claim_anything_cosmological(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "does not test any claim about cosmological expansion" in flat
        assert "this chart has none" in flat

    def test_it_says_a_negative_result_is_about_this_geometry_only(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "about this geometry only" in flat

    def test_it_says_it_is_not_the_universe(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "does not measure the actual universe" in flat

    def test_the_note_carries_the_same_bound(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "no expansion" in flat
        assert "does not speak to dark energy" in flat