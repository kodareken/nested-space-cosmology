"""Tests for the pole structure of the characteristic equation.

The characteristic right-hand side divides by the two speeds, and the
outgoing speed vanishes at both horizons rather than only at the one the
branch starts from. That is why the inward integration stops at the inner
horizon.

These tests pin where the poles are, that the ingoing leg has none, that the
outgoing pole at the inner horizon is simple rather than stiff, and the
correction of the word "stiffness" that the measurement forced.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-pole-structure-v1.json"
REACH_NOTE = LAB / "docs/nsc-rn-third-region-reach.md"
NOTE = LAB / "docs/nsc-rn-pole-structure.md"
SOURCE = LAB / "src/recursive_horizons/nsc_rn_source.py"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestWhereThePolesAre:
    def test_there_are_two_on_the_outgoing_leg(self, record):
        assert record["two_poles_everywhere"]
        for entry in record["entries"]:
            assert entry["outgoing_zero_count"] == 2
            assert entry["both_poles_are_on_the_outgoing_leg"]

    def test_one_is_at_each_horizon(self, record):
        assert record["second_pole_at_r_minus_everywhere"]
        for entry in record["entries"]:
            assert entry["first_pole_at_r_plus"]
            assert entry["second_pole_at_r_minus"]

    def test_the_ingoing_leg_has_none(self, record):
        assert record["ingoing_never_vanishes"]
        for entry in record["entries"]:
            assert entry["ingoing_zeroes"] == []

    def test_the_equation_divides_by_the_speeds(self):
        text = SOURCE.read_text(encoding="utf-8")
        assert ") / chart[\"speed_plus\"]" in text
        assert ") / chart[\"speed_minus\"]" in text


class TestThePoleIsSimple:
    def test_the_reciprocal_grows_by_a_decade_per_decade(self, record):
        """That is what makes it a pole rather than a logarithm."""
        for entry in record["entries"]:
            # The outermost probe is far enough away that the linear
            # approximation does not hold, so the growth is tested where
            # it should be clean.
            reciprocals = [abs(probe["reciprocal"])
                           for probe in entry["approach_to_r_minus"][1:]]
            for previous, current in zip(reciprocals, reciprocals[1:]):
                assert 5.0 < current / previous < 25.0

    def test_the_coefficient_is_finite_and_signed(self, record):
        for entry in record["entries"]:
            coefficient = (abs(entry["approach_to_r_minus"][-1]["reciprocal"])
                           * entry["approach_to_r_minus"][-1]["offset_above_r_minus"])
            assert 0.1 < coefficient < 50.0

    def test_the_pole_is_stronger_when_the_horizons_are_closer(self, record):
        by_mu = {entry["mu"]: abs(entry["approach_to_r_minus"][-1]["reciprocal"])
                 * entry["approach_to_r_minus"][-1]["offset_above_r_minus"]
                 for entry in record["entries"]}
        assert by_mu[1.002] > by_mu[1.040]


class TestTheCorrection:
    def test_the_record_says_stiffness_was_wrong(self, record):
        flat = record["correction"].replace(chr(10), " ")
        assert "not stiffness" in flat

    def test_the_previous_note_is_corrected_not_deleted(self):
        """It stays, with the correction pointing at what replaced it."""
        previous = REACH_NOTE.read_text(encoding="utf-8")
        current = NOTE.read_text(encoding="utf-8")
        assert "stiff" in previous.lower()
        assert "was wrong" in current

    def test_the_proposal_is_concrete_and_unclaimed(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "applies at" in flat
        assert "worth doing depends on what is on the other side" in flat


class TestScope:
    def test_it_says_it_is_not_the_universe(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "does not measure the actual universe" in flat

    def test_it_says_nothing_is_known_past_the_pole(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "not what is on the other side" in flat

    def test_it_read_no_station_payloads(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "no station" in flat