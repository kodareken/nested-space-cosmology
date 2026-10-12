"""Tests for the late-time behaviour of the compressed shell.

The question is whether the shell decays to nothing or leaves a remainder.
The answer decides whether the shell empties or persists, and it is a
measurement rather than a story.

These tests pin the three things that answer depends on: the decay turns
over, the rise after the trough is monotone rather than jitter, and the
amplitude never reaches zero at any recorded station.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-late-time-v1.json"
NOTE = LAB / "docs/nsc-rn-late-time.md"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestThereIsARemainder:
    def test_it_never_reaches_zero(self, record):
        assert record["never_reaches_zero_everywhere"]

    def test_the_decay_terminates_on_a_floor(self, record):
        for entry in record["entries"]:
            assert entry["trough_after_peak"]
            assert entry["trough_over_peak"] < 0.05

    def test_the_amplitude_grows_again_after_the_trough(self, record):
        assert record["turns_over_everywhere"]
        for entry in record["entries"]:
            assert entry["rises_after_the_trough"]
            assert entry["final_over_peak"] > entry["trough_over_peak"]

    def test_the_rise_is_monotone_not_jitter(self, record):
        for entry in record["entries"]:
            assert entry["rise_is_monotone"]
            assert entry["rise_steps"] >= 8


class TestTheFloorsAgree:
    def test_the_two_edges_give_the_same_floor(self, record):
        floors = [entry["trough_over_peak"] for entry in record["entries"]]
        assert abs(floors[0] - floors[1]) < 0.005

    def test_the_trough_is_at_the_same_time(self, record):
        times = {entry["trough_time"] for entry in record["entries"]}
        assert len(times) == 1

    def test_the_decay_is_more_than_a_factor_of_ten(self, record):
        for entry in record["entries"]:
            assert entry["decay_factor_peak_to_trough"] < 0.1
            assert entry["decay_steps"] >= 10


class TestWhatIsNotClaimed:
    def test_it_does_not_say_the_state_has_settled(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "long-time state is not known" in flat
        assert "still monotone at the last recorded station" in flat

    def test_it_does_not_say_the_remainder_is_physical(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "nothing here says it is physical" in flat

    def test_it_does_not_explain_the_horizon_coincidence(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "recorded and not explained" in flat

    def test_it_says_it_is_not_the_universe(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not measure the actual universe" in flat
