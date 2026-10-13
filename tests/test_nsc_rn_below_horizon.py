"""Tests for the region below the inner horizon.

The chart refuses to sample it, with the message "static modes stay strictly
above the inner horizon". These tests check that the refusal is a convention
in one function rather than a property of the geometry, that the model
evaluates below the inner horizon without trouble, and that the second zero
of the outgoing speed exists and sits on a horizon.

They also pin the decision not to relax the production guard, because that
is the thing most likely to be done carelessly next.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-below-horizon-v1.json"
NOTE = LAB / "docs/nsc-rn-geometry-map.md"
SOURCE = LAB / "src/recursive_horizons/nsc_rn_source.py"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestTheGuardIsAConvention:
    def test_the_guard_sits_beyond_the_real_domain(self, record):
        assert record["guard_is_beyond_the_real_domain_everywhere"]
        for entry in record["entries"]:
            assert entry["guard_is_beyond_the_real_domain"]

    def test_the_guard_is_not_at_the_geometry_limit(self, record):
        for entry in record["entries"]:
            assert entry["chart_guard_limit"] > entry["shift_real_domain_limit"]
            assert entry["guard_reaches_beyond_the_geometry"]

    def test_the_chart_really_does_refuse(self, record):
        assert "above the inner horizon" in record["chart_guard_message"]

    def test_the_guard_is_still_in_the_source(self):
        """It is not to be removed as a side effect of this finding."""
        text = SOURCE.read_text(encoding="utf-8")
        assert "static modes stay strictly above the inner horizon" in text


class TestTheModelEvaluatesThere:
    def test_it_does_for_every_entry(self, record):
        assert record["model_evaluates_below_inner_horizon_everywhere"]

    def test_every_probe_is_finite(self, record):
        for entry in record["entries"]:
            for probe in entry["model_probes"].values():
                assert probe["jets_ok"]
                assert probe["all_finite"]

    def test_the_light_constant_holds_down_there(self, record):
        for entry in record["entries"]:
            for probe in entry["model_probes"].values():
                assert probe["outgoing"] - probe["ingoing"] == pytest.approx(2.0)

    def test_the_lapse_is_still_one(self, record):
        assert record["lapse_is_one_everywhere"]
        for entry in record["entries"]:
            assert entry["lapse_is_one_throughout"]


class TestThereAreTwoZeroes:
    def test_two_for_every_entry(self, record):
        assert record["two_zeroes_everywhere"]

    def test_each_zero_sits_on_a_horizon(self, record):
        for entry in record["entries"]:
            for crossing in entry["outgoing_zeroes"]:
                on = (abs(crossing - entry["r_minus"]) < 1e-3
                      or abs(crossing - entry["r_plus"]) < 1e-3)
                assert on, (crossing, entry["r_minus"], entry["r_plus"])

    def test_outgoing_is_positive_again_below(self, record):
        for entry in record["entries"]:
            assert entry["positive_below_inner_horizon"]

    def test_this_is_why_every_reduction_reported_one(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "Only the outermost zero" in note
        assert "guarded domain" in note


class TestNothingWasRelaxed:
    def test_no_production_guard_was_changed(self, record):
        assert record["production_guard_changed"] is False

    def test_the_record_says_why_not(self, record):
        assert "wrong order of" in record["why_no_guard_was_relaxed"]

    def test_it_does_not_claim_the_third_region_contains_anything(self, record):
        assert "says nothing about whether it contains anything" in record["scope_note"]


class TestScope:
    def test_it_is_not_the_universe(self, record):
        flat = record["scope_note"].replace(chr(10), " ")
        assert "does not measure the actual universe" in flat

    def test_it_does_not_claim_the_nesting_iterates(self, record):
        assert "whether the nesting iterates" in record["scope_note"]

    def test_it_read_the_reference_model_only(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "no station payloads, no time step, no marching" in note.replace(
            chr(10), " ").lower()
