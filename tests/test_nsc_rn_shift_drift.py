"""Tests for the shift-drift reduction.

The claim is that everything causal in this metric reduces to one
radius-varying quantity: the lapse is pinned, the two characteristic speeds
sum to twice it, and the shift carries all of the spatial variation. These
tests check the identity exactly, check that the reference chart and the
evolved metric agree on where the outgoing speed vanishes, and pin the one
place the shift is not monotone.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-shift-drift-v2.json"
V1 = DEVELOPMENT / "nsc-rn-shift-drift-v1.json"
DRIVER = LAB / "scripts/derive_nsc_rn_shift_drift.py"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestTheLightConstant:
    def test_the_identity_holds_exactly(self, record):
        assert record["light_constant"]["holds"]

    def test_the_lapse_is_pinned_to_one(self, record):
        block = record["light_constant"]
        assert block["lapse_min"] == pytest.approx(1.0, abs=1e-9)
        assert block["lapse_max"] == pytest.approx(1.0, abs=1e-9)

    def test_the_speed_sum_is_pinned_to_two(self, record):
        block = record["light_constant"]
        assert block["constant_min"] == pytest.approx(2.0, abs=1e-9)
        assert block["constant_max"] == pytest.approx(2.0, abs=1e-9)

    def test_it_is_sampled_finely(self, record):
        assert record["chart"]["samples"] >= 1000


class TestTheShiftCarriesEverything:
    def test_the_shift_is_the_only_radius_varying_quantity(self, record):
        text = record["shift_drift"]["beta_is_the_only_radius_varying_quantity"]
        assert "carries all of the spatial variation" in text

    def test_the_shift_falls_across_the_chart(self, record):
        drift = record["shift_drift"]
        assert drift["beta_at_inner"] > drift["beta_at_outer"]

    def test_it_is_not_monotone_and_that_is_explained(self, record):
        """The turn is a property of the charge, not of the dynamics."""
        drift = record["shift_drift"]
        assert drift["beta_monotone_decreasing"] is False
        assert drift["beta_max"] > drift["beta_at_inner"]
        assert drift["beta_turn_is_above_default_edge"] is False

    def test_the_turn_is_where_the_radicand_extremises(self, record):
        """beta = sqrt(2M/r - Q^2/r^2) turns at r = Q^2/M."""
        drift = record["shift_drift"]
        mu = record["chart"]["mu"]
        assert drift["beta_max_radius"] == pytest.approx(1.0 / mu, abs=1e-3)

    def test_the_turn_lies_below_the_default_edge(self, record):
        drift = record["shift_drift"]
        assert drift["beta_max_radius"] < 0.5 * (record["chart"]["r_minus"]
                                                 + record["chart"]["r_plus"])


class TestTheTwoHorizonsAgree:
    def test_the_reference_zero_is_at_r_plus(self, record):
        assert record["reference_zero_of_outgoing"] == pytest.approx(
            record["chart"]["r_plus"], abs=1e-5)

    def test_the_evolved_horizon_was_located(self, record):
        assert record["evolved_trapping_horizon"] is not None

    def test_they_agree_within_a_grid_cell(self, record):
        assert record["reference_and_evolved_agree"] is True
        assert record["agreement_in_cells"] < 1.0

    def test_it_is_described_as_a_check_not_a_derivation(self, record):
        assert "not a derivation" in record["scope_note"]
        assert "consistency check" in record["scope_note"]


class TestScope:
    def test_it_says_it_is_not_the_universe(self, record):
        assert "does not measure the actual universe" in record["scope_note"]

    def test_it_says_the_identity_is_family_specific(self, record):
        assert "says nothing about any other" in record["scope_note"]

    def test_the_first_record_is_retained(self):
        first = _load(V1)
        assert "beta_max_radius" not in first["shift_drift"]


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_shift_drift.py"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr