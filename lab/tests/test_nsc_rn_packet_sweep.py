"""Tests for the inner-packet release-radius sweep.

The sweep exists because one run cannot distinguish a construction from an
accident. These tests pin the two halves of that: the quantities that hold
across release radii, and the ones that do not. A test that only checked
the invariants would let a non-invariant creep into the note as a claim.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-packet-sweep-v2.json"
RUN_SCRIPT = LAB / "scripts/run_nsc_rn_inner_packet.py"


@pytest.fixture(scope="module")
def sweep():
    return json.loads(RECORD.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def sweep_block(sweep):
    return sweep["sweep"]


class TestTheSweepIsASweep:
    def test_it_actually_varies_the_release_radius(self, sweep_block):
        radii = sweep_block["release_radii_over_rm"]
        assert len(radii) >= 3
        assert len(set(radii)) == len(radii)
        assert radii == sorted(radii)

    def test_everything_else_is_held_fixed(self, sweep):
        held = sweep["held_fixed"]
        for key in ("mu", "eps", "packet_width_over_rm", "points", "r_inner"):
            assert len(held[key]) == 1, (key, held[key])

    def test_every_run_reached_its_target(self, sweep):
        assert all(run["stop_reason"] == "target" for run in sweep["runs"])


class TestWhatHolds:
    def test_the_compression_gain_holds_across_release_radii(self, sweep_block):
        """The headline: x1.197 was not an accident of one run."""
        assert sweep_block["amplitude_gain_relative_spread"] < 0.05

    def test_every_gain_is_a_real_compression(self, sweep_block):
        assert all(gain > 1.1 for gain in sweep_block["amplitude_gain"])

    def test_the_rise_time_holds(self, sweep_block):
        assert sweep_block["rise_time_spread"] < 0.15

    def test_every_rise_is_measured_after_arrival(self, sweep):
        for run in sweep["runs"]:
            accumulation = run["accumulation_after_arrival"]
            assert accumulation["peak_t"] > accumulation["arrival_t"]
            assert accumulation["rise_time"] > 0.0


class TestWhatDoesNotHold:
    def test_the_relaxation_is_not_settled_and_the_record_says_so(self, sweep_block):
        spread = sweep_block["after_peak_fraction_spread"]
        assert spread > 0.1, "if this ever tightens, update the note as well"
        assert "does not" in sweep_block["what_does_not_hold"]

    def test_the_width_compression_grows_with_release_radius(self, sweep):
        """It is a property of the fall, not a fixed factor."""
        moving = [r for r in sweep["runs"] if r["width_compression"]]
        assert len(moving) >= 2
        compressions = [(r["release_radius_over_rm"], r["width_compression"]) for r in moving]
        compressions.sort()
        assert compressions[-1][1] > compressions[0][1]

    def test_no_relaxed_value_is_claimed_anywhere(self, sweep):
        text = json.dumps(sweep)
        assert "no relaxed value is claimed" in text


class TestThePinnedRun:
    def test_the_short_release_is_flagged_as_starting_pinned(self, sweep_block):
        assert sweep_block["runs_starting_pinned"] == [2.0]

    def test_it_is_excluded_from_the_falling_phase_count(self, sweep_block):
        assert sweep_block["runs_with_a_falling_phase"] == sweep_block["runs"] - 1

    def test_it_contributes_no_infall_rate(self, sweep):
        pinned = [r for r in sweep["runs"] if r["starts_pinned_at_inner_edge"]]
        assert pinned
        assert all(not r["motion"]["fit_available"] for r in pinned)

    def test_it_still_contributes_a_compression_measurement(self, sweep):
        """It measures the accumulation, just not a fall."""
        pinned = [r for r in sweep["runs"] if r["starts_pinned_at_inner_edge"]]
        assert all(r["accumulation_after_arrival"]["amplitude_gain"] > 1.1 for r in pinned)


class TestTheGuardOnTheReleaseRadius:
    @pytest.fixture(scope="class")
    def source(self):
        return RUN_SCRIPT.read_text(encoding="utf-8")

    def test_it_checks_the_left_edge_not_just_the_centre(self, source):
        assert "left_edge" in source
        assert "left edge" in source

    def test_it_explains_why_the_centre_alone_is_not_enough(self, source):
        assert "binding condition, not the centre" in source

    def test_it_still_checks_the_inner_horizon(self, source):
        assert "inner horizon" in source


class TestScope:
    def test_it_says_it_is_not_the_universe(self, sweep):
        assert "does not measure the actual universe" in sweep["scope_note"]

    def test_holding_a_result_is_not_claim_of_physics(self, sweep):
        assert "does not make it physical" in sweep["scope_note"]


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_packet_sweep.py"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr