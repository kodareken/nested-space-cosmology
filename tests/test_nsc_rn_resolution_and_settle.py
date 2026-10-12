"""Tests for the resolution control and the settled relaxation sweep.

The compression gain sits at the grid's left edge, where the one-sided SBP
closure lives, so it needed a resolution control before it could be claimed.
The relaxation needed a common clock, because the first sweep ended each run
at an arbitrary station and the resulting spread was measuring that rather
than the physics.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RESOLUTION = DEVELOPMENT / "nsc-rn-resolution-check-v1.json"
RELEASE = DEVELOPMENT / "nsc-rn-packet-sweep-v4.json"
SETTLED = DEVELOPMENT / "nsc-rn-packet-sweep-settle-v3.json"
MIXED = DEVELOPMENT / "nsc-rn-packet-sweep-v3.json"
SWEEP_DRIVER = LAB / "scripts/derive_nsc_rn_packet_sweep.py"


@pytest.fixture(scope="module")
def resolution():
    return json.loads(RESOLUTION.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def settled():
    return json.loads(SETTLED.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def release():
    return json.loads(RELEASE.read_text(encoding="utf-8"))


class TestResolutionControl:
    def test_it_actually_changes_the_resolution(self, resolution):
        points = resolution["points"]
        assert len(points) == 2
        assert points[1] > 1.9 * points[0]
        assert points[1] < 2.1 * points[0]

    def test_everything_else_is_held_fixed(self, resolution):
        held = resolution["held_fixed"]
        for key, values in held.items():
            assert len(set(values)) == 1, (key, values)

    def test_resolution_is_the_one_thing_allowed_to_differ(self, resolution):
        """The check must not require equal grids, or it proves nothing."""
        assert "points" not in resolution["held_fixed"]

    def test_both_runs_produced_a_compression(self, resolution):
        assert len(resolution["amplitude_gain"]) == 2
        assert all(gain > 1.1 for gain in resolution["amplitude_gain"])

    def test_the_gain_survives_the_resolution_change(self, resolution):
        assert resolution["amplitude_gain_relative_difference"] < 0.05
        assert resolution["verdict"] == "gain_holds_under_resolution_change"

    def test_the_arrival_time_is_resolution_independent(self, resolution):
        assert resolution["edge_arrival_t"][0] == pytest.approx(
            resolution["edge_arrival_t"][1], abs=0.05)

    def test_it_does_not_claim_a_convergence_order(self, resolution):
        assert "not a convergence order" in resolution["scope_note"]
        assert "not evidence of physics" in resolution["scope_note"]
        assert "Richardson" in resolution["what_this_does_not_establish"]

    def test_it_says_it_is_not_the_universe(self, resolution):
        assert "does not measure the actual universe" in resolution["scope_note"]


class TestTheSettledRelaxation:
    def test_the_settled_sweep_is_at_one_resolution(self, settled):
        assert len(set(settled["held_fixed"]["points"])) == 1

    def test_the_relaxation_spread_is_much_tighter_than_the_first_sweep(self, settled, release):
        first = release["sweep"]["after_peak_fraction_spread"]
        second = settled["sweep"]["after_peak_fraction_spread"]
        assert second < 0.25 * first

    def test_the_relaxation_agrees_to_within_a_tenth(self, settled):
        assert settled["sweep"]["after_peak_fraction_spread"] < 0.10

    def test_no_single_relaxed_value_is_claimed(self, settled):
        """Agreement at one lag bounds it; it does not settle it."""
        assert "bounds the relaxation" in settled["sweep"]["what_does_not_hold"]

    def test_the_compression_and_rise_time_still_hold(self, settled):
        assert settled["sweep"]["amplitude_gain_relative_spread"] < 0.05
        assert settled["sweep"]["rise_time_spread"] < 0.15

    def test_no_run_starts_pinned(self, settled):
        assert settled["sweep"]["runs_starting_pinned"] == []


class TestTheMixedResolutionRecord:
    def test_it_is_retained_and_self_evidently_invalid(self):
        """v3 was produced by a glob that mixed resolutions into a sweep."""
        mixed = json.loads(MIXED.read_text(encoding="utf-8"))
        assert len(set(mixed["held_fixed"]["points"])) > 1

    def test_the_driver_now_names_its_runs(self):
        source = SWEEP_DRIVER.read_text(encoding="utf-8")
        assert "--runs" in source
        assert "silently" not in source.split("--pattern")[0][:200]

    def test_the_driver_refuses_a_mixed_set(self):
        import subprocess
        import sys
        result = subprocess.run(
            [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_packet_sweep.py",
             "--runs", "nsc-rn-inner-packet-c2", "nsc-rn-inner-packet-v1",
             "nsc-rn-inner-packet-c3-fine", "--out", "_tmp-guard-check.json"],
            cwd=LAB.parent, capture_output=True, text=True)
        assert result.returncode == 2
        assert "not a release-radius sweep" in result.stderr


class TestScope:
    def test_the_release_sweep_still_holds(self, release):
        assert release["sweep"]["amplitude_gain_relative_spread"] < 0.05

    def test_the_first_sweep_still_records_that_relaxation_did_not_agree(self, release):
        assert release["sweep"]["after_peak_fraction_spread"] > 0.1
        assert "does not agree across runs" in release["sweep"]["what_does_not_hold"]

    def test_the_recorded_prose_is_data_driven(self):
        """Prose hardcoded from one run set misleads the next."""
        source = SWEEP_DRIVER.read_text(encoding="utf-8")
        assert "lines_hold" in source
        assert "No quantity in this record was found not to hold." in source


def test_sweep_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_packet_sweep.py",
         "--runs", "nsc-rn-inner-packet-v1", "nsc-rn-inner-packet-c2",
         "--out", "nsc-rn-packet-sweep-v4.json"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr