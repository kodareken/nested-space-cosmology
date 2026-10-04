"""Tests for the long deep-edge runs.

These settle two things the edge sweep left open: whether the compression
gain survived a longer run, and whether the relaxation completes. They also
pin a bug the longer runs exposed, where the curvature reduction read the
mass peak time off a different march.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
SHELL_077 = DEVELOPMENT / "nsc-rn-inner-packet-shell-exc077-long-v1.json"
SHELL_080 = DEVELOPMENT / "nsc-rn-inner-packet-shell-exc080-long-v1.json"
SHORT_077 = DEVELOPMENT / "nsc-rn-inner-packet-shell-exc085-v1.json"
CURV_077 = DEVELOPMENT / "nsc-rn-curvature-exc077-long-v2.json"
CURV_080 = DEVELOPMENT / "nsc-rn-curvature-exc080-long-v2.json"
CURV_077_BAD = DEVELOPMENT / "nsc-rn-curvature-exc077-long-v1.json"
EDGE_SWEEP = DEVELOPMENT / "nsc-rn-edge-sweep-v1.json"
NOTE = LAB / "docs/nsc-rn-deep-long.md"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def deep():
    return _load(SHELL_077)


@pytest.fixture(scope="module")
def slightly():
    return _load(SHELL_080)


class TestTheGainSurvivedTheLongerRun:
    def test_the_gain_is_unchanged_from_the_shorter_run(self, deep):
        sweep = _load(EDGE_SWEEP)
        before = dict(zip(sweep["sweep"]["r_inner"], sweep["sweep"]["amplitude_gain"]))
        assert deep["accumulation_after_arrival"]["amplitude_gain"] == pytest.approx(
            before[0.77], rel=1e-3)

    def test_the_peak_is_inside_the_run(self, deep, slightly):
        for record in (deep, slightly):
            peak = record["accumulation_after_arrival"]["peak_t"]
            assert 0.5 < peak < record["stations"][-1]["t_over_rm"] - 1.0


class TestTheRelaxationCompletes:
    def test_it_decays_to_a_few_percent(self, deep, slightly):
        for record in (deep, slightly):
            assert record["accumulation_after_arrival"]["amplitude_after_peak_fraction"] < 0.06

    def test_the_two_edges_agree(self, deep, slightly):
        a = deep["accumulation_after_arrival"]["amplitude_after_peak_fraction"]
        b = slightly["accumulation_after_arrival"]["amplitude_after_peak_fraction"]
        assert abs(a - b) < 0.01

    def test_the_deep_domain_retains_the_mass(self, deep, slightly):
        for record in (deep, slightly):
            assert record["mass_budget"]["fraction_below_inner_edge"] > 0.90

    def test_the_note_supersedes_the_too_short_claim(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "too short to measure a" in note
        assert "is now superseded" in note


class TestTheCurvatureLags:
    def test_the_global_maximum_is_the_formation_transient(self):
        for path in (CURV_077, CURV_080):
            record = _load(path)
            assert record["R4_edge_max_t"] == pytest.approx(0.1)

    def test_it_reads_this_runs_own_mass_peak(self):
        for path in (CURV_077, CURV_080):
            record = _load(path)
            assert record["mass_peak_t"] == pytest.approx(1.7)
            assert record["R4_edge_at_mass_peak_t"] == pytest.approx(1.7)

    def test_the_edge_curvature_at_the_mass_peak_matches_across_edges(self):
        a = _load(CURV_077)["R4_edge_at_mass_peak"]
        b = _load(CURV_080)["R4_edge_at_mass_peak"]
        assert abs(a - b) / a < 0.02

    def test_nothing_diverges(self):
        for path in (CURV_077, CURV_080):
            assert _load(path)["R4_edge_max_over_initial"] < 10.0


class TestTheWrongClockBug:
    def test_the_superseded_record_still_carries_the_wrong_column(self):
        bad = _load(CURV_077_BAD)
        assert "mass_peak_t" not in bad
        assert bad["R4_edge_at_mass_peak_t"] == pytest.approx(1.4)

    def test_the_wrong_column_is_off_by_a_lot(self):
        """It was a factor of five, so it could not have been subtle."""
        bad = _load(CURV_077_BAD)["R4_edge_at_mass_peak"]
        good = _load(CURV_077)["R4_edge_at_mass_peak"]
        assert good / bad > 3.0

    def test_the_driver_no_longer_hardcodes_a_run(self):
        source = (LAB / "scripts/derive_nsc_rn_curvature.py").read_text(encoding="utf-8")
        assert 'nsc-rn-inner-packet-shell-v2.json' not in source
        assert "--shell-record" in source

    def test_the_note_records_it(self):
        assert "hardcoded" in NOTE.read_text(encoding="utf-8")


class TestScope:
    def test_every_edge_is_still_a_boundary(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "still a **boundary**" in note
        assert "0.0157" in note

    def test_it_says_it_is_not_the_universe(self):
        assert "does not measure" in NOTE.read_text(encoding="utf-8")

    def test_it_does_not_claim_the_trend_turned_over(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "still rising with depth" in note or "has not turned over" in note

    def test_the_lag_is_marked_as_needing_its_own_sweep(self):
        assert "deserves its own" in NOTE.read_text(encoding="utf-8")


class TestStationCounts:
    @pytest.mark.parametrize("run", ["exc077-long", "exc080-long"])
    def test_the_long_runs_are_long(self, run):
        summary = _load(DEVELOPMENT / f"nsc-rn-inner-packet-{run}" / "summary.json")
        assert summary["target_over_rm"] == 6.0
        assert len(summary["recorded_stations"]) == 61
        assert summary["stop_reason"] == "target"