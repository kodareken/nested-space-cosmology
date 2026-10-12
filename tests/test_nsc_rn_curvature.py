"""Tests for the curvature reduction at the inner edge.

The question is whether mass compressing against the excision makes the
reconstructed curvature diverge there. The reduction's first version
answered with a statistic that moved for two reasons at once, so most of
these tests are about refusing to repeat that.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-curvature-v2.json"
RUN = DEVELOPMENT / "nsc-rn-inner-packet-curv"
DRIVER = LAB / "scripts/derive_nsc_rn_curvature.py"


@pytest.fixture(scope="module")
def curvature():
    """Named apart from the ``record`` method below, which would shadow it."""
    return json.loads(RECORD.read_text(encoding="utf-8"))


class TestTheCurvatureIsSourced:
    def test_every_station_carries_a_reading(self, curvature):
        assert curvature["stations_with_curvature"] == curvature["stations_total"] > 20

    def test_the_values_are_finite(self):
        files = sorted(RUN.glob("station-*.npz"))
        assert files
        for path in files[::3]:
            with np.load(path) as payload:
                assert "R4" in payload.files
                assert np.all(np.isfinite(payload["R4"]))
                assert np.all(np.isfinite(payload["Ricci2"]))

    def test_it_is_not_the_vacuum_reference(self, curvature):
        """A reference on the grid cannot know the shell arrived."""
        assert "does not know the shell is there" in curvature["what_this_measures"]


class TestTheEdgeBehaviour:
    def test_the_edge_value_is_reported_on_its_own(self, curvature):
        assert len(curvature["R4_edge_abs"]) == curvature["stations_with_curvature"]

    def test_the_maximum_is_found_and_located(self, curvature):
        assert curvature["R4_edge_max_abs"] == max(curvature["R4_edge_abs"])
        assert curvature["R4_edge_max_t"] is not None

    def test_the_maximum_is_bounded(self, curvature):
        assert curvature["R4_edge_max_over_initial"] < 10.0

    def test_it_reads_the_edge_at_the_mass_peak(self, curvature):
        """A curvature feature tracking the compression is the point."""
        assert curvature["R4_edge_at_mass_peak_t"] is not None
        assert curvature["R4_edge_at_mass_peak"] is not None

    def test_the_compression_peak_is_not_the_global_maximum(self, curvature):
        """The formation transient is larger; saying otherwise would be wrong."""
        assert curvature["R4_edge_max_t"] < curvature["R4_edge_at_mass_peak_t"]

    def test_the_edge_relaxes_below_where_it_started(self, curvature):
        assert curvature["R4_edge_final_over_initial"] < 1.0


class TestAgainstRepeatingTheEarlierMistake:
    def test_the_bulk_ratio_is_reported_with_its_caveat(self, curvature):
        assert "poor verdict on its own" in curvature["method"]["why_three_readings"]

    def test_the_neighbour_ratio_is_reported(self, curvature):
        """A local reading that does not move when the interior changes."""
        assert len(curvature["R4_edge_over_neighbour"]) == curvature["stations_with_curvature"]

    def test_the_verdict_is_not_an_endpoint_comparison(self, curvature):
        source = DRIVER.read_text(encoding="utf-8")
        assert "hides all three features" in curvature["reading"]
        assert "edge_overshoot" in source

    def test_the_first_verdict_is_retained_and_superseded(self, curvature):
        first = json.loads((DEVELOPMENT / "nsc-rn-curvature-v1.json").read_text(encoding="utf-8"))
        assert first["verdict"] == curvature["verdict"]
        # v1 reached the same label from a first-versus-last test that hid
        # the transient maximum. It must remain distinguishable.
        assert "R4_edge_max_abs" not in first
        assert "R4_edge_max_abs" in curvature


class TestScope:
    def test_it_says_it_is_not_the_universe(self, curvature):
        assert "does not measure the actual universe" in curvature["scope_note"]

    def test_it_does_not_call_the_excision_a_horizon(self, curvature):
        assert "not regularity at a horizon" in curvature["scope_note"]
        assert "horizon midpoint" in curvature["scope_note"]
        note = (LAB / "docs/nsc-rn-curvature.md").read_text(encoding="utf-8")
        assert "It is a boundary, not an event horizon" in note

    def test_it_names_the_excision_as_the_horizon_midpoint(self, curvature):
        assert "horizon midpoint" in curvature["scope_note"]

    def test_the_note_carries_the_same_limit(self):
        note = (LAB / "docs/nsc-rn-curvature.md").read_text(encoding="utf-8")
        assert "not regularity at a" in note.lower() or "is not regularity at a" in note
        assert "r_minus = 0.754" in note


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_curvature.py",
         "--run", "nsc-rn-inner-packet-curv",
         "--out", "nsc-rn-curvature-v2.json"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr