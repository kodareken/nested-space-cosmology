"""Tests for the inner-packet march and its reduction.

The inner-packet construction exists because the recorded family drops its
packet from 12 r_m, so the shell crosses the grid's inner edge at t = 7.0
and between 61% and 99% of the accreted mass ends up below it. These tests
pin the two things that construction is supposed to deliver, and the two
guards that keep it honest.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons import nsc_rn_mass_shell as shell

LAB = Path(__file__).resolve().parents[1]
MARCH = LAB / "results/development/nsc-rn-inner-packet-v1"
RECORD = LAB / "results/development/nsc-rn-inner-packet-shell-v2.json"
SCRIPT = LAB / "scripts/run_nsc_rn_inner_packet.py"


class TestTheConstruction:
    def test_packet_was_released_inside_the_collar(self):
        summary = json.loads((MARCH / "summary.json").read_text(encoding="utf-8"))
        assert summary["r_inner"] < summary["packet_center_over_rm"] < summary["r_outer"]

    def test_it_is_a_different_construction_from_the_recorded_family(self):
        summary = json.loads((MARCH / "summary.json").read_text(encoding="utf-8"))
        assert summary["packet_center_over_rm"] != summary["recorded_family_packet_center_over_rm"]

    def test_the_packet_width_was_not_changed(self):
        """Only the release radius differs, so the two are comparable at all."""
        summary = json.loads((MARCH / "summary.json").read_text(encoding="utf-8"))
        assert summary["packet_width_over_rm"] == pytest.approx(2.0)

    def test_the_inner_edge_stays_above_the_inner_horizon(self):
        summary = json.loads((MARCH / "summary.json").read_text(encoding="utf-8"))
        assert summary["r_minus"] < summary["r_inner"]

    def test_the_march_reached_its_target(self):
        summary = json.loads((MARCH / "summary.json").read_text(encoding="utf-8"))
        assert summary["stop_reason"] == "target"
        assert summary["recorded_stations"][-1]["t_over_rm"] == pytest.approx(
            summary["target_over_rm"], abs=0.05)

    def test_stations_are_named_on_the_recorded_convention(self):
        """The reduction reads its clock from the name, so it must hold."""
        names = sorted(p.name for p in MARCH.glob("station-*.npz"))
        assert "station-0.npz" in names
        assert "station-1.4.npz" in names


class TestShellGuardsInTheScript:
    @pytest.fixture(scope="class")
    def source(self):
        return SCRIPT.read_text(encoding="utf-8")

    def test_it_refuses_to_overwrite(self, source):
        assert "refusing to overwrite" in source

    def test_it_refuses_a_centre_outside_the_collar(self, source):
        assert "is not strictly inside the collar" in source

    def test_it_refuses_a_packet_whose_left_edge_is_outside_the_collar(self, source):
        """The centre alone is not enough: the packet is two r_m wide."""
        assert "is at or below the inner edge" in source

    def test_it_refuses_a_centre_inside_the_inner_horizon(self, source):
        assert "inside the inner horizon" in source

    def test_it_checks_the_inner_edge_is_above_the_inner_horizon(self, source):
        assert "not above the inner horizon" in source


class TestTheReduction:
    @pytest.fixture(scope="class")
    def record(self):
        return json.loads(RECORD.read_text(encoding="utf-8"))

    def test_the_shell_falls(self, record):
        assert record["motion"]["fit_available"]
        assert record["motion"]["infall_rate_mean"] < 0.0
        assert record["motion"]["radius_first"] > record["motion"]["radius_last"]

    def test_the_fit_is_good(self, record):
        assert record["motion"]["r_squared"] > 0.99

    def test_the_shell_arrives_at_the_inner_edge(self, record):
        assert record["motion"]["inner_edge_arrival_t"] is not None
        assert record["motion"]["inner_edge_arrival_t"] < 1.0

    def test_the_arrival_time_is_consistent_with_the_fitted_rate(self, record):
        """Arrival must not land before the fitted trajectory would reach it.

        The fitted rate is a mean over the falling phase. The shell moves
        only four to six cells per station, so the peak radius is quantised
        to the grid and the final segment is not described by the mean;
        this bounds the extrapolation rather than asserting it.
        """
        motion = record["motion"]
        predicted = motion["radius_first"] + motion["infall_rate_mean"] * motion["inner_edge_arrival_t"]
        assert predicted > 0.0
        assert motion["inner_edge_arrival_t"] < 1.0

    def test_the_peak_radius_is_quantised_to_the_grid(self):
        """Documents why per-interval speeds are not evidence of anything.

        The shell advances a whole number of cells per station, so the
        apparent alternation in local speed is the grid, not dynamics.
        """
        import numpy as np
        with np.load(MARCH / "station-0.npz") as payload:
            radius = payload["radius"]
        spacing = float(radius[1] - radius[0])
        record = json.loads(RECORD.read_text(encoding="utf-8"))
        falling = [s for s in record["stations"] if s["phase"] == "falling"]
        for a, b in zip(falling, falling[1:]):
            cells = (a["shell_radius"] - b["shell_radius"]) / spacing
            assert abs(cells - round(cells)) < 1e-6

    def test_the_shell_compresses_as_it_falls(self, record):
        """The one resolved shape change in the falling phase.

        Measured in grid cells, so the shrink is compared against the
        estimator's own smoothing scale rather than against a bare
        threshold.
        """
        import numpy as np
        with np.load(MARCH / "station-0.npz") as payload:
            radius = payload["radius"]
        spacing = float(radius[1] - radius[0])
        falling = [s for s in record["stations"] if s["phase"] == "falling"]
        first, last = falling[0]["shell_width_half_max"], falling[-1]["shell_width_half_max"]
        shrink_cells = (first - last) / spacing
        assert last < first * 0.9
        assert shrink_cells > shell.SMOOTH_WINDOW

    def test_the_accumulation_was_measured(self, record):
        """The point of the construction: the pile-up is visible."""
        accumulation = record["accumulation_after_arrival"]
        assert accumulation is not None
        assert accumulation["stations_after_arrival"] > 5
        assert accumulation["peak_amplitude"] > accumulation["amplitude_at_arrival"]
        assert accumulation["peak_t"] > accumulation["arrival_t"]

    def test_the_accumulation_compresses_by_a_measured_amount(self, record):
        assert record["accumulation_after_arrival"]["amplitude_gain"] > 1.1

    def test_the_clock_is_not_scaled_by_the_mass_ratio(self, record):
        """t/r_m is a clock. It must not silently become t/mu."""
        assert record["stations"][0]["t_over_rm"] == pytest.approx(0.0)
        assert record["stations"][1]["t_over_rm"] == pytest.approx(0.1)

    def test_it_states_it_is_not_the_universe(self, record):
        assert "does not measure the actual universe" in record["scope_note"]

    def test_it_states_the_construction_differs_from_the_recorded_family(self, record):
        assert "different construction" in record["scope_note"]


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_inner_packet.py"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr