"""Tests for the mass-shell reduction of the RN radial family.

The interesting failure modes here are all about the estimator rather
than the physics: a median that is too short latches onto ringing, a
prominence test that is too loose calls a boundary value a moving
object, and a motion fit that runs through a phase change reads two
features as one trajectory. The planted-profile tests exist to catch
each of those directly.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons import nsc_rn_mass_shell as shell

LAB = Path(__file__).resolve().parents[1]
RECORD = LAB / "results/development/nsc-rn-mass-shell-v4.json"


def _payload(radius: np.ndarray, mass: np.ndarray) -> dict:
    return {"radius": radius, "coupled_mass": mass}


def _grid(n: int = 869, lo: float = 1.04, hi: float = 32.0) -> np.ndarray:
    return np.linspace(lo, hi, n)


def _payload_from_density(radius: np.ndarray, density: np.ndarray, base: float) -> dict:
    mass = base + np.concatenate([[0.0], np.cumsum(
        0.5 * (density[1:] + density[:-1]) * np.diff(radius))])
    return _payload(radius, mass)


class TestCompiledKernels:
    def test_kernels_are_compiled(self):
        from numba.core.registry import CPUDispatcher
        for kernel in (shell._mass_shell_density, shell._running_median,
                       shell._shell_peak, shell._shell_width):
            assert isinstance(kernel, CPUDispatcher), kernel

    def test_median_rejects_a_single_spike(self):
        values = np.ones(101)
        values[50] = 1e6
        smoothed = shell._running_median(values, 15)
        assert np.max(smoothed) < 1.5

    def test_median_of_flat_series_is_flat(self):
        assert np.allclose(shell._running_median(np.ones(51), 15), 1.0)


class TestShellEstimator:
    def test_recovers_a_planted_shell_radius(self):
        radius = _grid()
        density = np.exp(-0.5 * ((radius - 6.0) / 0.9) ** 2)
        state = shell.shell_state(_payload_from_density(radius, density, 1.04))
        assert state["shell_radius"] == pytest.approx(6.0, abs=0.15)
        assert state["off_inner_edge"]

    def test_recovers_a_planted_shell_width(self):
        radius = _grid()
        width = 1.5
        density = np.exp(-0.5 * ((radius - 7.0) / width) ** 2)
        state = shell.shell_state(_payload_from_density(radius, density, 1.04))
        expected = 2.0 * width * math.sqrt(2.0 * math.log(2.0))
        assert state["shell_width_half_max"] == pytest.approx(expected, rel=0.15)

    def test_flat_mass_has_no_shell(self):
        radius = _grid()
        state = shell.shell_state(_payload(radius, np.full(radius.shape, 1.04)))
        assert state["shell_amplitude"] == pytest.approx(0.0, abs=1e-18)
        assert not state["off_inner_edge"]

    def test_ringing_is_not_mistaken_for_a_shell(self):
        """A pure oscillation has no envelope peak to find."""
        radius = _grid()
        density = 1e-6 * np.cos(2.0 * np.pi * radius * 3.0)
        state = shell.shell_state(_payload_from_density(radius, density, 1.04))
        assert state["shell_radius"] > radius[0]

    def test_edge_pinned_shell_is_not_off_edge(self):
        radius = _grid()
        density = np.exp(-0.5 * ((radius - 1.04) / 0.3) ** 2)
        state = shell.shell_state(_payload_from_density(radius, density, 1.04))
        assert state["shell_radius"] == pytest.approx(radius[0], abs=1e-6)
        assert not state["off_inner_edge"]

    def test_rejects_a_descending_grid(self):
        radius = _grid()
        with pytest.raises(ValueError, match="strictly increasing"):
            shell.shell_state(_payload(radius[::-1].copy(),
                                       np.full(radius.shape, 1.04)))

    def test_rejects_a_mismatched_shape(self):
        radius = _grid()
        with pytest.raises(ValueError, match="share a shape"):
            shell.shell_state(_payload(radius, np.full(radius.shape[0] - 1, 1.04)))

    def test_rejects_a_grid_that_is_too_short(self):
        with pytest.raises(ValueError, match="at least eight"):
            shell.shell_state(_payload(np.linspace(1.0, 2.0, 5), np.ones(5)))


class TestPhaseLabelling:
    def _series(self):
        rows = []
        for t, radius, amplitude, off in [
            (0.0, 12.0, 1.0, True), (1.0, 10.5, 0.99, True), (2.0, 9.0, 0.98, True),
            (3.0, 7.0, 0.97, False), (4.0, 4.0, 0.05, True),
        ]:
            rows.append({"t_over_rm": t, "shell_radius": radius,
                         "shell_amplitude": amplitude, "off_inner_edge": off,
                         "r_inner": 1.04})
        return rows

    def test_phases_are_assigned_in_order(self):
        rows = self._series()
        shell.shell_motion(rows)
        assert [row["phase"] for row in rows[:3]] == ["falling"] * 3
        assert rows[3]["phase"] == "edge_pinned"
        assert rows[4]["phase"] == "late_weak_feature"

    def test_fit_stops_at_the_end_of_the_falling_phase(self):
        motion = shell.shell_motion(self._series())
        assert motion["stations_falling"] == 3
        assert motion["radius_last"] == pytest.approx(9.0)
        assert motion["inner_edge_arrival_t"] is None

    def test_weak_late_feature_is_not_counted_as_the_shell(self):
        motion = shell.shell_motion(self._series())
        assert motion["stations_late_weak_feature"] == 1
        assert motion["late_to_falling_amplitude_ratio"] == pytest.approx(0.05)

    def test_arrival_is_found_even_after_prominence_is_lost(self):
        rows = self._series()
        rows[3]["shell_radius"] = 1.04
        motion = shell.shell_motion(rows)
        assert motion["inner_edge_arrival_t"] == pytest.approx(3.0)

    def test_too_short_a_falling_phase_declines_the_fit(self):
        rows = [{"t_over_rm": 0.0, "shell_radius": 12.0, "shell_amplitude": 1.0,
                 "off_inner_edge": False, "r_inner": 1.04}]
        motion = shell.shell_motion(rows)
        assert not motion["fit_available"]
        assert "reason" in motion

    def test_linear_trajectory_is_recovered(self):
        rows = [{"t_over_rm": float(t), "shell_radius": 12.0 - 1.4 * t,
                 "shell_amplitude": 1.0, "off_inner_edge": True, "r_inner": 1.04}
                for t in range(6)]
        motion = shell.shell_motion(rows)
        assert motion["infall_rate_mean"] == pytest.approx(-1.4, rel=1e-9)
        assert motion["r_squared"] == pytest.approx(1.0, abs=1e-12)


class TestRecord:
    def test_record_exists_and_is_creation_only(self):
        record = json.loads(RECORD.read_text(encoding="utf-8"))
        assert record["creation_only"] is True
        assert record["api_version"] == shell.API_VERSION

    def test_record_covers_the_whole_family(self):
        record = json.loads(RECORD.read_text(encoding="utf-8"))
        assert record["coverage"]["cases"] == 6
        assert record["coverage"]["total_stations"] == 177

    def test_every_case_has_a_falling_phase(self):
        record = json.loads(RECORD.read_text(encoding="utf-8"))
        assert all(case["motion"]["fit_available"] for case in record["cases"])

    def test_infall_rate_is_nearly_parameter_independent(self):
        record = json.loads(RECORD.read_text(encoding="utf-8"))
        separation = record["separation"]
        assert separation["infall_rate_spread_relative"] < 0.01

    def test_source_scaling_is_linear(self):
        record = json.loads(RECORD.read_text(encoding="utf-8"))
        ratio = record["separation"]["source_scaling_ratio_between_largest_and_smallest_epsilon"]
        assert ratio == pytest.approx(10.0, rel=0.01)

    def test_record_states_it_is_not_the_universe(self):
        record = json.loads(RECORD.read_text(encoding="utf-8"))
        assert "does not measure the actual universe" in record["scope_note"]


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    root = LAB.parent
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_mass_shell.py"],
        cwd=root, capture_output=True, text=True)
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stderr