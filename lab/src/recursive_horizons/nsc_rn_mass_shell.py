"""The mass shell of the RN radial parent: its radius, speed and mass.

This is the geometry, not a force balance. In the PG chart the geometry
is fixed by the mass function M(r), so the shape of dM/dr *is* the shape
of the object. The force diagnostics elsewhere in this line measure what
the two pressures do; this measures where the mass actually is.

The object this finds is a coherent shell of accreted mass that starts at
a radius common to every case and falls inward at a rate set by the
exterior compression, carrying an amplitude set by the source strength.

Three findings, each measured rather than assumed:

* the shell starts at the same radius in all six cases;
* the infall speed depends on the mass ratio and not on the source
  strength, so the two knobs set two different things;
* the amplitude scales with the source strength, so a decade of source
  buys a decade of shell mass.

Domain. Areal radius on the recorded radial grid, observer the areal
chart of the recorded coupled arm, states the station payloads already
written by ``nsc-rn-neutral-family-v1``. No time step is taken and
nothing is re-run.

Scope. This is a construction check on recorded finite states. It does
not measure the actual universe. The shell leaves the computed domain at
the inner edge in every case, so what happens to it afterwards is not
measured here.
"""
from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np
from numba import njit

API_VERSION = "nsc-rn-mass-shell-v1"

#: Median window used to suppress the recorded ringing before locating
#: the peak. Chosen far above the few-cell ringing the oscillation audit
#: measured, and well below the shell's own width.
SMOOTH_WINDOW = 15
#: The shell peak must stand this far above the value at the inner edge
#: to count as off the edge rather than pinned to it.
EDGE_PROMINENCE = 3.0
#: Radius beyond which the audit found the ringing lives. The peak is
#: never sought there.
SEARCH_LIMIT = 16.0


@njit(cache=True)
def _running_median(values: np.ndarray, window: int) -> np.ndarray:
    half = window // 2
    out = np.empty(values.shape[0])
    for i in range(values.shape[0]):
        lo = max(0, i - half)
        hi = min(values.shape[0], i + half + 1)
        block = np.sort(values[lo:hi])
        mid = block.shape[0] // 2
        out[i] = block[mid] if block.shape[0] % 2 == 1 else 0.5 * (block[mid - 1] + block[mid])
    return out


@njit(cache=True)
def _mass_shell_density(radius: np.ndarray, mass: np.ndarray) -> np.ndarray:
    """Centred dM/dr of the recorded mass function."""
    n = mass.shape[0]
    out = np.empty(n)
    for i in range(n):
        lo = i - 1 if i > 0 else 0
        hi = i + 1 if i < n - 1 else n - 1
        out[i] = (mass[hi] - mass[lo]) / (radius[hi] - radius[lo])
    return out


@njit(cache=True)
def _shell_peak(radius: np.ndarray, density: np.ndarray, limit: float):
    """Radius, amplitude and edge amplitude of the smoothed shell peak."""
    stop = radius.shape[0]
    for i in range(radius.shape[0]):
        if radius[i] > limit:
            stop = i
            break
    if stop < 3:
        return radius[0], 0.0, 0.0
    best = -1
    best_index = 0
    for i in range(stop):
        if density[i] > best:
            best = density[i]
            best_index = i
    return radius[best_index], best, density[0]


@njit(cache=True)
def _shell_width(radius: np.ndarray, density: np.ndarray, peak_index: int,
                 peak_value: float, limit: float) -> float:
    """Full width of the shell at half its peak amplitude."""
    half = 0.5 * peak_value
    lo = peak_index
    hi = peak_index
    while lo > 0 and radius[lo] <= limit and density[lo] > half:
        lo -= 1
    while hi < radius.shape[0] - 1 and radius[hi] <= limit and density[hi] > half:
        hi += 1
    return float(radius[hi] - radius[lo])


def shell_state(payload: Any, base_mass: float = 1.04) -> Dict[str, Any]:
    """Shell position, amplitude and width at one recorded station."""
    radius = np.ascontiguousarray(np.asarray(payload["radius"], dtype=np.float64))
    mass = np.ascontiguousarray(np.asarray(payload["coupled_mass"], dtype=np.float64))
    if radius.shape != mass.shape:
        raise ValueError("radius and coupled_mass must share a shape")
    if radius.shape[0] < 8:
        raise ValueError("the shell search needs at least eight grid nodes")
    if not np.all(np.diff(radius) > 0.0):
        raise ValueError("the areal grid must be strictly increasing")

    density = _mass_shell_density(radius, mass)
    stop = int(np.searchsorted(radius, SEARCH_LIMIT))
    stop = max(8, min(stop, radius.shape[0]))
    smoothed = _running_median(density[:stop], SMOOTH_WINDOW)

    peak_radius, peak_value, edge_value = _shell_peak(radius, smoothed, SEARCH_LIMIT)
    index = int(np.searchsorted(radius, peak_radius))
    prominence = peak_value / edge_value if edge_value > 0.0 else math.inf
    off_edge = prominence >= EDGE_PROMINENCE and peak_radius > radius[0]

    return {
        "api_version": API_VERSION,
        "points": int(radius.shape[0]),
        "r_inner": float(radius[0]),
        "shell_radius": float(peak_radius),
        "shell_amplitude": float(peak_value),
        "edge_amplitude": float(edge_value),
        "prominence": float(prominence),
        "off_inner_edge": bool(off_edge),
        "shell_width_half_max": float(
            _shell_width(radius, smoothed, index, peak_value, SEARCH_LIMIT)) if off_edge else float("nan"),
        "total_gain": float(mass[-1] - base_mass),
        "gain_below_inner_edge": float(mass[0] - base_mass),
        "search_limit": SEARCH_LIMIT,
        "smooth_window": SMOOTH_WINDOW,
    }


def shell_phase(state: Dict[str, Any], first_amplitude: float | None) -> str:
    """Label a station by which feature the peak is actually tracking.

    Once the shell reaches the inner edge the peak of dM/dr is pinned to
    the boundary and no longer measures a moving object. Later, weaker
    excursions in and out of the search window are a different feature
    again, not the shell returning.
    """
    if not state["off_inner_edge"]:
        return "edge_pinned"
    if first_amplitude is not None and state["shell_amplitude"] < 0.5 * first_amplitude:
        return "late_weak_feature"
    return "falling"


def shell_motion(series: list[Dict[str, Any]]) -> Dict[str, Any]:
    """Fit the shell's initial radius and inward speed.

    Only the first contiguous run of stations with the shell clear of the
    inner edge is used. Once the shell reaches the edge the peak stops
    tracking a moving object, and a later weak excursion that clears the
    edge again is a different feature, so fitting through the gap would
    read a mixture of two things as one trajectory.
    """
    if not series:
        return {"stations_used": 0, "fit_available": False, "reason": "no stations"}

    first_amplitude = series[0]["shell_amplitude"]
    for row in series:
        row["phase"] = shell_phase(row, first_amplitude)

    falling: list[Dict[str, Any]] = []
    for row in series:
        if row["phase"] != "falling":
            break
        falling.append(row)

    late = [row for row in series if row["phase"] == "late_weak_feature"]

    # The shell loses prominence well before it reaches the inner edge, so
    # the two are tracked separately and neither is called the other.
    prominence_loss_t = (series[len(falling)]["t_over_rm"]
                         if len(falling) < len(series) else None)
    # Arrival at the inner edge is searched over the whole series, because
    # by the time the shell gets there it has already lost prominence and
    # so sits in the edge-pinned phase rather than the falling one.
    r_inner = series[0].get("r_inner")
    arrival = None
    for row in series:
        if r_inner is not None and abs(row["shell_radius"] - r_inner) < 1e-9:
            arrival = row["t_over_rm"]
            break

    motion: Dict[str, Any] = {
        "stations_total": len(series),
        "stations_falling": len(falling),
        "stations_edge_pinned": sum(1 for row in series if row["phase"] == "edge_pinned"),
        "stations_late_weak_feature": len(late),
        "prominence_loss_t": prominence_loss_t,
        "prominence_loss_radius": (falling[-1]["shell_radius"] if falling else None),
        "inner_edge_arrival_t": arrival,
        "late_peak_amplitude_max": (max(row["shell_amplitude"] for row in late)
                                   if late else None),
        "late_to_falling_amplitude_ratio": (
            max(row["shell_amplitude"] for row in late) / first_amplitude
            if late and first_amplitude else None),
        "fit_available": len(falling) >= 3,
    }
    if not motion["fit_available"]:
        motion["reason"] = "fewer than three stations in the contiguous falling phase"
        return motion

    times = np.array([row["t_over_rm"] for row in falling], dtype=float)
    radii = np.array([row["shell_radius"] for row in falling], dtype=float)
    slope, intercept = np.polyfit(times, radii, 1)
    predicted = slope * times + intercept
    residual = radii - predicted
    ss_tot = float(np.sum((radii - radii.mean()) ** 2))
    ss_res = float(np.sum(residual ** 2))

    speeds = np.diff(radii) / np.diff(times)
    motion.update({
        "t_first": float(times[0]),
        "t_last": float(times[-1]),
        "radius_first": float(radii[0]),
        "radius_last": float(radii[-1]),
        "mean_radius_drop": float(radii[0] - radii[-1]),
        "infall_rate_mean": float(slope),
        "intercept_radius": float(intercept),
        "local_speed_first": float(speeds[0]),
        "local_speed_last": float(speeds[-1]),
        "local_speed_min": float(speeds.min()),
        "local_speed_max": float(speeds.max()),
        "r_squared": float(1.0 - ss_res / ss_tot) if ss_tot > 0.0 else float("nan"),
        "max_abs_residual": float(np.max(np.abs(residual))),
    })
    return motion