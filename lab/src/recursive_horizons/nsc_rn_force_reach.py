"""Effective reach of the two forces on the RN radial parent.

Douglas asked whether the larger force dominates a larger region and
reaches further, with nothing pushing back, and whether the opposing
force simply cancels away once other zones appear. This module turns
that into a finite measurement on the recorded areal grid.

Two quantities, both scale free and both computed per case:

* ``reach`` -- the areal radius out to which a force stays above a
  declared fraction of its own peak. The definition is the LAST crossing,
  not the first, so the ringing near the throat cannot make a force look
  like it died early. ``r`` counts as reached only while the force is
  still above the floor, so a later dip does not shorten the reach.
* ``dominant_fraction`` -- the fraction of grid nodes at which one force
  is strictly larger than the other.

A larger force with a larger reach supports the reading. A larger force
with a shorter reach contradicts it, and the record says which.

Scope. This is a construction check on recorded finite states. It does
not measure the actual universe, and a reach here is a property of one
computed collar on one grid, not a force law.

The oscillation in these profiles is characterised separately by
``nsc_rn_force_oscillation``; the last-crossing definition is chosen so
this module does not depend on that result.
"""
from __future__ import annotations

from typing import Any, Dict

import numpy as np
from numba import njit

API_VERSION = "nsc-rn-force-reach-v1"

#: A force counts as reaching a radius while it stays above this
#: fraction of its own peak. The same fraction as the balance wall's
#: live threshold, so both answers refer to the same supported region.
REACH_FLOOR = 1.0e-3


@njit(cache=True)
def _last_index_above(radius: np.ndarray, force: np.ndarray, floor: float) -> int:
    """Largest index whose magnitude is still above ``floor``.

    Scanning from the outside inward, the first node found above the
    floor is the last one, which is what makes the reach immune to the
    near-throat ringing.
    """
    for i in range(force.shape[0] - 1, -1, -1):
        if abs(force[i]) > floor:
            return i
    return -1


@njit(cache=True)
def _first_index_below(radius: np.ndarray, force: np.ndarray, floor: float) -> int:
    """First index at or inside which the force never returns above ``floor``.

    Returns the number of nodes in the supported prefix, so a force that
    dips below the floor and comes back is not counted as having died at
    the dip.
    """
    n = force.shape[0]
    last_below = -1
    for i in range(n):
        if abs(force[i]) <= floor:
            last_below = i
        elif last_below >= 0 and i > last_below:
            break
    return last_below + 1


@njit(cache=True)
def _dominant_fraction(inner: np.ndarray, outer: np.ndarray) -> float:
    """Fraction of nodes where the inner force is strictly larger."""
    if inner.shape[0] == 0:
        return 0.0
    count = 0
    for i in range(inner.shape[0]):
        if abs(inner[i]) > abs(outer[i]):
            count += 1
    return count / inner.shape[0]


@njit(cache=True)
def _crossover_radii(radius: np.ndarray, inner: np.ndarray, outer: np.ndarray) -> np.ndarray:
    """Radii where the ordering of the two magnitudes swaps, in order."""
    swaps = []
    for i in range(1, inner.shape[0]):
        before = abs(inner[i - 1]) > abs(outer[i - 1])
        after = abs(inner[i]) > abs(outer[i])
        if before != after:
            swaps.append(radius[i])
    return np.array(swaps, dtype=np.float64)


def force_reach(payload: Any, reach_floor: float = REACH_FLOOR) -> Dict[str, Any]:
    """Effective reach and dominance for one recorded station payload."""
    radius = np.ascontiguousarray(np.asarray(payload["radius"], dtype=np.float64))
    inner = np.ascontiguousarray(np.asarray(payload["radial_pressure"], dtype=np.float64))
    outer = np.ascontiguousarray(np.asarray(payload["angular_pressure"], dtype=np.float64))

    if radius.shape != inner.shape or radius.shape != outer.shape:
        raise ValueError("radius and both forces must share a shape")
    if radius.shape[0] < 3:
        raise ValueError("reach needs at least three grid nodes")
    if not np.all(np.diff(radius) > 0.0):
        raise ValueError("the areal grid must be strictly increasing")

    inner_peak = float(np.max(np.abs(inner)))
    outer_peak = float(np.max(np.abs(outer)))
    inner_floor = reach_floor * inner_peak
    outer_floor = reach_floor * outer_peak

    inner_last = _last_index_above(radius, inner, inner_floor)
    outer_last = _last_index_above(radius, outer, outer_floor)
    inner_prefix = _first_index_below(radius, inner, inner_floor)
    outer_prefix = _first_index_below(radius, outer, outer_floor)

    swappers = _crossover_radii(radius, inner, outer)
    inner_reach = float(radius[inner_last]) if inner_last >= 0 else float("nan")
    outer_reach = float(radius[outer_last]) if outer_last >= 0 else float("nan")
    larger = "radial" if inner_peak >= outer_peak else "angular"

    if np.isnan(inner_reach) or np.isnan(outer_reach):
        longer = None
    elif inner_reach > outer_reach:
        longer = "radial"
    elif outer_reach > inner_reach:
        longer = "angular"
    else:
        longer = "equal"

    return {
        "api_version": API_VERSION,
        "reach_floor": float(reach_floor),
        "points": int(radius.shape[0]),
        "r_inner": float(radius[0]),
        "r_outer": float(radius[-1]),
        "peak_inner": inner_peak,
        "peak_outer": outer_peak,
        "peak_ratio_inner_over_outer": inner_peak / (outer_peak or float("nan")),
        "larger_force": larger,
        "reach_inner": inner_reach,
        "reach_outer": outer_reach,
        "reach_longer": longer,
        "larger_force_reaches_further": bool(
            larger == longer and longer is not None),
        "supported_prefix_frac_inner": float(inner_prefix) / radius.shape[0],
        "supported_prefix_frac_outer": float(outer_prefix) / radius.shape[0],
        "dominant_fraction_inner": float(_dominant_fraction(inner, outer)),
        "crossover_count": int(swappers.size),
        "crossover_radii": [float(value) for value in swappers],
    }


def force_reach_verdict(profile: Dict[str, Any]) -> str:
    """One word for the reading, or ``mixed`` when the profile disagrees."""
    if profile["larger_force_reaches_further"]:
        return "larger_reaches_further"
    if profile["reach_longer"] == "equal":
        return "equal_reach"
    return "larger_reaches_less_far"