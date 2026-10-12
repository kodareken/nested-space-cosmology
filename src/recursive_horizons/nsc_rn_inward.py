"""The inward-march decision helper.

The recorded family stops at the horizon midpoint, outside the RN inner
horizon, so the Hayward regularity criterion (gr-qc/0506126: a static
region supported by finite stress behaves as r^g with g >= -1 inward)
could not be decided. The inward-march runner extends the grid to one
cell above r_minus. This module reads that extended payload and reports
the criterion on the region that is actually evaluated, without
extrapolating any further inward.

Inputs are the radial arrays already written by the march, so this is a
construction check on recorded finite states, not a new evolution.
"""
from __future__ import annotations

from typing import Any, Dict

import numpy as np
from numba import njit


@njit(cache=True)
def _log_slope(radius: np.ndarray, force: np.ndarray) -> np.ndarray:
    """Central log-log slope; nan where a stencil touches zero."""
    n = radius.shape[0]
    slope = np.full(n, np.nan)
    for i in range(1, n - 1):
        left = abs(force[i - 1])
        right = abs(force[i + 1])
        if left <= 0.0 or right <= 0.0:
            continue
        dlogr = np.log(radius[i + 1]) - np.log(radius[i - 1])
        if dlogr != 0.0:
            slope[i] = (np.log(right) - np.log(left)) / dlogr
    return slope


REGULAR_SLOPE_LIMIT = -1.0
REGULAR_SLOPE_TOLERANCE = 1.0e-9
INNER_FRACTION = 0.10
API_VERSION = "nsc-rn-inward-v1"


def inner_region_decision(radius: np.ndarray, radial_pressure: np.ndarray,
                          r_minus: float) -> Dict[str, Any]:
    """Steepness of the compressive force on the evaluated inner region.

    The grid now starts one cell above r_minus, so the inner ten percent
    of cells are evaluated, not extrapolated. The verdict is whether the
    measured slopes satisfy the Hayward g >= -1 bound there; the slope at
    the second cell is reported beside the minimum over the region.
    """
    radius = np.ascontiguousarray(np.asarray(radius, dtype=np.float64))
    pressure = np.ascontiguousarray(np.asarray(radial_pressure, dtype=np.float64))
    if radius.shape != pressure.shape:
        raise ValueError("radius and radial_pressure must share a shape")
    if radius.shape[0] < 3:
        raise ValueError("a logarithmic slope needs at least three grid points")
    slope = _log_slope(radius, pressure)
    inner_cells = max(2, int(round(INNER_FRACTION * radius.shape[0])))
    region = slope[1:inner_cells + 1]
    defined = region[np.isfinite(region)]
    if defined.size == 0:
        raise ValueError("no defined slope inside the inner region")
    region_min = float(np.min(defined))
    meets = bool(region_min >= REGULAR_SLOPE_LIMIT - REGULAR_SLOPE_TOLERANCE)
    gap = float(radius[0] - r_minus)
    return {
        "inner_edge": float(radius[0]),
        "r_minus": float(r_minus),
        "gap_to_horizon": gap,
        "gap_in_cells": float(gap / (radius[1] - radius[0])),
        "inner_cells_evaluated": int(inner_cells - 1),
        "throat_slope_inner": float(slope[1]),
        "inner_region_slope_min": region_min,
        "inner_region_slope_steepest_station": region_min,
        "meets_hayward_on_evaluated_domain": bool(meets),
        "regular_slope_limit": REGULAR_SLOPE_LIMIT,
        "decidable_on_this_domain": True,
        "extrapolation_below_inner_edge": "none performed",
    }