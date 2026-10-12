"""The two ends: one room's bottom and top, read from the same local space.

The programme reading is that Planck- and zero-like scales are not global
numbers but the two ends of the local space: the bottom end is the lowest
coherent frequency a room supports and the top end is its resolution. Both
are read from the same recorded state's metric, so their ratio is a
per-room dimensionless number, and whether that number changes from room
to room is a measured fact, not an assumption.

* Bottom end. The angular tower of the Dirac shell costs ``kappa/r`` in
  the local frame, so in chart units the lowest coherent frequency of the
  room at radius ``r`` is ``omega_zero(r) = N / r``. This is the zero-point
  of the spherical structure: the cheapest note the room can hold.
* Top end. The room resolves structure on its own lattice. In proper
  units the highest travel frequency is the characteristic-scale inverse
  ``omega_cutoff(r) = (N + |beta|) / spacing``, the CFL resolution the
  march itself respects.
* Room quality. ``Q(r) = omega_zero(r) / omega_cutoff(r)``, dimensionless
  and scale-invariant under a common rescaling of the room.

If the two ends were set globally, ``Q`` would be the same number in
every room. The record measures whether it is. For the recorded states it
is many orders of magnitude larger at the inner edge than at the outer
boundary, because ``r_inner << r_outer`` while the lapse is near one in
this flat-slice gauge: the same state holds a harder, sharper room at the
bottom than at the top.

Scale invariance. ``Q = N^2 spacing / (r (N + |beta|))`` is unchanged when
``(r, spacing)`` and the lapse save the same dimensionless profile, which
is the self-similarity the room-quality reading rests on.

Scope. A construction check on recorded finite states. It does not derive
Planck's constant or absolute zero, and it does not measure the actual
universe. It says what the recorded states' rooms say at their two ends.
"""
from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np

from . import nsc_rn_cone_map as cone

API_VERSION = "nsc-rn-two-ends-v1"


def local_room_scales(radius: np.ndarray, lapse: np.ndarray, shift: np.ndarray,
                      spacing: float) -> Dict[str, np.ndarray]:
    """Bottom end, top end and room quality as radial profiles."""
    if spacing <= 0.0 or not math.isfinite(spacing):
        raise ValueError("spacing must be positive and finite")
    lapse = np.asarray(lapse, dtype=float)
    shift = np.asarray(shift, dtype=float)
    radius = np.asarray(radius, dtype=float)
    if np.any(lapse <= 0.0) or np.any(radius <= 0.0):
        raise ValueError("lapse and radius must be positive")
    omega_zero = lapse / radius
    omega_cutoff = (lapse + np.abs(shift)) / spacing
    with np.errstate(divide="ignore", invalid="ignore"):
        quality = omega_zero / omega_cutoff
    return {
        "omega_zero": omega_zero,
        "omega_cutoff": omega_cutoff,
        "room_quality": quality,
    }


def _index_nearest(radius: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(np.asarray(radius, dtype=float) - target)))


def station_two_ends(payload: Dict[str, Any], mass_over_rm: float | None = None) -> Dict[str, Any]:
    """Room-quality profile and its two-end scalars for one recorded station."""
    grid = cone.station_grid(payload)
    metric = cone.coupled_metric(payload, grid)
    if not metric["reconstruction_is_the_state_metric"]:
        raise ValueError("reconstructed metric is not the state's metric")
    geometry = metric["geometry"]
    radius = grid["radius"]
    scales = local_room_scales(radius, geometry["lapse"], geometry["shift"], grid["spacing"])
    quality = scales["room_quality"]

    index_outer = int(radius.shape[0]) - 1
    index_ref = _index_nearest(radius, cone.reference_r_plus(mass_over_rm)) if mass_over_rm else int(radius.shape[0]) // 2
    index_mid = int(radius.shape[0]) // 2

    return {
        "api_version": API_VERSION,
        "t_over_rm": float(np.asarray(payload["coupled_t"])[0]),
        "trap_radius": cone.cone_profile(grid, geometry)["trap_radius"],
        "room_quality_inner": float(quality[0]),
        "room_quality_at_r_plus": float(quality[index_ref]),
        "room_quality_mid": float(quality[index_mid]),
        "room_quality_outer": float(quality[index_outer]),
        "quality_ratio_inner_over_outer": float(quality[0] / quality[index_outer]),
        "quality_ratio_inner_over_at_r_plus": float(quality[0] / quality[index_ref]),
        "omega_zero_inner": float(scales["omega_zero"][0]),
        "omega_cutoff_inner": float(scales["omega_cutoff"][0]),
        "omega_zero_outer": float(scales["omega_zero"][index_outer]),
        "omega_cutoff_outer": float(scales["omega_cutoff"][index_outer]),
        "dimensionless_invariant_is_per_room": bool(
            abs(float(quality[0] / quality[index_outer]) - 1.0) > 1.0e-6
        ),
        "scale_invariant_definition": "N^2 spacing / (r (N + |beta|))",
        "scope": (
            "the two ends of the recorded room, read from the same local "
            "space; not Planck's constant, not absolute zero, not the universe"
        ),
    }