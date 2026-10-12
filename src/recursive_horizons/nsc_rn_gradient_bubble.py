"""Radial/transverse pressure diagnostics on recorded finite RN states.

The input entries are orthonormal pressures from ``nsc_rn_source``, not
opposing forces from two regions. ``B = |p_r| / (|p_r| + |p_Omega|)``
locates equal-magnitude pressure crossings. It does not establish mechanical
equilibrium, a child boundary, exterior decoupling or a regular centre.
Slopes are local ``d ln|p| / d ln r`` measurements on a positive-radius
grid. A power-law extrapolation is an assumption, never a centre test.

Hayward gr-qc/0506126 supplies a finite-stress regular-core example, not
the former ``g >= -1`` test. Reaching an inner horizon does not reach the
centre. This consumer returns no regularity verdict.

Legacy force/wall names and diagnostic verdict strings remain for readers
of v1--v5. ``diagnostic_semantics`` and pressure-named aliases explain their
limited meaning. Old regularity fields are retained as unsupported/null.
The conventional flags concern only this diagnostic; they do not reject a
physical model. Frozen historical records are not rewritten.

``R4`` is reconstructed from the same recorded PG state. Analytic ``Ricci2``
is rejected as a different-provenance reference. The historical R4/pressure
ratio remains a reported convention in the record's units, not a new
invariant regularity or force-balance criterion.
"""
from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np
from numba import njit

API_VERSION = "nsc-rn-gradient-bubble-v6"

#: A local pressure slope steeper than this is reported, not interpreted
#: as collapse or central divergence.
THROAT_SLOPE_LIMIT = 8.0
#: Legacy threshold for a small transverse pressure relative to radial.
SINGLE_FORCE_BALANCE = 0.90
#: Legacy radius label: balance reaches this pressure-magnitude ratio.
EXTERIOR_DECOUPLED = 0.999
#: Conventional R4/pressure diagnostic threshold in the recorded units.
CURVATURE_FORCE_RATIO = 1.0e3
#: Fraction of the collar the throat test is judged on.
THROAT_FRACTION = 0.10
#: Deprecated v1--v5 value, retained only for import compatibility.
#: Unsupported by Hayward and never used in a v6 scientific assessment.
REGULAR_SLOPE_LIMIT = -1.0
#: The grid's inner edge counts as reaching the horizon only if it lies
#: this close to the inner root. This says nothing about the centre.
HORIZON_REACH_TOLERANCE = 1.0e-6
#: Deprecated companion to REGULAR_SLOPE_LIMIT; no longer used.
REGULAR_SLOPE_TOLERANCE = 1.0e-9
#: Both pressures must be at least this fraction of their own peak for a
#: crossing to be structure. Below it a crossing is a noise floor.
WALL_LIVE_FRACTION = 1.0e-3
#: A wall must hold the balance this far past one half before the next
#: crossing is believed, so a single-cell wobble is not a wall.
WALL_HYSTERESIS = 0.05
#: Minimum cells a wall must persist on each side.
WALL_MIN_CELLS = 2


def diagnostic_semantics() -> Dict[str, Any]:
    """Interpretation contract for live consumers and successor records."""
    return {
        "quantity": "orthonormal radial and transverse pressures",
        "balance": "abs(p_r) / (abs(p_r) + abs(p_Omega))",
        "wall": "equal-magnitude pressure crossing, not mechanical equilibrium",
        "integral": "4 pi integral abs(p) r^2 dr; volume-weighted pressure",
        "regularity_status": "unsupported_by_this_diagnostic",
        "regularity_criterion": None,
        "horizon_reach_determines_centre_regularity": False,
        "legacy_verdicts": "diagnostic labels only, not physical dead-end verdicts",
        "deprecated_fields": [
            "slope_meets_regularity_limit", "regularity_slope_limit",
            "regularity_slope_tolerance", "regularity_decidable_on_this_domain",
        ],
    }


@njit(cache=True)
def _log_slope(radius: np.ndarray, force: np.ndarray) -> np.ndarray:
    """Central log-log slope of ``force`` on a uniform areal grid.

    Returns ``nan`` where the stencil touches a zero, because the slope of
    a force that vanishes is not a number and must not be smoothed away.
    """
    n = radius.shape[0]
    slope = np.full(n, np.nan)
    for i in range(1, n - 1):
        a = abs(force[i - 1])
        b = abs(force[i + 1])
        if a <= 0.0 or b <= 0.0:
            continue
        dlogr = math.log(radius[i + 1]) - math.log(radius[i - 1])
        if dlogr == 0.0:
            continue
        slope[i] = (math.log(b) - math.log(a)) / dlogr
    return slope


@njit(cache=True)
def _sign_changes(force: np.ndarray) -> int:
    """Count radial sign changes of the signed force."""
    count = 0
    previous = 0
    for i in range(force.shape[0]):
        if force[i] > 0.0:
            sign = 1
        elif force[i] < 0.0:
            sign = -1
        else:
            continue
        if previous != 0 and sign != previous:
            count += 1
        previous = sign
    return count


@njit(cache=True)
def _integrated_force(radius: np.ndarray, force: np.ndarray) -> float:
    """Legacy name: volume-weighted pressure ``4 pi int |p| r^2 dr``."""
    total = 0.0
    for i in range(radius.shape[0] - 1):
        r0 = radius[i]
        r1 = radius[i + 1]
        f0 = abs(force[i])
        f1 = abs(force[i + 1])
        total += 0.5 * (f0 * r0 * r0 + f1 * r1 * r1) * (r1 - r0)
    return 4.0 * math.pi * total


@njit(cache=True)
def _balance_profile(radius: np.ndarray, inner: np.ndarray, outer: np.ndarray) -> np.ndarray:
    """``|F_in| / (|F_in| + |F_out|)`` pointwise."""
    n = radius.shape[0]
    balance = np.zeros(n)
    for i in range(n):
        a = abs(inner[i])
        b = abs(outer[i])
        denom = a + b
        if denom > 0.0:
            balance[i] = a / denom
    return balance


@njit(cache=True)
def _balance_walls(radius: np.ndarray, balance: np.ndarray, live: np.ndarray) -> np.ndarray:
    """Areal radii where the two forces are equal, inside the live collar.

    A crossing counts only when the balance genuinely holds beyond the
    hysteresis band on both sides and both forces are still above their
    own live floor. A single-cell wobble is not a wall, and a crossing
    out where one force has already died is not a wall either. The wall
    radius is interpolated where the balance passes one half, not taken
    at a grid node.
    """
    walls = []
    side = 0
    run = 0
    previous = -1
    candidate = -1.0
    for i in range(radius.shape[0]):
        if not live[i]:
            continue
        value = balance[i]
        if side == 0:
            if value >= 0.5 + WALL_HYSTERESIS:
                side = 1
                run = 1
                previous = i
            elif value <= 0.5 - WALL_HYSTERESIS:
                side = -1
                run = 1
                previous = i
            continue
        flipped = (side == 1 and value <= 0.5 - WALL_HYSTERESIS) or (
            side == -1 and value >= 0.5 + WALL_HYSTERESIS)
        if flipped:
            # The old side persisted; the new side still has to. Hold the
            # crossing as a candidate and only keep it if it survives.
            candidate = -1.0
            if run >= WALL_MIN_CELLS and previous >= 0 and previous < i:
                before = balance[previous]
                if value != before:
                    weight = (0.5 - before) / (value - before)
                    weight = min(max(weight, 0.0), 1.0)
                    candidate = radius[previous] + weight * (radius[i] - radius[previous])
            side = -side
            run = 1
            previous = i
        else:
            run += 1
            previous = i
            if run >= WALL_MIN_CELLS and candidate >= 0.0:
                walls.append(candidate)
                candidate = -1.0
    return np.array(walls, dtype=np.float64)


@njit(cache=True)
def _finite_extreme(values: np.ndarray, want_max: bool) -> float:
    """Largest or smallest finite entry, ignoring ``nan`` stencils."""
    best = math.nan
    for i in range(values.shape[0]):
        v = values[i]
        if math.isnan(v):
            continue
        if math.isnan(best):
            best = v
        elif want_max:
            if v > best:
                best = v
        else:
            if v < best:
                best = v
    return best


def gradient_bubble_profile(payload: Any) -> Dict[str, Any]:
    """Relative pressure magnitudes and slopes on a recorded station.

    ``payload`` is a mapping with ``radius``, ``radial_pressure`` and
    ``angular_pressure``. Nothing is rescaled and no value is fitted.
    """
    radius = np.ascontiguousarray(np.asarray(payload["radius"], dtype=np.float64))
    inner = np.ascontiguousarray(np.asarray(payload["radial_pressure"], dtype=np.float64))
    outer = np.ascontiguousarray(np.asarray(payload["angular_pressure"], dtype=np.float64))

    if radius.shape != inner.shape or radius.shape != outer.shape:
        raise ValueError("radius, radial_pressure and angular_pressure must share a shape")
    if radius.shape[0] < 3:
        raise ValueError("a logarithmic slope needs at least three grid points")
    dr = np.diff(radius)
    if not np.all(dr > 0.0):
        raise ValueError("the areal grid must be strictly increasing")
    if not np.all(np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("logarithmic pressure slopes require finite positive radii")

    balance = _balance_profile(radius, inner, outer)
    slope_inner = _log_slope(radius, inner)
    slope_outer = _log_slope(radius, outer)

    inner_peak = float(np.max(np.abs(inner)))
    outer_peak = float(np.max(np.abs(outer)))
    live = (np.abs(inner) >= WALL_LIVE_FRACTION * inner_peak) & (
        np.abs(outer) >= WALL_LIVE_FRACTION * outer_peak)
    walls = _balance_walls(radius, balance, live)

    throat = max(1, int(round(THROAT_FRACTION * radius.shape[0])))
    throat_slice = slice(1, throat + 1)
    mid = radius.shape[0] // 2

    defined = ~np.isnan(slope_inner)
    throat_defined = defined[throat_slice]
    throat_slope = _finite_extreme(np.ascontiguousarray(slope_inner[throat_slice]), False)

    decoupled_at = -1
    for i in range(mid):
        if balance[i] >= EXTERIOR_DECOUPLED:
            decoupled_at = i
            break

    return {
        "api_version": API_VERSION,
        "diagnostic_semantics": diagnostic_semantics(),
        "points": int(radius.shape[0]),
        "r_inner": float(radius[0]),
        "r_outer": float(radius[-1]),
        "balance": balance,
        "slope_inner": slope_inner,
        "slope_outer": slope_outer,
        "balance_min": float(np.min(balance)),
        "balance_max": float(np.max(balance)),
        "balance_at_throat": float(balance[1]),
        "balance_at_mid": float(balance[mid]),
        "balance_at_outer": float(balance[-1]),
        "throat_points": int(throat),
        "throat_slope_defined": int(np.count_nonzero(throat_defined)),
        "throat_slope_inner": float(throat_slope),
        "slope_outer_at_mid": float(slope_outer[mid]),
        "force_inner_integrated": float(_integrated_force(radius, inner)),
        "force_outer_integrated": float(_integrated_force(radius, outer)),
        "radial_pressure_volume_integral": float(_integrated_force(radius, inner)),
        "transverse_pressure_volume_integral": float(_integrated_force(radius, outer)),
        "force_ratio_integrated": float(
            _integrated_force(radius, outer) / _integrated_force(radius, inner)
        ),
        "sign_changes_inner": int(_sign_changes(inner)),
        "sign_changes_outer": int(_sign_changes(outer)),
        "exterior_decoupled_index": int(decoupled_at),
        "r_decouple": float(radius[decoupled_at]) if decoupled_at >= 0 else float(radius[-1]),
        "walls": walls,
        "wall_count": int(walls.size),
        "r_first_wall": float(walls[0]) if walls.size else float("nan"),
        "r_last_wall": float(walls[-1]) if walls.size else float("nan"),
        "balance_side_at_throat": int(np.sign(balance[1] - 0.5)) if balance[1] != 0.5 else 0,
    }


def inner_horizon(mass_over_rm: float, magnetic_r2: float = 1.0) -> float:
    """Smaller root of r^2 - 2 M r + magnetic_r2 = 0, the RN inner horizon.

    Recorded with `P^2 = magnetic_r2 = 1` for the classical flux-1
    sector, which is the sector every family case was run in.
    """
    discriminant = mass_over_rm * mass_over_rm - magnetic_r2
    if discriminant <= 0.0:
        raise ValueError("no inner horizon: the mass is below the charge")
    return mass_over_rm - math.sqrt(discriminant)


def horizon_extrapolation(r_inner: float, r_minus: float, slope: float) -> Dict[str, Any]:
    """What an assumed constant pressure power law does to the horizon.

    The computed domain stops at `r_inner`, the initial fixed inner mass,
    which lies outside `r_minus`. This reports the gap the computation
    never covers and what that power law would do across it. This is an
    extrapolation, not a measurement or a regularity test. Reaching the
    inner horizon does not establish the limit toward the centre.
    """
    gap = float(r_inner - r_minus)
    reaches = gap <= HORIZON_REACH_TOLERANCE
    extrapolated = 0.0
    if gap > 0.0 and not math.isnan(slope):
        extrapolated = float(slope * math.log(r_minus / r_inner) / math.log(10.0))
    return {
        "r_inner": float(r_inner),
        "r_minus": float(r_minus),
        "gap_to_horizon": gap,
        "domain_reaches_horizon": bool(reaches),
        "extrapolated_log10_force_ratio": extrapolated,
        "extrapolation_is_measurement": False,
        "centre_regularity_determined": False,
    }


def gradient_bubble_assessment(payload: Any, curvature: Dict[str, np.ndarray] | None = None,
                                mass_over_rm: float | None = None) -> Dict[str, Any]:
    """Pressure profile with legacy diagnostic flags and scoped semantics.

    ``curvature`` must carry only same-provenance entries. ``Ricci2`` is
    the analytic source-free reference and is rejected by name rather
    than silently averaged in.
    """
    profile = gradient_bubble_profile(payload)

    force_peak = float(np.max(np.abs(np.asarray(payload["radial_pressure"], dtype=np.float64))))

    rejected: list[str] = []
    curvature_max = 0.0
    if curvature:
        for key, block in curvature.items():
            if key != "R4":
                rejected.append(key)
                continue
            values = np.asarray(block, dtype=np.float64)
            curvature_max = max(curvature_max, float(np.max(np.abs(values))))

    # Preserve the historical diagnostic in the record's fixed units.
    # Same provenance alone does not establish force or regularity meaning.
    curvature_ratio = curvature_max / force_peak if force_peak > 0.0 else math.inf

    flags: list[str] = []

    extrapolation: Dict[str, Any] = {
        "r_inner": profile["r_inner"],
        "r_minus": None,
        "gap_to_horizon": None,
        "domain_reaches_horizon": False,
        "extrapolated_log10_force_ratio": None,
        "extrapolation_is_measurement": False,
        "centre_regularity_determined": False,
    }
    if mass_over_rm is not None:
        extrapolation = horizon_extrapolation(
            profile["r_inner"], inner_horizon(mass_over_rm), profile["throat_slope_inner"])

    # Local slope is reported without extrapolating a centre singularity.
    if (not math.isnan(profile["throat_slope_inner"])) and profile["throat_slope_inner"] < -THROAT_SLOPE_LIMIT:
        flags.append("throat_steep_reported")
    if profile["balance_min"] >= SINGLE_FORCE_BALANCE:
        flags.append("single_force")
    if curvature_ratio > CURVATURE_FORCE_RATIO:
        flags.append("curvature_dominates")

    hard = [f for f in flags if f != "throat_steep_reported"]

    if profile["wall_count"] > 0:
        verdict = "two_force_wall"
    elif hard:
        verdict = "dead_end"
    elif profile["balance_min"] <= 0.0 or profile["balance_max"] >= 1.0:
        verdict = "unresolved"
    else:
        verdict = "single_force_collar"

    return {
        "api_version": API_VERSION,
        "diagnostic_semantics": diagnostic_semantics(),
        "verdict": verdict,
        "dead_end_flags": flags,
        "slope_meets_regularity_limit": None,
        "regularity_decidable_on_this_domain": False,
        "regularity_extrapolation": extrapolation,
        "throat_slope_inner": float(profile["throat_slope_inner"]),
        "regularity_slope_limit": None,
        "regularity_slope_tolerance": None,
        "curvature_max": float(curvature_max),
        "curvature_keys_rejected": rejected,
        "force_peak_density": force_peak,
        "curvature_force_ratio": float(curvature_ratio),
        "limits": {
            "throat_slope_limit": THROAT_SLOPE_LIMIT,
            "single_force_balance": SINGLE_FORCE_BALANCE,
            "exterior_decoupled": EXTERIOR_DECOUPLED,
            "curvature_force_ratio": CURVATURE_FORCE_RATIO,
            "throat_fraction": THROAT_FRACTION,
            "wall_live_fraction": WALL_LIVE_FRACTION,
            "wall_hysteresis": WALL_HYSTERESIS,
            "wall_min_cells": WALL_MIN_CELLS,
            "regular_slope_limit": None,
        },
        "profile": profile,
    }
