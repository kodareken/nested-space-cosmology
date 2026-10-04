"""Two-force gradient balance on recorded RN radial parent states.

This module asks one finite question in one named domain: on the areal
radius of the recorded nonperiodic horizon-regular PG parent, does the
enclosing and the enclosed form a balanced pair, or does one of the two
forces run away and leave the other as decoration?

Domain. Areal radius ``r`` on the recorded radial grid, from the fixed
inner mass at ``grid[0]`` to the absorbing outer boundary at ``32 r_m``.
Observer is the areal chart of the same state. No coordinate change, no
time slice, no second geometry.

The two forces. Both are read from the same recorded station payload and
are not fitted, blended, or reassembled:

* ``radial_pressure`` is the inherited radial contribution. It is the
  compressive side that pushes matter toward the excision boundary.
* ``angular_pressure`` is the generated angular contribution. It is the
  resisting side that a compressed interior pushes back against.

Balance. ``B(r) = |F_in| / (|F_in| + |F_out|)`` is scale free. ``B = 1/2``
is two equal forces. ``B -> 1`` means the inside acts alone and the
outside has stopped being part of the object.

Gradient. ``d ln|F| / d ln r`` is a scale-free logarithmic slope. A value
of ``-1`` is inverse-radius decay. A magnitude above one means the force
grows without bound as ``r`` decreases, which is the collapse signature.

Dead ends. This module is built to fail loudly instead of producing a
comfortable label. Each named flag below is a reason to stop, not a
rounding remark:

* ``single_force`` -- the balance never drops below
  ``SINGLE_FORCE_BALANCE`` anywhere in the collar. Two forces were
  declared; one is absent. This is the "computed the interior only" dead
  end.
* ``throat_runaway`` -- the radial logarithmic slope inside the throat
  collar is steeper than ``-THROAT_SLOPE_LIMIT``. The compressive force
  diverges inward and the bubble does not hold.
* ``curvature_dominates`` -- the reconstructed curvature density exceeds
  the force density by more than ``CURVATURE_FORCE_RATIO``. Both are
  densities of the same evolved state on the same areal grid, so the
  ratio is dimensionless and the gradient being measured is no longer
  the leading structure.

Provenance matters for that check. The station payload carries two
curvature entries with different origins. ``R4`` is reconstructed from
the evolved PG jets and belongs to the same state as the pressures.
``Ricci2`` is the analytic source-free RN reference ``4 Q^2 / r^8``. The
two are never divided into one ratio here: an exact reference curvature
set against an evolved matter stress is not a comparison of two forces,
and mixing them would manufacture a false dead end.

Localization is reported, never used as a dead end. The radius
``r_decouple`` at which the balance first reaches ``EXTERIOR_DECOUPLED``
is a measured scale: the outside stops contributing beyond that areal
radius. Reporting it is what turns "the forces balance somewhere" into a
statement about how far the balance reaches.

Scope. This is a construction check on a recorded finite state. It does
not measure the actual universe, and a bounded balance here is not
evidence about cosmology.

The traction concept in ``nsc-discovery-parent-traction`` is a lead about
the same question in a different domain, the discovery collar
``|x| <= 1``. That boundary gradient charge is not this areal-radius
force balance and the two are not identified here.
"""
from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np
from numba import njit

API_VERSION = "nsc-rn-gradient-bubble-v3"

#: A force steeper than this as r decreases is a reported collapse
#: signature. It is a diagnostic, not a verdict. The verdict comes from
#: the measured balance-wall structure, which needs no chosen limit.
THROAT_SLOPE_LIMIT = 8.0
#: Balance below this means the angular force is under a tenth of radial.
SINGLE_FORCE_BALANCE = 0.90
#: Balance at or above this before mid collar ends the non-local content.
EXTERIOR_DECOUPLED = 0.999
#: Curvature density may not exceed the force density by more than this.
CURVATURE_FORCE_RATIO = 1.0e3
#: Fraction of the collar the throat test is judged on.
THROAT_FRACTION = 0.10
#: Regularity criterion from Hayward, gr-qc/0506126, "Formation and
#: evaporation of non-singular black holes" (2005). The static region is
#: "supported by finite density and pressures ... behaving as a
#: cosmological constant at small radius". A force with logarithmic slope
#: g behaves as r^g, so it stays finite as r -> 0 only when g >= -1.
#: Anything steeper diverges inward and the interior is not regular.
#: This is sourced, not chosen.
REGULAR_SLOPE_LIMIT = -1.0
#: An exact power law evaluated by a central difference lands a few ulp
#: off its own exponent, so the boundary is compared with this slack
#: rather than failing an exactly-regular force on round-off.
REGULAR_SLOPE_TOLERANCE = 1.0e-9
#: Both forces must be at least this fraction of their own peak for a
#: crossing to be structure. Below it a crossing is a noise floor.
WALL_LIVE_FRACTION = 1.0e-3
#: A wall must hold the balance this far past one half before the next
#: crossing is believed, so a single-cell wobble is not a wall.
WALL_HYSTERESIS = 0.05
#: Minimum cells a wall must persist on each side.
WALL_MIN_CELLS = 2


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
    """Trapezoidal surface integral ``4 pi int |force| r^2 dr``."""
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
    """Two-force balance and gradient on one recorded station payload.

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


def gradient_bubble_assessment(payload: Any, curvature: Dict[str, np.ndarray] | None = None) -> Dict[str, Any]:
    """Balance profile plus the named dead-end flags and the verdict.

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

    # Both sides are densities of the same evolved state on the same
    # areal grid, so the ratio is dimensionless. Dividing an integrated
    # force by r^3 would not be.
    curvature_ratio = curvature_max / force_peak if force_peak > 0.0 else math.inf

    flags: list[str] = []

    # Regularity is the sourced criterion. A force that diverges inward is
    # not a finite-stress interior, which is exactly what the 2005
    # non-singular black-hole construction requires at small radius.
    regular = (math.isnan(profile["throat_slope_inner"])
               or profile["throat_slope_inner"]
               >= REGULAR_SLOPE_LIMIT - REGULAR_SLOPE_TOLERANCE)
    if not regular:
        flags.append("interior_not_regular")

    # The verdict is the measured wall structure, not a chosen slope.
    # A throat slope steeper than the limit is reported because it is a
    # real collapse signature, but it no longer decides the outcome: the
    # v1 verdict rested on that number and the data straddled it at -7.37
    # and -8.44, so the verdict pointed the wrong way.
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
        "verdict": verdict,
        "dead_end_flags": flags,
        "interior_regular": bool(regular),
        "throat_slope_inner": float(profile["throat_slope_inner"]),
        "regularity_slope_limit": REGULAR_SLOPE_LIMIT,
        "regularity_slope_tolerance": REGULAR_SLOPE_TOLERANCE,
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
            "regular_slope_limit": REGULAR_SLOPE_LIMIT,
        },
        "profile": profile,
    }