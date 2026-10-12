"""The cone map: what light's value is after it rides the recorded space.

One finite question in one named domain. On the areal radius of the
recorded nonperiodic horizon-regular PG parent (the station payloads of
``nsc-rn-neutral-family-v1``), what does the state's own metric say the
light cone is, and how much has the value of light changed after riding
the space from the outer boundary down to each radius?

Domain. Areal radius ``r`` on the recorded radial grid. Observer is the
chart of the recorded coupled arm. The state is the payload already on
disk: no time step is taken and no file is rewritten.

The metric. The coupled arm's lapse ``N`` and shift ``beta`` are not
stored in the station payload (only the frozen control's are), but they
are deterministic: the march rebuilt them with
``nsc_rn_pg.reconstruct_geometry`` at weight 1 from the same
``coupled_phi``, ``occupations`` and ``coupled_inner_mass`` that the
payload does store. This module reuses exactly that route, like the
source-free ``_measure_case`` check of the parent driver, and it
refuses to report a station where the reconstructed mass does not
reproduce the stored ``coupled_mass``.

The cone facts that are recorded, never chosen:

* ``local_speed`` -- the null cone is one at every radius. This is the
  invariant statement: light always rides its own space at unit speed.
* ``chart_speeds`` -- ``v_out = -beta + N`` and ``v_in = -beta - N``,
  the values the chart gives. They change with radius because the
  recorded lapse and shift carry the room's stress. This is the "the
  value changes with the space it rides" reading, made into a list of
  measured numbers, not into a law.
* ``ray_frequency_ratio`` -- ``N(r_out)/N(r)`` for a conserved Killing
  frequency, the dilution of light's energy after riding the room. The
  archived Bronnikov-Fabris two-room comparison measured the same
  observable class as ``E/(1+beta)`` on the 2005 background; this is the
  sibling measurement on the live magnetic RN payloads.
* ``trap_radius`` -- the areal radius where the recorded metric's
  outgoing speed crosses zero. This is the chart's own one-way
  interface, measured from the same state as everything else. The
  source-free reference value is the RN outer horizon ``r_plus``.
* ``redshift_at_excision`` -- ``1/N(r_inner)``, how much the value of
  light has changed after riding the space all the way to the inner
  edge. The headline dimensionless invariant, reported per case and per
  station, on evaluated grids only.

Dead ends. No extrapolation below the recorded inner edge is performed:
the inner edge sits outside the RN inner horizon in every recorded case,
and the extrapolation lesson of the two-force balance applies here too.
The log-slope of the lapse inside the throat is reported as a flag, not
as a verdict, for the same reason.

Scope. A construction check on recorded finite states. It does not
measure the actual universe. The numbers below are what the recorded
state's metric says about light; "c is a constant whose value changes
with the space it rides" is a reading of them, and it stays a reading.
"""
from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np
from numba import njit

from . import nsc_rn_observables as observables
from . import nsc_rn_pg as pg

API_VERSION = "nsc-rn-cone-map-v1"

#: A reconstructed mass that does not reproduce the stored coupled mass
#: to this absolute tolerance is refused: the geometry would not be the
#: state's own metric. The march integrates with 1e-12 tolerances, so
#: anything above 1e-9 means the payload belongs to a different route.
MASS_RECONSTRUCTION_TOLERANCE = 1.0e-9
#: Both chart speeds below this at a node is the numerical floor of a
#: "speed zero" reading, not a verdict.
SPEED_FLOOR = 1.0e-12
#: Lapse slope window, in cells, judged from the inner edge. A reported
#: flag only; the source-free RN lapse diverges at the singularity, and
#: this grid never reaches it.
THROAT_FRACTION = 0.10
#: The source-free reference horizon comparison. Not a fitted constant.
MAGNETIC_R2 = 1.0


@njit(cache=True)
def _first_upcrossing(value: np.ndarray) -> int:
    """First index where ``value`` crosses zero to positive, or -1."""
    for i in range(1, value.shape[0]):
        if value[i - 1] < 0.0 and value[i] >= 0.0:
            return i
    return -1


@njit(cache=True)
def _interpolated_root(radius: np.ndarray, value: np.ndarray, index: int) -> float:
    """Linear root of ``value`` on the cell ending at ``index``."""
    left = value[index - 1]
    right = value[index]
    if right == left:
        return radius[index - 1]
    weight = (0.0 - left) / (right - left)
    return radius[index - 1] + weight * (radius[index] - radius[index - 1])


@njit(cache=True)
def _cumulative_inward(radius: np.ndarray, ingoing: np.ndarray) -> np.ndarray:
    """T(r) = int of dr/|v_in| accumulated inward from the outer edge."""
    count = radius.shape[0]
    crossed = np.zeros(count)
    accumulator = 0.0
    for i in range(count - 1):
        index = count - 2 - i
        step = radius[index + 1] - radius[index]
        accumulator += 0.5 * step * (1.0 / abs(ingoing[index + 1]) + 1.0 / abs(ingoing[index]))
        crossed[index] = accumulator
    return crossed


def station_grid(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Rebuild the production grid that the march used for this payload.

    The areal grid is uniform with the recorded spacing; the SBP weights
    and derivative come from the same ``nsc_rn_pg.sbp_operator`` that the
    march used, and the charge is the recorded P^2 = 1 magnetic radius.
    """
    radius = np.ascontiguousarray(np.asarray(payload["radius"], dtype=np.float64))
    if radius.ndim != 1 or radius.shape[0] < 8:
        raise ValueError("station grid needs a radial array of at least eight points")
    spacing_values = np.diff(radius)
    if np.max(np.abs(spacing_values - spacing_values[0])) > 1e-9 * spacing_values[0]:
        raise ValueError("station grid is not uniform in the areal coordinate")
    spacing = float(spacing_values[0])
    weights, derivative, _adjoint = pg.sbp_operator(int(radius.shape[0]), spacing, layout="csr")
    return {
        "radius": radius,
        "points": int(radius.shape[0]),
        "spacing": spacing,
        "weights": weights,
        "derivative": derivative,
        "charge": 1.0,
    }


def coupled_metric(payload: Dict[str, Any], grid: Dict[str, Any]) -> Dict[str, Any]:
    """The coupled arm's lapse and shift, reconstructed from the payload's own state.

    Uses the same route as ``nsc_rn_pg.stage_rates``: matter densities
    from the stored occupations at N=1, then ``reconstruct_geometry`` at
    weight 1 with the stored inner mass. Returns the geometry dict plus
    the reconstruction error against the stored ``coupled_mass``.
    """
    phi = observables.require_phi(payload["coupled_phi"], grid["points"])
    occupations = np.asarray(payload["occupations"], dtype=np.float64)
    inner_mass = float(payload["coupled_inner_mass"][0])
    densities = pg.matter_densities(phi, grid, occupations)
    geometry = pg.reconstruct_geometry(
        grid, inner_mass, densities["eps_N"], densities["eps_beta"],
        densities["angular"], weight=1.0,
        force_r_unit=densities.get("force_r_unit"),
    )
    stored_mass = np.ascontiguousarray(np.asarray(payload["coupled_mass"], dtype=np.float64))
    mass_residual = float(np.max(np.abs(geometry["mass"] - stored_mass)))
    return {
        "geometry": geometry,
        "mass_residual": mass_residual,
        "mass_residual_tolerance": MASS_RECONSTRUCTION_TOLERANCE,
        "reconstruction_is_the_state_metric": mass_residual <= MASS_RECONSTRUCTION_TOLERANCE,
    }


def cone_profile(grid: Dict[str, Any], geometry: Dict[str, Any]) -> Dict[str, Any]:
    """The local and chart light speeds and the accumulated riding value."""
    radius = grid["radius"]
    lapse = np.ascontiguousarray(np.asarray(geometry["lapse"], dtype=np.float64))
    shift = np.ascontiguousarray(np.asarray(geometry["shift"], dtype=np.float64))
    outgoing, ingoing = observables.chart_speeds(lapse, shift)

    crossing = _first_upcrossing(outgoing)
    trap_radius = float("nan")
    if crossing > 0:
        trap_radius = float(_interpolated_root(radius, outgoing, crossing))

    redshift = 1.0 / lapse
    crossing_coordinate = _cumulative_inward(radius, ingoing)
    delay = _cumulative_inward(radius, ingoing) - (radius[-1] - radius)

    throat = max(2, int(round(THROAT_FRACTION * radius.shape[0])))
    log_r = np.log(radius)
    dlogr = log_r[throat] - log_r[1]
    lapse_slope = (math.log(lapse[throat]) - math.log(lapse[1])) / dlogr if dlogr != 0.0 else math.nan

    return {
        "local_speed_everywhere": 1.0,
        "v_out": outgoing,
        "v_in": ingoing,
        "beta": shift,
        "lapse": lapse,
        "ray_frequency_ratio": redshift,
        "shift_factor": geometry["shift_factor"],
        "trapped_below_index": int(crossing),
        "trap_radius": trap_radius,
        "crossing_coordinate_time": crossing_coordinate,
        "light_delay_vs_flat": delay,
        "lapse_log_slope_at_throat": float(lapse_slope),
        "lapse_min": float(np.min(lapse)),
        "lapse_at_inner": float(lapse[0]),
        "redshift_at_excision": float(redshift[0]),
        "chart_speeds_at_outer": [float(outgoing[-1]), float(ingoing[-1])],
        "chart_speeds_at_inner": [float(outgoing[0]), float(ingoing[0])],
    }


def reference_r_plus(mass_over_rm: float) -> float:
    """Source-free RN outer horizon for P^2 = 1: mu + sqrt(mu^2 - 1)."""
    argument = float(mass_over_rm) ** 2 - MAGNETIC_R2
    if argument <= 0.0:
        raise ValueError("no RN horizon at or below the magnetic radius")
    return float(mass_over_rm) + math.sqrt(argument)


def station_assessment(payload: Dict[str, Any], mass_over_rm: float | None = None) -> Dict[str, Any]:
    """Cone facts plus the reconstruction check for one recorded station."""
    grid = station_grid(payload)
    metric = coupled_metric(payload, grid)
    if not metric["reconstruction_is_the_state_metric"]:
        raise ValueError(
            "reconstructed mass does not match the stored coupled mass "
            f"(residual {metric['mass_residual']:.3e} > {MASS_RECONSTRUCTION_TOLERANCE:.1e})"
        )
    geometry = metric["geometry"]
    profile = cone_profile(grid, geometry)

    flags: list[str] = []
    if profile["trapped_below_index"] <= 0:
        flags.append("no_trapped_strip")
    if not np.isfinite(profile["trap_radius"]):
        flags.append("trap_radius_unresolved")
    if profile["redshift_at_excision"] < 1.0 + 1e-12:
        flags.append("no_riding_dilution")

    comparison: Dict[str, Any] = {
        "reference_r_plus": None,
        "trap_minus_reference": None,
    }
    if mass_over_rm is not None:
        r_plus = reference_r_plus(mass_over_rm)
        comparison = {
            "reference_r_plus": float(r_plus),
            "trap_minus_reference": float(profile["trap_radius"] - r_plus),
            "reference_is_source_free": True,
        }

    return {
        "api_version": API_VERSION,
        "verdict": "cone_map_measured",
        "flags": flags,
        "mass_residual": metric["mass_residual"],
        "t_over_rm": float(np.asarray(payload["coupled_t"])[0]),
        "profile": profile,
        "horizon_comparison": comparison,
        "scope": (
            "What the recorded state's own metric says about light on the "
            "evaluated grid. It does not measure the actual universe, and it "
            "performs no extrapolation below the recorded inner edge."
        ),
    }