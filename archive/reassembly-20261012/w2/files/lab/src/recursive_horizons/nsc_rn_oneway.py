"""The one-way enclosure: transport direction of the recorded operator.

One finite question in one named domain. On the recorded nonperiodic
horizon-regular PG parent, do there exist modes of the state's own Dirac
operator that carry weight outward from inside the measured trap? The
2005 black-universe reading says the passage is one way: matter can be
drawn in through the horizon but the inside cannot send matter back out.
This module turns that reading into a spectral statement about the
recorded operator, station by station, with no time step taken.

Domain. Areal radius ``r`` on the recorded radial grid. The operator is
the instantaneous rate matrix the march uses in ``nsc_rn_pg.dirac_rates``
on the reconstructed coupled lapse and shift. A mode is interior when
more than half of its weight sits below the measured ``trap_radius`` of
the same state's metric (see ``nsc_rn_cone_map``).

Direction. For an instant of a mode ``(chi_plus, chi_minus)`` the
transport statistic is the signed probability current of the same chart,

    S = sum( v_out |chi_plus|^2 + v_in |chi_minus|^2 )
      / sum( |v_out| |chi_plus|^2 + |v_in| |chi_minus|^2 ),

with the speeds carrying their signs. Inside the trap BOTH characteristic
speeds are negative, so channel content there is inward weight; the
signed form is the one that reads "can the inside emit outward" instead
of merely "which channel carries the weight". ``S = +1`` is all-outward
current and ``S = -1`` all-inward. This is a statistic of the recorded
operator, not a conserved flux.

Verdict. An ``interior_cardinal_outward`` flag fires if any interior
mode carries outward weight above the live floor. Absence of that flag
plus the measured interior versus exterior direction split is the
one-way signature. A flag is a reported fact, not a proof about the
actual universe; the reverse retarded response of the trapped probe cell
is a separate, time-resolved statement and is not claimed here.

Scope. A construction check on recorded finite states. It does not
measure the actual universe, and no response across time is reconstructed
from a single instant.
"""
from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np
from numba import njit

from . import nsc_rn_cone_map as cone
from . import nsc_rn_observables as observables
from . import nsc_rn_pg as pg

API_VERSION = "nsc-rn-oneway-v1"
#: A direction statistic below this in magnitude is treated as balanced,
#: not as a directional claim. Ground in the chart's vacuum cleanliness.
DIRECTION_FLOOR = 0.05
#: A mode is interior-supported when this fraction of its weight sits
#: below the measured trap radius.
INTERIOR_FRACTION = 0.5
#: Bands for the energy dependence: edges of |Im(lambda)| in chart units.
FREQUENCY_BAND_EDGES = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, math.inf)


def operator_matrix(grid: Dict[str, Any], lapse: np.ndarray, shift: np.ndarray) -> np.ndarray:
    """Dense (2N, 2N) real rate matrix of ``nsc_rn_pg.dirac_rates``.

    The ordering is ``concat(plus_nodes, minus_nodes)`` and each column is
    the rates produced by a unit real basis vector, so ``dirac_rates`` is
    reproduced exactly by matrix-vector products on the real and imaginary
    parts separately.
    """
    points = int(grid["points"])
    size = 2 * points
    matrix = np.zeros((size, size))
    basis = np.zeros((2, points, 1))
    for column in range(size):
        basis[...] = 0.0
        component, node = divmod(column, points)
        basis[component, node, 0] = 1.0
        rates, _sat = pg.dirac_rates(grid, basis, lapse, shift, kappa=1)
        matrix[:points, column] = rates[0, :, 0]
        matrix[points:, column] = rates[1, :, 0]
    return matrix


@njit(cache=True)
def _direction_statistic(v_out: np.ndarray, v_in: np.ndarray, plus: np.ndarray,
                         minus: np.ndarray) -> float:
    signed = 0.0
    weight = 0.0
    for i in range(v_out.shape[0]):
        plus_power = plus[i].real ** 2 + plus[i].imag ** 2
        minus_power = minus[i].real ** 2 + minus[i].imag ** 2
        signed += v_out[i] * plus_power + v_in[i] * minus_power
        weight += abs(v_out[i]) * plus_power + abs(v_in[i]) * minus_power
    if weight <= 0.0:
        return 0.0
    return float(signed / weight)


@njit(cache=True)
def _interior_weight(radius: np.ndarray, trap_radius: float, plus: np.ndarray,
                     minus: np.ndarray) -> float:
    inside = 0.0
    total = 0.0
    for i in range(radius.shape[0]):
        weight = plus[i].real ** 2 + plus[i].imag ** 2 + minus[i].real ** 2 + minus[i].imag ** 2
        total += weight
        if radius[i] < trap_radius:
            inside += weight
    if total <= 0.0:
        return 0.0
    return float(inside / total)


def _mode_split(radius: np.ndarray, v_out: np.ndarray, v_in: np.ndarray,
                trap_radius: float, eigenvectors: np.ndarray) -> Dict[str, Any]:
    """Direction and interior weight of every eigenmode."""
    rows = []
    for column in range(eigenvectors.shape[1]):
        vector = eigenvectors[:, column]
        plus = vector[: radius.shape[0]]
        minus = vector[radius.shape[0]:]
        direction = _direction_statistic(v_out, v_in, plus, minus)
        interior = _interior_weight(radius, trap_radius, plus, minus)
        rows.append((direction, interior))
    return {"rows": rows}


def state_witness(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Eig-free witness: the signed current and interior weight of the recorded state.

    Cheap to compute, always available, and independent of the dense
    spectrum. ``station_oneway`` remains the spectral measure.
    """
    grid = cone.station_grid(payload)
    metric = cone.coupled_metric(payload, grid)
    if not metric["reconstruction_is_the_state_metric"]:
        raise ValueError("reconstructed metric is not the state's metric")
    geometry = metric["geometry"]
    profile = cone.cone_profile(grid, geometry)
    lapse = np.asarray(geometry["lapse"])
    shift = np.asarray(geometry["shift"])
    v_out, v_in = observables.chart_speeds(lapse, shift)
    phi = observables.require_phi(payload["coupled_phi"], grid["points"])
    plus = phi[0, :, 0]
    minus = phi[1, :, 0]
    return {
        "api_version": API_VERSION,
        "signed_current": float(_direction_statistic(v_out, v_in, plus, minus)),
        "interior_weight": float(_interior_weight(
            grid["radius"], profile["trap_radius"], plus, minus)),
        "trap_radius": profile["trap_radius"],
        "witness_is_eig_free": True,
    }


def station_oneway(payload: Dict[str, Any], mass_over_rm: float | None = None) -> Dict[str, Any]:
    """Operator, eigenmodes, and the transport-direction split of one station."""
    grid = cone.station_grid(payload)
    metric = cone.coupled_metric(payload, grid)
    if not metric["reconstruction_is_the_state_metric"]:
        raise ValueError("reconstructed metric is not the state's metric")
    geometry = metric["geometry"]
    profile = cone.cone_profile(grid, geometry)
    if not np.isfinite(profile["trap_radius"]):
        raise ValueError("no measured trap radius: the chart never flips outgoing light")

    lapse = np.ascontiguousarray(geometry["lapse"])
    shift = np.ascontiguousarray(geometry["shift"])
    radius = grid["radius"]
    v_out, v_in = observables.chart_speeds(lapse, shift)

    matrix = operator_matrix(grid, lapse, shift)
    eigenvalues, eigenvectors = np.linalg.eig(matrix)

    split = _mode_split(radius, v_out, v_in, profile["trap_radius"], eigenvectors)
    directions = np.asarray([s for s, _ in split["rows"]], dtype=float)
    interiors = np.asarray([w for _, w in split["rows"]], dtype=float)
    interior_mask = interiors >= INTERIOR_FRACTION

    interior_max = float(np.max(directions[interior_mask])) if np.any(interior_mask) else float("nan")
    exterior_min = float(np.min(directions[~interior_mask])) if np.any(~interior_mask) else float("nan")

    bands = []
    magnitudes = np.abs(eigenvalues.imag)
    for left, right in zip(FREQUENCY_BAND_EDGES[:-1], FREQUENCY_BAND_EDGES[1:]):
        band_mask = (magnitudes >= left) & (magnitudes < right)
        in_band = band_mask & interior_mask
        out_band = band_mask & ~interior_mask
        bands.append({
            "band": [float(left), float(right)],
            "interior_max_direction": float(np.max(directions[in_band])) if np.any(in_band) else None,
            "exterior_min_direction": float(np.min(directions[out_band])) if np.any(out_band) else None,
        })

    flags: list[str] = []
    if np.isfinite(interior_max) and interior_max > DIRECTION_FLOOR:
        flags.append("interior_cardinal_outward")
    if not np.isfinite(exterior_min):
        flags.append("no_exterior_modes")

    return {
        "api_version": API_VERSION,
        "verdict": "oneway_measured",
        "flags": flags,
        "t_over_rm": float(np.asarray(payload["coupled_t"])[0]),
        "trap_radius": profile["trap_radius"],
        "points": int(grid["points"]),
        "mode_count": int(matrix.shape[0]),
        "interior_mode_count": int(np.count_nonzero(interior_mask)),
        "interior_max_direction": interior_max,
        "exterior_min_direction": exterior_min,
        "direction_floor": DIRECTION_FLOOR,
        "interior_fraction": INTERIOR_FRACTION,
        "bands": bands,
        "eigenvalue_spectral_radius": float(np.max(np.abs(eigenvalues))),
        "max_real_eigenvalue": float(np.max(eigenvalues.real)),
        "scope": (
            "The instantaneous transport direction of the recorded operator. "
            "It is not a retarded response across time and it does not measure "
            "the actual universe."
        ),
    }