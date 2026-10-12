"""Per-level influence fee extracted from recorded NSC bytes, and its closure.

Saved-data consumer only: no evolution, no new solver. Three candidate
per-level coefficients f are extracted from already-recorded evidence and
the self-similar chain

    total = L (1 + f + f**2 + ...) = L / (1 - f)

is closed under each of them. local_fraction = 1 - f is the fraction of the
closed-chain budget carried by the innermost level.

Candidates (definitions pre-registered, none fitted):

    structural_resolvent : f = 1 / Omega, the multiplicative per-link scale
        of nsc_matryoshka.resolvent_levels, with Omega taken from the
        fixed-q scale binding record.
    recorded_transfer    : f = (Omega_times_transfer) / Omega, the recursive
        endpoint transfer already recorded in the scale binding record.
    measured_mass_fee    : f = (M_coupled_final - M_fixed_final) / M_fixed
        from the neutral continuation summary, the measured one-level
        metric-reconstruction fee per unit system mass.

Band: the pre-registered acceptance region for local_fraction is
(BAND_LO, BAND_HI) around the BBN-anchored baryon share 0.049. A candidate
is NEAR only if its local_fraction lands inside the band without any
parameter being tuned.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

API_VERSION = "nsc-influence-chain-v1"

TARGET_LOCAL_FRACTION = 0.049
BAND_LO = 0.03
BAND_HI = 0.10
ROUTE_RATIO = 32.0
NEAR_STATUS = "NEAR_LOCAL_BARYON_FRACTION"
NOT_NEAR_STATUS = "NOT_NEAR"
BAND_STATUS = "INDETERMINATE"
DEPTH_DEFAULT = None  # infinite self-similar chain


class ChainFeeError(RuntimeError):
    """A chain coefficient or record that cannot enter the closure."""


def path_sha256(path: Path) -> str:
    import hashlib
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_scale(value, name):
    value = float(value)
    if not math.isfinite(value):
        raise ChainFeeError(name + " must be finite")
    return value


def geometric_closure(f, *, depth=DEPTH_DEFAULT):
    """Sum 1 + f + ... + f**depth (depth None = infinite). 0 < f < 1."""
    f = _finite_scale(f, "f")
    if not 0.0 < f < 1.0:
        raise ChainFeeError("per-level coefficient must lie in (0, 1)")
    if depth is not None:
        depth = int(depth)
        if depth < 0:
            raise ChainFeeError("depth must be nonnegative")
        total = (1.0 - f ** (depth + 1)) / (1.0 - f)
    else:
        total = 1.0 / (1.0 - f)
    return {
        "f": f,
        "depth": depth,
        "total_multiplier": total,
        "local_fraction": 1.0 / total,
        "exterior_fraction": 1.0 - 1.0 / total,
    }


def verdict(local_fraction):
    """Pre-registered verdict taxonomy. No fitting is performed anywhere."""
    value = _finite_scale(local_fraction, "local_fraction")
    if BAND_LO <= value <= BAND_HI:
        return NEAR_STATUS
    if 0.0 < value < 1.0:
        return NOT_NEAR_STATUS
    return BAND_STATUS


def scale_binding_fees(record):
    """Structural and recorded-transfer coefficients from scale-binding.json."""
    try:
        branch = record["development_branch"]
        omega = _finite_scale(branch["omega"], "omega")
        zeta = _finite_scale(branch["zeta"], "zeta")
        product = _finite_scale(
            record["recursive_compatibility"]["Omega_times_transfer"],
            "Omega_times_transfer",
        )
    except (KeyError, TypeError) as error:
        raise ChainFeeError("scale binding record is not the expected shape: " + str(error))
    transfer = product / omega
    return {
        "omega": omega,
        "zeta": zeta,
        "zeta_minus_omega_squared": zeta - omega * omega,
        "structural_one_over_omega": 1.0 / omega,
        "recorded_transfer": transfer,
        "endpoint_condition_omega_q_le_one": bool(omega * transfer <= 1.0 + 1e-12),
    }


def continuation_mass_fee(summary):
    """Measured one-level mass renormalisation from the neutral continuation."""
    try:
        final = summary["final"]
        measured = _finite_scale(final["coupled_inner_mass"], "coupled_inner_mass")
        fixed = _finite_scale(final["fixed_inner_mass"], "fixed_inner_mass")
    except (KeyError, TypeError) as error:
        raise ChainFeeError("continuation summary is not the expected shape: " + str(error))
    if fixed <= 0.0:
        raise ChainFeeError("fixed inner mass must be positive")
    fee = (measured - fixed) / fixed
    route_reference = 1.0 / ROUTE_RATIO ** 2
    return {
        "route_ratio": ROUTE_RATIO,
        "measured_relative_mass_renormalisation": fee,
        "one_over_ratio_squared": route_reference,
        "fee_over_one_over_ratio_squared": fee / route_reference,
    }


def candidate_table(fees):
    """Closure and verdict for every candidate coefficient in ``fees``."""
    rows = []
    for name, f in fees.items():
        closure = geometric_closure(f)
        status = verdict(closure["local_fraction"])
        rows.append({
            "candidate": name,
            "f": closure["f"],
            "total_multiplier": closure["total_multiplier"],
            "local_fraction": closure["local_fraction"],
            "exterior_fraction": closure["exterior_fraction"],
            "verdict": status,
            "band": [BAND_LO, BAND_HI],
            "target": TARGET_LOCAL_FRACTION,
        })
    return rows


def transfer_baryon_proximity(fees):
    """How close the recorded transfer sits to the baryon share. Reported,
    not fitted: the comparison is made only after both numbers are fixed."""
    transfer = float(fees["recorded_transfer"])
    return {
        "recorded_transfer": transfer,
        "baryon_target": TARGET_LOCAL_FRACTION,
        "absolute_offset": transfer - TARGET_LOCAL_FRACTION,
        "relative_offset": (transfer - TARGET_LOCAL_FRACTION) / TARGET_LOCAL_FRACTION,
        "numerology_caveat": (
            "proximity is a recorded observation about two fixed numbers; "
            "no coefficient was moved to produce it"
        ),
    }


def matrix_coefficient(f):
    """The R a shallow-power chain would need: f = R^-2 -> R = 1/sqrt(f)."""
    f = _finite_scale(f, "f")
    if not 0.0 < f < 1.0:
        raise ChainFeeError("matrix coefficient must lie in (0, 1)")
    return {"f": f, "implied_layer_ratio": 1.0 / math.sqrt(f)}


def required_power(ratio, *, local_fraction=TARGET_LOCAL_FRACTION):
    """Power law index nu such that local = 1 - R^-nu holds."""
    ratio = _finite_scale(ratio, "ratio")
    local_fraction = _finite_scale(local_fraction, "local_fraction")
    if ratio <= 1.0:
        raise ChainFeeError("layer ratio must exceed 1")
    f = 1.0 - local_fraction
    if not 0.0 < f < 1.0:
        raise ChainFeeError("local fraction must lie in (0, 1)")
    return {
        "ratio": ratio,
        "local_fraction": local_fraction,
        "per_level_decay": f,
        "required_nu": -math.log(f) / math.log(ratio),
    }


def band_required_power(ratio):
    """nu interval that puts 1 - R^-nu inside the pre-registered band."""
    nu_min = required_power(ratio, local_fraction=BAND_LO)["required_nu"]
    nu_max = required_power(ratio, local_fraction=BAND_HI)["required_nu"]
    return {"ratio": ratio, "nu_min": nu_min, "nu_max": nu_max}


def lll_screen_power(g4, flux):
    """Power from the charged lowest-Landau screen mu^2 = |q| g4^2 / (4 pi^2 r^2).

    e^{-int mu dr} = r^-nu with nu = g4 sqrt(|q|) / (2 pi).
    """
    g4 = _finite_scale(g4, "g4")
    flux = _finite_scale(flux, "flux")
    if g4 < 0.0 or flux <= 0.0:
        raise ChainFeeError("positive coupling and flux required")
    nu = g4 * math.sqrt(flux) / (2.0 * math.pi)
    return {"g4": g4, "flux": flux, "nu": nu}


def coupling_for_power(nu, flux):
    """g4 = 2 pi nu / sqrt(|q|): the coupling a demanded power implies."""
    nu = _finite_scale(nu, "nu")
    flux = _finite_scale(flux, "flux")
    if nu < 0.0 or flux <= 0.0:
        raise ChainFeeError("positive power and flux required")
    return {"nu": nu, "flux": flux, "g4": 2.0 * math.pi * nu / math.sqrt(flux)}


def chain_fraction(nu, ratio):
    """Closure under f = R^-nu, without raising the layer ratio."""
    nu = _finite_scale(nu, "nu")
    ratio = _finite_scale(ratio, "ratio")
    if nu <= 0.0 or ratio <= 1.0:
        raise ChainFeeError("positive power and ratio greater than 1 required")
    f = ratio ** (-nu)
    closure = geometric_closure(f)
    closure["power"] = nu
    closure["ratio"] = ratio
    closure["verdict"] = verdict(closure["local_fraction"])
    return closure


def spatial_fee(radius, weights, response, seam):
    """L1 fraction of a response that lands inside the inner level.

    Pre-registered fee for the two-level tangent probe: with a direction
    localized in the outer level, the fraction of the mass-function
    response whose absolute weight sits inside ``seam`` is the per-level
    propagation fee f. The L1 form keeps the fraction in [0, 1] even when
    the signed mass response changes sign across the seam.
    """
    radius = np.asarray(radius, dtype=float)
    weights = np.asarray(weights, dtype=float)
    response = np.asarray(response, dtype=float)
    seam = _finite_scale(seam, "seam")
    if radius.shape != weights.shape or radius.shape != response.shape:
        raise ChainFeeError("radius, weights, and response must share a shape")
    if seam <= float(np.min(radius)) or seam >= float(np.max(radius)):
        raise ChainFeeError("the seam must cut the grid")
    inside = radius < seam
    magnitude = np.abs(response) * weights
    total = float(np.sum(magnitude))
    if total <= 0.0:
        raise ChainFeeError("the response vanished on this grid")
    fraction = float(np.sum(magnitude[inside])) / total
    return {
        "seam": seam,
        "inner_level": [float(np.min(radius)), seam],
        "outer_level": [seam, float(np.max(radius))],
        "l1_mass_response_fraction_inside": fraction,
        "same_fraction_definition_of_f": True,
    }


def interval_shares(radius, weights, response, edges):
    """L1 share of a response inside each interval of an edge partition.

    The edges must lie inside the radial domain and increase strictly.
    Shares sum to one. This is the atlas cell the tangent instrument and
    the infall instrument are compared through.
    """
    radius = np.asarray(radius, dtype=float)
    weights = np.asarray(weights, dtype=float)
    response = np.asarray(response, dtype=float)
    edges = [float(item) for item in edges]
    if radius.shape != weights.shape or radius.shape != response.shape:
        raise ChainFeeError("radius, weights, and response must share a shape")
    if len(edges) < 2:
        raise ChainFeeError("a partition needs at least two edges")
    for item in edges:
        if not np.isfinite(item):
            raise ChainFeeError("partition edges must be finite")
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise ChainFeeError("partition edges must increase strictly")
    if edges[0] < float(np.min(radius)) - 1e-12 or edges[-1] > float(np.max(radius)) + 1e-12:
        raise ChainFeeError("partition edges exceed the radial domain")
    magnitude = np.abs(response) * weights
    total = float(np.sum(magnitude))
    if total <= 0.0:
        raise ChainFeeError("the response vanished on this grid")
    shares = []
    for index, (left, right) in enumerate(zip(edges, edges[1:])):
        if index == len(edges) - 2:
            inside = radius >= left
        else:
            inside = (radius >= left) & (radius < right)
        shares.append({
            "interval": [left, right],
            "l1_share": float(np.sum(magnitude[inside])) / total,
        })
    if abs(sum(item["l1_share"] for item in shares) - 1.0) > 1e-9:
        raise ChainFeeError("interval shares did not close to one")
    return shares