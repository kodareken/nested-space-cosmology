"""Build a finite computational cap profile from supplied zero-case wall jets.

The caller supplies values and material rates evaluated at the same saved
stages. This adapter does not prepare a source, reconstruct a physical stage,
or evolve a field. Its second jets are interpolation curvature; they are not
an independently measured or computed physical second derivative.
"""
from collections.abc import Mapping
import math

import numpy as np

from .nsc_atlas_characteristic_background import CharacteristicProfile, MAX_POINTS


CURVATURE_METHOD = "numpy.gradient(F_prime, arrival_labels, edge_order=2)"
SCOPE = ("computational C2 quintic background from supplied zero-case wall values "
         "and same-stage first jets; finite boundary discrepancies remain evolved; "
         "interpolation curvature is not physical second-jet evidence")


def _rows(rows, name):
    if not isinstance(rows, Mapping) or not {"T", "r", "c"}.issubset(rows):
        raise ValueError(name+" requires T,r,c arrays")
    result = {}
    for key in ("T", "r", "c"):
        value = rows[key]
        if isinstance(value, (list, tuple)) and any(
                isinstance(v, (bool, np.bool_)) for v in np.asarray(value, dtype=object).flat):
            raise ValueError(name+" contains boolean data")
        raw = np.asarray(value)
        if (raw.ndim != 1 or raw.dtype.kind not in "iuf"
                or not 3 <= raw.size <= MAX_POINTS):
            raise ValueError(name+" requires 3..513 real one-dimensional rows")
        result[key] = np.array(raw, dtype=float, copy=True)
        if not np.all(np.isfinite(result[key])):
            raise ValueError(name+" must be finite")
    if len({a.shape for a in result.values()}) != 1:
        raise ValueError(name+" array shapes differ")
    return result


def build_background_profile(walls, wall_rates, *, alpha, budget):
    """Construct one bounded profile from wall arrays with matching row order.

    Both mappings contain ``T,r,c`` vectors. ``wall_rates`` are derivatives
    with respect to the same common evolution parameter, not derivatives of
    a saved-value interpolant. At alpha=0, before any reflected return,
    zeta=T+atanh(r), F=c and F'=c_x/zeta_x. Second jets are the declared
    second-order gradient of these supplied first jets in zeta. A single
    quintic Hermite interpolant then owns every runtime derivative.

    The first value and both resulting derivatives must be exactly quiet.
    Nothing is extrapolated beyond the last label. Later background queries
    must satisfy T+atanh(r)<=labels[-1]; the returned profile enforces this.
    Lower extension is only the certified zero branch. The caller retains
    responsibility for authenticating the supplied case and its same-stage
    provenance; numerical arrays alone cannot establish that provenance.

    Charges are six input scalar arrays, three derived arrays and one knot
    jet evaluation per row, all on the caller's shared budget.
    """
    if (isinstance(alpha, (bool, np.bool_)) or not np.isscalar(alpha)
            or not isinstance(alpha, (int, float, np.integer, np.floating))
            or not math.isfinite(float(alpha)) or float(alpha) != 0.):
        raise ValueError("background input must be the alpha=0 case")
    if not callable(getattr(budget, "check", None)):
        raise ValueError("one shared cumulative budget.check(count) is required")
    budget.check()
    before = budget.evaluations if hasattr(budget, "evaluations") else None
    walls, rates = _rows(walls, "walls"), _rows(wall_rates, "wall_rates")
    if walls["T"].shape != rates["T"].shape:
        raise ValueError("wall values and rates must have the same row count")
    n = len(walls["T"])
    budget.check(6*n)
    r = walls["r"]
    if np.any(r <= 0) or np.any(r >= 1):
        raise ValueError("wall radii must lie strictly between zero and one")
    with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
        try:
            labels = walls["T"]+np.arctanh(r)
            arrival_rate = rates["T"]+rates["r"]/(1-r*r)
            if (not np.all(np.isfinite(labels)) or not np.all(np.diff(labels) > 0)
                    or not np.all(np.isfinite(arrival_rate)) or np.any(arrival_rate <= 0)):
                raise ValueError("arrival labels must increase and their rates must be positive")
            first = rates["c"]/arrival_rate
            second = np.gradient(first, labels, edge_order=2)
        except (FloatingPointError, OverflowError) as error:
            raise ValueError("nonfinite profile chain or computational curvature") from error
    budget.check(3*n)
    values = walls["c"].copy()
    if not all(np.all(np.isfinite(a)) for a in (first, second)):
        raise ValueError("nonfinite profile derivatives")
    if values[0] != 0 or first[0] != 0 or second[0] != 0:
        raise ValueError("first value, first jet and computational curvature must be exactly quiet")
    no_return_margin = float(labels[0]-np.max(walls["T"]))
    if no_return_margin <= 0:
        raise ValueError("input must remain before the first central arrival")
    profile = CharacteristicProfile(labels, values, first, second)
    knots = profile.jet(labels, budget=budget)
    diagnostics = {
        "no_return_margin": no_return_margin,
        "minimum_arrival_rate": float(np.min(arrival_rate)),
        "minimum_label_spacing": float(np.min(np.diff(labels))),
        "max_abs_material_chain_defect": float(np.max(abs(first*arrival_rate-rates["c"]))),
        "max_abs_knot_value_defect": float(np.max(abs(knots["F"]-values))),
        "max_abs_knot_first_defect": float(np.max(abs(knots["F_prime"]-first))),
        "max_abs_knot_second_defect": float(np.max(abs(knots["F_second"]-second))),
        "quiet_anchor": {"value": float(values[0]), "first": float(first[0]),
                         "second": float(second[0])},
    }
    budget.check()
    return {
        "profile": profile, "labels": labels, "values": values, "first": first,
        "second": second, "walls": walls, "wall_rates": rates,
        "arrival_rate": arrival_rate, "label_range": [float(labels[0]), float(labels[-1])],
        "alpha": 0., "curvature_method": CURVATURE_METHOD, "scope": SCOPE,
        "diagnostics": diagnostics,
        "charged_evaluations": None if before is None else budget.evaluations-before,
    }
