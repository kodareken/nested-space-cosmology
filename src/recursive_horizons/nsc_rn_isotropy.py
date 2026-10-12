"""The isotropy trajectory: how far a recorded state sits from the isotropic point.

One finite question in one named domain. On the recorded states of the
Reissner-Nordstrom programme, the two forces that hold a region are the
compressive radial contribution and the resisting angular one. Their balance

    B(r) = |F_in(r)| / (|F_in(r)| + |F_out(r)|)

is 1/2 exactly when the two acts equally, which is the isotropic point of
that pair; B -> 1 means the radial force acts alone. This module reduces the
already-written balance profiles of the family and of the inward-march
campaign to a single scalar trajectory per case,

    I = <|B(r) - 1/2|>  over the recorded cells,

the departure from the isotropic point, and it registers the pass/fail
condition of the relaxation test that has not been run yet. The registration
is written before the run exists, so the test cannot be reinterpreted after
the fact.

Domain. The areal radius on the recorded grids, both the family layout and
the inward-march layout. The balance is read from the recorded pressures; no
time step is taken.

Two honest limits. First, I measures the *matter* balance, not the metric's
anisotropy; the metric reading is the tidal comparison in the discovery
family, and this module does not replace it. Second, the recorded cases are
driven: the shell is still falling through every recorded station. A decaying
or growing I here is therefore a statement about the driven phase only, which
is exactly why the relaxation test needs its own run and its own registered
thresholds.
"""
from __future__ import annotations

from typing import Any, Dict

import numpy as np
from numba import njit

API_VERSION = "nsc-rn-isotropy-v1"

#: A station whose balance profile never leaves this band around the
#: isotropic point is treated as isotropic within the recorded resolution.
ISOTROPIC_BAND = 0.05
#: Registered before the relaxation run: the tail anisotropy must fall below
#: this fraction of the driven-phase maximum to count as decay.
RELAXATION_DECAY_FRACTION = 0.5
#: Registered before the relaxation run: movement below this multiple of the
#: measured refinement indicator is not a reading.
REFINEMENT_MULTIPLE = 3.0


@njit(cache=True)
def departure_profile(balance: np.ndarray) -> np.ndarray:
    """|B - 1/2| on every recorded cell."""
    out = np.empty(balance.shape[0])
    for i in range(balance.shape[0]):
        out[i] = abs(balance[i] - 0.5)
    return out


@njit(cache=True)
def balance_from_pressures(radial: np.ndarray, angular: np.ndarray) -> np.ndarray:
    """The same pointwise balance as the owning module: |F_in|/(|F_in|+|F_out|).

    This mirrors ``nsc_rn_gradient_bubble._balance_profile`` exactly, so a
    reduction computed from pressures that were not written by the march (for
    example the frozen arm's, recomputed from its own field) uses the same
    statistic as the recorded one. A test pins the two against each other on a
    recorded station.
    """
    n = radial.shape[0]
    balance = np.zeros(n)
    for i in range(n):
        a = abs(radial[i])
        b = abs(angular[i])
        denom = a + b
        if denom > 0.0:
            balance[i] = a / denom
    return balance


def isotropy_index(balance: np.ndarray) -> float:
    """Mean departure from the isotropic point over the recorded cells."""
    profile = departure_profile(np.ascontiguousarray(np.asarray(balance, dtype=np.float64)))
    finite = profile[np.isfinite(profile)]
    if finite.size == 0:
        raise ValueError("the balance profile has no finite cells")
    return float(np.mean(finite))


def station_row(station_name: str, t_over_rm: float, balance: np.ndarray,
                refinement: float | None = None) -> Dict[str, Any]:
    index = isotropy_index(balance)
    return {
        "station": station_name,
        "t_over_rm": float(t_over_rm),
        "isotropy_index": index,
        "isotropic_within_band": bool(index <= ISOTROPIC_BAND),
        "refinement_indicator": None if refinement is None else float(refinement),
    }


def relaxation_verdict(rows: list[Dict[str, Any]]) -> Dict[str, Any]:
    """Registered pass/fail for the not-yet-run relaxation test.

    Driven-phase growth is expected and allowed. The test reads the tail: the
    late-time index must fall below RELAXATION_DECAY_FRACTION of the maximum
    departure seen in the driven phase, and the fall must exceed
    REFINEMENT_MULTIPLE times the largest recorded refinement indicator.
    Anything else is `inconclusive_on_this_domain`, which is a declared
    outcome rather than a failure to report.
    """
    if not rows:
        raise ValueError("the relaxation verdict needs at least one station")
    indices = np.array([row["isotropy_index"] for row in rows], dtype=float)
    times = np.array([row["t_over_rm"] for row in rows], dtype=float)
    order = np.argsort(times)
    indices = indices[order]
    times = times[order]
    driven_max = float(np.max(indices[:-1])) if indices.size > 1 else float(indices[0])
    tail = float(indices[-1])
    refinements = [row["refinement_indicator"] for row in rows
                   if row.get("refinement_indicator") is not None]
    refinement = float(max(refinements)) if refinements else 0.0
    fell_enough = tail <= RELAXATION_DECAY_FRACTION * driven_max
    exceeds_refinement = (driven_max - tail) > REFINEMENT_MULTIPLE * max(refinement, 1e-12)
    if fell_enough and exceeds_refinement:
        verdict = "relaxation_observed"
    elif not fell_enough:
        verdict = "no_relaxation_observed"
    else:
        verdict = "inconclusive_on_this_domain"
    return {
        "verdict": verdict,
        "driven_phase_max_departure": driven_max,
        "tail_departure": tail,
        "decay_fraction_registered": RELAXATION_DECAY_FRACTION,
        "refinement_multiple_registered": REFINEMENT_MULTIPLE,
        "largest_recorded_refinement": refinement,
        "tail_time": float(times[-1]),
        "registration": (
            "written before the relaxation run: pass requires the tail index at "
            "or below the registered fraction of the driven maximum AND a fall "
            "larger than the registered multiple of the refinement indicator; "
            "anything else is inconclusive on this domain"
        ),
    }