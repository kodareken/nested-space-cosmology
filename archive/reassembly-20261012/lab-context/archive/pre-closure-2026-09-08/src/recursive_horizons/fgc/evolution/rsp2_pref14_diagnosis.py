"""Independent checkpoint reduction for the RSP2/PREF14 result.

The RSP2 runner used :mod:`rsp2_constraint_runtime` to classify its endpoint.
PREF14 deliberately does not call that classifier.  It restores the three
terminal states, evaluates the unchanged owned-domain constraints, reconstructs
both adjacent-grid orders, and independently assembles the serialized endpoint
record consumed by the post-run certificate.

This module performs no I/O and cannot authorize PROTO13 or a candidate run.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite, log
from numbers import Real
from typing import Any, Mapping, Sequence

import numpy as np

from .cal4_common_event_diagnosis import (
    OwnedConstraintAdmission,
    owned_constraint_admission,
    owned_semidiscrete_constraint_snapshot,
)
from .numerical_engine import EvolutionState, UniformRadialGrid, array_content_sha256
from .proto4_admission import PROTO4_CONSTRAINT_ORDER


PREF14_METHOD = "SSPRK3"
PREF14_POINT_COUNTS = (4097, 8193, 16385)
PREF14_STAGE_COUNTS = (2524, 3948, 7084)
PREF14_TARGET_COMPONENT = "radial_momentum"
PREF14_TARGET_TIME = 23.0 / 16.0


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _pair_order(medium: float, fine: float, ratio: float) -> tuple[str, float | None]:
    if medium == 0.0 and fine == 0.0:
        return "accumulated_roundoff_enclosed_zero_pair", None
    if medium > 0.0 and fine == 0.0:
        return "resolved_to_accumulated_roundoff_enclosure", None
    if medium == 0.0 and fine > 0.0:
        return "nonzero_reappeared_after_accumulated_roundoff_enclosure", None
    return "finite_order", log(medium / fine) / log(ratio)


@dataclass(frozen=True, slots=True)
class RSP2PREF14Diagnosis:
    """Independent endpoint reconstruction and its additional pair evidence."""

    state_sha256: tuple[str, str, str]
    target_adjacent_pair_orders: tuple[float, float]
    target_adjacent_pair_margins: tuple[float, float]
    target_previous_pair_failed: bool
    target_finest_pair_passed: bool
    target_preasymptotic_on_tested_ladder: bool
    complete_constraint_admission_passed: bool
    minimum_finite_finest_pair_order: float
    endpoint_assessment: Mapping[str, Any]
    constraint_admission: OwnedConstraintAdmission

    def __post_init__(self) -> None:
        if len(self.state_sha256) != 3 or any(
            not isinstance(value, str) or len(value) != 64
            for value in self.state_sha256
        ):
            raise ValueError("PREF14 state hashes differ")
        if len(self.target_adjacent_pair_orders) != 2 or len(
            self.target_adjacent_pair_margins
        ) != 2:
            raise ValueError("PREF14 adjacent-pair evidence differs")
        for name in (
            "target_adjacent_pair_orders",
            "target_adjacent_pair_margins",
        ):
            if any(not isfinite(float(value)) for value in getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        for name in (
            "target_previous_pair_failed",
            "target_finest_pair_passed",
            "target_preasymptotic_on_tested_ladder",
            "complete_constraint_admission_passed",
        ):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")
        minimum = _finite(
            "minimum_finite_finest_pair_order",
            self.minimum_finite_finest_pair_order,
        )
        object.__setattr__(self, "minimum_finite_finest_pair_order", minimum)
        if not isinstance(self.constraint_admission, OwnedConstraintAdmission):
            raise TypeError("constraint_admission must be OwnedConstraintAdmission")


def diagnose_rsp2_pref14_endpoint(
    states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    *,
    accepted_stage_counts: Sequence[int],
    expected_state_sha256: Sequence[str],
    coordinate_time: Real,
    target_component: str = PREF14_TARGET_COMPONENT,
    minimum_finest_pair_order: Real = 1.5,
    fixed_outer_rows: int = 4,
    coarsest_guard_maximum: Real = 0.1,
    finest_guard_maximum: Real = 0.005,
    planck_mass: Real = 2.0,
    scalar_mass: Real = 3.0,
    quartic_coupling: Real = 0.5,
    length_unit: Real = 4.0,
) -> RSP2PREF14Diagnosis:
    """Recompute RSP2 from state arrays without trusting the runner result."""

    state_records = tuple(states)
    grid_records = tuple(grids)
    stage_records = tuple(accepted_stage_counts)
    expected_hashes = tuple(expected_state_sha256)
    if (
        len(state_records) != 3
        or len(grid_records) != 3
        or len(stage_records) != 3
        or len(expected_hashes) != 3
        or any(not isinstance(item, EvolutionState) for item in state_records)
        or any(not isinstance(item, UniformRadialGrid) for item in grid_records)
        or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in stage_records
        )
    ):
        raise ValueError("PREF14 requires three states, grids, stage counts, and hashes")
    point_counts = tuple(grid.point_count for grid in grid_records)
    if point_counts != PREF14_POINT_COUNTS:
        raise ValueError("PREF14 point-count ladder differs")
    if stage_records != PREF14_STAGE_COUNTS:
        raise ValueError("PREF14 accepted-stage counts differ")
    if target_component != PREF14_TARGET_COMPONENT:
        raise ValueError("PREF14 target component differs")
    time = _finite("coordinate_time", coordinate_time)
    if np.float64(time).tobytes() != np.float64(PREF14_TARGET_TIME).tobytes():
        raise ValueError("PREF14 requires the exact t=23/16 endpoint")
    minimum = _finite(
        "minimum_finest_pair_order", minimum_finest_pair_order, positive=True
    )
    observed_hashes = tuple(
        array_content_sha256(state.u, state.p, state.q) for state in state_records
    )
    if observed_hashes != expected_hashes:
        raise ValueError("PREF14 terminal state hash differs")

    snapshots = tuple(
        owned_semidiscrete_constraint_snapshot(
            state,
            grid,
            coordinate_time=time,
            diagnostic_spatial_order=2,
            fixed_outer_rows=fixed_outer_rows,
            accepted_stage_count=stage_count,
            planck_mass=planck_mass,
            scalar_mass=scalar_mass,
            quartic_coupling=quartic_coupling,
            length_unit=length_unit,
        )
        for state, grid, stage_count in zip(
            state_records, grid_records, stage_records, strict=True
        )
    )
    admission = owned_constraint_admission(
        snapshots,
        method=PREF14_METHOD,
        coarsest_guard_maximum=coarsest_guard_maximum,
        finest_guard_maximum=finest_guard_maximum,
        minimum_finest_pair_order=minimum,
    )

    component_effective: dict[str, tuple[float, float, float]] = {}
    component_adjacent: dict[str, tuple[tuple[str, float | None], ...]] = {}
    ratio = float(
        (point_counts[1] - 1) / (point_counts[0] - 1)
    )
    for index, name in enumerate(PROTO4_CONSTRAINT_ORDER):
        effective = tuple(
            float(
                snapshot.effective_owned_norms_for_order_only.component_infinity[
                    index
                ]
            )
            for snapshot in snapshots
        )
        component_effective[name] = effective
        component_adjacent[name] = tuple(
            _pair_order(left, right, ratio)
            for left, right in zip(effective[:-1], effective[1:], strict=True)
        )

    target_index = PROTO4_CONSTRAINT_ORDER.index(PREF14_TARGET_COMPONENT)
    target_raw = tuple(
        float(snapshot.owned_domain_norms.component_infinity[target_index])
        for snapshot in snapshots
    )
    target_enclosures = tuple(
        float(snapshot.accumulated_roundoff_enclosure) for snapshot in snapshots
    )
    target_effective = component_effective[PREF14_TARGET_COMPONENT]
    target_pairs = component_adjacent[PREF14_TARGET_COMPONENT]
    if any(status != "finite_order" or order is None for status, order in target_pairs):
        raise ValueError("PREF14 target adjacent pairs are not finite-order pairs")
    target_orders = tuple(float(order) for _status, order in target_pairs if order is not None)
    if len(target_orders) != 2:
        raise ValueError("PREF14 target adjacent orders differ")
    monotone = all(
        right <= left
        for left, right in zip(target_effective[:-1], target_effective[1:], strict=True)
    )
    target_status = admission.component_status[PREF14_TARGET_COMPONENT]
    target_finest_order = admission.component_finest_pair_orders[
        PREF14_TARGET_COMPONENT
    ]
    classification_premises = (
        admission.common_event_alignment_passed
        and admission.coarsest_guard_passed
        and admission.finest_guard_passed
        and monotone
        and target_status == "finite_order"
        and target_finest_order is not None
    )
    target_passed = bool(
        classification_premises and target_finest_order is not None
        and target_finest_order >= minimum
    )
    target_persistent = bool(
        classification_premises and target_finest_order is not None
        and target_finest_order < minimum
    )

    non_target: list[str] = []
    for name in PROTO4_CONSTRAINT_ORDER:
        if name == PREF14_TARGET_COMPONENT:
            continue
        status = admission.component_status[name]
        order = admission.component_finest_pair_orders[name]
        values = component_effective[name]
        if status == "nonzero_reappeared_after_accumulated_roundoff_enclosure":
            non_target.append(name)
        elif order is not None and order < minimum:
            non_target.append(name)
        elif any(
            right > left for left, right in zip(values[:-1], values[1:], strict=True)
        ):
            non_target.append(name)

    if not classification_premises:
        classification = "completed_target_inconclusive_premise_failure"
    elif target_passed and admission.admission_passed:
        classification = "completed_target_and_complete_constraint_pass"
    elif target_passed:
        classification = "completed_target_cleared_other_constraint_obstruction"
    elif non_target:
        classification = "completed_target_persistent_with_other_constraint_obstruction"
    else:
        classification = "completed_target_order_persistent_below_threshold"

    endpoint = {
        "method": PREF14_METHOD,
        "coordinate_time": time,
        "point_counts": list(point_counts),
        "target_component": PREF14_TARGET_COMPONENT,
        "target_minimum_order": minimum,
        "target_raw_owned_norms": list(target_raw),
        "target_accumulated_roundoff_enclosures": list(target_enclosures),
        "target_effective_norms_for_order_only": list(target_effective),
        "target_status": target_status,
        "target_finest_pair_order": target_finest_order,
        "target_order_margin": (
            None if target_finest_order is None else target_finest_order - minimum
        ),
        "target_profile_monotone": monotone,
        "target_classification_premises_passed": classification_premises,
        "target_order_passed": target_passed,
        "target_preasymptotic_on_tested_ladder": target_passed,
        "target_persistent_through_tested_ladder": target_persistent,
        "non_target_order_obstructions": non_target,
        "complete_constraint_admission_passed": admission.admission_passed,
        "classification": classification,
        "constraint_admission": asdict(admission),
    }
    finite_minimum = admission.minimum_finite_finest_pair_order
    if finite_minimum is None:
        raise ValueError("PREF14 finite-order minimum is absent")
    return RSP2PREF14Diagnosis(
        state_sha256=observed_hashes,
        target_adjacent_pair_orders=(target_orders[0], target_orders[1]),
        target_adjacent_pair_margins=(
            target_orders[0] - minimum,
            target_orders[1] - minimum,
        ),
        target_previous_pair_failed=target_orders[0] < minimum,
        target_finest_pair_passed=target_orders[1] >= minimum,
        target_preasymptotic_on_tested_ladder=(
            target_orders[0] < minimum <= target_orders[1]
        ),
        complete_constraint_admission_passed=admission.admission_passed,
        minimum_finite_finest_pair_order=finite_minimum,
        endpoint_assessment=endpoint,
        constraint_admission=admission,
    )


__all__ = [
    "PREF14_METHOD",
    "PREF14_POINT_COUNTS",
    "PREF14_STAGE_COUNTS",
    "PREF14_TARGET_COMPONENT",
    "PREF14_TARGET_TIME",
    "RSP2PREF14Diagnosis",
    "diagnose_rsp2_pref14_endpoint",
]
