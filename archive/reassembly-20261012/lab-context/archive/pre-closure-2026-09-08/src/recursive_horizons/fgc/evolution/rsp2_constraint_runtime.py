"""Outcome-neutral late-event constraint assessment for RSP2.

RSP2 consumes the immutable CAL9 SSPRK3 states at 4,097 and 8,193 points and
one independently evolved 16,385-point SSPRK3 state at the same coordinate
time.  It applies the unchanged evolution-owned constraint definition,
roundoff enclosure, absolute guards, and strict 3/2 finest-pair order gate.

The module performs no I/O, advances no state, and authorizes no candidate
branch.  A passing target says only that the former radial-momentum miss was
pre-asymptotic on the tested ladder.  A failing target says only that it
persists through the tested 8,193 -> 16,385 pair.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Sequence

import numpy as np

from .cal4_common_event_diagnosis import (
    OwnedConstraintAdmission,
    owned_constraint_admission,
    owned_semidiscrete_constraint_snapshot,
)
from .numerical_engine import EvolutionState, UniformRadialGrid, array_content_sha256
from .proto4_admission import PROTO4_CONSTRAINT_ORDER
from .rsp2_checkpoint import (
    RSP2_CAL9_STAGE_COUNTS,
    RSP2_CAL9_STATE_HASHES,
    RSP2_CAL9_TARGET_TIME,
)


RSP2_POINT_COUNTS = (4097, 8193, 16385)
RSP2_TARGET_COMPONENT = "radial_momentum"
RSP2_METHOD = "SSPRK3"


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


@dataclass(frozen=True, slots=True)
class RSP2ConstraintAssessment:
    """Complete RSP2 target decision plus every unchanged constraint gate."""

    method: str
    coordinate_time: float
    point_counts: tuple[int, int, int]
    target_component: str
    target_minimum_order: float
    target_raw_owned_norms: tuple[float, float, float]
    target_accumulated_roundoff_enclosures: tuple[float, float, float]
    target_effective_norms_for_order_only: tuple[float, float, float]
    target_status: str
    target_finest_pair_order: float | None
    target_order_margin: float | None
    target_profile_monotone: bool
    target_classification_premises_passed: bool
    target_order_passed: bool
    target_preasymptotic_on_tested_ladder: bool
    target_persistent_through_tested_ladder: bool
    non_target_order_obstructions: tuple[str, ...]
    complete_constraint_admission_passed: bool
    classification: str
    constraint_admission: OwnedConstraintAdmission

    def __post_init__(self) -> None:
        if self.method != RSP2_METHOD:
            raise ValueError("RSP2 method differs")
        if self.point_counts != RSP2_POINT_COUNTS:
            raise ValueError("RSP2 point-count ladder differs")
        if self.target_component != RSP2_TARGET_COMPONENT:
            raise ValueError("RSP2 target component differs")
        object.__setattr__(
            self,
            "coordinate_time",
            _finite("coordinate_time", self.coordinate_time),
        )
        minimum = _finite(
            "target_minimum_order", self.target_minimum_order, positive=True
        )
        object.__setattr__(self, "target_minimum_order", minimum)
        for name in (
            "target_raw_owned_norms",
            "target_accumulated_roundoff_enclosures",
            "target_effective_norms_for_order_only",
        ):
            values = tuple(
                _finite(f"{name} value", value) for value in getattr(self, name)
            )
            if len(values) != 3 or any(value < 0.0 for value in values):
                raise ValueError(f"{name} differs")
            object.__setattr__(self, name, values)
        if self.target_finest_pair_order is not None:
            object.__setattr__(
                self,
                "target_finest_pair_order",
                _finite(
                    "target_finest_pair_order", self.target_finest_pair_order
                ),
            )
        if self.target_order_margin is not None:
            object.__setattr__(
                self,
                "target_order_margin",
                _finite("target_order_margin", self.target_order_margin),
            )
        for name in (
            "target_profile_monotone",
            "target_classification_premises_passed",
            "target_order_passed",
            "target_preasymptotic_on_tested_ladder",
            "target_persistent_through_tested_ladder",
            "complete_constraint_admission_passed",
        ):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")
        if (
            self.target_preasymptotic_on_tested_ladder
            and self.target_persistent_through_tested_ladder
        ):
            raise ValueError("RSP2 target cannot both clear and persist")
        if not isinstance(self.constraint_admission, OwnedConstraintAdmission):
            raise TypeError("constraint_admission must be OwnedConstraintAdmission")
        allowed_classifications = {
            "completed_target_and_complete_constraint_pass",
            "completed_target_cleared_other_constraint_obstruction",
            "completed_target_order_persistent_below_threshold",
            "completed_target_persistent_with_other_constraint_obstruction",
            "completed_target_inconclusive_premise_failure",
        }
        if self.classification not in allowed_classifications:
            raise ValueError("RSP2 endpoint classification differs")
        if self.target_status not in {
            "finite_order",
            "accumulated_roundoff_enclosed_zero_pair",
            "resolved_to_accumulated_roundoff_enclosure",
            "nonzero_reappeared_after_accumulated_roundoff_enclosure",
        }:
            raise ValueError("RSP2 target status differs")


def rsp2_constraint_assessment(
    states: Sequence[EvolutionState],
    grids: Sequence[UniformRadialGrid],
    *,
    accepted_stage_counts: Sequence[int],
    coordinate_time: Real,
    fixed_outer_rows: int = 4,
    coarsest_guard_maximum: Real = 0.1,
    finest_guard_maximum: Real = 0.005,
    minimum_finest_pair_order: Real = 1.5,
    planck_mass: Real = 2.0,
    scalar_mass: Real = 3.0,
    quartic_coupling: Real = 0.5,
    length_unit: Real = 4.0,
    enforce_frozen_predecessors: bool = True,
) -> RSP2ConstraintAssessment:
    """Evaluate the frozen RSP2 three-grid SSPRK3 constraint question."""

    state_records = tuple(states)
    grid_records = tuple(grids)
    stage_records = tuple(accepted_stage_counts)
    if (
        len(state_records) != 3
        or len(grid_records) != 3
        or len(stage_records) != 3
        or any(not isinstance(item, EvolutionState) for item in state_records)
        or any(not isinstance(item, UniformRadialGrid) for item in grid_records)
        or any(
            isinstance(item, bool) or not isinstance(item, int) or item < 0
            for item in stage_records
        )
    ):
        raise ValueError("RSP2 requires three states, grids, and stage counts")
    counts = tuple(item.point_count for item in grid_records)
    if counts != RSP2_POINT_COUNTS:
        raise ValueError("RSP2 requires the 4097 -> 8193 -> 16385 ladder")
    time = _finite("coordinate_time", coordinate_time)
    if np.float64(time).tobytes() != np.float64(RSP2_CAL9_TARGET_TIME).tobytes():
        raise ValueError("RSP2 requires the exact frozen t=23/16 endpoint")
    if not isinstance(enforce_frozen_predecessors, bool):
        raise TypeError("enforce_frozen_predecessors must be bool")
    if enforce_frozen_predecessors:
        if stage_records[:2] != RSP2_CAL9_STAGE_COUNTS:
            raise ValueError("RSP2 CAL9 predecessor stage counts differ")
        observed_hashes = tuple(
            array_content_sha256(state.u, state.p, state.q)
            for state in state_records[:2]
        )
        if observed_hashes != RSP2_CAL9_STATE_HASHES:
            raise ValueError("RSP2 CAL9 predecessor states differ")
    minimum = _finite(
        "minimum_finest_pair_order", minimum_finest_pair_order, positive=True
    )
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
        method=RSP2_METHOD,
        coarsest_guard_maximum=coarsest_guard_maximum,
        finest_guard_maximum=finest_guard_maximum,
        minimum_finest_pair_order=minimum,
    )
    target_index = PROTO4_CONSTRAINT_ORDER.index(RSP2_TARGET_COMPONENT)
    raw = tuple(
        float(item.owned_domain_norms.component_infinity[target_index])
        for item in snapshots
    )
    enclosures = tuple(item.accumulated_roundoff_enclosure for item in snapshots)
    effective = tuple(
        float(
            item.effective_owned_norms_for_order_only.component_infinity[
                target_index
            ]
        )
        for item in snapshots
    )
    monotone = all(
        right <= left for left, right in zip(effective[:-1], effective[1:], strict=True)
    )
    status = admission.component_status[RSP2_TARGET_COMPONENT]
    order = admission.component_finest_pair_orders[RSP2_TARGET_COMPONENT]
    margin = None if order is None else order - minimum
    classification_premises = (
        admission.common_event_alignment_passed
        and admission.coarsest_guard_passed
        and admission.finest_guard_passed
        and monotone
        and status == "finite_order"
        and order is not None
    )
    target_passed = bool(classification_premises and order >= minimum)
    target_persistent = bool(classification_premises and order < minimum)

    non_target: list[str] = []
    for index, name in enumerate(PROTO4_CONSTRAINT_ORDER):
        if name == RSP2_TARGET_COMPONENT:
            continue
        component_status = admission.component_status[name]
        component_order = admission.component_finest_pair_orders[name]
        component_effective = tuple(
            float(item.effective_owned_norms_for_order_only.component_infinity[index])
            for item in snapshots
        )
        if component_status == "nonzero_reappeared_after_accumulated_roundoff_enclosure":
            non_target.append(name)
        elif component_order is not None and component_order < minimum:
            non_target.append(name)
        elif any(
            right > left
            for left, right in zip(
                component_effective[:-1], component_effective[1:], strict=True
            )
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

    return RSP2ConstraintAssessment(
        method=RSP2_METHOD,
        coordinate_time=time,
        point_counts=RSP2_POINT_COUNTS,
        target_component=RSP2_TARGET_COMPONENT,
        target_minimum_order=minimum,
        target_raw_owned_norms=raw,
        target_accumulated_roundoff_enclosures=enclosures,
        target_effective_norms_for_order_only=effective,
        target_status=status,
        target_finest_pair_order=order,
        target_order_margin=margin,
        target_profile_monotone=monotone,
        target_classification_premises_passed=classification_premises,
        target_order_passed=target_passed,
        target_preasymptotic_on_tested_ladder=target_passed,
        target_persistent_through_tested_ladder=target_persistent,
        non_target_order_obstructions=tuple(non_target),
        complete_constraint_admission_passed=admission.admission_passed,
        classification=classification,
        constraint_admission=admission,
    )


__all__ = [
    "RSP2_METHOD",
    "RSP2_POINT_COUNTS",
    "RSP2_TARGET_COMPONENT",
    "RSP2ConstraintAssessment",
    "rsp2_constraint_assessment",
]
