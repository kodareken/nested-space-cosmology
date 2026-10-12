"""Diagnose PROTO10's centre-adjacent binary64 source floor.

This module reads no candidate branch and advances no trajectory.  It reduces
the completed GR-0 CAL7 event ledger and reconstructs the unchanged initial
states in two numerical representations:

* the frozen PROTO10 representation, which differentiates the full ADM fields;
* a prospective well-balanced control, which differentiates deviations from
  the exact Minkowski spherical reference and restores its derivative
  analytically.

The second route is a diagnostic control only.  It is not silently substituted
into PROTO10 and does not authorize another campaign.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real
from typing import Any, Mapping, Sequence

import numpy as np

from .calibration_runtime import project_gr0_semidiscrete_state
from .gr0_calibration import GR0GridInitialData, Q_CENTER_PARITIES
from .gr0_direct_source import solve_gr0_grid_accelerations
from .health_monitor import FIELD_ORDER
from .numerical_engine import EvolutionState, SBPFirstDerivative
from .vectorized_source import ADM_CENTER_PARITIES


SOURCE_FAILURE = "newton_residual_limit"
EXPECTED_AMPLITUDES = ("5/2", "3")
EXPECTED_POINT_COUNTS = (2049, 4097, 8193)
EXPECTED_TERMINAL_MEMBER = "RK4-8193"
MINKOWSKI_U_ORDER = ("alpha", "shift", "lambda", "R", "phi", "chi")


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
class InitialSourceFloorRecord:
    amplitude: str
    method: str
    spatial_order: int
    point_count: int
    representation: str
    residual_infinity: float
    maximum_grid_index: int
    maximum_coordinate_radius: float
    maximum_field: str
    signed_maximum_residual: float
    refinement_iterations: int
    kinetic_condition_infinity: float
    roundoff_coefficient: float

    def __post_init__(self) -> None:
        if self.amplitude not in EXPECTED_AMPLITUDES:
            raise ValueError("source-floor amplitude differs")
        if self.method not in {"RK4", "SSPRK3"}:
            raise ValueError("source-floor method differs")
        if self.spatial_order not in {2, 4}:
            raise ValueError("source-floor spatial order differs")
        if self.point_count not in EXPECTED_POINT_COUNTS:
            raise ValueError("source-floor point count differs")
        if self.representation not in {"PROTO10_raw", "well_balanced_control"}:
            raise ValueError("source-floor representation differs")
        if not 1 <= self.maximum_grid_index < self.point_count:
            raise ValueError("source-floor maximum index differs")
        if self.maximum_field not in FIELD_ORDER:
            raise ValueError("source-floor maximum field differs")
        for name in (
            "residual_infinity",
            "maximum_coordinate_radius",
            "kinetic_condition_infinity",
            "roundoff_coefficient",
        ):
            object.__setattr__(
                self,
                name,
                _finite(name, getattr(self, name), positive=name != "roundoff_coefficient"),
            )
        object.__setattr__(
            self,
            "signed_maximum_residual",
            _finite("signed_maximum_residual", self.signed_maximum_residual),
        )
        if self.roundoff_coefficient < 0.0:
            raise ValueError("roundoff coefficient must be nonnegative")
        if not 0 <= self.refinement_iterations <= 16:
            raise ValueError("source-floor refinement count differs")

    def as_record(self) -> dict[str, Any]:
        return {
            name: getattr(self, name)
            for name in self.__dataclass_fields__
        }


def minkowski_reference_profiles(coordinates: object) -> tuple[np.ndarray, np.ndarray]:
    """Return exact binary64 ADM reference values and their radial derivative."""

    radius = np.asarray(coordinates, dtype=np.float64)
    if radius.ndim != 1 or radius.size < 9 or not np.all(np.isfinite(radius)):
        raise ValueError("reference coordinates must be a finite radial vector")
    if radius[0] != 0.0 or np.any(np.diff(radius) <= 0.0):
        raise ValueError("reference coordinates must increase from the centre")
    values = np.zeros((radius.size, len(FIELD_ORDER)), dtype=np.float64)
    derivative = np.zeros_like(values)
    values[:, 0] = 1.0
    values[:, 2] = 1.0
    values[:, 3] = radius
    derivative[:, 3] = 1.0
    return values, derivative


def well_balanced_projected_state(
    initial: GR0GridInitialData,
    *,
    spatial_order: int,
) -> EvolutionState:
    """Project ``q`` by differentiating only departure from Minkowski."""

    if not isinstance(initial, GR0GridInitialData):
        raise TypeError("initial must be GR0GridInitialData")
    reference_u, reference_q = minkowski_reference_profiles(
        initial.grid.coordinates
    )
    derivative = SBPFirstDerivative(initial.grid, spatial_order)
    q = derivative.differentiate(
        initial.state.u - reference_u,
        center_parities=ADM_CENTER_PARITIES,
    ) + reference_q
    return EvolutionState(initial.state.u, initial.state.p, q)


def initial_source_floor_record(
    initial: GR0GridInitialData,
    *,
    amplitude: str,
    method: str,
    spatial_order: int,
    representation: str,
    residual_tolerance: Real,
    condition_number_maximum: Real,
) -> InitialSourceFloorRecord:
    """Locate the largest initial GR-0 source residual on one grid."""

    if representation == "PROTO10_raw":
        state = project_gr0_semidiscrete_state(
            initial,
            spatial_order=spatial_order,
        )
        reference_q = np.zeros_like(state.q)
    elif representation == "well_balanced_control":
        state = well_balanced_projected_state(
            initial,
            spatial_order=spatial_order,
        )
        _reference_u, reference_q = minkowski_reference_profiles(
            initial.grid.coordinates
        )
    else:
        raise ValueError("unknown source-floor representation")
    derivative = SBPFirstDerivative(initial.grid, spatial_order)
    p_r = derivative.differentiate(
        state.p,
        center_parities=ADM_CENTER_PARITIES,
    )
    if representation == "well_balanced_control":
        q_r = derivative.differentiate(
            state.q - reference_q,
            center_parities=Q_CENTER_PARITIES,
        )
    else:
        q_r = derivative.differentiate(
            state.q,
            center_parities=Q_CENTER_PARITIES,
        )
    source = solve_gr0_grid_accelerations(
        state.u[1:],
        state.p[1:],
        state.q[1:],
        p_r[1:],
        q_r[1:],
        initial.grid.coordinates[1:],
        residual_tolerance=residual_tolerance,
        condition_number_maximum=condition_number_maximum,
    )
    flattened = int(np.argmax(np.abs(source.residuals)))
    row, field = np.unravel_index(flattened, source.residuals.shape)
    grid_index = int(row + 1)
    radius = float(initial.grid.coordinates[grid_index])
    residual = float(source.residuals[row, field])
    coefficient = (
        abs(residual) * radius**2 / np.finfo(np.float64).eps
    )
    return InitialSourceFloorRecord(
        amplitude=amplitude,
        method=method,
        spatial_order=spatial_order,
        point_count=initial.grid.point_count,
        representation=representation,
        residual_infinity=source.residual_infinity,
        maximum_grid_index=grid_index,
        maximum_coordinate_radius=radius,
        maximum_field=FIELD_ORDER[field],
        signed_maximum_residual=residual,
        refinement_iterations=source.refinement_iterations,
        kinetic_condition_infinity=source.kinetic_condition_infinity_maximum,
        roundoff_coefficient=coefficient,
    )


def _rejection_trace(
    events: Sequence[Mapping[str, Any]],
    amplitude: str,
) -> list[Mapping[str, Any]]:
    return [
        event
        for event in events
        if event.get("event_type") == "rejected_unaccepted_source_only_proposal"
        and event.get("amplitude") == amplitude
    ]


def diagnose_campaign_pair(
    events: Sequence[Mapping[str, Any]],
    campaign_result: Mapping[str, Any],
) -> dict[str, Any]:
    """Reduce the two amplitude stops without changing their classification."""

    if campaign_result.get("classification") != "calibration_failed_no_eligible_GR0_case":
        raise ValueError("PROTO10 campaign classification differs")
    amplitude_records = campaign_result.get("amplitude_records")
    if not isinstance(amplitude_records, list) or [
        item.get("amplitude") for item in amplitude_records
    ] != list(EXPECTED_AMPLITUDES):
        raise ValueError("PROTO10 amplitude record order differs")
    common = [event for event in events if event.get("event_type") == "common_event"]
    if [event.get("amplitude") for event in common] != list(EXPECTED_AMPLITUDES):
        raise ValueError("PROTO10 common-event order differs")
    for event in common:
        assessment = event.get("assessment", {})
        if (
            event.get("event_index") != 0
            or assessment.get("coordinate_time") != 0.0
            or assessment.get("qualified_trapped_common_event") is not False
            or assessment.get("primary_common_event", {}).get("admission_passed") is not True
            or assessment.get("comparator_common_event", {}).get("admission_passed") is not True
            or assessment.get("primary_raw_PROTO9_common_event", {}).get("admission_passed") is not False
            or assessment.get("comparator_raw_PROTO9_common_event", {}).get("admission_passed") is not False
        ):
            raise ValueError("PROTO10 initial common-event contract differs")

    traces = {amplitude: _rejection_trace(events, amplitude) for amplitude in EXPECTED_AMPLITUDES}
    if any(len(trace) != 131 for trace in traces.values()):
        raise ValueError("PROTO10 per-amplitude rejection count differs")
    for amplitude, trace in traces.items():
        for event in trace:
            if (
                event.get("retryable_source_only") is not True
                or event.get("non_source_failure_vetoed_retry") is not False
                or event.get("complete_failure_set") != [SOURCE_FAILURE]
                or event.get("fields_bitwise_preserved") is not True
                or event.get("external_transaction_state_preserved") is not True
                or event.get("time_advanced") is not False
                or event.get("accepted_stage_count_advanced") is not False
                or event.get("causal_debit_advanced") is not False
                or event.get("preaccept_called") is not False
            ):
                raise ValueError(f"PROTO10 rejection semantics differ for {amplitude}")
    paired_maximum_time_difference = 0.0
    paired_maximum_step_relative_difference = 0.0
    paired_maximum_residual_relative_difference = 0.0
    for left, right in zip(traces["5/2"], traces["3"], strict=True):
        if (
            left.get("member") != right.get("member")
            or left.get("retry_count_for_accepted_step")
            != right.get("retry_count_for_accepted_step")
        ):
            raise ValueError("PROTO10 paired rejection identity differs")
        paired_maximum_time_difference = max(
            paired_maximum_time_difference,
            abs(float(left["initial_time"]) - float(right["initial_time"])),
        )
        left_step = float(left["attempted_step_size"])
        right_step = float(right["attempted_step_size"])
        paired_maximum_step_relative_difference = max(
            paired_maximum_step_relative_difference,
            abs(left_step - right_step) / max(abs(left_step), abs(right_step)),
        )
        for left_failure, right_failure in zip(
            left.get("failed_evaluations", ()),
            right.get("failed_evaluations", ()),
            strict=True,
        ):
            left_residual = float(left_failure["source_residual_infinity"])
            right_residual = float(right_failure["source_residual_infinity"])
            paired_maximum_residual_relative_difference = max(
                paired_maximum_residual_relative_difference,
                abs(left_residual - right_residual)
                / max(abs(left_residual), abs(right_residual)),
            )

    terminal: dict[str, Any] = {}
    for item in amplitude_records:
        amplitude = item["amplitude"]
        stop = item.get("stop", {})
        evidence = stop.get("complete_failure_evidence", {})
        rejected = evidence.get("last_rejected_proposal", {})
        if (
            item.get("eligible") is not False
            or item.get("last_completed_common_event_index") != 0
            or item.get("last_completed_coordinate_time") != 0.0
            or stop.get("classification") != "scientific_source_retry_exhausted"
            or stop.get("reason") != SOURCE_FAILURE
            or stop.get("member") != EXPECTED_TERMINAL_MEMBER
            or stop.get("method") != "RK4"
            or stop.get("point_count") != 8193
            or stop.get("target_common_event_index") != 1
            or stop.get("cross_member_rollback_verified") is not True
            or stop.get("cross_member_state_rolled_back") is not True
            or evidence.get("exhaustion_kind") != "minimum_step_size"
            or rejected.get("retry_count_for_accepted_step") != 20
            or rejected.get("cumulative_member_source_retry_count") != 129
        ):
            raise ValueError(f"PROTO10 terminal stop differs for {amplitude}")
        terminal[amplitude] = {
            "last_accepted_member_time": rejected["initial_time"],
            "last_rejected_step_size": rejected["attempted_step_size"],
            "next_forbidden_step_size": evidence["attempted_step_size"],
            "minimum_step_size": evidence["minimum_step_size"],
            "last_candidate_residual": rejected["failed_evaluations"][-1][
                "source_residual_infinity"
            ],
            "accepted_step_index": rejected["accepted_step_index"],
            "accepted_state_sha256": rejected["accepted_state_sha256"],
        }
    terminal_time_difference = abs(
        float(terminal["5/2"]["last_accepted_member_time"])
        - float(terminal["3"]["last_accepted_member_time"])
    )
    terminal_residual_relative_difference = abs(
        float(terminal["5/2"]["last_candidate_residual"])
        - float(terminal["3"]["last_candidate_residual"])
    ) / max(
        float(terminal["5/2"]["last_candidate_residual"]),
        float(terminal["3"]["last_candidate_residual"]),
    )
    return {
        "event_count": len(events),
        "common_event_count": len(common),
        "source_only_rejection_count": sum(len(value) for value in traces.values()),
        "per_amplitude_source_only_rejection_count": {
            amplitude: len(value) for amplitude, value in traces.items()
        },
        "all_rejections_exactly_rolled_back": True,
        "no_non_source_failure_observed": True,
        "terminal": terminal,
        "terminal_time_absolute_difference": terminal_time_difference,
        "terminal_residual_relative_difference": terminal_residual_relative_difference,
        "paired_maximum_time_absolute_difference": paired_maximum_time_difference,
        "paired_maximum_step_relative_difference": paired_maximum_step_relative_difference,
        "paired_maximum_residual_relative_difference": paired_maximum_residual_relative_difference,
    }


__all__ = [
    "EXPECTED_AMPLITUDES",
    "EXPECTED_POINT_COUNTS",
    "EXPECTED_TERMINAL_MEMBER",
    "InitialSourceFloorRecord",
    "diagnose_campaign_pair",
    "initial_source_floor_record",
    "minkowski_reference_profiles",
    "well_balanced_projected_state",
]
