"""Outcome-neutral reduction of the frozen RSP1 resolution study.

The reducer is deliberately narrower than the generic PROTO11 admission
machinery.  It asks only the question prospectively frozen by
``FGC-1-RSP1-FRZ1``: whether both direct adjacent ``phi/Lambda`` derivative
tail ratios fall strictly below one quarter at the single evolved endpoint,
with the source, constraint, absolute-budget, and complete-profile premises
still intact for both independent numerical methods.

The generic all-field raw nested-tail diagnostic remains public.  It is not a
premise of this targeted study and is never relabelled as passing here.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence


EXPECTED_METHODS = ("RK4", "SSPRK3")
EXPECTED_POINT_COUNTS = (4097, 8193, 16385)
EXPECTED_MEMBER_KEYS = tuple(
    f"{method}-{point_count}"
    for method in EXPECTED_METHODS
    for point_count in EXPECTED_POINT_COUNTS
)
EXPECTED_CLASSIFICATION = "completed_direct_amplitude_three_spectrum_pass"
EXPECTED_RUNNER_ID = "FGC-1-RSP1-RUN1-RUNNER"
EXPECTED_TARGET_FIELD = "phi_over_Lambda"


def _require_mapping(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _require_bool(name: str, value: Any, expected: bool) -> None:
    if value is not expected:
        raise ValueError(f"{name} differs")


def _finite_pair(name: str, value: Any) -> list[float]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
    ):
        raise ValueError(f"{name} must contain two values")
    answer = [float(item) for item in value]
    if any(not math.isfinite(item) or item < 0.0 for item in answer):
        raise ValueError(f"{name} must be finite and nonnegative")
    return answer


def _method_summary(
    value: Any,
    *,
    method: str,
    coordinate_time: float,
    direct_tail_ceiling: float,
    endpoint: bool,
) -> dict[str, Any]:
    record = _require_mapping(f"{method} common event", value)
    if (
        record.get("method") != method
        or tuple(record.get("point_counts", ())) != EXPECTED_POINT_COUNTS
        or float(record.get("coordinate_time", math.nan)) != coordinate_time
        or record.get("target_field") != EXPECTED_TARGET_FIELD
    ):
        raise ValueError(f"{method} common-event identity differs")

    ratios = _finite_pair(
        f"{method} target derivative-tail ratios",
        record.get("target_derivative_tail_ratios"),
    )
    pair_passed = [item < direct_tail_ceiling for item in ratios]
    if record.get("target_direct_pair_passed") != pair_passed:
        raise ValueError(f"{method} target pair classification differs")

    constraint = _require_mapping(
        f"{method} constraint admission", record.get("constraint_admission")
    )
    profiles = _require_mapping(
        f"{method} profile convergence", record.get("profile_convergence")
    )
    raw_spectrum = _require_mapping(
        f"{method} raw spectrum", record.get("raw_spatial_spectral_admission")
    )
    _require_bool(
        f"{method} constraint admission", constraint.get("admission_passed"), True
    )
    _require_bool(
        f"{method} complete profile contraction",
        profiles.get("every_field_contracted"),
        True,
    )
    _require_bool(
        f"{method} target profile contraction",
        record.get("target_complete_profile_contraction_passed"),
        True,
    )
    _require_bool(
        f"{method} target absolute budgets",
        record.get("target_medium_and_fine_absolute_budgets_passed"),
        True,
    )
    _require_bool(
        f"{method} medium/fine absolute budgets",
        record.get("all_medium_and_fine_absolute_budgets_passed"),
        True,
    )
    _require_bool(
        f"{method} generic raw all-field spectral admission",
        raw_spectrum.get("admission_passed"),
        False,
    )
    _require_bool(
        f"{method} generic nested-tail aggregate",
        raw_spectrum.get("every_nested_tail_ratio_passed"),
        False,
    )
    _require_bool(
        f"{method} targeted study admission",
        record.get("study_admission_passed"),
        endpoint,
    )

    minimum_order = float(
        constraint.get("minimum_finite_finest_pair_order", math.nan)
    )
    if not math.isfinite(minimum_order) or minimum_order <= 0.0:
        raise ValueError(f"{method} minimum finite constraint order differs")

    return {
        "method": method,
        "target_derivative_tail_ratios": ratios,
        "strict_ratio_margins": [direct_tail_ceiling - item for item in ratios],
        "target_direct_pair_passed": pair_passed,
        "targeted_study_admission_passed": bool(
            record["study_admission_passed"]
        ),
        "constraint_admission_passed": True,
        "minimum_finite_finest_pair_order": minimum_order,
        "complete_profile_contraction_passed": True,
        "target_absolute_budgets_passed": True,
        "generic_raw_all_field_spectral_admission_passed": False,
        "generic_every_nested_tail_ratio_passed": False,
        "accepted_stage_counts": list(constraint["accepted_stage_counts"]),
        "owned_raw_constraint_norms": list(
            constraint["owned_raw_global_norms"]
        ),
    }


def _assessment_summary(
    value: Any,
    *,
    coordinate_time: float,
    direct_tail_ceiling: float,
    endpoint: bool,
) -> dict[str, Any]:
    assessment = _require_mapping("RSP1 assessment", value)
    if float(assessment.get("coordinate_time", math.nan)) != coordinate_time:
        raise ValueError("RSP1 assessment time differs")
    _require_bool(
        "RSP1 endpoint question",
        assessment.get("endpoint_question_passed"),
        endpoint,
    )
    _require_bool(
        "RSP1 collapse/trapped classification",
        assessment.get("collapse_or_trapped_outcome_classified"),
        False,
    )
    primary = _method_summary(
        assessment.get("primary_common_event"),
        method="RK4",
        coordinate_time=coordinate_time,
        direct_tail_ceiling=direct_tail_ceiling,
        endpoint=endpoint,
    )
    comparator = _method_summary(
        assessment.get("comparator_common_event"),
        method="SSPRK3",
        coordinate_time=coordinate_time,
        direct_tail_ceiling=direct_tail_ceiling,
        endpoint=endpoint,
    )
    return {
        "coordinate_time": coordinate_time,
        "endpoint_question_passed": endpoint,
        "collapse_or_trapped_outcome_classified": False,
        "methods": {"RK4": primary, "SSPRK3": comparator},
    }


def diagnose_rsp1_campaign(
    events: Sequence[Mapping[str, Any]],
    result: Mapping[str, Any],
    *,
    study_id: str,
    direct_tail_ceiling: float = 0.25,
) -> dict[str, Any]:
    """Reduce one terminal RSP1 run without promoting a physical claim."""

    if (
        not math.isfinite(direct_tail_ceiling)
        or direct_tail_ceiling <= 0.0
    ):
        raise ValueError("direct-tail ceiling must be positive and finite")
    if len(events) != 2:
        raise ValueError("RSP1 event log must contain exactly two records")
    initial_event = _require_mapping("RSP1 initial event", events[0])
    endpoint_event = _require_mapping("RSP1 endpoint event", events[1])
    if (
        initial_event.get("event_type") != "initial_premise_event"
        or initial_event.get("event_index") != 0
        or initial_event.get("amplitude") != "3"
        or initial_event.get(
            "t0_target_ratio_is_not_the_evolved_endpoint_outcome"
        )
        is not True
        or endpoint_event.get("event_type") != "evolved_resolution_endpoint"
        or endpoint_event.get("event_index") != 1
        or endpoint_event.get("amplitude") != "3"
        or endpoint_event.get("classification") != EXPECTED_CLASSIFICATION
    ):
        raise ValueError("RSP1 event sequence differs")

    if (
        result.get("schema_version") != 1
        or result.get("runner_id") != EXPECTED_RUNNER_ID
        or result.get("study_id") != study_id
        or result.get("classification") != EXPECTED_CLASSIFICATION
        or result.get("amplitude") != "3"
        or tuple(result.get("point_counts", ())) != EXPECTED_POINT_COUNTS
        or float(result.get("target_coordinate_time", math.nan)) != 0.0625
        or result.get("stop") is not None
        or result.get("endpoint_assessment") != endpoint_event.get("assessment")
    ):
        raise ValueError("RSP1 terminal result differs")

    for key in (
        "fresh_GR0_calibration_completed",
        "SGBL_outcome_read",
        "FGCQR_outcome_read",
        "holdout_execution_authorized",
        "retained_EFT_evolution_authorized",
        "physical_transition_claim_authorized",
        "mechanism_question_answered",
    ):
        _require_bool(f"RSP1 terminal nonclaim {key}", result.get(key), False)
    _require_bool(
        "RSP1 amplitude-three veto clearance",
        result.get("amplitude_three_spectral_veto_cleared_for_successor_design"),
        True,
    )

    initial = _assessment_summary(
        initial_event.get("assessment"),
        coordinate_time=0.0,
        direct_tail_ceiling=direct_tail_ceiling,
        endpoint=False,
    )
    endpoint = _assessment_summary(
        endpoint_event.get("assessment"),
        coordinate_time=0.0625,
        direct_tail_ceiling=direct_tail_ceiling,
        endpoint=True,
    )
    for method in EXPECTED_METHODS:
        if any(initial["methods"][method]["target_direct_pair_passed"]):
            raise ValueError(f"{method} t0 target veto was not preserved")
        if not all(endpoint["methods"][method]["target_direct_pair_passed"]):
            raise ValueError(f"{method} evolved target pair did not pass")

    endpoint_assessment = _require_mapping(
        "RSP1 endpoint assessment", result["endpoint_assessment"]
    )
    members = _require_mapping(
        "RSP1 member diagnostics", endpoint_assessment.get("member_diagnostics")
    )
    if tuple(sorted(members)) != tuple(sorted(EXPECTED_MEMBER_KEYS)):
        raise ValueError("RSP1 terminal member set differs")
    member_summary: dict[str, Any] = {}
    total_source_retries = 0
    minimum_boundary_margin = math.inf
    for key in EXPECTED_MEMBER_KEYS:
        item = _require_mapping(f"RSP1 member {key}", members[key])
        step_index = int(item.get("step_index", -1))
        transaction_serial = int(item.get("transaction_serial", -1))
        accepted_stage_count = int(item.get("accepted_stage_count", -1))
        cfl_retries = int(item.get("CFL_retry_count", -1))
        source_retries = int(item.get("source_retry_count", -1))
        boundary_margin = float(item.get("remaining_boundary_margin", math.nan))
        characteristic_distance = float(
            item.get("accumulated_characteristic_distance", math.nan)
        )
        state_sha = item.get("state_sha256")
        if (
            step_index <= 0
            or transaction_serial <= 0
            or accepted_stage_count <= 0
            or cfl_retries < 0
            or source_retries != 0
            or not math.isfinite(boundary_margin)
            or boundary_margin <= 0.0
            or not math.isfinite(characteristic_distance)
            or characteristic_distance <= 0.0
            or not isinstance(state_sha, str)
            or len(state_sha) != 64
        ):
            raise ValueError(f"RSP1 member diagnostics differ: {key}")
        total_source_retries += source_retries
        minimum_boundary_margin = min(minimum_boundary_margin, boundary_margin)
        member_summary[key] = {
            "step_index": step_index,
            "transaction_serial": transaction_serial,
            "accepted_stage_count": accepted_stage_count,
            "CFL_retry_count": cfl_retries,
            "source_retry_count": source_retries,
            "accumulated_characteristic_distance": characteristic_distance,
            "remaining_boundary_margin": boundary_margin,
            "state_sha256": state_sha,
        }

    return {
        "event_log_line_count": len(events),
        "study_id": study_id,
        "classification": EXPECTED_CLASSIFICATION,
        "amplitude": "3",
        "point_counts": list(EXPECTED_POINT_COUNTS),
        "direct_tail_ceiling": direct_tail_ceiling,
        "initial_premise": initial,
        "evolved_endpoint": endpoint,
        "members": member_summary,
        "all_six_members_reached_endpoint": True,
        "total_source_retry_count": total_source_retries,
        "minimum_remaining_boundary_margin": minimum_boundary_margin,
        "amplitude_three_spectral_veto_cleared_for_successor_design": True,
        "generic_all_field_raw_spectral_admission_passed": False,
        "GR0_case_eligible": False,
        "collapse_or_trapped_outcome_classified": False,
        "candidate_or_mechanism_outcome_read": False,
    }


__all__ = [
    "EXPECTED_CLASSIFICATION",
    "EXPECTED_MEMBER_KEYS",
    "EXPECTED_METHODS",
    "EXPECTED_POINT_COUNTS",
    "EXPECTED_RUNNER_ID",
    "EXPECTED_TARGET_FIELD",
    "diagnose_rsp1_campaign",
]
