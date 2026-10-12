"""Outcome-neutral diagnosis of the completed PROTO11 GR-0 campaign.

The module consumes only canonical campaign records.  It does not evolve a
state, alter a threshold, or authorize a successor.  Its two jobs are to
separate the amplitude-``5/2`` affine-source retry exhaustion from every
non-source monitor and to identify the exact direct spectral veto on the
amplitude-``3`` evolved common event.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import isfinite
from numbers import Real
from typing import Any

import numpy as np


EXPECTED_AMPLITUDES = ("5/2", "3")
EXPECTED_POINT_COUNTS = (2049, 4097, 8193)
EXPECTED_METHODS = ("RK4", "SSPRK3")
SOURCE_EVENT = "rejected_unaccepted_source_only_proposal"
COMMON_EVENT = "common_event"


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _common_event(
    events: Sequence[Mapping[str, Any]],
    *,
    amplitude: str,
    event_index: int,
    coordinate_time: float,
) -> Mapping[str, Any]:
    selected = [
        item
        for item in events
        if item.get("event_type") == COMMON_EVENT
        and item.get("amplitude") == amplitude
        and item.get("event_index") == event_index
    ]
    _require(len(selected) == 1, "PROTO11 common-event identity differs")
    event = selected[0]
    assessment = event.get("assessment", {})
    _require(
        assessment.get("coordinate_time") == coordinate_time,
        "PROTO11 common-event time differs",
    )
    return event


def _source_groups(
    records: Sequence[Mapping[str, Any]],
) -> list[list[Mapping[str, Any]]]:
    groups: list[list[Mapping[str, Any]]] = []
    for item in records:
        if not groups or groups[-1][0].get("initial_time") != item.get(
            "initial_time"
        ):
            groups.append([])
        groups[-1].append(item)
    return groups


def diagnose_amplitude_five_halves_source_obstruction(
    events: Sequence[Mapping[str, Any]],
    amplitude_record: Mapping[str, Any],
    *,
    raw_source_limit: Real,
    minimum_step_size: Real,
    accepted_increment_time_ulp_budget: int = 2,
) -> dict[str, Any]:
    """Prove the typed evolved-source obstruction for amplitude ``5/2``."""

    raw_limit = _finite("raw_source_limit", raw_source_limit, positive=True)
    minimum_step = _finite(
        "minimum_step_size", minimum_step_size, positive=True
    )
    if (
        isinstance(accepted_increment_time_ulp_budget, bool)
        or not isinstance(accepted_increment_time_ulp_budget, int)
        or accepted_increment_time_ulp_budget < 1
    ):
        raise ValueError("accepted increment ULP budget must be positive")

    rejected = [
        item for item in events if item.get("event_type") == SOURCE_EVENT
    ]
    _require(len(rejected) > 0, "PROTO11 source-retry trace is absent")
    _require(
        all(item.get("amplitude") == "5/2" for item in rejected),
        "PROTO11 source retries are not confined to amplitude 5/2",
    )
    groups = _source_groups(rejected)
    serialized_groups: list[dict[str, Any]] = []
    all_residuals: list[float] = []
    all_stage_times: list[float] = []

    for group_index, group in enumerate(groups):
        first = group[0]
        initial_time = _finite("initial_time", first.get("initial_time"))
        initial_step = _finite(
            "initial_step_size_for_accepted_step",
            first.get("initial_step_size_for_accepted_step"),
            positive=True,
        )
        accepted_hash = first.get("accepted_state_sha256")
        accepted_step_index = first.get("accepted_step_index")
        accepted_serial = first.get("accepted_transaction_serial")
        group_residuals: list[float] = []
        group_stage_times: list[float] = []

        for retry_index, item in enumerate(group, start=1):
            failed = item.get("failed_evaluations", [])
            expected_step = initial_step * (0.5 ** retry_index)
            _require(
                item.get("member") == "RK4-8193"
                and item.get("method") == "RK4"
                and item.get("point_count") == 8193
                and item.get("retry_count_for_accepted_step") == retry_index
                and item.get("attempted_step_size") == expected_step
                and item.get("initial_time") == initial_time
                and item.get("accepted_state_sha256") == accepted_hash
                and item.get("accepted_step_index") == accepted_step_index
                and item.get("accepted_transaction_serial") == accepted_serial
                and item.get("complete_failure_set")
                == ["newton_residual_limit"]
                and item.get("retryable_source_only") is True
                and item.get("non_source_failure_vetoed_retry") is False
                and item.get("fields_bitwise_preserved") is True
                and item.get("time_advanced") is False
                and item.get("accepted_stage_count_advanced") is False
                and item.get("causal_debit_advanced") is False
                and item.get("external_transaction_state_preserved") is True
                and item.get("preaccept_called") is False
                and item.get("complete_failure_evidence_serialized_before_retry")
                is True
                and isinstance(failed, list)
                and len(failed) > 0,
                "PROTO11 source-only rollback contract differs",
            )
            for evaluation in failed:
                residual = _finite(
                    "source_residual_infinity",
                    evaluation.get("source_residual_infinity"),
                    positive=True,
                )
                stage_time = _finite("source stage time", evaluation.get("stage_time"))
                _require(
                    evaluation.get("failures") == ["newton_residual_limit"]
                    and evaluation.get("source_residual_maximum") == raw_limit
                    and residual > raw_limit,
                    "PROTO11 failed source evaluation differs",
                )
                group_residuals.append(residual)
                group_stage_times.append(stage_time)
        all_residuals.extend(group_residuals)
        all_stage_times.extend(group_stage_times)

        accepted_advance: float | None = None
        accepted_matches_half_step: bool | None = None
        if group_index + 1 < len(groups):
            next_group = groups[group_index + 1]
            next_time = _finite("next accepted time", next_group[0]["initial_time"])
            accepted_advance = next_time - initial_time
            expected_advance = group[-1]["attempted_step_size"] * 0.5
            time_ulp = abs(float(np.spacing(np.float64(initial_time))))
            accepted_matches_half_step = (
                abs(accepted_advance - expected_advance)
                <= accepted_increment_time_ulp_budget * time_ulp
            )
            _require(
                accepted_advance > 0.0
                and accepted_matches_half_step
                and next_group[0].get("accepted_step_index")
                == accepted_step_index + 1
                and next_group[0].get("accepted_transaction_serial")
                == accepted_serial + 5,
                "PROTO11 accepted half-step sequence differs",
            )

        serialized_groups.append(
            {
                "group_index": group_index,
                "accepted_time": initial_time,
                "accepted_step_index": accepted_step_index,
                "accepted_transaction_serial": accepted_serial,
                "accepted_state_sha256": accepted_hash,
                "failed_retry_count": len(group),
                "first_failed_step_size": group[0]["attempted_step_size"],
                "last_failed_step_size": group[-1]["attempted_step_size"],
                "inferred_next_accepted_advance": accepted_advance,
                "next_accepted_advance_matches_exact_half_step_within_time_ULP_budget": (
                    accepted_matches_half_step
                ),
                "minimum_failed_residual": min(group_residuals),
                "maximum_failed_residual": max(group_residuals),
                "minimum_failed_stage_time": min(group_stage_times),
                "maximum_failed_stage_time": max(group_stage_times),
            }
        )

    stop = amplitude_record.get("stop", {})
    evidence = stop.get("complete_failure_evidence", {})
    last = rejected[-1]
    next_forbidden = _finite(
        "terminal attempted step", evidence.get("attempted_step_size"), positive=True
    )
    _require(
        amplitude_record.get("amplitude") == "5/2"
        and amplitude_record.get("eligible") is False
        and amplitude_record.get("last_completed_common_event_index") == 1
        and amplitude_record.get("last_completed_coordinate_time") == 0.0625
        and stop.get("classification") == "scientific_source_retry_exhausted"
        and stop.get("reason") == "newton_residual_limit"
        and stop.get("member") == "RK4-8193"
        and stop.get("method") == "RK4"
        and stop.get("point_count") == 8193
        and stop.get("target_common_event_index") == 2
        and stop.get("last_fully_completed_common_event_index") == 1
        and stop.get("cross_member_state_rolled_back") is True
        and stop.get("cross_member_rollback_verified") is True
        and evidence.get("exhaustion_kind") == "minimum_step_size"
        and evidence.get("minimum_step_size") == minimum_step
        and evidence.get("last_rejected_proposal") == last
        and evidence.get("source_retry_count_for_accepted_step")
        == len(groups[-1])
        and last.get("attempted_step_size") >= minimum_step
        and next_forbidden < minimum_step
        and next_forbidden == last.get("attempted_step_size") * 0.5,
        "PROTO11 terminal source-exhaustion record differs",
    )

    first_failed_step = rejected[0]["attempted_step_size"]
    last_failed_step = rejected[-1]["attempted_step_size"]
    return {
        "amplitude": "5/2",
        "rejected_source_only_proposal_count": len(rejected),
        "accepted_time_group_count": len(groups),
        "all_rejections_member": "RK4-8193",
        "all_rejections_are_raw_newton_residual_limit_only": True,
        "all_rejections_preserve_the_complete_accepted_boundary": True,
        "all_retry_sequences_are_exact_binary_halvings": True,
        "accepted_increment_time_ulp_budget": accepted_increment_time_ulp_budget,
        "all_inferred_accepted_advances_match_the_next_half_step": True,
        "first_failed_step_size": first_failed_step,
        "last_legal_failed_step_size": last_failed_step,
        "first_to_last_failed_step_ratio": first_failed_step / last_failed_step,
        "next_forbidden_step_size": next_forbidden,
        "frozen_minimum_step_size": minimum_step,
        "raw_source_residual_limit": raw_limit,
        "minimum_failed_source_residual": min(all_residuals),
        "maximum_failed_source_residual": max(all_residuals),
        "minimum_residual_over_raw_limit": min(all_residuals) / raw_limit,
        "first_failed_stage_time": min(
            evaluation["stage_time"]
            for evaluation in rejected[0]["failed_evaluations"]
        ),
        "last_failed_stage_time": max(
            evaluation["stage_time"]
            for evaluation in rejected[-1]["failed_evaluations"]
        ),
        "last_accepted_time": rejected[-1]["initial_time"],
        "last_accepted_step_index": rejected[-1]["accepted_step_index"],
        "last_accepted_transaction_serial": rejected[-1][
            "accepted_transaction_serial"
        ],
        "terminal_classification": stop["classification"],
        "terminal_exhaustion_kind": evidence["exhaustion_kind"],
        "cross_member_rollback_verified": True,
        "groups": serialized_groups,
    }


def diagnose_amplitude_three_direct_spectral_veto(
    event: Mapping[str, Any],
    amplitude_record: Mapping[str, Any],
    *,
    direct_tail_ceiling: Real,
) -> dict[str, Any]:
    """Identify the sole direct coarse-to-medium veto for amplitude ``3``."""

    ceiling = _finite("direct_tail_ceiling", direct_tail_ceiling, positive=True)
    assessment = event.get("assessment", {})
    primary = assessment.get("primary_common_event", {})
    comparator = assessment.get("comparator_common_event", {})
    trapped = assessment.get("trapped_assessment", {})
    temporal = assessment.get("temporal_spectral_admission", {})
    methods: dict[str, Any] = {}

    for method, common in zip(EXPECTED_METHODS, (primary, comparator), strict=True):
        constraint = common.get("constraint_admission", {})
        spatial = common.get("spatial_spectral_admission", {})
        direct_fields = spatial.get(
            "coarse_to_medium_direct_passed_by_field", {}
        )
        direct_derivatives = spatial.get(
            "coarse_to_medium_direct_passed_by_derivative", {}
        )
        field_ratios = spatial.get("field_power_tail_ratios", {})
        derivative_ratios = spatial.get("derivative_power_tail_ratios", {})
        failed_derivatives = sorted(
            key for key, passed in direct_derivatives.items() if passed is False
        )
        phi_ratios = derivative_ratios.get("phi_over_Lambda")
        _require(
            common.get("method") == method
            and common.get("point_counts") == list(EXPECTED_POINT_COUNTS)
            and common.get("accepted_stage_counts")
            and constraint.get("admission_passed") is True
            and constraint.get("minimum_finite_finest_pair_order") >= 1.5
            and common.get("admission_passed") is False
            and common.get("raw_PROTO9_admission_passed") is False
            and common.get("legacy_PROTO10_admission_passed") is False
            and common.get("reference_balanced_constraint_admission_passed")
            is True
            and common.get("reference_state_map_applied") is True
            and common.get("interior_q_reprojected") is False
            and spatial.get("admission_passed") is False
            and spatial.get("finest_pair_individual_budgets_passed") is True
            and spatial.get("every_profile_contracted") is True
            and all(value is True for value in direct_fields.values())
            and failed_derivatives == ["phi_over_Lambda"]
            and isinstance(phi_ratios, list)
            and len(phi_ratios) == 2
            and phi_ratios[0] >= ceiling
            and phi_ratios[1] < ceiling,
            "PROTO11 amplitude-three direct spectral veto differs",
        )
        methods[method] = {
            "point_counts": list(EXPECTED_POINT_COUNTS),
            "constraint_admission_passed": True,
            "minimum_finite_constraint_order": constraint[
                "minimum_finite_finest_pair_order"
            ],
            "all_coarse_to_medium_field_tails_passed_directly": True,
            "failed_coarse_to_medium_derivative_fields": failed_derivatives,
            "phi_derivative_tail_ratios": phi_ratios,
            "coarse_to_medium_phi_derivative_tail_passed": False,
            "medium_to_fine_phi_derivative_tail_passed": True,
            "direct_tail_ceiling": ceiling,
            "finest_pair_individual_budgets_passed": True,
            "every_whole_profile_contracted": True,
            "field_power_tail_ratios": field_ratios,
            "derivative_power_tail_ratios": derivative_ratios,
            "maximum_spatial_round_trip_interpolation_infinity": common[
                "maximum_spatial_round_trip_interpolation_infinity"
            ],
        }

    primary_scores = trapped.get("primary_trapped_scores", [])
    comparator_scores = trapped.get("comparator_trapped_scores", [])
    _require(
        amplitude_record.get("amplitude") == "3"
        and amplitude_record.get("eligible") is False
        and amplitude_record.get("last_completed_common_event_index") == 1
        and amplitude_record.get("last_completed_coordinate_time") == 0.0625
        and amplitude_record.get("stop")
        == {
            "classification": "common_event_constraint_or_spectral_stop",
            "reason": "common_event_constraint_or_spatial_spectral_admission",
            "event_index": 1,
            "coordinate_time": 0.0625,
        }
        and assessment.get("PROTO11_common_event_constraints_use_reference_map")
        is True
        and assessment.get("PROTO11_interior_q_reprojected") is False
        and assessment.get("PROTO11_complete_PROTO10_spectral_evidence_retained")
        is True
        and assessment.get("qualified_trapped_common_event") is False
        and trapped.get("both_common_event_admissions_passed") is False
        and trapped.get("trapped_sign_passed") is False
        and len(primary_scores) == 3
        and len(comparator_scores) == 3
        and all(value < 0.0 for value in (*primary_scores, *comparator_scores))
        and temporal
        == {
            "available": False,
            "admission_passed": False,
            "required_minimum_samples": 64,
        },
        "PROTO11 amplitude-three event classification differs",
    )
    return {
        "amplitude": "3",
        "coordinate_time": 0.0625,
        "event_index": 1,
        "source_only_rejection_count": 0,
        "both_constraint_admissions_passed": True,
        "both_spatial_admissions_passed": False,
        "sole_direct_coarse_to_medium_veto": "phi_over_Lambda derivative tail",
        "direct_tail_ceiling": ceiling,
        "methods": methods,
        "primary_trapped_scores": primary_scores,
        "comparator_trapped_scores": comparator_scores,
        "trapped_sign_passed": False,
        "temporal_spectrum_available": False,
        "interior_q_reprojected": False,
        "terminal_classification": amplitude_record["stop"]["classification"],
        "terminal_reason": amplitude_record["stop"]["reason"],
    }


def diagnose_proto11_campaign(
    events: Sequence[Mapping[str, Any]],
    campaign_result: Mapping[str, Any],
    *,
    raw_source_limit: Real = 1.0e-12,
    minimum_step_size: Real = 2.0 ** -30,
    direct_tail_ceiling: Real = 0.25,
    accepted_increment_time_ulp_budget: int = 2,
) -> dict[str, Any]:
    """Return the complete outcome-neutral CAL8 diagnosis."""

    records = tuple(events)
    amplitudes = campaign_result.get("amplitude_records", [])
    _require(
        campaign_result.get("runner_id") == "FGC-1-CAL8-RUN1-RUNNER"
        and campaign_result.get("classification")
        == "calibration_failed_no_eligible_GR0_case"
        and campaign_result.get("selected_amplitude") is None
        and campaign_result.get("holdout_execution_authorized") is False
        and campaign_result.get("SGBL_outcome_read") is False
        and campaign_result.get("FGCQR_outcome_read") is False
        and campaign_result.get("mechanism_question_answered") is False
        and campaign_result.get("retained_EFT_evolution_authorized") is False
        and campaign_result.get("PROTO11_reference_state_map_enabled") is True
        and campaign_result.get("PROTO11_interior_q_reprojection_enabled")
        is False
        and campaign_result.get(
            "PROTO11_complete_PROTO10_spectral_contract_retained"
        )
        is True
        and isinstance(amplitudes, list)
        and [item.get("amplitude") for item in amplitudes]
        == list(EXPECTED_AMPLITUDES),
        "PROTO11 campaign result boundary differs",
    )

    common = [item for item in records if item.get("event_type") == COMMON_EVENT]
    rejected = [item for item in records if item.get("event_type") == SOURCE_EVENT]
    _require(
        len(records) == len(common) + len(rejected)
        and [
            (item.get("amplitude"), item.get("event_index")) for item in common
        ]
        == [("5/2", 0), ("5/2", 1), ("3", 0), ("3", 1)],
        "PROTO11 event partition or common-event order differs",
    )
    five_t1 = _common_event(
        records, amplitude="5/2", event_index=1, coordinate_time=0.0625
    )
    three_t1 = _common_event(
        records, amplitude="3", event_index=1, coordinate_time=0.0625
    )
    for event in (five_t1, three_t1):
        assessment = event["assessment"]
        _require(
            assessment.get("trapped_assessment", {}).get("trapped_sign_passed")
            is False
            and assessment.get("temporal_spectral_admission", {}).get(
                "available"
            )
            is False,
            "PROTO11 first evolved event physical boundary differs",
        )
    five_assessment = five_t1["assessment"]
    _require(
        five_assessment.get("primary_common_event", {}).get("admission_passed")
        is True
        and five_assessment.get("comparator_common_event", {}).get(
            "admission_passed"
        )
        is True,
        "PROTO11 amplitude 5/2 first evolved event did not pass",
    )

    source = diagnose_amplitude_five_halves_source_obstruction(
        records,
        amplitudes[0],
        raw_source_limit=raw_source_limit,
        minimum_step_size=minimum_step_size,
        accepted_increment_time_ulp_budget=accepted_increment_time_ulp_budget,
    )
    spectral = diagnose_amplitude_three_direct_spectral_veto(
        three_t1,
        amplitudes[1],
        direct_tail_ceiling=direct_tail_ceiling,
    )
    return {
        "event_log_line_count": len(records),
        "common_event_count": len(common),
        "source_only_rejection_count": len(rejected),
        "both_amplitudes_reached_one_evolved_common_event": True,
        "amplitude_five_halves_first_evolved_event_fully_admitted": True,
        "amplitude_five_halves_source_obstruction": source,
        "amplitude_three_direct_spectral_obstruction": spectral,
        "no_GR0_amplitude_eligible": True,
        "candidate_or_mechanism_outcome_read": False,
    }


__all__ = [
    "COMMON_EVENT",
    "EXPECTED_AMPLITUDES",
    "EXPECTED_METHODS",
    "EXPECTED_POINT_COUNTS",
    "SOURCE_EVENT",
    "diagnose_amplitude_five_halves_source_obstruction",
    "diagnose_amplitude_three_direct_spectral_veto",
    "diagnose_proto11_campaign",
]
