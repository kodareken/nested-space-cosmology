"""Outcome-neutral diagnosis of the completed PROTO13 GR-0 campaign.

The reducer consumes canonical event/result records only.  It establishes the
exact gate at which the frozen calibration stopped and counts the failed
temporal spectral predicates without changing a threshold, evolving a state,
or assigning the stop to physics.  Terminal checkpoint restoration and the
independent numerical recomputation are owned by the PREF15 reproducer.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from math import isfinite
from numbers import Real
from typing import Any


EXPECTED_METHODS = ("RK4", "SSPRK3")
EXPECTED_POINT_COUNTS = {
    "RK4": (2049, 4097, 8193),
    "SSPRK3": (4097, 8193, 16385),
}
EXPECTED_FIELDS = (
    "alpha_minus_1",
    "shift_over_r",
    "lambda_minus_1",
    "R_over_r_minus_1",
    "phi_over_Lambda",
    "chi_over_Lambda",
)
RESTART_EVENT_INDEX = 23
TERMINAL_EVENT_INDEX = 63
TEMPORAL_MINIMUM_SAMPLES = 64
TEMPORAL_MAXIMUM_NESTED_TAIL_RATIO = 1.0 / 4.0


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _finite(name: str, value: object, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if nonnegative and answer < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


def _field_from_key(key: str) -> str:
    prefix, separator, field = key.partition(":")
    _require(
        separator == ":"
        and prefix.startswith("tracer_")
        and prefix[7:].isdigit()
        and field in EXPECTED_FIELDS,
        "PROTO13 temporal spectral key differs",
    )
    return field


def summarize_temporal_method(
    method: str,
    admission: Mapping[str, Any],
    *,
    maximum_nested_tail_ratio: Real = TEMPORAL_MAXIMUM_NESTED_TAIL_RATIO,
) -> dict[str, Any]:
    """Validate and summarize one serialized terminal temporal admission."""

    _require(method in EXPECTED_METHODS, "PROTO13 temporal method differs")
    threshold = _finite(
        "maximum_nested_tail_ratio", maximum_nested_tail_ratio, nonnegative=True
    )
    _require(
        threshold == TEMPORAL_MAXIMUM_NESTED_TAIL_RATIO,
        "PROTO13 temporal nested-tail threshold differs",
    )
    nested = admission.get("nested_admission")
    _require(
        admission.get("point_counts") == list(EXPECTED_POINT_COUNTS[method])
        and admission.get("sample_count") == TEMPORAL_MINIMUM_SAMPLES
        and admission.get("tracer_count") == 48
        and admission.get("field_names") == list(EXPECTED_FIELDS)
        and admission.get("admission_passed") is False
        and isinstance(admission.get("maximum_round_trip_interpolation_infinity"), Real)
        and isinstance(nested, Mapping),
        f"PROTO13 {method} temporal admission boundary differs",
    )
    interpolation = _finite(
        "maximum temporal round-trip interpolation",
        admission["maximum_round_trip_interpolation_infinity"],
        nonnegative=True,
    )
    _require(
        nested.get("grid_point_counts") == list(EXPECTED_POINT_COUNTS[method])
        and nested.get("every_individual_budget_passed") is False
        and nested.get("every_nested_tail_ratio_passed") is False
        and nested.get("admission_passed") is False,
        f"PROTO13 {method} nested temporal boundary differs",
    )

    expected_keys = {
        f"tracer_{tracer:04d}:{field}"
        for tracer in range(48)
        for field in EXPECTED_FIELDS
    }
    ratio_summaries: dict[str, Any] = {}
    for ratio_name in (
        "field_power_tail_ratios",
        "derivative_power_tail_ratios",
    ):
        ratios = nested.get(ratio_name)
        _require(
            isinstance(ratios, Mapping) and set(ratios) == expected_keys,
            f"PROTO13 {method} {ratio_name} key set differs",
        )
        failed_by_pair = [0, 0]
        failed_by_field: Counter[str] = Counter()
        none_count = 0
        largest: tuple[float, str, int] | None = None
        total = 0
        for key in sorted(ratios):
            values = ratios[key]
            _require(
                isinstance(values, list) and len(values) == 2,
                f"PROTO13 {method} {ratio_name} pair count differs",
            )
            field = _field_from_key(key)
            for pair_index, value in enumerate(values):
                total += 1
                if value is None:
                    none_count += 1
                    failed_by_pair[pair_index] += 1
                    failed_by_field[field] += 1
                    continue
                observed = _finite(
                    f"PROTO13 {method} {ratio_name}", value, nonnegative=True
                )
                if largest is None or observed > largest[0]:
                    largest = (observed, key, pair_index)
                if observed >= threshold:
                    failed_by_pair[pair_index] += 1
                    failed_by_field[field] += 1
        failed = sum(failed_by_pair)
        _require(
            total == 576 and failed > 0 and largest is not None,
            f"PROTO13 {method} {ratio_name} failure set differs",
        )
        ratio_summaries[ratio_name] = {
            "comparison_count": total,
            "failed_comparison_count": failed,
            "none_ratio_count": none_count,
            "failed_by_adjacent_pair": failed_by_pair,
            "failed_by_field": {
                field: failed_by_field[field] for field in EXPECTED_FIELDS
            },
            "largest_finite_ratio": largest[0],
            "largest_finite_ratio_key": largest[1],
            "largest_finite_ratio_pair_index": largest[2],
        }

    return {
        "method": method,
        "point_counts": list(EXPECTED_POINT_COUNTS[method]),
        "sample_count": TEMPORAL_MINIMUM_SAMPLES,
        "tracer_count": 48,
        "field_count": len(EXPECTED_FIELDS),
        "maximum_round_trip_interpolation_infinity": interpolation,
        "individual_budget_admission_passed": False,
        "nested_tail_ratio_admission_passed": False,
        "admission_passed": False,
        "ratio_summaries": ratio_summaries,
    }


def diagnose_proto13_campaign(
    events: Sequence[Mapping[str, Any]],
    campaign_result: Mapping[str, Any],
    *,
    maximum_nested_tail_ratio: Real = TEMPORAL_MAXIMUM_NESTED_TAIL_RATIO,
) -> dict[str, Any]:
    """Return the bounded CAL10/PREF15 terminal campaign diagnosis."""

    threshold = _finite(
        "maximum_nested_tail_ratio", maximum_nested_tail_ratio, nonnegative=True
    )
    _require(
        threshold == TEMPORAL_MAXIMUM_NESTED_TAIL_RATIO,
        "PROTO13 temporal threshold changed after outcome",
    )
    records = tuple(events)
    amplitude_records = campaign_result.get("amplitude_records")
    _require(
        campaign_result.get("runner_id") == "FGC-1-CAL10-RUN1-RUNNER"
        and campaign_result.get("classification")
        == "calibration_failed_no_eligible_GR0_case"
        and campaign_result.get("selected_amplitude") is None
        and campaign_result.get("holdout_execution_authorized") is False
        and campaign_result.get("mechanism_question_answered") is False
        and campaign_result.get("retained_EFT_evolution_authorized") is False
        and campaign_result.get("SGBL_outcome_read") is False
        and campaign_result.get("FGCQR_outcome_read") is False
        and campaign_result.get("method_owned_ladders_enabled") is True
        and campaign_result.get("common_physical_node_comparison_enabled") is True
        and isinstance(amplitude_records, list)
        and len(amplitude_records) == 1,
        "PROTO13 campaign result boundary differs",
    )
    amplitude = amplitude_records[0]
    expected_stop = {
        "classification": "common_event_constraint_or_spectral_stop",
        "reason": "causal_past_temporal_spectral_admission",
        "event_index": TERMINAL_EVENT_INDEX,
        "coordinate_time": TERMINAL_EVENT_INDEX / 16.0,
    }
    source_counts = amplitude.get("member_source_retry_counts")
    cfl_counts = amplitude.get("member_CFL_retry_counts")
    state_hashes = amplitude.get("member_final_state_hashes")
    expected_members = {
        f"{method}-{count}"
        for method, counts in EXPECTED_POINT_COUNTS.items()
        for count in counts
    }
    _require(
        amplitude.get("amplitude") == "3"
        and amplitude.get("eligible") is False
        and amplitude.get("last_completed_common_event_index")
        == TERMINAL_EVENT_INDEX
        and amplitude.get("last_completed_coordinate_time")
        == TERMINAL_EVENT_INDEX / 16.0
        and amplitude.get("consecutive_qualified_trapped_common_events") == 0
        and amplitude.get("stop") == expected_stop
        and isinstance(source_counts, Mapping)
        and set(source_counts) == expected_members
        and all(value == 0 for value in source_counts.values())
        and isinstance(cfl_counts, Mapping)
        and set(cfl_counts) == expected_members
        and all(isinstance(value, int) and not isinstance(value, bool) and value >= 0 for value in cfl_counts.values())
        and any(value > 0 for value in cfl_counts.values())
        and isinstance(state_hashes, Mapping)
        and set(state_hashes) == expected_members,
        "PROTO13 terminal amplitude, retry, or state-hash boundary differs",
    )

    _require(
        len(records) == TERMINAL_EVENT_INDEX - RESTART_EVENT_INDEX + 1
        and [record.get("event_index") for record in records]
        == list(range(RESTART_EVENT_INDEX, TERMINAL_EVENT_INDEX + 1)),
        "PROTO13 event sequence differs",
    )
    for position, event in enumerate(records):
        event_index = RESTART_EVENT_INDEX + position
        expected_type = "restart_common_event" if position == 0 else "common_event"
        assessment = event.get("assessment")
        _require(
            event.get("event_type") == expected_type
            and isinstance(assessment, Mapping)
            and assessment.get("coordinate_time") == event_index / 16.0
            and assessment.get("method_owned_ladders_enforced") is True
            and assessment.get("common_physical_node_comparison_enforced") is True
            and assessment.get("qualified_trapped_common_event") is False,
            f"PROTO13 event {event_index} identity or premise differs",
        )
        if position == 0:
            _require(
                event.get("trajectory_advanced_by_PROTO13") is False
                and event.get("counts_toward_new_consecutive_trapped_events") is False,
                "PROTO13 restart event semantics differ",
            )
        else:
            _require(
                event.get("amplitude") == "3"
                and event.get("consecutive_qualified_trapped_common_events") == 0,
                f"PROTO13 event {event_index} qualification history differs",
            )
        primary = assessment.get("primary_common_event")
        comparator = assessment.get("comparator_common_event")
        trapped = assessment.get("trapped_assessment")
        temporal = assessment.get("temporal_spectral_admission")
        _require(
            isinstance(primary, Mapping)
            and primary.get("method") == "RK4"
            and primary.get("point_counts") == list(EXPECTED_POINT_COUNTS["RK4"])
            and primary.get("admission_passed") is True
            and isinstance(comparator, Mapping)
            and comparator.get("method") == "SSPRK3"
            and comparator.get("point_counts")
            == list(EXPECTED_POINT_COUNTS["SSPRK3"])
            and comparator.get("admission_passed") is True
            and isinstance(trapped, Mapping)
            and trapped.get("both_method_owned_common_event_admissions_passed")
            is True
            and trapped.get("trapped_sign_passed") is False
            and isinstance(temporal, Mapping),
            f"PROTO13 event {event_index} common or trapped admission differs",
        )
        expected_samples = event_index + 1
        _require(
            temporal.get("sample_count") == expected_samples,
            f"PROTO13 event {event_index} temporal sample count differs",
        )
        if event_index < TERMINAL_EVENT_INDEX:
            _require(
                temporal
                == {
                    "available": False,
                    "required_minimum_samples": TEMPORAL_MINIMUM_SAMPLES,
                    "sample_count": expected_samples,
                    "admission_passed": False,
                },
                f"PROTO13 event {event_index} pre-temporal boundary differs",
            )

    terminal = records[-1]["assessment"]
    temporal = terminal["temporal_spectral_admission"]
    _require(
        temporal.get("available") is True
        and temporal.get("sample_count") == TEMPORAL_MINIMUM_SAMPLES
        and temporal.get("admission_passed") is False
        and isinstance(temporal.get("methods"), Mapping)
        and set(temporal["methods"]) == set(EXPECTED_METHODS),
        "PROTO13 terminal temporal composition differs",
    )
    methods = {
        method: summarize_temporal_method(
            method,
            temporal["methods"][method],
            maximum_nested_tail_ratio=threshold,
        )
        for method in EXPECTED_METHODS
    }
    return {
        "event_log_line_count": len(records),
        "restart_event_index": RESTART_EVENT_INDEX,
        "new_common_event_count": TERMINAL_EVENT_INDEX - RESTART_EVENT_INDEX,
        "terminal_event_index": TERMINAL_EVENT_INDEX,
        "terminal_coordinate_time": TERMINAL_EVENT_INDEX / 16.0,
        "terminal_classification": campaign_result["classification"],
        "terminal_stop": expected_stop,
        "serialized_source_retry_count": 0,
        "terminal_CFL_retry_count_by_member": dict(cfl_counts),
        "some_adaptive_CFL_step_reductions_occurred": True,
        "campaign_stop_was_not_CFL_retry_exhaustion": True,
        "terminal_common_spatial_and_constraint_admissions_passed": True,
        "terminal_trapped_sign_passed": False,
        "temporal_minimum_samples_reached": True,
        "temporal_maximum_nested_tail_ratio": threshold,
        "temporal_methods": methods,
        "both_temporal_method_admissions_failed": True,
        "both_temporal_individual_budget_layers_failed": True,
        "both_temporal_nested_ratio_layers_failed": True,
        "temporal_failure_cause_derived": False,
        "GR0_case_eligible": False,
        "candidate_or_mechanism_outcome_read": False,
    }


__all__ = [
    "EXPECTED_FIELDS",
    "EXPECTED_METHODS",
    "EXPECTED_POINT_COUNTS",
    "RESTART_EVENT_INDEX",
    "TEMPORAL_MAXIMUM_NESTED_TAIL_RATIO",
    "TEMPORAL_MINIMUM_SAMPLES",
    "TERMINAL_EVENT_INDEX",
    "diagnose_proto13_campaign",
    "summarize_temporal_method",
]
