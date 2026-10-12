"""Outcome-neutral diagnosis of the completed PROTO12 GR-0 campaign.

The reducer consumes only canonical campaign records.  It does not evolve a
state, change a threshold, infer a continuum cause from two amplitudes, or open
the SGB-L/FGC-QR holdout.  It proves where the frozen calibration contract
stopped and preserves the distinction between source retries and ordinary CFL
step reductions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import isfinite, log
from numbers import Real
from typing import Any


EXPECTED_AMPLITUDES = ("5/2", "3")
EXPECTED_METHODS = ("RK4", "SSPRK3")
EXPECTED_POINT_COUNTS = (2049, 4097, 8193)
COMMON_EVENT = "common_event"
CONSTRAINT_COMPONENTS = (
    "Hamiltonian",
    "radial_momentum",
    "gauge_t",
    "gauge_r",
    "reduction_alpha",
    "reduction_shift",
    "reduction_lambda",
    "reduction_R",
    "reduction_phi",
    "reduction_chi",
)
FINITE_CONSTRAINT_COMPONENTS = (
    "Hamiltonian",
    "radial_momentum",
    "gauge_t",
    "gauge_r",
)
TERMINAL_EVENT = {
    "5/2": (24, 1.5),
    "3": (23, 1.4375),
}


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


def _method_common_event(
    assessment: Mapping[str, Any], method: str
) -> Mapping[str, Any]:
    key = "primary_common_event" if method == "RK4" else "comparator_common_event"
    value = assessment.get(key)
    _require(isinstance(value, Mapping), f"PROTO12 {method} common event is absent")
    return value


def _recompute_constraint_orders(
    constraint: Mapping[str, Any], *, minimum_order: float
) -> dict[str, Any]:
    """Recompute the frozen order decision including roundoff enclosures."""

    counts = constraint.get("point_counts")
    raw = constraint.get("owned_domain_component_infinity")
    enclosures = constraint.get("accumulated_roundoff_enclosures")
    observed_orders = constraint.get("component_finest_pair_orders")
    observed_status = constraint.get("component_status")
    _require(
        counts == list(EXPECTED_POINT_COUNTS)
        and isinstance(raw, Mapping)
        and isinstance(enclosures, list)
        and len(enclosures) == 3
        and isinstance(observed_orders, Mapping)
        and isinstance(observed_status, Mapping)
        and set(raw) == set(CONSTRAINT_COMPONENTS)
        and set(observed_orders) == set(CONSTRAINT_COMPONENTS)
        and set(observed_status) == set(CONSTRAINT_COMPONENTS),
        "PROTO12 constraint-order inputs differ",
    )
    intervals = tuple(count - 1 for count in EXPECTED_POINT_COUNTS)
    refinement = intervals[-1] / intervals[-2]
    _require(refinement == 2.0, "PROTO12 refinement ratio differs")
    enclosure_values = tuple(
        _finite("accumulated roundoff enclosure", value, positive=True)
        for value in enclosures
    )

    effective_by_component: dict[str, list[float]] = {}
    recomputed_orders: dict[str, float | None] = {}
    recomputed_status: dict[str, str] = {}
    finite_orders: list[float] = []
    below_threshold: list[str] = []
    monotone = True
    for component in CONSTRAINT_COMPONENTS:
        values = raw[component]
        _require(
            isinstance(values, list) and len(values) == 3,
            f"PROTO12 {component} constraint norms differ",
        )
        raw_values = tuple(
            _finite(f"{component} constraint norm", value) for value in values
        )
        effective = tuple(
            max(value - enclosure, 0.0)
            for value, enclosure in zip(raw_values, enclosure_values, strict=True)
        )
        effective_by_component[component] = list(effective)
        if any(right > left for left, right in zip(effective[:-1], effective[1:], strict=True)):
            monotone = False
        medium, fine = effective[-2:]
        if medium == 0.0 and fine == 0.0:
            order = None
            status = "accumulated_roundoff_enclosed_zero_pair"
        elif medium > 0.0 and fine == 0.0:
            order = None
            status = "resolved_to_accumulated_roundoff_enclosure"
        elif medium == 0.0 and fine > 0.0:
            order = None
            status = "nonzero_reappeared_after_accumulated_roundoff_enclosure"
        else:
            order = log(medium / fine) / log(refinement)
            status = "finite_order"
            finite_orders.append(order)
            if order < minimum_order:
                below_threshold.append(component)
        _require(
            observed_orders[component] == order
            and observed_status[component] == status,
            f"PROTO12 {component} serialized constraint order differs",
        )
        recomputed_orders[component] = order
        recomputed_status[component] = status

    minimum_finite = min(finite_orders) if finite_orders else None
    _require(
        constraint.get("minimum_finite_finest_pair_order") == minimum_finite
        and constraint.get("monotone_refinement_passed") is monotone
        and constraint.get("finest_pair_order_passed")
        is (not below_threshold),
        "PROTO12 aggregate constraint-order decision differs",
    )
    return {
        "raw_owned_norms": {key: list(raw[key]) for key in CONSTRAINT_COMPONENTS},
        "accumulated_roundoff_enclosures": list(enclosure_values),
        "effective_norms_for_order_only": effective_by_component,
        "component_finest_pair_orders": recomputed_orders,
        "component_status": recomputed_status,
        "minimum_finite_finest_pair_order": minimum_finite,
        "components_below_frozen_minimum": below_threshold,
        "monotone_refinement_passed": monotone,
    }


def _validate_common_event_contract(
    event: Mapping[str, Any],
    *,
    amplitude: str,
    event_index: int,
    minimum_order: float,
    terminal: bool,
) -> dict[str, Any]:
    expected_time = event_index / 16.0
    assessment = event.get("assessment")
    _require(
        event.get("event_type") == COMMON_EVENT
        and event.get("amplitude") == amplitude
        and event.get("event_index") == event_index
        and isinstance(assessment, Mapping)
        and assessment.get("coordinate_time") == expected_time,
        "PROTO12 common-event identity or schedule differs",
    )
    temporal = assessment.get("temporal_spectral_admission")
    trapped = assessment.get("trapped_assessment")
    members = assessment.get("member_diagnostics")
    _require(
        assessment.get("PROTO12_SRC4_source_backend_applied") is True
        and assessment.get("PROTO12_pairwise_spectral_classifier_applied") is True
        and assessment.get("PROTO12_all_raw_ratios_retained") is True
        and assessment.get("PROTO12_resolution_ladder_enforced") is True
        and assessment.get("PROTO11_common_event_constraints_use_reference_map")
        is True
        and assessment.get("PROTO11_interior_q_reprojected") is False
        and assessment.get("qualified_trapped_common_event") is False
        and temporal
        == {
            "available": False,
            "admission_passed": False,
            "required_minimum_samples": 64,
        }
        and isinstance(trapped, Mapping)
        and trapped.get("trapped_sign_passed") is False
        and isinstance(members, Mapping)
        and set(members)
        == {
            f"{method}-{count}"
            for method in EXPECTED_METHODS
            for count in EXPECTED_POINT_COUNTS
        },
        "PROTO12 common-event premise or physical boundary differs",
    )
    scores = tuple(trapped.get("primary_trapped_scores", ())) + tuple(
        trapped.get("comparator_trapped_scores", ())
    )
    _require(
        len(scores) == 6
        and all(_finite("trapped score", value) < 0.0 for value in scores),
        "PROTO12 common-event trapped-score boundary differs",
    )
    cfl_counts: dict[str, int] = {}
    for member, diagnostic in members.items():
        _require(isinstance(diagnostic, Mapping), "PROTO12 member diagnostic differs")
        source_retries = diagnostic.get("source_retry_count")
        cfl_retries = diagnostic.get("CFL_retry_count")
        _require(
            source_retries == 0
            and isinstance(cfl_retries, int)
            and not isinstance(cfl_retries, bool)
            and cfl_retries >= 0,
            "PROTO12 source/CFL retry provenance differs",
        )
        cfl_counts[member] = cfl_retries

    method_records: dict[str, Any] = {}
    for method in EXPECTED_METHODS:
        common = _method_common_event(assessment, method)
        constraint = common.get("constraint_admission")
        spatial = common.get("spatial_spectral_admission")
        expected_admission = not (terminal and method == "SSPRK3")
        _require(
            common.get("method") == method
            and common.get("point_counts") == list(EXPECTED_POINT_COUNTS)
            and common.get("coordinate_time") == expected_time
            and common.get("SRC4_source_backend_bound") is True
            and common.get("reference_state_map_applied") is True
            and common.get("interior_q_reprojected") is False
            and common.get("admission_passed") is expected_admission
            and isinstance(constraint, Mapping)
            and constraint.get("admission_passed") is expected_admission
            and isinstance(spatial, Mapping)
            and spatial.get("admission_passed") is True
            and spatial.get("medium_and_fine_individual_budgets_passed") is True
            and spatial.get("every_profile_contracted") is True
            and spatial.get("every_tail_resolved_or_saturated") is True,
            f"PROTO12 {method} common-event admission differs",
        )
        order = _recompute_constraint_orders(
            constraint, minimum_order=minimum_order
        )
        method_records[method] = {
            "admission_passed": common["admission_passed"],
            "constraint_admission_passed": constraint["admission_passed"],
            "spatial_spectral_admission_passed": True,
            "constraint_order": order,
        }

    _require(
        trapped.get("both_common_event_admissions_passed") is (not terminal),
        "PROTO12 trapped-assessment admission composition differs",
    )
    return {
        "event_index": event_index,
        "coordinate_time": expected_time,
        "qualified_trapped_common_event": False,
        "finest_primary_trapped_score": trapped["primary_trapped_scores"][-1],
        "finest_comparator_trapped_score": trapped[
            "comparator_trapped_scores"
        ][-1],
        "CFL_retry_count_by_member": cfl_counts,
        "methods": method_records,
    }


def _diagnose_amplitude(
    events: Sequence[Mapping[str, Any]],
    amplitude_record: Mapping[str, Any],
    *,
    amplitude: str,
    minimum_order: float,
    trailing_event_count: int,
) -> dict[str, Any]:
    terminal_index, terminal_time = TERMINAL_EVENT[amplitude]
    _require(
        len(events) == terminal_index + 1
        and [item.get("event_index") for item in events]
        == list(range(terminal_index + 1)),
        f"PROTO12 amplitude {amplitude} event sequence differs",
    )
    records: list[dict[str, Any]] = []
    for index, event in enumerate(events):
        records.append(
            _validate_common_event_contract(
                event,
                amplitude=amplitude,
                event_index=index,
                minimum_order=minimum_order,
                terminal=index == terminal_index,
            )
        )

    stop = amplitude_record.get("stop")
    expected_stop = {
        "classification": "common_event_constraint_or_spectral_stop",
        "reason": "common_event_constraint_or_spatial_spectral_admission",
        "event_index": terminal_index,
        "coordinate_time": terminal_time,
    }
    source_counts = amplitude_record.get("member_source_retry_counts")
    _require(
        amplitude_record.get("amplitude") == amplitude
        and amplitude_record.get("eligible") is False
        and amplitude_record.get("last_completed_common_event_index")
        == terminal_index
        and amplitude_record.get("last_completed_coordinate_time") == terminal_time
        and amplitude_record.get("consecutive_qualified_trapped_common_events") == 0
        and stop == expected_stop
        and isinstance(source_counts, Mapping)
        and set(source_counts) == set(records[-1]["CFL_retry_count_by_member"])
        and all(value == 0 for value in source_counts.values()),
        f"PROTO12 amplitude {amplitude} terminal record differs",
    )

    terminal = records[-1]
    previous = records[-2]
    primary_constraint = terminal["methods"]["RK4"]["constraint_order"]
    comparator_constraint = terminal["methods"]["SSPRK3"]["constraint_order"]
    previous_comparator = previous["methods"]["SSPRK3"]["constraint_order"]
    radial_order = comparator_constraint["component_finest_pair_orders"][
        "radial_momentum"
    ]
    previous_radial_order = previous_comparator["component_finest_pair_orders"][
        "radial_momentum"
    ]
    _require(
        terminal["methods"]["RK4"]["admission_passed"] is True
        and terminal["methods"]["RK4"]["constraint_admission_passed"] is True
        and primary_constraint["minimum_finite_finest_pair_order"] >= minimum_order
        and terminal["methods"]["SSPRK3"]["admission_passed"] is False
        and terminal["methods"]["SSPRK3"]["constraint_admission_passed"] is False
        and terminal["methods"]["SSPRK3"]["spatial_spectral_admission_passed"]
        is True
        and comparator_constraint["components_below_frozen_minimum"]
        == ["radial_momentum"]
        and radial_order
        == comparator_constraint["minimum_finite_finest_pair_order"]
        and radial_order < minimum_order
        and previous["methods"]["SSPRK3"]["admission_passed"] is True
        and previous_radial_order >= minimum_order,
        f"PROTO12 amplitude {amplitude} terminal comparator veto differs",
    )

    terminal_constraint = _method_common_event(
        events[-1]["assessment"], "SSPRK3"
    )["constraint_admission"]
    _require(
        terminal_constraint.get("common_event_alignment_passed") is True
        and terminal_constraint.get("coarsest_guard_passed") is True
        and terminal_constraint.get("finest_guard_passed") is True
        and terminal_constraint.get("monotone_refinement_passed") is True
        and terminal_constraint.get("finest_pair_order_passed") is False,
        f"PROTO12 amplitude {amplitude} sole constraint gate differs",
    )

    trailing = records[-trailing_event_count:]
    trailing_orders = [
        item["methods"]["SSPRK3"]["constraint_order"]
        ["component_finest_pair_orders"]["radial_momentum"]
        for item in trailing
    ]
    _require(
        len(trailing_orders) == trailing_event_count
        and all(isinstance(value, float) for value in trailing_orders)
        and all(
            right < left
            for left, right in zip(trailing_orders[:-1], trailing_orders[1:], strict=True)
        ),
        f"PROTO12 amplitude {amplitude} late radial-order sequence differs",
    )
    terminal_cfl = terminal["CFL_retry_count_by_member"]
    _require(
        any(value > 0 for value in terminal_cfl.values()),
        f"PROTO12 amplitude {amplitude} CFL-retry provenance differs",
    )
    return {
        "amplitude": amplitude,
        "common_event_count": len(records),
        "last_fully_admitted_event_index": terminal_index - 1,
        "last_fully_admitted_coordinate_time": (terminal_index - 1) / 16.0,
        "terminal_event_index": terminal_index,
        "terminal_coordinate_time": terminal_time,
        "terminal_classification": stop["classification"],
        "terminal_reason": stop["reason"],
        "serialized_source_retry_count": 0,
        "terminal_CFL_retry_count_by_member": terminal_cfl,
        "some_adaptive_CFL_step_reductions_occurred": True,
        "campaign_stop_was_not_CFL_retry_exhaustion": True,
        "primary_constraint_minimum_order": primary_constraint[
            "minimum_finite_finest_pair_order"
        ],
        "comparator_frozen_minimum_order": minimum_order,
        "comparator_terminal_radial_momentum_order": radial_order,
        "comparator_order_shortfall": minimum_order - radial_order,
        "comparator_previous_radial_momentum_order": previous_radial_order,
        "comparator_only_component_below_minimum": "radial_momentum",
        "comparator_constraint_alignment_passed": True,
        "comparator_constraint_absolute_guards_passed": True,
        "comparator_constraint_monotone_refinement_passed": True,
        "both_PROTO12_spatial_spectral_admissions_passed": True,
        "terminal_comparator_radial_momentum_raw_owned_norms": (
            comparator_constraint["raw_owned_norms"]["radial_momentum"]
        ),
        "terminal_comparator_accumulated_roundoff_enclosures": (
            comparator_constraint["accumulated_roundoff_enclosures"]
        ),
        "terminal_comparator_radial_momentum_effective_norms_for_order_only": (
            comparator_constraint["effective_norms_for_order_only"][
                "radial_momentum"
            ]
        ),
        "terminal_comparator_component_orders": comparator_constraint[
            "component_finest_pair_orders"
        ],
        "trailing_comparator_radial_momentum_orders": trailing_orders,
        "trailing_event_indices": [item["event_index"] for item in trailing],
        "terminal_finest_primary_trapped_score": terminal[
            "finest_primary_trapped_score"
        ],
        "terminal_finest_comparator_trapped_score": terminal[
            "finest_comparator_trapped_score"
        ],
        "qualified_trapped_event_observed": False,
    }


def diagnose_proto12_campaign(
    events: Sequence[Mapping[str, Any]],
    campaign_result: Mapping[str, Any],
    *,
    minimum_constraint_order: Real = 1.5,
    trailing_event_count: int = 6,
) -> dict[str, Any]:
    """Return the complete outcome-neutral CAL9/PREF13 diagnosis."""

    minimum_order = _finite(
        "minimum_constraint_order", minimum_constraint_order, positive=True
    )
    if (
        isinstance(trailing_event_count, bool)
        or not isinstance(trailing_event_count, int)
        or trailing_event_count < 2
    ):
        raise ValueError("trailing event count must be an integer of at least two")
    records = tuple(events)
    amplitude_records = campaign_result.get("amplitude_records")
    _require(
        campaign_result.get("runner_id") == "FGC-1-CAL9-RUN1-RUNNER"
        and campaign_result.get("classification")
        == "calibration_failed_no_eligible_GR0_case"
        and campaign_result.get("selected_amplitude") is None
        and campaign_result.get("holdout_execution_authorized") is False
        and campaign_result.get("mechanism_question_answered") is False
        and campaign_result.get("retained_EFT_evolution_authorized") is False
        and campaign_result.get("SGBL_outcome_read") is False
        and campaign_result.get("FGCQR_outcome_read") is False
        and campaign_result.get("PROTO12_SRC4_source_backend_enabled") is True
        and campaign_result.get("PROTO12_pairwise_spectral_classifier_enabled")
        is True
        and campaign_result.get("PROTO12_all_raw_ratios_retained") is True
        and campaign_result.get("PROTO11_reference_state_map_enabled") is True
        and campaign_result.get("PROTO11_interior_q_reprojection_enabled") is False
        and isinstance(amplitude_records, list)
        and [item.get("amplitude") for item in amplitude_records]
        == list(EXPECTED_AMPLITUDES),
        "PROTO12 campaign result boundary differs",
    )
    _require(
        len(records) == 49
        and all(item.get("event_type") == COMMON_EVENT for item in records)
        and [item.get("amplitude") for item in records]
        == ["5/2"] * 25 + ["3"] * 24,
        "PROTO12 event partition or amplitude order differs",
    )

    amplitude_diagnoses: dict[str, Any] = {}
    offset = 0
    for amplitude, amplitude_record in zip(
        EXPECTED_AMPLITUDES, amplitude_records, strict=True
    ):
        count = TERMINAL_EVENT[amplitude][0] + 1
        amplitude_diagnoses[amplitude] = _diagnose_amplitude(
            records[offset : offset + count],
            amplitude_record,
            amplitude=amplitude,
            minimum_order=minimum_order,
            trailing_event_count=trailing_event_count,
        )
        offset += count
    _require(offset == len(records), "PROTO12 event consumption differs")

    cfl_counts = {
        amplitude: diagnosis["terminal_CFL_retry_count_by_member"]
        for amplitude, diagnosis in amplitude_diagnoses.items()
    }
    return {
        "event_log_line_count": len(records),
        "common_event_count": len(records),
        "serialized_source_retry_count": 0,
        "some_adaptive_CFL_step_reductions_occurred": True,
        "campaign_stop_was_not_CFL_retry_exhaustion": True,
        "terminal_CFL_retry_count_by_amplitude_and_member": cfl_counts,
        "both_amplitudes_reached_late_evolved_common_events": True,
        "amplitudes": amplitude_diagnoses,
        "shared_declared_veto": {
            "method": "SSPRK3",
            "constraint_component": "radial_momentum",
            "gate": "finest_pair_order_passed",
            "frozen_minimum_order": minimum_order,
            "same_configured_gate_failed_for_both_amplitudes": True,
            "all_other_terminal_comparator_constraint_subgates_passed": True,
            "both_terminal_PROTO12_spatial_admissions_passed": True,
            "common_physical_cause_derived": False,
            "common_numerical_cause_derived": False,
            "preasymptotic_or_persistent_decided": False,
        },
        "SRC4_arithmetic_wall_cleared_for_both_amplitudes": True,
        "no_GR0_amplitude_eligible": True,
        "candidate_or_mechanism_outcome_read": False,
    }


__all__ = [
    "COMMON_EVENT",
    "CONSTRAINT_COMPONENTS",
    "EXPECTED_AMPLITUDES",
    "EXPECTED_METHODS",
    "EXPECTED_POINT_COUNTS",
    "FINITE_CONSTRAINT_COMPONENTS",
    "TERMINAL_EVENT",
    "diagnose_proto12_campaign",
]
