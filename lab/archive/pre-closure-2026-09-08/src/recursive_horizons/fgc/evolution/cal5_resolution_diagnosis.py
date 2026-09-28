"""Outcome-neutral diagnosis of the PROTO8 GR-0 calibration stops.

The PROTO8 campaign is calibration data, not a candidate or physical result.
This module reduces already committed common-event records and asks one
bounded numerical question: is a one-level refinement of the *entire* nested
grid ladder a discriminating successor experiment while every admission
threshold remains fixed?

No trajectory is advanced here.  Conditional Richardson projections are
reported only as planning diagnostics; they never count as admission.
"""

from __future__ import annotations

from math import isfinite
from numbers import Real
from typing import Any, Mapping, Sequence


PROTO8_POINT_COUNTS = (1025, 2049, 4097)
PROSPECTIVE_POINT_COUNTS = (2049, 4097, 8193)
MINIMUM_FINEST_PAIR_ORDER = 1.5
MAXIMUM_NESTED_TAIL_RATIO = 0.25
METHOD_GUARDS = {
    "RK4": {"coarsest": 1.0 / 50.0, "finest": 1.0 / 1000.0},
    "SSPRK3": {"coarsest": 1.0 / 10.0, "finest": 1.0 / 200.0},
}
METHOD_EVENT_KEYS = {
    "RK4": "primary_common_event",
    "SSPRK3": "comparator_common_event",
}


def _finite(name: str, value: Real, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _bool(name: str, value: object) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be Boolean")
    return value


def _mapping(name: str, value: object) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


def _tail_failures(
    spectral: Mapping[str, Any],
) -> tuple[dict[str, object], ...]:
    failures: list[dict[str, object]] = []
    for family, key in (
        ("field", "field_power_tail_ratios"),
        ("derivative", "derivative_power_tail_ratios"),
    ):
        ratios = _mapping(key, spectral.get(key))
        for field in sorted(ratios):
            values = tuple(ratios[field])
            if len(values) != 2:
                raise ValueError(f"{key}.{field} must contain two adjacent ratios")
            numeric = tuple(
                _finite(f"{key}.{field}[{index}]", value, positive=True)
                for index, value in enumerate(values)
            )
            for pair_index, value in enumerate(numeric):
                if value > MAXIMUM_NESTED_TAIL_RATIO:
                    failures.append(
                        {
                            "family": family,
                            "field": field,
                            "pair": (
                                "1025_to_2049"
                                if pair_index == 0
                                else "2049_to_4097"
                            ),
                            "ratio": value,
                            "maximum": MAXIMUM_NESTED_TAIL_RATIO,
                        }
                    )
    return tuple(failures)


def diagnose_method_common_event(
    assessment: Mapping[str, Any],
    *,
    method: str,
) -> dict[str, Any]:
    """Reduce one method at one committed PROTO8 common event."""

    if method not in METHOD_EVENT_KEYS:
        raise ValueError("method must be RK4 or SSPRK3")
    event = _mapping(METHOD_EVENT_KEYS[method], assessment.get(METHOD_EVENT_KEYS[method]))
    if event.get("method") != method:
        raise ValueError("common-event method label differs")
    point_counts = tuple(event.get("point_counts", ()))
    if point_counts != PROTO8_POINT_COUNTS:
        raise ValueError("PROTO8 point-count ladder differs")
    constraint = _mapping("constraint_admission", event.get("constraint_admission"))
    if tuple(constraint.get("point_counts", ())) != PROTO8_POINT_COUNTS:
        raise ValueError("constraint point-count ladder differs")
    norms = tuple(
        _finite(f"owned_raw_global_norms[{index}]", value, positive=True)
        for index, value in enumerate(constraint.get("owned_raw_global_norms", ()))
    )
    if len(norms) != 3:
        raise ValueError("owned raw global norms must contain three grids")
    if not (norms[0] > norms[1] > norms[2]):
        raise ValueError("owned raw constraint norms do not contract monotonically")
    order = _finite(
        "minimum_finite_finest_pair_order",
        constraint.get("minimum_finite_finest_pair_order"),
        positive=True,
    )
    guards = METHOD_GUARDS[method]
    prospective_coarse = norms[1]
    current_fine = norms[2]
    conditional_next = current_fine / (2.0**order)
    spectral = _mapping(
        "spatial_spectral_admission", event.get("spatial_spectral_admission")
    )
    failures = _tail_failures(spectral)
    budgets = _mapping(
        "individual_budget_passed_by_grid",
        spectral.get("individual_budget_passed_by_grid"),
    )
    medium = _mapping("medium spectral budget", budgets.get("2049"))
    fine = _mapping("fine spectral budget", budgets.get("4097"))
    medium_all = bool(medium) and all(_bool("medium budget", item) for item in medium.values())
    fine_all = bool(fine) and all(_bool("fine budget", item) for item in fine.values())

    return {
        "method": method,
        "coordinate_time": _finite("coordinate_time", event.get("coordinate_time")),
        "PROTO8_admission_passed": _bool(
            "admission_passed", event.get("admission_passed")
        ),
        "constraint_admission_passed": _bool(
            "constraint admission", constraint.get("admission_passed")
        ),
        "constraint_owned_raw_global_norms": list(norms),
        "constraint_minimum_finite_finest_pair_order": order,
        "constraint_order_threshold": MINIMUM_FINEST_PAIR_ORDER,
        "constraint_order_passed": _bool(
            "finest_pair_order_passed", constraint.get("finest_pair_order_passed")
        ),
        "current_coarsest_guard": guards["coarsest"],
        "current_finest_guard": guards["finest"],
        "current_coarsest_guard_passed": _bool(
            "coarsest_guard_passed", constraint.get("coarsest_guard_passed")
        ),
        "current_finest_guard_passed": _bool(
            "finest_guard_passed", constraint.get("finest_guard_passed")
        ),
        "prospective_2049_coarsest_norm": prospective_coarse,
        "prospective_2049_coarsest_guard_passed": (
            prospective_coarse <= guards["coarsest"]
        ),
        "conditional_8193_norm_if_observed_order_persists": conditional_next,
        "conditional_8193_finest_guard_passed": conditional_next <= guards["finest"],
        "conditional_projection_is_admission_evidence": False,
        "medium_and_fine_absolute_spectral_budgets_passed": medium_all and fine_all,
        "current_nested_tail_failures": [dict(item) for item in failures],
        "new_4097_to_8193_tail_is_unmeasured": True,
    }


def diagnose_stop_event(event: Mapping[str, Any]) -> dict[str, Any]:
    """Reduce one failed amplitude stop without changing its decision."""

    if event.get("event_type") != "common_event":
        raise ValueError("stop event must be a common event")
    amplitude = event.get("amplitude")
    if amplitude not in {"5/2", "3"}:
        raise ValueError("unexpected calibration amplitude")
    assessment = _mapping("assessment", event.get("assessment"))
    methods = {
        method: diagnose_method_common_event(assessment, method=method)
        for method in METHOD_EVENT_KEYS
    }
    if any(item["PROTO8_admission_passed"] for item in methods.values()):
        raise ValueError("diagnosed stop unexpectedly passed a method")
    if any(
        not item["prospective_2049_coarsest_guard_passed"]
        or not item["conditional_8193_finest_guard_passed"]
        or not item["medium_and_fine_absolute_spectral_budgets_passed"]
        for item in methods.values()
    ):
        raise ValueError("one-level resolution shift lacks the declared premise evidence")
    return {
        "amplitude": amplitude,
        "event_index": event.get("event_index"),
        "coordinate_time": methods["RK4"]["coordinate_time"],
        "methods": methods,
        "one_level_resolution_shift_is_a_discriminating_test": True,
        "one_level_resolution_shift_is_already_admitted": False,
    }


def diagnose_proto8_campaign_events(
    events: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Validate the public event sequence and reduce both terminal stops."""

    records = tuple(events)
    if len(records) != 9:
        raise ValueError("PROTO8 event log must contain exactly nine records")
    common = [item for item in records if item.get("event_type") == "common_event"]
    retries = [
        item
        for item in records
        if item.get("event_type") == "rejected_unaccepted_source_only_proposal"
    ]
    observed = [
        (
            item.get("amplitude"),
            item.get("event_index"),
            _mapping("assessment", item.get("assessment"))
            .get("primary_common_event", {})
            .get("coordinate_time"),
        )
        for item in common
    ]
    expected = [
        ("5/2", 0, 0.0),
        ("5/2", 1, 1.0 / 16.0),
        ("5/2", 2, 1.0 / 8.0),
        ("3", 0, 0.0),
        ("3", 1, 1.0 / 16.0),
    ]
    if observed != expected:
        raise ValueError("PROTO8 common-event sequence differs")
    if len(retries) != 4:
        raise ValueError("PROTO8 source-only retry count differs")
    for amplitude in ("5/2", "3"):
        owned = [item for item in retries if item.get("amplitude") == amplitude]
        if (
            len(owned) != 2
            or any(item.get("method") != "RK4" for item in owned)
            or any(item.get("point_count") != 4097 for item in owned)
            or [item.get("retry_count_for_accepted_step") for item in owned] != [1, 2]
            or any(item.get("retryable_source_only") is not True for item in owned)
            or any(item.get("non_source_failure_vetoed_retry") is not False for item in owned)
            or any(item.get("fields_bitwise_preserved") is not True for item in owned)
        ):
            raise ValueError("PROTO8 source-only retry evidence differs")

    stop_5 = diagnose_stop_event(common[2])
    stop_3 = diagnose_stop_event(common[4])
    if stop_5["methods"]["RK4"]["current_nested_tail_failures"]:
        raise ValueError("amplitude 5/2 unexpectedly has a nested-tail failure")
    if stop_5["methods"]["SSPRK3"]["current_nested_tail_failures"]:
        raise ValueError("amplitude 5/2 unexpectedly has a nested-tail failure")
    for method in METHOD_EVENT_KEYS:
        failures = stop_3["methods"][method]["current_nested_tail_failures"]
        if failures != [
            {
                "family": "derivative",
                "field": "phi_over_Lambda",
                "pair": "2049_to_4097",
                "ratio": failures[0]["ratio"] if failures else None,
                "maximum": MAXIMUM_NESTED_TAIL_RATIO,
            }
        ]:
            raise ValueError("amplitude 3 nested-tail failure differs")

    return {
        "public_event_count": len(records),
        "common_event_count": len(common),
        "source_only_retry_count": len(retries),
        "amplitude_5over2_stop": stop_5,
        "amplitude_3_stop": stop_3,
        "current_point_counts": list(PROTO8_POINT_COUNTS),
        "smallest_prospective_point_counts": list(PROSPECTIVE_POINT_COUNTS),
        "threshold_changes_required": False,
        "new_8193_data_required_before_admission": True,
        "candidate_or_physical_trajectory_read": False,
    }


__all__ = [
    "MAXIMUM_NESTED_TAIL_RATIO",
    "METHOD_GUARDS",
    "MINIMUM_FINEST_PAIR_ORDER",
    "PROTO8_POINT_COUNTS",
    "PROSPECTIVE_POINT_COUNTS",
    "diagnose_method_common_event",
    "diagnose_proto8_campaign_events",
    "diagnose_stop_event",
]
