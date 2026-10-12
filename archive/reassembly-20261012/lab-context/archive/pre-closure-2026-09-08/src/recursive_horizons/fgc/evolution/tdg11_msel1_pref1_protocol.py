"""Independent, I/O-free validation of TDG11-MSEL1 evidence records.

The wire format is shared with the frozen diagnostic, not its implementation.
This module imports no runner, reconstruction, localizer adapter, NumPy, Git,
store, or publication decision. It recomputes the rational predicates and
selection from saved enclosures. Actual reconstruction belongs to the live
PREF1 binder; validating these compact records alone is not a replay.
"""

from __future__ import annotations

from dataclasses import fields, is_dataclass
from fractions import Fraction
from hashlib import sha256
import math
import re

from . import tdg11_msel1_contract as freeze


RAW_SCHEMA = "FGC-1-TDG11-MSEL1-raw-v1"
RUNNER_ID = "FGC-1-TDG11-MSEL1-RUN1"
SELECTED = "qualified_candidate_selected_no_state_advance"
NONPASS = "all_three_candidates_nonpass_no_state_advance"
INCONCLUSIVE = "bounded_selection_inconclusive_no_state_advance"
CHANNELS = freeze.TDG6_COMPLETE_STATE_CHANNELS
CANDIDATES = freeze.CANDIDATES
GROUP_STATUSES = frozenset(
    {"not_attempted", "complete", "premise_stop", "resource_stop"}
)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_INTEGER_TEXT = re.compile(r"-?(?:0|[1-9][0-9]*)\Z", re.ASCII)
_MAX_INTEGER_BITS = 32768
_MAX_INTEGER_TEXT = 10000
_MAX_DEPTH = 48


class PREF1ProtocolError(ValueError):
    """Malformed or internally inconsistent evidence, never a nonpass."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PREF1ProtocolError(message)


def exact_keys(value: object, names: set[str], label: str) -> dict:
    require(type(value) is dict and set(value) == names, f"{label}: field set differs")
    return value


def strict_equal(actual: object, expected: object, label: str = "evidence") -> None:
    """Compare the typed JSON tree without bool/int or list/tuple aliases."""
    require(type(actual) is type(expected), f"{label}: type differs")
    if type(expected) is dict:
        exact_keys(actual, set(expected), label)
        for key in expected:
            strict_equal(actual[key], expected[key], f"{label}.{key}")
    elif type(expected) is list:
        require(len(actual) == len(expected), f"{label}: length differs")
        for index, (left, right) in enumerate(zip(actual, expected, strict=True)):
            strict_equal(left, right, f"{label}[{index}]")
    else:
        require(actual == expected, f"{label}: value differs")


def integer(value: object, label: str, *, minimum: int = 0) -> int:
    require(type(value) is int and value >= minimum, f"{label}: integer required")
    return value


def digest(value: object, label: str) -> str:
    require(
        type(value) is str and _SHA256.fullmatch(value) is not None,
        f"{label}: SHA-256 required",
    )
    return value


def integer_text(value: int) -> str:
    """Encode bounded integers without changing Python's global digit limit."""
    require(
        type(value) is int and abs(value).bit_length() <= _MAX_INTEGER_BITS,
        "integer exceeds exact serialization domain",
    )
    if value == 0:
        return "0"
    remaining = abs(value)
    groups = []
    while remaining:
        remaining, group = divmod(remaining, 1_000_000_000)
        groups.append(group)
    digits = str(groups.pop())
    digits += "".join(format(group, "09d") for group in reversed(groups))
    return ("-" if value < 0 else "") + digits


def parse_integer_text(value: object) -> int:
    require(
        type(value) is str
        and len(value) <= _MAX_INTEGER_TEXT
        and value != "-0"
        and _INTEGER_TEXT.fullmatch(value) is not None,
        "noncanonical exact integer text",
    )
    negative = value.startswith("-")
    digits = value[1:] if negative else value
    result = 0
    # Independent left-to-right radix reduction; no unbounded int(text).
    for start in range(0, len(digits), 8):
        part = digits[start : start + 8]
        result = result * (10 ** len(part)) + int(part)
        require(result.bit_length() <= _MAX_INTEGER_BITS, "exact integer bit ceiling")
    return -result if negative else result


def fraction(value: object, *, nonnegative: bool = True) -> Fraction:
    item = exact_keys(value, {"numerator", "denominator"}, "rational")
    numerator = parse_integer_text(item["numerator"])
    denominator = parse_integer_text(item["denominator"])
    require(denominator > 0, "nonpositive rational denominator")
    result = Fraction(numerator, denominator)
    require(
        (result.numerator, result.denominator) == (numerator, denominator),
        "unreduced rational",
    )
    require(not nonnegative or result >= 0, "negative magnitude")
    return result


def wire(value: object, depth: int = 0) -> object:
    """Serialize binder-owned values to the frozen exact-rational wire shape."""
    require(depth <= _MAX_DEPTH, "exact record nesting ceiling")
    if type(value) is Fraction:
        return {
            "numerator": integer_text(value.numerator),
            "denominator": integer_text(value.denominator),
        }
    if is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: wire(getattr(value, field.name), depth + 1)
            for field in fields(value)
        }
    if type(value) is dict:
        require(all(type(key) is str for key in value), "non-string record key")
        return {key: wire(item, depth + 1) for key, item in value.items()}
    if type(value) in (tuple, list):
        return [wire(item, depth + 1) for item in value]
    if type(value) is float:
        require(math.isfinite(value), "nonfinite binary64 fingerprint")
        return {"binary64_hex": value.hex()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise PREF1ProtocolError(f"unsupported evidence type: {type(value).__name__}")


def decision_from_bounds(
    l01: Fraction, u01: Fraction, l12: Fraction, u12: Fraction
) -> dict:
    """Independent transcription of the unchanged, five-case TDG6 law."""
    require(
        all(type(item) is Fraction and item >= 0 for item in (l01, u01, l12, u12)),
        "decision needs exact magnitudes",
    )
    require(l01 <= u01 and l12 <= u12, "reversed magnitude interval")
    resolved = False
    order_passed = None
    if u01 == 0 and u12 == 0:
        classification = "exact_zero"
        admitted = True
    elif l01 == 0 and l12 == 0 and u12 <= u01:
        classification = "enclosure_dominated_debit_only"
        admitted = True
    elif l01 > 0:
        resolved = True
        order_passed = 8 * u12 * u12 <= l01 * l01
        classification = (
            "resolved_order_pass" if order_passed else "resolved_order_failure"
        )
        admitted = order_passed
    else:
        classification = "order_inconclusive"
        admitted = False
    return {
        "classification": classification,
        "admission_passed": admitted,
        "temporal_retry_permitted": not admitted,
        "minimum_observed_order": Fraction(3, 2),
        "outer_difference": {"lower": l01, "upper": u01},
        "finest_difference": {"lower": l12, "upper": u12},
        "finest_pair_debit": u12,
        "order_threshold_resolved": resolved,
        "order_threshold_passed": order_passed,
        "absolute_state_tolerance_used": False,
        "physical_signal_used_for_normalization": False,
    }


def validate_decision(value: object, bounds: tuple[Fraction, ...]) -> dict:
    expected = decision_from_bounds(*bounds)
    strict_equal(value, wire(expected), "TDG6 decision")
    return expected


def three_magnitudes(value: object, label: str) -> tuple[Fraction, Fraction, Fraction]:
    require(type(value) is list and len(value) == 3, f"{label}: three levels required")
    return tuple(fraction(item) for item in value)


def typed_stop(value: object, label: str = "stop") -> dict:
    stop = exact_keys(value, {"owner", "code", "detail"}, label)
    require(
        all(type(item) is str and 0 < len(item) <= 640 for item in stop.values()),
        f"{label}: bounded nonempty text required",
    )
    return stop


def validate_difference(
    value: object, *, polynomials: int, ceiling: int, label: str
) -> tuple[Fraction, Fraction]:
    item = exact_keys(
        value,
        {
            "lower",
            "upper",
            "polynomial_count",
            "candidate_count",
            "survivor_count",
            "localization_classification",
            "coefficient_stream_sha256",
            "survivor_key_stream_sha256",
            "primary_stationary_count_stream_sha256",
            "independent_stationary_count_stream_sha256",
            "primary_evaluator_id",
            "independent_evaluator_id",
            "maximum_candidates",
            "routes_agree",
        },
        label,
    )
    lower, upper = fraction(item["lower"]), fraction(item["upper"])
    require(lower <= upper, f"{label}: reversed enclosure")
    strict_equal(item["polynomial_count"], polynomials, f"{label}.polynomials")
    strict_equal(item["maximum_candidates"], ceiling, f"{label}.ceiling")
    count = integer(item["candidate_count"], f"{label}.candidates", minimum=1)
    survivors = integer(item["survivor_count"], f"{label}.survivors", minimum=1)
    require(
        2 * polynomials <= count <= min(4 * polynomials, ceiling),
        f"{label}: impossible cubic candidate count",
    )
    require(survivors <= count, f"{label}: excess survivors")
    if item["localization_classification"] == "unique_maximum":
        require(survivors == 1, f"{label}: unique maximum has multiple survivors")
    else:
        strict_equal(
            item["localization_classification"],
            "nonunique_or_interval_inconclusive",
            f"{label}.localization",
        )
        require(survivors >= 2, f"{label}: missing co-survivors")
    for key in (
        "coefficient_stream_sha256",
        "survivor_key_stream_sha256",
        "primary_stationary_count_stream_sha256",
        "independent_stationary_count_stream_sha256",
    ):
        digest(item[key], f"{label}.{key}")
    strict_equal(
        item["primary_stationary_count_stream_sha256"],
        item["independent_stationary_count_stream_sha256"],
        f"{label}.stationary-count agreement",
    )
    strict_equal(
        item["primary_evaluator_id"],
        "tdg9_loc1_derivative_monotone_bisection_v1",
        f"{label}.primary",
    )
    strict_equal(
        item["independent_evaluator_id"],
        "tdg9_loc2_endpoint_deflated_adaptive_discriminant_v2",
        f"{label}.independent",
    )
    strict_equal(item["routes_agree"], True, f"{label}.routes_agree")
    return lower, upper


def validate_complete_c(value: object, *, expected_rows: int) -> dict:
    integer(expected_rows, "expected rows", minimum=1)
    data = exact_keys(
        value,
        {
            "d01",
            "d12",
            "decision",
            "row_count",
            "expected_row_count",
            "row_stream_sha256",
            "combined_coefficient_stream_sha256",
            "sufficient_pass_left",
            "sufficient_pass_right",
            "sufficient_contraction_pass",
            "sufficient_contraction_failure",
            "threshold_inconclusive",
            "maximum_candidates_D01",
            "maximum_candidates_D12",
            "refinement_depth",
            "schema_version",
            "evaluator_id",
            "absolute_tolerance_used",
            "physical_signal_used_for_normalization",
            "declared_cubic_is_exact_PDE_history",
            "production_authority",
        },
        "complete-C",
    )
    for key in ("row_count", "expected_row_count"):
        strict_equal(data[key], expected_rows, key)
    ceilings = (
        freeze.RESOURCES["per_channel_maximum_candidates_D01"],
        freeze.RESOURCES["per_channel_maximum_candidates_D12"],
    )
    l01, u01 = validate_difference(
        data["d01"], polynomials=2 * expected_rows, ceiling=ceilings[0], label="D01"
    )
    l12, u12 = validate_difference(
        data["d12"], polynomials=4 * expected_rows, ceiling=ceilings[1], label="D12"
    )
    strict_equal(data["maximum_candidates_D01"], ceilings[0], "D01 ceiling")
    strict_equal(data["maximum_candidates_D12"], ceilings[1], "D12 ceiling")
    strict_equal(
        data["refinement_depth"],
        freeze.RESOURCES["primary_refinement_depth"],
        "root isolation depth",
    )
    digest(data["row_stream_sha256"], "row stream")
    combined = sha256(
        b"TDG11-RATIONAL-COMPLETE-C-COMBINED-v1\n"
        + (
            data["d01"]["coefficient_stream_sha256"]
            + "\n"
            + data["d12"]["coefficient_stream_sha256"]
            + "\n"
        ).encode("ascii")
    )
    strict_equal(
        data["combined_coefficient_stream_sha256"],
        combined.hexdigest(),
        "combined coefficient stream",
    )
    left, right = 8 * u12 * u12, l01 * l01
    require(
        fraction(data["sufficient_pass_left"]) == left
        and fraction(data["sufficient_pass_right"]) == right,
        "complete-C inequality witness differs",
    )
    zero = u01 == 0 and u12 == 0
    passed = not zero and left <= right
    failed = not zero and not passed and 8 * l12 * l12 > u01 * u01
    unresolved = not zero and not passed and not failed
    strict_equal(data["sufficient_contraction_pass"], passed, "sufficient pass")
    strict_equal(data["sufficient_contraction_failure"], failed, "sufficient failure")
    strict_equal(data["threshold_inconclusive"], unresolved, "threshold inconclusive")
    decision = validate_decision(data["decision"], (l01, u01, l12, u12))
    strict_equal(data["schema_version"], 1, "complete-C schema")
    strict_equal(
        data["evaluator_id"],
        "tdg11_rational_complete_c_dual_route_v1",
        "complete-C evaluator",
    )
    for key in (
        "absolute_tolerance_used",
        "physical_signal_used_for_normalization",
        "declared_cubic_is_exact_PDE_history",
        "production_authority",
    ):
        strict_equal(data[key], False, f"complete-C.{key}")
    return {
        "d01": (l01, u01),
        "d12": (l12, u12),
        "decision": decision,
        "sufficient_failure": failed,
        "threshold_inconclusive": unresolved,
    }


def classify_channel(value: object, candidate: str, *, expected_rows: int) -> str:
    """Recompute a channel outcome; a failed sufficient pass is not a failure."""
    require(
        type(candidate) is str and candidate in ("baseline", *CANDIDATES),
        "unknown candidate",
    )
    if candidate == CANDIDATES[2]:
        data = exact_keys(
            value,
            {
                "channel",
                "estimators",
                "roundoff_bounds",
                "contraction_01",
                "contraction_12",
                "public_fine_debit",
                "global_PDE_enclosure",
            },
            "embedded channel",
        )
        e0, e1, e2 = three_magnitudes(data["estimators"], "embedded estimators")
        rounds = three_magnitudes(data["roundoff_bounds"], "accumulation bounds")
        first = validate_decision(data["contraction_01"], (e0, e0, e1, e1))
        second = validate_decision(data["contraction_12"], (e1, e1, e2, e2))
        strict_equal(data["global_PDE_enclosure"], False, "embedded nonclaim")
        require(
            fraction(data["public_fine_debit"]) == e2 + rounds[2],
            "embedded debit omits an error contribution",
        )
        if first["admission_passed"] and second["admission_passed"]:
            outcome = "pass"
        else:
            outcome = (
                "nonpass"
                if 8 * e1 * e1 > e0 * e0 or 8 * e2 * e2 > e1 * e1
                else "inconclusive"
            )
    else:
        names = {"channel", "exact_complete_C", "public_fine_debit"}
        if candidate == CANDIDATES[0]:
            names.add("roundoff_bounds")
        data = exact_keys(value, names, "complete-C channel")
        evidence = validate_complete_c(
            data["exact_complete_C"], expected_rows=expected_rows
        )
        debit = evidence["d12"][1]
        if candidate == CANDIDATES[0]:
            debit += three_magnitudes(data["roundoff_bounds"], "accumulation bounds")[2]
        require(
            fraction(data["public_fine_debit"]) == debit, "complete-C debit differs"
        )
        if evidence["decision"]["admission_passed"]:
            outcome = "pass"
        else:
            outcome = "nonpass" if evidence["sufficient_failure"] else "inconclusive"
    require(
        type(data["channel"]) is str and data["channel"] in CHANNELS,
        "unknown complete-state channel",
    )
    return outcome


def classify_group(value: object, candidate: str, *, expected_rows: int) -> list[str]:
    group = exact_keys(value, {"status", "channels", "stop"}, "group")
    require(
        type(group["status"]) is str and group["status"] in GROUP_STATUSES,
        "unknown group status",
    )
    records = group["channels"]
    require(
        type(records) is list and len(records) <= len(CHANNELS), "invalid channel count"
    )
    if group["status"] == "complete":
        require(len(records) == len(CHANNELS), "completed group is partial")
    if group["status"] == "not_attempted":
        require(not records, "unattempted group contains measurements")
    if group["status"] in {"complete", "not_attempted"}:
        strict_equal(group["stop"], None, "group stop")
    else:
        typed_stop(group["stop"])
    outcomes = []
    for index, record in enumerate(records):
        require(type(record) is dict, "channel is not a record")
        strict_equal(record.get("channel"), CHANNELS[index], "channel order")
        outcomes.append(
            classify_channel(record, candidate, expected_rows=expected_rows)
        )
    return outcomes


_WIDTH_FIELDS = {
    "retry",
    "generation",
    "width_hex",
    "restored",
    "families",
    "baseline",
    "candidates",
}
_COUNT_FIELDS = {"accepted_source_prechecks", "stage_and_endpoint_records", "rhs_calls"}
_FAMILY_LIMITS = {
    "accepted_source_prechecks": 7,
    "stage_and_endpoint_records": 35,
    "rhs_calls": 42,
}


def reduce_widths(value: object, *, expected_rows: int) -> dict:
    """Reduce authenticated slots logically; this function grants no authority."""
    integer(expected_rows, "expected rows", minimum=1)
    require(
        type(value) is list and len(value) == 3, "exactly three width slots required"
    )
    rows = {candidate: [] for candidate in CANDIDATES}
    completed = {candidate: True for candidate in CANDIDATES}
    for width, spec in zip(value, freeze.replay_specs(), strict=True):
        exact_keys(width, _WIDTH_FIELDS, "width")
        for key in ("retry", "generation", "width_hex"):
            strict_equal(width[key], spec[key], f"width.{key}")
        exact_keys(width["candidates"], set(CANDIDATES), "candidate slots")
        classify_group(width["baseline"], "baseline", expected_rows=expected_rows)
        for candidate in CANDIDATES:
            group = width["candidates"][candidate]
            rows[candidate].extend(
                classify_group(group, candidate, expected_rows=expected_rows)
            )
            completed[candidate] = (
                completed[candidate] and group["status"] == "complete"
            )
    summaries = []
    for candidate in CANDIDATES:
        outcomes = rows[candidate]
        if completed[candidate] and len(outcomes) == 54 and set(outcomes) == {"pass"}:
            status = "pass"
        elif "nonpass" in outcomes:
            # One valid counterexample refutes an all-width/all-channel claim.
            status = "nonpass"
        else:
            status = "inconclusive"
        summaries.append(
            {
                "candidate": candidate,
                "status": status,
                "classified_channels": len(outcomes),
                "complete_all_widths": completed[candidate],
            }
        )
    passing = [item["candidate"] for item in summaries if item["status"] == "pass"]
    selected = passing[0] if passing else None
    if selected is not None:
        classification = SELECTED
    elif all(item["status"] == "nonpass" for item in summaries):
        classification = NONPASS
    else:
        classification = INCONCLUSIVE
    return {
        "classification": classification,
        "selected_candidate": selected,
        "candidate_summaries": summaries,
        "licenses_only_separate_TDG11_IMP1": False,
        "independent_PREF1_required_before_IMP1": True,
    }


def _family_groups(width: dict, arithmetic: str) -> list[dict]:
    if arithmetic == freeze.ORIGINAL_ARITHMETIC_ID:
        return [
            width["baseline"],
            width["candidates"][CANDIDATES[0]],
            width["candidates"][CANDIDATES[2]],
        ]
    return [width["candidates"][CANDIDATES[1]]]


def _validate_families(
    widths: list, accounting: object, *, global_stop: object
) -> None:
    totals = dict.fromkeys(_COUNT_FIELDS, 0)
    attempted_family = False
    for width, spec in zip(widths, freeze.replay_specs(), strict=True):
        exact_keys(width, _WIDTH_FIELDS, "width")
        families = exact_keys(
            width["families"],
            {freeze.ORIGINAL_ARITHMETIC_ID, freeze.COMPENSATED_ARITHMETIC_ID},
            "arithmetic families",
        )
        restored = width["restored"]
        if restored is not None:
            exact_keys(
                restored,
                {
                    "checkpoint_sha256",
                    "descriptor_sha256",
                    "state_sha256",
                    "fingerprint_sha256",
                },
                "restoration",
            )
            strict_equal(
                restored["checkpoint_sha256"], spec["checkpoint_sha256"], "checkpoint"
            )
            strict_equal(
                restored["descriptor_sha256"], freeze.DESCRIPTOR_SHA256, "descriptor"
            )
            strict_equal(
                restored["state_sha256"], freeze.PHYSICAL_STATE_SHA256, "state"
            )
            digest(restored["fingerprint_sha256"], "restored fingerprint")
        for arithmetic, family in families.items():
            groups = _family_groups(width, arithmetic)
            if family is None:
                require(
                    all(group["status"] == "not_attempted" for group in groups),
                    "measurements without a guarded family",
                )
                continue
            require(
                restored is not None and type(family) is dict,
                "guarded family lacks a restored predecessor",
            )
            attempted_family = True
            status = family.get("status")
            require(
                type(status) is str
                and status in {"complete", "premise_stop", "resource_stop"},
                "guarded family status differs",
            )
            complete = status == "complete"
            exact_keys(
                family,
                _COUNT_FIELDS | {"status", "family_sha256" if complete else "stop"},
                "guarded family",
            )
            for key, limit in _FAMILY_LIMITS.items():
                count = integer(family[key], f"family.{key}")
                require(count <= limit, "guarded family exceeded its budget")
                if complete:
                    strict_equal(count, limit, f"complete family.{key}")
                totals[key] += count
            if complete:
                digest(family["family_sha256"], "recorded family")
                candidate_order = (
                    ("baseline", CANDIDATES[0], CANDIDATES[2])
                    if arithmetic == freeze.ORIGINAL_ARITHMETIC_ID
                    else (CANDIDATES[1],)
                )
                for candidate, group in zip(candidate_order, groups, strict=True):
                    require(
                        group["status"] != "premise_stop",
                        "completed family cannot contain a premise-stop group",
                    )
                    if group["status"] == "resource_stop":
                        stop = typed_stop(group["stop"])
                        require(
                            stop["owner"] in {"exact_localizer", "resource"},
                            "completed-family stop has an impossible owner",
                        )
                        if stop["owner"] == "exact_localizer":
                            require(
                                candidate != CANDIDATES[2],
                                "embedded estimator never invokes an exact localizer",
                            )
                            require(
                                stop["code"]
                                in {
                                    "candidate_ceiling_exhausted",
                                    "root_isolation_inconclusive",
                                    "route_disagreement",
                                },
                                "unknown exact-localizer stop",
                            )
                            require(
                                len(group["channels"]) < len(CHANNELS),
                                "localizer stop cannot follow all eighteen records",
                            )
                        else:
                            require(
                                global_stop is not None,
                                "resource group lacks global stop",
                            )
                            strict_equal(
                                stop, global_stop, "group/global resource stop"
                            )
            else:
                stop = typed_stop(family["stop"], "family stop")
                if status == "premise_stop":
                    strict_equal(
                        stop["owner"], "guarded_shadow", "family premise owner"
                    )
                    require(
                        stop["code"] != "MSEL1ResourceStop",
                        "global resource exception mislabeled as a continuable premise",
                    )
                else:
                    strict_equal(stop["owner"], "resource", "family resource owner")
                    strict_equal(
                        stop["code"], "MSEL1ResourceStop", "family resource cause"
                    )
                    require(
                        global_stop is not None, "resource family lacks global stop"
                    )
                prechecks = family["accepted_source_prechecks"]
                records = family["stage_and_endpoint_records"]
                rhs_calls = family["rhs_calls"]
                require(
                    records % 5 == 0, "partial RK4 proposal cannot retain stage records"
                )
                proposals = records // 5
                require(
                    proposals <= prechecks <= proposals + 1,
                    "prechecks do not describe a one-proposal prefix",
                )
                upper_calls = 6 * proposals + (1 if prechecks == proposals else 6)
                require(
                    prechecks + records <= rhs_calls <= upper_calls,
                    "RHS/precheck/retained-record prefix is impossible",
                )
                if stop["code"] == "source_retry":
                    require(
                        1 <= prechecks <= 7
                        and records == 5 * prechecks
                        and rhs_calls == 6 * prechecks,
                        "source retry lacks its complete proposal evidence",
                    )
                for group in groups:
                    strict_equal(group["status"], status, "stopped family/group status")
                    strict_equal(group["channels"], [], "stopped family channels")
                    strict_equal(group["stop"], stop, "stopped family/group reason")
    counts = exact_keys(accounting, _COUNT_FIELDS | {"static_shells"}, "accounting")
    for key in _COUNT_FIELDS:
        strict_equal(counts[key], totals[key], f"total.{key}")
        require(totals[key] <= 6 * _FAMILY_LIMITS[key], "total shadow budget exceeded")
    shells = integer(counts["static_shells"], "static shell count")
    require(shells in (0, 6), "partial static shell set")
    if attempted_family:
        require(
            shells == 6 and all(width["restored"] is not None for width in widths),
            "shadow records precede complete static/predecessor preparation",
        )

    # A global stop can leave only a suffix unattempted in the fixed schedule.
    halted = False
    for arithmetic in (freeze.ORIGINAL_ARITHMETIC_ID, freeze.COMPENSATED_ARITHMETIC_ID):
        for width in widths:
            family = width["families"][arithmetic]
            if halted:
                require(family is None, "family appears after a global stop")
                continue
            if family is None:
                require(
                    global_stop is not None, "unfinished schedule has no global stop"
                )
                halted = True
                continue
            if family["status"] == "resource_stop":
                # Failed family construction marks all its groups together.
                halted = True
                continue
            if family["status"] == "premise_stop":
                continue
            for group in _family_groups(width, arithmetic):
                if halted:
                    strict_equal(
                        group["status"], "not_attempted", "group after global stop"
                    )
                    continue
                if group["status"] == "not_attempted":
                    require(
                        global_stop is not None, "unattempted group has no global stop"
                    )
                    halted = True
                elif (
                    group["status"] == "resource_stop"
                    and group["stop"]["owner"] == "resource"
                ):
                    halted = True


def validate_raw_terminal(
    value: object,
    *,
    authority_commit: str,
    config_sha256: str,
    freeze_sha256: str,
    environment: dict,
    expected_rows: int = freeze.OWNED_ROW_COUNT,
) -> dict:
    """Authenticate the shape and predicates, without I/O or a shadow replay."""
    require(
        type(authority_commit) is str
        and _COMMIT.fullmatch(authority_commit) is not None,
        "authority must be a full commit identity",
    )
    require(
        type(environment) is dict
        and all(
            type(key) is str and type(item) is str for key, item in environment.items()
        ),
        "environment must be an exact text mapping",
    )
    terminal = exact_keys(
        value,
        {
            "schema",
            "artifact_id",
            "runner_id",
            "authority_commit",
            "config_sha256",
            "freeze_sha256",
            "environment",
            "store_snapshot_before",
            "store_snapshot_after",
            "widths",
            "accounting",
            "global_stop",
            "elapsed_seconds_hex",
            "classification",
            "selected_candidate",
            "candidate_summaries",
            "licenses_only_separate_TDG11_IMP1",
            "independent_PREF1_required_before_IMP1",
            "nonclaims",
        },
        "raw terminal",
    )
    digest(config_sha256, "freeze config hash")
    digest(freeze_sha256, "freeze compact hash")
    fixed = {
        "schema": RAW_SCHEMA,
        "artifact_id": freeze.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "config_sha256": config_sha256,
        "freeze_sha256": freeze_sha256,
        "environment": environment,
        "nonclaims": freeze.nonclaims(),
        "store_snapshot_before": [freeze.STORE_LEAF_COUNT, freeze.STORE_SHA256],
        "store_snapshot_after": [freeze.STORE_LEAF_COUNT, freeze.STORE_SHA256],
    }
    for key, expected in fixed.items():
        strict_equal(terminal[key], expected, f"terminal.{key}")
    reduction = reduce_widths(terminal["widths"], expected_rows=expected_rows)
    for key, expected in reduction.items():
        strict_equal(terminal[key], expected, f"terminal.{key}")
    if terminal["global_stop"] is not None:
        stop = typed_stop(terminal["global_stop"], "global stop")
        strict_equal(stop["owner"], "resource", "global stop owner")
    _validate_families(
        terminal["widths"], terminal["accounting"], global_stop=terminal["global_stop"]
    )
    elapsed = terminal["elapsed_seconds_hex"]
    require(type(elapsed) is str and len(elapsed) <= 64, "elapsed time text differs")
    try:
        clock = float.fromhex(elapsed)
    except (ValueError, OverflowError) as error:
        raise PREF1ProtocolError("elapsed time is not binary64") from error
    require(
        math.isfinite(clock) and clock >= 0 and clock.hex() == elapsed,
        "elapsed time is not canonical/nonnegative",
    )
    return terminal


def validate_raw_manifest(
    value: object,
    *,
    authority_commit: str,
    config_sha256: str,
    freeze_sha256: str,
    terminal_sha256: str,
) -> dict:
    expected = {
        "schema": RAW_SCHEMA,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "config_sha256": config_sha256,
        "freeze_sha256": freeze_sha256,
        "terminal_sha256": terminal_sha256,
        "leaf_names": ["manifest.json", "terminal.json"],
        "endpoint_serialized": False,
        "historical_store_written": False,
    }
    digest(terminal_sha256, "terminal hash")
    strict_equal(value, expected, "raw manifest")
    return value
