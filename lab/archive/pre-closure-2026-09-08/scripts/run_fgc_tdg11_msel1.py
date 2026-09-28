#!/usr/bin/env python3
"""One prospectively authorized TDG11 diagnostic; never a production runner."""

from __future__ import annotations

import argparse
from dataclasses import fields, is_dataclass
from fractions import Fraction as Q
from hashlib import sha256
import math
from pathlib import Path
import re
import signal
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons import evidence_io as io  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_msel1_contract as contract  # noqa: E402
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    CertifiedMagnitudeInterval,
    TDG6ChannelAdmission,
    classify_tdg6_channel,
)
from recursive_horizons.fgc.evolution.tdg11_rational_complete_c import (  # noqa: E402
    RationalCompleteCAdmissionEvidence,
    RationalCompleteCDifferenceEvidence,
)


RUNNER_ID = "FGC-1-TDG11-MSEL1-RUN1"
RAW_SCHEMA = "FGC-1-TDG11-MSEL1-raw-v1"
SELECTED_CLASS = "qualified_candidate_selected_no_state_advance"
NONPASS_CLASS = "all_three_candidates_nonpass_no_state_advance"
INCONCLUSIVE_CLASS = "bounded_selection_inconclusive_no_state_advance"
GROUP_STATUSES = frozenset(
    {"not_attempted", "complete", "premise_stop", "resource_stop"}
)
CHANNELS = contract.TDG6_COMPLETE_STATE_CHANNELS
_HEX = re.compile(r"[0-9a-f]{64}\Z")
_DECIMAL = re.compile(r"-?(?:0|[1-9][0-9]*)\Z")


class MSEL1RunnerError(RuntimeError):
    """Invalid implementation/provenance, not a scientific nonpass."""


class MSEL1ResourceStop(RuntimeError):
    """The prospectively bounded invocation cannot complete further work."""


def _keys(value: object, expected: set[str], label: str) -> dict:
    if type(value) is not dict or set(value) != expected:
        raise MSEL1RunnerError(f"{label} keys differ")
    return value


def _integer(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        raise MSEL1RunnerError(f"{label} must be a nonnegative integer")
    return value


def _digest(value: object, label: str) -> str:
    if type(value) is not str or not _HEX.fullmatch(value):
        raise MSEL1RunnerError(f"{label} is not a SHA-256")
    return value


def _decimal(value: int) -> str:
    """Bounded integer text without changing Python's process-wide digit cap."""
    if abs(value).bit_length() > contract.RESOURCES["maximum_rational_bits"]:
        raise MSEL1ResourceStop("rational_bit_ceiling")
    if value == 0:
        return "0"
    sign, remaining = ("-" if value < 0 else ""), abs(value)
    chunks = []
    while remaining:
        remaining, remainder = divmod(remaining, 10**9)
        chunks.append(remainder)
    return (
        sign
        + str(chunks[-1])
        + "".join(f"{chunk:09d}" for chunk in reversed(chunks[:-1]))
    )


def _parse_decimal(value: object) -> int:
    if (
        type(value) is not str
        or len(value) > 10000
        or not _DECIMAL.fullmatch(value)
        or value == "-0"
    ):
        raise MSEL1RunnerError("invalid exact integer text")
    sign = -1 if value.startswith("-") else 1
    digits = value.removeprefix("-")
    number = 0
    for index in range(0, len(digits), 9):
        chunk = digits[index : index + 9]
        number = number * 10 ** len(chunk) + int(chunk)
    if number.bit_length() > contract.RESOURCES["maximum_rational_bits"]:
        raise MSEL1RunnerError("exact integer exceeds the frozen bit ceiling")
    return sign * number


def exact_object(value: object, depth: int = 0) -> object:
    if depth > 48:
        raise MSEL1RunnerError("exact record nesting exceeds its bound")
    if type(value) is Q:
        return {
            "numerator": _decimal(value.numerator),
            "denominator": _decimal(value.denominator),
        }
    if is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: exact_object(getattr(value, field.name), depth + 1)
            for field in fields(value)
        }
    if type(value) is dict:
        if any(type(key) is not str for key in value):
            raise MSEL1RunnerError("record keys must be strings")
        return {key: exact_object(item, depth + 1) for key, item in value.items()}
    if type(value) in (tuple, list):
        return [exact_object(item, depth + 1) for item in value]
    if type(value) is float:
        if not math.isfinite(value):
            raise MSEL1RunnerError("nonfinite fingerprint scalar")
        return {"binary64_hex": value.hex()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise MSEL1RunnerError(f"unsupported exact record type: {type(value).__name__}")


def exact_fraction(value: object, *, nonnegative: bool = True) -> Q:
    record = _keys(value, {"numerator", "denominator"}, "rational")
    numerator, denominator = (
        _parse_decimal(record["numerator"]),
        _parse_decimal(record["denominator"]),
    )
    if denominator <= 0:
        raise MSEL1RunnerError("rational denominator must be positive")
    answer = Q(numerator, denominator)
    if (answer.numerator, answer.denominator) != (numerator, denominator) or (
        nonnegative and answer < 0
    ):
        raise MSEL1RunnerError("rational is not canonical/nonnegative")
    return answer


def _decode(value: object) -> object:
    if type(value) is dict:
        if set(value) == {"numerator", "denominator"}:
            return exact_fraction(value, nonnegative=False)
        return {key: _decode(item) for key, item in value.items()}
    if type(value) is list:
        return [_decode(item) for item in value]
    return value


def _decision(record: object) -> TDG6ChannelAdmission:
    values = _keys(
        record, {field.name for field in fields(TDG6ChannelAdmission)}, "decision"
    )
    decoded = _decode(values)
    for key in ("outer_difference", "finest_difference"):
        decoded[key] = CertifiedMagnitudeInterval(**decoded[key])
    decision = TDG6ChannelAdmission(**decoded)
    expected = classify_tdg6_channel(
        decision.outer_difference, decision.finest_difference
    )
    contract._strict_equal(exact_object(decision), exact_object(expected), "decision")
    return decision


def validate_exact_evidence(
    record: object, *, expected_rows: int
) -> RationalCompleteCAdmissionEvidence:
    values = _keys(
        record,
        {field.name for field in fields(RationalCompleteCAdmissionEvidence)},
        "exact evidence",
    )
    decoded = _decode(values)
    for key in ("d01", "d12"):
        _keys(
            values[key],
            {field.name for field in fields(RationalCompleteCDifferenceEvidence)},
            key,
        )
        decoded[key] = RationalCompleteCDifferenceEvidence(**decoded[key])
    decoded["decision"] = _decision(values["decision"])
    result = RationalCompleteCAdmissionEvidence(**decoded)
    if (
        result.row_count != expected_rows
        or result.maximum_candidates_D01
        != contract.RESOURCES["per_channel_maximum_candidates_D01"]
        or result.maximum_candidates_D12
        != contract.RESOURCES["per_channel_maximum_candidates_D12"]
        or result.refinement_depth != contract.RESOURCES["primary_refinement_depth"]
    ):
        raise MSEL1RunnerError("exact evidence differs from frozen resources/rows")
    return result


def _roundoff_bounds(value: object) -> tuple[Q, Q, Q]:
    if type(value) is not list or len(value) != 3:
        raise MSEL1RunnerError("roundoff bounds must contain three levels")
    return tuple(exact_fraction(item) for item in value)


def channel_outcome(
    record: object, candidate: str, *, expected_rows: int = contract.OWNED_ROW_COUNT
) -> str:
    """Recompute a saved channel's gate, not merely trust its pass label."""
    if candidate in ("baseline", contract.CANDIDATES[0], contract.CANDIDATES[1]):
        expected_keys = {"channel", "exact_complete_C", "public_fine_debit"}
        if candidate == contract.CANDIDATES[0]:
            expected_keys.add("roundoff_bounds")
        data = _keys(record, expected_keys, "complete-C channel")
        evidence = validate_exact_evidence(
            data["exact_complete_C"], expected_rows=expected_rows
        )
        debit = evidence.d12.upper
        if candidate == contract.CANDIDATES[0]:
            debit += _roundoff_bounds(data["roundoff_bounds"])[2]
        if exact_fraction(data["public_fine_debit"]) != debit:
            raise MSEL1RunnerError("public fine debit omits or changes a contribution")
        if evidence.decision.admission_passed:
            return "pass"
        return "nonpass" if evidence.sufficient_contraction_failure else "inconclusive"
    if candidate != contract.CANDIDATES[2]:
        raise MSEL1RunnerError("unknown candidate")
    data = _keys(
        record,
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
    e0, e1, e2 = _roundoff_bounds(data["estimators"])
    first, second = _decision(data["contraction_01"]), _decision(data["contraction_12"])
    for decision, left, right in ((first, e0, e1), (second, e1, e2)):
        if decision.outer_difference != CertifiedMagnitudeInterval(
            left, left
        ) or decision.finest_difference != CertifiedMagnitudeInterval(right, right):
            raise MSEL1RunnerError(
                "embedded contractions do not bind their estimator levels"
            )
    if data["global_PDE_enclosure"] is not False:
        raise MSEL1RunnerError("embedded defect is not a global PDE enclosure")
    if (
        exact_fraction(data["public_fine_debit"])
        != e2 + _roundoff_bounds(data["roundoff_bounds"])[2]
    ):
        raise MSEL1RunnerError("embedded debit differs")
    if first.admission_passed and second.admission_passed:
        return "pass"
    return "nonpass" if (8 * e1**2 > e0**2 or 8 * e2**2 > e1**2) else "inconclusive"


def empty_group() -> dict[str, object]:
    return {"status": "not_attempted", "channels": [], "stop": None}


def _group_outcomes(group: object, candidate: str, *, expected_rows: int) -> list[str]:
    data = _keys(group, {"status", "channels", "stop"}, "channel group")
    if (
        type(data["status"]) is not str
        or data["status"] not in GROUP_STATUSES
        or type(data["channels"]) is not list
    ):
        raise MSEL1RunnerError("channel group status/list differs")
    values = data["channels"]
    if len(values) > len(CHANNELS) or (
        data["status"] == "complete" and len(values) != len(CHANNELS)
    ):
        raise MSEL1RunnerError("channel group is incomplete or oversized")
    if data["status"] == "not_attempted" and values:
        raise MSEL1RunnerError("unattempted group contains measurements")
    if data["status"] in ("complete", "not_attempted"):
        if data["stop"] is not None:
            raise MSEL1RunnerError("completed/unattempted group carries a stop")
    else:
        stop = _keys(data["stop"], {"owner", "code", "detail"}, "typed stop")
        if any(
            type(value) is not str or not value or len(value) > 640
            for value in stop.values()
        ):
            raise MSEL1RunnerError("invalid typed stop")
    outcomes = []
    for index, record in enumerate(values):
        if type(record) is not dict or record.get("channel") != CHANNELS[index]:
            raise MSEL1RunnerError("channel order differs")
        outcomes.append(channel_outcome(record, candidate, expected_rows=expected_rows))
    return outcomes


def reduce_selection(
    width_records: object, *, expected_rows: int = contract.OWNED_ROW_COUNT
) -> dict[str, object]:
    if type(width_records) is not list or len(width_records) != 3:
        raise MSEL1RunnerError("selection needs exactly three width records")
    outcomes = {candidate: [] for candidate in contract.CANDIDATES}
    complete = {candidate: True for candidate in contract.CANDIDATES}
    for width, spec in zip(width_records, contract.replay_specs(), strict=True):
        data = _keys(
            width,
            {
                "retry",
                "generation",
                "width_hex",
                "restored",
                "families",
                "baseline",
                "candidates",
            },
            "width",
        )
        for key in ("retry", "generation", "width_hex"):
            contract._strict_equal(data[key], spec[key], f"width.{key}")
        _keys(data["candidates"], set(contract.CANDIDATES), "candidate groups")
        _group_outcomes(data["baseline"], "baseline", expected_rows=expected_rows)
        for candidate in contract.CANDIDATES:
            group = data["candidates"][candidate]
            outcomes[candidate].extend(
                _group_outcomes(group, candidate, expected_rows=expected_rows)
            )
            complete[candidate] &= group["status"] == "complete"
    summaries = []
    for candidate in contract.CANDIDATES:
        values = outcomes[candidate]
        if (
            complete[candidate]
            and len(values) == 54
            and all(value == "pass" for value in values)
        ):
            status = "pass"
        elif "nonpass" in values:
            status = "nonpass"
        else:
            status = "inconclusive"
        summaries.append(
            {
                "candidate": candidate,
                "status": status,
                "classified_channels": len(values),
                "complete_all_widths": complete[candidate],
            }
        )
    selected = next(
        (item["candidate"] for item in summaries if item["status"] == "pass"), None
    )
    classification = (
        SELECTED_CLASS
        if selected is not None
        else NONPASS_CLASS
        if all(item["status"] == "nonpass" for item in summaries)
        else INCONCLUSIVE_CLASS
    )
    return {
        "classification": classification,
        "selected_candidate": selected,
        "candidate_summaries": summaries,
        "licenses_only_separate_TDG11_IMP1": False,
        "independent_PREF1_required_before_IMP1": True,
    }


def _new_width(spec: dict) -> dict:
    return {
        "retry": spec["retry"],
        "generation": spec["generation"],
        "width_hex": spec["width_hex"],
        "restored": None,
        "families": {
            contract.ORIGINAL_ARITHMETIC_ID: None,
            contract.COMPENSATED_ARITHMETIC_ID: None,
        },
        "baseline": empty_group(),
        "candidates": {candidate: empty_group() for candidate in contract.CANDIDATES},
    }


def _stop(owner: str, code: str, error: object) -> dict[str, str]:
    return {
        "owner": owner,
        "code": code,
        "detail": (" ".join(str(error).replace(str(ROOT), "<repo>").split()) or code)[
            :640
        ],
    }


def _exact_channel(
    channel: str, evidence, *, roundoff: tuple[Q, Q, Q] | None = None
) -> dict:
    debit = evidence.d12.upper
    answer = {"channel": channel, "exact_complete_C": exact_object(evidence)}
    if roundoff is not None:
        answer["roundoff_bounds"] = exact_object(roundoff)
        debit += roundoff[2]
    answer["public_fine_debit"] = exact_object(debit)
    return answer


def _embedded_channel(channel: str, data) -> dict:
    e0, e1, e2 = data.embedded_defect_bounds
    return {
        "channel": channel,
        "estimators": exact_object((e0, e1, e2)),
        "roundoff_bounds": exact_object(data.accumulation_bounds),
        "contraction_01": exact_object(
            classify_tdg6_channel(
                CertifiedMagnitudeInterval(e0, e0), CertifiedMagnitudeInterval(e1, e1)
            )
        ),
        "contraction_12": exact_object(
            classify_tdg6_channel(
                CertifiedMagnitudeInterval(e1, e1), CertifiedMagnitudeInterval(e2, e2)
            )
        ),
        "public_fine_debit": exact_object(e2 + data.accumulation_bounds[2]),
        "global_PDE_enclosure": False,
    }


def _assess_group(family, group: dict, candidate: str, progress) -> None:
    from recursive_horizons.fgc.evolution.tdg11_msel1_reconstruction import (
        reconstruct_channel,
    )
    from recursive_horizons.fgc.evolution.tdg11_rational_complete_c import (
        RationalCompleteCResourceExhausted,
        RationalCompleteCRouteDisagreement,
        assess_rational_complete_c_rows,
    )

    try:
        for channel in CHANNELS:
            data = reconstruct_channel(family.recorded, channel)
            if candidate == contract.CANDIDATES[2]:
                record = _embedded_channel(channel, data)
            else:
                evidence = assess_rational_complete_c_rows(
                    data.corrected_rows
                    if candidate == contract.CANDIDATES[0]
                    else data.raw_rows,
                    expected_row_count=contract.OWNED_ROW_COUNT,
                    maximum_candidates_D01=contract.RESOURCES[
                        "per_channel_maximum_candidates_D01"
                    ],
                    maximum_candidates_D12=contract.RESOURCES[
                        "per_channel_maximum_candidates_D12"
                    ],
                    refinement_depth=contract.RESOURCES["primary_refinement_depth"],
                )
                record = _exact_channel(
                    channel,
                    evidence,
                    roundoff=data.accumulation_bounds
                    if candidate == contract.CANDIDATES[0]
                    else None,
                )
            group["channels"].append(record)
            progress({"candidate": candidate, "channel": channel})
        group["status"] = "complete"
    except RationalCompleteCResourceExhausted as error:
        group.update(
            status="resource_stop",
            stop=_stop("exact_localizer", error.evidence.reason, error),
        )
    except RationalCompleteCRouteDisagreement as error:
        group.update(
            status="resource_stop",
            stop=_stop("exact_localizer", "route_disagreement", error),
        )


def _execution_authority(root: Path, commit: str):
    return _authority_module().authorize_execution(root, authority_commit=commit)


def _authority_module():
    from recursive_horizons.fgc.evolution import tdg11_msel1_authority

    return tdg11_msel1_authority


def validate_terminal(
    value: object,
    *,
    authority_commit: str,
    config_sha256: str,
    freeze_sha256: str,
    environment: dict,
    expected_rows: int = contract.OWNED_ROW_COUNT,
) -> dict:
    expected_keys = {
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
    }
    terminal = _keys(value, expected_keys, "terminal")
    for key, expected in (
        ("schema", RAW_SCHEMA),
        ("artifact_id", contract.ARTIFACT_ID),
        ("runner_id", RUNNER_ID),
        ("authority_commit", authority_commit),
        ("config_sha256", config_sha256),
        ("freeze_sha256", freeze_sha256),
        ("environment", environment),
        ("nonclaims", contract.nonclaims()),
    ):
        contract._strict_equal(terminal[key], expected, f"terminal.{key}")
    _digest(config_sha256, "config hash")
    _digest(freeze_sha256, "freeze hash")
    for key in ("store_snapshot_before", "store_snapshot_after"):
        contract._strict_equal(
            terminal[key], [contract.STORE_LEAF_COUNT, contract.STORE_SHA256], key
        )
    if type(terminal["widths"]) is not list or len(terminal["widths"]) != 3:
        raise MSEL1RunnerError("terminal does not contain all width slots")
    counted = {
        "accepted_source_prechecks": 0,
        "stage_and_endpoint_records": 0,
        "rhs_calls": 0,
    }
    for width, spec in zip(terminal["widths"], contract.replay_specs(), strict=True):
        _keys(
            width,
            {
                "retry",
                "generation",
                "width_hex",
                "restored",
                "families",
                "baseline",
                "candidates",
            },
            "width",
        )
        families = _keys(
            width["families"],
            {contract.ORIGINAL_ARITHMETIC_ID, contract.COMPENSATED_ARITHMETIC_ID},
            "families",
        )
        _keys(width["candidates"], set(contract.CANDIDATES), "candidate groups")
        for group in [width["baseline"], *width["candidates"].values()]:
            _keys(group, {"status", "channels", "stop"}, "group")
        if width["restored"] is not None:
            restored = _keys(
                width["restored"],
                {
                    "checkpoint_sha256",
                    "descriptor_sha256",
                    "state_sha256",
                    "fingerprint_sha256",
                },
                "restored",
            )
            for key, expected in (
                ("checkpoint_sha256", spec["checkpoint_sha256"]),
                ("descriptor_sha256", contract.DESCRIPTOR_SHA256),
                ("state_sha256", contract.PHYSICAL_STATE_SHA256),
            ):
                contract._strict_equal(restored[key], expected, f"restored.{key}")
            _digest(restored["fingerprint_sha256"], "restored fingerprint")
        for arithmetic, family in families.items():
            groups = (
                [
                    width["baseline"],
                    width["candidates"][contract.CANDIDATES[0]],
                    width["candidates"][contract.CANDIDATES[2]],
                ]
                if arithmetic == contract.ORIGINAL_ARITHMETIC_ID
                else [width["candidates"][contract.CANDIDATES[1]]]
            )
            if family is None:
                if any(group["status"] != "not_attempted" for group in groups):
                    raise MSEL1RunnerError("measurements lack a guarded family")
                continue
            if width["restored"] is None:
                raise MSEL1RunnerError("family lacks an authenticated predecessor")
            if type(family) is not dict:
                raise MSEL1RunnerError("family record is malformed")
            complete = family.get("status") == "complete"
            required = {
                "status",
                "accepted_source_prechecks",
                "stage_and_endpoint_records",
                "rhs_calls",
                "family_sha256" if complete else "stop",
            }
            _keys(family, required, "family record")
            if family["status"] not in {"complete", "premise_stop", "resource_stop"}:
                raise MSEL1RunnerError("family status differs")
            for key in counted:
                counted[key] += _integer(family[key], f"family.{key}")
            if complete:
                _digest(family["family_sha256"], "family hash")
                if (
                    family["accepted_source_prechecks"],
                    family["stage_and_endpoint_records"],
                    family["rhs_calls"],
                ) != (7, 35, 42):
                    raise MSEL1RunnerError("complete RK4 family counts differ")
            else:
                _keys(family["stop"], {"owner", "code", "detail"}, "family stop")
                if any(
                    group["status"] != family["status"] or group["channels"]
                    for group in groups
                ):
                    raise MSEL1RunnerError(
                        "failed family contains completed measurements"
                    )
            if (
                family["accepted_source_prechecks"] > 7
                or family["stage_and_endpoint_records"] > 35
                or family["rhs_calls"] > 42
            ):
                raise MSEL1RunnerError("family exceeded its frozen work bound")
    accounting = _keys(
        terminal["accounting"], set(counted) | {"static_shells"}, "accounting"
    )
    for key, expected in counted.items():
        contract._strict_equal(accounting[key], expected, f"accounting.{key}")
    if type(accounting["static_shells"]) is not int or accounting[
        "static_shells"
    ] not in (0, 6):
        raise MSEL1RunnerError("static shell return count differs")
    if counted["rhs_calls"] and accounting["static_shells"] != 6:
        raise MSEL1RunnerError("shadows lack the complete static factory")
    if terminal["global_stop"] is not None:
        stop = _keys(
            terminal["global_stop"], {"owner", "code", "detail"}, "global stop"
        )
        if stop["owner"] != "resource":
            raise MSEL1RunnerError("global stop is not a typed resource stop")
    elapsed = terminal["elapsed_seconds_hex"]
    if type(elapsed) is not str:
        raise MSEL1RunnerError("elapsed time must be binary64 text")
    clock = float.fromhex(elapsed)
    if not math.isfinite(clock) or clock < 0 or clock.hex() != elapsed:
        raise MSEL1RunnerError("elapsed time is invalid")
    reduction = reduce_selection(terminal["widths"], expected_rows=expected_rows)
    for key, expected in reduction.items():
        contract._strict_equal(terminal[key], expected, f"terminal.{key}")
    return terminal


def status(root: Path, *, authority_commit: str) -> dict:
    try:
        receipt = _execution_authority(root, authority_commit)
    except Exception as error:
        return {
            "safe_to_run": False,
            "authority_commit": authority_commit,
            "reason": _stop("preflight", type(error).__name__, error),
        }
    return {
        "safe_to_run": True,
        "authority_commit": receipt.authority_commit,
        "output_namespace": contract.OUTPUT_NAMESPACE,
        "scientific_run_started": False,
    }


def run(root: Path, *, authority_commit: str, progress=None) -> dict:
    """Authorize first, construct each frozen family once, publish only two leaves."""
    authority = _authority_module()
    receipt = _execution_authority(root, authority_commit)
    # Heavy/runtime imports occur only after the authority gate.
    from recursive_horizons.fgc.evolution.numerical_engine import (
        PRIMARY_METHOD,
        array_content_sha256,
    )
    from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (
        build_static_gr0_shells,
    )
    from recursive_horizons.fgc.evolution.tdg11_msel1_runtime import (
        build_guarded_shadow_family,
        TDG11ShadowPremiseStop,
    )

    started = time.monotonic()
    callback = (lambda _event: None) if progress is None else progress
    widths = [_new_width(spec) for spec in contract.replay_specs()]
    totals = {
        "accepted_source_prechecks": 0,
        "stage_and_endpoint_records": 0,
        "rhs_calls": 0,
        "static_shells": 0,
    }
    global_stop = None
    active_group = None

    def report(event):
        if time.monotonic() - started > contract.RESOURCES["maximum_wall_seconds"]:
            raise MSEL1ResourceStop("wall_time_ceiling")
        callback(event)

    def alarm(_signum, _frame):
        raise MSEL1ResourceStop("wall_time_ceiling")

    if signal.getitimer(signal.ITIMER_REAL) != (0.0, 0.0):
        raise MSEL1RunnerError(
            "an existing process alarm prevents exclusive budget ownership"
        )
    previous_handler = signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, contract.RESOURCES["maximum_wall_seconds"])
    try:
        templates = build_static_gr0_shells(root)
        totals["static_shells"] = len(templates)
        if totals["static_shells"] != 6:
            raise MSEL1RunnerError("static shell count differs")
        restored = {}
        for width in widths:
            item = authority.restore_predecessor(
                root, width["retry"], templates=templates
            )
            restored[width["retry"]] = item
            width["restored"] = {
                "checkpoint_sha256": item.checkpoint_sha256,
                "descriptor_sha256": item.descriptor_sha256,
                "state_sha256": array_content_sha256(
                    item.member.state.u, item.member.state.p, item.member.state.q
                ),
                "fingerprint_sha256": sha256(
                    io.canonical_json_bytes(exact_object(item.fingerprint))
                ).hexdigest(),
            }
        for arithmetic_id in (
            contract.ORIGINAL_ARITHMETIC_ID,
            contract.COMPENSATED_ARITHMETIC_ID,
        ):
            for width in widths:
                report(
                    {
                        "retry": width["retry"],
                        "arithmetic": arithmetic_id,
                        "phase": "guarded_family",
                    }
                )
                member = restored[width["retry"]].member
                candidates = (
                    ("baseline", contract.CANDIDATES[0], contract.CANDIDATES[2])
                    if arithmetic_id == contract.ORIGINAL_ARITHMETIC_ID
                    else (contract.CANDIDATES[1],)
                )
                groups = [
                    width["baseline"]
                    if candidate == "baseline"
                    else width["candidates"][candidate]
                    for candidate in candidates
                ]
                try:
                    family = build_guarded_shadow_family(
                        method=PRIMARY_METHOD,
                        arithmetic_id=arithmetic_id,
                        time=member.time,
                        step_size=float.fromhex(width["width_hex"]),
                        state=member.state,
                        rhs=member.operator,
                        projector=member.projector,
                        transaction=member.transaction,
                        tracers=member.tracers,
                        coordinates=member.initial.grid.coordinates,
                        previous_step_index=member.step_index,
                        previous_transaction_serial=member.transaction_serial,
                    )
                except TDG11ShadowPremiseStop as error:
                    totals["accepted_source_prechecks"] += (
                        error.accepted_source_precheck_count
                    )
                    totals["stage_and_endpoint_records"] += error.stage_record_count
                    totals["rhs_calls"] += error.rhs_call_count
                    resource_failure = isinstance(error.cause, MSEL1ResourceStop)
                    stop = _stop(
                        "resource" if resource_failure else "guarded_shadow",
                        "source_retry"
                        if error.source_retry is not None
                        else type(error.cause).__name__,
                        error,
                    )
                    width["families"][arithmetic_id] = {
                        "status": "resource_stop"
                        if resource_failure
                        else "premise_stop",
                        "stop": stop,
                        "accepted_source_prechecks": error.accepted_source_precheck_count,
                        "stage_and_endpoint_records": error.stage_record_count,
                        "rhs_calls": error.rhs_call_count,
                    }
                    for group in groups:
                        group.update(
                            status="resource_stop"
                            if resource_failure
                            else "premise_stop",
                            stop=stop,
                        )
                    if resource_failure:
                        raise error.cause
                    continue
                totals["accepted_source_prechecks"] += (
                    family.accepted_source_precheck_count
                )
                totals["stage_and_endpoint_records"] += family.stage_record_count
                totals["rhs_calls"] += family.rhs_call_count
                width["families"][arithmetic_id] = {
                    "status": "complete",
                    "family_sha256": family.recorded.family_sha256,
                    "accepted_source_prechecks": family.accepted_source_precheck_count,
                    "stage_and_endpoint_records": family.stage_record_count,
                    "rhs_calls": family.rhs_call_count,
                }
                for candidate, group in zip(candidates, groups, strict=True):
                    active_group = group
                    _assess_group(family, group, candidate, report)
                    active_group = None
    except MSEL1ResourceStop as error:
        global_stop = _stop("resource", str(error), error)
        if active_group is not None and active_group["status"] != "complete":
            active_group.update(status="resource_stop", stop=global_stop)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, previous_handler)

    # Reauthenticate before publication; source/implementation drift is invalid,
    # not a publishable scientific nonpass or a permission to retry.
    final_authority = _execution_authority(root, authority_commit)
    if final_authority.store_snapshot != receipt.store_snapshot:
        raise MSEL1RunnerError("historical store drifted")
    if (
        totals["accepted_source_prechecks"] > 42
        or totals["stage_and_endpoint_records"] > 210
        or totals["rhs_calls"] > 252
    ):
        raise MSEL1RunnerError("frozen shadow budget was exceeded")
    reduction = reduce_selection(widths, expected_rows=contract.OWNED_ROW_COUNT)
    terminal = {
        "schema": RAW_SCHEMA,
        "artifact_id": contract.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "config_sha256": receipt.config_sha256,
        "freeze_sha256": receipt.result_sha256,
        "environment": receipt.environment,
        "store_snapshot_before": list(receipt.store_snapshot),
        "store_snapshot_after": list(final_authority.store_snapshot),
        "widths": widths,
        "accounting": totals,
        "global_stop": global_stop,
        "elapsed_seconds_hex": float(time.monotonic() - started).hex(),
        **reduction,
        "nonclaims": contract.nonclaims(),
    }
    validate_terminal(
        terminal,
        authority_commit=authority_commit,
        config_sha256=receipt.config_sha256,
        freeze_sha256=receipt.result_sha256,
        environment=receipt.environment,
        expected_rows=contract.OWNED_ROW_COUNT,
    )
    raw = io.canonical_json_bytes(terminal)
    manifest = {
        "schema": RAW_SCHEMA,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "config_sha256": receipt.config_sha256,
        "freeze_sha256": receipt.result_sha256,
        "terminal_sha256": sha256(raw).hexdigest(),
        "leaf_names": ["manifest.json", "terminal.json"],
        "endpoint_serialized": False,
        "historical_store_written": False,
    }
    io.publish_exclusive_directory(
        root,
        contract.OUTPUT_NAMESPACE,
        {"manifest.json": io.canonical_json_bytes(manifest), "terminal.json": raw},
    )
    return terminal


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--status", action="store_true")
    action.add_argument("--run", action="store_true")
    parser.add_argument("--authority-commit", required=True)
    args = parser.parse_args(argv)
    if args.status:
        answer = status(ROOT, authority_commit=args.authority_commit)
    else:
        answer = run(
            ROOT,
            authority_commit=args.authority_commit,
            progress=lambda event: print(
                io.canonical_json_bytes(event).decode(), file=sys.stderr, flush=True
            ),
        )
        answer = {
            key: answer[key]
            for key in (
                "classification",
                "selected_candidate",
                "candidate_summaries",
                "accounting",
            )
        }
    print(io.canonical_json_bytes(answer).decode())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
