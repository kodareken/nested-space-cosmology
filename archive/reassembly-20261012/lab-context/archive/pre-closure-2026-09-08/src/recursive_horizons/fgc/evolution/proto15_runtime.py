"""Synthetic-only durable cursor runtime frozen by FGC-2-SF1-PROTO15.

This is deliberately a persistence/state-machine instrument.  It owns neither
the PDE nor a production runner: TDG6 supplies temporal evidence and typed
outcomes, TDG7 supplies safe plans.  The on-disk format is canonical JSON and
is intended for synthetic fault qualification only.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from typing import Any, Callable, Final, Mapping, NoReturn

import numpy as np

from .tdg6_temporal_admission_design import (
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
    TDG6_MINIMUM_MACRO_STEP,
)
from .tdg6_temporal_admission_runtime import (
    TDG6TemporalRetryEvidence,
    TDG6TemporalLedger,
    TDG6TemporalRetryExhausted,
    TDG6TemporalRetryRequired,
    _magnitude_interval,
    classify_tdg6_runtime_intervals,
)
from .tdg6_temporal_admission_design import TDG6_COMPLETE_STATE_CHANNELS
from .tdg7_binary64_subdivision_lattice import (
    TDG7Binary64SubdivisionPlan,
    plan_forward_proto14_subdivision,
    validate_tdg7_binary64_subdivision_plan,
)
from .numerical_engine import COMPARATOR_METHOD, PRIMARY_METHOD


ROOT_SHA256: Final[str] = "0" * 64
MEMBER_KEYS: Final[tuple[str, ...]] = (
    "RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
)
CURSOR_MODES: Final[frozenset[str]] = frozenset({"FRESH_READY", "RETRY_PENDING"})
RECORD_KINDS: Final[frozenset[str]] = frozenset({
    "TDG6_REJECTION", "CURSOR_TRANSITION", "ACCEPTED_FINE_COMMIT",
    "TERMINAL_INVALID", "COMMON_EVENT_ROLLBACK",
})
FaultHook = Callable[[str, str], None]
FAULT_PHASES: Final[tuple[str, ...]] = (
    "before_write", "after_write", "before_flush", "after_flush",
    "before_fsync", "after_fsync", "before_replace", "after_replace",
    "before_dir_fsync", "after_dir_fsync",
)
PLAN_FIELDS: Final[tuple[str, ...]] = (
    "current_hex", "event_target_hex", "requested_cap_hex", "quantum_hex",
    "selection_budget_hex", "macro_width_hex", "fine_width_hex",
    "conservative_reduction_hex", "minimum_width_hex_or_none",
    "boundaries_hex", "target_limited",
)
RUNTIME_METHOD_BY_LABEL: Final[Mapping[str, str]] = {
    "RK4": PRIMARY_METHOD,
    "SSPRK3": COMPARATOR_METHOD,
}
ROLLBACK_PROTOCOL_JOURNAL_FIELD_LOCATIONS: Final[Mapping[str, str]] = {
    "new_journal_parent_sha256": "journal_envelope.journal_parent_sha256",
    "new_journal_record_sha256": "journal_envelope.record_sha256",
}
_AMBIGUOUS_REPLACE_CRASH_PHASES: Final[frozenset[str]] = frozenset({
    "after_replace", "before_dir_fsync",
})


class Proto15DurabilityError(RuntimeError):
    """A write boundary did not become durable; no adoption may follow."""


class Proto15CursorContractError(ValueError):
    """A cursor or transition differs from the frozen lifecycle contract."""


class Proto15JournalIntegrityError(ValueError):
    """A journal frame is malformed, noncanonical, or not contiguous."""


class Proto15CheckpointIntegrityError(ValueError):
    """A campaign checkpoint is incomplete, forked, or hash-invalid."""


class Proto15RecoveryInvalidStop(RuntimeError):
    """Recovery found a tail other than PROTO15's exact permitted branches."""

    physical_classification = False


class Proto15TerminalLocked(RuntimeError):
    """The terminal lock forbids both entrypoints for this campaign identity."""

    physical_classification = False


def _canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False).encode("utf-8")


def _hash(value: object) -> str:
    return sha256(_canonical(value)).hexdigest()


def _require_sha(name: str, value: object) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise Proto15CursorContractError(f"{name} must be a SHA-256 hex digest")
    try:
        int(value, 16)
    except ValueError as error:
        raise Proto15CursorContractError(f"{name} must be hexadecimal") from error
    if value.lower() != value:
        raise Proto15CursorContractError(f"{name} must use lowercase hexadecimal")
    return value


def _nonnegative_integer(name: str, value: object, *, error: type[Exception] = Proto15CursorContractError) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise error(f"{name} must be a nonnegative integer")
    return value


def _finite_hex(name: str, value: object, *, positive: bool = False,
                error: type[Exception] = Proto15CursorContractError) -> float:
    if not isinstance(value, str):
        raise error(f"{name} must be canonical binary64 hexadecimal text")
    try:
        answer = float.fromhex(value)
    except ValueError as cause:
        raise error(f"{name} is malformed") from cause
    if not math.isfinite(answer) or answer.hex() != value or (positive and answer <= 0.0):
        raise error(f"{name} is not canonical finite binary64 hexadecimal text")
    return answer


def _same_binary64(left: object, right: object) -> bool:
    try:
        return np.float64(float(left)).tobytes() == np.float64(float(right)).tobytes()
    except (TypeError, ValueError, OverflowError):
        return False


def _decode_canonical_json(encoded: bytes, *, context: str,
                           error: type[Exception]) -> object:
    def reject_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        answer: dict[str, object] = {}
        for key, value in pairs:
            if key in answer:
                raise ValueError(f"duplicate key {key}")
            answer[key] = value
        return answer

    try:
        value = json.loads(encoded.decode("utf-8"), object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as cause:
        raise error(f"{context} is malformed") from cause
    if _canonical(value) != encoded:
        raise error(f"{context} is noncanonical")
    return value


def _json_safe(value: object) -> object:
    """Return the canonical JSON value, normalizing tuples to arrays."""
    return _decode_canonical_json(
        _canonical(value), context="canonical JSON value", error=Proto15CursorContractError
    )


def _plan_mapping(value: TDG7Binary64SubdivisionPlan | Mapping[str, object]) -> dict[str, object]:
    if isinstance(value, TDG7Binary64SubdivisionPlan):
        plan = value
    elif isinstance(value, Mapping) and set(value) == set(PLAN_FIELDS):
        boundaries = value["boundaries_hex"]
        if not isinstance(boundaries, list) or len(boundaries) != 5:
            raise Proto15CursorContractError("TDG7 plan boundaries differ")
        minimum = value["minimum_width_hex_or_none"]
        try:
            plan = TDG7Binary64SubdivisionPlan(
                current=_finite_hex("plan current", value["current_hex"]),
                event_target=_finite_hex("plan event target", value["event_target_hex"]),
                requested_cap=_finite_hex("plan requested cap", value["requested_cap_hex"], positive=True),
                quantum=_finite_hex("plan quantum", value["quantum_hex"], positive=True),
                selection_budget=_finite_hex("plan selection budget", value["selection_budget_hex"], positive=True),
                macro_width=_finite_hex("plan macro width", value["macro_width_hex"], positive=True),
                fine_width=_finite_hex("plan fine width", value["fine_width_hex"], positive=True),
                conservative_reduction=_finite_hex("plan conservative reduction", value["conservative_reduction_hex"]),
                minimum_width=None if minimum is None else _finite_hex("plan minimum width", minimum, positive=True),
                boundaries=tuple(_finite_hex("plan boundary", item) for item in boundaries),
                target_limited=value["target_limited"],
            )
        except (TypeError, ValueError) as cause:
            raise Proto15CursorContractError("TDG7 plan mapping differs") from cause
    else:
        raise Proto15CursorContractError("retry plan must be a TDG7 plan or exact persisted plan mapping")
    validate_tdg7_binary64_subdivision_plan(plan)
    return {
        "current_hex": plan.current.hex(),
        "event_target_hex": plan.event_target.hex(),
        "requested_cap_hex": plan.requested_cap.hex(),
        "quantum_hex": plan.quantum.hex(),
        "selection_budget_hex": plan.selection_budget.hex(),
        "macro_width_hex": plan.macro_width.hex(),
        "fine_width_hex": plan.fine_width.hex(),
        "conservative_reduction_hex": plan.conservative_reduction.hex(),
        "minimum_width_hex_or_none": None if plan.minimum_width is None else plan.minimum_width.hex(),
        "boundaries_hex": [item.hex() for item in plan.boundaries],
        "target_limited": plan.target_limited,
    }


def _validated_accepted_plan(
    cursor: "Proto15Cursor",
    executed_plan: TDG7Binary64SubdivisionPlan | Mapping[str, object],
    accepted_boundary_time: Mapping[str, object],
) -> dict[str, object]:
    """Bind an accepted commit to the exact TDG7 plan that reached it."""
    plan = _plan_mapping(executed_plan)
    accepted = float.fromhex(str(accepted_boundary_time["binary64_hex"]))
    current = float.fromhex(str(cursor.payload["accepted_boundary_time"]["binary64_hex"]))
    target = float.fromhex(str(cursor.payload["event_target_time"]["binary64_hex"]))
    endpoint = float.fromhex(str(plan["boundaries_hex"][-1]))
    if (not _same_binary64(float.fromhex(str(plan["current_hex"])), current)
            or not _same_binary64(float.fromhex(str(plan["event_target_hex"])), target)
            or not _same_binary64(endpoint, accepted)):
        raise Proto15CursorContractError(
            "accepted fine commit does not equal its exact executed TDG7 plan"
        )
    if cursor.mode == "RETRY_PENDING":
        pending = cursor.payload["retry_successor_payload_or_none"]
        if not isinstance(pending, Mapping) or plan != pending["successor_plan"]:
            raise Proto15CursorContractError(
                "accepted pending retry did not execute its exact durable TDG7 successor plan"
            )
    return plan


def _time_identity(value: object) -> dict[str, str]:
    """Return PROTO15's rational plus exact binary64 identity."""
    if isinstance(value, Mapping):
        if set(value) != {"rational", "binary64_hex"}:
            raise Proto15CursorContractError("time identity keys differ")
        rational, binary = value["rational"], value["binary64_hex"]
        if not isinstance(rational, str) or not isinstance(binary, str):
            raise Proto15CursorContractError("time identity values must be text")
        try:
            fraction = Fraction(rational)
        except (ValueError, ZeroDivisionError) as error:
            raise Proto15CursorContractError("time identity is malformed") from error
        observed = _finite_hex("time binary64 identity", binary)
        if float(fraction).hex() != observed.hex() or str(fraction) != rational:
            raise Proto15CursorContractError("rational and binary64 time identities differ")
        return {"rational": str(fraction), "binary64_hex": binary}
    if isinstance(value, bool) or not isinstance(value, (int, float, Fraction)):
        raise Proto15CursorContractError("time must be a finite numeric identity")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise Proto15CursorContractError("time must be finite")
    rational = value if isinstance(value, Fraction) else Fraction.from_float(numeric)
    if float(rational).hex() != numeric.hex():
        raise Proto15CursorContractError("time rational is not represented by its binary64 identity")
    return {"rational": str(rational), "binary64_hex": numeric.hex()}


def _ledger_mapping(ledger: TDG6TemporalLedger) -> dict[str, object]:
    if not isinstance(ledger, TDG6TemporalLedger):
        raise TypeError("ledger must be TDG6TemporalLedger")
    return {
        "last_accepted_time_hex": ledger.last_accepted_time.hex(),
        "accepted_macro_step_count": ledger.accepted_macro_step_count,
        "cumulative_temporal_retry_count": ledger.cumulative_temporal_retry_count,
        "current_macro_step_temporal_retry_count": ledger.current_macro_step_temporal_retry_count,
        "last_accepted_macro_step_temporal_retry_count": ledger.last_accepted_macro_step_temporal_retry_count,
        "accumulated_debit_vector_hex": [item.hex() for item in ledger.accumulated_debit_vector],
        "serialized_temporal_rejections": list(ledger.serialized_temporal_rejections),
    }


def _ledger_from_mapping(value: Mapping[str, object]) -> TDG6TemporalLedger:
    required = {"last_accepted_time_hex", "accepted_macro_step_count", "cumulative_temporal_retry_count",
                "current_macro_step_temporal_retry_count", "last_accepted_macro_step_temporal_retry_count",
                "accumulated_debit_vector_hex", "serialized_temporal_rejections"}
    if set(value) != required:
        raise Proto15CheckpointIntegrityError("complete TDG6 ledger mapping differs")
    debit = value["accumulated_debit_vector_hex"]
    serialized = value["serialized_temporal_rejections"]
    if not isinstance(debit, list) or len(debit) != len(TDG6_COMPLETE_STATE_CHANNELS):
        raise Proto15CheckpointIntegrityError("TDG6 ledger debit vector differs")
    if not isinstance(serialized, list) or any(not isinstance(item, str) for item in serialized):
        raise Proto15CheckpointIntegrityError("TDG6 rejection ledger differs")
    try:
        return TDG6TemporalLedger(
            last_accepted_time=_finite_hex("ledger accepted time", value["last_accepted_time_hex"], error=Proto15CheckpointIntegrityError),
            accepted_macro_step_count=_nonnegative_integer("accepted step count", value["accepted_macro_step_count"], error=Proto15CheckpointIntegrityError),
            cumulative_temporal_retry_count=_nonnegative_integer("cumulative retry count", value["cumulative_temporal_retry_count"], error=Proto15CheckpointIntegrityError),
            current_macro_step_temporal_retry_count=_nonnegative_integer("current retry count", value["current_macro_step_temporal_retry_count"], error=Proto15CheckpointIntegrityError),
            last_accepted_macro_step_temporal_retry_count=_nonnegative_integer("last accepted retry count", value["last_accepted_macro_step_temporal_retry_count"], error=Proto15CheckpointIntegrityError),
            accumulated_debit_vector=tuple(_finite_hex("accumulated debit", item, error=Proto15CheckpointIntegrityError) for item in debit),
            serialized_temporal_rejections=tuple(serialized),
        )
    except (TypeError, ValueError) as error:
        if isinstance(error, Proto15CheckpointIntegrityError):
            raise
        raise Proto15CheckpointIntegrityError("TDG6 ledger mapping is invalid") from error


_RETRY_EVIDENCE_FIELDS: Final[frozenset[str]] = frozenset({
    "event_type", "method", "initial_time", "attempted_macro_step_size",
    "retry_step_size", "initial_state_sha256", "previous_step_index",
    "previous_transaction_serial", "retry_count_for_current_macro_step",
    "cumulative_temporal_retry_count", "failed_channels", "channel_admissions",
    "accepted_state_bitwise_preserved", "real_monitor_and_causal_ledgers_preserved",
    "real_tracer_ledger_preserved", "accumulated_debit_unchanged",
    "classification_is_numerical_not_physical",
})


def _validate_retry_evidence_mapping(
    value: object, *, method: str, state_sha256: str,
    accepted_time: Mapping[str, object], step_index: int,
    transaction_serial: int, prior_ledger: TDG6TemporalLedger,
) -> dict[str, object]:
    if not isinstance(value, Mapping) or set(value) != _RETRY_EVIDENCE_FIELDS:
        raise Proto15CursorContractError("TDG6 retry evidence fields differ")
    evidence = dict(value)
    if (evidence["event_type"] != "rejected_TDG6_temporal_admission"
            or evidence["method"] != RUNTIME_METHOD_BY_LABEL.get(method)):
        raise Proto15CursorContractError("TDG6 retry evidence identity differs")
    _require_sha("TDG6 initial state", evidence["initial_state_sha256"])
    if evidence["initial_state_sha256"] != state_sha256:
        raise Proto15CursorContractError("TDG6 retry state differs from cursor")
    if not _same_binary64(evidence["initial_time"], float.fromhex(str(accepted_time["binary64_hex"]))):
        raise Proto15CursorContractError("TDG6 retry initial time differs")
    attempted = evidence["attempted_macro_step_size"]
    retry = evidence["retry_step_size"]
    if (isinstance(attempted, bool) or not isinstance(attempted, (int, float))
            or isinstance(retry, bool) or not isinstance(retry, (int, float))
            or not math.isfinite(float(attempted)) or not math.isfinite(float(retry))
            or float(attempted) <= 0.0 or float(retry) <= 0.0
            or not _same_binary64(float(retry), float(attempted) * 0.5)):
        raise Proto15CursorContractError("TDG6 retry size is not the exact binary64 half step")
    if evidence["previous_step_index"] != step_index or evidence["previous_transaction_serial"] != transaction_serial:
        raise Proto15CursorContractError("TDG6 retry accepted serials differ")
    retry_count = _nonnegative_integer("TDG6 step retry count", evidence["retry_count_for_current_macro_step"])
    cumulative = _nonnegative_integer("TDG6 cumulative retry count", evidence["cumulative_temporal_retry_count"])
    if retry_count != prior_ledger.current_macro_step_temporal_retry_count + 1 or cumulative != prior_ledger.cumulative_temporal_retry_count + 1:
        raise Proto15CursorContractError("TDG6 retry counts do not extend the ledger")
    failed = evidence["failed_channels"]
    admissions = evidence["channel_admissions"]
    if not isinstance(failed, list) or not failed or any(item not in TDG6_COMPLETE_STATE_CHANNELS for item in failed):
        raise Proto15CursorContractError("TDG6 failed-channel evidence differs")
    if not isinstance(admissions, list) or len(admissions) != len(TDG6_COMPLETE_STATE_CHANNELS):
        raise Proto15CursorContractError("TDG6 channel-admission evidence differs")
    observed_channels: list[str] = []
    observed_failed: list[str] = []
    for item in admissions:
        if not isinstance(item, Mapping) or item.get("channel") not in TDG6_COMPLETE_STATE_CHANNELS:
            raise Proto15CursorContractError("TDG6 channel admission is malformed")
        channel = str(item["channel"])
        observed_channels.append(channel)
        if item.get("admission_passed") is False:
            observed_failed.append(channel)
    if tuple(observed_channels) != TDG6_COMPLETE_STATE_CHANNELS or observed_failed != failed:
        raise Proto15CursorContractError("TDG6 failed channels do not match channel admissions")
    for name in (
        "accepted_state_bitwise_preserved", "real_monitor_and_causal_ledgers_preserved",
        "real_tracer_ledger_preserved", "accumulated_debit_unchanged",
        "classification_is_numerical_not_physical",
    ):
        if evidence[name] is not True:
            raise Proto15CursorContractError("TDG6 retry evidence crossed its numerical boundary")
    return evidence


def _updated_retry_ledger(prior: TDG6TemporalLedger, evidence: Mapping[str, object]) -> TDG6TemporalLedger:
    serialized = _canonical(dict(evidence)).decode("utf-8")
    return TDG6TemporalLedger(
        last_accepted_time=prior.last_accepted_time,
        accepted_macro_step_count=prior.accepted_macro_step_count,
        cumulative_temporal_retry_count=int(evidence["cumulative_temporal_retry_count"]),
        current_macro_step_temporal_retry_count=int(evidence["retry_count_for_current_macro_step"]),
        last_accepted_macro_step_temporal_retry_count=prior.last_accepted_macro_step_temporal_retry_count,
        accumulated_debit_vector=prior.accumulated_debit_vector,
        serialized_temporal_rejections=(*prior.serialized_temporal_rejections, serialized),
    )


def _resolved_rollback_journal_fields(
    record: Mapping[str, object],
) -> dict[str, str]:
    """Resolve PROTO15's rollback journal names from the hash envelope.

    The record self-hash cannot be embedded inside the payload it hashes.  The
    sealed protocol's two ``new_journal_*`` names are therefore implemented by
    the immutable journal envelope and are checked as such during recovery.
    """
    try:
        parent = _require_sha(
            "rollback journal parent", record["journal_parent_sha256"]
        )
        digest = _require_sha("rollback journal record", record["record_sha256"])
    except (KeyError, Proto15CursorContractError) as cause:
        raise Proto15CheckpointIntegrityError(
            "rollback journal envelope mapping differs"
        ) from cause
    return {
        "new_journal_parent_sha256": parent,
        "new_journal_record_sha256": digest,
    }


def _prior_retry_ledger(
    updated: TDG6TemporalLedger, evidence: Mapping[str, object],
) -> TDG6TemporalLedger:
    """Invert exactly one TDG6 rejection append for persisted-payload checks."""
    retry_count = _nonnegative_integer(
        "TDG6 persisted retry count", evidence.get("retry_count_for_current_macro_step")
    )
    cumulative = _nonnegative_integer(
        "TDG6 persisted cumulative retry count", evidence.get("cumulative_temporal_retry_count")
    )
    if retry_count < 1 or cumulative < 1 or not updated.serialized_temporal_rejections:
        raise Proto15CursorContractError("persisted TDG6 rejection cannot have a predecessor ledger")
    return TDG6TemporalLedger(
        last_accepted_time=updated.last_accepted_time,
        accepted_macro_step_count=updated.accepted_macro_step_count,
        cumulative_temporal_retry_count=cumulative - 1,
        current_macro_step_temporal_retry_count=retry_count - 1,
        last_accepted_macro_step_temporal_retry_count=(
            updated.last_accepted_macro_step_temporal_retry_count
        ),
        accumulated_debit_vector=updated.accumulated_debit_vector,
        serialized_temporal_rejections=updated.serialized_temporal_rejections[:-1],
    )


def _validate_retry_payload(cursor: "Proto15Cursor", ledger: TDG6TemporalLedger) -> None:
    payload = cursor.payload["retry_successor_payload_or_none"]
    if not isinstance(payload, Mapping):
        raise Proto15CursorContractError("pending cursor omits retry payload")
    predecessor = _plan_mapping(payload["predecessor_plan"])
    successor = _plan_mapping(payload["successor_plan"])
    evidence = payload["TDG6_rejection_evidence"]
    if not isinstance(evidence, Mapping):
        raise Proto15CursorContractError("pending cursor retry evidence differs")
    prior = _prior_retry_ledger(ledger, evidence)
    evidence = _validate_retry_evidence_mapping(
        evidence,
        method=str(cursor.payload["method"]),
        state_sha256=str(cursor.payload["accepted_state_sha256"]),
        accepted_time=cursor.payload["accepted_boundary_time"],
        step_index=int(cursor.payload["previous_step_index"]),
        transaction_serial=int(cursor.payload["previous_transaction_serial"]),
        prior_ledger=prior,
    )
    if _updated_retry_ledger(prior, evidence) != ledger:
        raise Proto15CursorContractError("pending cursor evidence does not produce its ledger")
    prefix = payload["complete_rejection_prefix"]
    if not isinstance(prefix, list) or tuple(prefix) != ledger.serialized_temporal_rejections or not prefix:
        raise Proto15CursorContractError("pending cursor rejection prefix differs")
    if prefix[-1] != _canonical(dict(evidence)).decode("utf-8"):
        raise Proto15CursorContractError("pending cursor rejection tail differs")
    if (evidence.get("retry_count_for_current_macro_step") != ledger.current_macro_step_temporal_retry_count
            or evidence.get("cumulative_temporal_retry_count") != ledger.cumulative_temporal_retry_count):
        raise Proto15CursorContractError("pending cursor TDG6 evidence counts differ from ledger")
    if payload["member_key"] != cursor.payload["member_key"] or payload["method"] != cursor.payload["method"] or payload["point_count"] != cursor.payload["point_count"]:
        raise Proto15CursorContractError("pending cursor member identity differs")
    if payload["accepted_state_sha256"] != cursor.payload["accepted_state_sha256"] or payload["accepted_boundary_time"] != cursor.payload["accepted_boundary_time"]:
        raise Proto15CursorContractError("pending cursor accepted identity differs")
    if payload["previous_step_index"] != cursor.payload["previous_step_index"] or payload["previous_transaction_serial"] != cursor.payload["previous_transaction_serial"]:
        raise Proto15CursorContractError("pending cursor accepted serials differ")
    if payload["retry_count"] != ledger.current_macro_step_temporal_retry_count or ledger.current_macro_step_temporal_retry_count < 1:
        raise Proto15CursorContractError("pending cursor retry count differs")
    if payload["successor_attempt_id"] != cursor.payload["attempt_id"] or cursor.payload["attempt_serial"] < 1:
        raise Proto15CursorContractError("pending cursor successor attempt identity differs")
    predecessor_attempt = int(cursor.payload["attempt_serial"]) - 1
    expected_predecessor_id = f"{cursor.payload['campaign_id']}:{cursor.payload['member_key']}:{cursor.payload['committed_common_event_index']}:{predecessor_attempt}"
    if payload["predecessor_attempt_id"] != expected_predecessor_id:
        raise Proto15CursorContractError("pending cursor predecessor attempt identity differs")
    rejection = _require_sha("durable rejection record", payload["durable_rejection_record_sha256"])
    if rejection != cursor.payload["journal_tip_sha256"]:
        raise Proto15CursorContractError("pending cursor journal tip differs from rejection")
    _require_sha("predecessor journal tip", payload["predecessor_journal_tip_sha256"])
    half_cap = _finite_hex("pending half cap", payload["half_cap"], positive=True)
    minimum = _finite_hex("pending minimum step", payload["minimum_macro_step"], positive=True)
    if not _same_binary64(half_cap, evidence["retry_step_size"]) or not _same_binary64(minimum, float(TDG6_MINIMUM_MACRO_STEP)):
        raise Proto15CursorContractError("pending cursor cap or minimum differs")
    if not (_same_binary64(float.fromhex(predecessor["current_hex"]), evidence["initial_time"])
            and _same_binary64(float.fromhex(predecessor["macro_width_hex"]), evidence["attempted_macro_step_size"])
            and _same_binary64(float.fromhex(successor["current_hex"]), evidence["initial_time"])
            and _same_binary64(float.fromhex(successor["requested_cap_hex"]), evidence["retry_step_size"])
            and float.fromhex(successor["macro_width_hex"]) < float.fromhex(predecessor["macro_width_hex"])):
        raise Proto15CursorContractError("pending cursor TDG7 plan relation differs")
    expected_successor = _plan_mapping(plan_forward_proto14_subdivision(
        float(evidence["initial_time"]),
        float.fromhex(str(cursor.payload["event_target_time"]["binary64_hex"])),
        float(evidence["retry_step_size"]),
        minimum_width=float(TDG6_MINIMUM_MACRO_STEP),
    ))
    if successor != expected_successor:
        raise Proto15CursorContractError("pending cursor TDG7 successor is not the exact frozen plan")


@dataclass(frozen=True, slots=True)
class Proto15Cursor:
    payload: Mapping[str, object]

    @property
    def sha256(self) -> str:
        return str(self.payload["cursor_chain_sha256"])

    @property
    def mode(self) -> str:
        return str(self.payload["mode"])

    def validate(self) -> None:
        p = self.payload
        required = {"protocol_artifact_id", "campaign_id", "member_key", "method", "point_count",
                    "committed_common_event_index", "accepted_boundary_time", "accepted_state_sha256",
                    "previous_step_index", "previous_transaction_serial", "TDG6_ledger_sha256",
                    "cursor_generation", "event_target_time", "attempt_serial", "attempt_id", "mode",
                    "retry_successor_payload_or_none", "journal_tip_sha256", "cursor_chain_parent_sha256",
                    "cursor_chain_sha256"}
        if set(p) != required or p.get("member_key") not in MEMBER_KEYS or p.get("mode") not in CURSOR_MODES:
            raise Proto15CursorContractError("cursor fields or mode differ")
        if not isinstance(p["protocol_artifact_id"], str) or not p["protocol_artifact_id"]:
            raise Proto15CursorContractError("cursor protocol identity differs")
        if not isinstance(p["campaign_id"], str) or not p["campaign_id"]:
            raise Proto15CursorContractError("cursor campaign identity differs")
        expected_method, expected_points = str(p["member_key"]).split("-", 1)
        if p["method"] != expected_method or p["point_count"] != int(expected_points):
            raise Proto15CursorContractError("cursor method or point count differs from member key")
        _time_identity(p["accepted_boundary_time"]); _time_identity(p["event_target_time"])
        for name in ("accepted_state_sha256", "TDG6_ledger_sha256", "journal_tip_sha256",
                     "cursor_chain_parent_sha256", "cursor_chain_sha256"):
            _require_sha(name, p[name])
        for name in ("point_count", "committed_common_event_index", "previous_step_index",
                     "previous_transaction_serial", "cursor_generation", "attempt_serial"):
            _nonnegative_integer(f"cursor {name}", p[name])
        expected_id = f"{p['campaign_id']}:{p['member_key']}:{p['committed_common_event_index']}:{p['attempt_serial']}"
        if p["attempt_id"] != expected_id:
            raise Proto15CursorContractError("cursor attempt identity differs")
        retry = p["retry_successor_payload_or_none"]
        if p["mode"] == "FRESH_READY" and retry is not None:
            raise Proto15CursorContractError("fresh cursor has retry payload")
        if p["mode"] == "RETRY_PENDING" and not isinstance(retry, Mapping):
            raise Proto15CursorContractError("pending cursor omits retry payload")
        if p["mode"] == "RETRY_PENDING" and set(retry) != {
            "predecessor_plan", "successor_plan", "TDG6_rejection_evidence",
            "complete_rejection_prefix", "durable_rejection_record_sha256",
            "predecessor_journal_tip_sha256", "predecessor_attempt_id",
            "successor_attempt_id", "member_key", "method", "point_count",
            "accepted_state_sha256", "accepted_boundary_time", "previous_step_index",
            "previous_transaction_serial", "retry_count", "half_cap", "minimum_macro_step",
        }:
            raise Proto15CursorContractError("pending cursor retry payload fields differ")
        bare = dict(p); observed = bare.pop("cursor_chain_sha256")
        if _hash(bare) != observed:
            raise Proto15CursorContractError("cursor hash differs")


def _cursor(payload: Mapping[str, object]) -> Proto15Cursor:
    value = dict(payload)
    value.pop("cursor_chain_sha256", None)
    value["cursor_chain_sha256"] = _hash(value)
    answer = Proto15Cursor(value); answer.validate(); return answer


class Proto15SyntheticCampaign:
    """A synthetic campaign store with durable, fail-closed cursor transitions."""

    def __init__(self, root: Path, checkpoint: Mapping[str, object], *, fault_hook: FaultHook | None = None) -> None:
        self.root, self._checkpoint, self.fault_hook = Path(root), dict(checkpoint), fault_hook
        self._load_checkpoint(self._checkpoint)

    @classmethod
    def initialize(cls, root: Path, *, protocol_artifact_id: str, campaign_id: str,
                   members: Mapping[str, Mapping[str, object]], common_event_index: int = 23,
                   event_target_time: object = Fraction(3, 2), fault_hook: FaultHook | None = None) -> "Proto15SyntheticCampaign":
        if not isinstance(protocol_artifact_id, str) or not protocol_artifact_id:
            raise Proto15CursorContractError("protocol identity must be nonempty text")
        if not isinstance(campaign_id, str) or not campaign_id or ":" in campaign_id:
            raise Proto15CursorContractError("campaign identity must be nonempty text without colon")
        _nonnegative_integer("common event index", common_event_index)
        if set(members) != set(MEMBER_KEYS):
            raise Proto15CursorContractError("initialization requires exactly PROTO15's six members")
        root = Path(root)
        if root.exists() and any(root.iterdir()):
            raise FileExistsError(root)
        for name in ("states", "journal", "checkpoints"):
            (root / name).mkdir(parents=True, exist_ok=True)
        ledgers: dict[str, dict[str, object]] = {}; states: dict[str, str] = {}; cursors: dict[str, Mapping[str, object]] = {}
        targets = _time_identity(event_target_time)
        for key in MEMBER_KEYS:
            item = members[key]
            state = item.get("state")
            ledger = item.get("ledger")
            if not isinstance(ledger, TDG6TemporalLedger) or not isinstance(state, Mapping):
                raise TypeError("synthetic member requires mapping state and TDG6 ledger")
            method, points = key.split("-", 1)
            if item.get("method") != method or item.get("point_count") != int(points):
                raise Proto15CursorContractError("initial member method/grid identity differs")
            if ledger.current_macro_step_temporal_retry_count != 0:
                raise Proto15CursorContractError("initial FRESH_READY ledger has a current retry")
            step_index = _nonnegative_integer("initial step index", item.get("step_index", 0))
            transaction_serial = _nonnegative_integer("initial transaction serial", item.get("transaction_serial", 0))
            state_sha = cls._store_state_static(root, state, fault_hook)
            states[key] = state_sha; ledgers[key] = _ledger_mapping(ledger)
            accepted = _time_identity(item.get("accepted_boundary_time", ledger.last_accepted_time))
            if accepted["binary64_hex"] != ledger.last_accepted_time.hex():
                raise Proto15CursorContractError("initial accepted time differs from TDG6 ledger")
            cursor = _cursor({"protocol_artifact_id": protocol_artifact_id, "campaign_id": campaign_id,
                "member_key": key, "method": item["method"], "point_count": item["point_count"],
                "committed_common_event_index": common_event_index, "accepted_boundary_time": accepted,
                "accepted_state_sha256": state_sha, "previous_step_index": step_index,
                "previous_transaction_serial": transaction_serial,
                "TDG6_ledger_sha256": _hash(ledgers[key]), "cursor_generation": 0,
                "event_target_time": targets, "attempt_serial": 0,
                "attempt_id": f"{campaign_id}:{key}:{common_event_index}:0", "mode": "FRESH_READY",
                "retry_successor_payload_or_none": None, "journal_tip_sha256": ROOT_SHA256,
                "cursor_chain_parent_sha256": ROOT_SHA256})
            cursors[key] = cursor.payload
        checkpoint = cls._make_checkpoint(protocol_artifact_id, campaign_id, 0, common_event_index,
                                          targets, cursors, ledgers, states, ROOT_SHA256, False,
                                          {key: 0 for key in MEMBER_KEYS}, {key: 0 for key in MEMBER_KEYS}, ROOT_SHA256)
        cls._publish_checkpoint_static(root, checkpoint, fault_hook)
        return cls(root, checkpoint, fault_hook=fault_hook)

    @staticmethod
    def _hook(hook: FaultHook | None, phase: str, target: str) -> None:
        if hook is not None:
            hook(phase, target)

    @classmethod
    def _atomic_bytes(cls, path: Path, payload: bytes, hook: FaultHook | None, target: str) -> None:
        temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
        try:
            cls._hook(hook, "before_write", target)
            with temporary.open("xb") as handle:
                handle.write(payload); cls._hook(hook, "after_write", target)
                cls._hook(hook, "before_flush", target); handle.flush(); cls._hook(hook, "after_flush", target)
                cls._hook(hook, "before_fsync", target); os.fsync(handle.fileno()); cls._hook(hook, "after_fsync", target)
            cls._hook(hook, "before_replace", target); os.replace(temporary, path); cls._hook(hook, "after_replace", target)
            descriptor = os.open(path.parent, os.O_RDONLY)
            try:
                cls._hook(hook, "before_dir_fsync", target); os.fsync(descriptor); cls._hook(hook, "after_dir_fsync", target)
            finally: os.close(descriptor)
        except OSError as error:
            raise Proto15DurabilityError(f"durability failure for {target}") from error
        finally:
            temporary.unlink(missing_ok=True)

    @classmethod
    def _store_state_static(cls, root: Path, state: Mapping[str, object], hook: FaultHook | None) -> str:
        payload = _canonical(dict(state)); digest = sha256(payload).hexdigest(); path = root / "states" / f"{digest}.json"
        if path.exists():
            if path.read_bytes() != payload: raise Proto15CheckpointIntegrityError("content-addressed state collision")
        else: cls._atomic_bytes(path, payload, hook, "state_object")
        return digest

    @staticmethod
    def _sets(cursors: Mapping[str, Mapping[str, object]], ledgers: Mapping[str, Mapping[str, object]], states: Mapping[str, str]) -> tuple[str, str, str]:
        ordered = [{"member_key": key, "cursor_chain_sha256": cursors[key]["cursor_chain_sha256"],
                    "TDG6_ledger_sha256": cursors[key]["TDG6_ledger_sha256"], "accepted_state_sha256": states[key],
                    "accepted_boundary_time": cursors[key]["accepted_boundary_time"], "cursor_generation": cursors[key]["cursor_generation"],
                    "committed_common_event_index": cursors[key]["committed_common_event_index"]} for key in MEMBER_KEYS]
        return _hash(ordered), _hash({key: _hash(ledgers[key]) for key in MEMBER_KEYS}), _hash({key: states[key] for key in MEMBER_KEYS})

    @classmethod
    def _make_checkpoint(cls, protocol: str, campaign: str, generation: int, event: int, target: Mapping[str, str],
                         cursors: Mapping[str, Mapping[str, object]], ledgers: Mapping[str, Mapping[str, object]], states: Mapping[str, str],
                         journal_tip: str, terminal_lock: bool, attempt_high: Mapping[str, int], cursor_high: Mapping[str, int], parent: str,
                         disposition: str = "nonterminal", terminal_closure: Mapping[str, object] | None = None) -> dict[str, object]:
        cursor_set, ledger_set, state_set = cls._sets(cursors, ledgers, states)
        terminal_member = None if terminal_closure is None else terminal_closure.get("member_key")
        terminal_record = None if terminal_closure is None else journal_tip
        bare: dict[str, object] = {"protocol_artifact_id": protocol, "campaign_id": campaign, "campaign_generation": generation,
            "committed_common_event_index": event, "active_event_target_time": dict(target), "disposition": disposition,
            "terminal_member_key": terminal_member, "terminal_record_sha256": terminal_record,
            "member_keys": list(MEMBER_KEYS), "cursors": {k: dict(cursors[k]) for k in MEMBER_KEYS},
            "ledgers": {k: dict(ledgers[k]) for k in MEMBER_KEYS}, "states": dict(states), "cursor_set_sha256": cursor_set,
            "ledger_set_sha256": ledger_set, "accepted_state_set_sha256": state_set, "journal_tip_sha256": journal_tip,
            "terminal_lock": terminal_lock, "attempt_high_water": dict(attempt_high), "cursor_generation_high_water": dict(cursor_high),
            "parent_checkpoint_sha256": parent, "terminal_closure": None if terminal_closure is None else dict(terminal_closure)}
        bare["checkpoint_sha256"] = _hash(bare); return bare

    @classmethod
    def _publish_checkpoint_static(cls, root: Path, checkpoint: Mapping[str, object], hook: FaultHook | None) -> Path:
        payload = _canonical(dict(checkpoint)); digest = str(checkpoint["checkpoint_sha256"])
        path = root / "checkpoints" / f"{int(checkpoint['campaign_generation']):020d}-{digest}.json"
        if path.exists():
            if path.read_bytes() != payload: raise Proto15CheckpointIntegrityError("same generation checkpoint fork")
            return path
        # A second valid file with this generation is a fork, even if its name differs.
        if list((root / "checkpoints").glob(f"{int(checkpoint['campaign_generation']):020d}-*.json")):
            raise Proto15CheckpointIntegrityError("same generation checkpoint fork")
        cls._atomic_bytes(path, payload, hook, "checkpoint"); return path

    @staticmethod
    def _validate_state_object(root: Path, digest: str) -> None:
        digest = _require_sha("state object", digest)
        path = root / "states" / f"{digest}.json"
        if not path.is_file():
            raise Proto15CheckpointIntegrityError("referenced state object is missing")
        payload = path.read_bytes()
        if sha256(payload).hexdigest() != digest:
            raise Proto15CheckpointIntegrityError("state object content hash differs")
        value = _decode_canonical_json(
            payload, context="state object", error=Proto15CheckpointIntegrityError
        )
        if not isinstance(value, Mapping):
            raise Proto15CheckpointIntegrityError("state object is not a mapping")

    def _load_checkpoint(self, checkpoint: Mapping[str, object]) -> None:
        required = {
            "protocol_artifact_id", "campaign_id", "campaign_generation",
            "committed_common_event_index", "active_event_target_time", "disposition",
            "terminal_member_key", "terminal_record_sha256", "member_keys", "cursors",
            "ledgers", "states", "cursor_set_sha256", "ledger_set_sha256",
            "accepted_state_set_sha256", "journal_tip_sha256", "terminal_lock",
            "attempt_high_water", "cursor_generation_high_water",
            "parent_checkpoint_sha256", "terminal_closure", "checkpoint_sha256",
        }
        if not isinstance(checkpoint, Mapping) or set(checkpoint) != required:
            raise Proto15CheckpointIntegrityError("checkpoint fields differ")
        bare = dict(checkpoint); observed = bare.pop("checkpoint_sha256")
        try:
            observed = _require_sha("checkpoint", observed)
            parent = _require_sha("parent checkpoint", checkpoint["parent_checkpoint_sha256"])
        except Proto15CursorContractError as cause:
            raise Proto15CheckpointIntegrityError(str(cause)) from cause
        if _hash(bare) != observed or checkpoint["member_keys"] != list(MEMBER_KEYS):
            raise Proto15CheckpointIntegrityError("checkpoint hash or member order differs")
        protocol = checkpoint["protocol_artifact_id"]
        campaign = checkpoint["campaign_id"]
        if not isinstance(protocol, str) or not protocol or not isinstance(campaign, str) or not campaign:
            raise Proto15CheckpointIntegrityError("checkpoint protocol or campaign identity differs")
        generation = _nonnegative_integer("campaign generation", checkpoint["campaign_generation"], error=Proto15CheckpointIntegrityError)
        event = _nonnegative_integer("common event index", checkpoint["committed_common_event_index"], error=Proto15CheckpointIntegrityError)
        if (generation == 0) is not (parent == ROOT_SHA256):
            raise Proto15CheckpointIntegrityError("checkpoint root-parent relation differs")
        try:
            target = _time_identity(checkpoint["active_event_target_time"])
            journal_tip = _require_sha("journal tip", checkpoint["journal_tip_sha256"])
        except Proto15CursorContractError as cause:
            raise Proto15CheckpointIntegrityError(str(cause)) from cause
        disposition = checkpoint["disposition"]
        terminal_lock = checkpoint["terminal_lock"]
        terminal_member = checkpoint["terminal_member_key"]
        terminal_record = checkpoint["terminal_record_sha256"]
        closure = checkpoint["terminal_closure"]
        if disposition not in {"nonterminal", "rolled_back", "terminal_invalid"} or not isinstance(terminal_lock, bool):
            raise Proto15CheckpointIntegrityError("checkpoint disposition differs")
        if disposition == "terminal_invalid":
            if terminal_lock is not True or terminal_member not in MEMBER_KEYS or terminal_record != journal_tip or not isinstance(closure, Mapping):
                raise Proto15CheckpointIntegrityError("terminal checkpoint identity differs")
        elif terminal_lock or terminal_member is not None or terminal_record is not None or closure is not None:
            raise Proto15CheckpointIntegrityError("nonterminal checkpoint carries terminal state")
        cursors, ledgers, states = checkpoint["cursors"], checkpoint["ledgers"], checkpoint["states"]
        high_attempt, high_cursor = checkpoint["attempt_high_water"], checkpoint["cursor_generation_high_water"]
        if not all(isinstance(x, Mapping) and set(x) == set(MEMBER_KEYS) for x in (cursors, ledgers, states, high_attempt, high_cursor)):
            raise Proto15CheckpointIntegrityError("checkpoint does not contain all six member mappings")
        actual = self._sets(cursors, ledgers, states)
        if actual != (checkpoint["cursor_set_sha256"], checkpoint["ledger_set_sha256"], checkpoint["accepted_state_set_sha256"]):
            raise Proto15CheckpointIntegrityError("checkpoint member-set hashes differ")
        parsed_cursors: dict[str, Proto15Cursor] = {}
        parsed_ledgers: dict[str, TDG6TemporalLedger] = {}
        parsed_states: dict[str, str] = {}
        attempts: dict[str, int] = {}
        generations: dict[str, int] = {}
        for key in MEMBER_KEYS:
            try:
                cursor = Proto15Cursor(dict(cursors[key])); cursor.validate()
                state = _require_sha("state", states[key])
            except (TypeError, Proto15CursorContractError) as cause:
                raise Proto15CheckpointIntegrityError("checkpoint cursor or state differs") from cause
            ledger = _ledger_from_mapping(ledgers[key])
            if (cursor.payload["protocol_artifact_id"] != protocol
                    or cursor.payload["campaign_id"] != campaign
                    or cursor.payload["member_key"] != key
                    or cursor.payload["committed_common_event_index"] != event
                    or cursor.payload["event_target_time"] != target
                    or cursor.payload["accepted_state_sha256"] != state
                    or cursor.payload["accepted_boundary_time"]["binary64_hex"] != ledger.last_accepted_time.hex()):
                raise Proto15CheckpointIntegrityError("checkpoint cursor/member cross-binding differs")
            ledger_hash = _hash(_ledger_mapping(ledger))
            if key == terminal_member:
                if not isinstance(closure, Mapping) or closure.get("TDG6_ledger_sha256") != ledger_hash:
                    raise Proto15CheckpointIntegrityError("terminal member ledger binding differs")
            elif cursor.payload["TDG6_ledger_sha256"] != ledger_hash:
                raise Proto15CheckpointIntegrityError("cursor TDG6 ledger binding differs")
            if (cursor.mode == "FRESH_READY" and key != terminal_member
                    and ledger.current_macro_step_temporal_retry_count != 0):
                raise Proto15CheckpointIntegrityError("fresh cursor carries a pending TDG6 retry")
            if cursor.mode == "RETRY_PENDING" and key != terminal_member:
                try:
                    _validate_retry_payload(cursor, ledger)
                except Proto15CursorContractError as cause:
                    raise Proto15CheckpointIntegrityError("pending cursor payload differs") from cause
            attempt = _nonnegative_integer("attempt high-water", high_attempt[key], error=Proto15CheckpointIntegrityError)
            cursor_generation = _nonnegative_integer("cursor high-water", high_cursor[key], error=Proto15CheckpointIntegrityError)
            if attempt != cursor.payload["attempt_serial"] or cursor_generation != cursor.payload["cursor_generation"]:
                raise Proto15CheckpointIntegrityError("cursor high-water marks differ")
            self._validate_state_object(self.root, state)
            parsed_cursors[key] = cursor; parsed_ledgers[key] = ledger; parsed_states[key] = state
            attempts[key] = attempt; generations[key] = cursor_generation
        if disposition == "terminal_invalid":
            self._validate_terminal_closure(
                closure, protocol=protocol, campaign=campaign, generation=generation,
                event=event, member_key=str(terminal_member), journal_tip=journal_tip,
                cursor=parsed_cursors[str(terminal_member)], ledger=parsed_ledgers[str(terminal_member)],
                state_sha=parsed_states[str(terminal_member)],
            )
        self.protocol_artifact_id, self.campaign_id = protocol, campaign
        self.generation, self.event_index = generation, event
        self.event_target_time, self.journal_tip = target, journal_tip
        self.terminal_lock, self.cursors = terminal_lock, parsed_cursors
        self.ledgers, self.states = parsed_ledgers, parsed_states
        self.attempt_high, self.cursor_high = attempts, generations
        self.checkpoint_sha256 = observed

    @staticmethod
    def _validate_terminal_closure(
        closure: Mapping[str, object], *, protocol: str, campaign: str,
        generation: int, event: int, member_key: str, journal_tip: str,
        cursor: Proto15Cursor, ledger: TDG6TemporalLedger, state_sha: str,
    ) -> None:
        required = {
            "protocol_artifact_id", "campaign_id", "committed_common_event_index",
            "campaign_generation", "attempt_serial", "attempt_id", "member_key",
            "method", "point_count", "accepted_boundary_time", "accepted_state_sha256",
            "previous_step_index", "previous_transaction_serial", "TDG6_ledger_sha256",
            "exhaustion_updated_ledger", "durable_rejection_record_sha256",
            "predecessor_cursor_chain_sha256", "journal_tip_sha256", "exhaustion_reason",
            "classification", "terminal_lock",
        }
        if set(closure) != required:
            raise Proto15CheckpointIntegrityError("terminal closure fields differ")
        if (closure["protocol_artifact_id"] != protocol or closure["campaign_id"] != campaign
                or closure["committed_common_event_index"] != event
                or closure["campaign_generation"] != generation - 1
                or closure["member_key"] != member_key
                or closure["method"] != cursor.payload["method"]
                or closure["point_count"] != cursor.payload["point_count"]
                or closure["attempt_serial"] != cursor.payload["attempt_serial"]
                or closure["attempt_id"] != cursor.payload["attempt_id"]
                or closure["accepted_boundary_time"] != cursor.payload["accepted_boundary_time"]
                or closure["accepted_state_sha256"] != state_sha
                or closure["previous_step_index"] != cursor.payload["previous_step_index"]
                or closure["previous_transaction_serial"] != cursor.payload["previous_transaction_serial"]
                or closure["predecessor_cursor_chain_sha256"] != cursor.sha256
                or closure["classification"] != "invalid_implementation_or_nonconverged_run"
                or closure["terminal_lock"] is not True
                or closure["exhaustion_reason"] not in {"maximum_temporal_retries", "minimum_macro_step"}):
            raise Proto15CheckpointIntegrityError("terminal closure identity differs")
        try:
            rejection = _require_sha("terminal rejection", closure["durable_rejection_record_sha256"])
        except Proto15CursorContractError as cause:
            raise Proto15CheckpointIntegrityError("terminal rejection hash differs") from cause
        if rejection == ROOT_SHA256 or closure["journal_tip_sha256"] != rejection:
            raise Proto15CheckpointIntegrityError("terminal closure omits its rejection")
        if not isinstance(closure["exhaustion_updated_ledger"], Mapping):
            raise Proto15CheckpointIntegrityError("terminal closure ledger differs")
        embedded = _ledger_from_mapping(closure["exhaustion_updated_ledger"])
        if embedded != ledger or closure["TDG6_ledger_sha256"] != _hash(_ledger_mapping(ledger)):
            raise Proto15CheckpointIntegrityError("terminal closure updated ledger differs")

    def _append_record(self, kind: str, payload: Mapping[str, object]) -> str:
        if kind not in RECORD_KINDS: raise Proto15JournalIntegrityError("unknown journal kind")
        records = self._read_journal_chain(self.root)
        observed_tip = ROOT_SHA256 if not records else str(records[-1]["record_sha256"])
        if observed_tip != self.journal_tip:
            raise Proto15JournalIntegrityError("in-memory journal tip differs from durable chain")
        sequence = len(records) + 1
        inner = dict(payload); inner["sequence"] = sequence
        bare = {"record_kind": kind, "journal_parent_sha256": self.journal_tip, "payload": inner}
        record = dict(bare); record["record_sha256"] = _hash(bare)
        encoded = _canonical(record); frame = f"{len(encoded):08x} ".encode("ascii") + encoded + b"\n"
        path = self.root / "journal" / f"{sequence:020d}-{record['record_sha256']}.journal"
        self._atomic_bytes(path, frame, self.fault_hook, "journal_record")
        self.journal_tip = str(record["record_sha256"]); return self.journal_tip

    def _publish(self, *, disposition: str = "nonterminal", terminal_closure: Mapping[str, object] | None = None) -> None:
        checkpoint = self._make_checkpoint(self.protocol_artifact_id, self.campaign_id, self.generation + 1, self.event_index,
            self.event_target_time, {k: v.payload for k, v in self.cursors.items()}, {k: _ledger_mapping(v) for k, v in self.ledgers.items()}, self.states,
            self.journal_tip, self.terminal_lock, self.attempt_high, self.cursor_high, self.checkpoint_sha256, disposition, terminal_closure)
        self._publish_checkpoint_static(self.root, checkpoint, self.fault_hook); self._checkpoint = checkpoint; self._load_checkpoint(checkpoint)

    def guard_fresh_entrypoint(self, member_key: str) -> Proto15Cursor:
        if self.terminal_lock: raise Proto15TerminalLocked("terminal lock blocks fresh entrypoint")
        cursor = self.cursors[member_key]
        if cursor.mode != "FRESH_READY": raise Proto15CursorContractError("fresh entrypoint forbids RETRY_PENDING cursor")
        return cursor

    def guard_retry_entrypoint(self, member_key: str) -> Proto15Cursor:
        if self.terminal_lock: raise Proto15TerminalLocked("terminal lock blocks retry entrypoint")
        cursor = self.cursors[member_key]
        if cursor.mode != "RETRY_PENDING": raise Proto15CursorContractError("retry entrypoint requires RETRY_PENDING cursor")
        return cursor

    def _require_clean_checkpoint_tip(self) -> None:
        if self.journal_tip != self._checkpoint["journal_tip_sha256"]:
            raise Proto15JournalIntegrityError("campaign has an unresolved durable journal suffix")

    def durable_rejection_sink(
        self, member_key: str, *, predecessor_plan: TDG7Binary64SubdivisionPlan | Mapping[str, object],
    ) -> Callable[[Mapping[str, object]], None]:
        self.guard_fresh_entrypoint(member_key) if self.cursors[member_key].mode == "FRESH_READY" else self.guard_retry_entrypoint(member_key)
        self._require_clean_checkpoint_tip()
        persisted_plan = _plan_mapping(predecessor_plan)
        active = self.cursors[member_key]
        if active.mode == "RETRY_PENDING":
            pending = active.payload["retry_successor_payload_or_none"]
            if not isinstance(pending, Mapping) or persisted_plan != pending["successor_plan"]:
                raise Proto15CursorContractError(
                    "pending retry must execute its exact durable TDG7 successor plan"
                )
        def sink(evidence: Mapping[str, object]) -> None:
            cursor = self.cursors[member_key]
            normalized = _json_safe(dict(evidence))
            if not isinstance(normalized, Mapping):
                raise Proto15CursorContractError("TDG6 rejection evidence is not a mapping")
            _validate_retry_evidence_mapping(
                normalized, method=str(cursor.payload["method"]), state_sha256=self.states[member_key],
                accepted_time=cursor.payload["accepted_boundary_time"],
                step_index=int(cursor.payload["previous_step_index"]),
                transaction_serial=int(cursor.payload["previous_transaction_serial"]),
                prior_ledger=self.ledgers[member_key],
            )
            if (not _same_binary64(float.fromhex(persisted_plan["current_hex"]), normalized["initial_time"])
                    or not _same_binary64(float.fromhex(persisted_plan["macro_width_hex"]), normalized["attempted_macro_step_size"])):
                raise Proto15CursorContractError("TDG7 predecessor plan differs from TDG6 rejection")
            self._append_record("TDG6_REJECTION", {"member_key": member_key, "TDG6_rejection_evidence": dict(normalized),
                "predecessor_plan": persisted_plan,
                "predecessor_cursor_chain_sha256": self.cursors[member_key].sha256})
        return sink

    def _validated_tdg6_outcome(
        self, member_key: str,
        outcome: TDG6TemporalRetryRequired | TDG6TemporalRetryExhausted,
    ) -> tuple[Proto15Cursor, dict[str, object], TDG6TemporalLedger, str]:
        if not isinstance(outcome, (TDG6TemporalRetryRequired, TDG6TemporalRetryExhausted)):
            raise TypeError("outcome must be TDG6 typed retry outcome")
        if not isinstance(outcome.evidence, TDG6TemporalRetryEvidence):
            raise Proto15CursorContractError("typed TDG6 outcome omits canonical evidence type")
        cursor = self.cursors[member_key]
        evidence_value = _json_safe(outcome.evidence.as_mapping())
        if not isinstance(evidence_value, Mapping):
            raise Proto15CursorContractError("typed TDG6 evidence is not a mapping")
        evidence = _validate_retry_evidence_mapping(
            evidence_value, method=str(cursor.payload["method"]), state_sha256=self.states[member_key],
            accepted_time=cursor.payload["accepted_boundary_time"],
            step_index=int(cursor.payload["previous_step_index"]),
            transaction_serial=int(cursor.payload["previous_transaction_serial"]),
            prior_ledger=self.ledgers[member_key],
        )
        expected_ledger = _updated_retry_ledger(self.ledgers[member_key], evidence)
        if outcome.updated_ledger != expected_ledger:
            raise Proto15CursorContractError("typed TDG6 outcome ledger differs from durable evidence")
        exhausted = (int(evidence["retry_count_for_current_macro_step"]) > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
                     or float(evidence["retry_step_size"]) < float(TDG6_MINIMUM_MACRO_STEP))
        if isinstance(outcome, TDG6TemporalRetryRequired) and exhausted:
            raise Proto15CursorContractError("retry-required outcome is already exhausted")
        if isinstance(outcome, TDG6TemporalRetryExhausted):
            expected_reason = ("maximum_temporal_retries"
                               if int(evidence["retry_count_for_current_macro_step"]) > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
                               else "minimum_macro_step")
            if not exhausted or outcome.reason != expected_reason:
                raise Proto15CursorContractError("TDG6 exhaustion reason differs from evidence")
        rejection = self._find_tip_rejection(member_key, evidence)
        return cursor, evidence, expected_ledger, rejection

    def close_tdg6_outcome(self, member_key: str, outcome: TDG6TemporalRetryRequired | TDG6TemporalRetryExhausted,
                           *, retry_successor: Mapping[str, object] | None = None) -> Proto15Cursor:
        cursor, evidence, ledger, rejection = self._validated_tdg6_outcome(member_key, outcome)
        if isinstance(outcome, TDG6TemporalRetryExhausted):
            return self._commit_terminal_closure(member_key, outcome, rejection=rejection)
        if not isinstance(retry_successor, Mapping) or set(retry_successor) != {"predecessor_plan", "successor_plan"}:
            raise Proto15CursorContractError("retry-required outcome needs exact predecessor and successor TDG7 plans")
        predecessor_plan = _plan_mapping(retry_successor["predecessor_plan"])
        successor_plan = _plan_mapping(retry_successor["successor_plan"])
        rejection_record = self._read_record(sorted((self.root / "journal").glob("*.journal"))[-1])
        if rejection_record["payload"].get("predecessor_plan") != predecessor_plan:
            raise Proto15CursorContractError("retry successor predecessor plan differs from durable rejection")
        if not (_same_binary64(float.fromhex(predecessor_plan["current_hex"]), evidence["initial_time"])
                and _same_binary64(float.fromhex(predecessor_plan["macro_width_hex"]), evidence["attempted_macro_step_size"])
                and _same_binary64(float.fromhex(successor_plan["current_hex"]), evidence["initial_time"])
                and _same_binary64(float.fromhex(successor_plan["requested_cap_hex"]), evidence["retry_step_size"])
                and float.fromhex(successor_plan["macro_width_hex"]) < float.fromhex(predecessor_plan["macro_width_hex"])):
            raise Proto15CursorContractError("TDG7 retry plans differ from typed TDG6 evidence")
        expected_successor = _plan_mapping(plan_forward_proto14_subdivision(
            float(evidence["initial_time"]),
            float.fromhex(str(cursor.payload["event_target_time"]["binary64_hex"])),
            float(evidence["retry_step_size"]),
            minimum_width=float(TDG6_MINIMUM_MACRO_STEP),
        ))
        if successor_plan != expected_successor:
            raise Proto15CursorContractError("TDG7 retry successor is not the exact frozen plan")
        predecessor_journal_tip = str(rejection_record["journal_parent_sha256"])
        successor: dict[str, object] = {"predecessor_plan": predecessor_plan, "successor_plan": successor_plan}
        successor.update({"TDG6_rejection_evidence": evidence, "complete_rejection_prefix": list(ledger.serialized_temporal_rejections),
            "durable_rejection_record_sha256": rejection, "predecessor_journal_tip_sha256": predecessor_journal_tip,
            "predecessor_attempt_id": cursor.payload["attempt_id"], "member_key": member_key, "method": cursor.payload["method"],
            "point_count": cursor.payload["point_count"], "accepted_state_sha256": self.states[member_key],
            "accepted_boundary_time": cursor.payload["accepted_boundary_time"], "previous_step_index": cursor.payload["previous_step_index"],
            "previous_transaction_serial": cursor.payload["previous_transaction_serial"], "retry_count": ledger.current_macro_step_temporal_retry_count,
            "half_cap": outcome.retry_step_size.hex(), "minimum_macro_step": float(TDG6_MINIMUM_MACRO_STEP).hex()})
        next_attempt = self.attempt_high[member_key] + 1; successor["successor_attempt_id"] = f"{self.campaign_id}:{member_key}:{self.event_index}:{next_attempt}"
        new = _cursor({**cursor.payload, "TDG6_ledger_sha256": _hash(_ledger_mapping(ledger)), "cursor_generation": self.cursor_high[member_key] + 1,
            "attempt_serial": next_attempt, "attempt_id": successor["successor_attempt_id"], "mode": "RETRY_PENDING", "retry_successor_payload_or_none": successor,
            "journal_tip_sha256": rejection, "cursor_chain_parent_sha256": cursor.sha256})
        _validate_retry_payload(new, ledger)
        self._append_record("CURSOR_TRANSITION", {"member_key": member_key, "predecessor_cursor_chain_sha256": cursor.sha256,
            "durable_rejection_record_sha256": rejection, "successor_cursor_chain_sha256": new.sha256, "successor_attempt_id": new.payload["attempt_id"], "cursor": dict(new.payload)})
        self.cursors[member_key] = new; self.ledgers[member_key] = ledger; self.attempt_high[member_key] = next_attempt; self.cursor_high[member_key] += 1; self._publish(); return new

    def _find_tip_rejection(self, member_key: str, evidence: Mapping[str, object]) -> str:
        files = sorted((self.root / "journal").glob("*.journal"))
        if not files: raise Proto15JournalIntegrityError("typed outcome lacks durable rejection")
        record = self._read_record(files[-1])
        if (record["record_kind"] != "TDG6_REJECTION"
                or record["record_sha256"] != self.journal_tip
                or record["payload"].get("member_key") != member_key
                or _canonical(record["payload"].get("TDG6_rejection_evidence")) != _canonical(dict(evidence))):
            raise Proto15JournalIntegrityError("typed outcome differs from durable rejection")
        return str(record["record_sha256"])

    def accepted_fine_commit(
        self,
        member_key: str,
        *,
        state: Mapping[str, object],
        ledger: TDG6TemporalLedger,
        executed_plan: TDG7Binary64SubdivisionPlan | Mapping[str, object],
        accepted_boundary_time: object,
        step_index: int,
        transaction_serial: int,
    ) -> Proto15Cursor:
        if self.terminal_lock: raise Proto15TerminalLocked("terminal lock blocks accepted commit")
        self._require_clean_checkpoint_tip()
        current = self.cursors[member_key]
        if not isinstance(state, Mapping) or not isinstance(ledger, TDG6TemporalLedger):
            raise TypeError("accepted fine commit requires mapping state and TDG6 ledger")
        old_ledger = self.ledgers[member_key]
        accepted = _time_identity(accepted_boundary_time)
        old_time = float.fromhex(str(current.payload["accepted_boundary_time"]["binary64_hex"]))
        new_time = float.fromhex(accepted["binary64_hex"])
        target_time = float.fromhex(str(current.payload["event_target_time"]["binary64_hex"]))
        plan = _validated_accepted_plan(current, executed_plan, accepted)
        step_index = _nonnegative_integer("accepted step index", step_index)
        transaction_serial = _nonnegative_integer("accepted transaction serial", transaction_serial)
        if not old_time < new_time <= target_time:
            raise Proto15CursorContractError("accepted fine time does not advance within the event target")
        if (step_index <= current.payload["previous_step_index"]
                or transaction_serial <= current.payload["previous_transaction_serial"]):
            raise Proto15CursorContractError("accepted fine serials do not advance")
        if (ledger.last_accepted_time.hex() != accepted["binary64_hex"]
                or ledger.current_macro_step_temporal_retry_count != 0
                or ledger.accepted_macro_step_count != old_ledger.accepted_macro_step_count + 1
                or ledger.cumulative_temporal_retry_count != old_ledger.cumulative_temporal_retry_count
                or ledger.last_accepted_macro_step_temporal_retry_count != old_ledger.current_macro_step_temporal_retry_count
                or ledger.serialized_temporal_rejections != old_ledger.serialized_temporal_rejections
                or any(new < old for new, old in zip(ledger.accumulated_debit_vector, old_ledger.accumulated_debit_vector))):
            raise Proto15CursorContractError("accepted fine TDG6 ledger is not the exact monotone commit successor")
        state_sha = self._store_state_static(self.root, state, self.fault_hook); next_attempt = self.attempt_high[member_key] + 1
        new = _cursor({**current.payload, "accepted_boundary_time": accepted, "accepted_state_sha256": state_sha,
            "previous_step_index": step_index, "previous_transaction_serial": transaction_serial, "TDG6_ledger_sha256": _hash(_ledger_mapping(ledger)),
            "cursor_generation": self.cursor_high[member_key] + 1, "attempt_serial": next_attempt,
            "attempt_id": f"{self.campaign_id}:{member_key}:{self.event_index}:{next_attempt}", "mode": "FRESH_READY", "retry_successor_payload_or_none": None,
            "journal_tip_sha256": self.journal_tip, "cursor_chain_parent_sha256": current.sha256})
        self._append_record("ACCEPTED_FINE_COMMIT", {
            "member_key": member_key, "predecessor_cursor_chain_sha256": current.sha256,
            "successor_cursor_chain_sha256": new.sha256, "accepted_state_sha256": state_sha,
            "TDG6_ledger_sha256": new.payload["TDG6_ledger_sha256"],
            "executed_plan": plan, "cursor": dict(new.payload),
        })
        self.cursors[member_key] = new; self.ledgers[member_key] = ledger; self.states[member_key] = state_sha; self.attempt_high[member_key] = next_attempt; self.cursor_high[member_key] += 1; self._publish(); return new

    def terminal_closure(self, member_key: str, outcome: TDG6TemporalRetryExhausted) -> Proto15Cursor:
        if not isinstance(outcome, TDG6TemporalRetryExhausted):
            raise TypeError("terminal closure requires TDG6 exhaustion")
        _cursor_value, _evidence, _ledger, rejection = self._validated_tdg6_outcome(member_key, outcome)
        return self._commit_terminal_closure(member_key, outcome, rejection=rejection)

    def _commit_terminal_closure(
        self, member_key: str, outcome: TDG6TemporalRetryExhausted, *, rejection: str,
    ) -> Proto15Cursor:
        current = self.cursors[member_key]
        closure = {"protocol_artifact_id": self.protocol_artifact_id, "campaign_id": self.campaign_id, "committed_common_event_index": self.event_index,
            "campaign_generation": self.generation, "attempt_serial": current.payload["attempt_serial"], "attempt_id": current.payload["attempt_id"], "member_key": member_key,
            "method": current.payload["method"], "point_count": current.payload["point_count"], "accepted_boundary_time": current.payload["accepted_boundary_time"],
            "accepted_state_sha256": self.states[member_key], "previous_step_index": current.payload["previous_step_index"], "previous_transaction_serial": current.payload["previous_transaction_serial"],
            "TDG6_ledger_sha256": _hash(_ledger_mapping(outcome.updated_ledger)), "exhaustion_updated_ledger": _ledger_mapping(outcome.updated_ledger),
            "durable_rejection_record_sha256": rejection, "predecessor_cursor_chain_sha256": current.sha256, "journal_tip_sha256": rejection,
            "exhaustion_reason": outcome.reason, "classification": "invalid_implementation_or_nonconverged_run", "terminal_lock": True}
        self._append_record("TERMINAL_INVALID", closure)
        self.ledgers[member_key] = outcome.updated_ledger
        self.terminal_lock = True
        self._publish(disposition="terminal_invalid", terminal_closure=closure)
        return current

    def common_event_snapshot(self) -> Mapping[str, object]:
        self._require_clean_checkpoint_tip()
        event_times = {cursor.payload["accepted_boundary_time"]["binary64_hex"] for cursor in self.cursors.values()}
        if len(event_times) != 1:
            raise Proto15CursorContractError("common-event snapshot requires one shared accepted time")
        return {
            "campaign_id": self.campaign_id,
            "campaign_generation": self.generation,
            "committed_common_event_index": self.event_index,
            "event_time": next(iter(self.cursors.values())).payload["accepted_boundary_time"],
            "member_keys": list(MEMBER_KEYS),
            "cursor_set_sha256": self._checkpoint["cursor_set_sha256"],
            "ledger_set_sha256": self._checkpoint["ledger_set_sha256"],
            "accepted_state_set_sha256": self._checkpoint["accepted_state_set_sha256"],
            "journal_tip_sha256": self.journal_tip,
            "checkpoint_sha256": self.checkpoint_sha256,
            "checkpoint_payload": dict(self._checkpoint),
        }

    def common_event_rollback(self, snapshot: Mapping[str, object], *, aborted_attempt_id: str) -> None:
        if self.terminal_lock: raise Proto15TerminalLocked("terminal lock blocks rollback")
        self._require_clean_checkpoint_tip()
        required = {
            "campaign_id", "campaign_generation", "committed_common_event_index", "event_time",
            "member_keys", "cursor_set_sha256", "ledger_set_sha256", "accepted_state_set_sha256",
            "journal_tip_sha256", "checkpoint_sha256", "checkpoint_payload",
        }
        if not isinstance(snapshot, Mapping) or set(snapshot) != required or not isinstance(snapshot["checkpoint_payload"], Mapping):
            raise Proto15CheckpointIntegrityError("rollback snapshot fields differ")
        if not isinstance(aborted_attempt_id, str) or not aborted_attempt_id:
            raise Proto15CursorContractError("aborted attempt identity must be nonempty text")
        if aborted_attempt_id not in {
            cursor.payload["attempt_id"] for cursor in self.cursors.values()
        }:
            raise Proto15CursorContractError(
                "aborted attempt identity must name one active member attempt"
            )
        chain, _records = self._validated_persisted_chain(self.root)
        if chain[-1].checkpoint_sha256 != self.checkpoint_sha256:
            raise Proto15CheckpointIntegrityError("live campaign is not the highest persisted checkpoint")
        matches = [item for item in chain if item.checkpoint_sha256 == snapshot["checkpoint_sha256"]]
        if len(matches) != 1:
            raise Proto15CheckpointIntegrityError("rollback snapshot is not persisted in the campaign chain")
        restored = matches[0]
        if restored.generation >= self.generation or dict(snapshot["checkpoint_payload"]) != dict(restored._checkpoint):
            raise Proto15CheckpointIntegrityError("rollback snapshot is not an ancestral checkpoint")
        if (snapshot["campaign_id"] != restored.campaign_id
                or snapshot["campaign_generation"] != restored.generation
                or snapshot["committed_common_event_index"] != restored.event_index
                or snapshot["member_keys"] != list(MEMBER_KEYS)
                or snapshot["cursor_set_sha256"] != restored._checkpoint["cursor_set_sha256"]
                or snapshot["ledger_set_sha256"] != restored._checkpoint["ledger_set_sha256"]
                or snapshot["accepted_state_set_sha256"] != restored._checkpoint["accepted_state_set_sha256"]
                or snapshot["journal_tip_sha256"] != restored.journal_tip):
            raise Proto15CheckpointIntegrityError("rollback snapshot envelope differs from checkpoint")
        event_time = _time_identity(snapshot["event_time"])
        restored_times = {cursor.payload["accepted_boundary_time"]["binary64_hex"] for cursor in restored.cursors.values()}
        if len(restored_times) != 1 or event_time["binary64_hex"] not in restored_times:
            raise Proto15CheckpointIntegrityError("rollback snapshot event time differs")
        aborted_cursors = dict(self.cursors)
        old_sets = self._sets({k: v.payload for k, v in self.cursors.items()}, {k: _ledger_mapping(v) for k, v in self.ledgers.items()}, self.states)
        self.cursors, self.ledgers, self.states = restored.cursors, restored.ledgers, restored.states; self.event_index, self.event_target_time = restored.event_index, restored.event_target_time
        self.attempt_high = {k: max(self.attempt_high[k], restored.attempt_high[k]) + 1 for k in MEMBER_KEYS}; self.cursor_high = {k: max(self.cursor_high[k], restored.cursor_high[k]) + 1 for k in MEMBER_KEYS}
        for key in MEMBER_KEYS:
            base = self.cursors[key]; self.cursors[key] = _cursor({**base.payload, "cursor_generation": self.cursor_high[key], "attempt_serial": self.attempt_high[key],
                "attempt_id": f"{self.campaign_id}:{key}:{self.event_index}:{self.attempt_high[key]}", "cursor_chain_parent_sha256": aborted_cursors[key].sha256, "journal_tip_sha256": self.journal_tip})
        new_sets = self._sets({k: v.payload for k, v in self.cursors.items()}, {k: _ledger_mapping(v) for k, v in self.ledgers.items()}, self.states)
        if aborted_attempt_id in {cursor.payload["attempt_id"] for cursor in self.cursors.values()}:
            raise Proto15CursorContractError("rollback reused the aborted attempt identity")
        self._append_record("COMMON_EVENT_ROLLBACK", {"aborted_attempt_id": aborted_attempt_id, "aborted_campaign_generation": self.generation,
            "aborted_cursor_set_sha256": old_sets[0], "aborted_ledger_set_sha256": old_sets[1], "aborted_accepted_state_set_sha256": old_sets[2],
            "restored_snapshot_campaign_generation": restored.generation, "restored_snapshot_checkpoint_sha256": restored.checkpoint_sha256,
            "restored_snapshot_cursor_set_sha256": restored._checkpoint["cursor_set_sha256"], "restored_snapshot_ledger_set_sha256": restored._checkpoint["ledger_set_sha256"], "restored_snapshot_accepted_state_set_sha256": restored._checkpoint["accepted_state_set_sha256"],
            "new_campaign_generation": self.generation + 1, "new_cursor_set_sha256": new_sets[0], "new_ledger_set_sha256": new_sets[1], "new_accepted_state_set_sha256": new_sets[2]})
        self._publish(disposition="rolled_back")

    @staticmethod
    def _read_record(path: Path) -> Mapping[str, object]:
        raw = path.read_bytes()
        try:
            length_text, payload = raw.split(b" ", 1); encoded, newline = payload[:-1], payload[-1:]
            if (newline != b"\n" or len(length_text) != 8
                    or length_text.decode("ascii") != f"{len(encoded):08x}"):
                raise ValueError
            record = _decode_canonical_json(
                encoded, context="journal record", error=Proto15JournalIntegrityError
            )
        except (ValueError, UnicodeDecodeError) as error:
            raise Proto15JournalIntegrityError("journal frame is malformed") from error
        if not isinstance(record, Mapping) or set(record) != {"record_kind", "journal_parent_sha256", "payload", "record_sha256"}:
            raise Proto15JournalIntegrityError("journal record is noncanonical")
        bare = dict(record); observed = bare.pop("record_sha256")
        try:
            observed = _require_sha("journal record", observed)
            _require_sha("journal parent", record["journal_parent_sha256"])
        except Proto15CursorContractError as cause:
            raise Proto15JournalIntegrityError("journal hash field differs") from cause
        if record["record_kind"] not in RECORD_KINDS or _hash(bare) != observed or not isinstance(record["payload"], Mapping):
            raise Proto15JournalIntegrityError("journal record hash differs")
        name_prefix = path.name.split("-", 1)
        if len(name_prefix) != 2 or not name_prefix[0].isdigit() or name_prefix[1] != f"{observed}.journal":
            raise Proto15JournalIntegrityError("journal filename differs from record identity")
        sequence = _nonnegative_integer("journal sequence", record["payload"].get("sequence"), error=Proto15JournalIntegrityError)
        if sequence < 1 or int(name_prefix[0]) != sequence or name_prefix[0] != f"{sequence:020d}":
            raise Proto15JournalIntegrityError("journal sequence differs")
        return record

    @classmethod
    def _read_journal_chain(cls, root: Path) -> list[Mapping[str, object]]:
        records = [cls._read_record(path) for path in sorted((root / "journal").glob("*.journal"))]
        parent = ROOT_SHA256
        for expected, record in enumerate(records, 1):
            if record["payload"]["sequence"] != expected or record["journal_parent_sha256"] != parent:
                raise Proto15JournalIntegrityError("journal chain is noncontiguous")
            parent = str(record["record_sha256"])
        return records

    @classmethod
    def _read_checkpoint_file(cls, path: Path) -> Mapping[str, object]:
        try:
            encoded = path.read_bytes()
        except OSError as cause:
            raise Proto15CheckpointIntegrityError("checkpoint unreadable") from cause
        value = _decode_canonical_json(
            encoded, context="checkpoint", error=Proto15CheckpointIntegrityError
        )
        if not isinstance(value, Mapping):
            raise Proto15CheckpointIntegrityError("checkpoint is not a mapping")
        generation = value.get("campaign_generation")
        digest = value.get("checkpoint_sha256")
        if isinstance(generation, bool) or not isinstance(generation, int) or generation < 0 or not isinstance(digest, str):
            raise Proto15CheckpointIntegrityError("checkpoint filename identity differs")
        if path.name != f"{generation:020d}-{digest}.json":
            raise Proto15CheckpointIntegrityError("checkpoint filename differs from payload")
        return value

    @classmethod
    def _validate_persisted_rejection(
        cls,
        previous: "Proto15SyntheticCampaign",
        member: str,
        rejection: Mapping[str, object],
    ) -> tuple[dict[str, object], TDG6TemporalLedger, dict[str, object]]:
        payload = rejection["payload"]
        required = {
            "member_key", "TDG6_rejection_evidence", "predecessor_plan",
            "predecessor_cursor_chain_sha256", "sequence",
        }
        if (set(payload) != required or payload["member_key"] != member
                or payload["predecessor_cursor_chain_sha256"] != previous.cursors[member].sha256):
            raise Proto15CheckpointIntegrityError("persisted TDG6 rejection identity differs")
        cursor = previous.cursors[member]
        try:
            plan = _plan_mapping(payload["predecessor_plan"])
            evidence = _validate_retry_evidence_mapping(
                payload["TDG6_rejection_evidence"],
                method=str(cursor.payload["method"]),
                state_sha256=previous.states[member],
                accepted_time=cursor.payload["accepted_boundary_time"],
                step_index=int(cursor.payload["previous_step_index"]),
                transaction_serial=int(cursor.payload["previous_transaction_serial"]),
                prior_ledger=previous.ledgers[member],
            )
        except (KeyError, TypeError, ValueError, Proto15CursorContractError) as cause:
            raise Proto15CheckpointIntegrityError(
                "persisted TDG6 rejection evidence is invalid"
            ) from cause
        if (not _same_binary64(float.fromhex(plan["current_hex"]), evidence["initial_time"])
                or not _same_binary64(
                    float.fromhex(plan["macro_width_hex"]),
                    evidence["attempted_macro_step_size"],
                )):
            raise Proto15CheckpointIntegrityError(
                "persisted TDG6 rejection plan differs from its evidence"
            )
        if cursor.mode == "RETRY_PENDING":
            pending = cursor.payload["retry_successor_payload_or_none"]
            if not isinstance(pending, Mapping) or plan != pending["successor_plan"]:
                raise Proto15CheckpointIntegrityError(
                    "persisted retry did not execute the pending durable plan"
                )
        return evidence, _updated_retry_ledger(previous.ledgers[member], evidence), plan

    @staticmethod
    def _validate_accepted_successor(
        previous: "Proto15SyntheticCampaign",
        current: "Proto15SyntheticCampaign",
        member: str,
        accepted_record: Mapping[str, object],
    ) -> None:
        old_cursor = previous.cursors[member]
        new_cursor = current.cursors[member]
        old_ledger = previous.ledgers[member]
        new_ledger = current.ledgers[member]
        try:
            _validated_accepted_plan(
                old_cursor,
                accepted_record["payload"]["executed_plan"],
                new_cursor.payload["accepted_boundary_time"],
            )
        except (KeyError, TypeError, ValueError, Proto15CursorContractError) as cause:
            raise Proto15CheckpointIntegrityError(
                "persisted accepted TDG7 plan differs"
            ) from cause
        old_time = float.fromhex(str(old_cursor.payload["accepted_boundary_time"]["binary64_hex"]))
        new_time = float.fromhex(str(new_cursor.payload["accepted_boundary_time"]["binary64_hex"]))
        target = float.fromhex(str(new_cursor.payload["event_target_time"]["binary64_hex"]))
        if not old_time < new_time <= target:
            raise Proto15CheckpointIntegrityError("persisted accepted time is not monotone")
        if (new_cursor.payload["previous_step_index"] <= old_cursor.payload["previous_step_index"]
                or new_cursor.payload["previous_transaction_serial"] <= old_cursor.payload["previous_transaction_serial"]):
            raise Proto15CheckpointIntegrityError("persisted accepted serials are not monotone")
        if (new_cursor.mode != "FRESH_READY"
                or new_cursor.payload["retry_successor_payload_or_none"] is not None
                or new_cursor.payload["cursor_generation"] != old_cursor.payload["cursor_generation"] + 1
                or new_cursor.payload["attempt_serial"] != old_cursor.payload["attempt_serial"] + 1
                or new_cursor.payload["cursor_chain_parent_sha256"] != old_cursor.sha256
                or new_cursor.payload["journal_tip_sha256"] != accepted_record["journal_parent_sha256"]):
            raise Proto15CheckpointIntegrityError("persisted accepted cursor transition differs")
        if (new_ledger.last_accepted_time.hex()
                    != new_cursor.payload["accepted_boundary_time"]["binary64_hex"]
                or new_ledger.current_macro_step_temporal_retry_count != 0
                or new_ledger.accepted_macro_step_count != old_ledger.accepted_macro_step_count + 1
                or new_ledger.cumulative_temporal_retry_count != old_ledger.cumulative_temporal_retry_count
                or new_ledger.last_accepted_macro_step_temporal_retry_count
                    != old_ledger.current_macro_step_temporal_retry_count
                or new_ledger.serialized_temporal_rejections != old_ledger.serialized_temporal_rejections
                or any(new < old for new, old in zip(
                    new_ledger.accumulated_debit_vector,
                    old_ledger.accumulated_debit_vector,
                ))):
            raise Proto15CheckpointIntegrityError("persisted accepted TDG6 ledger differs")

    @classmethod
    def _validate_generation_transition(
        cls, previous: "Proto15SyntheticCampaign",
        current: "Proto15SyntheticCampaign",
        records: list[Mapping[str, object]],
        ancestry: Mapping[str, "Proto15SyntheticCampaign"],
    ) -> None:
        kinds = [record["record_kind"] for record in records]
        changed_cursors = [key for key in MEMBER_KEYS if previous.cursors[key].sha256 != current.cursors[key].sha256]
        changed_ledgers = [key for key in MEMBER_KEYS if previous.ledgers[key] != current.ledgers[key]]
        changed_states = [key for key in MEMBER_KEYS if previous.states[key] != current.states[key]]
        if current._checkpoint["disposition"] == "terminal_invalid":
            if kinds != ["TDG6_REJECTION", "TERMINAL_INVALID"] or changed_cursors or changed_states:
                raise Proto15CheckpointIntegrityError("terminal checkpoint transition differs")
            member = current._checkpoint["terminal_member_key"]
            if changed_ledgers != [member]:
                raise Proto15CheckpointIntegrityError("terminal checkpoint changed the wrong ledger set")
            rejection, terminal = records
            evidence, updated, _plan = cls._validate_persisted_rejection(
                previous, str(member), rejection
            )
            retry_count = int(evidence["retry_count_for_current_macro_step"])
            retry_size = float(evidence["retry_step_size"])
            exhausted = (
                retry_count > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
                or retry_size < float(TDG6_MINIMUM_MACRO_STEP)
            )
            reason = (
                "maximum_temporal_retries"
                if retry_count > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
                else "minimum_macro_step"
            )
            if (not exhausted or current.ledgers[str(member)] != updated
                    or terminal["journal_parent_sha256"] != rejection["record_sha256"]):
                raise Proto15CheckpointIntegrityError("terminal journal transition differs")
            terminal_payload = dict(terminal["payload"]); terminal_payload.pop("sequence", None)
            if (terminal_payload != current._checkpoint["terminal_closure"]
                    or terminal_payload.get("durable_rejection_record_sha256") != rejection["record_sha256"]
                    or terminal_payload.get("journal_tip_sha256") != rejection["record_sha256"]
                    or terminal_payload.get("exhaustion_reason") != reason):
                raise Proto15CheckpointIntegrityError("terminal record differs from checkpoint closure")
            return
        if current._checkpoint["disposition"] == "rolled_back":
            if kinds != ["COMMON_EVENT_ROLLBACK"] or set(changed_cursors) != set(MEMBER_KEYS):
                raise Proto15CheckpointIntegrityError("rollback checkpoint transition differs")
            record = records[0]
            payload = record["payload"]
            required = {
                "aborted_attempt_id", "aborted_campaign_generation", "aborted_cursor_set_sha256",
                "aborted_ledger_set_sha256", "aborted_accepted_state_set_sha256",
                "restored_snapshot_campaign_generation", "restored_snapshot_checkpoint_sha256",
                "restored_snapshot_cursor_set_sha256", "restored_snapshot_ledger_set_sha256",
                "restored_snapshot_accepted_state_set_sha256", "new_campaign_generation",
                "new_cursor_set_sha256", "new_ledger_set_sha256", "new_accepted_state_set_sha256",
                "sequence",
            }
            if set(payload) != required or payload["aborted_campaign_generation"] != previous.generation:
                raise Proto15CheckpointIntegrityError("rollback record fields differ")
            resolved_journal = _resolved_rollback_journal_fields(record)
            if (resolved_journal["new_journal_parent_sha256"] != previous.journal_tip
                    or resolved_journal["new_journal_record_sha256"] != current.journal_tip):
                raise Proto15CheckpointIntegrityError(
                    "rollback journal envelope mapping differs"
                )
            ancestor_sha = payload["restored_snapshot_checkpoint_sha256"]
            restored = ancestry.get(str(ancestor_sha))
            if restored is None or restored.generation >= previous.generation:
                raise Proto15CheckpointIntegrityError("rollback record does not name a strict ancestor")
            if (payload["aborted_cursor_set_sha256"] != previous._checkpoint["cursor_set_sha256"]
                    or payload["aborted_ledger_set_sha256"] != previous._checkpoint["ledger_set_sha256"]
                    or payload["aborted_accepted_state_set_sha256"] != previous._checkpoint["accepted_state_set_sha256"]
                    or payload["new_campaign_generation"] != current.generation
                    or payload["new_cursor_set_sha256"] != current._checkpoint["cursor_set_sha256"]
                    or payload["new_ledger_set_sha256"] != current._checkpoint["ledger_set_sha256"]
                    or payload["new_accepted_state_set_sha256"] != current._checkpoint["accepted_state_set_sha256"]
                    or payload["restored_snapshot_campaign_generation"] != restored.generation
                    or payload["restored_snapshot_cursor_set_sha256"] != restored._checkpoint["cursor_set_sha256"]
                    or payload["restored_snapshot_ledger_set_sha256"] != restored._checkpoint["ledger_set_sha256"]
                    or payload["restored_snapshot_accepted_state_set_sha256"] != restored._checkpoint["accepted_state_set_sha256"]
                    or current.event_index != restored.event_index
                    or current.event_target_time != restored.event_target_time
                    or current.states != restored.states
                    or current.ledgers != restored.ledgers):
                raise Proto15CheckpointIntegrityError("rollback record set binding differs")
            aborted_ids = {cursor.payload["attempt_id"] for cursor in previous.cursors.values()}
            if payload["aborted_attempt_id"] not in aborted_ids:
                raise Proto15CheckpointIntegrityError("rollback aborted attempt identity differs")
            for key in MEMBER_KEYS:
                base = restored.cursors[key].payload
                new = current.cursors[key].payload
                ignored = {
                    "cursor_generation", "attempt_serial", "attempt_id",
                    "cursor_chain_parent_sha256", "cursor_chain_sha256", "journal_tip_sha256",
                }
                if ({name: value for name, value in new.items() if name not in ignored}
                        != {name: value for name, value in base.items() if name not in ignored}
                        or new["cursor_generation"] != previous.cursor_high[key] + 1
                        or new["attempt_serial"] != previous.attempt_high[key] + 1
                        or new["cursor_chain_parent_sha256"] != previous.cursors[key].sha256
                        or new["journal_tip_sha256"] != record["journal_parent_sha256"]
                        or new["attempt_id"] in aborted_ids):
                    raise Proto15CheckpointIntegrityError("rollback cursor restoration differs")
            return
        if current._checkpoint["disposition"] != "nonterminal":
            raise Proto15CheckpointIntegrityError("unknown checkpoint transition disposition")
        if kinds == ["ACCEPTED_FINE_COMMIT"]:
            if (len(changed_cursors) != 1 or changed_ledgers != changed_cursors
                    or changed_states not in ([], changed_cursors)):
                raise Proto15CheckpointIntegrityError("accepted fine checkpoint did not change exactly one member")
            member = changed_cursors[0]
            payload = records[0]["payload"]
            required = {
                "member_key", "predecessor_cursor_chain_sha256", "successor_cursor_chain_sha256",
                "accepted_state_sha256", "TDG6_ledger_sha256", "executed_plan",
                "cursor", "sequence",
            }
            if (set(payload) != required or payload["member_key"] != member
                    or payload["predecessor_cursor_chain_sha256"] != previous.cursors[member].sha256
                    or payload["successor_cursor_chain_sha256"] != current.cursors[member].sha256
                    or payload["accepted_state_sha256"] != current.states[member]
                    or payload["TDG6_ledger_sha256"] != current.cursors[member].payload["TDG6_ledger_sha256"]
                    or payload["cursor"] != current.cursors[member].payload):
                raise Proto15CheckpointIntegrityError("accepted fine record differs from checkpoint")
            cls._validate_accepted_successor(previous, current, member, records[0])
            return
        if kinds == ["TDG6_REJECTION", "CURSOR_TRANSITION"]:
            if len(changed_cursors) != 1 or changed_ledgers != changed_cursors or changed_states:
                raise Proto15CheckpointIntegrityError("retry checkpoint changed the wrong member sets")
            member = changed_cursors[0]
            rejection, transition = records
            evidence, updated, predecessor_plan = cls._validate_persisted_rejection(
                previous, member, rejection
            )
            payload = transition["payload"]
            required = {
                "member_key", "predecessor_cursor_chain_sha256", "durable_rejection_record_sha256",
                "successor_cursor_chain_sha256", "successor_attempt_id", "cursor", "sequence",
            }
            if (set(payload) != required or rejection["payload"].get("member_key") != member
                    or payload["member_key"] != member
                    or payload["predecessor_cursor_chain_sha256"] != previous.cursors[member].sha256
                    or payload["durable_rejection_record_sha256"] != rejection["record_sha256"]
                    or transition["journal_parent_sha256"] != rejection["record_sha256"]
                    or payload["successor_cursor_chain_sha256"] != current.cursors[member].sha256
                    or payload["successor_attempt_id"] != current.cursors[member].payload["attempt_id"]
                    or payload["cursor"] != current.cursors[member].payload
                    or current.cursors[member].payload["cursor_generation"] != previous.cursors[member].payload["cursor_generation"] + 1
                    or current.cursors[member].payload["attempt_serial"] != previous.cursors[member].payload["attempt_serial"] + 1):
                raise Proto15CheckpointIntegrityError("retry transition record differs from checkpoint")
            pending = current.cursors[member].payload["retry_successor_payload_or_none"]
            if (not isinstance(pending, Mapping)
                    or pending["predecessor_plan"] != predecessor_plan
                    or pending["TDG6_rejection_evidence"] != evidence
                    or pending["predecessor_journal_tip_sha256"] != rejection["journal_parent_sha256"]
                    or pending["durable_rejection_record_sha256"] != rejection["record_sha256"]
                    or current.ledgers[member] != updated):
                raise Proto15CheckpointIntegrityError("retry evidence differs across durable records")
            return
        raise Proto15CheckpointIntegrityError("checkpoint publication has an unauthorized journal transition")

    @classmethod
    def _validated_persisted_chain(
        cls, root: Path, *, fault_hook: FaultHook | None = None,
    ) -> tuple[list["Proto15SyntheticCampaign"], list[Mapping[str, object]]]:
        root = Path(root)
        paths = sorted((root / "checkpoints").glob("*.json"))
        if not paths:
            raise Proto15CheckpointIntegrityError("no checkpoint exists")
        candidates = [cls._read_checkpoint_file(path) for path in paths]
        by_generation: dict[int, list[Mapping[str, object]]] = {}
        for candidate in candidates:
            by_generation.setdefault(int(candidate["campaign_generation"]), []).append(candidate)
        if any(len(values) != 1 for values in by_generation.values()):
            raise Proto15CheckpointIntegrityError("same generation checkpoint fork")
        if sorted(by_generation) != list(range(max(by_generation) + 1)):
            raise Proto15CheckpointIntegrityError("checkpoint generation chain is incomplete")
        records = cls._read_journal_chain(root)
        record_positions = {str(record["record_sha256"]): index for index, record in enumerate(records)}
        chain: list[Proto15SyntheticCampaign] = []
        previous: Proto15SyntheticCampaign | None = None
        prior_tip_position = -1
        for generation in range(max(by_generation) + 1):
            campaign = cls(root, by_generation[generation][0], fault_hook=fault_hook)
            if previous is None:
                if campaign.generation != 0 or campaign._checkpoint["parent_checkpoint_sha256"] != ROOT_SHA256:
                    raise Proto15CheckpointIntegrityError("checkpoint genesis differs")
            else:
                if (campaign.generation != previous.generation + 1
                        or campaign._checkpoint["parent_checkpoint_sha256"] != previous.checkpoint_sha256
                        or campaign.protocol_artifact_id != previous.protocol_artifact_id
                        or campaign.campaign_id != previous.campaign_id):
                    raise Proto15CheckpointIntegrityError("checkpoint parent chain differs")
                for key in MEMBER_KEYS:
                    old = previous.cursors[key]
                    new = campaign.cursors[key]
                    if new.sha256 == old.sha256:
                        if (new.payload["cursor_generation"] != old.payload["cursor_generation"]
                                or new.payload["attempt_serial"] != old.payload["attempt_serial"]):
                            raise Proto15CheckpointIntegrityError("unchanged cursor changed its high-water identity")
                    elif (new.payload["cursor_chain_parent_sha256"] != old.sha256
                          or new.payload["cursor_generation"] <= old.payload["cursor_generation"]
                          or new.payload["attempt_serial"] <= old.payload["attempt_serial"]):
                        raise Proto15CheckpointIntegrityError("cursor transition is not monotone from its parent")
                    if (campaign.attempt_high[key] < previous.attempt_high[key]
                            or campaign.cursor_high[key] < previous.cursor_high[key]):
                        raise Proto15CheckpointIntegrityError("checkpoint high-water mark regressed")
            tip = campaign.journal_tip
            if tip == ROOT_SHA256:
                tip_position = -1
            elif tip in record_positions:
                tip_position = record_positions[tip]
            else:
                raise Proto15CheckpointIntegrityError("checkpoint journal tip is absent")
            if tip_position < prior_tip_position:
                raise Proto15CheckpointIntegrityError("checkpoint journal tip regressed")
            if previous is None:
                if tip_position != -1:
                    raise Proto15CheckpointIntegrityError("genesis checkpoint has a journal history")
            else:
                cls._validate_generation_transition(
                    previous,
                    campaign,
                    records[prior_tip_position + 1:tip_position + 1],
                    {item.checkpoint_sha256: item for item in chain},
                )
            if campaign._checkpoint["disposition"] == "terminal_invalid":
                if tip_position < 0 or records[tip_position]["record_kind"] != "TERMINAL_INVALID":
                    raise Proto15CheckpointIntegrityError("terminal checkpoint tip differs")
            prior_tip_position = tip_position
            chain.append(campaign); previous = campaign
        return chain, records

    @classmethod
    def recover_campaign(cls, root: Path, *, fault_hook: FaultHook | None = None) -> "Proto15SyntheticCampaign":
        root = Path(root)
        chain, records = cls._validated_persisted_chain(root, fault_hook=fault_hook)
        campaign = chain[-1]
        if campaign.journal_tip == ROOT_SHA256:
            start = 0
        else:
            positions = [index for index, record in enumerate(records) if record["record_sha256"] == campaign.journal_tip]
            if len(positions) != 1:
                raise Proto15RecoveryInvalidStop("checkpoint journal tip is ambiguous")
            start = positions[0] + 1
        suffix = records[start:]
        if not suffix:
            return campaign
        if campaign.terminal_lock:
            raise Proto15RecoveryInvalidStop("terminal checkpoint has an uncheckpointed suffix")
        if suffix[0]["record_kind"] != "TDG6_REJECTION" or len(suffix) not in {1, 2}:
            raise Proto15RecoveryInvalidStop("uncheckpointed journal suffix is not an exact rejection branch")
        rejection = suffix[0]
        payload = rejection["payload"]
        if set(payload) != {"member_key", "TDG6_rejection_evidence", "predecessor_plan", "predecessor_cursor_chain_sha256", "sequence"}:
            raise Proto15RecoveryInvalidStop("durable rejection record fields differ")
        member = payload["member_key"]
        if member not in MEMBER_KEYS:
            raise Proto15RecoveryInvalidStop("rejection member differs")
        member = str(member)
        cursor = campaign.cursors[member]
        if payload["predecessor_cursor_chain_sha256"] != cursor.sha256:
            raise Proto15RecoveryInvalidStop("rejection cursor parent differs")
        try:
            predecessor_plan = _plan_mapping(payload["predecessor_plan"])
            evidence = _validate_retry_evidence_mapping(
                payload["TDG6_rejection_evidence"], method=str(cursor.payload["method"]),
                state_sha256=campaign.states[member], accepted_time=cursor.payload["accepted_boundary_time"],
                step_index=int(cursor.payload["previous_step_index"]),
                transaction_serial=int(cursor.payload["previous_transaction_serial"]),
                prior_ledger=campaign.ledgers[member],
            )
            if (not _same_binary64(float.fromhex(predecessor_plan["current_hex"]), evidence["initial_time"])
                    or not _same_binary64(float.fromhex(predecessor_plan["macro_width_hex"]), evidence["attempted_macro_step_size"])):
                raise Proto15CursorContractError("rejection predecessor plan differs")
            updated = _updated_retry_ledger(campaign.ledgers[member], evidence)
        except (TypeError, ValueError, Proto15CursorContractError) as cause:
            raise Proto15RecoveryInvalidStop("durable rejection evidence cannot reconstruct its branch") from cause
        retry_count = int(evidence["retry_count_for_current_macro_step"])
        retry_size = float(evidence["retry_step_size"])
        exhausted = retry_count > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP or retry_size < float(TDG6_MINIMUM_MACRO_STEP)
        reason = ("maximum_temporal_retries" if retry_count > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
                  else "minimum_macro_step")

        def terminal_payload() -> dict[str, object]:
            return {
                "protocol_artifact_id": campaign.protocol_artifact_id,
                "campaign_id": campaign.campaign_id,
                "committed_common_event_index": campaign.event_index,
                "campaign_generation": campaign.generation,
                "attempt_serial": cursor.payload["attempt_serial"],
                "attempt_id": cursor.payload["attempt_id"],
                "member_key": member,
                "method": cursor.payload["method"],
                "point_count": cursor.payload["point_count"],
                "accepted_boundary_time": cursor.payload["accepted_boundary_time"],
                "accepted_state_sha256": campaign.states[member],
                "previous_step_index": cursor.payload["previous_step_index"],
                "previous_transaction_serial": cursor.payload["previous_transaction_serial"],
                "TDG6_ledger_sha256": _hash(_ledger_mapping(updated)),
                "exhaustion_updated_ledger": _ledger_mapping(updated),
                "durable_rejection_record_sha256": rejection["record_sha256"],
                "predecessor_cursor_chain_sha256": cursor.sha256,
                "journal_tip_sha256": rejection["record_sha256"],
                "exhaustion_reason": reason,
                "classification": "invalid_implementation_or_nonconverged_run",
                "terminal_lock": True,
            }

        if exhausted:
            closure = terminal_payload()
            if len(suffix) == 2:
                terminal = suffix[1]
                terminal_without_sequence = dict(terminal["payload"]); terminal_without_sequence.pop("sequence", None)
                if terminal["record_kind"] != "TERMINAL_INVALID" or terminal_without_sequence != closure:
                    raise Proto15RecoveryInvalidStop("terminal suffix differs from exact exhausted branch")
                campaign.journal_tip = str(terminal["record_sha256"])
            else:
                campaign.journal_tip = str(rejection["record_sha256"])
                campaign._append_record("TERMINAL_INVALID", closure)
            campaign.ledgers[member] = updated
            campaign.terminal_lock = True
            campaign._publish(disposition="terminal_invalid", terminal_closure=closure)
            return campaign

        try:
            successor_plan = _plan_mapping(plan_forward_proto14_subdivision(
                float(evidence["initial_time"]),
                float.fromhex(str(cursor.payload["event_target_time"]["binary64_hex"])),
                retry_size,
                minimum_width=float(TDG6_MINIMUM_MACRO_STEP),
            ))
        except (TypeError, ValueError) as cause:
            raise Proto15RecoveryInvalidStop("TDG7 retry successor cannot be reconstructed") from cause
        next_attempt = campaign.attempt_high[member] + 1
        successor = {
            "predecessor_plan": predecessor_plan,
            "successor_plan": successor_plan,
            "TDG6_rejection_evidence": dict(evidence),
            "complete_rejection_prefix": list(updated.serialized_temporal_rejections),
            "durable_rejection_record_sha256": rejection["record_sha256"],
            "predecessor_journal_tip_sha256": rejection["journal_parent_sha256"],
            "predecessor_attempt_id": cursor.payload["attempt_id"],
            "successor_attempt_id": f"{campaign.campaign_id}:{member}:{campaign.event_index}:{next_attempt}",
            "member_key": member,
            "method": cursor.payload["method"],
            "point_count": cursor.payload["point_count"],
            "accepted_state_sha256": campaign.states[member],
            "accepted_boundary_time": cursor.payload["accepted_boundary_time"],
            "previous_step_index": cursor.payload["previous_step_index"],
            "previous_transaction_serial": cursor.payload["previous_transaction_serial"],
            "retry_count": retry_count,
            "half_cap": retry_size.hex(),
            "minimum_macro_step": float(TDG6_MINIMUM_MACRO_STEP).hex(),
        }
        pending = _cursor({
            **cursor.payload,
            "TDG6_ledger_sha256": _hash(_ledger_mapping(updated)),
            "cursor_generation": campaign.cursor_high[member] + 1,
            "attempt_serial": next_attempt,
            "attempt_id": successor["successor_attempt_id"],
            "mode": "RETRY_PENDING",
            "retry_successor_payload_or_none": successor,
            "journal_tip_sha256": rejection["record_sha256"],
            "cursor_chain_parent_sha256": cursor.sha256,
        })
        try:
            _validate_retry_payload(pending, updated)
        except Proto15CursorContractError as cause:
            raise Proto15RecoveryInvalidStop("reconstructed pending cursor differs") from cause
        expected_transition = {
            "member_key": member,
            "predecessor_cursor_chain_sha256": cursor.sha256,
            "durable_rejection_record_sha256": rejection["record_sha256"],
            "successor_cursor_chain_sha256": pending.sha256,
            "successor_attempt_id": pending.payload["attempt_id"],
            "cursor": dict(pending.payload),
        }
        if len(suffix) == 2:
            transition = suffix[1]
            transition_without_sequence = dict(transition["payload"]); transition_without_sequence.pop("sequence", None)
            if transition["record_kind"] != "CURSOR_TRANSITION" or transition_without_sequence != expected_transition:
                raise Proto15RecoveryInvalidStop("cursor transition suffix differs from exact retry branch")
            campaign.journal_tip = str(transition["record_sha256"])
        else:
            campaign.journal_tip = str(rejection["record_sha256"])
            campaign._append_record("CURSOR_TRANSITION", expected_transition)
        campaign.cursors[member] = pending
        campaign.ledgers[member] = updated
        campaign.attempt_high[member] = next_attempt
        campaign.cursor_high[member] += 1
        campaign._publish()
        return campaign


def recover_campaign(root: Path, *, fault_hook: FaultHook | None = None) -> Proto15SyntheticCampaign:
    return Proto15SyntheticCampaign.recover_campaign(root, fault_hook=fault_hook)


_VISIBLE_AFTER_FAULT: Final[frozenset[str]] = frozenset({
    "after_replace", "before_dir_fsync", "after_dir_fsync",
})


def _synthetic_members() -> dict[str, dict[str, object]]:
    return {
        key: {
            "method": key.split("-", 1)[0],
            "point_count": int(key.split("-", 1)[1]),
            "state": {"synthetic_member": key, "version": 0},
            "ledger": TDG6TemporalLedger.zero(initial_time=23 / 16),
            "step_index": 1,
            "transaction_serial": 2,
        }
        for key in MEMBER_KEYS
    }


def _synthetic_failed_admission(method: str) -> object:
    outer = _magnitude_interval(
        np.ones((2, 4), dtype=np.float64),
        np.zeros((2, 4), dtype=np.float64),
        subinterval_count=2,
        owned_row_count=1,
    )
    finest = _magnitude_interval(
        np.ones((4, 4), dtype=np.float64),
        np.zeros((4, 4), dtype=np.float64),
        subinterval_count=4,
        owned_row_count=1,
    )
    return classify_tdg6_runtime_intervals(
        method=method,
        outer_intervals={channel: outer for channel in TDG6_COMPLETE_STATE_CHANNELS},
        finest_intervals={channel: finest for channel in TDG6_COMPLETE_STATE_CHANNELS},
    )


def _synthetic_plan(
    campaign: Proto15SyntheticCampaign,
    member: str,
    requested_cap: float = 0.25,
) -> TDG7Binary64SubdivisionPlan:
    cursor = campaign.cursors[member]
    return plan_forward_proto14_subdivision(
        float.fromhex(str(cursor.payload["accepted_boundary_time"]["binary64_hex"])),
        float.fromhex(str(cursor.payload["event_target_time"]["binary64_hex"])),
        requested_cap,
        minimum_width=float(TDG6_MINIMUM_MACRO_STEP),
    )


def _synthetic_retry_outcome(
    campaign: Proto15SyntheticCampaign,
    member: str,
    predecessor_plan: TDG7Binary64SubdivisionPlan,
) -> TDG6TemporalRetryRequired | TDG6TemporalRetryExhausted:
    cursor = campaign.cursors[member]
    prior = campaign.ledgers[member]
    admission = _synthetic_failed_admission(
        RUNTIME_METHOD_BY_LABEL[str(cursor.payload["method"])]
    )
    evidence = TDG6TemporalRetryEvidence(
        event_type="rejected_TDG6_temporal_admission",
        method=RUNTIME_METHOD_BY_LABEL[str(cursor.payload["method"])],
        initial_time=float.fromhex(str(predecessor_plan.current.hex())),
        attempted_macro_step_size=predecessor_plan.macro_width,
        retry_step_size=predecessor_plan.macro_width * 0.5,
        initial_state_sha256=campaign.states[member],
        previous_step_index=int(cursor.payload["previous_step_index"]),
        previous_transaction_serial=int(cursor.payload["previous_transaction_serial"]),
        retry_count_for_current_macro_step=(
            prior.current_macro_step_temporal_retry_count + 1
        ),
        cumulative_temporal_retry_count=prior.cumulative_temporal_retry_count + 1,
        failed_channels=admission.failed_channels,
        channel_admissions=admission.channel_admissions,
    )
    updated = _updated_retry_ledger(prior, evidence.as_mapping())
    if evidence.retry_count_for_current_macro_step > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP:
        return TDG6TemporalRetryExhausted(
            reason="maximum_temporal_retries", evidence=evidence, updated_ledger=updated
        )
    if evidence.retry_step_size < float(TDG6_MINIMUM_MACRO_STEP):
        return TDG6TemporalRetryExhausted(
            reason="minimum_macro_step", evidence=evidence, updated_ledger=updated
        )
    return TDG6TemporalRetryRequired(evidence, updated)


def _synthetic_accept(campaign: Proto15SyntheticCampaign, member: str = "RK4-2049") -> None:
    prior = campaign.ledgers[member]
    cursor = campaign.cursors[member]
    if cursor.mode == "RETRY_PENDING":
        retry = cursor.payload["retry_successor_payload_or_none"]
        if not isinstance(retry, Mapping):
            raise AssertionError("synthetic pending cursor lost its retry plan")
        executed_plan: TDG7Binary64SubdivisionPlan | Mapping[str, object] = retry[
            "successor_plan"
        ]
        new_time = float.fromhex(str(executed_plan["boundaries_hex"][-1]))
    else:
        executed_plan = _synthetic_plan(campaign, member, 0.25)
        new_time = executed_plan.endpoint
    ledger = TDG6TemporalLedger(
        last_accepted_time=new_time,
        accepted_macro_step_count=prior.accepted_macro_step_count + 1,
        cumulative_temporal_retry_count=prior.cumulative_temporal_retry_count,
        current_macro_step_temporal_retry_count=0,
        last_accepted_macro_step_temporal_retry_count=(
            prior.current_macro_step_temporal_retry_count
        ),
        accumulated_debit_vector=prior.accumulated_debit_vector,
        serialized_temporal_rejections=prior.serialized_temporal_rejections,
    )
    campaign.accepted_fine_commit(
        member,
        state={"synthetic_member": member, "version": prior.accepted_macro_step_count + 1},
        ledger=ledger,
        executed_plan=executed_plan,
        accepted_boundary_time=new_time,
        step_index=int(campaign.cursors[member].payload["previous_step_index"]) + 1,
        transaction_serial=(
            int(campaign.cursors[member].payload["previous_transaction_serial"]) + 1
        ),
    )


def _synthetic_reject(
    campaign: Proto15SyntheticCampaign,
    member: str = "RK4-2049",
    *,
    requested_cap: float = 0.25,
) -> tuple[TDG6TemporalRetryRequired | TDG6TemporalRetryExhausted, TDG7Binary64SubdivisionPlan]:
    predecessor = _synthetic_plan(campaign, member, requested_cap)
    outcome = _synthetic_retry_outcome(campaign, member, predecessor)
    campaign.durable_rejection_sink(member, predecessor_plan=predecessor)(
        outcome.evidence.as_mapping()
    )
    return outcome, predecessor


def _synthetic_close_retry(
    campaign: Proto15SyntheticCampaign,
    member: str,
    outcome: TDG6TemporalRetryRequired,
    predecessor: TDG7Binary64SubdivisionPlan,
) -> None:
    successor = plan_forward_proto14_subdivision(
        outcome.evidence.initial_time,
        float.fromhex(str(campaign.cursors[member].payload["event_target_time"]["binary64_hex"])),
        outcome.retry_step_size,
        minimum_width=float(TDG6_MINIMUM_MACRO_STEP),
    )
    campaign.close_tdg6_outcome(
        member,
        outcome,
        retry_successor={"predecessor_plan": predecessor, "successor_plan": successor},
    )


def _synthetic_recovery_label(root: Path) -> str:
    try:
        campaign = recover_campaign(root)
    except Proto15RecoveryInvalidStop:
        return "explicit_recovery_invalid_stop"
    except (Proto15CheckpointIntegrityError, Proto15JournalIntegrityError):
        return "explicit_integrity_invalid_stop"
    if campaign.terminal_lock:
        return "terminal_invalid"
    if campaign._checkpoint["disposition"] == "rolled_back":
        return "rolled_back"
    if any(cursor.mode == "RETRY_PENDING" for cursor in campaign.cursors.values()):
        return "retry_pending"
    return f"generation_{campaign.generation}_fresh"


def _fault_once(phase: str, target: str) -> tuple[FaultHook, dict[str, bool]]:
    fired = {"value": False}

    def hook(observed_phase: str, observed_target: str) -> None:
        if (not fired["value"] and observed_phase == phase and observed_target == target):
            fired["value"] = True
            raise OSError("deterministic HLT13 durability fault")

    return hook, fired


def _directory_image(path: Path) -> dict[str, bytes]:
    return {item.name: item.read_bytes() for item in sorted(path.iterdir())}


def _restore_directory_image(path: Path, image: Mapping[str, bytes]) -> None:
    for item in path.iterdir():
        if not item.is_file():
            raise AssertionError("synthetic durability image contains a directory")
        item.unlink()
    for name, payload in image.items():
        (path / name).write_bytes(payload)


def _fault_expected(scenario: str, *, visible: bool) -> str:
    return {
        "state_object": "generation_0_fresh",
        "accepted_journal": (
            "explicit_recovery_invalid_stop" if visible else "generation_0_fresh"
        ),
        "accepted_checkpoint": (
            "generation_1_fresh" if visible else "explicit_recovery_invalid_stop"
        ),
        "rejection_journal": "retry_pending" if visible else "generation_0_fresh",
        "cursor_transition": "retry_pending",
        "terminal_record": "terminal_invalid",
        "terminal_checkpoint": "terminal_invalid",
        "rollback_record": (
            "explicit_recovery_invalid_stop" if visible else "generation_1_fresh"
        ),
        "rollback_checkpoint": (
            "rolled_back" if visible else "explicit_recovery_invalid_stop"
        ),
    }[scenario]


def _run_fault_case(scenario: str, phase: str) -> dict[str, object]:
    target = "state_object" if scenario == "state_object" else (
        "checkpoint" if scenario.endswith("checkpoint") else "journal_record"
    )
    visible = phase in _VISIBLE_AFTER_FAULT
    with TemporaryDirectory(prefix="hlt13-fault-") as temporary:
        root = Path(temporary) / "campaign"
        campaign = Proto15SyntheticCampaign.initialize(
            root,
            protocol_artifact_id="FGC-2-SF1-PROTO15",
            campaign_id="qualification",
            members=_synthetic_members(),
        )
        snapshot: Mapping[str, object] | None = None
        outcome: TDG6TemporalRetryRequired | TDG6TemporalRetryExhausted | None = None
        predecessor: TDG7Binary64SubdivisionPlan | None = None
        if scenario in {"cursor_transition", "terminal_record", "terminal_checkpoint"}:
            cap = float(TDG6_MINIMUM_MACRO_STEP) if scenario.startswith("terminal") else 0.25
            outcome, predecessor = _synthetic_reject(campaign, requested_cap=cap)
        elif scenario in {"rollback_record", "rollback_checkpoint"}:
            snapshot = campaign.common_event_snapshot()
            _synthetic_accept(campaign)
        target_directory = {
            "state_object": "states",
            "journal_record": "journal",
            "checkpoint": "checkpoints",
        }[target]
        old_image = _directory_image(root / target_directory)
        hook, fired = _fault_once(phase, target)
        campaign.fault_hook = hook
        try:
            if scenario in {"state_object", "accepted_journal", "accepted_checkpoint"}:
                _synthetic_accept(campaign)
            elif scenario == "rejection_journal":
                _synthetic_reject(campaign)
            elif scenario == "cursor_transition":
                if not isinstance(outcome, TDG6TemporalRetryRequired) or predecessor is None:
                    raise AssertionError("synthetic retry setup differs")
                _synthetic_close_retry(campaign, "RK4-2049", outcome, predecessor)
            elif scenario in {"terminal_record", "terminal_checkpoint"}:
                if not isinstance(outcome, TDG6TemporalRetryExhausted):
                    raise AssertionError("synthetic terminal setup differs")
                campaign.close_tdg6_outcome("RK4-2049", outcome)
            elif scenario in {"rollback_record", "rollback_checkpoint"}:
                if snapshot is None:
                    raise AssertionError("synthetic rollback setup differs")
                aborted = str(campaign.cursors["RK4-2049"].payload["attempt_id"])
                campaign.common_event_rollback(snapshot, aborted_attempt_id=aborted)
            else:
                raise AssertionError(f"unknown synthetic fault scenario {scenario}")
        except Proto15DurabilityError:
            pass
        else:
            raise AssertionError("synthetic fault hook did not interrupt publication")
        if not fired["value"]:
            raise AssertionError("synthetic fault hook was not reached")
        if phase in _AMBIGUOUS_REPLACE_CRASH_PHASES:
            old_root = Path(temporary) / "old-image"
            new_root = Path(temporary) / "new-image"
            shutil.copytree(root, old_root)
            shutil.copytree(root, new_root)
            _restore_directory_image(old_root / target_directory, old_image)
            observed_old = _synthetic_recovery_label(old_root)
            observed_new = _synthetic_recovery_label(new_root)
            expected_old = _fault_expected(scenario, visible=False)
            expected_new = _fault_expected(scenario, visible=True)
            return {
                "case": f"fault:{scenario}:{phase}",
                "phase": phase,
                "target": target,
                "visible_after_fault": "old_or_new_directory_image",
                "crash_images_tested": ["old", "new"],
                "observed": {"old": observed_old, "new": observed_new},
                "expected": {"old": expected_old, "new": expected_new},
                "passed": observed_old == expected_old and observed_new == expected_new,
            }
        observed = _synthetic_recovery_label(root)
    expected = _fault_expected(scenario, visible=visible)
    return {
        "case": f"fault:{scenario}:{phase}",
        "phase": phase,
        "target": target,
        "visible_after_fault": visible,
        "observed": observed,
        "expected": expected,
        "passed": observed == expected,
    }


def _run_nominal_cases() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []

    def record(name: str, passed: bool, observed: str) -> None:
        records.append({"case": name, "observed": observed, "passed": bool(passed)})

    with TemporaryDirectory(prefix="hlt13-nominal-") as temporary:
        root = Path(temporary) / "campaign"
        campaign = Proto15SyntheticCampaign.initialize(
            root, protocol_artifact_id="FGC-2-SF1-PROTO15",
            campaign_id="nominal-init", members=_synthetic_members(),
        )
        recovered = recover_campaign(root)
        record("nominal:initialize", recovered.generation == 0 and len(recovered.cursors) == 6,
               _synthetic_recovery_label(root))
    with TemporaryDirectory(prefix="hlt13-nominal-") as temporary:
        root = Path(temporary) / "campaign"
        campaign = Proto15SyntheticCampaign.initialize(
            root, protocol_artifact_id="FGC-2-SF1-PROTO15",
            campaign_id="nominal-accept", members=_synthetic_members(),
        )
        _synthetic_accept(campaign)
        record("nominal:accepted_fine", _synthetic_recovery_label(root) == "generation_1_fresh",
               _synthetic_recovery_label(root))
    with TemporaryDirectory(prefix="hlt13-nominal-") as temporary:
        root = Path(temporary) / "campaign"
        campaign = Proto15SyntheticCampaign.initialize(
            root, protocol_artifact_id="FGC-2-SF1-PROTO15",
            campaign_id="nominal-rejection", members=_synthetic_members(),
        )
        _synthetic_reject(campaign)
        label = _synthetic_recovery_label(root)
        record("nominal:rejection_suffix", label == "retry_pending", label)
    with TemporaryDirectory(prefix="hlt13-nominal-") as temporary:
        root = Path(temporary) / "campaign"
        campaign = Proto15SyntheticCampaign.initialize(
            root, protocol_artifact_id="FGC-2-SF1-PROTO15",
            campaign_id="nominal-retry", members=_synthetic_members(),
        )
        outcome, predecessor = _synthetic_reject(campaign)
        if not isinstance(outcome, TDG6TemporalRetryRequired):
            raise AssertionError("nominal retry unexpectedly exhausted")
        _synthetic_close_retry(campaign, "RK4-2049", outcome, predecessor)
        pending_label = _synthetic_recovery_label(root)
        recovered = recover_campaign(root)
        _synthetic_accept(recovered)
        accepted_label = _synthetic_recovery_label(root)
        record(
            "nominal:retry_pending",
            pending_label == "retry_pending" and accepted_label == "generation_2_fresh",
            f"{pending_label}->{accepted_label}",
        )
    with TemporaryDirectory(prefix="hlt13-nominal-") as temporary:
        root = Path(temporary) / "campaign"
        campaign = Proto15SyntheticCampaign.initialize(
            root, protocol_artifact_id="FGC-2-SF1-PROTO15",
            campaign_id="nominal-terminal", members=_synthetic_members(),
        )
        outcome, _predecessor = _synthetic_reject(
            campaign, requested_cap=float(TDG6_MINIMUM_MACRO_STEP)
        )
        if not isinstance(outcome, TDG6TemporalRetryExhausted):
            raise AssertionError("nominal terminal did not exhaust")
        campaign.close_tdg6_outcome("RK4-2049", outcome)
        label = _synthetic_recovery_label(root)
        record("nominal:terminal", label == "terminal_invalid", label)
    with TemporaryDirectory(prefix="hlt13-nominal-") as temporary:
        root = Path(temporary) / "campaign"
        campaign = Proto15SyntheticCampaign.initialize(
            root, protocol_artifact_id="FGC-2-SF1-PROTO15",
            campaign_id="nominal-rollback", members=_synthetic_members(),
        )
        snapshot = campaign.common_event_snapshot()
        _synthetic_accept(campaign)
        rollback_parent = campaign.journal_tip
        aborted = str(campaign.cursors["RK4-2049"].payload["attempt_id"])
        campaign.common_event_rollback(snapshot, aborted_attempt_id=aborted)
        label = _synthetic_recovery_label(root)
        rollback_record = Proto15SyntheticCampaign._read_record(
            sorted((root / "journal").glob("*.journal"))[-1]
        )
        rollback_fields = _resolved_rollback_journal_fields(rollback_record)
        mapping_valid = (
            rollback_fields["new_journal_parent_sha256"] == rollback_parent
            and rollback_fields["new_journal_record_sha256"]
            == rollback_record["record_sha256"]
        )
        record("nominal:rollback", label == "rolled_back" and mapping_valid, label)
    return records


def _synthetic_record_at(root: Path, sequence: int) -> dict[str, object]:
    path = next((root / "journal").glob(f"{sequence:020d}-*.journal"))
    return dict(Proto15SyntheticCampaign._read_record(path))


def _rewrite_synthetic_record(
    root: Path, record: Mapping[str, object],
) -> dict[str, object]:
    bare = {
        "record_kind": record["record_kind"],
        "journal_parent_sha256": record["journal_parent_sha256"],
        "payload": record["payload"],
    }
    rewritten = {**bare, "record_sha256": _hash(bare)}
    encoded = _canonical(rewritten)
    sequence = int(rewritten["payload"]["sequence"])
    for path in (root / "journal").glob(f"{sequence:020d}-*.journal"):
        path.unlink()
    path = root / "journal" / f"{sequence:020d}-{rewritten['record_sha256']}.journal"
    path.write_bytes(f"{len(encoded):08x} ".encode("ascii") + encoded + b"\n")
    return rewritten


def _synthetic_checkpoint_at(root: Path, generation: int) -> dict[str, object]:
    path = next((root / "checkpoints").glob(f"{generation:020d}-*.json"))
    value = _decode_canonical_json(
        path.read_bytes(), context="synthetic checkpoint",
        error=Proto15CheckpointIntegrityError,
    )
    if not isinstance(value, Mapping):
        raise AssertionError("synthetic checkpoint is not a mapping")
    return dict(value)


def _rewrite_synthetic_checkpoint(
    root: Path, checkpoint: dict[str, object],
) -> None:
    cursor_set, ledger_set, state_set = Proto15SyntheticCampaign._sets(
        checkpoint["cursors"], checkpoint["ledgers"], checkpoint["states"]
    )
    checkpoint["cursor_set_sha256"] = cursor_set
    checkpoint["ledger_set_sha256"] = ledger_set
    checkpoint["accepted_state_set_sha256"] = state_set
    bare = dict(checkpoint)
    bare.pop("checkpoint_sha256", None)
    checkpoint["checkpoint_sha256"] = _hash(bare)
    generation = int(checkpoint["campaign_generation"])
    for path in (root / "checkpoints").glob(f"{generation:020d}-*.json"):
        path.unlink()
    path = root / "checkpoints" / f"{generation:020d}-{checkpoint['checkpoint_sha256']}.json"
    path.write_bytes(_canonical(checkpoint))


def _run_mutation_cases() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []

    def invalid(
        name: str,
        mutate: Callable[[Path, Proto15SyntheticCampaign], None],
    ) -> None:
        with TemporaryDirectory(prefix="hlt13-mutation-") as temporary:
            root = Path(temporary) / "campaign"
            campaign = Proto15SyntheticCampaign.initialize(
                root, protocol_artifact_id="FGC-2-SF1-PROTO15",
                campaign_id=f"mutation-{name}", members=_synthetic_members(),
            )
            mutate(root, campaign)
            label = _synthetic_recovery_label(root)
            records.append({
                "case": f"mutation:{name}", "observed": label,
                "passed": label.startswith("explicit_"),
            })

    def remove_state(root: Path, _campaign: Proto15SyntheticCampaign) -> None:
        next((root / "states").glob("*.json")).unlink()

    def corrupt_state(root: Path, _campaign: Proto15SyntheticCampaign) -> None:
        next((root / "states").glob("*.json")).write_bytes(b"{\"corrupt\":true}")

    def fork_checkpoint(root: Path, _campaign: Proto15SyntheticCampaign) -> None:
        original = next((root / "checkpoints").glob("*.json"))
        duplicate = original.with_name(f"{0:020d}-{'f' * 64}.json")
        duplicate.write_bytes(original.read_bytes())

    def pending_plan_substitution(
        root: Path, campaign: Proto15SyntheticCampaign,
    ) -> None:
        outcome, predecessor = _synthetic_reject(campaign)
        if not isinstance(outcome, TDG6TemporalRetryRequired):
            raise AssertionError("synthetic pending-plan case exhausted")
        _synthetic_close_retry(campaign, "RK4-2049", outcome, predecessor)
        successor = campaign.cursors["RK4-2049"].payload[
            "retry_successor_payload_or_none"
        ]["successor_plan"]
        substitute = plan_forward_proto14_subdivision(
            outcome.evidence.initial_time,
            float.fromhex(str(campaign.cursors["RK4-2049"].payload[
                "event_target_time"
            ]["binary64_hex"])),
            np.nextafter(outcome.retry_step_size, float("inf")),
            minimum_width=float(TDG6_MINIMUM_MACRO_STEP),
        )
        if substitute.endpoint.hex() != successor["boundaries_hex"][-1]:
            raise AssertionError("synthetic substitute plan endpoint differs")
        _synthetic_accept(campaign)
        accepted = _synthetic_record_at(root, 3)
        accepted["payload"]["executed_plan"] = _plan_mapping(substitute)
        accepted = _rewrite_synthetic_record(root, accepted)
        checkpoint = _synthetic_checkpoint_at(root, 2)
        checkpoint["journal_tip_sha256"] = accepted["record_sha256"]
        _rewrite_synthetic_checkpoint(root, checkpoint)

    def retry_evidence_fork(
        root: Path, campaign: Proto15SyntheticCampaign,
    ) -> None:
        alternative_plan = _synthetic_plan(campaign, "RK4-2049", 0.03125)
        alternative = _synthetic_retry_outcome(
            campaign, "RK4-2049", alternative_plan
        )
        outcome, predecessor = _synthetic_reject(campaign)
        if not isinstance(outcome, TDG6TemporalRetryRequired):
            raise AssertionError("synthetic retry-evidence case exhausted")
        _synthetic_close_retry(campaign, "RK4-2049", outcome, predecessor)
        rejection = _synthetic_record_at(root, 1)
        rejection["payload"]["TDG6_rejection_evidence"] = (
            alternative.evidence.as_mapping()
        )
        rejection["payload"]["predecessor_plan"] = _plan_mapping(alternative_plan)
        rejection = _rewrite_synthetic_record(root, rejection)
        transition = _synthetic_record_at(root, 2)
        transition["journal_parent_sha256"] = rejection["record_sha256"]
        transition["payload"]["durable_rejection_record_sha256"] = (
            rejection["record_sha256"]
        )
        transition = _rewrite_synthetic_record(root, transition)
        checkpoint = _synthetic_checkpoint_at(root, 1)
        checkpoint["journal_tip_sha256"] = transition["record_sha256"]
        _rewrite_synthetic_checkpoint(root, checkpoint)

    def accepted_forgery(
        root: Path, campaign: Proto15SyntheticCampaign, *, ledger: bool,
    ) -> None:
        _synthetic_accept(campaign)
        checkpoint = _synthetic_checkpoint_at(root, 1)
        accepted = _synthetic_record_at(root, 1)
        cursor = dict(checkpoint["cursors"]["RK4-2049"])
        if ledger:
            checkpoint["ledgers"]["RK4-2049"]["accepted_macro_step_count"] = 2
            cursor["TDG6_ledger_sha256"] = _hash(
                checkpoint["ledgers"]["RK4-2049"]
            )
            accepted["payload"]["TDG6_ledger_sha256"] = cursor[
                "TDG6_ledger_sha256"
            ]
        else:
            cursor["previous_step_index"] = 1
        cursor = dict(_cursor(cursor).payload)
        checkpoint["cursors"]["RK4-2049"] = cursor
        accepted["payload"]["cursor"] = cursor
        accepted["payload"]["successor_cursor_chain_sha256"] = cursor[
            "cursor_chain_sha256"
        ]
        accepted = _rewrite_synthetic_record(root, accepted)
        checkpoint["journal_tip_sha256"] = accepted["record_sha256"]
        _rewrite_synthetic_checkpoint(root, checkpoint)

    def accepted_time_forgery(
        root: Path, campaign: Proto15SyntheticCampaign, *, target_overrun: bool,
    ) -> None:
        """Rehash an accepted record with an impossible time relation.

        The exact plan must reject either relation before a forged accepted
        cursor can become a new durable FRESH_READY state.  Rehashing the
        cursor, record, and checkpoint makes this a persisted-graph attack,
        not merely an entrypoint misuse.
        """
        _synthetic_accept(campaign)
        checkpoint = _synthetic_checkpoint_at(root, 1)
        accepted = _synthetic_record_at(root, 1)
        cursor = dict(checkpoint["cursors"]["RK4-2049"])
        # Obtain the original boundary from the immutable predecessor cursor,
        # rather than selecting a new timing value from the live accepted state.
        previous = _synthetic_checkpoint_at(root, 0)
        old_identity = previous["cursors"]["RK4-2049"]["accepted_boundary_time"]
        target_identity = cursor["event_target_time"]
        forged_identity = target_identity if target_overrun else old_identity
        if target_overrun:
            target = float.fromhex(str(target_identity["binary64_hex"]))
            forged_identity = _time_identity(np.nextafter(target, float("inf")))
        cursor["accepted_boundary_time"] = forged_identity
        checkpoint["ledgers"]["RK4-2049"]["last_accepted_time_hex"] = (
            forged_identity["binary64_hex"]
        )
        cursor["TDG6_ledger_sha256"] = _hash(checkpoint["ledgers"]["RK4-2049"])
        cursor = dict(_cursor(cursor).payload)
        checkpoint["cursors"]["RK4-2049"] = cursor
        accepted["payload"]["cursor"] = cursor
        accepted["payload"]["successor_cursor_chain_sha256"] = cursor[
            "cursor_chain_sha256"
        ]
        accepted["payload"]["TDG6_ledger_sha256"] = cursor["TDG6_ledger_sha256"]
        accepted = _rewrite_synthetic_record(root, accepted)
        checkpoint["journal_tip_sha256"] = accepted["record_sha256"]
        _rewrite_synthetic_checkpoint(root, checkpoint)

    def terminal_reason_forgery(
        root: Path, campaign: Proto15SyntheticCampaign,
    ) -> None:
        outcome, _predecessor = _synthetic_reject(
            campaign, requested_cap=float(TDG6_MINIMUM_MACRO_STEP)
        )
        if not isinstance(outcome, TDG6TemporalRetryExhausted):
            raise AssertionError("synthetic terminal mutation did not exhaust")
        campaign.close_tdg6_outcome("RK4-2049", outcome)
        terminal = _synthetic_record_at(root, 2)
        terminal["payload"]["exhaustion_reason"] = "maximum_temporal_retries"
        terminal = _rewrite_synthetic_record(root, terminal)
        checkpoint = _synthetic_checkpoint_at(root, 1)
        checkpoint["terminal_closure"]["exhaustion_reason"] = (
            "maximum_temporal_retries"
        )
        checkpoint["journal_tip_sha256"] = terminal["record_sha256"]
        checkpoint["terminal_record_sha256"] = terminal["record_sha256"]
        _rewrite_synthetic_checkpoint(root, checkpoint)

    def terminal_older_rejection_pointer_forgery(
        root: Path, campaign: Proto15SyntheticCampaign,
    ) -> None:
        """Point a terminal closure at an older valid rejection after rehashing."""
        first, first_plan = _synthetic_reject(campaign)
        if not isinstance(first, TDG6TemporalRetryRequired):
            raise AssertionError("first terminal-pointer rejection exhausted")
        _synthetic_close_retry(campaign, "RK4-2049", first, first_plan)
        terminal: TDG6TemporalRetryExhausted | None = None
        for _ in range(TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP + 2):
            pending = campaign.cursors["RK4-2049"].payload[
                "retry_successor_payload_or_none"
            ]
            if not isinstance(pending, Mapping):
                raise AssertionError("terminal-pointer retry cursor is not pending")
            predecessor = plan_forward_proto14_subdivision(
                float.fromhex(str(campaign.cursors["RK4-2049"].payload[
                    "accepted_boundary_time"
                ]["binary64_hex"])),
                float.fromhex(str(campaign.cursors["RK4-2049"].payload[
                    "event_target_time"
                ]["binary64_hex"])),
                float.fromhex(str(pending["successor_plan"]["requested_cap_hex"])),
                minimum_width=float(TDG6_MINIMUM_MACRO_STEP),
            )
            if _plan_mapping(predecessor) != pending["successor_plan"]:
                raise AssertionError("terminal-pointer retry plan differs from cursor")
            outcome = _synthetic_retry_outcome(campaign, "RK4-2049", predecessor)
            campaign.durable_rejection_sink(
                "RK4-2049", predecessor_plan=predecessor
            )(outcome.evidence.as_mapping())
            if isinstance(outcome, TDG6TemporalRetryExhausted):
                campaign.close_tdg6_outcome("RK4-2049", outcome)
                terminal = outcome
                break
            _synthetic_close_retry(campaign, "RK4-2049", outcome, predecessor)
        if terminal is None:
            raise AssertionError("terminal-pointer retry sequence did not exhaust")
        first_rejection = _synthetic_record_at(root, 1)
        terminal = Proto15SyntheticCampaign._read_record(
            sorted((root / "journal").glob("*.journal"))[-1]
        )
        terminal = dict(terminal)
        terminal["payload"]["durable_rejection_record_sha256"] = first_rejection[
            "record_sha256"
        ]
        terminal["payload"]["journal_tip_sha256"] = first_rejection["record_sha256"]
        terminal = _rewrite_synthetic_record(root, terminal)
        checkpoint = _synthetic_checkpoint_at(root, campaign.generation)
        checkpoint["terminal_closure"]["durable_rejection_record_sha256"] = (
            first_rejection["record_sha256"]
        )
        checkpoint["terminal_closure"]["journal_tip_sha256"] = first_rejection[
            "record_sha256"
        ]
        checkpoint["journal_tip_sha256"] = terminal["record_sha256"]
        checkpoint["terminal_record_sha256"] = terminal["record_sha256"]
        _rewrite_synthetic_checkpoint(root, checkpoint)

    def rollback_forgery(
        root: Path, campaign: Proto15SyntheticCampaign, *, ancestor: bool,
    ) -> None:
        snapshot = campaign.common_event_snapshot()
        _synthetic_accept(campaign)
        previous = _synthetic_checkpoint_at(root, 1)
        aborted = str(campaign.cursors["RK4-2049"].payload["attempt_id"])
        campaign.common_event_rollback(snapshot, aborted_attempt_id=aborted)
        rollback = _synthetic_record_at(root, 2)
        if ancestor:
            rollback["payload"]["restored_snapshot_checkpoint_sha256"] = (
                previous["checkpoint_sha256"]
            )
            rollback["payload"]["restored_snapshot_campaign_generation"] = 1
            rollback["payload"]["restored_snapshot_cursor_set_sha256"] = previous[
                "cursor_set_sha256"
            ]
            rollback["payload"]["restored_snapshot_ledger_set_sha256"] = previous[
                "ledger_set_sha256"
            ]
            rollback["payload"]["restored_snapshot_accepted_state_set_sha256"] = (
                previous["accepted_state_set_sha256"]
            )
        else:
            rollback["payload"]["new_cursor_set_sha256"] = "f" * 64
        rollback = _rewrite_synthetic_record(root, rollback)
        checkpoint = _synthetic_checkpoint_at(root, 2)
        checkpoint["journal_tip_sha256"] = rollback["record_sha256"]
        _rewrite_synthetic_checkpoint(root, checkpoint)

    def rollback_cursor_content_forgery(
        root: Path, campaign: Proto15SyntheticCampaign,
    ) -> None:
        """Rehash a rollback with one restored cursor field replaced."""
        snapshot = campaign.common_event_snapshot()
        _synthetic_accept(campaign)
        aborted = str(campaign.cursors["RK4-2049"].payload["attempt_id"])
        campaign.common_event_rollback(snapshot, aborted_attempt_id=aborted)
        checkpoint = _synthetic_checkpoint_at(root, 2)
        cursor = dict(checkpoint["cursors"]["RK4-2049"])
        cursor["previous_transaction_serial"] = (
            int(cursor["previous_transaction_serial"]) + 1
        )
        cursor = dict(_cursor(cursor).payload)
        checkpoint["cursors"]["RK4-2049"] = cursor
        _rewrite_synthetic_checkpoint(root, checkpoint)
        rollback = _synthetic_record_at(root, 2)
        rollback["payload"]["new_cursor_set_sha256"] = checkpoint[
            "cursor_set_sha256"
        ]
        rollback = _rewrite_synthetic_record(root, rollback)
        checkpoint["journal_tip_sha256"] = rollback["record_sha256"]
        _rewrite_synthetic_checkpoint(root, checkpoint)

    invalid("missing_state", remove_state)
    invalid("mutated_state", corrupt_state)
    invalid("same_generation_fork", fork_checkpoint)
    invalid("pending_plan_substitution", pending_plan_substitution)
    invalid("retry_evidence_fork", retry_evidence_fork)
    invalid("accepted_serial_forgery", lambda root, campaign: accepted_forgery(root, campaign, ledger=False))
    invalid("accepted_ledger_forgery", lambda root, campaign: accepted_forgery(root, campaign, ledger=True))
    invalid("accepted_no_advance_forgery", lambda root, campaign: accepted_time_forgery(root, campaign, target_overrun=False))
    invalid("accepted_target_overrun_forgery", lambda root, campaign: accepted_time_forgery(root, campaign, target_overrun=True))
    invalid("terminal_reason_forgery", terminal_reason_forgery)
    invalid("terminal_older_rejection_pointer_forgery", terminal_older_rejection_pointer_forgery)
    invalid("rollback_ancestor_forgery", lambda root, campaign: rollback_forgery(root, campaign, ancestor=True))
    invalid("rollback_set_forgery", lambda root, campaign: rollback_forgery(root, campaign, ancestor=False))
    invalid("rollback_cursor_content_forgery", rollback_cursor_content_forgery)
    return records


def synthetic_qualification() -> dict[str, object]:
    """Execute and summarize HLT13's deterministic synthetic qualification corpus."""
    scenarios = (
        "state_object", "accepted_journal", "accepted_checkpoint",
        "rejection_journal", "cursor_transition", "terminal_record",
        "terminal_checkpoint", "rollback_record", "rollback_checkpoint",
    )
    nominal = _run_nominal_cases()
    mutation = _run_mutation_cases()
    fault = [
        _run_fault_case(scenario, phase)
        for scenario in scenarios
        for phase in FAULT_PHASES
    ]
    cases = [*nominal, *mutation, *fault]
    failed = [str(case["case"]) for case in cases if case["passed"] is not True]
    if failed:
        raise AssertionError(f"HLT13 synthetic qualification failed: {failed}")
    phase_coverage = sorted({str(case["phase"]) for case in fault})
    target_coverage = sorted({str(case["target"]) for case in fault})
    ambiguous_crash_cases = [
        case for case in fault if case.get("crash_images_tested") == ["old", "new"]
    ]
    return {
        "artifact_id": "FGC-1-HLT13-MON13",
        "runtime_module": "src/recursive_horizons/fgc/evolution/proto15_runtime.py",
        "synthetic_only": True,
        "six_member_contract": len(_synthetic_members()) == len(MEMBER_KEYS),
        "six_member_order": list(MEMBER_KEYS),
        "cursor_modes": ["FRESH_READY", "RETRY_PENDING"],
        "journal_record_kinds": sorted(RECORD_KINDS),
        "canonical_encoding": "utf8_sorted_compact_json_ensure_ascii_true_allow_nan_false",
        "fault_hook_windows": list(FAULT_PHASES),
        "qualification_case_count": len(cases),
        "nominal_case_count": len(nominal),
        "mutation_case_count": len(mutation),
        "fault_case_count": len(fault),
        "fault_scenario_count": len(scenarios),
        "fault_scenario_coverage": list(scenarios),
        "fault_phase_coverage": phase_coverage,
        "fault_target_coverage": target_coverage,
        "ambiguous_replace_crash_case_count": len(ambiguous_crash_cases),
        "fault_recovery_image_count": len(fault) + len(ambiguous_crash_cases),
        "qualification_case_digest": _hash(cases),
        "all_qualification_cases_passed": not failed,
        "content_addressed_state_objects": any(
            case["case"] == "mutation:mutated_state" and case["passed"] for case in mutation
        ),
        "immutable_one_record_journal_framing": all(
            case["passed"] for case in fault if case["target"] == "journal_record"
        ),
        "rational_and_hex_time_identity": nominal[0]["passed"],
        "high_water_marks_persisted": nominal[-1]["passed"],
        "checkpoint_fork_rejected": any(
            case["case"] == "mutation:same_generation_fork" and case["passed"] for case in mutation
        ),
        "terminal_cursor_set_preserved": any(
            case["case"] == "nominal:terminal" and case["passed"] for case in nominal
        ),
        "terminal_ledger_embedded": any(
            case["case"] == "nominal:terminal" and case["passed"] for case in nominal
        ),
        "exact_rejection_suffix_recovery_only": any(
            case["case"] == "nominal:rejection_suffix" and case["passed"] for case in nominal
        ),
        "uncheckpointed_accepted_or_rollback_suffix_invalid": all(
            case["passed"] for case in fault
            if case["case"].startswith("fault:accepted_")
            or case["case"].startswith("fault:rollback_")
        ),
        "six_member_fine_checkpoint": any(
            case["case"] == "nominal:accepted_fine" and case["passed"] for case in nominal
        ),
        "accepted_pending_retry_plan_bound": any(
            case["case"] == "nominal:retry_pending" and case["passed"] for case in nominal
        ),
        "rollback_new_identities": any(
            case["case"] == "nominal:rollback" and case["passed"] for case in nominal
        ),
        "rollback_protocol_journal_field_locations": dict(
            ROLLBACK_PROTOCOL_JOURNAL_FIELD_LOCATIONS
        ),
        "rollback_journal_hashes_bound_by_record_envelope": any(
            case["case"] == "nominal:rollback" and case["passed"] for case in nominal
        ),
        "all_fault_hooks_exercised": phase_coverage == sorted(FAULT_PHASES),
        "old_or_new_replace_crash_images_exercised": (
            len(ambiguous_crash_cases)
            == len(scenarios) * len(_AMBIGUOUS_REPLACE_CRASH_PHASES)
        ),
        "hash_integrity_requires_external_trusted_genesis": True,
        "wholesale_campaign_replacement_authenticated": False,
        "production_runner_implemented": False,
        "output_namespace_created": False,
        "trajectory_authorized": False,
        "fresh_GR0_calibration_authorized": False,
        "physical_classification": False,
    }


__all__ = ["MEMBER_KEYS", "ROLLBACK_PROTOCOL_JOURNAL_FIELD_LOCATIONS", "Proto15CheckpointIntegrityError", "Proto15Cursor", "Proto15CursorContractError", "Proto15DurabilityError", "Proto15JournalIntegrityError", "Proto15RecoveryInvalidStop", "Proto15SyntheticCampaign", "Proto15TerminalLocked", "recover_campaign", "synthetic_qualification"]
