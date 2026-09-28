"""One pure-in-process GR-0 attempt for the authenticated HLT16 campaign.

This adapter makes exactly one TDG7/TDG6 proposal.  It has no loop, no file
I/O, no checkpoint ownership, and no candidate-branch capability.  Its caller
either persists the returned evidence or discards the in-memory attempt.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass
import json
import math
from typing import Any, Callable, Literal, Mapping

from . import hlt16_lifecycle as lifecycle
from . import proto15_runtime as p15
from . import tdg6_temporal_admission_runtime as tdg6
from . import tdg7_stage_safe_runtime as tdg7
from . import tdg8_persisted_retry_replay as tdg8_replay
from .proto14_runtime import (
    Proto14InvalidTemporalRuntime,
    Proto14RunMember,
    Proto14SourceRetryExhausted,
)
from .proto5_runtime import CFLRetryRequired, GR0RuntimeStop
from .proto7_runtime import Proto7TerminalStop
from .numerical_engine import array_content_sha256


Disposition = Literal[
    "accepted_fine",
    "source_retry_required",
    "source_retry_exhausted",
    "cfl_retry_required",
    "cfl_retry_exhausted",
    "temporal_retry_required",
    "temporal_retry_exhausted",
    "scientific_stop",
    "invalid",
]
DISPOSITIONS: frozenset[str] = frozenset(Disposition.__args__)


@dataclass(frozen=True, slots=True)
class RetryCounters:
    """Campaign-owned retry counters supplied without mutating durable state."""

    source_current: int = 0
    source_total: int = 0
    cfl_current: int = 0
    cfl_total: int = 0
    source_cap: int = 32
    cfl_cap: int = 32

    def __post_init__(self) -> None:
        values = (
            self.source_current, self.source_total, self.cfl_current,
            self.cfl_total, self.source_cap, self.cfl_cap,
        )
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in values):
            raise ValueError("retry counters must be nonnegative integers")
        if self.source_current > self.source_total or self.cfl_current > self.cfl_total:
            raise ValueError("retry current count exceeds its total")
        if self.source_cap != 32 or self.cfl_cap != 32:
            raise ValueError("retry caps differ from the frozen HLT16 contract")
        if self.source_current > self.source_cap + 1 or self.cfl_current > self.cfl_cap + 1:
            raise ValueError("retry current count exceeds the terminal exhaustion edge")


@dataclass(frozen=True, slots=True)
class MemberAttemptOutcome:
    """Immutable evidence for exactly one proposal; never a durable record."""

    disposition: Disposition
    executed_plan: object | None
    accepted_time: float | None
    accepted_state: object | None
    accepted_step_index: int | None
    accepted_transaction_serial: int | None
    accepted_ledger: tdg6.TDG6TemporalLedger | None
    retry_outcome: object | None
    retry_counters: RetryCounters
    next_cap: float | None
    evidence: Mapping[str, object]
    accepted_boundary_restored: bool

    def __post_init__(self) -> None:
        if self.disposition not in DISPOSITIONS:
            raise ValueError("unknown HLT16 attempt disposition")
        if not isinstance(self.retry_counters, RetryCounters):
            raise TypeError("retry counters are absent")
        if not isinstance(self.evidence, Mapping):
            raise TypeError("attempt evidence is absent")
        try:
            json.dumps(
                self.evidence,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
                allow_nan=False,
            )
        except (TypeError, ValueError) as error:
            raise ValueError("attempt evidence is not canonical-JSON compatible") from error
        if self.disposition == "accepted_fine":
            if (
                self.executed_plan is None or self.accepted_time is None
                or self.accepted_state is None or self.accepted_step_index is None
                or self.accepted_transaction_serial is None or self.accepted_ledger is None
                or self.retry_outcome is not None or self.next_cap is not None
                or self.accepted_boundary_restored
            ):
                raise ValueError("accepted attempt evidence is incomplete")
        elif self.accepted_boundary_restored is not True:
            raise ValueError("nonaccepted attempt must restore the accepted boundary")


class HLT16AttemptError(ValueError):
    """The in-memory input lineage is malformed before a proposal is trusted."""


def _same_binary64(left: float, right: float) -> bool:
    return float(left).hex() == float(right).hex()


def _positive(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise HLT16AttemptError(f"{label} is not a finite positive scalar")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise HLT16AttemptError(f"{label} is not a finite positive scalar")
    return result


def _jsonable(value: object) -> object:
    """Normalize retained evidence without preserving live exception objects."""
    if is_dataclass(value) and not isinstance(value, type):
        value = asdict(value)
    if isinstance(value, Mapping):
        result: dict[str, object] = {}
        for key, item in value.items():
            normalized = str(key)
            if normalized in result:
                raise HLT16AttemptError("attempt evidence keys collide after normalization")
            result[normalized] = _jsonable(item)
        return result
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        if isinstance(value, float) and not math.isfinite(value):
            raise HLT16AttemptError("attempt evidence contains a nonfinite scalar")
        return value
    raise HLT16AttemptError(
        f"attempt evidence contains unsupported {type(value).__name__}"
    )


def _safe_restore(member: Proto14RunMember, boundary: Mapping[str, Any]) -> bool:
    """Restore every nonaccepted route and verify the inherited boundary hook."""
    member.restore(boundary)
    checker = getattr(member, "_accepted_boundary_preserved", None)
    if callable(checker) and checker(boundary) is not True:
        raise HLT16AttemptError("nonaccepted attempt did not restore the accepted boundary")
    return True


def _cursor_and_member_preflight(
    member: Proto14RunMember,
    cursor: p15.Proto15Cursor,
    target: float,
    descriptor_sha256: str,
) -> None:
    if not isinstance(member, Proto14RunMember) or not isinstance(cursor, p15.Proto15Cursor):
        raise HLT16AttemptError("attempt requires a PROTO14 member and PROTO18 cursor")
    cursor.validate()
    if cursor.payload["protocol_artifact_id"] != "FGC-2-SF1-PROTO18":
        raise HLT16AttemptError("attempt cursor is not PROTO18")
    if cursor.payload["member_key"] != member.key or cursor.payload["method"] != member.method_label:
        raise HLT16AttemptError("cursor and GR-0 member identity differs")
    if cursor.payload["accepted_state_sha256"] != descriptor_sha256:
        raise HLT16AttemptError("restored descriptor differs from cursor state identity")
    if not _same_binary64(member.time, float.fromhex(cursor.payload["accepted_boundary_time"]["binary64_hex"])):
        raise HLT16AttemptError("member time differs from cursor accepted boundary")
    if not _same_binary64(target, float.fromhex(cursor.payload["event_target_time"]["binary64_hex"])):
        raise HLT16AttemptError("requested target differs from cursor event target")
    if member.temporal_ledger is None:
        raise HLT16AttemptError("member omits TDG6 ledger")
    if cursor.payload["TDG6_ledger_sha256"] != p15._hash(p15._ledger_mapping(member.temporal_ledger)):
        raise HLT16AttemptError("member TDG6 ledger differs from cursor")


def _prepare_initial(member: Proto14RunMember, target: float, ledger: tdg6.TDG6TemporalLedger, cap: float) -> object:
    return tdg7.prepare_tdg7_initial_stage_safe_runtime(
        method=member.integrator_id,
        time=member.time,
        event_target=target,
        requested_cap=cap,
        state=member.state,
        rhs=member.operator,
        projector=member.projector,
        transaction=member.transaction,
        tracers=member.tracers,
        coordinates=member.initial.grid.coordinates,
        temporal_ledger=ledger,
        previous_step_index=member.step_index,
        previous_transaction_serial=member.transaction_serial,
    )


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Decode persisted evidence without accepting duplicate object keys."""
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise HLT16AttemptError("pending TDG6 rejection prefix has duplicate keys")
        result[key] = value
    return result


def _active_rejection_chain(
    cursor: p15.Proto15Cursor,
    ledger: tdg6.TDG6TemporalLedger,
    *,
    physical_state_sha256: str,
) -> tuple[tdg6.TDG6TemporalLedger, tuple[dict[str, object], ...]]:
    """Return the latest predecessor ledger and exact active rejection suffix.

    ``serialized_temporal_rejections`` is cumulative across accepted macro
    steps.  Only its final ``current_macro_step_temporal_retry_count`` records
    belong to the retry currently represented by the cursor.  Roll that suffix
    back exactly, validate it forward, and cross-bind every adjacent TDG7 plan.
    Older PDE proposals are deliberately not rerun because their full plans
    are not retained; only the latest exact persisted predecessor is replayed.
    """
    pending = cursor.payload.get("retry_successor_payload_or_none")
    if not isinstance(pending, Mapping):
        raise HLT16AttemptError("pending cursor omits its replay payload")
    prefix = pending.get("complete_rejection_prefix")
    retry_count = ledger.current_macro_step_temporal_retry_count
    if (
        not isinstance(prefix, list)
        or tuple(prefix) != ledger.serialized_temporal_rejections
        or retry_count < 1
        or retry_count > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
        or len(prefix) < retry_count
        or pending.get("retry_count") != retry_count
    ):
        raise HLT16AttemptError("pending TDG6 rejection prefix differs")

    decoded: list[dict[str, object]] = []
    for serialized in prefix:
        if not isinstance(serialized, str):
            raise HLT16AttemptError("pending TDG6 rejection prefix is not text")
        try:
            value = json.loads(serialized, object_pairs_hook=_unique_json_object)
        except (TypeError, ValueError, json.JSONDecodeError) as error:
            if isinstance(error, HLT16AttemptError):
                raise
            raise HLT16AttemptError(
                "pending TDG6 rejection prefix is not canonical JSON"
            ) from error
        if (
            not isinstance(value, Mapping)
            or p15._canonical(dict(value)).decode("utf-8") != serialized
        ):
            raise HLT16AttemptError(
                "pending TDG6 rejection prefix is not canonical JSON"
            )
        decoded.append(dict(value))

    active = tuple(decoded[-retry_count:])
    durable_tail = pending.get("TDG6_rejection_evidence")
    if (
        not isinstance(durable_tail, Mapping)
        or active[-1] != dict(durable_tail)
        or p15._canonical(active[-1]) != p15._canonical(dict(durable_tail))
    ):
        raise HLT16AttemptError("pending TDG6 rejection tail differs")

    base = ledger
    active_serialized = tuple(prefix[-retry_count:])
    try:
        for serialized, evidence in reversed(tuple(zip(active_serialized, active))):
            if (
                not base.serialized_temporal_rejections
                or base.serialized_temporal_rejections[-1] != serialized
                or evidence.get("retry_count_for_current_macro_step")
                != base.current_macro_step_temporal_retry_count
                or evidence.get("cumulative_temporal_retry_count")
                != base.cumulative_temporal_retry_count
            ):
                raise HLT16AttemptError(
                    "pending TDG6 rejection prefix is noncontiguous"
                )
            base = p15._prior_retry_ledger(base, evidence)
    except (TypeError, ValueError, KeyError, p15.Proto15CursorContractError) as error:
        if isinstance(error, HLT16AttemptError):
            raise
        raise HLT16AttemptError(
            "pending TDG6 rejection prefix is noncontiguous"
        ) from error
    if base.current_macro_step_temporal_retry_count != 0:
        raise HLT16AttemptError("pending TDG6 rejection prefix is noncontiguous")

    prior = base
    latest_prior = base
    checked_active: list[dict[str, object]] = []
    try:
        for index, evidence in enumerate(active):
            checked = dict(
                p15._validate_retry_evidence_mapping(
                    evidence,
                    method=str(cursor.payload["method"]),
                    state_sha256=physical_state_sha256,
                    accepted_time=cursor.payload["accepted_boundary_time"],
                    step_index=int(cursor.payload["previous_step_index"]),
                    transaction_serial=int(
                        cursor.payload["previous_transaction_serial"]
                    ),
                    prior_ledger=prior,
                )
            )
            if (
                checked != evidence
                or p15._canonical(checked) != p15._canonical(evidence)
                or int(checked["retry_count_for_current_macro_step"])
                > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
                or float(checked["retry_step_size"])
                < float(p15.TDG6_MINIMUM_MACRO_STEP)
            ):
                raise HLT16AttemptError(
                    "pending TDG6 rejection prefix is not retryable"
                )
            if index:
                previous = checked_active[-1]
                expected_plan = p15._plan_mapping(
                    p15.plan_forward_proto14_subdivision(
                        float(previous["initial_time"]),
                        float.fromhex(
                            str(cursor.payload["event_target_time"]["binary64_hex"])
                        ),
                        float(previous["retry_step_size"]),
                        minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
                    )
                )
                if not _same_binary64(
                    float.fromhex(expected_plan["macro_width_hex"]),
                    float(checked["attempted_macro_step_size"]),
                ):
                    raise HLT16AttemptError(
                        "pending TDG6 rejection prefix has a broken plan adjacency"
                    )
                if index == len(active) - 1 and expected_plan != p15._plan_mapping(
                    pending["predecessor_plan"]
                ):
                    raise HLT16AttemptError(
                        "pending TDG6 rejection prefix has a broken plan adjacency"
                    )
            if index == len(active) - 1:
                latest_prior = prior
            prior = p15._updated_retry_ledger(prior, checked)
            checked_active.append(checked)
    except (TypeError, ValueError, KeyError, p15.Proto15CursorContractError) as error:
        if isinstance(error, HLT16AttemptError):
            raise
        raise HLT16AttemptError(
            "pending TDG6 rejection prefix is noncontiguous"
        ) from error
    if prior != ledger:
        raise HLT16AttemptError("pending TDG6 rejection prefix is noncontiguous")
    return latest_prior, tuple(checked_active)


def _replay_pending_rejection(
    member: Proto14RunMember,
    cursor: p15.Proto15Cursor,
    target: float,
) -> object:
    """Recompute the exact latest rejected predecessor before retry prepare.

    The complete active evidence suffix is authenticated, rolled back, and
    cross-bound first.  Only the latest proposal is shadow-rerun: persistence
    retains its complete predecessor plan, but not the full plans of older
    rejected proposals.  At retry depth one the original fresh path is
    unchanged; deeper histories use the disjoint TDG8 persisted-replay overlay.
    """
    ledger = member.temporal_ledger
    assert ledger is not None
    physical_state_sha256 = array_content_sha256(
        member.state.u, member.state.p, member.state.q
    )
    try:
        lifecycle.validate_retry_payload(
            cursor,
            ledger,
            evolution_state_sha256=physical_state_sha256,
        )
    except (TypeError, ValueError, KeyError) as error:
        raise HLT16AttemptError("pending cursor replay payload differs") from error
    pending = cursor.payload["retry_successor_payload_or_none"]
    if not isinstance(pending, Mapping):
        raise HLT16AttemptError("pending cursor omits its replay payload")
    prior, evidence_chain = _active_rejection_chain(
        cursor,
        ledger,
        physical_state_sha256=physical_state_sha256,
    )
    evidence = evidence_chain[-1]
    predecessor = p15._plan_mapping(pending["predecessor_plan"])
    successor_plan = p15._plan_mapping(pending["successor_plan"])
    try:
        replay_plan = p15.plan_forward_proto14_subdivision(
            float.fromhex(predecessor["current_hex"]),
            float.fromhex(predecessor["event_target_hex"]),
            float.fromhex(predecessor["requested_cap_hex"]),
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        )
        if p15._plan_mapping(replay_plan) != predecessor:
            raise HLT16AttemptError(
                "pending TDG6 predecessor plan cannot be reconstructed"
            )
        if len(evidence_chain) == 1:
            prepared = _prepare_initial(
                member,
                target,
                prior,
                float.fromhex(predecessor["requested_cap_hex"]),
            )
        else:
            prepared = tdg8_replay.prepare_tdg8_persisted_retry_replay(
                replay_plan=replay_plan,
                expected_prior_retry_count=(
                    prior.current_macro_step_temporal_retry_count
                ),
                method=member.integrator_id,
                state=member.state,
                rhs=member.operator,
                projector=member.projector,
                transaction=member.transaction,
                tracers=member.tracers,
                coordinates=member.initial.grid.coordinates,
                temporal_ledger=prior,
                previous_step_index=member.step_index,
                previous_transaction_serial=member.transaction_serial,
            )
    except (tdg6.TDG6RefinementPathStop, TypeError, ValueError) as error:
        raise HLT16AttemptError("pending TDG6 predecessor no longer replays") from error
    if (
        p15._plan_mapping(prepared.plan) != predecessor
        or not _same_binary64(
            float.fromhex(predecessor["macro_width_hex"]),
            float(evidence["attempted_macro_step_size"]),
        )
    ):
        raise HLT16AttemptError(
            "pending TDG6 predecessor plan differs from its durable chain"
        )

    seen: list[Mapping[str, object]] = []
    require_admission = (
        tdg7.require_tdg7_stage_safe_admission
        if len(evidence_chain) == 1
        else tdg8_replay.require_tdg8_persisted_retry_replay_admission
    )
    replan_after_rejection = (
        tdg7.replan_tdg7_after_tdg6_retry
        if len(evidence_chain) == 1
        else tdg8_replay.replan_tdg8_after_persisted_retry_rejection
    )
    try:
        require_admission(
            prepared,
            transaction=member.transaction,
            tracers=member.tracers,
            temporal_ledger=prior,
            current_time=member.time,
            current_state=member.state,
            current_step_index=member.step_index,
            current_transaction_serial=member.transaction_serial,
            durable_rejection_sink=lambda item: seen.append(dict(item)),
        )
    except tdg6.TDG6TemporalRetryRequired as replayed:
        if len(seen) != 1:
            raise HLT16AttemptError(
                "pending TDG6 replay evidence differs from the durable cursor"
            )
        normalized_seen = p15._json_safe(seen[0])
        if (
            normalized_seen != evidence
            or p15._canonical(normalized_seen) != p15._canonical(evidence)
        ):
            # Python considers ``1`` and ``1.0`` equal.  Canonical persistence
            # bytes retain that distinction, so tuple/list normalization cannot
            # weaken scalar or ordered-array identity.
            raise HLT16AttemptError(
                "pending TDG6 replay evidence differs from the durable cursor"
            )
        expected_ledger = p15._updated_retry_ledger(prior, evidence)
        if replayed.updated_ledger != expected_ledger or expected_ledger != ledger:
            raise HLT16AttemptError(
                "pending TDG6 replay ledger differs from the durable cursor"
            )
        try:
            successor = replan_after_rejection(prepared, replayed)
        except (TypeError, ValueError) as error:
            raise HLT16AttemptError(
                "pending TDG7 replay successor cannot be reconstructed"
            ) from error
        if p15._plan_mapping(successor.plan) != successor_plan:
            raise HLT16AttemptError(
                "pending TDG7 replay successor differs from the durable cursor"
            )
        return successor
    except tdg6.TDG6TemporalRetryExhausted as error:
        raise HLT16AttemptError("pending cursor replays as temporal exhaustion") from error
    raise HLT16AttemptError("pending cursor replays as an admitted predecessor")


def _scientific(error: BaseException) -> bool:
    return isinstance(error, (GR0RuntimeStop, Proto7TerminalStop))


def _source_evidence(stop: tdg6.TDG6RefinementPathStop, plan: object) -> dict[str, object]:
    if stop.source_retry is None:
        raise HLT16AttemptError("source retry omits proposal evidence")
    return {
        "kind": "source_retry", "path": stop.path,
        "source_retry": _jsonable(stop.source_retry),
        "attempted_cap_hex": _plan_width(plan).hex(),
        "accepted_boundary_preserved": True, "classification_is_numerical_not_physical": True,
    }


def _cfl_evidence(stop: tdg6.TDG6RefinementPathStop, plan: object) -> dict[str, object]:
    cause = stop.cause
    assert isinstance(cause, CFLRetryRequired)
    return {
        "kind": "cfl_retry", "path": stop.path, "observed_ratio": cause.observed_ratio,
        "maximum_ratio": cause.maximum_ratio, "attempted_cap_hex": _plan_width(plan).hex(),
        "accepted_boundary_preserved": True, "classification_is_numerical_not_physical": True,
    }


def _plan_width(plan: object) -> float:
    if isinstance(plan, Mapping):
        return float.fromhex(str(plan["macro_width_hex"]))
    width = getattr(plan, "macro_width", None)
    if (
        isinstance(width, (int, float))
        and not isinstance(width, bool)
        and math.isfinite(float(width))
        and width > 0.0
    ):
        return float(width)
    raise HLT16AttemptError("attempted TDG7 plan has no finite macro width")


def attempt_once(
    member: Proto14RunMember,
    cursor: p15.Proto15Cursor,
    target: float,
    *,
    descriptor_sha256: str,
    requested_cap: float,
    retry_counters: RetryCounters = RetryCounters(),
    temporal_rejection_sink: Callable[[Mapping[str, object]], None] | None = None,
) -> MemberAttemptOutcome:
    """Run one exact fresh or pending GR-0 proposal, with no retry loop.

    A temporal rejection invokes ``temporal_rejection_sink`` before it is
    returned.  The campaign writer supplies the durable implementation; this
    adapter does not alter a checkpoint or cursor itself.
    """
    boundary: Mapping[str, Any] | None = None
    planned: object | None = None
    try:
        _cursor_and_member_preflight(member, cursor, target, descriptor_sha256)
        requested = _positive(requested_cap, "requested cap")
        remaining = target - member.time
        if requested > remaining or not math.isfinite(remaining) or remaining <= 0.0:
            raise HLT16AttemptError("requested cap exceeds the positive event remainder")
        boundary = member.snapshot()
        ledger = member.temporal_ledger
        assert ledger is not None
        if cursor.mode == "FRESH_READY":
            planned = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
                member.time, target, requested,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            ))
            prepared = _prepare_initial(member, target, ledger, requested)
            planned = prepared.plan
        elif cursor.mode == "RETRY_PENDING":
            pending = cursor.payload.get("retry_successor_payload_or_none")
            if not isinstance(pending, Mapping) or not _same_binary64(
                requested, float.fromhex(str(pending.get("half_cap")))
            ):
                raise HLT16AttemptError(
                    "requested cap differs from the durable temporal successor"
                )
            successor = _replay_pending_rejection(member, cursor, target)
            planned = successor.plan
            prepared = tdg7.prepare_tdg7_retry_successor(
                successor=successor,
                method=member.integrator_id,
                state=member.state,
                rhs=member.operator,
                projector=member.projector,
                transaction=member.transaction,
                tracers=member.tracers,
                coordinates=member.initial.grid.coordinates,
                previous_step_index=member.step_index,
                previous_transaction_serial=member.transaction_serial,
            )
            planned = prepared.plan
        else:
            raise HLT16AttemptError("cursor mode is not executable")
        try:
            tdg7.require_tdg7_stage_safe_admission(
                prepared,
                transaction=member.transaction,
                tracers=member.tracers,
                temporal_ledger=ledger,
                current_time=member.time,
                current_state=member.state,
                current_step_index=member.step_index,
                current_transaction_serial=member.transaction_serial,
                durable_rejection_sink=(temporal_rejection_sink or (lambda _: None)),
            )
        except tdg6.TDG6TemporalRetryRequired as retry:
            restored = _safe_restore(member, boundary)
            return MemberAttemptOutcome(
                "temporal_retry_required", prepared.plan, None, None, None, None, None,
                retry, retry_counters, retry.retry_step_size,
                {"kind": "TDG6_temporal_retry", "evidence": _jsonable(retry.evidence.as_mapping())}, restored,
            )
        except tdg6.TDG6TemporalRetryExhausted as retry:
            restored = _safe_restore(member, boundary)
            return MemberAttemptOutcome(
                "temporal_retry_exhausted", prepared.plan, None, None, None, None, None,
                retry, retry_counters, None,
                {"kind": "TDG6_temporal_exhaustion", "reason": retry.reason, "evidence": _jsonable(retry.evidence.as_mapping())}, restored,
            )
        committed = tdg7.commit_tdg7_stage_safe_runtime(
            prepared,
            transaction=member.transaction,
            tracers=member.tracers,
            temporal_ledger=ledger,
            current_time=member.time,
            current_state=member.state,
            current_step_index=member.step_index,
            current_transaction_serial=member.transaction_serial,
        )
        member.state = committed.state
        member.time = committed.time
        member.step_index = committed.step_index
        member.transaction_serial = committed.transaction_serial
        member.temporal_ledger = committed.temporal_ledger
        return MemberAttemptOutcome(
            "accepted_fine", prepared.plan, committed.time, committed.state,
            committed.step_index, committed.transaction_serial, committed.temporal_ledger,
            None, retry_counters, None,
            {"kind": "accepted_fine", "continuous_admission": _jsonable(committed.continuous_admission)}, False,
        )
    except tdg6.TDG6RefinementPathStop as stop:
        if boundary is None:
            return MemberAttemptOutcome(
                "invalid", None, None, None, None, None, None, stop,
                retry_counters, None,
                {
                    "kind": "preflight_refinement_stop",
                    "exception_type": type(stop).__name__,
                    "cause_type": type(stop.cause).__name__
                    if stop.cause is not None else "None",
                },
                True,
            )
        restored = _safe_restore(member, boundary)
        plan = planned
        if stop.source_retry is not None:
            next_count = retry_counters.source_current + 1
            counters = RetryCounters(next_count, retry_counters.source_total + 1, retry_counters.cfl_current, retry_counters.cfl_total, retry_counters.source_cap, retry_counters.cfl_cap)
            return MemberAttemptOutcome(
                "source_retry_exhausted" if next_count > retry_counters.source_cap else "source_retry_required",
                plan, None, None, None, None, None, stop, counters,
                None if next_count > retry_counters.source_cap else _plan_width(plan) * .5,
                _source_evidence(stop, plan), restored,
            )
        if isinstance(stop.cause, CFLRetryRequired):
            next_count = retry_counters.cfl_current + 1
            counters = RetryCounters(retry_counters.source_current, retry_counters.source_total, next_count, retry_counters.cfl_total + 1, retry_counters.source_cap, retry_counters.cfl_cap)
            return MemberAttemptOutcome(
                "cfl_retry_exhausted" if next_count > retry_counters.cfl_cap else "cfl_retry_required",
                plan, None, None, None, None, None, stop, counters,
                None if next_count > retry_counters.cfl_cap else _plan_width(plan) * .5,
                _cfl_evidence(stop, plan), restored,
            )
        if stop.cause is not None and _scientific(stop.cause):
            return MemberAttemptOutcome(
                "scientific_stop", plan, None, None, None, None, None,
                stop.cause, retry_counters, None,
                {
                    "kind": "scientific_stop",
                    "exception_type": type(stop.cause).__name__,
                    "reason": str(getattr(stop.cause, "reason", None)),
                    "failures": list(getattr(stop.cause, "failures", ())),
                }, restored,
            )
        return MemberAttemptOutcome(
            "invalid", plan, None, None, None, None, None, stop,
            retry_counters, None,
            {
                "kind": "unowned_refinement_stop",
                "exception_type": type(stop.cause).__name__
                if stop.cause is not None else "None",
            }, restored,
        )
    except Proto14SourceRetryExhausted as error:
        if boundary is None:
            return MemberAttemptOutcome("invalid", None, None, None, None, None, None, error, retry_counters, None, {"kind": "preflight_source_exhaustion"}, True)
        restored = _safe_restore(member, boundary)
        return MemberAttemptOutcome("source_retry_exhausted", None, None, None, None, None, None, error, retry_counters, None, {"kind": "source_retry_exhaustion", "evidence": _jsonable(error.evidence)}, restored)
    except (GR0RuntimeStop, Proto7TerminalStop) as error:
        restored = True if boundary is None else _safe_restore(member, boundary)
        return MemberAttemptOutcome(
            "scientific_stop", None, None, None, None, None, None, error,
            retry_counters, None,
            {
                "kind": "scientific_stop",
                "exception_type": type(error).__name__,
                "reason": str(getattr(error, "reason", None)),
                "failures": list(getattr(error, "failures", ())),
            }, restored,
        )
    except Exception as error:
        restored = True if boundary is None else _safe_restore(member, boundary)
        return MemberAttemptOutcome("invalid", None, None, None, None, None, None, error, retry_counters, None, {"kind": "invalid", "exception_type": type(error).__name__}, restored)


__all__ = ["DISPOSITIONS", "Disposition", "HLT16AttemptError", "MemberAttemptOutcome", "RetryCounters", "attempt_once"]
