"""One-attempt coordinator for the authenticated HLT16 GR-0 campaign.

This is deliberately a *coordinator*, not another integrator.  It joins the
pure PRO19 freeze, the HLT16 store, codec, lifecycle, and one-attempt adapter
into one durable transition.  In particular it has no loop, no candidate
branch import, and no authority to alter an event target or retry limit.

All public functions operate on a complete six-member checkpoint.  A caller
may invoke :func:`execute_one_attempt` repeatedly, but each invocation makes
at most one proposal for one named member and either publishes its exact
successor checkpoint or raises before publishing anything new.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import TYPE_CHECKING, Any, Callable, Mapping

from . import proto15_runtime as p15
from . import hlt16_lifecycle as lifecycle
from .hlt16_campaign_schema import (
    CampaignCheckpoint, MemberCampaignState, PROTOCOL, journal_envelope,
)
from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError
from .hlt16_member_codec import (
    GEN0_PROTOCOL, RUNTIME_PROTOCOL, HLT16MemberSnapshot, encode_member,
    restore_member as restore_snapshot,
)
from .hlt16_progression_attempt import MemberAttemptOutcome, RetryCounters, attempt_once
from .hlt16_state_store import authorize_gen0_import, semantic_state_sha256
from .numerical_engine import array_content_sha256
from .proto17_pure_construction import MEMBER_KEYS
from .proto19_progression_contract import ProgressionPlan, TimeIdentity

if TYPE_CHECKING:
    from .proto19_progression_inputs import ReconstructedGR0MemberSet


class HLT16CampaignRuntimeError(ValueError):
    """A coordinator input cannot produce a lawful durable transition."""


@dataclass(frozen=True, slots=True)
class DurableAttempt:
    """The exact checkpoint and record suffix created by one call."""

    disposition: str
    checkpoint: CampaignCheckpoint
    record_sha256s: tuple[str, ...]


def _time(value: object) -> dict[str, str]:
    # The prospective plan deliberately carries a typed rational/binary64 pair,
    # while accepted numerical outcomes carry scalars.  Normalize both through
    # PROTO15's one convention-locked validator instead of asking the runtime
    # caller to erase the plan's exact rational identity first.
    if isinstance(value, TimeIdentity):
        value = {
            "rational": value.rational,
            "binary64_hex": value.binary64_hex,
        }
    try:
        return p15._time_identity(value)
    except Exception as exc:
        raise HLT16CampaignRuntimeError("time identity differs") from exc


def _cursor(value: Mapping[str, Any]) -> p15.Proto15Cursor:
    try:
        result = p15.Proto15Cursor(dict(value)); result.validate()
    except Exception as exc:
        raise HLT16CampaignRuntimeError("cursor differs") from exc
    if result.payload.get("protocol_artifact_id") != PROTOCOL:
        raise HLT16CampaignRuntimeError("cursor protocol differs")
    return result


def _ledger(value: Mapping[str, Any]):
    try:
        return p15._ledger_from_mapping(value)
    except Exception as exc:
        raise HLT16CampaignRuntimeError("TDG6 ledger differs") from exc


def _successor(previous: CampaignCheckpoint, *, members: Mapping[str, MemberCampaignState],
               records: tuple[Mapping[str, Any], ...], disposition: str = "nonterminal",
               terminal: Mapping[str, Any] | None = None, event: int | None = None,
               target: Mapping[str, Any] | None = None) -> CampaignCheckpoint:
    if not records:
        raise HLT16CampaignRuntimeError("successor record suffix is absent")
    current = CampaignCheckpoint(
        plan_sha256=previous.plan_sha256, authorization_commit=previous.authorization_commit,
        protocol=previous.protocol, campaign_id=previous.campaign_id,
        event=previous.event if event is None else event,
        target=dict(previous.target if target is None else target), members=dict(members),
        journal_tip_sha256=str(records[-1]["record_sha256"]), parent_sha256=previous.sha256,
        generation=previous.generation + 1, journal_sequence=int(records[-1]["sequence"]),
        disposition=disposition, terminal=None if terminal is None else dict(terminal),
    )
    return current


def _record(previous: CampaignCheckpoint, *, kind: str, payload: Mapping[str, Any]) -> dict[str, Any]:
    return journal_envelope(
        sequence=previous.journal_sequence + 1,
        campaign_id=previous.campaign_id, generation=previous.generation + 1,
        event=previous.event, previous_record_sha256=previous.journal_tip_sha256,
        kind=kind, payload=payload,
    )


def _record_after(previous: CampaignCheckpoint, first: Mapping[str, Any], *, kind: str,
                  payload: Mapping[str, Any]) -> dict[str, Any]:
    return journal_envelope(
        sequence=int(first["sequence"]) + 1, campaign_id=previous.campaign_id,
        generation=previous.generation + 1, event=previous.event,
        previous_record_sha256=str(first["record_sha256"]), kind=kind, payload=payload,
    )


def initialize_or_recover_generation_one(
    plan: ProgressionPlan,
    store: HLT16CampaignStore,
    reconstructed: "ReconstructedGR0MemberSet",
) -> CampaignCheckpoint:
    """Bridge six authenticated GEN0 imports into one complete HLT16 checkpoint.

    This function only serializes already reconstructed GEN0 arrays.  It never
    makes a numerical proposal and is idempotent across an interruption at any
    store publication boundary.
    """
    # The historical GEN0 bridge alone owns this broad reconstruction module.
    # Persisted resume paths must be able to import the campaign runtime
    # without also importing the retired raw-source/calibration factory.
    from .proto19_progression_inputs import (
        ReconstructedGR0MemberSet as _ReconstructedGR0MemberSet,
    )

    if not isinstance(plan, ProgressionPlan) or not isinstance(store, HLT16CampaignStore):
        raise HLT16CampaignRuntimeError("bridge types differ")
    if not isinstance(reconstructed, _ReconstructedGR0MemberSet):
        raise HLT16CampaignRuntimeError("GEN0 reconstruction differs")
    if (plan.branch, plan.amplitude, plan.common_event_index) != ("GR-0", "3", 23):
        raise HLT16CampaignRuntimeError("candidate branch or first-event identity differs")
    if tuple(reconstructed.members) != MEMBER_KEYS:
        raise HLT16CampaignRuntimeError("GEN0 members differ")
    snapshots: dict[str, HLT16MemberSnapshot] = {}
    authorities = {}
    cursors = reconstructed.evidence.checkpoint["cursors"]
    for key in MEMBER_KEYS:
        member = reconstructed.members[key]
        snapshot = encode_member(member, cursors[key], protocol_artifact_id=PROTOCOL,
                                 campaign_id=plan.campaign_id, generation=0,
                                 provenance_source_protocol=GEN0_PROTOCOL)
        snapshots[key] = snapshot
        authorities[key] = authorize_gen0_import(snapshot)

    def member_factory(key: str, descriptor: str) -> MemberCampaignState:
        ledger = reconstructed.members[key].temporal_ledger
        cursor = lifecycle.bridge_gen0(
            cursors[key],
            ledger,
            descriptor,
            successor_campaign_id=plan.campaign_id,
        )
        return MemberCampaignState(descriptor_sha256=descriptor, cursor=dict(cursor.payload),
                                   ledger=p15._ledger_mapping(ledger))

    try:
        return store.initialize_or_recover_bridge(
            plan_sha256=plan.sha256, authorization_commit=plan.authorization_commit,
            campaign_id=plan.campaign_id, snapshots=snapshots, authorities=authorities,
            member_factory=member_factory, target=_time(plan.target_time), event=plan.common_event_index,
        )
    except (HLT16CampaignStoreError, ValueError, KeyError) as exc:
        raise HLT16CampaignRuntimeError("generation-one bridge differs") from exc


def recover_persisted_generation_one(
    plan: ProgressionPlan, store: HLT16CampaignStore,
) -> CampaignCheckpoint:
    """Return a verified HLT16 successor without reopening GEN0 source arrays.

    This is intentionally the only runtime recovery path used after the first
    bridge checkpoint exists.  It validates the campaign/plan identity from
    the authenticated HLT16 store; array decoding happens later through
    :func:`restore_member_with_overlay` and never through the historical
    PREF26 raw-container importer.
    """
    if not isinstance(plan, ProgressionPlan) or not isinstance(store, HLT16CampaignStore):
        raise HLT16CampaignRuntimeError("persisted recovery types differ")
    try:
        checkpoint = store.authenticated_snapshot().checkpoint
    except (HLT16CampaignStoreError, ValueError, KeyError) as exc:
        raise HLT16CampaignRuntimeError("persisted generation-one checkpoint differs") from exc
    if (
        checkpoint.generation < 1
        or checkpoint.protocol != PROTOCOL
        or checkpoint.plan_sha256 != plan.sha256
        or checkpoint.authorization_commit != plan.authorization_commit
        or checkpoint.campaign_id != plan.campaign_id
        or tuple(checkpoint.members) != MEMBER_KEYS
    ):
        raise HLT16CampaignRuntimeError("persisted generation-one authority differs")
    return checkpoint


def restore_member_with_overlay(
    store: HLT16CampaignStore, checkpoint: CampaignCheckpoint, member: Any,
    *, key: str,
) -> Any:
    """Restore one persisted accepted state and apply its checkpoint retry overlay."""
    if key not in MEMBER_KEYS or checkpoint.members[key].descriptor_sha256 != _cursor(checkpoint.members[key].cursor).payload["accepted_state_sha256"]:
        raise HLT16CampaignRuntimeError("member/checkpoint identity differs")
    state = checkpoint.members[key]
    try:
        persisted = store.load_state(state.descriptor_sha256)
        snapshot = HLT16MemberSnapshot(persisted.arrays, persisted.descriptor["metadata"])
        restore_snapshot(member, snapshot)
    except Exception as exc:
        raise HLT16CampaignRuntimeError("persisted member cannot be restored") from exc
    # Persisted arrays own the physical boundary; the checkpoint owns retry
    # state.  GEN0 descriptors contain the sealed historical counter baseline,
    # so a pre-acceptance retry suffix must be overlaid once.  Every evolved
    # descriptor is written *after* that overlay and consequently already
    # contains baseline plus the campaign total; adding again would double
    # count after a restart.
    generation = persisted.descriptor["metadata"]["runtime_identity"]["generation"]
    if generation == 0:
        member.source_retry_count += state.source_total
        member.CFL_retry_count += state.cfl_total
    elif (member.source_retry_count < state.source_total
          or member.CFL_retry_count < state.cfl_total):
        raise HLT16CampaignRuntimeError("retry overlay exceeds persisted member counter")
    if member.time.hex() != _cursor(state.cursor).payload["accepted_boundary_time"]["binary64_hex"]:
        raise HLT16CampaignRuntimeError("persisted member time differs")
    member.temporal_ledger = _ledger(state.ledger)
    return member


def choose_requested_cap(checkpoint: CampaignCheckpoint, key: str, member: Any) -> float:
    """Choose the one lawful cap: pending owner first, otherwise frozen CFL cap."""
    state = checkpoint.members[key]; cursor = _cursor(state.cursor)
    target = float.fromhex(str(checkpoint.target["binary64_hex"]))
    remaining = target - float(member.time)
    if not math.isfinite(remaining) or remaining <= 0.0:
        raise HLT16CampaignRuntimeError("member has no positive event remainder")
    if state.pending_owner in {"source", "cfl"}:
        if cursor.mode != "FRESH_READY" or state.pending_cap_hex is None:
            raise HLT16CampaignRuntimeError("source/CFL retry ownership differs")
        cap = float.fromhex(state.pending_cap_hex)
    elif state.pending_owner == "temporal":
        pending = cursor.payload.get("retry_successor_payload_or_none")
        if cursor.mode != "RETRY_PENDING" or not isinstance(pending, Mapping):
            raise HLT16CampaignRuntimeError("temporal retry ownership differs")
        cap = float.fromhex(str(pending["half_cap"]))
        if state.pending_cap_hex != cap.hex():
            raise HLT16CampaignRuntimeError("temporal retry cap differs")
    elif state.pending_owner is None:
        if cursor.mode != "FRESH_READY":
            raise HLT16CampaignRuntimeError("fresh retry cursor differs")
        speed = max(float(member.transaction.causal_state.previous_speed_upper), math.ulp(0.0))
        cap = min(remaining, float(member.transaction.cfl_maximum) * float(member.initial.grid.spacing) / speed)
    else:
        raise HLT16CampaignRuntimeError("unknown retry owner")
    if not math.isfinite(cap) or cap <= 0.0 or cap > remaining:
        raise HLT16CampaignRuntimeError("requested cap differs from frozen event remainder")
    return cap


def _counters(state: MemberCampaignState) -> RetryCounters:
    return RetryCounters(state.source_current, state.source_total, state.cfl_current, state.cfl_total)


def _retry_state(previous: MemberCampaignState, outcome: MemberAttemptOutcome, owner: str) -> MemberCampaignState:
    counters = outcome.retry_counters
    if owner == "source":
        return replace(previous, source_current=counters.source_current, source_total=counters.source_total,
                       pending_owner="source", pending_cap_hex=None if outcome.next_cap is None else outcome.next_cap.hex())
    return replace(previous, cfl_current=counters.cfl_current, cfl_total=counters.cfl_total,
                   pending_owner="cfl", pending_cap_hex=None if outcome.next_cap is None else outcome.next_cap.hex())


def _terminal_records(previous: CampaignCheckpoint, key: str, *, reason: str, owner: str,
                      classification: str, evidence: Mapping[str, Any], rejection: str | None = None,
                      pre_record: Mapping[str, Any] | None = None) -> tuple[DurableAttempt, tuple[dict[str, Any], ...]]:
    """Build one terminal checkpoint and its exact canonical journal suffix.

    The terminal checkpoint derives its journal tip from these exact records.
    Returning them to the caller prevents a second independently assembled
    terminal-lock envelope from becoming an accidental second source of truth.
    """
    state = previous.members[key]; cursor = _cursor(state.cursor)
    payload = {"member_key": key, "reason": reason, "classification": classification,
               "descriptor_sha256": state.descriptor_sha256, "cursor_sha256": cursor.sha256,
               "owner": owner, "rejection_record_sha256": rejection, "evidence": dict(evidence)}
    record = (_record_after(previous, pre_record, kind="terminal_lock", payload=payload)
              if pre_record is not None else _record(previous, kind="terminal_lock", payload=payload))
    records: tuple[dict[str, Any], ...] = (record,) if pre_record is None else (dict(pre_record), record)
    terminal = {"classification": classification, "event": previous.event, "target": dict(previous.target),
                "member_key": key, "reason": reason, "owner": owner}
    current = _successor(previous, members=previous.members, records=records,
                         disposition=f"{classification}_terminal", terminal=terminal)
    return (
        DurableAttempt(f"{classification}_terminal", current, tuple(str(item["record_sha256"]) for item in records)),
        records,
    )


def _source_or_cfl_rejection_payload(
    *,
    owner: str,
    key: str,
    state: MemberCampaignState,
    cursor: p15.Proto15Cursor,
    outcome: MemberAttemptOutcome,
    exhausted: bool,
) -> dict[str, Any]:
    """Build the complete, schema-owned source/CFL rejection evidence.

    The attempt adapter owns the numerical evidence and counters; this
    coordinator owns the durable relation to the accepted descriptor and the
    frozen halving contract.  Keeping the record construction in one helper
    prevents the terminal path from silently dropping a field that the retry
    path preserves.
    """
    if owner not in {"source", "cfl"}:
        raise HLT16CampaignRuntimeError("source/CFL rejection owner differs")
    plan = outcome.executed_plan
    if not isinstance(plan, Mapping):
        raise HLT16CampaignRuntimeError("source/CFL retry omits its executed plan")
    macro_hex = plan.get("macro_width_hex")
    if not isinstance(macro_hex, str):
        raise HLT16CampaignRuntimeError("source/CFL executed plan omits macro width")
    try:
        attempted_cap = float.fromhex(macro_hex)
    except ValueError as exc:
        raise HLT16CampaignRuntimeError("source/CFL executed plan macro width differs") from exc
    if (
        not math.isfinite(attempted_cap)
        or attempted_cap <= 0.0
        or attempted_cap.hex() != macro_hex
    ):
        raise HLT16CampaignRuntimeError("source/CFL executed plan macro width differs")
    current = getattr(outcome.retry_counters, f"{owner}_current")
    total = getattr(outcome.retry_counters, f"{owner}_total")
    if exhausted:
        if outcome.next_cap is not None:
            raise HLT16CampaignRuntimeError("exhausted source/CFL retry retained a successor cap")
        next_cap_hex: str | None = None
    else:
        if outcome.next_cap is None or not math.isfinite(outcome.next_cap) or outcome.next_cap <= 0.0:
            raise HLT16CampaignRuntimeError("source/CFL retry omits its exact successor cap")
        if outcome.next_cap.hex() != (attempted_cap * 0.5).hex():
            raise HLT16CampaignRuntimeError("source/CFL successor cap does not halve executed macro width")
        next_cap_hex = outcome.next_cap.hex()
    return {
        "member_key": key,
        "evidence": dict(outcome.evidence),
        "predecessor_descriptor_sha256": state.descriptor_sha256,
        "predecessor_cursor_sha256": cursor.sha256,
        "retry_current": current,
        "retry_total": total,
        "retry_limit": 32,
        "attempted_cap_hex": attempted_cap.hex(),
        "next_cap_hex": next_cap_hex,
        "exhausted": exhausted,
        "no_commit": True,
    }


def execute_one_attempt(
    plan: ProgressionPlan, store: HLT16CampaignStore, checkpoint: CampaignCheckpoint,
    member: Any, *, key: str,
    runner: Callable[..., MemberAttemptOutcome] = attempt_once,
    expected_orphan_semantic_sha256: str | None = None,
) -> DurableAttempt:
    """Execute and durably close exactly one lawfully selected member proposal.

    ``expected_orphan_semantic_sha256`` is used only after an interrupted
    accepted-state publication left a verified payload without its descriptor.
    The replay may adopt that payload only if the independently recomputed
    semantic state identity is byte-for-byte identical.  It cannot turn an
    unrelated rerun into the missing descriptor or authorize another outcome.
    """
    if not isinstance(plan, ProgressionPlan) or checkpoint.plan_sha256 != plan.sha256:
        raise HLT16CampaignRuntimeError("plan/checkpoint authority differs")
    if key not in MEMBER_KEYS or checkpoint.terminal is not None:
        raise HLT16CampaignRuntimeError("member execution is not authorized")
    try:
        # This is not a convenience lookup: it proves that this coordinator
        # owns the active on-disk lease for the exact current predecessor
        # before it can create an orphan payload or append a rejection.
        writer_capability = store.active_writer_capability(checkpoint)
    except HLT16CampaignStoreError as exc:
        raise HLT16CampaignRuntimeError("member execution lacks active writer capability") from exc
    state = checkpoint.members[key]; cursor = _cursor(state.cursor)
    restore_member_with_overlay(store, checkpoint, member, key=key)
    evolution_state_sha256 = array_content_sha256(
        member.state.u, member.state.p, member.state.q
    )
    cap = choose_requested_cap(checkpoint, key, member)
    # TDG6 invokes this sink before it installs a retry ledger.  The rejection
    # must therefore be a durable append-only record at that exact point, not
    # merely a Python value carried until `attempt_once` returns.  A crash
    # after this append but before the complete successor checkpoint is a
    # recognized uncheckpointed suffix, never an invitation to rerun it.
    captured: list[dict[str, Any]] = []

    def durable_temporal_rejection(evidence: Mapping[str, object]) -> None:
        if captured:
            raise HLT16CampaignRuntimeError("TDG6 attempt emitted multiple rejection records")
        if evidence.get("initial_state_sha256") != evolution_state_sha256:
            raise HLT16CampaignRuntimeError(
                "TDG6 rejection differs from the restored evolution state"
            )
        record = _record(checkpoint, kind="tdg6_rejection", payload={
            "member_key": key, "predecessor_descriptor_sha256": state.descriptor_sha256,
            "predecessor_cursor_sha256": cursor.sha256, "evidence": dict(evidence),
        })
        store.publish_journal(record, capability=writer_capability)
        captured.append(record)

    outcome = runner(member, cursor, float.fromhex(str(checkpoint.target["binary64_hex"])),
                     descriptor_sha256=state.descriptor_sha256, requested_cap=cap,
                     retry_counters=_counters(state), temporal_rejection_sink=durable_temporal_rejection)
    if not isinstance(outcome, MemberAttemptOutcome):
        raise HLT16CampaignRuntimeError("attempt adapter returned no typed outcome")
    # A TDG6 sink is allowed to fire only for the two temporal outcomes.  If a
    # malformed adapter emits a durable rejection and then reports another
    # disposition, preserve the append-only suffix and stop fail-closed rather
    # than attaching a contradictory successor checkpoint to it.
    if captured and outcome.disposition not in {"temporal_retry_required", "temporal_retry_exhausted"}:
        raise HLT16CampaignRuntimeError("attempt disposition contradicts its durable TDG6 rejection")
    if expected_orphan_semantic_sha256 is not None:
        if (not isinstance(expected_orphan_semantic_sha256, str)
                or len(expected_orphan_semantic_sha256) != 64
                or expected_orphan_semantic_sha256.lower() != expected_orphan_semantic_sha256):
            raise HLT16CampaignRuntimeError("orphan semantic identity differs")
        try:
            int(expected_orphan_semantic_sha256, 16)
        except ValueError as exc:
            raise HLT16CampaignRuntimeError("orphan semantic identity differs") from exc
        if outcome.disposition != "accepted_fine":
            raise HLT16CampaignRuntimeError("orphan replay did not reproduce an accepted fine state")
    if outcome.disposition == "accepted_fine":
        if outcome.accepted_ledger is None or outcome.executed_plan is None or outcome.accepted_time is None:
            raise HLT16CampaignRuntimeError("accepted attempt omits its state")
        # The descriptor is created before its successor checkpoint is
        # published, but it belongs to that successor generation.  The store
        # independently verifies this relation when it binds the payload to
        # the accepted-fine journal record.
        snapshot = encode_member(member, dict(cursor.payload), protocol_artifact_id=PROTOCOL,
                                 campaign_id=checkpoint.campaign_id, generation=checkpoint.generation + 1,
                                 provenance_source_protocol=RUNTIME_PROTOCOL, executed_plan=outcome.executed_plan)
        if (expected_orphan_semantic_sha256 is not None
                and semantic_state_sha256(snapshot.arrays) != expected_orphan_semantic_sha256):
            raise HLT16CampaignRuntimeError(
                "replayed accepted state differs from the authenticated orphan payload"
            )
        persisted = store.persist_snapshot(
            snapshot, checkpoint=checkpoint, capability=writer_capability,
        )
        successor_cursor = lifecycle.close_fine(cursor, _ledger(state.ledger), outcome.accepted_ledger,
                                                 persisted.descriptor_sha256, outcome.executed_plan,
                                                 _time(outcome.accepted_time), outcome.accepted_step_index,
                                                 outcome.accepted_transaction_serial,
                                                 checkpoint.journal_tip_sha256)
        after = MemberCampaignState(persisted.descriptor_sha256, dict(successor_cursor.payload),
                                    p15._ledger_mapping(outcome.accepted_ledger), 0, state.source_total,
                                    0, state.cfl_total, None, None)
        members = dict(checkpoint.members); members[key] = after
        record = _record(checkpoint, kind="accepted_fine", payload={
            "member_key": key, "predecessor_cursor_sha256": cursor.sha256,
            "successor_cursor": dict(successor_cursor.payload), "successor_cursor_sha256": successor_cursor.sha256,
            "descriptor_sha256": persisted.descriptor_sha256,
            "TDG6_ledger_sha256": successor_cursor.payload["TDG6_ledger_sha256"],
            "executed_plan": outcome.executed_plan,
        })
        current = _successor(checkpoint, members=members, records=(record,))
        _publish(store, checkpoint, current, (record,), capability=writer_capability)
        return DurableAttempt(outcome.disposition, current, (record["record_sha256"],))
    if outcome.disposition in {"source_retry_required", "cfl_retry_required"}:
        owner = "source" if outcome.disposition.startswith("source") else "cfl"
        after = _retry_state(state, outcome, owner); members = dict(checkpoint.members); members[key] = after
        record = _record(checkpoint, kind=f"{owner}_rejection", payload=_source_or_cfl_rejection_payload(
            owner=owner, key=key, state=state, cursor=cursor, outcome=outcome,
            exhausted=False,
        ))
        current = _successor(checkpoint, members=members, records=(record,))
        _publish(store, checkpoint, current, (record,), capability=writer_capability)
        return DurableAttempt(outcome.disposition, current, (record["record_sha256"],))
    if outcome.disposition in {"source_retry_exhausted", "cfl_retry_exhausted"}:
        owner = "source" if outcome.disposition.startswith("source") else "cfl"
        rejection = _record(checkpoint, kind=f"{owner}_rejection", payload=_source_or_cfl_rejection_payload(
            owner=owner, key=key, state=state, cursor=cursor, outcome=outcome,
            exhausted=True,
        ))
        durable, records = _terminal_records(
            checkpoint, key, reason=outcome.disposition, owner=owner, classification="invalid",
            evidence=outcome.evidence, rejection=rejection["record_sha256"], pre_record=rejection,
        )
        _publish(store, checkpoint, durable.checkpoint, records, capability=writer_capability)
        store.publish_terminal_lock(durable.checkpoint, capability=writer_capability)
        return durable
    if outcome.disposition == "temporal_retry_required":
        if (len(captured) != 1 or outcome.retry_outcome is None or outcome.executed_plan is None
                or outcome.next_cap is None or not math.isfinite(outcome.next_cap) or outcome.next_cap <= 0.0):
            raise HLT16CampaignRuntimeError("TDG6 retry evidence was not durably captured")
        rejection = captured[0]
        successor_cursor = lifecycle.close_retry(cursor, _ledger(state.ledger), outcome.retry_outcome,
                                                  outcome.executed_plan,
                                                  p15.plan_forward_proto14_subdivision(
                                                      member.time,
                                                      float.fromhex(str(checkpoint.target["binary64_hex"])),
                                                      outcome.next_cap,
                                                      minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
                                                  ),
                                                  rejection["record_sha256"], checkpoint.journal_tip_sha256,
                                                  evolution_state_sha256=evolution_state_sha256)
        after = replace(state, cursor=dict(successor_cursor.payload), ledger=p15._ledger_mapping(outcome.retry_outcome.updated_ledger),
                        pending_owner="temporal", pending_cap_hex=outcome.next_cap.hex())
        members = dict(checkpoint.members); members[key] = after
        transition = _record_after(checkpoint, rejection, kind="cursor_transition", payload={"member_key": key,
            "predecessor_cursor_sha256": cursor.sha256, "rejection_record_sha256": rejection["record_sha256"],
            "successor_cursor": dict(successor_cursor.payload), "successor_cursor_sha256": successor_cursor.sha256})
        current = _successor(checkpoint, members=members, records=(rejection, transition))
        _publish(store, checkpoint, current, (rejection, transition), capability=writer_capability)
        return DurableAttempt(outcome.disposition, current, (rejection["record_sha256"], transition["record_sha256"]))
    if outcome.disposition == "temporal_retry_exhausted":
        if len(captured) != 1:
            raise HLT16CampaignRuntimeError("TDG6 exhaustion lacks captured rejection")
        rejection = captured[0]
        # The lifecycle owner validates both fresh and pending executable
        # cursor edges.  The checkpoint still retains the last accepted member
        # state: the reconstructed rejected ledger belongs only to terminal
        # evidence, never to a restartable physical state.
        try:
            lifecycle.retry_exhaustion(
                cursor, _ledger(state.ledger), outcome.retry_outcome,
                rejection["record_sha256"], checkpoint.journal_tip_sha256,
                evolution_state_sha256=evolution_state_sha256,
            )
        except Exception as exc:
            raise HLT16CampaignRuntimeError(
                "TDG6 exhaustion has no validated lifecycle owner"
            ) from exc
        durable, records = _terminal_records(
            checkpoint, key, reason=outcome.disposition, owner="temporal", classification="invalid",
            evidence=rejection["payload"]["evidence"], rejection=rejection["record_sha256"], pre_record=rejection,
        )
        _publish(store, checkpoint, durable.checkpoint, records, capability=writer_capability); store.publish_terminal_lock(durable.checkpoint, capability=writer_capability)
        return durable
    if outcome.disposition in {"scientific_stop", "invalid"}:
        classification = "scientific" if outcome.disposition == "scientific_stop" else "invalid"
        durable, records = _terminal_records(
            checkpoint, key, reason=outcome.disposition, owner=classification,
            classification=classification, evidence=outcome.evidence,
        )
        _publish(store, checkpoint, durable.checkpoint, records, capability=writer_capability); store.publish_terminal_lock(durable.checkpoint, capability=writer_capability)
        return durable
    raise HLT16CampaignRuntimeError("attempt disposition is unowned")


def _publish(store: HLT16CampaignStore, previous: CampaignCheckpoint, current: CampaignCheckpoint,
             records: tuple[Mapping[str, Any], ...], *, capability: object) -> None:
    """Publish the journal suffix before its complete checkpoint, never vice versa."""
    for record in records:
        store.publish_journal(record, capability=capability)
    store.publish_checkpoint(current, previous=previous, records=records, capability=capability)


def commit_common_event(checkpoint: CampaignCheckpoint, store: HLT16CampaignStore) -> DurableAttempt:
    """Commit event 24 only after all six members are fresh at the exact target."""
    if checkpoint.event != 23 or checkpoint.target != _time(3.0 / 2.0):
        raise HLT16CampaignRuntimeError("this coordinator owns only the first common event")
    if checkpoint.terminal is not None:
        raise HLT16CampaignRuntimeError("common event cannot follow a terminal checkpoint")
    try:
        writer_capability = store.active_writer_capability(checkpoint)
    except HLT16CampaignStoreError as exc:
        raise HLT16CampaignRuntimeError("common event lacks active writer capability") from exc
    if any(_cursor(checkpoint.members[key].cursor).mode != "FRESH_READY" for key in MEMBER_KEYS):
        raise HLT16CampaignRuntimeError("common event members are not all fresh")
    if any(
        checkpoint.members[key].pending_owner is not None
        or checkpoint.members[key].pending_cap_hex is not None
        or _cursor(checkpoint.members[key].cursor).payload.get("retry_successor_payload_or_none") is not None
        for key in MEMBER_KEYS
    ):
        raise HLT16CampaignRuntimeError("common event members retain pending retry state")
    if any(_cursor(checkpoint.members[key].cursor).payload["accepted_boundary_time"] != checkpoint.target for key in MEMBER_KEYS):
        raise HLT16CampaignRuntimeError("common event target is incomplete")
    next_target = _time(25.0 / 16.0)
    record = _record(checkpoint, kind="common_event_commit", payload={"member_keys": list(MEMBER_KEYS),
        "predecessor_event": checkpoint.event, "completed_event": checkpoint.event + 1,
        "completed_target": dict(checkpoint.target), "next_target": next_target,
        "predecessor_checkpoint_sha256": checkpoint.sha256})
    members: dict[str, MemberCampaignState] = {}
    for key in MEMBER_KEYS:
        before = checkpoint.members[key]; old = _cursor(before.cursor); payload = dict(old.payload)
        payload.update({"committed_common_event_index": checkpoint.event + 1, "event_target_time": next_target,
                        "cursor_generation": old.payload["cursor_generation"] + 1, "attempt_serial": 0,
                        "attempt_id": f"{checkpoint.campaign_id}:{key}:{checkpoint.event + 1}:0",
                        "journal_tip_sha256": record["record_sha256"],
                        "cursor_chain_parent_sha256": old.sha256, "mode": "FRESH_READY", "retry_successor_payload_or_none": None})
        payload.pop("cursor_chain_sha256", None)
        # Event commit changes cursor identity fields.  Rebuild through the
        # canonical PROTO15 constructor so ``cursor_chain_sha256`` is derived
        # from the new payload rather than validating a stale predecessor hash.
        fresh = p15._cursor(payload); fresh.validate()
        members[key] = replace(before, cursor=dict(fresh.payload), pending_owner=None, pending_cap_hex=None)
    current = _successor(checkpoint, members=members, records=(record,), disposition="event_complete",
                         event=checkpoint.event + 1, target=next_target)
    _publish(store, checkpoint, current, (record,), capability=writer_capability)
    return DurableAttempt("event_complete", current, (record["record_sha256"],))


__all__ = [
    "DurableAttempt", "HLT16CampaignRuntimeError", "choose_requested_cap",
    "commit_common_event", "execute_one_attempt", "initialize_or_recover_generation_one",
    "recover_persisted_generation_one", "restore_member_with_overlay",
]
