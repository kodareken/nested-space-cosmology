"""Deterministic replay of authenticated HLT16 publication suffixes.

This module never advances a member, replays a PDE proposal, removes a file,
or chooses an alternative branch.  It only converts an already validated store
snapshot into the sole checkpoint that its append-only evidence proves.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Literal, Mapping

from . import hlt16_lifecycle as lifecycle
from . import proto15_runtime as p15
from .hlt16_campaign_schema import CampaignCheckpoint, MemberCampaignState, journal_envelope
from .hlt16_campaign_store import AuthenticatedCampaignSnapshot, HLT16CampaignStore
from .hlt16_member_codec import validate_codec_metadata
from .numerical_engine import array_content_sha256
from .proto17_pure_construction import MEMBER_KEYS


class HLT16RecoveryError(ValueError):
    """Authenticated facts do not determine one safe publication action."""


@dataclass(frozen=True, slots=True)
class RecoveryDecision:
    disposition: Literal[
        "clean_checkpoint", "reconciled_checkpoint", "rerun_from_checkpoint_required",
        "publication_replay_required", "invalid_ambiguous_suffix",
    ]
    checkpoint_sha256: str
    published: bool
    detail: str


@dataclass(frozen=True, slots=True)
class TDG6RecoveryPreview:
    """The sole checkpoint implied by one authenticated TDG6 suffix."""

    predecessor_checkpoint_sha256: str
    checkpoint: CampaignCheckpoint
    records: tuple[Mapping[str, Any], ...]
    evolution_state_sha256: str
    physical_state_advanced: bool = False


def _cursor(value: Mapping[str, Any]) -> p15.Proto15Cursor:
    result = p15.Proto15Cursor(dict(value)); result.validate()
    return result


def _ledger(metadata: Mapping[str, Any]) -> p15.TDG6TemporalLedger:
    checked = validate_codec_metadata(metadata)
    extension = checked["tdg6"]["extension"]
    debit = [float.fromhex(item) for item in checked["tdg6"]["accumulated_debit_hex"]]
    from .tdg6_temporal_admission_runtime import restore_tdg6_checkpoint_extension
    return restore_tdg6_checkpoint_extension(extension, debit)


def _record(previous: CampaignCheckpoint, kind: str, payload: Mapping[str, Any]) -> dict[str, Any]:
    return journal_envelope(
        sequence=previous.journal_sequence + 1, campaign_id=previous.campaign_id,
        generation=previous.generation + 1, event=previous.event,
        previous_record_sha256=previous.journal_tip_sha256, kind=kind, payload=payload,
    )


def _after(
    previous: CampaignCheckpoint,
    records: tuple[Mapping[str, Any], ...],
    members: Mapping[str, MemberCampaignState],
    *,
    disposition: str = "nonterminal",
    terminal: Mapping[str, Any] | None = None,
    event: int | None = None,
    target: Mapping[str, Any] | None = None,
) -> CampaignCheckpoint:
    return CampaignCheckpoint(
        plan_sha256=previous.plan_sha256, authorization_commit=previous.authorization_commit,
        protocol=previous.protocol, campaign_id=previous.campaign_id,
        event=previous.event if event is None else event,
        target=dict(previous.target if target is None else target), members=dict(members),
        journal_tip_sha256=str(records[-1]["record_sha256"]), parent_sha256=previous.sha256,
        generation=previous.generation + 1, journal_sequence=int(records[-1]["sequence"]),
        disposition=disposition, terminal=None if terminal is None else dict(terminal),
    )


def _descriptor_acceptance(
    previous: CampaignCheckpoint,
    descriptor: str,
    metadata: Mapping[str, Any],
) -> tuple[str, MemberCampaignState, dict[str, Any]]:
    checked = validate_codec_metadata(metadata)
    identity = checked["runtime_identity"]
    key = str(identity["member_key"])
    before = previous.members[key]
    predecessor = _cursor(before.cursor)
    if checked["provenance_cursor"]["payload"] != dict(predecessor.payload):
        raise HLT16RecoveryError("orphan descriptor predecessor differs")
    ledger = _ledger(metadata)
    successor = lifecycle.close_fine(
        predecessor, p15._ledger_from_mapping(before.ledger), ledger, descriptor,
        checked["accepted_plan"], identity["accepted_time"],
        checked["counters"]["step_index"], checked["counters"]["transaction_serial"],
        previous.journal_tip_sha256,
    )
    after = MemberCampaignState(
        descriptor_sha256=descriptor, cursor=dict(successor.payload),
        ledger=p15._ledger_mapping(ledger), source_current=0,
        source_total=before.source_total, cfl_current=0, cfl_total=before.cfl_total,
    )
    payload = {
        "member_key": key, "predecessor_cursor_sha256": predecessor.sha256,
        "successor_cursor": dict(successor.payload), "successor_cursor_sha256": successor.sha256,
        "descriptor_sha256": descriptor, "TDG6_ledger_sha256": successor.payload["TDG6_ledger_sha256"],
        "executed_plan": checked["accepted_plan"],
    }
    return key, after, payload


def _terminal(previous: CampaignCheckpoint, records: tuple[Mapping[str, Any], ...]) -> CampaignCheckpoint:
    lock = records[-1]["payload"]
    classification = str(lock["classification"])
    terminal = {"classification": classification, "event": previous.event, "target": dict(previous.target),
                "member_key": lock["member_key"], "reason": lock["reason"], "owner": lock["owner"]}
    return _after(previous, records, previous.members, disposition=f"{classification}_terminal", terminal=terminal)


def _terminal_lock_record(
    previous: CampaignCheckpoint,
    rejection: Mapping[str, Any],
    *,
    owner: str,
    reason: str,
) -> dict[str, Any]:
    payload = rejection["payload"]
    key = payload["member_key"]
    before = previous.members[key]
    cursor = _cursor(before.cursor)
    lock_payload = {
        "member_key": key,
        "reason": reason,
        "classification": "invalid",
        "descriptor_sha256": before.descriptor_sha256,
        "cursor_sha256": cursor.sha256,
        "owner": owner,
        "rejection_record_sha256": rejection["record_sha256"],
        "evidence": payload["evidence"],
    }
    return journal_envelope(
        sequence=int(rejection["sequence"]) + 1,
        campaign_id=previous.campaign_id,
        generation=previous.generation + 1,
        event=previous.event,
        previous_record_sha256=rejection["record_sha256"],
        kind="terminal_lock",
        payload=lock_payload,
    )


def _resolve_tdg6_rejection(
    previous: CampaignCheckpoint,
    rejection: Mapping[str, Any],
    arrays: Mapping[str, Any],
) -> TDG6RecoveryPreview:
    """Derive, but do not publish, the unique TDG6 suffix successor."""
    record = dict(rejection)
    payload = record["payload"]
    key = payload["member_key"]
    before = previous.members[key]
    cursor = _cursor(before.cursor)
    ledger = p15._ledger_from_mapping(before.ledger)
    try:
        physical = array_content_sha256(arrays["u"], arrays["p"], arrays["q"])
        typed = lifecycle.typed_retry_from_evidence(
            cursor,
            ledger,
            payload["evidence"],
            evolution_state_sha256=physical,
        )
        if isinstance(typed, p15.TDG6TemporalRetryExhausted):
            lock_payload = {
                "member_key": key,
                "reason": "temporal_retry_exhausted",
                "classification": "invalid",
                "descriptor_sha256": before.descriptor_sha256,
                "cursor_sha256": cursor.sha256,
                "owner": "temporal",
                "rejection_record_sha256": record["record_sha256"],
                "evidence": payload["evidence"],
            }
            lock = journal_envelope(
                sequence=record["sequence"] + 1,
                campaign_id=previous.campaign_id,
                generation=previous.generation + 1,
                event=previous.event,
                previous_record_sha256=record["record_sha256"],
                kind="terminal_lock",
                payload=lock_payload,
            )
            records = (record, lock)
            current = _terminal(previous, records)
        else:
            accepted_time = float.fromhex(
                cursor.payload["accepted_boundary_time"]["binary64_hex"]
            )
            target_time = float.fromhex(
                cursor.payload["event_target_time"]["binary64_hex"]
            )
            if cursor.mode == "RETRY_PENDING":
                # The pending cursor is the durable owner of the plan that was
                # actually executed.  Its requested cap may differ from the
                # selected macro width by one or more binary64 quanta after the
                # TDG7 lattice applies a conservative reduction.  Replanning
                # from ``attempted_macro_step_size`` would erase that fact and
                # make an otherwise exact second rejection unrecoverable.
                pending = cursor.payload["retry_successor_payload_or_none"]
                if not isinstance(pending, Mapping):
                    raise HLT16RecoveryError(
                        "pending TDG6 cursor omits its durable successor plan"
                    )
                predecessor_plan = pending["successor_plan"]
            else:
                predecessor_plan = p15.plan_forward_proto14_subdivision(
                    accepted_time,
                    target_time,
                    typed.evidence.attempted_macro_step_size,
                    minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
                )
            successor_plan = p15.plan_forward_proto14_subdivision(
                accepted_time,
                target_time,
                typed.evidence.retry_step_size,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
            successor = lifecycle.close_retry(
                cursor,
                ledger,
                typed,
                predecessor_plan,
                successor_plan,
                record["record_sha256"],
                previous.journal_tip_sha256,
                evolution_state_sha256=physical,
            )
            members = dict(previous.members)
            members[key] = replace(
                before,
                cursor=dict(successor.payload),
                ledger=p15._ledger_mapping(typed.updated_ledger),
                pending_owner="temporal",
                pending_cap_hex=typed.retry_step_size.hex(),
            )
            transition_payload = {
                "member_key": key,
                "predecessor_cursor_sha256": cursor.sha256,
                "rejection_record_sha256": record["record_sha256"],
                "successor_cursor": dict(successor.payload),
                "successor_cursor_sha256": successor.sha256,
            }
            transition = journal_envelope(
                sequence=record["sequence"] + 1,
                campaign_id=previous.campaign_id,
                generation=previous.generation + 1,
                event=previous.event,
                previous_record_sha256=record["record_sha256"],
                kind="cursor_transition",
                payload=transition_payload,
            )
            records = (record, transition)
            current = _after(previous, records, members)
    except Exception as exc:
        raise HLT16RecoveryError(
            "authenticated TDG6 suffix did not determine one checkpoint"
        ) from exc
    return TDG6RecoveryPreview(
        predecessor_checkpoint_sha256=previous.sha256,
        checkpoint=current,
        records=records,
        evolution_state_sha256=physical,
    )


def preview_tdg6_recovery(store: HLT16CampaignStore) -> TDG6RecoveryPreview:
    """Read and derive the exact TDG6 recovery edge without publication."""
    snapshot = store.authenticated_snapshot()
    if (
        snapshot.suffix_kinds != ("tdg6_rejection",)
        or snapshot.staging_paths
        or snapshot.orphan_descriptor_sha256 is not None
        or snapshot.orphan_payload_semantic_sha256 is not None
    ):
        raise HLT16RecoveryError("store does not expose one clean TDG6 suffix")
    rejection = snapshot.suffix_records[0]
    key = rejection["payload"]["member_key"]
    descriptor = snapshot.checkpoint.members[key].descriptor_sha256
    persisted = store.load_state(descriptor)
    return _resolve_tdg6_rejection(
        snapshot.checkpoint,
        rejection,
        persisted.arrays,
    )


def reconcile(store: HLT16CampaignStore) -> RecoveryDecision:
    """Publish only the unique checkpoint proved by the current snapshot."""
    snapshot = store.authenticated_snapshot()
    previous = snapshot.checkpoint
    # A clean/suffix classification remains readable without a lease, but any
    # branch below that publishes a checkpoint or terminal lock must be owned
    # by the active writer bound to this exact predecessor.
    capability = None
    if snapshot.suffix_classification != "clean" or (
        previous.disposition in {"scientific_terminal", "invalid_terminal"}
        and not snapshot.terminal_lock_present
    ):
        try:
            capability = store.active_writer_capability(previous)
        except Exception as exc:
            raise HLT16RecoveryError("recovery lacks active writer capability") from exc
    if snapshot.staging_paths:
        return RecoveryDecision("publication_replay_required", previous.sha256, False, "staging suffix requires explicit replay")
    if (
        snapshot.suffix_classification == "clean"
        and previous.disposition in {"scientific_terminal", "invalid_terminal"}
        and not snapshot.terminal_lock_present
    ):
        # The terminal checkpoint itself is already authenticated.  The only
        # lawful completion is its store-owned exact lock, not a new journal
        # record or a replacement checkpoint.
        store.publish_terminal_lock(previous, capability=capability)
        return RecoveryDecision(
            "reconciled_checkpoint", previous.sha256, True, "terminal lock"
        )
    if snapshot.suffix_classification == "clean":
        return RecoveryDecision("clean_checkpoint", previous.sha256, False, "no suffix")
    if snapshot.orphan_payload_semantic_sha256 is not None and snapshot.orphan_descriptor_sha256 is None:
        return RecoveryDecision("rerun_from_checkpoint_required", previous.sha256, False, "payload has no descriptor")

    records = tuple(dict(item) for item in snapshot.suffix_records)
    kinds = snapshot.suffix_kinds
    try:
        if kinds == ("accepted_fine",) or (not kinds and snapshot.orphan_descriptor_sha256 is not None):
            descriptor = (str(records[0]["payload"]["descriptor_sha256"]) if records else snapshot.orphan_descriptor_sha256)
            if descriptor is None:
                raise HLT16RecoveryError("accepted suffix omits descriptor")
            persisted = store.load_state(descriptor)
            key, after, payload = _descriptor_acceptance(previous, descriptor, persisted.descriptor["metadata"])
            expected = _record(previous, "accepted_fine", payload)
            if records and records != (expected,):
                raise HLT16RecoveryError("accepted suffix does not match descriptor")
            members = dict(previous.members); members[key] = after
            current = _after(previous, (expected,), members)
            if not records:
                store.publish_journal(expected, capability=capability)
            store.publish_checkpoint(current, previous=previous, records=(expected,), capability=capability)
            return RecoveryDecision("reconciled_checkpoint", current.sha256, True, "accepted fine")

        if kinds in (("source_rejection",), ("cfl_rejection",)):
            record = records[0]; p = record["payload"]; key = p["member_key"]
            before = previous.members[key]
            owner = "source" if kinds[0] == "source_rejection" else "cfl"
            if p["exhausted"] is True:
                lock = _terminal_lock_record(
                    previous, record, owner=owner,
                    reason=f"{owner}_retry_exhausted",
                )
                records = (record, lock)
                current = _terminal(previous, records)
                store.publish_journal(lock, capability=capability)
            else:
                changes = (
                    {
                        "source_current": p["retry_current"],
                        "source_total": p["retry_total"],
                        "pending_owner": owner,
                        "pending_cap_hex": p["next_cap_hex"],
                    }
                    if owner == "source"
                    else {
                        "cfl_current": p["retry_current"],
                        "cfl_total": p["retry_total"],
                        "pending_owner": owner,
                        "pending_cap_hex": p["next_cap_hex"],
                    }
                )
                members = dict(previous.members)
                members[key] = replace(before, **changes)
                current = _after(previous, records, members)
        elif kinds == ("tdg6_rejection",):
            key = records[0]["payload"]["member_key"]
            persisted = store.load_state(previous.members[key].descriptor_sha256)
            preview = _resolve_tdg6_rejection(
                previous,
                records[0],
                persisted.arrays,
            )
            records = tuple(dict(item) for item in preview.records)
            current = preview.checkpoint
            store.publish_journal(records[1], capability=capability)
        elif kinds == ("tdg6_rejection", "cursor_transition"):
            transition = records[1]["payload"]; key = transition["member_key"]; before = previous.members[key]
            members = dict(previous.members); members[key] = replace(before, cursor=dict(transition["successor_cursor"]), ledger=p15._ledger_mapping(p15._updated_retry_ledger(p15._ledger_from_mapping(before.ledger), records[0]["payload"]["evidence"])), pending_owner="temporal", pending_cap_hex=transition["successor_cursor"]["retry_successor_payload_or_none"]["half_cap"])
            current = _after(previous, records, members)
        elif kinds in (("terminal_lock",), ("source_rejection", "terminal_lock"), ("cfl_rejection", "terminal_lock"), ("tdg6_rejection", "terminal_lock")):
            current = _terminal(previous, records)
        elif kinds == ("common_event_commit",):
            p = records[0]["payload"]; members = {}
            for key in MEMBER_KEYS:
                before = previous.members[key]; old = _cursor(before.cursor); payload = dict(old.payload)
                payload.update({"committed_common_event_index": p["completed_event"], "event_target_time": p["next_target"], "cursor_generation": old.payload["cursor_generation"] + 1, "attempt_serial": 0, "attempt_id": f"{previous.campaign_id}:{key}:{p['completed_event']}:0", "journal_tip_sha256": records[0]["record_sha256"], "cursor_chain_parent_sha256": old.sha256, "mode": "FRESH_READY", "retry_successor_payload_or_none": None}); payload.pop("cursor_chain_sha256", None)
                fresh = p15._cursor(payload)
                fresh.validate()
                members[key] = replace(
                    before,
                    cursor=dict(fresh.payload),
                    pending_owner=None,
                    pending_cap_hex=None,
                )
            current = _after(previous, records, members, disposition="event_complete", event=p["completed_event"], target=p["next_target"])
        else:
            return RecoveryDecision("invalid_ambiguous_suffix", previous.sha256, False, "unrecognized authenticated suffix")
        store.publish_checkpoint(current, previous=previous, records=records, capability=capability)
        if current.disposition in {"scientific_terminal", "invalid_terminal"}:
            store.publish_terminal_lock(current, capability=capability)
        return RecoveryDecision("reconciled_checkpoint", current.sha256, True, "+".join(kinds))
    except Exception as exc:
        raise HLT16RecoveryError("authenticated suffix did not determine one checkpoint") from exc


__all__ = [
    "HLT16RecoveryError",
    "RecoveryDecision",
    "TDG6RecoveryPreview",
    "preview_tdg6_recovery",
    "reconcile",
]
