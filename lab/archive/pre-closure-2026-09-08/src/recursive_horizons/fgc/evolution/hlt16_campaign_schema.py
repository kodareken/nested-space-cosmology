"""Pure, fail-closed HLT16 campaign journal and checkpoint schema.

PROTO17 generation zero remains an external immutable source.  This module
starts at the HLT16 generation-one bridge and has no filesystem, PDE, or
candidate-branch capability.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

from . import hlt16_lifecycle as lifecycle
from . import proto15_runtime as p15
from .proto17_pure_construction import MEMBER_KEYS, ROOT_SHA256, digest

PROTOCOL = "FGC-2-SF1-PROTO18"
JOURNAL_SCHEMA = "FGC-1-HLT16-campaign-journal-v1"
RECORD_KINDS = (
    "genesis_bridge", "source_rejection", "cfl_rejection", "tdg6_rejection",
    "cursor_transition", "accepted_fine", "terminal_lock", "common_event_commit",
)
TERMINAL_DISPOSITIONS = frozenset({"scientific_terminal", "invalid_terminal"})
NONTERMINAL_DISPOSITIONS = frozenset({"nonterminal", "event_complete"})


class HLT16CampaignSchemaError(ValueError):
    """A durable HLT16 claim cannot be reconstructed exactly."""


def _sha(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise HLT16CampaignSchemaError(f"{label} differs")
    try:
        int(value, 16)
    except ValueError as error:
        raise HLT16CampaignSchemaError(f"{label} differs") from error
    return value


def _retry_evolution_state_sha256(value: object) -> str:
    if not isinstance(value, Mapping):
        raise HLT16CampaignSchemaError("TDG6 retry evidence differs")
    return _sha(value.get("initial_state_sha256"), "TDG6 evolution state")


def _commit(value: object) -> str:
    if not isinstance(value, str) or len(value) != 40 or value.lower() != value:
        raise HLT16CampaignSchemaError("authorization commit differs")
    try:
        int(value, 16)
    except ValueError as error:
        raise HLT16CampaignSchemaError("authorization commit differs") from error
    return value


def _integer(value: object, label: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise HLT16CampaignSchemaError(f"{label} differs")
    return value


def _time(value: object) -> dict[str, str]:
    try:
        return p15._time_identity(value)
    except Exception as error:
        raise HLT16CampaignSchemaError("time identity differs") from error


def _cursor(value: object) -> p15.Proto15Cursor:
    if not isinstance(value, Mapping):
        raise HLT16CampaignSchemaError("cursor differs")
    try:
        # ``_cursor`` is PROTO15's constructor: it deliberately recomputes a
        # chain hash.  Durable HLT16 evidence must instead verify the observed
        # hash verbatim, so a one-bit cursor mutation is never normalized away.
        answer = p15.Proto15Cursor(dict(value))
        answer.validate()
    except Exception as error:
        raise HLT16CampaignSchemaError("cursor differs") from error
    if answer.payload["protocol_artifact_id"] != PROTOCOL:
        raise HLT16CampaignSchemaError("cursor protocol differs")
    return answer


def _ledger(value: object) -> p15.TDG6TemporalLedger:
    if not isinstance(value, Mapping):
        raise HLT16CampaignSchemaError("TDG6 ledger differs")
    try:
        return p15._ledger_from_mapping(value)
    except Exception as error:
        raise HLT16CampaignSchemaError("TDG6 ledger differs") from error


def _keys(value: object, expected: frozenset[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != expected:
        raise HLT16CampaignSchemaError(f"{label} differs")
    return value


def _member(value: object) -> str:
    if not isinstance(value, str) or value not in MEMBER_KEYS:
        raise HLT16CampaignSchemaError("member key differs")
    return value


def _cursor_mapping(value: object, key: str | None = None) -> dict[str, Any]:
    cursor = _cursor(value)
    if key is not None and cursor.payload["member_key"] != key:
        raise HLT16CampaignSchemaError("cursor member differs")
    return dict(cursor.payload)


def _hex_positive(value: object, label: str) -> str:
    if not isinstance(value, str):
        raise HLT16CampaignSchemaError(f"{label} differs")
    try:
        number = float.fromhex(value)
    except ValueError as error:
        raise HLT16CampaignSchemaError(f"{label} differs") from error
    if not math.isfinite(number) or not number > 0.0 or number.hex() != value:
        raise HLT16CampaignSchemaError(f"{label} differs")
    return value


def _record_payload(kind: str, value: object) -> dict[str, Any]:
    """Validate a standalone record; predecessor relations are checked later."""
    if kind == "genesis_bridge":
        p = _keys(value, frozenset({
            "external_protocol", "external_checkpoint_sha256", "plan_sha256",
            "member_descriptors", "bridge_cursors",
        }), "genesis bridge payload")
        if p["external_protocol"] != "FGC-2-SF1-PROTO17":
            raise HLT16CampaignSchemaError("external protocol differs")
        _sha(p["external_checkpoint_sha256"], "external checkpoint")
        _sha(p["plan_sha256"], "plan")
        if not isinstance(p["member_descriptors"], Mapping) or set(p["member_descriptors"]) != set(MEMBER_KEYS):
            raise HLT16CampaignSchemaError("genesis descriptors differ")
        if not isinstance(p["bridge_cursors"], Mapping) or set(p["bridge_cursors"]) != set(MEMBER_KEYS):
            raise HLT16CampaignSchemaError("genesis cursors differ")
        for key in MEMBER_KEYS:
            _sha(p["member_descriptors"][key], "genesis descriptor")
            _cursor_mapping(p["bridge_cursors"][key], key)
        return dict(p)

    if kind in {"source_rejection", "cfl_rejection"}:
        p = _keys(value, frozenset({
            "member_key", "evidence", "predecessor_descriptor_sha256",
            "predecessor_cursor_sha256", "retry_current", "retry_total", "retry_limit",
            "attempted_cap_hex", "next_cap_hex", "exhausted", "no_commit",
        }), f"{kind} payload")
        _member(p["member_key"])
        # These owners intentionally preserve their native canonical numerical
        # evidence rather than reimplementing a stale source/CFL field list.
        # The campaign writer must supply the exact asdict/canonical evidence.
        if not isinstance(p["evidence"], Mapping) or p["no_commit"] is not True:
            raise HLT16CampaignSchemaError(f"{kind} evidence/no-commit differs")
        _sha(p["predecessor_descriptor_sha256"], "predecessor descriptor")
        _sha(p["predecessor_cursor_sha256"], "predecessor cursor")
        current = _integer(p["retry_current"], "retry current", 1)
        total = _integer(p["retry_total"], "retry total", 1)
        limit = _integer(p["retry_limit"], "retry limit", 1)
        attempted = float.fromhex(_hex_positive(p["attempted_cap_hex"], "attempted cap"))
        if not isinstance(p["exhausted"], bool):
            raise HLT16CampaignSchemaError(f"{kind} exhaustion flag differs")
        if p["exhausted"]:
            if p["next_cap_hex"] is not None or current != limit + 1:
                raise HLT16CampaignSchemaError(f"{kind} exhaustion relation differs")
        else:
            next_cap = float.fromhex(_hex_positive(p["next_cap_hex"], "next cap"))
            if next_cap.hex() != (attempted * 0.5).hex() or current > limit:
                raise HLT16CampaignSchemaError(f"{kind} cap reduction differs")
        if current > total or current > limit + 1 or limit != 32:
            raise HLT16CampaignSchemaError(f"{kind} retry relation differs")
        return {**p, "evidence": dict(p["evidence"])}

    if kind == "tdg6_rejection":
        p = _keys(value, frozenset({
            "member_key", "predecessor_descriptor_sha256",
            "predecessor_cursor_sha256", "evidence",
        }), "TDG6 rejection payload")
        _member(p["member_key"])
        _sha(p["predecessor_descriptor_sha256"], "predecessor descriptor")
        _sha(p["predecessor_cursor_sha256"], "predecessor cursor")
        if not isinstance(p["evidence"], Mapping):
            raise HLT16CampaignSchemaError("TDG6 evidence differs")
        # TDG6's native dataclass mapping deliberately retains tuples, while
        # the durable journal is canonical JSON and therefore represents
        # those sequences as arrays.  Normalize at this single record-
        # construction boundary so the in-memory record validated for its
        # immediate checkpoint is byte-identical to the record recovered
        # from disk.  The lifecycle's downstream evidence validator remains
        # strict and continues to require concrete lists in durable evidence.
        try:
            evidence = p15._json_safe(dict(p["evidence"]))
        except (TypeError, ValueError) as error:
            raise HLT16CampaignSchemaError("TDG6 evidence differs") from error
        if not isinstance(evidence, Mapping):
            raise HLT16CampaignSchemaError("TDG6 evidence differs")
        return {**p, "evidence": dict(evidence)}

    if kind == "cursor_transition":
        p = _keys(value, frozenset({
            "member_key", "predecessor_cursor_sha256", "rejection_record_sha256",
            "successor_cursor", "successor_cursor_sha256",
        }), "cursor transition payload")
        key = _member(p["member_key"])
        _sha(p["predecessor_cursor_sha256"], "predecessor cursor")
        _sha(p["rejection_record_sha256"], "rejection record")
        successor = _cursor_mapping(p["successor_cursor"], key)
        if _sha(p["successor_cursor_sha256"], "successor cursor") != successor["cursor_chain_sha256"]:
            raise HLT16CampaignSchemaError("successor cursor hash differs")
        return {**p, "successor_cursor": successor}

    if kind == "accepted_fine":
        p = _keys(value, frozenset({
            "member_key", "predecessor_cursor_sha256", "successor_cursor",
            "successor_cursor_sha256", "descriptor_sha256", "TDG6_ledger_sha256",
            "executed_plan",
        }), "accepted fine payload")
        key = _member(p["member_key"])
        _sha(p["predecessor_cursor_sha256"], "predecessor cursor")
        successor = _cursor_mapping(p["successor_cursor"], key)
        if _sha(p["successor_cursor_sha256"], "successor cursor") != successor["cursor_chain_sha256"]:
            raise HLT16CampaignSchemaError("successor cursor hash differs")
        if _sha(p["descriptor_sha256"], "accepted descriptor") != successor["accepted_state_sha256"]:
            raise HLT16CampaignSchemaError("accepted descriptor differs from cursor")
        if _sha(p["TDG6_ledger_sha256"], "accepted ledger") != successor["TDG6_ledger_sha256"]:
            raise HLT16CampaignSchemaError("accepted ledger differs from cursor")
        try:
            plan = p15._plan_mapping(p["executed_plan"])
        except Exception as error:
            raise HLT16CampaignSchemaError("accepted TDG7 plan differs") from error
        return {**p, "successor_cursor": successor, "executed_plan": plan}

    if kind == "terminal_lock":
        p = _keys(value, frozenset({
            "member_key", "reason", "classification", "descriptor_sha256", "cursor_sha256",
            "owner", "rejection_record_sha256", "evidence",
        }), "terminal lock payload")
        _member(p["member_key"])
        if (p["classification"] not in {"scientific", "invalid"} or not isinstance(p["reason"], str)
                or not p["reason"] or p["owner"] not in {"source", "cfl", "temporal", "scientific", "invalid"}
                or not isinstance(p["evidence"], Mapping)):
            raise HLT16CampaignSchemaError("terminal lock differs")
        _sha(p["descriptor_sha256"], "terminal descriptor")
        _sha(p["cursor_sha256"], "terminal cursor")
        if p["rejection_record_sha256"] is not None:
            _sha(p["rejection_record_sha256"], "terminal rejection record")
        return {**p, "evidence": dict(p["evidence"])}

    if kind == "common_event_commit":
        p = _keys(value, frozenset({
            "member_keys", "predecessor_event", "completed_event", "completed_target",
            "next_target", "predecessor_checkpoint_sha256",
        }), "common event payload")
        if tuple(p["member_keys"]) != MEMBER_KEYS:
            raise HLT16CampaignSchemaError("common event members differ")
        predecessor = _integer(p["predecessor_event"], "predecessor event", 23)
        if _integer(p["completed_event"], "completed event", 24) != predecessor + 1:
            raise HLT16CampaignSchemaError("common event index differs")
        prior, following = _time(p["completed_target"]), _time(p["next_target"])
        if not float.fromhex(following["binary64_hex"]) > float.fromhex(prior["binary64_hex"]):
            raise HLT16CampaignSchemaError("common target does not advance")
        _sha(p["predecessor_checkpoint_sha256"], "common predecessor checkpoint")
        return {**p, "completed_target": prior, "next_target": following}

    raise HLT16CampaignSchemaError("journal kind differs")


def journal_envelope(*, sequence: int, campaign_id: str, generation: int, event: int,
                     previous_record_sha256: str, kind: str, payload: Mapping[str, Any]) -> dict[str, Any]:
    if kind not in RECORD_KINDS or not isinstance(campaign_id, str) or not campaign_id:
        raise HLT16CampaignSchemaError("journal identity differs")
    body = {
        "schema": JOURNAL_SCHEMA, "sequence": _integer(sequence, "journal sequence"),
        "campaign_id": campaign_id, "generation": _integer(generation, "journal generation", 1),
        "event": _integer(event, "journal event", 23),
        "previous_record_sha256": _sha(previous_record_sha256, "journal parent"),
        "kind": kind, "payload": _record_payload(kind, payload),
    }
    return {**body, "record_sha256": digest(body)}


def validate_journal(value: object) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise HLT16CampaignSchemaError("journal differs")
    body = dict(value)
    observed = body.pop("record_sha256", None)
    expected = journal_envelope(
        sequence=body.get("sequence"), campaign_id=body.get("campaign_id"),
        generation=body.get("generation"), event=body.get("event"),
        previous_record_sha256=body.get("previous_record_sha256"),
        kind=body.get("kind"), payload=body.get("payload"),
    )
    if dict(value) != expected or observed != expected["record_sha256"]:
        raise HLT16CampaignSchemaError("journal hash differs")
    return expected


@dataclass(frozen=True, slots=True)
class MemberCampaignState:
    descriptor_sha256: str
    cursor: Mapping[str, Any]
    ledger: Mapping[str, Any]
    source_current: int = 0
    source_total: int = 0
    cfl_current: int = 0
    cfl_total: int = 0
    pending_owner: str | None = None
    pending_cap_hex: str | None = None

    def validate(self) -> "MemberCampaignState":
        descriptor = _sha(self.descriptor_sha256, "descriptor")
        cursor, ledger = _cursor(self.cursor), _ledger(self.ledger)
        if cursor.payload["accepted_state_sha256"] != descriptor:
            raise HLT16CampaignSchemaError("cursor descriptor differs")
        if cursor.payload["TDG6_ledger_sha256"] != p15._hash(p15._ledger_mapping(ledger)):
            raise HLT16CampaignSchemaError("cursor ledger differs")
        if ledger.last_accepted_time.hex() != cursor.payload["accepted_boundary_time"]["binary64_hex"]:
            raise HLT16CampaignSchemaError("cursor accepted time differs from ledger")
        for value, label in ((self.source_current, "source current"), (self.source_total, "source total"),
                             (self.cfl_current, "CFL current"), (self.cfl_total, "CFL total")):
            _integer(value, label)
        if self.source_current > self.source_total or self.cfl_current > self.cfl_total or self.source_current > 32 or self.cfl_current > 32:
            raise HLT16CampaignSchemaError("retry counters differ")
        if self.pending_owner not in {None, "source", "cfl", "temporal"}:
            raise HLT16CampaignSchemaError("pending owner differs")
        if self.pending_cap_hex is not None:
            _hex_positive(self.pending_cap_hex, "pending cap")
        if cursor.mode == "FRESH_READY":
            if self.pending_owner == "temporal" or ledger.current_macro_step_temporal_retry_count:
                raise HLT16CampaignSchemaError("fresh retry state differs")
            if self.pending_owner is None and self.pending_cap_hex is not None:
                raise HLT16CampaignSchemaError("fresh pending cap differs")
            if self.pending_owner in {"source", "cfl"} and self.pending_cap_hex is None:
                raise HLT16CampaignSchemaError("fresh owner cap differs")
        else:
            try:
                pending = cursor.payload["retry_successor_payload_or_none"]
                if not isinstance(pending, Mapping):
                    raise HLT16CampaignSchemaError("pending retry relation differs")
                lifecycle.validate_retry_payload(
                    cursor,
                    ledger,
                    evolution_state_sha256=_retry_evolution_state_sha256(
                        pending.get("TDG6_rejection_evidence")
                    ),
                )
            except Exception as error:
                raise HLT16CampaignSchemaError("pending retry relation differs") from error
            pending = cursor.payload["retry_successor_payload_or_none"]
            assert isinstance(pending, Mapping)
            if self.pending_owner != "temporal" or self.pending_cap_hex != pending["half_cap"]:
                raise HLT16CampaignSchemaError("pending cap differs")
        return self

    def mapping(self) -> dict[str, Any]:
        self.validate()
        return {
            "descriptor_sha256": self.descriptor_sha256, "cursor": dict(self.cursor), "ledger": dict(self.ledger),
            "source_current": self.source_current, "source_total": self.source_total,
            "cfl_current": self.cfl_current, "cfl_total": self.cfl_total,
            "pending_owner": self.pending_owner, "pending_cap_hex": self.pending_cap_hex,
        }

    @classmethod
    def from_mapping(cls, value: object) -> "MemberCampaignState":
        mapping = _keys(value, frozenset({
            "descriptor_sha256", "cursor", "ledger", "source_current", "source_total",
            "cfl_current", "cfl_total", "pending_owner", "pending_cap_hex",
        }), "member state")
        return cls(**dict(mapping)).validate()


@dataclass(frozen=True, slots=True)
class CampaignCheckpoint:
    plan_sha256: str
    authorization_commit: str
    protocol: str
    campaign_id: str
    event: int
    target: Mapping[str, Any]
    members: Mapping[str, MemberCampaignState]
    journal_tip_sha256: str
    parent_sha256: str
    generation: int
    journal_sequence: int = -1
    disposition: str = "nonterminal"
    terminal: Mapping[str, Any] | None = None

    def validate(self) -> "CampaignCheckpoint":
        _sha(self.plan_sha256, "plan"); _commit(self.authorization_commit)
        _sha(self.journal_tip_sha256, "journal tip"); _sha(self.parent_sha256, "checkpoint parent")
        target = _time(self.target)
        if self.protocol != PROTOCOL or not isinstance(self.campaign_id, str) or not self.campaign_id:
            raise HLT16CampaignSchemaError("checkpoint identity differs")
        _integer(self.event, "checkpoint event", 23); _integer(self.generation, "checkpoint generation", 1)
        if isinstance(self.journal_sequence, bool) or not isinstance(self.journal_sequence, int) or self.journal_sequence < 0:
            raise HLT16CampaignSchemaError("checkpoint sequence differs")
        if set(self.members) != set(MEMBER_KEYS):
            raise HLT16CampaignSchemaError("checkpoint must contain six members")
        for key in MEMBER_KEYS:
            state = self.members[key]
            if not isinstance(state, MemberCampaignState):
                raise HLT16CampaignSchemaError("member state type differs")
            state.validate()
            cursor = _cursor(state.cursor)
            if (cursor.payload["campaign_id"] != self.campaign_id or cursor.payload["member_key"] != key
                    or cursor.payload["committed_common_event_index"] != self.event
                    or cursor.payload["event_target_time"] != target):
                raise HLT16CampaignSchemaError("member checkpoint relation differs")
        if self.disposition not in TERMINAL_DISPOSITIONS | NONTERMINAL_DISPOSITIONS:
            raise HLT16CampaignSchemaError("checkpoint disposition differs")
        if (self.disposition in TERMINAL_DISPOSITIONS) != (self.terminal is not None):
            raise HLT16CampaignSchemaError("terminal field differs")
        if self.terminal is not None:
            terminal = _keys(self.terminal, frozenset({"classification", "event", "target", "member_key", "reason", "owner"}), "terminal")
            expected = "scientific" if self.disposition == "scientific_terminal" else "invalid"
            owners = ({"scientific"} if expected == "scientific"
                      else {"source", "cfl", "temporal", "invalid"})
            if (terminal["classification"] != expected or terminal["event"] != self.event
                    or terminal["target"] != target or not isinstance(terminal["reason"], str)
                    or not terminal["reason"] or terminal["owner"] not in owners):
                raise HLT16CampaignSchemaError("terminal relation differs")
            _member(terminal["member_key"])
        if self.generation == 1 and (self.journal_sequence != 0 or self.journal_tip_sha256 == ROOT_SHA256
                                     or self.parent_sha256 == ROOT_SHA256):
            raise HLT16CampaignSchemaError("generation-one bridge differs")
        return self

    def body(self) -> dict[str, Any]:
        self.validate()
        return {
            "plan_sha256": self.plan_sha256, "authorization_commit": self.authorization_commit,
            "protocol": self.protocol, "campaign_id": self.campaign_id, "event": self.event,
            "target": _time(self.target), "members": {key: self.members[key].mapping() for key in MEMBER_KEYS},
            "journal_tip_sha256": self.journal_tip_sha256, "parent_sha256": self.parent_sha256,
            "generation": self.generation, "journal_sequence": self.journal_sequence,
            "disposition": self.disposition, "terminal": None if self.terminal is None else dict(self.terminal),
        }

    @property
    def sha256(self) -> str:
        return digest(self.body())

    def mapping(self) -> dict[str, Any]:
        return {**self.body(), "checkpoint_sha256": self.sha256}

    @classmethod
    def from_mapping(cls, value: object) -> "CampaignCheckpoint":
        mapping = _keys(value, frozenset({
            "plan_sha256", "authorization_commit", "protocol", "campaign_id", "event", "target", "members",
            "journal_tip_sha256", "parent_sha256", "generation", "journal_sequence", "disposition",
            "terminal", "checkpoint_sha256",
        }), "checkpoint")
        members = mapping["members"]
        if not isinstance(members, Mapping) or set(members) != set(MEMBER_KEYS):
            raise HLT16CampaignSchemaError("checkpoint members differ")
        result = cls(
            **{key: mapping[key] for key in mapping if key not in {"members", "checkpoint_sha256"}},
            members={key: MemberCampaignState.from_mapping(members[key]) for key in MEMBER_KEYS},
        )
        if mapping["checkpoint_sha256"] != result.sha256:
            raise HLT16CampaignSchemaError("checkpoint hash differs")
        return result


def _changed(previous: CampaignCheckpoint, current: CampaignCheckpoint) -> tuple[str, ...]:
    # Checkpoints were validated by the caller.  Compare their serializable
    # fields directly; rejected proposals never mutate restartable members.
    def raw(state: MemberCampaignState) -> dict[str, Any]:
        return {
            "descriptor_sha256": state.descriptor_sha256, "cursor": dict(state.cursor), "ledger": dict(state.ledger),
            "source_current": state.source_current, "source_total": state.source_total,
            "cfl_current": state.cfl_current, "cfl_total": state.cfl_total,
            "pending_owner": state.pending_owner, "pending_cap_hex": state.pending_cap_hex,
        }
    return tuple(key for key in MEMBER_KEYS if raw(previous.members[key]) != raw(current.members[key]))


def _record_chain(previous: CampaignCheckpoint, current: CampaignCheckpoint,
                  records: Sequence[Mapping[str, Any]]) -> tuple[dict[str, Any], ...]:
    if not records:
        raise HLT16CampaignSchemaError("checkpoint successor lacks records")
    parsed = tuple(validate_journal(record) for record in records)
    sequence, parent = previous.journal_sequence + 1, previous.journal_tip_sha256
    for record in parsed:
        if (record["sequence"] != sequence or record["previous_record_sha256"] != parent
                or record["campaign_id"] != previous.campaign_id or record["generation"] != current.generation):
            raise HLT16CampaignSchemaError("journal successor chain differs")
        sequence, parent = sequence + 1, record["record_sha256"]
    if current.journal_sequence != parsed[-1]["sequence"] or current.journal_tip_sha256 != parent:
        raise HLT16CampaignSchemaError("checkpoint journal tip differs")
    # All ordinary records describe the event being progressed.  A common-event
    # commit is written at the completed event, then advances the checkpoint to
    # the next event.
    if any(record["event"] != previous.event for record in parsed):
        raise HLT16CampaignSchemaError("journal event differs from predecessor")
    return parsed


def validate_checkpoint_successor(previous: CampaignCheckpoint, current: CampaignCheckpoint,
                                  records: Sequence[Mapping[str, Any]]) -> None:
    """Prove a complete checkpoint successor and its exact journal suffix."""
    previous.validate(); current.validate()
    if current.generation != previous.generation + 1 or current.parent_sha256 != previous.sha256:
        raise HLT16CampaignSchemaError("checkpoint parent/generation differs")
    if (current.plan_sha256 != previous.plan_sha256 or current.authorization_commit != previous.authorization_commit
            or current.protocol != previous.protocol or current.campaign_id != previous.campaign_id):
        raise HLT16CampaignSchemaError("checkpoint authority differs")
    parsed, kinds, changed = _record_chain(previous, current, records), (), _changed(previous, current)
    kinds = tuple(record["kind"] for record in parsed)
    if "genesis_bridge" in kinds:
        raise HLT16CampaignSchemaError("genesis bridge cannot be successor")

    if kinds == ("accepted_fine",):
        p = parsed[0]["payload"]; key = p["member_key"]; before, after = previous.members[key], current.members[key]
        if (current.event != previous.event or current.target != previous.target
                or current.disposition != "nonterminal" or changed != (key,)):
            raise HLT16CampaignSchemaError("accepted scope differs")
        if (p["predecessor_cursor_sha256"] != _cursor(before.cursor).sha256
                or p["successor_cursor_sha256"] != _cursor(after.cursor).sha256
                or p["descriptor_sha256"] != after.descriptor_sha256):
            raise HLT16CampaignSchemaError("accepted identity differs")
        try:
            p15._validated_accepted_plan(_cursor(before.cursor), p["executed_plan"], _cursor(after.cursor).payload["accepted_boundary_time"])
        except Exception as error:
            raise HLT16CampaignSchemaError("accepted plan differs") from error
        old, new = _ledger(before.ledger), _ledger(after.ledger)
        if (after.source_current != 0 or after.cfl_current != 0
                or after.source_total != before.source_total or after.cfl_total != before.cfl_total
                or after.pending_owner is not None or after.pending_cap_hex is not None
                # The accepted record embeds this successor cursor, so the
                # cursor binds the already-known predecessor checkpoint tip;
                # binding its own accepted-record hash would be cyclic.
                or _cursor(after.cursor).payload["journal_tip_sha256"]
                != previous.journal_tip_sha256):
            raise HLT16CampaignSchemaError("accepted retry/journal reset differs")
        if (new.accepted_macro_step_count != old.accepted_macro_step_count + 1
                or new.last_accepted_time.hex() != _cursor(after.cursor).payload["accepted_boundary_time"]["binary64_hex"]
                or new.current_macro_step_temporal_retry_count != 0
                or new.cumulative_temporal_retry_count != old.cumulative_temporal_retry_count
                or new.last_accepted_macro_step_temporal_retry_count != old.current_macro_step_temporal_retry_count
                or new.serialized_temporal_rejections != old.serialized_temporal_rejections
                or any(item < prior for item, prior in zip(new.accumulated_debit_vector, old.accumulated_debit_vector))):
            raise HLT16CampaignSchemaError("accepted ledger differs")
        old_cursor, new_cursor = _cursor(before.cursor), _cursor(after.cursor)
        if (new_cursor.mode != "FRESH_READY" or new_cursor.payload["cursor_generation"] != old_cursor.payload["cursor_generation"] + 1
                or new_cursor.payload["attempt_serial"] != old_cursor.payload["attempt_serial"] + 1
                or new_cursor.payload["cursor_chain_parent_sha256"] != old_cursor.sha256
                or new_cursor.payload["previous_step_index"] <= old_cursor.payload["previous_step_index"]
                or new_cursor.payload["previous_transaction_serial"] <= old_cursor.payload["previous_transaction_serial"]):
            raise HLT16CampaignSchemaError("accepted cursor lifecycle differs")
        return

    if kinds == ("tdg6_rejection", "cursor_transition"):
        rejection, transition = parsed; key = rejection["payload"]["member_key"]
        before, after = previous.members[key], current.members[key]
        if (transition["payload"]["member_key"] != key or changed != (key,)
                or current.event != previous.event or current.target != previous.target
                or rejection["payload"]["predecessor_descriptor_sha256"] != before.descriptor_sha256
                or rejection["payload"]["predecessor_cursor_sha256"] != _cursor(before.cursor).sha256
                or transition["payload"]["rejection_record_sha256"] != rejection["record_sha256"]
                or transition["payload"]["predecessor_cursor_sha256"] != _cursor(before.cursor).sha256
                or transition["payload"]["successor_cursor_sha256"] != _cursor(after.cursor).sha256
                or after.descriptor_sha256 != before.descriptor_sha256
                or current.disposition != "nonterminal"
                or after.source_current != before.source_current
                or after.source_total != before.source_total
                or after.cfl_current != before.cfl_current
                or after.cfl_total != before.cfl_total):
            raise HLT16CampaignSchemaError("TDG6 transition differs")
        try:
            retry_evidence = rejection["payload"]["evidence"]
            lifecycle.validate_retry_payload(
                _cursor(after.cursor),
                _ledger(after.ledger),
                evolution_state_sha256=_retry_evolution_state_sha256(retry_evidence),
            )
        except Exception as error:
            raise HLT16CampaignSchemaError("pending cursor differs") from error
        old_ledger, new_ledger = _ledger(before.ledger), _ledger(after.ledger)
        try:
            evidence = lifecycle.validate_retry_evidence(
                _cursor(before.cursor),
                old_ledger,
                retry_evidence,
                evolution_state_sha256=_retry_evolution_state_sha256(retry_evidence),
            )
        except Exception as error:
            raise HLT16CampaignSchemaError("TDG6 rejection evidence differs") from error
        if p15._ledger_mapping(p15._updated_retry_ledger(old_ledger, evidence)) != p15._ledger_mapping(new_ledger):
            raise HLT16CampaignSchemaError("TDG6 updated ledger differs")
        if _cursor(after.cursor).payload["journal_tip_sha256"] != rejection["record_sha256"]:
            raise HLT16CampaignSchemaError("pending cursor journal differs")
        return

    if kinds in (("source_rejection",), ("cfl_rejection",)):
        p = parsed[0]["payload"]; key = p["member_key"]; before, after = previous.members[key], current.members[key]
        if (changed != (key,) or current.event != previous.event or current.target != previous.target
                or current.disposition != "nonterminal" or p["exhausted"] is not False
                or after.descriptor_sha256 != before.descriptor_sha256
                or _cursor(after.cursor).sha256 != _cursor(before.cursor).sha256
                or after.ledger != before.ledger
                or p["predecessor_descriptor_sha256"] != before.descriptor_sha256
                or p["predecessor_cursor_sha256"] != _cursor(before.cursor).sha256):
            raise HLT16CampaignSchemaError("retry predecessor differs")
        actual = (after.source_current, after.source_total) if kinds[0] == "source_rejection" else (after.cfl_current, after.cfl_total)
        if actual != (p["retry_current"], p["retry_total"]):
            raise HLT16CampaignSchemaError("retry counter differs")
        expected_owner = "source" if kinds[0] == "source_rejection" else "cfl"
        old_current, old_total = ((before.source_current, before.source_total)
                                  if expected_owner == "source" else (before.cfl_current, before.cfl_total))
        other_before = ((before.cfl_current, before.cfl_total)
                        if expected_owner == "source" else (before.source_current, before.source_total))
        other_after = ((after.cfl_current, after.cfl_total)
                       if expected_owner == "source" else (after.source_current, after.source_total))
        if (actual != (old_current + 1, old_total + 1) or other_after != other_before
                or after.pending_owner != expected_owner or after.pending_cap_hex != p["next_cap_hex"]):
            raise HLT16CampaignSchemaError("retry pending owner/cap differs")
        return

    if kinds in (("terminal_lock",), ("source_rejection", "terminal_lock"),
                 ("cfl_rejection", "terminal_lock"), ("tdg6_rejection", "terminal_lock")):
        lock = parsed[-1]; p = lock["payload"]; member = previous.members[p["member_key"]]
        expected_classification = "scientific" if current.disposition == "scientific_terminal" else "invalid"
        if (current.disposition not in TERMINAL_DISPOSITIONS or changed != ()
                or current.terminal is None or current.terminal["classification"] != expected_classification
                or current.terminal["member_key"] != p["member_key"]
                or current.terminal["reason"] != p["reason"]
                or p["classification"] != expected_classification
                or current.event != previous.event or current.target != previous.target):
            raise HLT16CampaignSchemaError("terminal transition differs")
        if p["descriptor_sha256"] != member.descriptor_sha256 or p["cursor_sha256"] != _cursor(member.cursor).sha256:
            raise HLT16CampaignSchemaError("terminal identity differs")
        if kinds == ("terminal_lock",):
            expected_owner = "scientific" if expected_classification == "scientific" else "invalid"
            if p["owner"] != expected_owner or p["rejection_record_sha256"] is not None:
                raise HLT16CampaignSchemaError("terminal owner differs")
        else:
            rejection = parsed[0]
            if (p["rejection_record_sha256"] != rejection["record_sha256"]
                    or p["member_key"] != rejection["payload"]["member_key"]):
                raise HLT16CampaignSchemaError("terminal rejection linkage differs")
        if kinds in (("source_rejection", "terminal_lock"), ("cfl_rejection", "terminal_lock")):
            rejection = parsed[0]["payload"]
            owner = "source" if kinds[0] == "source_rejection" else "cfl"
            if (p["owner"] != owner or p["evidence"] != rejection["evidence"]
                    or expected_classification != "invalid" or rejection["exhausted"] is not True
                    or rejection["next_cap_hex"] is not None):
                raise HLT16CampaignSchemaError("terminal owner exhaustion differs")
            old_pair = ((member.source_current, member.source_total) if owner == "source"
                        else (member.cfl_current, member.cfl_total))
            if (rejection["retry_current"], rejection["retry_total"]) != (old_pair[0] + 1, old_pair[1] + 1):
                raise HLT16CampaignSchemaError("terminal exhaustion counters differ")
        if kinds == ("tdg6_rejection", "terminal_lock"):
            rejection = parsed[0]["payload"]
            if (p["owner"] != "temporal" or p["evidence"] != rejection["evidence"]
                    or expected_classification != "invalid"):
                raise HLT16CampaignSchemaError("terminal TDG6 exhaustion differs")
            try:
                evidence = lifecycle.validate_retry_evidence(
                    _cursor(member.cursor),
                    _ledger(member.ledger),
                    rejection["evidence"],
                    evolution_state_sha256=_retry_evolution_state_sha256(
                        rejection["evidence"]
                    ),
                )
            except Exception as error:
                raise HLT16CampaignSchemaError("terminal TDG6 evidence differs") from error
            exhausted = (
                int(evidence["retry_count_for_current_macro_step"])
                > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
                or float(evidence["retry_step_size"])
                < float(p15.TDG6_MINIMUM_MACRO_STEP)
            )
            if not exhausted or p["reason"] != "temporal_retry_exhausted":
                raise HLT16CampaignSchemaError("terminal TDG6 exhaustion differs")
            cursor = _cursor(member.cursor)
            if cursor.mode == "RETRY_PENDING":
                try:
                    pending_evidence = cursor.payload[
                        "retry_successor_payload_or_none"
                    ]["TDG6_rejection_evidence"]
                    lifecycle.validate_retry_payload(
                        cursor,
                        _ledger(member.ledger),
                        evolution_state_sha256=_retry_evolution_state_sha256(
                            pending_evidence
                        ),
                    )
                    pending = cursor.payload["retry_successor_payload_or_none"]
                    assert isinstance(pending, Mapping)
                    successor = p15._plan_mapping(pending["successor_plan"])
                except Exception as error:
                    raise HLT16CampaignSchemaError(
                        "terminal pending TDG6 predecessor differs"
                    ) from error
                if (
                    not p15._same_binary64(
                        evidence["initial_time"],
                        float.fromhex(str(successor["current_hex"])),
                    )
                    or not p15._same_binary64(
                        evidence["attempted_macro_step_size"],
                        float.fromhex(str(successor["macro_width_hex"])),
                    )
                ):
                    raise HLT16CampaignSchemaError(
                        "terminal pending TDG6 successor differs"
                    )
            # The deterministic rejected ledger is reconstructed here but is
            # deliberately not installed as the restartable accepted member.
            p15._updated_retry_ledger(_ledger(member.ledger), evidence)
        return

    if kinds == ("common_event_commit",):
        p = parsed[0]["payload"]
        if (current.event != previous.event + 1 or p["predecessor_event"] != previous.event
                or p["completed_event"] != current.event or p["completed_target"] != previous.target
                or p["next_target"] != current.target or p["predecessor_checkpoint_sha256"] != previous.sha256
                or current.disposition != "event_complete" or changed != MEMBER_KEYS):
            raise HLT16CampaignSchemaError("common event successor differs")
        for key in MEMBER_KEYS:
            before, after = previous.members[key], current.members[key]
            old_cursor, new_cursor = _cursor(before.cursor), _cursor(after.cursor)
            if (old_cursor.mode != "FRESH_READY" or before.pending_owner is not None
                    or before.pending_cap_hex is not None
                    or old_cursor.payload["accepted_boundary_time"] != previous.target
                    or after.descriptor_sha256 != before.descriptor_sha256
                    or after.ledger != before.ledger
                    or after.source_current != before.source_current
                    or after.source_total != before.source_total
                    or after.cfl_current != before.cfl_current
                    or after.cfl_total != before.cfl_total
                    or after.pending_owner is not None or after.pending_cap_hex is not None
                    or new_cursor.mode != "FRESH_READY"
                    or new_cursor.payload["accepted_boundary_time"] != previous.target
                    or new_cursor.payload["committed_common_event_index"] != current.event
                    or new_cursor.payload["event_target_time"] != current.target
                    or new_cursor.payload["attempt_serial"] != 0
                    or new_cursor.payload["cursor_generation"] != old_cursor.payload["cursor_generation"] + 1
                    or new_cursor.payload["journal_tip_sha256"] != parsed[0]["record_sha256"]
                    or new_cursor.payload["cursor_chain_parent_sha256"] != old_cursor.sha256):
                raise HLT16CampaignSchemaError("common event member successor differs")
        return
    raise HLT16CampaignSchemaError("unsupported checkpoint transition")


def validate_generation_one_bridge(checkpoint: CampaignCheckpoint, record: Mapping[str, Any],
                                   *, external_checkpoint_sha256: str) -> None:
    """Bind one generation-one HLT16 checkpoint to immutable PROTO17 GEN0."""
    checkpoint.validate()
    parsed = validate_journal(record)
    if (checkpoint.generation != 1 or parsed["kind"] != "genesis_bridge" or parsed["sequence"] != 0
            or parsed["previous_record_sha256"] != ROOT_SHA256 or checkpoint.journal_tip_sha256 != parsed["record_sha256"]
            or checkpoint.journal_sequence != 0):
        raise HLT16CampaignSchemaError("generation-one bridge differs")
    p = parsed["payload"]
    if (p["external_checkpoint_sha256"] != _sha(external_checkpoint_sha256, "external checkpoint")
            or checkpoint.parent_sha256 != external_checkpoint_sha256 or p["plan_sha256"] != checkpoint.plan_sha256):
        raise HLT16CampaignSchemaError("bridge authority differs")
    for key in MEMBER_KEYS:
        if p["member_descriptors"][key] != checkpoint.members[key].descriptor_sha256 or p["bridge_cursors"][key] != dict(checkpoint.members[key].cursor):
            raise HLT16CampaignSchemaError("bridge member differs")


__all__ = [
    "CampaignCheckpoint", "HLT16CampaignSchemaError", "JOURNAL_SCHEMA", "MemberCampaignState",
    "NONTERMINAL_DISPOSITIONS", "PROTOCOL", "RECORD_KINDS", "TERMINAL_DISPOSITIONS",
    "journal_envelope", "validate_checkpoint_successor", "validate_generation_one_bridge", "validate_journal",
]
