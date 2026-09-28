from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import proto15_runtime as p15
from recursive_horizons.fgc.evolution.hlt16_campaign_schema import (
    CampaignCheckpoint, HLT16CampaignSchemaError, MemberCampaignState, PROTOCOL,
    journal_envelope, validate_checkpoint_successor, validate_generation_one_bridge,
    validate_journal,
)
from recursive_horizons.fgc.evolution.hlt16_lifecycle import (
    close_fine, close_retry as _close_retry, retry_exhaustion as _retry_exhaustion,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import MEMBER_KEYS, ROOT_SHA256

SHA = "a" * 64
COMMIT = "b" * 40
CAMPAIGN = "FGC-2-SF1-PROTO17-GR0-A3-CAL-test"
TARGET = {"rational": "3/2", "binary64_hex": float(1.5).hex()}
START = {"rational": "23/16", "binary64_hex": float(23.0 / 16.0).hex()}
EVOLUTION_STATE = "9" * 64


def close_retry(*args, **kwargs):
    kwargs["evolution_state_sha256"] = EVOLUTION_STATE
    return _close_retry(*args, **kwargs)


def retry_exhaustion(*args, **kwargs):
    kwargs["evolution_state_sha256"] = EVOLUTION_STATE
    return _retry_exhaustion(*args, **kwargs)


def ledger() -> p15.TDG6TemporalLedger:
    return p15.TDG6TemporalLedger(
        last_accepted_time=23.0 / 16.0,
        accepted_macro_step_count=0,
        cumulative_temporal_retry_count=0,
        current_macro_step_temporal_retry_count=0,
        last_accepted_macro_step_temporal_retry_count=0,
        accumulated_debit_vector=(0.0,) * 18,
        serialized_temporal_rejections=(),
    )


def cursor(key: str, descriptor: str, value: p15.TDG6TemporalLedger | None = None) -> dict[str, object]:
    value = ledger() if value is None else value
    method, points = key.split("-", 1)
    raw = {
        "protocol_artifact_id": PROTOCOL, "campaign_id": CAMPAIGN, "member_key": key,
        "method": method, "point_count": int(points), "committed_common_event_index": 23,
        "accepted_boundary_time": START, "accepted_state_sha256": descriptor,
        "previous_step_index": 9, "previous_transaction_serial": 13,
        "TDG6_ledger_sha256": p15._hash(p15._ledger_mapping(value)),
        "cursor_generation": 0, "event_target_time": TARGET, "attempt_serial": 0,
        "attempt_id": f"{CAMPAIGN}:{key}:23:0", "mode": "FRESH_READY",
        "retry_successor_payload_or_none": None, "journal_tip_sha256": ROOT_SHA256,
        "cursor_chain_parent_sha256": "c" * 64,
    }
    return dict(p15._cursor(raw).payload)


def member(key: str, descriptor: str) -> MemberCampaignState:
    value = ledger()
    return MemberCampaignState(descriptor, cursor(key, descriptor, value), p15._ledger_mapping(value))


def bridge_checkpoint() -> tuple[CampaignCheckpoint, dict[str, object]]:
    descriptors = {key: f"{index:x}" * 64 for index, key in enumerate(MEMBER_KEYS, 1)}
    members = {key: member(key, descriptors[key]) for key in MEMBER_KEYS}
    record = journal_envelope(
        sequence=0, campaign_id=CAMPAIGN, generation=1, event=23,
        previous_record_sha256=ROOT_SHA256, kind="genesis_bridge",
        payload={
            "external_protocol": "FGC-2-SF1-PROTO17",
            "external_checkpoint_sha256": "d" * 64, "plan_sha256": SHA,
            "member_descriptors": descriptors,
            "bridge_cursors": {key: dict(members[key].cursor) for key in MEMBER_KEYS},
        },
    )
    checkpoint = CampaignCheckpoint(
        plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
        campaign_id=CAMPAIGN, event=23, target=TARGET, members=members,
        journal_tip_sha256=record["record_sha256"], parent_sha256="d" * 64,
        generation=1, journal_sequence=0,
    )
    return checkpoint, record


def event_plan(value: p15.Proto15Cursor, cap: float = 1 / 16) -> dict[str, object]:
    return p15._plan_mapping(p15.plan_forward_proto14_subdivision(
        float.fromhex(str(value.payload["accepted_boundary_time"]["binary64_hex"])),
        float.fromhex(str(value.payload["event_target_time"]["binary64_hex"])), cap,
        minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
    ))


def accepted_ledger(prior: p15.TDG6TemporalLedger, endpoint: float) -> p15.TDG6TemporalLedger:
    return p15.TDG6TemporalLedger(
        last_accepted_time=endpoint,
        accepted_macro_step_count=prior.accepted_macro_step_count + 1,
        cumulative_temporal_retry_count=prior.cumulative_temporal_retry_count,
        current_macro_step_temporal_retry_count=0,
        last_accepted_macro_step_temporal_retry_count=prior.current_macro_step_temporal_retry_count,
        accumulated_debit_vector=prior.accumulated_debit_vector,
        serialized_temporal_rejections=prior.serialized_temporal_rejections,
    )


def temporal_outcome(
    cursor_value: p15.Proto15Cursor,
    prior: p15.TDG6TemporalLedger,
    plan: dict[str, object],
    *,
    exhausted: bool = False,
) -> p15.TDG6TemporalRetryRequired | p15.TDG6TemporalRetryExhausted:
    admission = p15._synthetic_failed_admission(
        p15.RUNTIME_METHOD_BY_LABEL[cursor_value.payload["method"]]
    )
    count = prior.current_macro_step_temporal_retry_count + 1
    evidence = p15.TDG6TemporalRetryEvidence(
        event_type="rejected_TDG6_temporal_admission",
        method=p15.RUNTIME_METHOD_BY_LABEL[cursor_value.payload["method"]],
        initial_time=float.fromhex(str(plan["current_hex"])),
        attempted_macro_step_size=float.fromhex(str(plan["macro_width_hex"])),
        retry_step_size=float.fromhex(str(plan["macro_width_hex"])) * 0.5,
        initial_state_sha256=EVOLUTION_STATE,
        previous_step_index=int(cursor_value.payload["previous_step_index"]),
        previous_transaction_serial=int(cursor_value.payload["previous_transaction_serial"]),
        retry_count_for_current_macro_step=count,
        cumulative_temporal_retry_count=prior.cumulative_temporal_retry_count + 1,
        failed_channels=admission.failed_channels,
        channel_admissions=admission.channel_admissions,
    )
    updated = p15._updated_retry_ledger(prior, evidence.as_mapping())
    if exhausted:
        reason = (
            "maximum_temporal_retries"
            if count > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
            else "minimum_macro_step"
        )
        return p15.TDG6TemporalRetryExhausted(
            reason=reason, evidence=evidence, updated_ledger=updated,
        )
    return p15.TDG6TemporalRetryRequired(evidence, updated)


def accepted_successor(
    previous: CampaignCheckpoint,
    key: str,
    descriptor: str,
) -> tuple[CampaignCheckpoint, dict[str, object]]:
    """Build a one-member accepted edge exclusively through HLT16 lifecycle."""
    before = previous.members[key]
    cursor_value = p15.Proto15Cursor(dict(before.cursor))
    prior = p15._ledger_from_mapping(before.ledger)
    plan = event_plan(cursor_value)
    endpoint = float.fromhex(str(plan["boundaries_hex"][-1]))
    evolved_ledger = accepted_ledger(prior, endpoint)
    successor = close_fine(
        cursor_value, prior, evolved_ledger, descriptor, plan,
        {"rational": str(Fraction.from_float(endpoint)), "binary64_hex": endpoint.hex()},
        int(cursor_value.payload["previous_step_index"]) + 1,
        int(cursor_value.payload["previous_transaction_serial"]) + 1,
        previous.journal_tip_sha256,
    )
    states = dict(previous.members)
    states[key] = MemberCampaignState(descriptor, dict(successor.payload), p15._ledger_mapping(evolved_ledger))
    record = journal_envelope(
        sequence=previous.journal_sequence + 1, campaign_id=CAMPAIGN,
        generation=previous.generation + 1, event=previous.event,
        previous_record_sha256=previous.journal_tip_sha256, kind="accepted_fine",
        payload={
            "member_key": key,
            "predecessor_cursor_sha256": cursor_value.sha256,
            "successor_cursor": dict(successor.payload),
            "successor_cursor_sha256": successor.sha256,
            "descriptor_sha256": descriptor,
            "TDG6_ledger_sha256": successor.payload["TDG6_ledger_sha256"],
            "executed_plan": plan,
        },
    )
    current = CampaignCheckpoint(
        plan_sha256=previous.plan_sha256, authorization_commit=previous.authorization_commit,
        protocol=PROTOCOL, campaign_id=CAMPAIGN, event=previous.event,
        target=previous.target, members=states, journal_tip_sha256=record["record_sha256"],
        parent_sha256=previous.sha256, generation=previous.generation + 1,
        journal_sequence=record["sequence"],
    )
    validate_checkpoint_successor(previous, current, [record])
    return current, record


class HLT16CampaignSchemaTests(unittest.TestCase):
    def test_generation_one_bridge_binds_all_six_members(self) -> None:
        checkpoint, record = bridge_checkpoint()
        validate_generation_one_bridge(checkpoint, record, external_checkpoint_sha256="d" * 64)
        self.assertEqual(CampaignCheckpoint.from_mapping(checkpoint.mapping()).sha256, checkpoint.sha256)

    def test_canonical_journal_is_independent_of_mapping_insertion_order(self) -> None:
        checkpoint, record = bridge_checkpoint()
        payload = dict(record["payload"])
        reordered = {name: payload[name] for name in reversed(tuple(payload))}
        same = journal_envelope(
            sequence=0, campaign_id=CAMPAIGN, generation=1, event=23,
            previous_record_sha256=ROOT_SHA256, kind="genesis_bridge", payload=reordered,
        )
        self.assertEqual(record, same)
        self.assertEqual(validate_journal(same)["record_sha256"], record["record_sha256"])
        self.assertEqual(checkpoint.members["RK4-2049"].cursor["protocol_artifact_id"], PROTOCOL)

    def test_tdg6_record_boundary_normalizes_native_tuple_evidence(self) -> None:
        previous, _bridge = bridge_checkpoint()
        key = "RK4-2049"
        state = previous.members[key]
        cursor_value = p15.Proto15Cursor(dict(state.cursor))
        outcome = temporal_outcome(
            cursor_value,
            p15._ledger_from_mapping(state.ledger),
            event_plan(cursor_value),
        )
        assert isinstance(outcome, p15.TDG6TemporalRetryRequired)
        native = outcome.evidence.as_mapping()
        self.assertIsInstance(native["failed_channels"], tuple)
        self.assertIsInstance(native["channel_admissions"], tuple)

        record = journal_envelope(
            sequence=previous.journal_sequence + 1,
            campaign_id=CAMPAIGN,
            generation=previous.generation + 1,
            event=previous.event,
            previous_record_sha256=previous.journal_tip_sha256,
            kind="tdg6_rejection",
            payload={
                "member_key": key,
                "predecessor_descriptor_sha256": state.descriptor_sha256,
                "predecessor_cursor_sha256": cursor_value.sha256,
                "evidence": native,
            },
        )
        durable = record["payload"]["evidence"]
        self.assertIsInstance(durable["failed_channels"], list)
        self.assertIsInstance(durable["channel_admissions"], list)
        self.assertEqual(durable, p15._json_safe(native))
        self.assertEqual(validate_journal(record), record)

        malformed = dict(native)
        malformed["unsupported"] = {"not", "canonical", "JSON"}
        with self.assertRaisesRegex(HLT16CampaignSchemaError, "TDG6 evidence differs"):
            journal_envelope(
                sequence=previous.journal_sequence + 1,
                campaign_id=CAMPAIGN,
                generation=previous.generation + 1,
                event=previous.event,
                previous_record_sha256=previous.journal_tip_sha256,
                kind="tdg6_rejection",
                payload={
                    "member_key": key,
                    "predecessor_descriptor_sha256": state.descriptor_sha256,
                    "predecessor_cursor_sha256": cursor_value.sha256,
                    "evidence": malformed,
                },
            )

    def test_source_retry_successor_changes_one_member_and_binds_journal(self) -> None:
        previous, _bridge = bridge_checkpoint()
        key = "RK4-2049"
        states = dict(previous.members)
        old = states[key]
        states[key] = MemberCampaignState(
            old.descriptor_sha256, old.cursor, old.ledger,
            source_current=1, source_total=1, cfl_current=0, cfl_total=0,
            pending_owner="source", pending_cap_hex="0x1.0000000000000p-5",
        )
        record = journal_envelope(
            sequence=1, campaign_id=CAMPAIGN, generation=2, event=23,
            previous_record_sha256=previous.journal_tip_sha256, kind="source_rejection",
            payload={
                "member_key": key, "evidence": {"owner": "source", "residual": "0x1.0p+0"}, "no_commit": True,
                "predecessor_descriptor_sha256": old.descriptor_sha256,
                "predecessor_cursor_sha256": old.cursor["cursor_chain_sha256"],
                "retry_current": 1, "retry_total": 1, "retry_limit": 32,
                "attempted_cap_hex": "0x1.0000000000000p-4",
                "next_cap_hex": "0x1.0000000000000p-5", "exhausted": False,
            },
        )
        current = CampaignCheckpoint(
            plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
            campaign_id=CAMPAIGN, event=23, target=TARGET, members=states,
            journal_tip_sha256=record["record_sha256"], parent_sha256=previous.sha256,
            generation=2, journal_sequence=1,
        )
        validate_checkpoint_successor(previous, current, [record])

    def test_lifecycle_accepted_fine_binds_record_to_current_checkpoint_tip(self) -> None:
        previous, _bridge = bridge_checkpoint()
        current, record = accepted_successor(previous, "RK4-2049", "e" * 64)
        self.assertEqual(record["previous_record_sha256"], previous.journal_tip_sha256)
        self.assertEqual(current.members["RK4-2049"].cursor["accepted_boundary_time"], TARGET)
        # The cursor binds the current *predecessor checkpoint* tip, not its
        # own accepted-record hash (which would create a self-hash cycle).
        self.assertEqual(
            current.members["RK4-2049"].cursor["journal_tip_sha256"],
            previous.journal_tip_sha256,
        )
        mutated = deepcopy(record)
        mutated["previous_record_sha256"] = ROOT_SHA256
        with self.assertRaises(HLT16CampaignSchemaError):
            validate_checkpoint_successor(previous, current, [mutated])

    def test_lifecycle_tdg6_retry_pending_fixture_and_lineage_mutation(self) -> None:
        previous, _bridge = bridge_checkpoint()
        key = "RK4-2049"
        state = previous.members[key]
        before = p15.Proto15Cursor(dict(state.cursor))
        prior = p15._ledger_from_mapping(state.ledger)
        old_plan = event_plan(before, 1 / 16)
        outcome = temporal_outcome(before, prior, old_plan)
        assert isinstance(outcome, p15.TDG6TemporalRetryRequired)
        successor_plan = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            outcome.evidence.initial_time, 1.5, outcome.retry_step_size,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ))
        rejection = journal_envelope(
            sequence=previous.journal_sequence + 1, campaign_id=CAMPAIGN,
            generation=previous.generation + 1, event=23,
            previous_record_sha256=previous.journal_tip_sha256, kind="tdg6_rejection",
            payload={
                "member_key": key,
                "predecessor_descriptor_sha256": state.descriptor_sha256,
                "predecessor_cursor_sha256": before.sha256,
                "evidence": p15._json_safe(outcome.evidence.as_mapping()),
            },
        )
        pending = close_retry(
            before, prior, outcome, old_plan, successor_plan,
            rejection["record_sha256"], previous.journal_tip_sha256,
        )
        states = dict(previous.members)
        states[key] = MemberCampaignState(
            state.descriptor_sha256, dict(pending.payload),
            p15._ledger_mapping(outcome.updated_ledger), pending_owner="temporal",
            pending_cap_hex=outcome.retry_step_size.hex(),
        )
        transition = journal_envelope(
            sequence=previous.journal_sequence + 2, campaign_id=CAMPAIGN,
            generation=previous.generation + 1, event=23,
            previous_record_sha256=rejection["record_sha256"], kind="cursor_transition",
            payload={
                "member_key": key, "predecessor_cursor_sha256": before.sha256,
                "rejection_record_sha256": rejection["record_sha256"],
                "successor_cursor": dict(pending.payload),
                "successor_cursor_sha256": pending.sha256,
            },
        )
        current = CampaignCheckpoint(
            plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
            campaign_id=CAMPAIGN, event=23, target=TARGET, members=states,
            journal_tip_sha256=transition["record_sha256"], parent_sha256=previous.sha256,
            generation=2, journal_sequence=transition["sequence"],
        )
        validate_checkpoint_successor(previous, current, [rejection, transition])
        bad_transition = deepcopy(transition)
        bad_transition["payload"]["predecessor_cursor_sha256"] = "f" * 64
        with self.assertRaises(HLT16CampaignSchemaError):
            validate_checkpoint_successor(previous, current, [rejection, bad_transition])

    def test_source_and_cfl_terminal_exhaustion_fixtures(self) -> None:
        for kind, owner in (("source_rejection", "source"), ("cfl_rejection", "cfl")):
            with self.subTest(owner=owner):
                previous, _bridge = bridge_checkpoint()
                key = "RK4-2049"; state = previous.members[key]
                exhausted_state = replace(
                    state,
                    source_current=32 if owner == "source" else 0,
                    source_total=32 if owner == "source" else 0,
                    cfl_current=32 if owner == "cfl" else 0,
                    cfl_total=32 if owner == "cfl" else 0,
                    pending_owner=owner,
                    pending_cap_hex="0x1.0000000000000p-36",
                )
                members = dict(previous.members); members[key] = exhausted_state
                # This is a valid, already-pending predecessor snapshot; the
                # terminal suffix itself has no restart-state mutation.
                prepared = replace(
                    previous, members=members, generation=2,
                    journal_sequence=1, parent_sha256="e" * 64,
                    journal_tip_sha256="f" * 64,
                )
                prepared.validate()
                rejection = journal_envelope(
                    sequence=2, campaign_id=CAMPAIGN, generation=3, event=23,
                    previous_record_sha256=prepared.journal_tip_sha256, kind=kind,
                    payload={
                        "member_key": key, "evidence": {"owner": owner, "exact": True},
                        "no_commit": True,
                        "predecessor_descriptor_sha256": state.descriptor_sha256,
                        "predecessor_cursor_sha256": state.cursor["cursor_chain_sha256"],
                        "retry_current": 33, "retry_total": 33, "retry_limit": 32,
                        "attempted_cap_hex": "0x1.0000000000000p-35",
                        "next_cap_hex": None, "exhausted": True,
                    },
                )
                lock = journal_envelope(
                    sequence=3, campaign_id=CAMPAIGN, generation=3, event=23,
                    previous_record_sha256=rejection["record_sha256"], kind="terminal_lock",
                    payload={
                        "member_key": key, "reason": f"{owner}_retry_exhausted",
                        "classification": "invalid", "descriptor_sha256": state.descriptor_sha256,
                        "cursor_sha256": state.cursor["cursor_chain_sha256"], "owner": owner,
                        "rejection_record_sha256": rejection["record_sha256"],
                        "evidence": {"owner": owner, "exact": True},
                    },
                )
                terminal = CampaignCheckpoint(
                    plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
                    campaign_id=CAMPAIGN, event=23, target=TARGET, members=prepared.members,
                    journal_tip_sha256=lock["record_sha256"], parent_sha256=prepared.sha256,
                    generation=3, journal_sequence=3, disposition="invalid_terminal",
                    terminal={"classification": "invalid", "event": 23, "target": TARGET,
                              "member_key": key, "reason": f"{owner}_retry_exhausted", "owner": owner},
                )
                validate_checkpoint_successor(prepared, terminal, [rejection, lock])

    def test_pending_tdg6_terminal_exhaustion_has_exact_successor_lineage(self) -> None:
        previous, _bridge = bridge_checkpoint(); key = "RK4-2049"
        state = previous.members[key]; before = p15.Proto15Cursor(dict(state.cursor))
        prior = p15._ledger_from_mapping(state.ledger)
        old_plan = event_plan(before, 2.0 * float(p15.TDG6_MINIMUM_MACRO_STEP))
        first = temporal_outcome(before, prior, old_plan)
        assert isinstance(first, p15.TDG6TemporalRetryRequired)
        successor_plan = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            first.evidence.initial_time, 1.5, first.retry_step_size,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ))
        first_record = journal_envelope(
            sequence=1, campaign_id=CAMPAIGN, generation=2, event=23,
            previous_record_sha256=previous.journal_tip_sha256, kind="tdg6_rejection",
            payload={"member_key": key, "predecessor_descriptor_sha256": state.descriptor_sha256,
                     "predecessor_cursor_sha256": before.sha256,
                     "evidence": p15._json_safe(first.evidence.as_mapping())},
        )
        pending = close_retry(before, prior, first, old_plan, successor_plan,
                              first_record["record_sha256"], previous.journal_tip_sha256)
        states = dict(previous.members)
        states[key] = MemberCampaignState(state.descriptor_sha256, dict(pending.payload),
                                           p15._ledger_mapping(first.updated_ledger),
                                           pending_owner="temporal", pending_cap_hex=first.retry_step_size.hex())
        transition = journal_envelope(
            sequence=2, campaign_id=CAMPAIGN, generation=2, event=23,
            previous_record_sha256=first_record["record_sha256"], kind="cursor_transition",
            payload={"member_key": key, "predecessor_cursor_sha256": before.sha256,
                     "rejection_record_sha256": first_record["record_sha256"],
                     "successor_cursor": dict(pending.payload), "successor_cursor_sha256": pending.sha256},
        )
        pending_checkpoint = CampaignCheckpoint(
            plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
            campaign_id=CAMPAIGN, event=23, target=TARGET, members=states,
            journal_tip_sha256=transition["record_sha256"], parent_sha256=previous.sha256,
            generation=2, journal_sequence=2,
        )
        validate_checkpoint_successor(previous, pending_checkpoint, [first_record, transition])
        exhausted = temporal_outcome(pending, first.updated_ledger, successor_plan, exhausted=True)
        assert isinstance(exhausted, p15.TDG6TemporalRetryExhausted)
        # Exercise the lifecycle owner before presenting the schema suffix.
        retry_exhaustion(pending, first.updated_ledger, exhausted, "e" * 64,
                         pending_checkpoint.journal_tip_sha256)
        rejection = journal_envelope(
            sequence=3, campaign_id=CAMPAIGN, generation=3, event=23,
            previous_record_sha256=pending_checkpoint.journal_tip_sha256, kind="tdg6_rejection",
            payload={"member_key": key, "predecessor_descriptor_sha256": state.descriptor_sha256,
                     "predecessor_cursor_sha256": pending.sha256,
                     "evidence": p15._json_safe(exhausted.evidence.as_mapping())},
        )
        lock = journal_envelope(
            sequence=4, campaign_id=CAMPAIGN, generation=3, event=23,
            previous_record_sha256=rejection["record_sha256"], kind="terminal_lock",
            payload={"member_key": key, "reason": "temporal_retry_exhausted", "classification": "invalid",
                     "descriptor_sha256": state.descriptor_sha256, "cursor_sha256": pending.sha256,
                     "owner": "temporal", "rejection_record_sha256": rejection["record_sha256"],
                     "evidence": p15._json_safe(exhausted.evidence.as_mapping())},
        )
        terminal = CampaignCheckpoint(
            plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
            campaign_id=CAMPAIGN, event=23, target=TARGET, members=pending_checkpoint.members,
            journal_tip_sha256=lock["record_sha256"], parent_sha256=pending_checkpoint.sha256,
            generation=3, journal_sequence=4, disposition="invalid_terminal",
            terminal={"classification": "invalid", "event": 23, "target": TARGET,
                      "member_key": key, "reason": "temporal_retry_exhausted", "owner": "temporal"},
        )
        validate_checkpoint_successor(pending_checkpoint, terminal, [rejection, lock])
        bad = deepcopy(rejection)
        bad["payload"]["evidence"]["attempted_macro_step_size"] *= 2.0
        with self.assertRaises(HLT16CampaignSchemaError):
            validate_checkpoint_successor(pending_checkpoint, terminal, [bad, lock])

    def test_common_event_transition_after_six_lifecycle_accepted_fines(self) -> None:
        """Exercise the event edge only after all six real lifecycle commits."""
        previous, _bridge = bridge_checkpoint()
        for index, key in enumerate(MEMBER_KEYS, 1):
            previous, _record = accepted_successor(previous, key, f"{index + 8:x}" * 64)
        self.assertTrue(all(
            state.cursor["accepted_boundary_time"] == TARGET
            for state in previous.members.values()
        ))
        next_target = {"rational": "25/16", "binary64_hex": float(25 / 16).hex()}
        record = journal_envelope(
            sequence=previous.journal_sequence + 1, campaign_id=CAMPAIGN,
            generation=previous.generation + 1, event=previous.event,
            previous_record_sha256=previous.journal_tip_sha256, kind="common_event_commit",
            payload={
                "member_keys": list(MEMBER_KEYS), "predecessor_event": previous.event,
                "completed_event": previous.event + 1, "completed_target": TARGET,
                "next_target": next_target,
                "predecessor_checkpoint_sha256": previous.sha256,
            },
        )
        members: dict[str, MemberCampaignState] = {}
        for key in MEMBER_KEYS:
            state = previous.members[key]
            old = p15.Proto15Cursor(dict(state.cursor))
            payload = dict(old.payload)
            payload.update({
                "committed_common_event_index": previous.event + 1,
                "event_target_time": next_target,
                "cursor_generation": int(old.payload["cursor_generation"]) + 1,
                "attempt_serial": 0,
                "attempt_id": f"{CAMPAIGN}:{key}:{previous.event + 1}:0",
                "journal_tip_sha256": record["record_sha256"],
                "cursor_chain_parent_sha256": old.sha256,
            })
            payload.pop("cursor_chain_sha256", None)
            successor = p15._cursor(payload)
            members[key] = MemberCampaignState(
                state.descriptor_sha256, dict(successor.payload), state.ledger,
            )
        current = CampaignCheckpoint(
            plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
            campaign_id=CAMPAIGN, event=previous.event + 1, target=next_target,
            members=members, journal_tip_sha256=record["record_sha256"],
            parent_sha256=previous.sha256, generation=previous.generation + 1,
            journal_sequence=record["sequence"], disposition="event_complete",
        )
        validate_checkpoint_successor(previous, current, [record])
        broken_members = dict(members)
        broken = dict(broken_members["RK4-2049"].cursor)
        broken["event_target_time"] = TARGET
        broken_members["RK4-2049"] = MemberCampaignState(
            members["RK4-2049"].descriptor_sha256, p15._cursor(broken).payload,
            members["RK4-2049"].ledger,
        )
        with self.assertRaises(HLT16CampaignSchemaError):
            validate_checkpoint_successor(previous, replace(current, members=broken_members), [record])

    def test_record_and_checkpoint_mutations_fail_closed(self) -> None:
        previous, bridge = bridge_checkpoint()
        mutated = deepcopy(bridge)
        mutated["payload"]["external_protocol"] = PROTOCOL
        with self.assertRaises(HLT16CampaignSchemaError):
            validate_journal(mutated)
        bad_members = dict(previous.members)
        bad_members.pop("RK4-2049")
        with self.assertRaises(HLT16CampaignSchemaError):
            CampaignCheckpoint(
                plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
                campaign_id=CAMPAIGN, event=23, target=TARGET, members=bad_members,
                journal_tip_sha256=previous.journal_tip_sha256, parent_sha256="e" * 64,
                generation=1, journal_sequence=0,
            ).validate()
        old = previous.members["RK4-2049"]
        bad_cursor = dict(old.cursor)
        bad_cursor["cursor_chain_sha256"] = "0" * 64
        with self.assertRaises(HLT16CampaignSchemaError):
            MemberCampaignState(old.descriptor_sha256, bad_cursor, old.ledger).validate()

    def test_successor_rejects_changed_authority_or_wrong_record_parent(self) -> None:
        previous, _bridge = bridge_checkpoint()
        key = "RK4-2049"
        states = dict(previous.members)
        old = states[key]
        states[key] = MemberCampaignState(old.descriptor_sha256, old.cursor, old.ledger, source_current=1, source_total=1, pending_owner="source", pending_cap_hex="0x1.0000000000000p-5")
        record = journal_envelope(
            sequence=1, campaign_id=CAMPAIGN, generation=2, event=23,
            previous_record_sha256="f" * 64, kind="source_rejection",
            payload={
                "member_key": key, "evidence": {"owner": "source", "residual": "0x1.0p+0"}, "no_commit": True,
                "predecessor_descriptor_sha256": old.descriptor_sha256,
                "predecessor_cursor_sha256": old.cursor["cursor_chain_sha256"],
                "retry_current": 1, "retry_total": 1, "retry_limit": 32,
                "attempted_cap_hex": "0x1.0000000000000p-4",
                "next_cap_hex": "0x1.0000000000000p-5", "exhausted": False,
            },
        )
        current = CampaignCheckpoint(
            plan_sha256=SHA, authorization_commit=COMMIT, protocol=PROTOCOL,
            campaign_id=CAMPAIGN, event=23, target=TARGET, members=states,
            journal_tip_sha256=record["record_sha256"], parent_sha256=previous.sha256,
            generation=2, journal_sequence=1,
        )
        with self.assertRaises(HLT16CampaignSchemaError):
            validate_checkpoint_successor(previous, current, [record])

    def test_pending_owner_and_canonical_finite_cap_are_fail_closed(self) -> None:
        checkpoint, _bridge = bridge_checkpoint()
        old = checkpoint.members["RK4-2049"]
        with self.assertRaises(HLT16CampaignSchemaError):
            MemberCampaignState(
                old.descriptor_sha256, old.cursor, old.ledger,
                pending_owner="temporal", pending_cap_hex="0x1.0000000000000p-5",
            ).validate()
        with self.assertRaises(HLT16CampaignSchemaError):
            journal_envelope(
                sequence=1, campaign_id=CAMPAIGN, generation=2, event=23,
                previous_record_sha256=checkpoint.journal_tip_sha256, kind="source_rejection",
                payload={
                    "member_key": "RK4-2049", "evidence": {"owner": "source"}, "no_commit": True,
                    "predecessor_descriptor_sha256": old.descriptor_sha256,
                    "predecessor_cursor_sha256": old.cursor["cursor_chain_sha256"],
                    "retry_current": 1, "retry_total": 1, "retry_limit": 32,
                    "attempted_cap_hex": "0x1.0000000000000p-4",
                    "next_cap_hex": "inf", "exhausted": False,
                },
            )

    def test_singleton_invalid_terminal_has_an_exact_invalid_owner(self) -> None:
        previous, _bridge = bridge_checkpoint()
        key = "RK4-2049"
        state = previous.members[key]
        record = journal_envelope(
            sequence=1,
            campaign_id=CAMPAIGN,
            generation=2,
            event=23,
            previous_record_sha256=previous.journal_tip_sha256,
            kind="terminal_lock",
            payload={
                "member_key": key,
                "reason": "invalid_runtime_contract",
                "classification": "invalid",
                "descriptor_sha256": state.descriptor_sha256,
                "cursor_sha256": state.cursor["cursor_chain_sha256"],
                "owner": "invalid",
                "rejection_record_sha256": None,
                "evidence": {"exception": "typed-invalid"},
            },
        )
        terminal = {
            "classification": "invalid",
            "event": 23,
            "target": TARGET,
            "member_key": key,
            "reason": "invalid_runtime_contract",
            "owner": "invalid",
        }
        current = CampaignCheckpoint(
            plan_sha256=SHA,
            authorization_commit=COMMIT,
            protocol=PROTOCOL,
            campaign_id=CAMPAIGN,
            event=23,
            target=TARGET,
            members=previous.members,
            journal_tip_sha256=record["record_sha256"],
            parent_sha256=previous.sha256,
            generation=2,
            journal_sequence=1,
            disposition="invalid_terminal",
            terminal=terminal,
        )
        validate_checkpoint_successor(previous, current, [record])
        self.assertEqual(
            CampaignCheckpoint.from_mapping(current.mapping()).sha256,
            current.sha256,
        )
        bad = dict(terminal)
        bad["owner"] = "scientific"
        with self.assertRaises(HLT16CampaignSchemaError):
            replace(current, terminal=bad).validate()


if __name__ == "__main__":
    unittest.main()
