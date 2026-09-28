"""Behavioral synthetic controls for the PROTO15 durable campaign runtime."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution import proto15_runtime as runtime
from recursive_horizons.fgc.evolution.numerical_engine import PRIMARY_METHOD
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import TDG6_COMPLETE_STATE_CHANNELS, TDG6_MINIMUM_MACRO_STEP
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import TDG6TemporalLedger, TDG6TemporalRetryEvidence, TDG6TemporalRetryExhausted, TDG6TemporalRetryRequired, _magnitude_interval, classify_tdg6_runtime_intervals
from recursive_horizons.fgc.evolution.tdg7_binary64_subdivision_lattice import plan_forward_proto14_subdivision


def members():
    return {k: {"method": k.split("-", 1)[0], "point_count": int(k.split("-", 1)[1]), "state": {"member": k, "v": 0}, "ledger": TDG6TemporalLedger.zero(initial_time=23 / 16), "step_index": 1, "transaction_serial": 2} for k in runtime.MEMBER_KEYS}


def plan(cap=0.25):
    return plan_forward_proto14_subdivision(23 / 16, 3 / 2, cap, minimum_width=float(TDG6_MINIMUM_MACRO_STEP))


def failure_admission():
    outer = _magnitude_interval(np.ones((2, 4)), np.zeros((2, 4)), subinterval_count=2, owned_row_count=1)
    fine = _magnitude_interval(np.ones((4, 4)), np.zeros((4, 4)), subinterval_count=4, owned_row_count=1)
    return classify_tdg6_runtime_intervals(method=PRIMARY_METHOD, outer_intervals={c: outer for c in TDG6_COMPLETE_STATE_CHANNELS}, finest_intervals={c: fine for c in TDG6_COMPLETE_STATE_CHANNELS})


def outcome(campaign, member, prior_plan, *, exhausted=False):
    cursor, ledger, admission = campaign.cursors[member], campaign.ledgers[member], failure_admission()
    evidence = TDG6TemporalRetryEvidence(event_type="rejected_TDG6_temporal_admission", method=PRIMARY_METHOD, initial_time=23 / 16, attempted_macro_step_size=prior_plan.macro_width, retry_step_size=prior_plan.macro_width * .5, initial_state_sha256=campaign.states[member], previous_step_index=cursor.payload["previous_step_index"], previous_transaction_serial=cursor.payload["previous_transaction_serial"], retry_count_for_current_macro_step=ledger.current_macro_step_temporal_retry_count + 1, cumulative_temporal_retry_count=ledger.cumulative_temporal_retry_count + 1, failed_channels=admission.failed_channels, channel_admissions=admission.channel_admissions)
    updated = runtime._updated_retry_ledger(ledger, evidence.as_mapping())
    return TDG6TemporalRetryExhausted(reason="minimum_macro_step", evidence=evidence, updated_ledger=updated) if exhausted else TDG6TemporalRetryRequired(evidence, updated)


def record_at(root, sequence):
    path = next((root / "journal").glob(f"{sequence:020d}-*.journal"))
    return dict(runtime.Proto15SyntheticCampaign._read_record(path))


def rewrite_record(root, record):
    bare = {
        "record_kind": record["record_kind"],
        "journal_parent_sha256": record["journal_parent_sha256"],
        "payload": record["payload"],
    }
    rewritten = {**bare, "record_sha256": runtime._hash(bare)}
    encoded = runtime._canonical(rewritten)
    sequence = int(rewritten["payload"]["sequence"])
    for path in (root / "journal").glob(f"{sequence:020d}-*.journal"):
        path.unlink()
    path = root / "journal" / f"{sequence:020d}-{rewritten['record_sha256']}.journal"
    path.write_bytes(f"{len(encoded):08x} ".encode("ascii") + encoded + b"\n")
    return rewritten


def checkpoint_at(root, generation):
    path = next((root / "checkpoints").glob(f"{generation:020d}-*.json"))
    return json.loads(path.read_text())


def rewrite_checkpoint(root, checkpoint):
    cursor_set, ledger_set, state_set = runtime.Proto15SyntheticCampaign._sets(
        checkpoint["cursors"], checkpoint["ledgers"], checkpoint["states"]
    )
    checkpoint["cursor_set_sha256"] = cursor_set
    checkpoint["ledger_set_sha256"] = ledger_set
    checkpoint["accepted_state_set_sha256"] = state_set
    bare = dict(checkpoint)
    bare.pop("checkpoint_sha256", None)
    checkpoint["checkpoint_sha256"] = runtime._hash(bare)
    generation = int(checkpoint["campaign_generation"])
    for path in (root / "checkpoints").glob(f"{generation:020d}-*.json"):
        path.unlink()
    path = root / "checkpoints" / f"{generation:020d}-{checkpoint['checkpoint_sha256']}.json"
    path.write_bytes(runtime._canonical(checkpoint))
    return checkpoint


class Proto15RuntimeTests(unittest.TestCase):
    def campaign(self, root, **kwargs):
        return runtime.Proto15SyntheticCampaign.initialize(root, protocol_artifact_id="FGC-2-SF1-PROTO15", campaign_id="synthetic", members=members(), **kwargs)

    def reject(self, campaign, member="RK4-2049", prior_plan=None, exhausted=False):
        prior_plan = plan() if prior_plan is None else prior_plan
        typed = outcome(campaign, member, prior_plan, exhausted=exhausted)
        campaign.durable_rejection_sink(member, predecessor_plan=prior_plan)(typed.evidence.as_mapping())
        return typed, prior_plan

    def test_initialization_is_complete_content_addressed_and_recoverable(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root)
            self.assertEqual(tuple(campaign.cursors), runtime.MEMBER_KEYS)
            self.assertEqual(len(list((root / "states").glob("*.json"))), 6)
            self.assertEqual(campaign._checkpoint["member_keys"], list(runtime.MEMBER_KEYS))
            for key in runtime.MEMBER_KEYS:
                self.assertEqual(campaign.guard_fresh_entrypoint(key).payload["accepted_boundary_time"], {"rational": "23/16", "binary64_hex": (23 / 16).hex()})
            self.assertEqual(runtime.recover_campaign(root).checkpoint_sha256, campaign.checkpoint_sha256)

    def test_accepted_commit_advances_one_member_and_publishes_complete_checkpoint(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root); prior = dict(campaign.states)
            ledger = replace(campaign.ledgers["RK4-2049"], last_accepted_time=3 / 2, accepted_macro_step_count=1)
            cursor = campaign.accepted_fine_commit("RK4-2049", state={"member": "RK4-2049", "v": 1}, ledger=ledger, executed_plan=plan(), accepted_boundary_time=3 / 2, step_index=3, transaction_serial=4)
            self.assertEqual(cursor.mode, "FRESH_READY"); self.assertEqual(campaign.generation, 1)
            self.assertNotEqual(campaign.states["RK4-2049"], prior["RK4-2049"])
            self.assertTrue(all(campaign.states[k] == prior[k] for k in runtime.MEMBER_KEYS[1:]))
            self.assertEqual(runtime.recover_campaign(root).ledgers["RK4-2049"], ledger)

    def test_retry_required_is_durable_then_pending_and_exact_suffix_recovers(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root); typed, previous = self.reject(campaign)
            with self.assertRaises(runtime.Proto15JournalIntegrityError): campaign.accepted_fine_commit("RK4-4097", state={"x": 1}, ledger=campaign.ledgers["RK4-4097"], executed_plan=plan(), accepted_boundary_time=3 / 2, step_index=3, transaction_serial=4)
            recovered = runtime.recover_campaign(root)
            self.assertEqual(recovered.guard_retry_entrypoint("RK4-2049").mode, "RETRY_PENDING")
            self.assertEqual(recovered.ledgers["RK4-2049"], typed.updated_ledger)
        with TemporaryDirectory() as td:
            campaign = self.campaign(Path(td) / "campaign"); typed, previous = self.reject(campaign)
            pending = campaign.close_tdg6_outcome("RK4-2049", typed, retry_successor={"predecessor_plan": previous, "successor_plan": plan(typed.retry_step_size)})
            self.assertEqual(pending.mode, "RETRY_PENDING")
            with self.assertRaises(runtime.Proto15CursorContractError): campaign.guard_fresh_entrypoint("RK4-2049")
            wrong_ledger = replace(
                campaign.ledgers["RK4-2049"],
                last_accepted_time=3 / 2,
                accepted_macro_step_count=1,
                current_macro_step_temporal_retry_count=0,
                last_accepted_macro_step_temporal_retry_count=1,
            )
            with self.assertRaisesRegex(
                runtime.Proto15CursorContractError,
                "exact durable TDG7 successor plan",
            ):
                campaign.accepted_fine_commit(
                    "RK4-2049", state={"member": "RK4-2049", "v": 1},
                    ledger=wrong_ledger, executed_plan=plan(),
                    accepted_boundary_time=3 / 2, step_index=3,
                    transaction_serial=4,
                )
            successor = pending.payload["retry_successor_payload_or_none"]["successor_plan"]
            endpoint = float.fromhex(successor["boundaries_hex"][-1])
            accepted_ledger = replace(
                campaign.ledgers["RK4-2049"],
                last_accepted_time=endpoint,
                accepted_macro_step_count=1,
                current_macro_step_temporal_retry_count=0,
                last_accepted_macro_step_temporal_retry_count=1,
            )
            accepted = campaign.accepted_fine_commit(
                "RK4-2049", state={"member": "RK4-2049", "v": 1},
                ledger=accepted_ledger, executed_plan=successor,
                accepted_boundary_time=endpoint, step_index=3,
                transaction_serial=4,
            )
            self.assertEqual(accepted.mode, "FRESH_READY")
            self.assertEqual(runtime.recover_campaign(campaign.root).cursors["RK4-2049"].mode, "FRESH_READY")

    def test_exhaustion_is_terminal_and_locks_both_entrypoints(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root); tiny = plan(float(TDG6_MINIMUM_MACRO_STEP)); typed, _ = self.reject(campaign, prior_plan=tiny, exhausted=True)
            campaign.close_tdg6_outcome("RK4-2049", typed)
            self.assertTrue(campaign.terminal_lock); self.assertEqual(campaign._checkpoint["terminal_closure"]["classification"], "invalid_implementation_or_nonconverged_run")
            self.assertEqual(campaign.ledgers["RK4-2049"], typed.updated_ledger)
            with self.assertRaises(runtime.Proto15TerminalLocked): campaign.guard_fresh_entrypoint("RK4-4097")
            with self.assertRaises(runtime.Proto15TerminalLocked): campaign.guard_retry_entrypoint("RK4-2049")
            self.assertTrue(runtime.recover_campaign(root).terminal_lock)

    def test_recovery_rejects_hash_valid_pending_plan_substitution(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root)
            typed, previous = self.reject(campaign)
            pending = campaign.close_tdg6_outcome(
                "RK4-2049", typed,
                retry_successor={"predecessor_plan": previous, "successor_plan": plan(typed.retry_step_size)},
            )
            successor = pending.payload["retry_successor_payload_or_none"]["successor_plan"]
            endpoint = float.fromhex(successor["boundaries_hex"][-1])
            accepted_ledger = replace(
                campaign.ledgers["RK4-2049"], last_accepted_time=endpoint,
                accepted_macro_step_count=1, current_macro_step_temporal_retry_count=0,
                last_accepted_macro_step_temporal_retry_count=1,
            )
            campaign.accepted_fine_commit(
                "RK4-2049", state={"member": "RK4-2049", "v": 1},
                ledger=accepted_ledger, executed_plan=successor,
                accepted_boundary_time=endpoint, step_index=3, transaction_serial=4,
            )
            accepted = record_at(root, 3)
            substitute = plan(np.nextafter(typed.retry_step_size, float("inf")))
            self.assertEqual(substitute.endpoint.hex(), endpoint.hex())
            self.assertNotEqual(runtime._plan_mapping(substitute), successor)
            accepted["payload"]["executed_plan"] = runtime._plan_mapping(substitute)
            accepted = rewrite_record(root, accepted)
            checkpoint = checkpoint_at(root, 2)
            checkpoint["journal_tip_sha256"] = accepted["record_sha256"]
            rewrite_checkpoint(root, checkpoint)
            with self.assertRaisesRegex(
                runtime.Proto15CheckpointIntegrityError,
                "persisted accepted TDG7 plan differs",
            ):
                runtime.recover_campaign(root)

    def test_recovery_rejects_hash_valid_accepted_serial_and_ledger_forgery(self):
        for kind in ("serial", "ledger"):
            with self.subTest(kind=kind), TemporaryDirectory() as td:
                root = Path(td) / "campaign"; campaign = self.campaign(root)
                ledger = replace(campaign.ledgers["RK4-2049"], last_accepted_time=3 / 2, accepted_macro_step_count=1)
                campaign.accepted_fine_commit(
                    "RK4-2049", state={"member": "RK4-2049", "v": 1},
                    ledger=ledger, executed_plan=plan(), accepted_boundary_time=3 / 2,
                    step_index=3, transaction_serial=4,
                )
                checkpoint = checkpoint_at(root, 1)
                accepted = record_at(root, 1)
                cursor = dict(checkpoint["cursors"]["RK4-2049"])
                if kind == "serial":
                    cursor["previous_step_index"] = 1
                else:
                    checkpoint["ledgers"]["RK4-2049"]["accepted_macro_step_count"] = 2
                    cursor["TDG6_ledger_sha256"] = runtime._hash(checkpoint["ledgers"]["RK4-2049"])
                    accepted["payload"]["TDG6_ledger_sha256"] = cursor["TDG6_ledger_sha256"]
                cursor = dict(runtime._cursor(cursor).payload)
                checkpoint["cursors"]["RK4-2049"] = cursor
                accepted["payload"]["cursor"] = cursor
                accepted["payload"]["successor_cursor_chain_sha256"] = cursor["cursor_chain_sha256"]
                accepted = rewrite_record(root, accepted)
                checkpoint["journal_tip_sha256"] = accepted["record_sha256"]
                rewrite_checkpoint(root, checkpoint)
                message = "persisted accepted serials are not monotone" if kind == "serial" else "persisted accepted TDG6 ledger differs"
                with self.assertRaisesRegex(runtime.Proto15CheckpointIntegrityError, message):
                    runtime.recover_campaign(root)

    def test_semantic_rehash_mutation_corpus_rejects_time_terminal_and_rollback_edges(self):
        cases = {
            item["case"]: item
            for item in runtime._run_mutation_cases()
        }
        required = {
            "mutation:accepted_no_advance_forgery",
            "mutation:accepted_target_overrun_forgery",
            "mutation:terminal_older_rejection_pointer_forgery",
            "mutation:rollback_cursor_content_forgery",
        }
        self.assertEqual(len(cases), 14)
        self.assertTrue(required <= set(cases))
        for name in required:
            self.assertTrue(cases[name]["passed"], name)
            self.assertTrue(cases[name]["observed"].startswith("explicit_"), name)

    def test_rollback_restores_ancestor_and_issues_new_identities(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root); snapshot = campaign.common_event_snapshot(); prior_states = dict(campaign.states); prior_ids = {k: v.payload["attempt_id"] for k, v in campaign.cursors.items()}
            ledger = replace(campaign.ledgers["RK4-2049"], last_accepted_time=3 / 2, accepted_macro_step_count=1)
            campaign.accepted_fine_commit("RK4-2049", state={"member": "RK4-2049", "v": 1}, ledger=ledger, executed_plan=plan(), accepted_boundary_time=3 / 2, step_index=3, transaction_serial=4)
            aborted = campaign.cursors["RK4-2049"].payload["attempt_id"]
            campaign.common_event_rollback(snapshot, aborted_attempt_id=aborted)
            self.assertEqual(campaign._checkpoint["disposition"], "rolled_back"); self.assertEqual(campaign.states, prior_states)
            self.assertTrue(all(campaign.cursors[k].payload["attempt_id"] != prior_ids[k] for k in runtime.MEMBER_KEYS))
            rollback = record_at(root, 2)
            self.assertEqual(
                runtime._resolved_rollback_journal_fields(rollback),
                {
                    "new_journal_parent_sha256": rollback["journal_parent_sha256"],
                    "new_journal_record_sha256": rollback["record_sha256"],
                },
            )
            self.assertEqual(
                runtime.ROLLBACK_PROTOCOL_JOURNAL_FIELD_LOCATIONS,
                {
                    "new_journal_parent_sha256": "journal_envelope.journal_parent_sha256",
                    "new_journal_record_sha256": "journal_envelope.record_sha256",
                },
            )
            self.assertEqual(runtime.recover_campaign(root).checkpoint_sha256, campaign.checkpoint_sha256)
            self.assertNotIn(aborted, {c.payload["attempt_id"] for c in campaign.cursors.values()})

    def test_invalid_suffix_and_same_generation_fork_fail_closed(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root); self.reject(campaign)
            rejection = runtime.Proto15SyntheticCampaign._read_record(sorted((root / "journal").glob("*.journal"))[-1])
            bare = {"record_kind": "ACCEPTED_FINE_COMMIT", "journal_parent_sha256": rejection["record_sha256"], "payload": {"sequence": 2}}; record = {**bare, "record_sha256": runtime._hash(bare)}; encoded = runtime._canonical(record)
            (root / "journal" / f"{2:020d}-{record['record_sha256']}.journal").write_bytes(f"{len(encoded):08x} ".encode() + encoded + b"\n")
            with self.assertRaises(runtime.Proto15RecoveryInvalidStop): runtime.recover_campaign(root)
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; self.campaign(root); original = next((root / "checkpoints").glob("*.json")); value = json.loads(original.read_text()); value["checkpoint_sha256"] = "f" * 64
            (root / "checkpoints" / f"{0:020d}-{'f' * 64}.json").write_bytes(runtime._canonical(value))
            with self.assertRaises(runtime.Proto15CheckpointIntegrityError): runtime.recover_campaign(root)

    def test_recovery_rejects_hash_valid_retry_evidence_fork(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root)
            alternative_plan = plan(0.03125)
            alternative = outcome(campaign, "RK4-2049", alternative_plan)
            typed, previous = self.reject(campaign)
            campaign.close_tdg6_outcome(
                "RK4-2049", typed,
                retry_successor={"predecessor_plan": previous, "successor_plan": plan(typed.retry_step_size)},
            )
            rejection = record_at(root, 1)
            rejection["payload"]["TDG6_rejection_evidence"] = alternative.evidence.as_mapping()
            rejection["payload"]["predecessor_plan"] = runtime._plan_mapping(alternative_plan)
            rejection = rewrite_record(root, rejection)
            transition = record_at(root, 2)
            transition["journal_parent_sha256"] = rejection["record_sha256"]
            transition["payload"]["durable_rejection_record_sha256"] = rejection["record_sha256"]
            transition = rewrite_record(root, transition)
            checkpoint = checkpoint_at(root, 1)
            checkpoint["journal_tip_sha256"] = transition["record_sha256"]
            rewrite_checkpoint(root, checkpoint)
            with self.assertRaisesRegex(
                runtime.Proto15CheckpointIntegrityError,
                "retry evidence differs across durable records",
            ):
                runtime.recover_campaign(root)

    def test_recovery_derives_terminal_reason_from_rejection(self):
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root)
            tiny = plan(float(TDG6_MINIMUM_MACRO_STEP))
            typed, _ = self.reject(campaign, prior_plan=tiny, exhausted=True)
            campaign.close_tdg6_outcome("RK4-2049", typed)
            terminal = record_at(root, 2)
            terminal["payload"]["exhaustion_reason"] = "maximum_temporal_retries"
            terminal = rewrite_record(root, terminal)
            checkpoint = checkpoint_at(root, 1)
            checkpoint["terminal_closure"]["exhaustion_reason"] = "maximum_temporal_retries"
            checkpoint["journal_tip_sha256"] = terminal["record_sha256"]
            checkpoint["terminal_record_sha256"] = terminal["record_sha256"]
            rewrite_checkpoint(root, checkpoint)
            with self.assertRaisesRegex(
                runtime.Proto15CheckpointIntegrityError,
                "terminal record differs from checkpoint closure",
            ):
                runtime.recover_campaign(root)

    def test_recovery_rejects_hash_valid_rollback_ancestor_and_set_forgery(self):
        for kind in ("ancestor", "set"):
            with self.subTest(kind=kind), TemporaryDirectory() as td:
                root = Path(td) / "campaign"; campaign = self.campaign(root)
                snapshot = campaign.common_event_snapshot()
                ledger = replace(campaign.ledgers["RK4-2049"], last_accepted_time=3 / 2, accepted_macro_step_count=1)
                campaign.accepted_fine_commit(
                    "RK4-2049", state={"member": "RK4-2049", "v": 1},
                    ledger=ledger, executed_plan=plan(), accepted_boundary_time=3 / 2,
                    step_index=3, transaction_serial=4,
                )
                previous_checkpoint = checkpoint_at(root, 1)
                aborted = campaign.cursors["RK4-2049"].payload["attempt_id"]
                campaign.common_event_rollback(snapshot, aborted_attempt_id=aborted)
                rollback = record_at(root, 2)
                if kind == "ancestor":
                    rollback["payload"]["restored_snapshot_checkpoint_sha256"] = previous_checkpoint["checkpoint_sha256"]
                    rollback["payload"]["restored_snapshot_campaign_generation"] = 1
                    rollback["payload"]["restored_snapshot_cursor_set_sha256"] = previous_checkpoint["cursor_set_sha256"]
                    rollback["payload"]["restored_snapshot_ledger_set_sha256"] = previous_checkpoint["ledger_set_sha256"]
                    rollback["payload"]["restored_snapshot_accepted_state_set_sha256"] = previous_checkpoint["accepted_state_set_sha256"]
                else:
                    rollback["payload"]["new_cursor_set_sha256"] = "f" * 64
                rollback = rewrite_record(root, rollback)
                checkpoint = checkpoint_at(root, 2)
                checkpoint["journal_tip_sha256"] = rollback["record_sha256"]
                rewrite_checkpoint(root, checkpoint)
                message = "strict ancestor" if kind == "ancestor" else "set binding differs"
                with self.assertRaisesRegex(runtime.Proto15CheckpointIntegrityError, message):
                    runtime.recover_campaign(root)

    def test_missing_or_mutated_state_cursor_ledger_and_high_water_fail_closed(self):
        for kind in ("missing_state", "mutated_state", "cursor", "ledger", "high_water"):
            with self.subTest(kind=kind), TemporaryDirectory() as td:
                root = Path(td) / "campaign"; self.campaign(root); checkpoint_path = next((root / "checkpoints").glob("*.json")); value = json.loads(checkpoint_path.read_text())
                if kind == "missing_state": (root / "states" / f"{value['states']['RK4-2049']}.json").unlink()
                elif kind == "mutated_state": (root / "states" / f"{value['states']['RK4-2049']}.json").write_text('{"changed":true}')
                else:
                    if kind == "cursor": value["cursors"]["RK4-2049"]["attempt_id"] = "forged"
                    elif kind == "ledger": value["ledgers"]["RK4-2049"]["accepted_macro_step_count"] = 99
                    else: value["attempt_high_water"]["RK4-2049"] = -1
                    bare = dict(value); bare.pop("checkpoint_sha256"); value["checkpoint_sha256"] = runtime._hash(bare); checkpoint_path.unlink(); checkpoint_path = checkpoint_path.with_name(f"{0:020d}-{value['checkpoint_sha256']}.json"); checkpoint_path.write_bytes(runtime._canonical(value))
                with self.assertRaises((runtime.Proto15CheckpointIntegrityError, runtime.Proto15CursorContractError)): runtime.recover_campaign(root)

    def test_named_fault_hooks_cover_state_journal_checkpoint_and_abort_durably(self):
        seen = []
        with TemporaryDirectory() as td:
            root = Path(td) / "campaign"; campaign = self.campaign(root, fault_hook=lambda phase, target: seen.append((phase, target)))
            ledger = replace(campaign.ledgers["RK4-2049"], last_accepted_time=3 / 2, accepted_macro_step_count=1)
            campaign.accepted_fine_commit("RK4-2049", state={"member": "RK4-2049", "v": 1}, ledger=ledger, executed_plan=plan(), accepted_boundary_time=3 / 2, step_index=3, transaction_serial=4)
        self.assertEqual({phase for phase, _ in seen}, set(runtime.FAULT_PHASES)); self.assertTrue(all(any(target == expected for _, target in seen) for expected in ("state_object", "journal_record", "checkpoint")))
        with TemporaryDirectory() as td:
            with self.assertRaises(runtime.Proto15DurabilityError): self.campaign(Path(td) / "campaign", fault_hook=lambda phase, target: (_ for _ in ()).throw(OSError("fault")) if (phase, target) == ("before_replace", "checkpoint") else None)


if __name__ == "__main__": unittest.main()
