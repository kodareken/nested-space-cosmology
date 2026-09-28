"""Bounded coordinator tests for the first PRO19 event; no PDE is run."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as campaign_runtime
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore
from recursive_horizons.fgc.evolution.proto19_progression_contract import construct_first_event
from recursive_horizons.fgc.evolution.proto19_progression_inputs import reconstruct_gr0_members
from tests.fgc_gen0_fixture import (
    copy_sealed_gen0_store,
    reconstruct_sealed_gen0_members,
)


def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "fgc_pro19_event1_runner_test_module",
        ROOT / "scripts" / "run_fgc_pro19_event1.py",
    )
    if spec is None or spec.loader is None:
        raise AssertionError("event-one runner cannot be imported")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


runner = _load_runner()


class _State:
    def __init__(self, key: str, time: float, *, mode: str = "FRESH_READY", owner=None):
        self.cursor = {
            "accepted_boundary_time": {
                "rational": "23/16" if time == 23 / 16 else "3/2",
                "binary64_hex": time.hex(),
            },
            "mode": mode,
        }
        self.pending_owner = owner
        self.pending_cap_hex = None
        self.descriptor_sha256 = "6" * 64
        self.cfl_current = 0
        self.cfl_total = 0
        self.ledger = {"cumulative_temporal_retry_count": 0}


class _Checkpoint:
    def __init__(self, times, *, event: int = 23, disposition: str = "nonterminal"):
        self.authorization_commit = "2" * 40
        self.plan_sha256 = "3" * 64
        self.campaign_id = "campaign"
        self.protocol = "FGC-2-SF1-PROTO18"
        self.generation = 1
        self.event = event
        self.target = {"rational": "3/2", "binary64_hex": (3 / 2).hex()}
        self.members = {
            key: _State(key, times[key]) for key in runner.MEMBER_KEYS
        }
        self.journal_sequence = 0
        self.journal_tip_sha256 = "4" * 64
        self.disposition = disposition
        self.terminal = None
        self.parent_sha256 = "0" * 64
        self.sha256 = "5" * 64


def _sid1_receipt(checkpoint):
    return SimpleNamespace(
        continuation_commit="9" * 40,
        manifest_sha256="a" * 64,
        campaign_id=runner.sid1.CAMPAIGN_ID,
        recovery_checkpoint_generation=7,
        recovery_checkpoint_sha256=runner.sid1.CHECKPOINT_SHA256,
        suffix_sequence=7,
        suffix_sha256=runner.sid1.TDG6_SUFFIX_SHA256,
        descriptor_sha256=runner.sid1.DESCRIPTOR_SHA256,
        evolution_state_sha256=runner.sid1.EVOLUTION_STATE_SHA256,
        expected_transition_sha256=runner.sid1.RECOVERY_TRANSITION_SHA256,
        expected_checkpoint_sha256=checkpoint.sha256,
        expected_cursor_sha256=runner.sid1.RECOVERY_CURSOR_SHA256,
        expected_pending_cap_hex=runner.sid1.RECOVERY_PENDING_CAP_HEX,
        original_plan_sha256=runner.sid1.PLAN_SHA256,
        progression_plan=SimpleNamespace(branch="GR-0", amplitude="3"),
    )


def _sid1_successor():
    checkpoint = _Checkpoint({key: 23 / 16 for key in runner.MEMBER_KEYS})
    checkpoint.authorization_commit = runner.sid1.ORIGINAL_AUTHORIZATION_COMMIT
    checkpoint.plan_sha256 = runner.sid1.PLAN_SHA256
    checkpoint.campaign_id = runner.sid1.CAMPAIGN_ID
    checkpoint.generation = 8
    checkpoint.parent_sha256 = runner.sid1.CHECKPOINT_SHA256
    checkpoint.journal_sequence = 8
    checkpoint.journal_tip_sha256 = runner.sid1.RECOVERY_TRANSITION_SHA256
    checkpoint.sha256 = runner.sid1.RECOVERY_CHECKPOINT_SHA256
    member = checkpoint.members[runner.sid1.MEMBER_KEY]
    member.descriptor_sha256 = runner.sid1.DESCRIPTOR_SHA256
    member.cursor.update({
        "accepted_state_sha256": runner.sid1.DESCRIPTOR_SHA256,
        "cursor_chain_sha256": runner.sid1.RECOVERY_CURSOR_SHA256,
        "mode": "RETRY_PENDING",
    })
    member.cursor["accepted_boundary_time"] = {
        "rational": "206891375579527/140737488355328",
        "binary64_hex": runner.sid1.MEMBER_TIME_HEX,
    }
    member.cfl_current = 1
    member.cfl_total = 1
    member.ledger["cumulative_temporal_retry_count"] = 1
    member.pending_owner = "temporal"
    member.pending_cap_hex = runner.sid1.RECOVERY_PENDING_CAP_HEX
    return checkpoint


def _sid1_preview(checkpoint):
    return SimpleNamespace(
        predecessor_checkpoint_sha256=runner.sid1.CHECKPOINT_SHA256,
        checkpoint=checkpoint,
        records=(
            {
                "kind": "tdg6_rejection",
                "record_sha256": runner.sid1.TDG6_SUFFIX_SHA256,
            },
            {
                "kind": "cursor_transition",
                "record_sha256": runner.sid1.RECOVERY_TRANSITION_SHA256,
                "payload": {
                    "successor_cursor_sha256": runner.sid1.RECOVERY_CURSOR_SHA256,
                },
            },
        ),
        evolution_state_sha256=runner.sid1.EVOLUTION_STATE_SHA256,
        physical_state_advanced=False,
    )


class Proto19Event1RunnerTests(unittest.TestCase):
    def test_member_schedule_is_canonical_and_pending_owner_has_priority(self) -> None:
        times = {key: 23 / 16 for key in runner.MEMBER_KEYS}
        checkpoint = _Checkpoint(times)
        self.assertEqual(runner._member_to_advance(checkpoint), runner.MEMBER_KEYS[0])
        checkpoint.members[runner.MEMBER_KEYS[-1]].pending_owner = "source"
        self.assertEqual(runner._member_to_advance(checkpoint), runner.MEMBER_KEYS[-1])
        checkpoint.members[runner.MEMBER_KEYS[-2]].pending_owner = "cfl"
        with self.assertRaisesRegex(runner.Proto19Event1RunnerError, "more than one"):
            runner._member_to_advance(checkpoint)

    def test_all_target_members_commit_only_the_common_event(self) -> None:
        before = _Checkpoint({key: 3 / 2 for key in runner.MEMBER_KEYS})
        after = _Checkpoint(
            {key: 3 / 2 for key in runner.MEMBER_KEYS},
            event=24,
            disposition="event_complete",
        )
        after.target = {"rational": "25/16", "binary64_hex": (25 / 16).hex()}
        receipt = SimpleNamespace(progression_plan=SimpleNamespace())
        reconstructed = SimpleNamespace(members={key: object() for key in runner.MEMBER_KEYS})
        store = SimpleNamespace(
            runtime_checkpoint_present=lambda: False,
            active_writer_capability=lambda _checkpoint: object(),
        )
        event = SimpleNamespace(
            checkpoint=after,
            record_sha256s=("6" * 64,),
        )
        with patch.object(runner, "authorize_and_reconstruct", return_value=(receipt, reconstructed)), \
             patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner, "_recoverable_bootstrap_prefix", return_value=True), \
             patch.object(runner.runtime, "initialize_or_recover_generation_one", return_value=before), \
             patch.object(runner, "_acquire_writer", return_value="7" * 64), \
             patch.object(runner, "_reconcile_before_attempt", return_value=(before, None)), \
             patch.object(runner.runtime, "commit_common_event", return_value=event) as commit, \
             patch.object(runner.runtime, "execute_one_attempt") as execute, \
             patch.object(runner, "_release_writer_if_closed", return_value=True):
            result = runner.run_first_event(
                ROOT,
                authorization_commit="2" * 40,
                manifest_path=ROOT / "manifest.toml",
                store_root=ROOT / "unused-store",
            )
        commit.assert_called_once_with(before, store)
        execute.assert_not_called()
        self.assertEqual(result["state"], "event_one_complete")
        self.assertFalse(result["candidate_branch_opened"])

    def test_payload_orphan_identity_is_forwarded_to_the_exact_member_replay(self) -> None:
        before = _Checkpoint({key: 23 / 16 for key in runner.MEMBER_KEYS})
        terminal = _Checkpoint(
            {key: 23 / 16 for key in runner.MEMBER_KEYS},
            disposition="invalid_terminal",
        )
        receipt = SimpleNamespace(progression_plan=SimpleNamespace())
        reconstructed = SimpleNamespace(members={key: object() for key in runner.MEMBER_KEYS})
        store = SimpleNamespace(
            runtime_checkpoint_present=lambda: False,
            authenticated_snapshot=lambda: SimpleNamespace(terminal_lock_present=True),
            active_writer_capability=lambda _checkpoint: object(),
        )
        attempt = SimpleNamespace(
            checkpoint=terminal,
            record_sha256s=("8" * 64,),
        )
        expected = "9" * 64
        with patch.object(runner, "authorize_and_reconstruct", return_value=(receipt, reconstructed)), \
             patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner, "_recoverable_bootstrap_prefix", return_value=True), \
             patch.object(runner.runtime, "initialize_or_recover_generation_one", return_value=before), \
             patch.object(runner, "_acquire_writer", return_value="7" * 64), \
             patch.object(runner, "_reconcile_before_attempt", return_value=(before, expected)), \
             patch.object(runner.copy, "deepcopy", return_value=object()), \
             patch.object(runner.runtime, "execute_one_attempt", return_value=attempt) as execute, \
             patch.object(runner, "_release_writer_if_closed", return_value=True):
            result = runner.run_first_event(
                ROOT,
                authorization_commit="2" * 40,
                manifest_path=ROOT / "manifest.toml",
                store_root=ROOT / "unused-store",
            )
        self.assertEqual(
            execute.call_args.kwargs["expected_orphan_semantic_sha256"], expected
        )
        self.assertEqual(execute.call_args.kwargs["key"], runner.MEMBER_KEYS[0])
        self.assertEqual(result["state"], "invalid_terminal")

    def test_ambiguous_exception_deliberately_leaves_recovery_lease(self) -> None:
        before = _Checkpoint({key: 23 / 16 for key in runner.MEMBER_KEYS})
        receipt = SimpleNamespace(progression_plan=SimpleNamespace())
        reconstructed = SimpleNamespace(members={key: object() for key in runner.MEMBER_KEYS})
        with patch.object(runner, "authorize_and_reconstruct", return_value=(receipt, reconstructed)), \
             patch.object(
                 runner,
                 "HLT16CampaignStore",
                 return_value=SimpleNamespace(
                     runtime_checkpoint_present=lambda: False,
                     active_writer_capability=lambda _checkpoint: object(),
                 ),
             ), \
             patch.object(runner, "_recoverable_bootstrap_prefix", return_value=True), \
             patch.object(runner.runtime, "initialize_or_recover_generation_one", return_value=before), \
             patch.object(runner, "_acquire_writer", return_value="7" * 64), \
             patch.object(runner, "_reconcile_before_attempt", side_effect=runner.Proto19Event1RunnerError("ambiguous")), \
             patch.object(runner, "_release_writer_if_closed") as release:
            with self.assertRaisesRegex(runner.Proto19Event1RunnerError, "ambiguous"):
                runner.run_first_event(
                    ROOT,
                    authorization_commit="2" * 40,
                    manifest_path=ROOT / "manifest.toml",
                    store_root=ROOT / "unused-store",
                )
        release.assert_not_called()

    def test_persisted_generation_uses_static_shells_without_reopening_gen0(self) -> None:
        terminal = _Checkpoint(
            {key: 23 / 16 for key in runner.MEMBER_KEYS},
            disposition="invalid_terminal",
        )
        receipt = SimpleNamespace(progression_plan=SimpleNamespace())
        store = SimpleNamespace(
            runtime_checkpoint_present=lambda: True,
            authenticated_snapshot=lambda: SimpleNamespace(terminal_lock_present=True),
            active_writer_capability=lambda _checkpoint: object(),
        )
        templates = {key: object() for key in runner.MEMBER_KEYS}
        with patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner, "_recoverable_bootstrap_prefix", return_value=False), \
             patch.object(runner, "authorize_event_authority", return_value=receipt), \
             patch.object(runner, "reconstruct_gr0_members", side_effect=AssertionError("historical raw reopen")) as raw, \
             patch.object(runner.runtime, "recover_persisted_generation_one", return_value=terminal), \
             patch.object(runner, "_acquire_writer", return_value="7" * 64), \
             patch.object(runner, "_reconcile_before_attempt", return_value=(terminal, None)), \
             patch.object(runner, "_release_writer_if_closed", return_value=True):
            result = runner.run_first_event(
                ROOT,
                authorization_commit="2" * 40,
                manifest_path=ROOT / "manifest.toml",
                store_root=ROOT / "unused-store",
                reconstruct=runner.reconstruct_gr0_members,
                static_shells=lambda _root: templates,
            )
        raw.assert_not_called()
        self.assertEqual(result["state"], "invalid_terminal")

    def test_continuation_preflight_authenticates_exact_checkpoint_without_writer(self) -> None:
        checkpoint = _Checkpoint({key: 23 / 16 for key in runner.MEMBER_KEYS})
        checkpoint.generation = 6
        receipt = SimpleNamespace(
            continuation_commit="9" * 40,
            original_authorization_commit="2" * 40,
            original_plan_sha256="3" * 64,
            campaign_id="campaign",
            recovery_checkpoint_generation=6,
            recovery_checkpoint_sha256=checkpoint.sha256,
            recovery_journal_sequence=checkpoint.journal_sequence,
            recovery_journal_tip_sha256=checkpoint.journal_tip_sha256,
            recovery_member_key="RK4-2049",
            recovery_member_descriptor_sha256=(
                checkpoint.members["RK4-2049"].descriptor_sha256
            ),
            progression_plan=SimpleNamespace(branch="GR-0"),
        )
        store = object()
        with patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner, "_recoverable_bootstrap_prefix", return_value=False), \
             patch.object(runner, "authorize_continuation_event_authority", return_value=receipt), \
             patch.object(runner.runtime, "recover_persisted_generation_one", return_value=checkpoint), \
             patch.object(runner, "_acquire_writer") as acquire:
            result = runner.continuation_preflight(
                ROOT,
                continuation_commit="9" * 40,
                manifest_path=ROOT / "manifest.toml",
                store_root=ROOT / "unused-store",
            )
        acquire.assert_not_called()
        self.assertTrue(result["safe_to_resume"])
        self.assertFalse(result["state_advanced"])
        self.assertFalse(result["candidate_branch_opened"])

    def test_corrected_resume_uses_persisted_shells_and_original_plan_only(self) -> None:
        checkpoint = _Checkpoint({key: 23 / 16 for key in runner.MEMBER_KEYS})
        checkpoint.generation = 6
        original_plan = SimpleNamespace(branch="GR-0")
        receipt = SimpleNamespace(
            continuation_commit="9" * 40,
            continuation_config_sha256="a" * 64,
            continuation_result_sha256="b" * 64,
            recovery_checkpoint_generation=6,
            recovery_checkpoint_sha256=checkpoint.sha256,
            recovery_journal_sequence=checkpoint.journal_sequence,
            recovery_journal_tip_sha256=checkpoint.journal_tip_sha256,
            recovery_member_key="RK4-2049",
            recovery_member_descriptor_sha256=(
                checkpoint.members["RK4-2049"].descriptor_sha256
            ),
            progression_plan=original_plan,
        )
        templates = {key: object() for key in runner.MEMBER_KEYS}
        store = object()
        terminal = {
            "state": "event_one_complete",
            "candidate_branch_opened": False,
            "physical_result_earned": False,
        }
        with patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner, "_recoverable_bootstrap_prefix", return_value=False), \
             patch.object(runner, "authorize_continuation_event_authority", return_value=receipt), \
             patch.object(runner.runtime, "recover_persisted_generation_one", return_value=checkpoint) as recover, \
             patch.object(runner, "reconstruct_gr0_members", side_effect=AssertionError("GEN0 reopened")) as raw, \
             patch.object(runner, "_advance_authenticated_event", return_value=terminal) as advance:
            result = runner.resume_corrected_first_event(
                ROOT,
                continuation_commit="9" * 40,
                manifest_path=ROOT / "manifest.toml",
                store_root=ROOT / "unused-store",
                static_shells=lambda _root: templates,
            )
        raw.assert_not_called()
        recover.assert_called_once_with(original_plan, store)
        advance.assert_called_once_with(
            original_plan,
            store,
            checkpoint,
            templates,
            attempt_runner=None,
        )
        self.assertEqual(result["continuation_commit"], "9" * 40)
        self.assertFalse(result["earlier_stop_physical_inference"])

    def test_continuation_refuses_store_movement_after_anchor_capture(self) -> None:
        checkpoint = _Checkpoint({key: 23 / 16 for key in runner.MEMBER_KEYS})
        checkpoint.generation = 7
        receipt = SimpleNamespace(
            continuation_commit="9" * 40,
            original_authorization_commit="2" * 40,
            original_plan_sha256="3" * 64,
            campaign_id="campaign",
            recovery_checkpoint_generation=6,
            recovery_checkpoint_sha256="c" * 64,
            recovery_journal_sequence=5,
            recovery_journal_tip_sha256="d" * 64,
            recovery_member_key="RK4-2049",
            recovery_member_descriptor_sha256="e" * 64,
            progression_plan=SimpleNamespace(branch="GR-0"),
        )
        store = object()
        with patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner, "_recoverable_bootstrap_prefix", return_value=False), \
             patch.object(runner, "authorize_continuation_event_authority", return_value=receipt), \
             patch.object(runner.runtime, "recover_persisted_generation_one", return_value=checkpoint), \
             patch.object(runner, "_acquire_writer") as acquire:
            with self.assertRaisesRegex(
                runner.Proto19Event1RunnerError,
                "checkpoint moved after authority capture",
            ):
                runner.continuation_preflight(
                    ROOT,
                    continuation_commit="9" * 40,
                    manifest_path=ROOT / "manifest.toml",
                    store_root=ROOT / "unused-store",
                )
        acquire.assert_not_called()

    def test_corrected_resume_never_reopens_generation_zero(self) -> None:
        store = object()
        with patch.object(runner, "HLT16CampaignStore", return_value=store), \
             patch.object(runner, "_recoverable_bootstrap_prefix", return_value=True), \
             patch.object(runner, "authorize_continuation_event_authority") as authority:
            with self.assertRaisesRegex(
                runner.Proto19Event1RunnerError,
                "cannot reopen generation-zero",
            ):
                runner.resume_corrected_first_event(
                    ROOT,
                    continuation_commit="9" * 40,
                    manifest_path=ROOT / "manifest.toml",
                    store_root=ROOT / "unused-store",
                )
        authority.assert_not_called()

    def test_partial_runtime_suffix_refuses_preflight_before_raw_reconstruction(self) -> None:
        with TemporaryDirectory() as temporary:
            copied = Path(temporary) / "calibration"
            copy_sealed_gen0_store(ROOT, copied)
            (copied / "payloads").mkdir()
            (copied / "payloads" / ("a" * 64 + ".npz")).write_bytes(b"partial")
            with patch.object(runner, "authorize_and_reconstruct") as raw:
                with self.assertRaisesRegex(
                    runner.Proto19Event1RunnerError,
                    "post-GEN0 runtime, staging, lock, or foreign suffix",
                ):
                    runner.preflight(
                        ROOT,
                        authorization_commit="2" * 40,
                        manifest_path=ROOT / "manifest.toml",
                        store_root=copied,
                    )
        raw.assert_not_called()

    def test_foreign_pre_generation_prefix_refuses_run_before_raw_reconstruction(self) -> None:
        with TemporaryDirectory() as temporary:
            copied = Path(temporary) / "calibration"
            copy_sealed_gen0_store(ROOT, copied)
            (copied / "journal").mkdir()
            (copied / "journal" / "foreign.journal").write_bytes(b"partial")
            with patch.object(runner, "authorize_and_reconstruct") as raw:
                with self.assertRaisesRegex(
                    runner.Proto19Event1RunnerError,
                    "bootstrap prefix has a foreign runtime leaf",
                ):
                    runner.run_first_event(
                        ROOT,
                        authorization_commit="2" * 40,
                        manifest_path=ROOT / "manifest.toml",
                        store_root=copied,
                    )
        raw.assert_not_called()

    def test_unlinked_bootstrap_stage_refuses_run_before_raw_reconstruction(self) -> None:
        reconstructed = reconstruct_sealed_gen0_members(ROOT)
        plan = construct_first_event(ROOT, authorization_commit="2" * 40)
        with TemporaryDirectory() as temporary:
            copied = Path(temporary) / "calibration"
            copy_sealed_gen0_store(ROOT, copied)
            store = HLT16CampaignStore(
                copied,
                fault_hook=lambda observed: (
                    (_ for _ in ()).throw(RuntimeError(observed))
                    if observed == "state:after_write" else None
                ),
            )
            with self.assertRaisesRegex(RuntimeError, "state:after_write"):
                campaign_runtime.initialize_or_recover_generation_one(
                    plan, store, reconstructed,
                )
            before = runner._walk_tree(copied)[0]
            with patch.object(runner, "authorize_and_reconstruct") as raw:
                with self.assertRaisesRegex(
                    runner.Proto19Event1RunnerError,
                    "unlinked bootstrap stage ownership is unprovable",
                ):
                    runner.run_first_event(
                        ROOT,
                        authorization_commit="2" * 40,
                        manifest_path=ROOT / "manifest.toml",
                        store_root=copied,
                    )
            raw.assert_not_called()
            after = runner._walk_tree(copied)[0]
            self.assertEqual(after, before)

    def test_production_cli_rejects_alternate_paths_before_preflight(self) -> None:
        alternate = ROOT / "alternate-store-forbidden"
        with patch.object(sys, "argv", [
            str(ROOT / "scripts/run_fgc_pro19_event1.py"),
            "--preflight", "--store-root", str(alternate),
        ]), patch.object(runner, "preflight") as preflight:
            self.assertEqual(runner.main(), 2)
        preflight.assert_not_called()

    def test_bootstrap_publication_cuts_select_exact_replay_or_persisted_recovery(self) -> None:
        """A real copied GEN0 bridge never falls back after an interrupted cut."""
        reconstructed = reconstruct_sealed_gen0_members(ROOT)
        plan = construct_first_event(ROOT, authorization_commit="2" * 40)
        for phase, expected_bootstrap in (
            ("state:after_link", True),
            ("journal:after_link", True),
            ("checkpoint:after_link", False),
        ):
            with self.subTest(phase=phase), TemporaryDirectory() as temporary:
                copied = Path(temporary) / "calibration"
                copy_sealed_gen0_store(ROOT, copied)
                store = HLT16CampaignStore(
                    copied,
                    fault_hook=lambda observed, cut=phase: (
                        (_ for _ in ()).throw(RuntimeError(cut))
                        if observed == cut else None
                    ),
                )
                with self.assertRaisesRegex(RuntimeError, phase):
                    campaign_runtime.initialize_or_recover_generation_one(
                        plan, store, reconstructed,
                    )
                self.assertEqual(
                    runner._recoverable_bootstrap_prefix(store),
                    expected_bootstrap,
                )
                store.fault_hook = None
                checkpoint = campaign_runtime.initialize_or_recover_generation_one(
                    plan, store, reconstructed,
                )
                self.assertEqual(checkpoint.generation, 1)

    def test_sid1_preflight_reports_exact_no_advance_edge_without_writing(self) -> None:
        checkpoint = _sid1_successor()
        receipt = _sid1_receipt(checkpoint)
        preview = _sid1_preview(checkpoint)
        with patch.object(
            runner,
            "_sid1_preflight_context",
            return_value=(object(), receipt, preview),
        ), patch.object(runner, "_acquire_writer") as acquire, patch.object(
            runner.recovery, "reconcile"
        ) as reconcile, patch.object(
            runner, "_advance_authenticated_event"
        ) as advance:
            result = runner.sid1_preflight(
                ROOT,
                continuation_commit="9" * 40,
                manifest_path=ROOT / "manifest.toml",
                store_root=ROOT / "unused-store",
            )
        acquire.assert_not_called()
        reconcile.assert_not_called()
        advance.assert_not_called()
        self.assertTrue(result["safe_to_recover"])
        self.assertFalse(result["safe_to_resume_trajectory"])
        self.assertFalse(result["physical_state_advanced"])

    def test_sid1_successor_rejects_each_frozen_identity_or_retry_mutation(self) -> None:
        fields = (
            "checkpoint",
            "transition",
            "cursor",
            "cfl",
            "temporal",
            "cap",
        )
        for field in fields:
            checkpoint = _sid1_successor()
            receipt = _sid1_receipt(checkpoint)
            member = checkpoint.members[runner.sid1.MEMBER_KEY]
            if field == "checkpoint":
                checkpoint.sha256 = "0" * 64
            elif field == "transition":
                checkpoint.journal_tip_sha256 = "0" * 64
            elif field == "cursor":
                member.cursor["cursor_chain_sha256"] = "0" * 64
            elif field == "cfl":
                member.cfl_total = 2
            elif field == "temporal":
                member.ledger["cumulative_temporal_retry_count"] = 2
            else:
                member.pending_cap_hex = "0x1.0p-9"
            with self.subTest(field=field), self.assertRaises(
                runner.Proto19Event1RunnerError
            ):
                runner._require_exact_sid1_successor(checkpoint, receipt)

    def test_sid1_recovery_publishes_only_one_metadata_successor(self) -> None:
        predecessor = _Checkpoint({key: 23 / 16 for key in runner.MEMBER_KEYS})
        predecessor.sha256 = runner.sid1.CHECKPOINT_SHA256
        successor = _sid1_successor()
        receipt = _sid1_receipt(successor)
        preview = _sid1_preview(successor)
        store = SimpleNamespace(
            root=ROOT / "unused-store",
            authenticated_snapshot=lambda: SimpleNamespace(checkpoint=predecessor),
        )
        post_store = SimpleNamespace(root=store.root)
        final_store = object()
        decision = SimpleNamespace(
            disposition="reconciled_checkpoint",
            published=True,
            checkpoint_sha256=successor.sha256,
        )
        with patch.object(
            runner,
            "_sid1_preflight_context",
            side_effect=[(store, receipt, preview), (post_store, receipt, preview)],
        ), patch.object(runner, "_acquire_writer", return_value="7" * 64) as acquire, \
             patch.object(runner.recovery, "reconcile", return_value=decision) as reconcile, \
             patch.object(
                 runner,
                 "_require_exact_sid1_post_recovery",
                 side_effect=[successor, successor],
             ) as postcheck, patch.object(
                 runner, "_release_writer_if_closed", return_value=True
             ) as release, patch.object(
                 runner, "HLT16CampaignStore", return_value=final_store
             ), patch.object(
                 runner, "_advance_authenticated_event"
             ) as advance, patch.object(
                 runner.runtime, "execute_one_attempt"
             ) as attempt:
            result = runner.sid1_recover(
                ROOT,
                continuation_commit="9" * 40,
                manifest_path=ROOT / "manifest.toml",
                store_root=store.root,
            )
        acquire.assert_called_once_with(store, predecessor)
        reconcile.assert_called_once_with(store)
        release.assert_called_once_with(store, owner_token="7" * 64)
        postcheck.assert_any_call(store, receipt, expected_owner_token="7" * 64)
        postcheck.assert_any_call(final_store, receipt, expected_owner_token=None)
        advance.assert_not_called()
        attempt.assert_not_called()
        self.assertTrue(result["recovery_checkpoint_published"])
        self.assertFalse(result["physical_state_advanced"])
        self.assertFalse(result["rejected_proposal_replayed"])

    def test_sid1_recovery_ambiguity_leaves_writer_as_recovery_marker(self) -> None:
        predecessor = _Checkpoint({key: 23 / 16 for key in runner.MEMBER_KEYS})
        predecessor.sha256 = runner.sid1.CHECKPOINT_SHA256
        successor = _sid1_successor()
        receipt = _sid1_receipt(successor)
        preview = _sid1_preview(successor)
        store = SimpleNamespace(
            root=ROOT / "unused-store",
            authenticated_snapshot=lambda: SimpleNamespace(checkpoint=predecessor),
        )
        with patch.object(
            runner,
            "_sid1_preflight_context",
            side_effect=[(store, receipt, preview), (store, receipt, preview)],
        ), patch.object(runner, "_acquire_writer", return_value="7" * 64), \
             patch.object(
                 runner.recovery,
                 "reconcile",
                 side_effect=runner.recovery.HLT16RecoveryError("ambiguous"),
             ), patch.object(runner, "_release_writer_if_closed") as release:
            with self.assertRaisesRegex(
                runner.Proto19Event1RunnerError, "publication failed"
            ):
                runner.sid1_recover(
                    ROOT,
                    continuation_commit="9" * 40,
                    manifest_path=ROOT / "manifest.toml",
                    store_root=store.root,
                )
        release.assert_not_called()

    def test_sid2_binder_does_not_unlock_sid1_numerical_resume(self) -> None:
        with patch.object(runner, "HLT16CampaignStore") as store, patch.object(
            runner, "authorize_sid1_event_authority"
        ) as authority, patch.object(
            runner, "_advance_authenticated_event"
        ) as advance, patch.object(
            runner, "reconstruct_gr0_members"
        ) as gen0:
            with self.assertRaisesRegex(
                runner.Proto19Event1RunnerError, "resume authority and read-only preflight"
            ):
                runner.sid1_resume(
                    ROOT,
                    continuation_commit="9" * 40,
                    manifest_path=ROOT / "manifest.toml",
                    store_root=ROOT / "unused-store",
                )
        store.assert_not_called()
        authority.assert_not_called()
        advance.assert_not_called()
        gen0.assert_not_called()

    def test_real_copied_sid1_suffix_recovers_exactly_without_advancing_arrays(self) -> None:
        with TemporaryDirectory() as temporary:
            copied = Path(temporary) / "calibration"
            shutil.copytree(ROOT / runner.sid1.STORE_PATH, copied)
            # The live campaign has lawfully crossed SID1's one-time temporal
            # boundary. Reconstruct the exact pre-recovery prefix only inside
            # this disposable copy by removing SID1's two metadata-only
            # successor leaves; every descriptor and payload remains intact.
            (copied / "checkpoints" / (
                f"{8:020d}-{runner.sid1.RECOVERY_CHECKPOINT_SHA256}.json"
            )).unlink()
            (copied / "journal" / (
                f"{8:020d}-{runner.sid1.RECOVERY_TRANSITION_SHA256}.journal"
            )).unlink()
            # PREF28 later sealed a generic generation-nine terminal over the
            # recovered generation-eight boundary. Remove that exact later
            # metadata suffix from this disposable copy as well; the test owns
            # SID1 recovery, not the consumed SID3 continuation attempt.
            generation9 = list((copied / "checkpoints").glob(f"{9:020d}-*.json"))
            sequence9 = list((copied / "journal").glob(f"{9:020d}-*.journal"))
            self.assertEqual(len(generation9), 1)
            self.assertEqual(len(sequence9), 1)
            generation9[0].unlink()
            sequence9[0].unlink()
            (copied / "locks" / "terminal.lock").unlink()
            for path in (copied / "locks").glob(
                ".active-write.lock.hlt16-quarantine-*"
            ):
                if json.loads(path.read_text(encoding="ascii")).get(
                    "checkpoint_sha256"
                ) == runner.sid1.RECOVERY_CHECKPOINT_SHA256:
                    path.unlink()
            retired = [
                path
                for path in (copied / "locks").glob(
                    ".active-write.lock.hlt16-quarantine-*"
                )
                if json.loads(path.read_text(encoding="ascii")).get(
                    "checkpoint_sha256"
                )
                == runner.sid1.CHECKPOINT_SHA256
            ]
            self.assertEqual(len(retired), 1)
            retired[0].rename(copied / "locks" / "active-write.lock")
            store = runner.HLT16CampaignStore(copied)
            predecessor = store.authenticated_snapshot().checkpoint
            self.assertEqual(predecessor.sha256, runner.sid1.CHECKPOINT_SHA256)
            preview = runner.recovery.preview_tdg6_recovery(store)
            receipt = _sid1_receipt(preview.checkpoint)
            runner._require_exact_sid1_preview(preview, receipt)
            token = runner._acquire_writer(store, predecessor)
            decision = runner.recovery.reconcile(store)
            self.assertEqual(decision.checkpoint_sha256, receipt.expected_checkpoint_sha256)
            live = runner._require_exact_sid1_post_recovery(
                store, receipt, expected_owner_token=token
            )
            self.assertEqual(live.sha256, receipt.expected_checkpoint_sha256)
            with self.assertRaisesRegex(
                runner.Proto19Event1RunnerError, "writer ownership"
            ):
                runner._require_exact_sid1_post_recovery(
                    store, receipt, expected_owner_token="0" * 64
                )
            self.assertTrue(
                runner._release_writer_if_closed(store, owner_token=token)
            )
            restarted = runner._require_exact_sid1_post_recovery(
                runner.HLT16CampaignStore(copied),
                receipt,
                expected_owner_token=None,
            )
            self.assertEqual(restarted.sha256, receipt.expected_checkpoint_sha256)
            self.assertEqual(
                restarted.members[runner.sid1.MEMBER_KEY].descriptor_sha256,
                runner.sid1.DESCRIPTOR_SHA256,
            )


if __name__ == "__main__":
    unittest.main()
