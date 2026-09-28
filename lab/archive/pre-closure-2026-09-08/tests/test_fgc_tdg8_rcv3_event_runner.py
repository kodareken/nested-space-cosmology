"""Focused orchestration tests for the bounded RCV3 TDG8 event runner."""

from __future__ import annotations

import ast
from dataclasses import replace
import inspect
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from scripts import run_fgc_tdg8_rcv3_event as runner


ROOT = Path(__file__).resolve().parents[1]
AUTHORITY_COMMIT = "a" * 40
CONFIG_RAW = b"frozen RCV3 authority config\n"
RESULT_RAW = b"frozen RCV3 authority result\n"


def _receipt(**changes):
    values = {
        "authority_commit": AUTHORITY_COMMIT,
        "pref1_authority_commit": runner.authority.PREF1_AUTHORITY_COMMIT,
        "pref1_result_sha256": runner.authority.PREF1_RESULT_SHA256,
        "projection_id": runner.authority.PROJECTION_ID,
        "receipt_sha256": runner.authority.RECEIPT_SHA256,
        "original_execution_commit": runner.ORIGINAL_EXECUTION_COMMIT,
        "original_plan_sha256": runner.ORIGINAL_PLAN_SHA256,
        "config_sha256": "b" * 64,
        "result_sha256": "c" * 64,
        "campaign_id": runner.authority.CAMPAIGN_ID,
        "destination_path": runner.authority.DESTINATION_PATH,
        "entry_checkpoint_sha256": runner.authority.ENTRY_CHECKPOINT_SHA256,
        "entry_journal_sha256": runner.authority.ENTRY_JOURNAL_SHA256,
        "temporal_retry_cap": runner.authority.TEMPORAL_RETRY_CAP,
        "minimum_macro_step_hex": runner.authority.MINIMUM_MACRO_STEP_HEX,
        "implementation_inventory": (),
        "environment": {},
        "candidate_branches_forbidden": True,
    }
    values.update(changes)
    return runner.authority.TDG8RCV3ExecutionAuthority(**values)


def _plan(**changes):
    values = {
        "sha256": runner.ORIGINAL_PLAN_SHA256,
        "authorization_commit": runner.ORIGINAL_EXECUTION_COMMIT,
        "protocol_artifact_id": runner.EXPECTED_PROTOCOL,
        "campaign_id": runner.authority.CAMPAIGN_ID,
        "branch": "GR-0",
        "amplitude": "3",
        "common_event_index": runner.authority.EVENT,
        "start_time": SimpleNamespace(rational="23/16"),
        "target_time": SimpleNamespace(rational="3/2"),
        "members": tuple(
            SimpleNamespace(member_key=key) for key in runner.MEMBER_KEYS
        ),
    }
    values.update(changes)
    return SimpleNamespace(**values)


def _boundary(**changes):
    values = {
        "checkpoint_generation": 9,
        "checkpoint_sha256": runner.authority.ENTRY_CHECKPOINT_SHA256,
        "journal_sequence": 10,
        "journal_tip_sha256": runner.authority.ENTRY_JOURNAL_SHA256,
        "event": runner.authority.EVENT,
        "target": dict(runner.authority.TARGET),
        "disposition": "nonterminal",
        "suffix_kinds": (),
        "suffix_classification": "clean",
        "active_write": False,
        "writer_state": None,
        "external_receipt_authenticated": True,
        "entry_ancestor_present": True,
        "temporal_pending_members": ("RK4-2049",),
        "terminal": False,
        "terminal_lock_present": False,
        "safe_to_restart": True,
    }
    values.update(changes)
    return runner.authority.RCV3BoundedBoundary(**values)


def _completed_result(**changes):
    values = {
        "state": "event_one_complete",
        "authorization_commit": runner.ORIGINAL_EXECUTION_COMMIT,
        "campaign_id": runner.authority.CAMPAIGN_ID,
        "plan_sha256": runner.ORIGINAL_PLAN_SHA256,
        "event": runner.authority.SUCCESSOR_EVENT,
        "target": dict(runner.authority.SUCCESSOR_TARGET),
        "candidate_branch_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
    }
    values.update(changes)
    return values


def _runtime_checkpoint(**changes):
    values = {
        "event": runner.authority.EVENT,
        "target": dict(runner.authority.TARGET),
        "disposition": "nonterminal",
    }
    values.update(changes)
    return SimpleNamespace(**values)


class TDG8RCV3EventRunnerTests(unittest.TestCase):
    def test_execution_authority_reads_only_fixed_bytes_and_checks_typed_scope(self) -> None:
        receipt = _receipt()
        with (
            patch.object(
                runner,
                "_read_authority_leaf",
                side_effect=(CONFIG_RAW, RESULT_RAW),
            ) as read,
            patch.object(
                runner.authority,
                "authorize_execution",
                return_value=receipt,
            ) as authorize,
        ):
            observed = runner._execution_authority(ROOT, AUTHORITY_COMMIT)
        self.assertIs(observed, receipt)
        self.assertEqual(read.call_count, 2)
        authorize.assert_called_once_with(
            ROOT, CONFIG_RAW, RESULT_RAW, AUTHORITY_COMMIT
        )

        with (
            patch.object(
                runner,
                "_read_authority_leaf",
                side_effect=(CONFIG_RAW, RESULT_RAW),
            ),
            patch.object(
                runner.authority,
                "authorize_execution",
                return_value=replace(receipt, projection_id="foreign"),
            ),
        ):
            with self.assertRaisesRegex(
                runner.TDG8RCV3EventRunnerError,
                "authority/typed_scope_rejected",
            ):
                runner._execution_authority(ROOT, AUTHORITY_COMMIT)

    def test_original_plan_always_uses_immutable_execution_commit(self) -> None:
        base = object()
        plan = _plan()
        with (
            patch.object(
                runner,
                "construct_first_event",
                return_value=base,
            ) as construct,
            patch.object(runner, "replace", return_value=plan) as bind_campaign,
        ):
            self.assertIs(runner._original_plan(ROOT), plan)
        construct.assert_called_once_with(
            ROOT,
            authorization_commit=runner.ORIGINAL_EXECUTION_COMMIT,
        )
        bind_campaign.assert_called_once_with(
            base,
            campaign_id=runner.authority.CAMPAIGN_ID,
        )
        self.assertEqual(
            runner.ORIGINAL_EXECUTION_COMMIT,
            "6df67967e6e0bc63eca4dc6951a60d3787368f26",
        )
        self.assertNotEqual(runner.ORIGINAL_EXECUTION_COMMIT, AUTHORITY_COMMIT)

    def test_live_original_plan_reconstructs_the_frozen_digest(self) -> None:
        plan = runner._original_plan(ROOT.resolve())
        self.assertEqual(plan.sha256, runner.ORIGINAL_PLAN_SHA256)
        self.assertEqual(
            plan.authorization_commit,
            "6df67967e6e0bc63eca4dc6951a60d3787368f26",
        )
        self.assertEqual(plan.campaign_id, runner.authority.CAMPAIGN_ID)

    def test_destination_is_the_fixed_rcv3_projection(self) -> None:
        self.assertEqual(
            runner.EXPECTED_DESTINATION_PATH.as_posix(),
            "runs/fgc-2-sf1/tdg8-rcv3/calibration",
        )
        with patch.object(
            runner.authority,
            "DESTINATION_PATH",
            "runs/foreign/calibration",
        ):
            with self.assertRaisesRegex(
                runner.TDG8RCV3EventRunnerError,
                "authority/destination_scope",
            ):
                runner._destination(ROOT.resolve())

    def test_status_authenticates_descendant_without_opening_a_writer_or_store(self) -> None:
        boundary = _boundary(
            suffix_kinds=("tdg6_rejection",),
            suffix_classification="tdg6_rejection",
            safe_to_restart=False,
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(
                runner.authority,
                "inspect_descendant",
                return_value=boundary,
            ) as inspect_descendant,
            patch.object(runner, "HLT16CampaignStore") as store,
            patch.object(
                runner.runtime,
                "recover_persisted_generation_one",
            ) as restore,
            patch.object(runner, "build_static_gr0_shells") as shells,
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            result = runner.status(ROOT, authority_commit=AUTHORITY_COMMIT)
        inspect_descendant.assert_called_once()
        store.assert_not_called()
        restore.assert_not_called()
        shells.assert_not_called()
        advance.assert_not_called()
        self.assertEqual(result["state"], "temporal_retry_pending")
        self.assertFalse(result["boundary"]["safe_to_restart"])
        self.assertFalse(result["candidate_branch_opened"])

    def test_run_restores_persisted_state_and_advances_only_the_existing_event(self) -> None:
        plan = _plan()
        checkpoint = _runtime_checkpoint()
        store = SimpleNamespace()
        templates = {key: object() for key in runner.MEMBER_KEYS}
        completed = _completed_result()
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=plan),
            patch.object(
                runner.authority,
                "inspect_descendant",
                return_value=_boundary(),
            ),
            patch.object(runner, "_destination", return_value=ROOT / "projection"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(
                runner.runtime,
                "recover_persisted_generation_one",
                return_value=checkpoint,
            ) as restore,
            patch.object(
                runner,
                "build_static_gr0_shells",
                return_value=templates,
            ) as shells,
            patch.object(
                runner.event1,
                "_advance_authenticated_event",
                return_value=completed,
            ) as advance,
        ):
            result = runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        restore.assert_called_once_with(plan, store)
        shells.assert_called_once_with(ROOT)
        advance.assert_called_once_with(plan, store, checkpoint, templates)
        self.assertEqual(result["state"], "event_one_complete")
        self.assertEqual(result["result"]["event"], 24)
        self.assertFalse(result["candidate_branch_opened"])
        self.assertFalse(result["calibration_result_earned"])
        self.assertFalse(result["physical_result_earned"])

    def test_completed_boundary_returns_without_store_or_pde_path(self) -> None:
        boundary = _boundary(
            event=runner.authority.SUCCESSOR_EVENT,
            target=dict(runner.authority.SUCCESSOR_TARGET),
            disposition="event_complete",
            temporal_pending_members=(),
            safe_to_restart=False,
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(
                runner.authority,
                "inspect_descendant",
                return_value=boundary,
            ),
            patch.object(runner, "HLT16CampaignStore") as store,
            patch.object(
                runner.runtime,
                "recover_persisted_generation_one",
            ) as restore,
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            result = runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        self.assertEqual(result["state"], "event_one_complete")
        self.assertFalse(result["output_created"])
        store.assert_not_called()
        restore.assert_not_called()
        advance.assert_not_called()

    def test_unsafe_restart_fails_before_store_or_event_engine(self) -> None:
        boundary = _boundary(
            active_write=True,
            writer_state="live_lock",
            safe_to_restart=False,
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(
                runner.authority,
                "inspect_descendant",
                return_value=boundary,
            ),
            patch.object(runner, "HLT16CampaignStore") as store,
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            with self.assertRaisesRegex(
                runner.TDG8RCV3EventRunnerError,
                "writer/restart_not_authorized",
            ):
                runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        store.assert_not_called()
        advance.assert_not_called()

    def test_store_movement_after_inspection_cannot_escape_one_event_scope(self) -> None:
        moved = SimpleNamespace(
            event=25,
            target={"rational": "13/8", "binary64_hex": (13 / 8).hex()},
            disposition="nonterminal",
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner.authority, "inspect_descendant", return_value=_boundary()),
            patch.object(runner, "_destination", return_value=ROOT / "projection"),
            patch.object(runner, "HLT16CampaignStore", return_value=SimpleNamespace()),
            patch.object(
                runner.runtime, "recover_persisted_generation_one", return_value=moved,
            ),
            patch.object(runner, "build_static_gr0_shells") as shells,
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            with self.assertRaisesRegex(
                runner.TDG8RCV3EventRunnerError,
                "scope/runtime_checkpoint_event",
            ):
                runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        shells.assert_not_called()
        advance.assert_not_called()

    def test_post_recovery_scope_accepts_only_event23_or_exact_event24(self) -> None:
        runner._require_runtime_checkpoint_scope(
            SimpleNamespace(
                event=23, target=dict(runner.authority.TARGET), disposition="nonterminal",
            )
        )
        runner._require_runtime_checkpoint_scope(
            SimpleNamespace(
                event=24,
                target=dict(runner.authority.SUCCESSOR_TARGET),
                disposition="event_complete",
            )
        )
        with self.assertRaises(runner.TDG8RCV3EventRunnerError):
            runner._require_runtime_checkpoint_scope(
                SimpleNamespace(
                    event=24,
                    target=dict(runner.authority.SUCCESSOR_TARGET),
                    disposition="nonterminal",
                )
            )

    def test_recovery_failure_has_typed_owner_and_sanitized_detail(self) -> None:
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(
                runner.authority,
                "inspect_descendant",
                return_value=_boundary(),
            ),
            patch.object(runner, "_destination", return_value=ROOT / "projection"),
            patch.object(
                runner,
                "HLT16CampaignStore",
                return_value=SimpleNamespace(),
            ),
            patch.object(
                runner.runtime,
                "recover_persisted_generation_one",
                return_value=_runtime_checkpoint(),
            ),
            patch.object(
                runner,
                "build_static_gr0_shells",
                return_value={key: object() for key in runner.MEMBER_KEYS},
            ),
            patch.object(
                runner.event1,
                "_advance_authenticated_event",
                side_effect=runner.recovery.HLT16RecoveryError(
                    f"suffix failed at {ROOT}/private-detail"
                ),
            ),
        ):
            with self.assertRaises(runner.TDG8RCV3EventRunnerError) as caught:
                runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        record = caught.exception.as_record()
        self.assertEqual(record["error_owner"], "recovery")
        self.assertEqual(record["error_code"], "metadata_reconciliation_rejected")
        self.assertIn("<repo>/private-detail", record["error_detail"])
        self.assertNotIn(str(ROOT), record["error_detail"])
        self.assertFalse(record["physical_result_earned"])

    def test_result_scope_rejects_promotion_new_authority_or_extra_event(self) -> None:
        variants = (
            {**_completed_result(), "candidate_branch_opened": True},
            {**_completed_result(), "calibration_result_earned": True},
            {**_completed_result(), "physical_result_earned": True},
            {**_completed_result(), "authorization_commit": AUTHORITY_COMMIT},
            {**_completed_result(), "event": 25},
        )
        for result in variants:
            with self.subTest(result=result):
                with self.assertRaises(runner.TDG8RCV3EventRunnerError):
                    runner._validate_run_result(result)

        terminal = _completed_result(
            state="scientific_terminal",
            event=runner.authority.EVENT,
            target=dict(runner.authority.TARGET),
        )
        self.assertIs(runner._validate_run_result(terminal), terminal)

    def test_runner_has_no_gen0_bootstrap_candidate_calibration_or_direct_pde_route(self) -> None:
        source = inspect.getsource(runner)
        tree = ast.parse(source)
        imported = {
            node.module or ""
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        }
        for name in imported:
            lowered = name.lower()
            self.assertNotIn("bootstrap", lowered)
            self.assertNotIn("candidate", lowered)
            self.assertNotIn("calibration", lowered)
            self.assertNotIn("sgbl", lowered)
            self.assertNotIn("fgcqr", lowered)
        for forbidden in (
            "reconstruct_gr0_members",
            "initialize_or_recover_generation_one",
            "install_generation_zero_prefix",
            "authenticate_and_reconstruct",
            "run_fgc_gr0_calibration",
            "execute_one_attempt",
        ):
            self.assertNotIn(forbidden, source)
        self.assertIn("recover_persisted_generation_one", source)
        self.assertIn("_advance_authenticated_event", source)
        main_source = inspect.getsource(runner.main)
        for forbidden_option in ("--root", "--destination", "--store-root"):
            self.assertNotIn(forbidden_option, main_source)


if __name__ == "__main__":
    unittest.main()
