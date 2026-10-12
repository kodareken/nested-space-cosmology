"""Focused orchestration tests for the generic bounded TDG8 retry runner."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from scripts import run_fgc_tdg8_bounded_retry as runner


ROOT = Path(__file__).resolve().parents[1]
AUTHORITY_COMMIT = "a" * 40


def _receipt():
    return SimpleNamespace(
        authority_commit=AUTHORITY_COMMIT,
        original_execution_commit=runner.authority.ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=runner.authority.ORIGINAL_PLAN_SHA256,
        campaign_id=runner.authority.CAMPAIGN_ID,
        destination_path=runner.authority.DESTINATION_PATH,
    )


def _plan():
    return SimpleNamespace(
        sha256=runner.authority.ORIGINAL_PLAN_SHA256,
        authorization_commit=runner.authority.ORIGINAL_EXECUTION_COMMIT,
        protocol_artifact_id=runner.authority.TARGET_PROTOCOL,
        campaign_id=runner.authority.CAMPAIGN_ID,
        branch="GR-0",
        amplitude="3",
        common_event_index=23,
        start_time=SimpleNamespace(rational="23/16"),
        target_time=SimpleNamespace(rational="3/2"),
        members=tuple(SimpleNamespace(member_key=key) for key in runner.MEMBER_KEYS),
    )


def _boundary(**changes):
    value = dict(
        checkpoint_generation=9,
        checkpoint_sha256=runner.authority.FIRST_RECOVERY_CHECKPOINT_SHA256,
        journal_sequence=10,
        journal_tip_sha256=runner.authority.ENTRY_TRANSITION_SHA256,
        event=23,
        target=dict(runner.authority.TARGET),
        disposition="nonterminal",
        suffix_kinds=(),
        suffix_classification="clean",
        active_write=False,
        writer_state=None,
        entry_ancestor_present=True,
        first_recovery_ancestor_present=True,
        temporal_pending_members=(runner.authority.ENTRY_MEMBER_KEY,),
        terminal=False,
        terminal_lock_present=False,
        safe_to_restart=True,
        recoverable_under_rcv2=True,
    )
    value.update(changes)
    return runner.authority.BoundedRetryBoundary(**value)


class TDG8BoundedRetryRunnerTests(unittest.TestCase):
    def test_status_authenticates_without_recovery_or_advance(self) -> None:
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
                runner.authority, "inspect_descendant", return_value=boundary,
            ),
            patch.object(runner.recovery, "reconcile") as reconcile,
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            result = runner.status(ROOT, authority_commit=AUTHORITY_COMMIT)
        self.assertEqual(result["state"], "temporal_retry_pending")
        self.assertFalse(result["boundary"]["safe_to_restart"])
        reconcile.assert_not_called()
        advance.assert_not_called()

    def test_run_uses_only_persisted_gr0_and_the_original_event(self) -> None:
        boundary = _boundary()
        plan = _plan()
        checkpoint = SimpleNamespace(event=23)
        store = SimpleNamespace()
        templates = {key: object() for key in runner.MEMBER_KEYS}
        completed = {
            "state": "event_one_complete",
            "campaign_id": runner.authority.CAMPAIGN_ID,
            "plan_sha256": runner.authority.ORIGINAL_PLAN_SHA256,
            "event": 24,
            "target": dict(runner.authority.SUCCESSOR_TARGET),
            "candidate_branch_opened": False,
            "calibration_result_earned": False,
            "physical_result_earned": False,
        }
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=plan),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(
                runner.authority, "inspect_descendant", return_value=boundary,
            ),
            patch.object(
                runner.runtime, "recover_persisted_generation_one",
                return_value=checkpoint,
            ) as restore,
            patch.object(
                runner, "build_static_gr0_shells", return_value=templates,
            ) as shells,
            patch.object(
                runner.event1, "_advance_authenticated_event",
                return_value=completed,
            ) as advance,
        ):
            result = runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        restore.assert_called_once_with(plan, store)
        shells.assert_called_once_with(ROOT)
        advance.assert_called_once_with(plan, store, checkpoint, templates)
        self.assertEqual(result["state"], "event_one_complete")
        self.assertFalse(result["candidate_branch_opened"])
        self.assertFalse(result["physical_result_earned"])

    def test_completed_boundary_returns_without_taking_a_writer(self) -> None:
        boundary = _boundary(
            event=24,
            target=dict(runner.authority.SUCCESSOR_TARGET),
            disposition="event_complete",
            temporal_pending_members=(),
            recoverable_under_rcv2=False,
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=SimpleNamespace()),
            patch.object(
                runner.authority, "inspect_descendant", return_value=boundary,
            ),
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            result = runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        self.assertEqual(result["state"], "event_one_complete")
        self.assertFalse(result["output_created"])
        advance.assert_not_called()

    def test_recovery_failure_retains_owning_layer_and_sanitized_detail(self) -> None:
        boundary = _boundary()
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=SimpleNamespace()),
            patch.object(
                runner.authority, "inspect_descendant", return_value=boundary,
            ),
            patch.object(
                runner.runtime, "recover_persisted_generation_one",
                return_value=SimpleNamespace(),
            ),
            patch.object(
                runner, "build_static_gr0_shells",
                return_value={key: object() for key in runner.MEMBER_KEYS},
            ),
            patch.object(
                runner.event1, "_advance_authenticated_event",
                side_effect=runner.recovery.HLT16RecoveryError(
                    f"exact suffix failed at {ROOT}/secret-like-detail"
                ),
            ),
        ):
            with self.assertRaises(runner.TDG8BoundedRetryRunnerError) as caught:
                runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        record = caught.exception.as_record()
        self.assertEqual(record["error_owner"], "recovery")
        self.assertEqual(record["error_code"], "metadata_reconciliation_rejected")
        self.assertIn("<repo>/secret-like-detail", record["error_detail"])
        self.assertNotIn(str(ROOT), record["error_detail"])
        self.assertFalse(record["physical_result_earned"])

    def test_runner_has_no_gen0_candidate_or_calibration_route(self) -> None:
        source = (ROOT / "scripts/run_fgc_tdg8_bounded_retry.py").read_text(
            encoding="utf-8",
        )
        self.assertNotIn("reconstruct_gr0_members", source)
        self.assertNotIn("initialize_or_recover_generation_one", source)
        self.assertNotIn("install_generation_zero_prefix", source)
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertNotIn("SGB-L", source)
        self.assertNotIn("FGC-QR", source)

    def test_result_scope_rejects_hidden_promotion(self) -> None:
        promoted = {
            "state": "event_one_complete",
            "campaign_id": runner.authority.CAMPAIGN_ID,
            "plan_sha256": runner.authority.ORIGINAL_PLAN_SHA256,
            "event": 24,
            "target": dict(runner.authority.SUCCESSOR_TARGET),
            "candidate_branch_opened": True,
            "calibration_result_earned": False,
            "physical_result_earned": False,
        }
        with self.assertRaisesRegex(
            runner.TDG8BoundedRetryRunnerError, "scope/result_promotion",
        ):
            runner._validate_run_result(promoted)

        incomplete = {
            **promoted,
            "candidate_branch_opened": False,
            "state": "event_one_complete",
            "event": 23,
            "target": dict(runner.authority.TARGET),
        }
        with self.assertRaisesRegex(
            runner.TDG8BoundedRetryRunnerError, "scope/event_boundary",
        ):
            runner._validate_run_result(incomplete)


if __name__ == "__main__":
    unittest.main()
