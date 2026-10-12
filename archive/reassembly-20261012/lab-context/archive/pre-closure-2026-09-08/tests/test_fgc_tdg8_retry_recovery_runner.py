"""Focused orchestration tests for exact TDG8 recovery and resume."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from scripts import run_fgc_tdg8_retry_recovery as runner


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
        members=tuple(
            SimpleNamespace(member_key=key) for key in runner.MEMBER_KEYS
        ),
    )


def _checkpoint(sha256: str, generation: int):
    return SimpleNamespace(
        sha256=sha256,
        generation=generation,
        event=23,
        target={"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"},
        disposition="nonterminal",
        journal_sequence=runner.authority.EXPECTED_JOURNAL_SEQUENCE,
        journal_tip_sha256=runner.authority.EXPECTED_JOURNAL_SHA256,
    )


class TDG8RetryRecoveryRunnerTests(unittest.TestCase):
    def test_recover_publishes_only_the_unique_metadata_edge(self) -> None:
        predecessor = _checkpoint(
            runner.authority.ANCHOR_CHECKPOINT_SHA256,
            runner.authority.ANCHOR_CHECKPOINT_GENERATION,
        )
        recovered = _checkpoint(
            runner.authority.EXPECTED_CHECKPOINT_SHA256,
            runner.authority.EXPECTED_CHECKPOINT_GENERATION,
        )
        store = SimpleNamespace(
            authenticated_snapshot=Mock(
                return_value=SimpleNamespace(checkpoint=predecessor)
            )
        )
        decision = runner.recovery.RecoveryDecision(
            disposition="reconciled_checkpoint",
            checkpoint_sha256=runner.authority.EXPECTED_CHECKPOINT_SHA256,
            published=True,
            detail="tdg6_rejection",
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(
                runner.authority,
                "inspect_recovery_boundary",
                return_value={"state": "generation7_rejection_only"},
            ) as inspect,
            patch.object(
                runner.event1, "_acquire_writer", return_value="writer-token"
            ) as acquire,
            patch.object(runner.recovery, "reconcile", return_value=decision) as reconcile,
            patch.object(
                runner,
                "_expected_generation_eight",
                side_effect=(recovered, recovered),
            ) as expected,
            patch.object(
                runner.event1, "_release_writer_if_closed", return_value=True
            ) as release,
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            result = runner.recover(ROOT, authority_commit=AUTHORITY_COMMIT)

        self.assertEqual(result["state"], "metadata_recovery_complete")
        self.assertFalse(result["physical_state_advanced"])
        self.assertFalse(result["candidate_branch_opened"])
        self.assertEqual(inspect.call_count, 2)
        acquire.assert_called_once_with(store, predecessor)
        reconcile.assert_called_once_with(store)
        self.assertEqual(expected.call_count, 2)
        release.assert_called_once_with(store, owner_token="writer-token")
        advance.assert_not_called()

    def test_resume_starts_from_exact_generation_eight_and_original_plan(self) -> None:
        recovered = _checkpoint(
            runner.authority.EXPECTED_CHECKPOINT_SHA256,
            runner.authority.EXPECTED_CHECKPOINT_GENERATION,
        )
        plan = _plan()
        store = SimpleNamespace(
            authenticated_snapshot=Mock(
                return_value=SimpleNamespace(checkpoint=recovered)
            )
        )
        completed = {
            "state": "event_one_complete",
            "campaign_id": runner.authority.CAMPAIGN_ID,
            "plan_sha256": runner.authority.ORIGINAL_PLAN_SHA256,
            "candidate_branch_opened": False,
            "calibration_result_earned": False,
            "physical_result_earned": False,
        }
        templates = {key: object() for key in runner.MEMBER_KEYS}
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=plan),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(runner, "_require_generation_eight_ancestor") as ancestor,
            patch.object(
                runner, "_expected_generation_eight", return_value=recovered
            ) as expected,
            patch.object(
                runner.runtime,
                "recover_persisted_generation_one",
                return_value=recovered,
            ) as persisted,
            patch.object(
                runner, "build_static_gr0_shells", return_value=templates
            ) as shells,
            patch.object(
                runner.event1,
                "_advance_authenticated_event",
                return_value=completed,
            ) as advance,
            patch.object(runner.recovery, "reconcile") as reconcile,
            patch.object(runner.authority, "inspect_recovery_boundary") as inspect,
        ):
            result = runner.resume(ROOT, authority_commit=AUTHORITY_COMMIT)

        ancestor.assert_called_once_with(store)
        expected.assert_called_once_with(store, require_no_writer=True)
        persisted.assert_called_once_with(plan, store)
        shells.assert_called_once_with(ROOT)
        advance.assert_called_once_with(plan, store, recovered, templates)
        reconcile.assert_not_called()
        inspect.assert_not_called()
        self.assertEqual(result["state"], "event_one_complete")
        self.assertFalse(result["candidate_branch_opened"])
        self.assertFalse(result["calibration_result_earned"])
        self.assertFalse(result["physical_result_earned"])

    def test_recovery_releases_a_clean_exact_publication_when_postcheck_rejects(
        self,
    ) -> None:
        predecessor = _checkpoint(
            runner.authority.ANCHOR_CHECKPOINT_SHA256,
            runner.authority.ANCHOR_CHECKPOINT_GENERATION,
        )
        store = SimpleNamespace(
            authenticated_snapshot=Mock(
                return_value=SimpleNamespace(checkpoint=predecessor)
            )
        )
        decision = runner.recovery.RecoveryDecision(
            disposition="reconciled_checkpoint",
            checkpoint_sha256=runner.authority.EXPECTED_CHECKPOINT_SHA256,
            published=True,
            detail="tdg6_rejection",
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(
                runner.authority,
                "inspect_recovery_boundary",
                return_value={"state": "generation7_rejection_only"},
            ),
            patch.object(
                runner.event1, "_acquire_writer", return_value="writer-token"
            ),
            patch.object(runner.recovery, "reconcile", return_value=decision),
            patch.object(
                runner,
                "_expected_generation_eight",
                side_effect=runner.TDG8RetryRecoveryRunnerError("postcheck rejected"),
            ),
            patch.object(
                runner.event1, "_release_writer_if_closed", return_value=True
            ) as release,
        ):
            with self.assertRaisesRegex(
                runner.TDG8RetryRecoveryRunnerError, "postcheck rejected"
            ):
                runner.recover(ROOT, authority_commit=AUTHORITY_COMMIT)
        release.assert_called_once_with(store, owner_token="writer-token")

    def test_recovery_completes_an_exact_sequence_eight_publication_cut(self) -> None:
        predecessor = _checkpoint(
            runner.authority.ANCHOR_CHECKPOINT_SHA256,
            runner.authority.ANCHOR_CHECKPOINT_GENERATION,
        )
        recovered = _checkpoint(
            runner.authority.EXPECTED_CHECKPOINT_SHA256,
            runner.authority.EXPECTED_CHECKPOINT_GENERATION,
        )
        store = SimpleNamespace(
            authenticated_snapshot=Mock(
                return_value=SimpleNamespace(checkpoint=predecessor)
            )
        )
        decision = runner.recovery.RecoveryDecision(
            disposition="reconciled_checkpoint",
            checkpoint_sha256=runner.authority.EXPECTED_CHECKPOINT_SHA256,
            published=True,
            detail="tdg6_rejection+cursor_transition",
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(
                runner.authority,
                "inspect_recovery_boundary",
                return_value={"state": "generation7_transition_published"},
            ),
            patch.object(
                runner.event1, "_acquire_writer", return_value="writer-token"
            ),
            patch.object(runner.recovery, "reconcile", return_value=decision) as reconcile,
            patch.object(
                runner,
                "_expected_generation_eight",
                side_effect=(recovered, recovered),
            ),
            patch.object(
                runner.event1, "_release_writer_if_closed", return_value=True
            ),
        ):
            result = runner.recover(ROOT, authority_commit=AUTHORITY_COMMIT)
        reconcile.assert_called_once_with(store)
        self.assertEqual(result["checkpoint_sha256"], recovered.sha256)
        self.assertFalse(result["physical_state_advanced"])

    def test_recovery_closes_a_stale_writer_after_generation_eight_cut(self) -> None:
        recovered = _checkpoint(
            runner.authority.EXPECTED_CHECKPOINT_SHA256,
            runner.authority.EXPECTED_CHECKPOINT_GENERATION,
        )
        store = SimpleNamespace(
            authenticated_snapshot=Mock(
                return_value=SimpleNamespace(checkpoint=recovered)
            )
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(
                runner.authority,
                "inspect_recovery_boundary",
                return_value={"state": "generation8_checkpoint_published"},
            ),
            patch.object(
                runner.event1, "_acquire_writer", return_value="writer-token"
            ),
            patch.object(runner.recovery, "reconcile") as reconcile,
            patch.object(
                runner,
                "_expected_generation_eight",
                side_effect=(recovered, recovered),
            ),
            patch.object(
                runner.event1, "_release_writer_if_closed", return_value=True
            ) as release,
        ):
            result = runner.recover(ROOT, authority_commit=AUTHORITY_COMMIT)
        reconcile.assert_not_called()
        release.assert_called_once_with(store, owner_token="writer-token")
        self.assertEqual(result["checkpoint_sha256"], recovered.sha256)
        self.assertFalse(result["physical_state_advanced"])

    def test_resume_rejects_generation_seven_before_any_event_advance(self) -> None:
        predecessor = _checkpoint(
            runner.authority.ANCHOR_CHECKPOINT_SHA256,
            runner.authority.ANCHOR_CHECKPOINT_GENERATION,
        )
        store = SimpleNamespace(
            authenticated_snapshot=Mock(
                return_value=SimpleNamespace(checkpoint=predecessor)
            )
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(runner, "_require_generation_eight_ancestor") as ancestor,
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            with self.assertRaisesRegex(
                runner.TDG8RetryRecoveryRunnerError,
                "metadata recovery must complete",
            ):
                runner.resume(ROOT, authority_commit=AUTHORITY_COMMIT)
        ancestor.assert_not_called()
        advance.assert_not_called()

    def test_runner_has_no_generation_zero_or_destination_absence_route(self) -> None:
        source = (ROOT / "scripts/run_fgc_tdg8_retry_recovery.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("reconstruct_gr0_members", source)
        self.assertNotIn("install_generation_zero_prefix", source)
        self.assertNotIn("destination_state", source)
        self.assertNotIn("SGB-L", source.replace("No mode can import a candidate branch, reopen", ""))
        self.assertNotIn("FGC-QR", source.replace("No mode can import a candidate branch, reopen", ""))

    def test_resume_refuses_a_descendant_beyond_the_one_event_scope(self) -> None:
        foreign = _checkpoint("b" * 64, 12)
        foreign.event = 25
        foreign.target = {
            "rational": "13/8",
            "binary64_hex": "0x1.a000000000000p+0",
        }
        store = SimpleNamespace(
            authenticated_snapshot=Mock(
                return_value=SimpleNamespace(checkpoint=foreign)
            )
        )
        with (
            patch.object(runner, "_fixed_root", return_value=ROOT),
            patch.object(runner, "_execution_authority", return_value=_receipt()),
            patch.object(runner, "_original_plan", return_value=_plan()),
            patch.object(runner, "_destination", return_value=ROOT / "runs"),
            patch.object(runner, "HLT16CampaignStore", return_value=store),
            patch.object(runner, "_require_generation_eight_ancestor") as ancestor,
            patch.object(runner.event1, "_advance_authenticated_event") as advance,
        ):
            with self.assertRaisesRegex(
                runner.TDG8RetryRecoveryRunnerError,
                "escaped the one-event scope",
            ):
                runner.resume(ROOT, authority_commit=AUTHORITY_COMMIT)
        ancestor.assert_not_called()
        advance.assert_not_called()


if __name__ == "__main__":
    unittest.main()
