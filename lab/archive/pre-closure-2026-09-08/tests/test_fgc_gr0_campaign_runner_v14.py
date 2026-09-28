from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import run_fgc_gr0_calibration_v14 as runner  # noqa: E402


class FGCGR0CampaignRunnerV14Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # Restoration validates the immutable ignored payloads and is the one
        # comparatively expensive fixture in this focused runner suite.
        cls.members = runner._restore_frozen_members({})

    def test_frozen_run_plan_is_strict_and_candidate_closed(self) -> None:
        plan = runner._validate_run_plan(runner.DEFAULT_PLAN)
        self.assertEqual(plan["scope"]["target_protocol"], "FGC-2-SF1-PROTO14")
        self.assertTrue(
            plan["temporal_admission"]["every_accepted_post_restart_macro_step_must_pass"]
        )
        self.assertTrue(plan["temporal_admission"]["sampled_64_history_diagnostic_only"])
        self.assertTrue(all(value is False for value in plan["claims"].values()))
        self.assertFalse((REPOSITORY / plan["provenance"]["output_root"]).exists())

    def test_restart_wraps_all_six_members_with_zero_tdg6_ledgers(self) -> None:
        members = self.members
        self.assertEqual(set(members), set(runner._MEMBER_KEYS))
        for member in members.values():
            self.assertEqual(member.time, runner.RESTART_TIME)
            self.assertEqual(member.temporal_ledger.last_accepted_time, runner.RESTART_TIME)
            self.assertEqual(member.temporal_ledger.accepted_macro_step_count, 0)
            self.assertEqual(member.temporal_ledger.cumulative_temporal_retry_count, 0)
            self.assertFalse(member.temporal_ledger.serialized_temporal_rejections)
            self.assertTrue(all(value == 0.0 for value in member.temporal_ledger.accumulated_debit_vector))

    def test_event_assessment_makes_sampled_history_diagnostic_only(self) -> None:
        members = self.members
        assessment = runner._event_assessment(
            runner._validate_run_plan(runner.DEFAULT_PLAN), members, runner.RESTART_TIME
        )
        self.assertFalse(assessment["sampled_temporal_history_diagnostic"]["available"])
        self.assertFalse(assessment["temporal_spectral_admission"]["sampled_history_rule_used_as_veto"])
        self.assertEqual(set(assessment["member_temporal_ledgers"]), set(members))
        self.assertFalse(assessment["accumulated_temporal_debit_is_Raychaudhuri_error_bound"])
        # The restart event has no post-restart admitted TDG6 macro step and
        # must therefore be unable to start a qualified-trapped streak.
        self.assertFalse(assessment["TDG6_ledger_consistent_at_common_event"])
        self.assertFalse(assessment["qualified_trapped_common_event"])

    def test_durable_retry_entry_has_common_event_attempt_identity(self) -> None:
        member = next(iter(self.members.values()))
        with tempfile.TemporaryDirectory() as directory:
            event_log = Path(directory) / "events.jsonl"
            sink = runner._durable_attempt_sink(
                event_log, attempt_id="proto14-common-0024", member=member,
                kind="TDG6_temporal",
            )
            sink({"event_type": "rejected_TDG6_temporal_admission"})
            record = json.loads(event_log.read_text(encoding="utf-8"))
        self.assertEqual(record["common_event_attempt_id"], "proto14-common-0024")
        self.assertEqual(record["member"], member.key)
        self.assertEqual(record["retry_kind"], "TDG6_temporal")
        self.assertTrue(record["classification_is_numerical_not_physical"])

    def test_atomic_checkpoint_round_trips_all_tdg6_extensions(self) -> None:
        manifest = {"campaign_id": "test-proto14", "schema_version": 1}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            event_log = root / "events.jsonl"
            event_log.write_bytes(b"")
            checkpoint = root / "latest-checkpoint.npz"
            metadata = runner._checkpoint_metadata(
                manifest=manifest, event_index=23, consecutive=0, event_log=event_log,
                members=self.members, terminal=False, terminal_result=None,
            )
            runner._write_checkpoint(
                checkpoint, metadata=metadata, members=self.members, event_log=event_log
            )
            restored_metadata, restored, payload = runner._restore_campaign_checkpoint(
                checkpoint, manifest=manifest
            )
        self.assertEqual(payload, b"")
        self.assertEqual(restored_metadata["TDG6_checkpoint_extension_present"], True)
        for key in self.members:
            self.assertEqual(
                restored[key].temporal_ledger, self.members[key].temporal_ledger
            )

    def test_terminal_checkpoint_publication_is_idempotently_recoverable(self) -> None:
        manifest = {"campaign_id": "test-proto14", "schema_version": 1}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            event_log = root / "events.jsonl"
            event_log.write_bytes(b"{}\n")
            result_path = root / "campaign-result.json"
            final_hashes = {
                key: runner.array_content_sha256(
                    member.state.u, member.state.p, member.state.q
                )
                for key, member in self.members.items()
            }
            result = {
                "schema_version": 2,
                "runner_id": runner.RUNNER_ID,
                "campaign_id": manifest["campaign_id"],
                "classification": "calibration_failed_no_eligible_GR0_case",
                "amplitude_records": [{
                    "last_completed_common_event_index": 23,
                    "consecutive_qualified_trapped_common_events": 0,
                    "member_final_state_hashes": final_hashes,
                }],
                "event_log": "runs/test/events.jsonl",
                "event_log_sha256": runner.sha256(event_log.read_bytes()).hexdigest(),
                "FGCQR_outcome_read": False,
                "SGBL_outcome_read": False,
                "holdout_execution_authorized": False,
                "retained_EFT_evolution_authorized": False,
                "mechanism_question_answered": False,
            }
            metadata = runner._checkpoint_metadata(
                manifest=manifest,
                event_index=23,
                consecutive=0,
                event_log=event_log,
                members=self.members,
                terminal=True,
                terminal_result=result,
            )
            checkpoint = root / "latest-checkpoint.npz"
            runner._write_checkpoint(
                checkpoint,
                metadata=metadata,
                members=self.members,
                event_log=event_log,
            )
            restored_metadata, _, restored_log = runner._restore_campaign_checkpoint(
                checkpoint,
                manifest=manifest,
            )

            # Crash window one: the terminal checkpoint exists but publication
            # has not happened. Resume must publish exactly its embedded result.
            published = runner._reconcile_terminal_publication(
                metadata=restored_metadata,
                manifest=manifest,
                event_log_payload=restored_log,
                event_log_relative="runs/test/events.jsonl",
                result_path=result_path,
            )
            self.assertEqual(published, runner.inherited._serial(result))
            self.assertTrue(result_path.is_file())

            # Crash window two: publication succeeded but the process did not
            # return. A second resume verifies and returns the same evidence.
            repeated = runner._reconcile_terminal_publication(
                metadata=restored_metadata,
                manifest=manifest,
                event_log_payload=restored_log,
                event_log_relative="runs/test/events.jsonl",
                result_path=result_path,
            )
            self.assertEqual(repeated, published)

    def test_terminal_publication_mismatch_and_nonterminal_result_fail_closed(self) -> None:
        manifest = {"campaign_id": "test-proto14", "schema_version": 1}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            event_log = root / "events.jsonl"
            event_log.write_bytes(b"{}\n")
            result_path = root / "campaign-result.json"
            result_path.write_text("{}\n", encoding="utf-8")
            nonterminal = runner._checkpoint_metadata(
                manifest=manifest,
                event_index=23,
                consecutive=0,
                event_log=event_log,
                members=self.members,
                terminal=False,
                terminal_result=None,
            )
            with self.assertRaisesRegex(ValueError, "nonterminal checkpoint"):
                runner._reconcile_terminal_publication(
                    metadata=nonterminal,
                    manifest=manifest,
                    event_log_payload=event_log.read_bytes(),
                    event_log_relative="runs/test/events.jsonl",
                    result_path=result_path,
                )

            final_hashes = {
                key: runner.array_content_sha256(
                    member.state.u, member.state.p, member.state.q
                )
                for key, member in self.members.items()
            }
            expected = {
                "schema_version": 2,
                "runner_id": runner.RUNNER_ID,
                "campaign_id": manifest["campaign_id"],
                "classification": "calibration_failed_no_eligible_GR0_case",
                "amplitude_records": [{
                    "last_completed_common_event_index": 23,
                    "consecutive_qualified_trapped_common_events": 0,
                    "member_final_state_hashes": final_hashes,
                }],
                "event_log": "runs/test/events.jsonl",
                "event_log_sha256": runner.sha256(event_log.read_bytes()).hexdigest(),
                "FGCQR_outcome_read": False,
                "SGBL_outcome_read": False,
                "holdout_execution_authorized": False,
                "retained_EFT_evolution_authorized": False,
                "mechanism_question_answered": False,
            }
            terminal = runner._checkpoint_metadata(
                manifest=manifest,
                event_index=23,
                consecutive=0,
                event_log=event_log,
                members=self.members,
                terminal=True,
                terminal_result=expected,
            )
            with self.assertRaisesRegex(ValueError, "published result differs"):
                runner._reconcile_terminal_publication(
                    metadata=terminal,
                    manifest=manifest,
                    event_log_payload=event_log.read_bytes(),
                    event_log_relative="runs/test/events.jsonl",
                    result_path=result_path,
                )

    def test_failed_common_constraint_cannot_become_eighth_qualified_event(self) -> None:
        originals = {key: member.temporal_ledger for key, member in self.members.items()}
        try:
            for member in self.members.values():
                member.temporal_ledger = replace(
                    member.temporal_ledger, accepted_macro_step_count=1
                )
            base = {
                "primary_common_event": SimpleNamespace(admission_passed=False),
                "comparator_common_event": SimpleNamespace(admission_passed=True),
                "temporal_spectral_admission": {"available": True},
                "trapped_assessment": SimpleNamespace(trapped_sign_passed=True),
            }
            with patch.object(runner.proto13_runner, "_event_assessment", return_value=base):
                assessment = runner._event_assessment({}, self.members, runner.RESTART_TIME)
            self.assertTrue(assessment["TDG6_ledger_consistent_at_common_event"])
            self.assertFalse(assessment["qualified_trapped_common_event"])
        finally:
            for key, member in self.members.items():
                member.temporal_ledger = originals[key]


if __name__ == "__main__":
    unittest.main()
