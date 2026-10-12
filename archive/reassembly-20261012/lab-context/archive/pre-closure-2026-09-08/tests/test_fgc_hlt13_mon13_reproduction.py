from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import reproduce_fgc_hlt13_mon13 as hlt13  # noqa: E402


class HLT13MON13ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = hlt13.load_canonical_result()

    def test_canonical_record_reproduces_as_synthetic_only(self) -> None:
        hlt13.verify_canonical()
        self.assertEqual(self.record, hlt13.record())
        self.assertEqual(self.record["artifact_id"], "FGC-1-HLT13-MON13")
        self.assertTrue(self.record["gate_status"]["PROTO15_successor_runtime_implemented"])
        self.assertTrue(self.record["gate_status"]["PROTO15_runtime_synthetic_qualification_passed"])
        self.assertFalse(self.record["gate_status"]["PROTO15_pretrajectory_runtime_authorized"])
        self.assertFalse(self.record["gate_status"]["PROTO15_fresh_GR0_dynamic_calibration_authorized"])

    def test_compact_lineage_is_bound_at_f381_and_live(self) -> None:
        lineage = self.record["artifact_payload"]["immutable_lineage"]
        self.assertEqual(lineage["evidence_checkpoint_commit"], "f3812d0b93e318890e973b4b21f681dbde8f437e")
        self.assertTrue(lineage["checkpoint_is_ancestor_of_HEAD"])
        self.assertEqual(
            [item["path"] for item in lineage["tracked_compact_records"]],
            list(hlt13.COMPACT_LINEAGE),
        )
        self.assertTrue(all(item["matches_checkpoint_and_worktree"] for item in lineage["tracked_compact_records"]))

    def test_runtime_summary_covers_the_frozen_durability_contract(self) -> None:
        payload = self.record["artifact_payload"]
        summary = payload["synthetic_qualification_summary"]
        self.assertTrue(summary["synthetic_only"])
        self.assertTrue(summary["six_member_contract"])
        self.assertTrue(summary["checkpoint_fork_rejected"])
        self.assertTrue(summary["terminal_cursor_set_preserved"])
        self.assertTrue(summary["terminal_ledger_embedded"])
        self.assertTrue(summary["exact_rejection_suffix_recovery_only"])
        self.assertTrue(summary["uncheckpointed_accepted_or_rollback_suffix_invalid"])
        self.assertTrue(summary["accepted_pending_retry_plan_bound"])
        self.assertTrue(summary["rollback_new_identities"])
        self.assertTrue(summary["rollback_journal_hashes_bound_by_record_envelope"])
        self.assertTrue(summary["old_or_new_replace_crash_images_exercised"])
        self.assertTrue(summary["hash_integrity_requires_external_trusted_genesis"])
        self.assertFalse(summary["wholesale_campaign_replacement_authenticated"])
        self.assertEqual(summary["fault_hook_windows"], payload["runtime_binding"]["named_fault_hooks"])
        self.assertTrue(summary["all_fault_hooks_exercised"])
        for key, value in hlt13.EXECUTED_QUALIFICATION.items():
            self.assertEqual(summary[key], value, key)
        self.assertEqual(
            summary["qualification_case_count"],
            summary["nominal_case_count"] + summary["mutation_case_count"] + summary["fault_case_count"],
        )
        self.assertEqual(
            summary["fault_recovery_image_count"],
            summary["fault_case_count"] + summary["ambiguous_replace_crash_case_count"],
        )

    def test_namespace_and_historical_nonclaim_boundary_are_public(self) -> None:
        payload = self.record["artifact_payload"]
        namespace = payload["namespace_precondition"]
        self.assertTrue(namespace["both_new_namespaces_absent"])
        self.assertTrue(namespace["qualification_created_no_namespace"])
        self.assertTrue(all(not (ROOT / item["path"]).exists() for item in namespace["records"]))
        historical = payload["historical_boundary"]
        self.assertFalse(historical["historical_PROTO14_raw_bundle_opened"])
        self.assertTrue(historical["PROTO14_invalid_runtime_history_preserved"])
        self.assertTrue(historical["PROTO14_terminal_checkpoint_may_not_resume"])
        self.assertEqual(self.record["nonclaims"], {key: False for key in hlt13.FALSE_CLAIMS})

    def test_claim_and_configuration_shapes_are_closed(self) -> None:
        config = hlt13.load_config()
        self.assertEqual(set(config["claims"]), hlt13.TRUE_CLAIMS | hlt13.FALSE_CLAIMS)
        self.assertEqual(config["synthetic_qualification"]["member_keys_in_canonical_order"], hlt13.MEMBERS)
        self.assertEqual(config["runtime_binding"]["journal_record_envelope"], ["record_kind", "journal_parent_sha256", "payload", "record_sha256"])
        self.assertEqual(config["runtime_binding"]["accepted_fine_commit_record_payload_fields"], hlt13.ACCEPTED_RECORD_FIELDS)
        self.assertEqual(config["runtime_binding"]["rollback_protocol_journal_field_locations"], hlt13.ROLLBACK_JOURNAL_FIELD_LOCATIONS)
        self.assertTrue(config["runtime_binding"]["hash_integrity_requires_external_trusted_genesis"])
        self.assertFalse(config["runtime_binding"]["wholesale_campaign_replacement_authenticated"])
        self.assertEqual(len(config["runtime_binding"]["named_fault_hooks"]), 10)
        self.assertEqual(
            {key: config["synthetic_qualification"][key] for key in hlt13.EXECUTED_QUALIFICATION},
            hlt13.EXECUTED_QUALIFICATION,
        )
        altered = deepcopy(config)
        altered["claims"]["PROTO15_fresh_GR0_dynamic_calibration_authorized"] = True
        self.assertNotEqual(altered["claims"], self.record["gate_status"])
        altered = deepcopy(config)
        altered["runtime_binding"]["same_generation_checkpoint_fork_is_invalid"] = False
        self.assertNotEqual(altered["runtime_binding"], self.record["artifact_payload"]["runtime_binding"])

    def test_verify_entrypoint_does_not_write_or_launch(self) -> None:
        completed = subprocess.run(
            [sys.executable, "scripts/reproduce_fgc_hlt13_mon13.py", "--verify"],
            cwd=ROOT, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("verified", completed.stdout)
        self.assertIn("results/fgc-1-hlt13-mon13.json", completed.stdout)

    def test_focused_repository_checker_does_not_traverse_historical_campaigns(self) -> None:
        completed = subprocess.run(
            [sys.executable, "scripts/check_repo.py", "--only-hlt13-mon13"],
            cwd=ROOT, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("PASS FGC, NGS, DSF, and epistemic alignment", completed.stdout)
        self.assertIn("PASS FGC-1-HLT13-MON13 compact certificate", completed.stdout)
        self.assertNotIn("CAL11-PREF22", completed.stdout)
        self.assertNotIn("machine-readable result", completed.stdout)
        for relative in (
            "runs/fgc-2-sf1/proto15/calibration",
            "runs/fgc-2-sf1/proto15/holdout",
        ):
            self.assertFalse((ROOT / relative).exists())


if __name__ == "__main__":
    unittest.main()
