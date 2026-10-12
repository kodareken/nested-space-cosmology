from __future__ import annotations

from copy import deepcopy
import unittest

from scripts import reproduce_fgc_pro15_frz1 as pro15


class PRO15FreezeReproductionTests(unittest.TestCase):
    def test_canonical_record_reproduces_without_authorizing_execution(self) -> None:
        record = pro15.load_canonical_result()
        pro15.verify_canonical()
        self.assertEqual(record, pro15.record())
        self.assertEqual(record["artifact_id"], pro15.ARTIFACT_ID)
        self.assertTrue(record["gate_status"]["PROTO15_frozen"])
        self.assertFalse(record["gate_status"]["PROTO15_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(record["gate_status"]["FGCQR_holdout_execution_authorized"])

    def test_compact_restart_lineage_creates_cursor_origin_templates_only(self) -> None:
        members = pro15.load_canonical_result()["artifact_payload"]["restart_members"]
        self.assertEqual([member["key"] for member in members], ["RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"])
        for member in members:
            template = member["cursor_origin_template"]
            self.assertTrue(template["not_a_serialized_complete_cursor"])
            self.assertEqual(template["mode"], "FRESH_READY")
            self.assertIsNone(template["retry_successor_payload_or_none"])
            self.assertEqual(template["current_retry_count"], 0)
            self.assertEqual(template["accepted_boundary_time_rational"], "23/16")
            self.assertEqual(template["accepted_boundary_time_binary64_hex"], float(23 / 16).hex())
            self.assertEqual(template["accepted_state_sha256"], member["state_sha256"])
            self.assertEqual(template["journal_tip_sha256"], "0" * 64)
        self.assertTrue(all(member["TDG6_ledger_initialization"]["all_debits_exact_zero"] for member in members))

    def test_historical_invalid_boundary_and_namespace_absence_are_public(self) -> None:
        payload = pro15.load_canonical_result()["artifact_payload"]
        history = payload["historical_boundary"]
        namespaces = payload["namespace_precondition"]
        self.assertTrue(history["PROTO14_invalid_runtime_or_nonconverged_not_physics"])
        self.assertTrue(history["TDG7_stateless_adapter_is_not_campaign_cursor"])
        self.assertTrue(history["campaign_owned_cursor_required_before_calibration"])
        self.assertTrue(namespaces["both_new_namespaces_absent"])
        self.assertTrue(namespaces["freeze_created_no_namespace"])

    def test_retry_and_exhaustion_branches_are_distinct(self) -> None:
        lifecycle = pro15.load_canonical_result()["artifact_payload"]["lifecycle_contract"]
        self.assertEqual(lifecycle["nonexhausted_retry_transition_order"][:2], ["durable_TDG6_rejection_journal_record", "typed_TDG6_retry_required_raised"])
        self.assertEqual(lifecycle["exhaustion_transition_order"], ["durable_TDG6_rejection_journal_record", "typed_TDG6_retry_exhaustion_raised", "append_and_fsync_terminal_invalid_nonconverged_record", "checkpoint_terminal_state"])
        self.assertTrue(lifecycle["TDG6_retry_exhaustion_must_not_create_RETRY_PENDING"])
        self.assertTrue(lifecycle["TDG6_retry_exhaustion_is_invalid_nonconverged_not_physical"])
        self.assertIn("journal_tip_sha256", lifecycle["terminal_exhaustion_required_fields"])
        self.assertTrue(lifecycle["durable_journal_append_requires_one_canonical_record_flush_and_fsync_before_return"])
        self.assertTrue(lifecycle["atomic_checkpoint_requires_file_fsync_atomic_replace_and_parent_directory_fsync"])
        self.assertTrue(lifecycle["rollback_must_append_new_generation_binding_restored_snapshot_not_reuse_old_cursor_bytes"])
        self.assertEqual(
            lifecycle["journal_record_envelope_required_fields"],
            ["record_kind", "journal_parent_sha256", "payload", "record_sha256"],
        )
        self.assertTrue(lifecycle["cursor_transition_record_sha256_becomes_checkpoint_journal_tip"])
        self.assertTrue(
            lifecycle[
                "accepted_fine_state_TDG6_ledger_and_FRESH_READY_cursor_must_publish_as_one_atomic_checkpoint_generation"
            ]
        )
        self.assertTrue(lifecycle["crash_before_new_durable_record_restores_exact_last_durable_cursor"])
        self.assertTrue(lifecycle["checkpoint_must_bind_complete_six_member_cursor_ledger_state_sets_and_journal_tip"])
        self.assertTrue(lifecycle["common_event_rollback_is_one_campaign_record_over_complete_six_member_snapshot"])
        self.assertEqual(lifecycle["terminal_exhaustion_classification_exact"], "invalid_implementation_or_nonconverged_run")
        self.assertTrue(lifecycle["terminal_exhaustion_requires_terminal_lock_true"])
        self.assertTrue(lifecycle["campaign_owner_must_complete_retry_or_terminal_persistence_before_typed_outcome_escapes"])

    def test_embedded_mutations_all_fail_closed(self) -> None:
        controls = pro15.load_canonical_result()["artifact_payload"]["mutation_controls"]
        self.assertTrue(all(controls.values()))
        self.assertIn("exhaustion_cannot_retry", controls)
        self.assertIn("crash_recovery", controls)
        self.assertIn("terminal_exhaustion_fields", controls)
        self.assertIn("rollback_generation", controls)
        self.assertIn("journal_envelope", controls)
        self.assertIn("pre_record_crash_restore", controls)
        self.assertIn("member_set", controls)
        self.assertIn("retry_semantics", controls)
        self.assertIn("terminal_lock", controls)
        self.assertIn("recovery_order", controls)

    def test_config_mutations_differ_from_frozen_contract(self) -> None:
        config = pro15.load_config()
        for name, mutate in {
            "state": lambda v: v["immutable_lineage"].__setitem__("proto14_freeze_result_sha256", "0" * 64),
            "namespace": lambda v: v["namespace_precondition"].__setitem__("calibration_output_root", "runs/fgc-2-sf1/proto14/calibration"),
            "claim": lambda v: v["claims"].__setitem__("FGCQR_holdout_execution_authorized", True),
            "extra_root": lambda v: v.__setitem__("unexpected", True),
            "removed_proof_key": lambda v: v["proof_contract"].pop("adversarial_mutations_must_fail_closed"),
            "false_proof_flag": lambda v: v["proof_contract"].__setitem__("adversarial_mutations_must_fail_closed", False),
            "protocol_core_path": lambda v: v.__setitem__("protocol_core", "src/recursive_horizons/fgc/evolution/protocol_v14.py"),
        }.items():
            with self.subTest(name=name):
                attacked = deepcopy(config)
                mutate(attacked)
                with self.assertRaises(ValueError):
                    pro15.load_config_from_value(attacked)


if __name__ == "__main__":
    unittest.main()
