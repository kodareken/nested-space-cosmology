from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution.protocol_v15 import (
    PROTO15_LIFECYCLE,
    PROTO15_RESTART_INPUTS,
    validate_sf1_protocol_v15,
)


CONFIG = ROOT / "configs/fgc/fgc-2-sf1-protocol-v15.toml"


class FGCProtocolV15Tests(unittest.TestCase):
    def setUp(self) -> None:
        with CONFIG.open("rb") as handle:
            self.protocol = tomllib.load(handle)

    def test_freeze_validates_and_does_not_authorize_execution(self) -> None:
        result = validate_sf1_protocol_v15(self.protocol)
        self.assertEqual(result["artifact_id"], "FGC-2-SF1-PROTO15")
        self.assertEqual(result["protocol_version"], 15)
        self.assertTrue(result["campaign_owned_cursor_contract_defined"])
        self.assertFalse(result["production_trajectory_authorized"])
        self.assertFalse(result["historical_PROTO14_result_reclassified"])
        self.assertFalse(result["claims"]["PROTO15_fresh_GR0_dynamic_calibration_authorized"])

    def test_restart_inputs_and_lifecycle_are_exact(self) -> None:
        self.assertEqual(self.protocol["restart_inputs"], PROTO15_RESTART_INPUTS)
        self.assertEqual(self.protocol["lifecycle"], PROTO15_LIFECYCLE)
        self.assertEqual(len(self.protocol["restart_inputs"]["members"]), 6)
        self.assertEqual(self.protocol["restart_inputs"]["members"][0]["method"], "RK4")
        self.assertEqual(self.protocol["restart_inputs"]["members"][0]["point_count"], 2049)
        self.assertEqual(self.protocol["restart_inputs"]["members"][0]["step_index"], 342)
        self.assertEqual(self.protocol["restart_inputs"]["members"][0]["transaction_serial"], 1710)
        self.assertEqual(self.protocol["lifecycle"]["cursor_modes"], ["FRESH_READY", "RETRY_PENDING"])
        self.assertTrue(self.protocol["lifecycle"]["retry_must_be_durable_before_execution"])
        self.assertTrue(self.protocol["lifecycle"]["common_event_abort_must_restore_complete_six_member_snapshot_and_append_explicit_rollback"])
        self.assertEqual(
            self.protocol["lifecycle"]["nonexhausted_retry_transition_order"][:2],
            ["durable_TDG6_rejection_journal_record", "typed_TDG6_retry_required_raised"],
        )
        self.assertEqual(
            self.protocol["lifecycle"]["exhaustion_transition_order"][:2],
            ["durable_TDG6_rejection_journal_record", "typed_TDG6_retry_exhaustion_raised"],
        )
        self.assertTrue(self.protocol["lifecycle"]["TDG6_retry_exhaustion_must_not_create_RETRY_PENDING"])
        self.assertTrue(self.protocol["lifecycle"]["TDG6_retry_exhaustion_is_invalid_nonconverged_not_physical"])
        self.assertIn("journal_tip_sha256", self.protocol["lifecycle"]["cursor_required_fields"])
        self.assertTrue(self.protocol["lifecycle"]["time_identity_fields_require_rational_and_binary64_hex"])
        self.assertEqual(
            self.protocol["lifecycle"]["cursor_chain_root_parent_sha256"],
            "0" * 64,
        )
        self.assertTrue(
            self.protocol["lifecycle"][
                "rollback_must_append_new_generation_binding_restored_snapshot_not_reuse_old_cursor_bytes"
            ]
        )
        self.assertTrue(
            self.protocol["lifecycle"][
                "atomic_checkpoint_requires_file_fsync_atomic_replace_and_parent_directory_fsync"
            ]
        )
        self.assertEqual(
            self.protocol["lifecycle"]["journal_record_envelope_required_fields"],
            ["record_kind", "journal_parent_sha256", "payload", "record_sha256"],
        )
        self.assertEqual(
            self.protocol["lifecycle"]["cursor_journal_tip_semantics"],
            "hash_of_latest_durable_journal_record_preceding_cursor_transition",
        )
        self.assertTrue(
            self.protocol["lifecycle"][
                "cursor_transition_record_sha256_becomes_checkpoint_journal_tip"
            ]
        )
        self.assertTrue(
            self.protocol["lifecycle"][
                "accepted_fine_state_TDG6_ledger_and_FRESH_READY_cursor_must_publish_as_one_atomic_checkpoint_generation"
            ]
        )
        self.assertTrue(
            self.protocol["lifecycle"][
                "crash_before_new_durable_record_restores_exact_last_durable_cursor"
            ]
        )
        self.assertTrue(self.protocol["lifecycle"]["FRESH_READY_requires_retry_successor_payload_json_null"])
        self.assertTrue(self.protocol["lifecycle"]["FRESH_READY_requires_TDG6_current_retry_count_zero"])
        self.assertTrue(self.protocol["lifecycle"]["RETRY_PENDING_requires_TDG6_current_retry_count_positive"])
        self.assertEqual(
            self.protocol["lifecycle"]["member_keys_in_canonical_order"],
            ["RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"],
        )
        self.assertTrue(
            self.protocol["lifecycle"][
                "checkpoint_must_bind_complete_six_member_cursor_ledger_state_sets_and_journal_tip"
            ]
        )
        self.assertEqual(
            self.protocol["lifecycle"]["terminal_exhaustion_classification_exact"],
            "invalid_implementation_or_nonconverged_run",
        )
        self.assertTrue(
            self.protocol["lifecycle"][
                "terminal_lock_forbids_fresh_and_retry_entrypoints_until_separately_authorized_new_campaign_id"
            ]
        )
        self.assertEqual(
            self.protocol["lifecycle"]["TDG6_typed_outcome_delivery"],
            "typed_exception_caught_by_campaign_owner_before_escape",
        )
        self.assertEqual(len(self.protocol["lifecycle"]["recovery_authority_order"]), 5)

    def test_adversarial_contract_mutations_fail_closed(self) -> None:
        attacks = {
            "threshold_inheritance": lambda value: value["inheritance"].__setitem__("PROTO14_TDG6_channel_order_threshold_and_minimum_step_inherited_unchanged", False),
            "lattice_history": lambda value: value["amendment"].__setitem__("stage_safe_repair_result_sha256", "0" * 64),
            "invalid_mode": lambda value: value["lifecycle"].__setitem__("cursor_modes", ["FRESH_READY"]),
            "successor_prefix": lambda value: value["lifecycle"]["retry_successor_required_fields"].remove("complete_rejection_prefix"),
            "retry_reset": lambda value: value["lifecycle"].__setitem__("RETRY_PENDING_may_not_reset_to_FRESH_READY_without_accepted_fine_commit_or_explicit_common_event_rollback", False),
            "checkpoint": lambda value: value["lifecycle"].__setitem__("checkpoint_must_bind_complete_six_member_cursor_ledger_state_sets_and_journal_tip", False),
            "unmatched_rejection": lambda value: value["lifecycle"].__setitem__("unmatched_durable_rejection_recovery_requires_exact_bound_evidence_and_typed_branch_reconstruction_or_invalid_stop", False),
            "mismatch_stop": lambda value: value["lifecycle"].__setitem__("mismatched_durable_rejection_or_successor_is_invalid_not_fresh", False),
            "fresh_retry_count": lambda value: value["lifecycle"].__setitem__("FRESH_READY_requires_TDG6_current_retry_count_zero", False),
            "retry_order": lambda value: value["lifecycle"]["nonexhausted_retry_transition_order"].reverse(),
            "exhaustion_order": lambda value: value["lifecycle"]["exhaustion_transition_order"].reverse(),
            "exhaustion_retries": lambda value: value["lifecycle"].__setitem__("TDG6_retry_exhaustion_must_not_create_RETRY_PENDING", False),
            "journal_tip": lambda value: value["lifecycle"]["cursor_required_fields"].remove("journal_tip_sha256"),
            "journal_envelope": lambda value: value["lifecycle"]["journal_record_envelope_required_fields"].remove("journal_parent_sha256"),
            "journal_tip_semantics": lambda value: value["lifecycle"].__setitem__("cursor_transition_record_sha256_becomes_checkpoint_journal_tip", False),
            "ledger_hash_rule": lambda value: value["lifecycle"].__setitem__("TDG6_ledger_sha256_rule", "unspecified"),
            "rollback_generation": lambda value: value["lifecycle"].__setitem__("rollback_must_append_new_generation_binding_restored_snapshot_not_reuse_old_cursor_bytes", False),
            "checkpoint_directory_fsync": lambda value: value["lifecycle"].__setitem__("atomic_checkpoint_requires_file_fsync_atomic_replace_and_parent_directory_fsync", False),
            "accepted_commit_generation": lambda value: value["lifecycle"].__setitem__("accepted_fine_state_TDG6_ledger_and_FRESH_READY_cursor_must_publish_as_one_atomic_checkpoint_generation", False),
            "pre_record_crash_restore": lambda value: value["lifecycle"].__setitem__("crash_before_new_durable_record_restores_exact_last_durable_cursor", False),
            "member_set": lambda value: value["lifecycle"]["member_keys_in_canonical_order"].pop(),
            "retry_relation": lambda value: value["lifecycle"]["retry_successor_semantic_relations"].remove("successor_cap_is_exact_binary64_half_of_predecessor_cap"),
            "terminal_lock": lambda value: value["lifecycle"].__setitem__("terminal_exhaustion_requires_terminal_lock_true", False),
            "common_event_set": lambda value: value["lifecycle"].__setitem__("common_event_rollback_is_one_campaign_record_over_complete_six_member_snapshot", False),
            "recovery_order": lambda value: value["lifecycle"]["recovery_authority_order"].reverse(),
            "attempt_semantics": lambda value: value["lifecycle"].__setitem__("attempt_serial_strictly_increases_for_every_fresh_or_retry_proposal", False),
            "typed_outcome_owner": lambda value: value["lifecycle"].__setitem__("campaign_owner_must_complete_retry_or_terminal_persistence_before_typed_outcome_escapes", False),
            "restart_payload": lambda value: value["restart_inputs"]["members"][0].__setitem__("state_sha256", "0" * 64),
            "claim": lambda value: value["claims"].__setitem__("FGCQR_holdout_execution_authorized", True),
        }
        for name, mutate in attacks.items():
            with self.subTest(name=name):
                attacked = deepcopy(self.protocol)
                mutate(attacked)
                with self.assertRaises(ValueError):
                    validate_sf1_protocol_v15(attacked)


if __name__ == "__main__":
    unittest.main()
