"""Pure validation for the pre-trajectory PROTO15 lifecycle replacement.

PROTO15 freezes a new *campaign-owned* temporal cursor contract after the
historical PROTO14 runtime was classified invalid.  It deliberately implements
no cursor, journal, checkpoint, runner, namespace, or trajectory: those are
separate future artifacts.  This module therefore validates only the immutable
protocol data and cannot authorize a numerical or physical result.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V15_ARTIFACT_ID = "FGC-2-SF1-PROTO15"

PROTO15_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO14",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v14.toml",
    "predecessor_protocol_sha256": "3043cbc23e7a78af764b6e74397512cb3f7bc78eace64485ab188963ba35f3a3",
    "predecessor_freeze_artifact_id": "FGC-1-PRO14-FRZ1",
    "predecessor_freeze_result": "results/fgc-1-pro14-frz1.json",
    "predecessor_freeze_result_sha256": "a0bac04e78fe9ed5af851694057e51a0de62e8682a2fbfb734e63c98537c305b",
    "historical_terminal_artifact_id": "FGC-1-CAL11-PREF22",
    "historical_terminal_result": "results/fgc-1-cal11-pref22.json",
    "historical_terminal_result_sha256": "1339792198ec78ad83bc3e11a2805f09f961607ab6c1f6f06e79bb88ba6b58b8",
    "historical_runtime_authorization_artifact_id": "FGC-1-HLT12-MON12",
    "historical_runtime_authorization_result": "results/fgc-1-hlt12-mon12.json",
    "historical_runtime_authorization_result_sha256": "7f4068a3116ff1838ce7f0165ae18e72d2eea68933b87f056b26d8149c935731",
    "stage_safe_repair_artifact_id": "FGC-1-TDG7-IMP3",
    "stage_safe_repair_config": "configs/fgc/fgc-1-tdg7-imp3.toml",
    "stage_safe_repair_config_sha256": "1c989c6e698e40a3237931508b558fe1c8528316e35e8c29c8fc0e2e320c85e8",
    "stage_safe_repair_result": "results/fgc-1-tdg7-imp3.json",
    "stage_safe_repair_result_sha256": "364223f0f0f93c446b9d5b6cf09d3aa3cbe5b48bd1c175faf96a34459fd5b779",
    "evidence_checkpoint_commit": "64316d8f80f8532524c306817239f823a73b14be",
    "revision_classification": "premise_only_campaign_owned_temporal_lifecycle_replacement_freeze",
    "revision_reasons": [
        "CAL11_preserves_PROTO14_as_invalid_runtime_or_nonconverged_not_physics",
        "TDG7_IMP3_qualifies_the_stage_safe_lattice_repair_only_synthetically",
        "TDG7_IMP3_exposes_that_a_stateless_adapter_cannot_detect_a_discarded_external_retry",
        "a_fresh_campaign_must_record_retry_intent_durably_before_retry_execution",
    ],
    "revision_changes": [
        "replace_stateless_entrypoint_inference_with_campaign_owned_FRESH_READY_and_RETRY_PENDING_cursor",
        "require_exact_TDG6_retry_successor_and_complete_rejection_prefix_in_every_retry_transition",
        "require_fsynced_retry_transition_before_retry_execution",
        "require_cursor_chain_and_checkpoint_binding_with_explicit_common_event_rollback",
        "require_retry_exhaustion_to_commit_terminal_invalid_record_without_retry_cursor",
        "freeze_canonical_cursor_ledger_journal_and_rollback_hashing_rules",
    ],
    "all_physical_numerical_and_non_lifecycle_PROTO14_contract_fields_inherited_by_exact_predecessor_hash": True,
    "PROTO14_invalid_runtime_history_preserved": True,
    "PROTO14_terminal_checkpoint_may_not_resume": True,
    "historical_raw_campaign_bundle_must_not_be_read_by_this_freeze": True,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO15_RESTART_INPUTS = {
    "restart_coordinate_time": "23/16",
    "branch": "GR-0",
    "amplitude": "3",
    "members": [
        {"key": "RK4-2049", "source": "PROTO12", "method": "RK4", "point_count": 2049, "state_sha256": "d119a3b0ab450b0f0916b801254bbbfa964c8e5e95223cb8e51d9034f62c665e", "restart_payload_sha256": "4bb04ae97c09f66881084b5816d14ce300c54152ae11c182f53b64d5e887c072", "step_index": 342, "transaction_serial": 1710},
        {"key": "RK4-4097", "source": "PROTO12", "method": "RK4", "point_count": 4097, "state_sha256": "9d27ca221fbac2a5b3b86a8d9aef49a8b8ac451b255b80529563b8cb7f8b73a9", "restart_payload_sha256": "cc6d3b6fbb3a5b6abd93f69f7f59ba67c38e9b7d332fac991a88847970b57aac", "step_index": 673, "transaction_serial": 3365},
        {"key": "RK4-8193", "source": "PROTO12", "method": "RK4", "point_count": 8193, "state_sha256": "d5b4cfe25199862fff5be7d524cc1c3a6030cc634d268cedb626c3ced7220f4a", "restart_payload_sha256": "b95ea0be0779f712fb0c354f9c10f28b9ccc6c737fb403d5710dbd07e0e506f9", "step_index": 976, "transaction_serial": 4880},
        {"key": "SSPRK3-4097", "source": "PROTO12", "method": "SSPRK3", "point_count": 4097, "state_sha256": "ab96602c0f05633375710efee2be01be84e24c32a31d07d6ef23390788a72162", "restart_payload_sha256": "dc056db806e5d983852135138dd48dbb9b995a8585f56900258bece7e16e7016", "step_index": 631, "transaction_serial": 2524},
        {"key": "SSPRK3-8193", "source": "PROTO12", "method": "SSPRK3", "point_count": 8193, "state_sha256": "c46abc3f0f72fa2a13c1fee7bb329b94698cc2c26003406ae3eb1592f11b648c", "restart_payload_sha256": "b262bf4bd180717d42a7fa4d502f006bac0950669bc9e86c474ad6e4a21ce7ad", "step_index": 987, "transaction_serial": 3948},
        {"key": "SSPRK3-16385", "source": "RSP2", "method": "SSPRK3", "point_count": 16385, "state_sha256": "91e64d0864390cb4ac63cfdbdbbcb8ca5ad6c2f47cdd4d6741e05ef560adad0d", "restart_payload_sha256": "62784c8cc63eedc210366f146f28d535343b966e1f988a4298894e78ec560884", "step_index": 1771, "transaction_serial": 7084},
    ],
    "no_state_reinitialization_reprojection_interpolation_or_refitting_is_permitted": True,
    "PROTO14_post_restart_or_terminal_state_may_not_be_used": True,
}

PROTO15_LIFECYCLE = {
    "owner": "campaign_runtime_not_TDG7_entrypoint_inference",
    "canonical_record_encoding": "utf8_sorted_compact_json_ensure_ascii_true_allow_nan_false",
    "binary64_identity_encoding": "lowercase_float_hex",
    "member_keys_in_canonical_order": ["RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"],
    "cursor_modes": ["FRESH_READY", "RETRY_PENDING"],
    "cursor_required_fields": ["protocol_artifact_id", "campaign_id", "member_key", "method", "point_count", "committed_common_event_index", "accepted_boundary_time", "accepted_state_sha256", "previous_step_index", "previous_transaction_serial", "TDG6_ledger_sha256", "cursor_generation", "event_target_time", "attempt_serial", "attempt_id", "mode", "retry_successor_payload_or_none", "journal_tip_sha256", "cursor_chain_parent_sha256", "cursor_chain_sha256"],
    "time_identity_fields_require_rational_and_binary64_hex": True,
    "attempt_id_rule": "canonical_campaign_id_member_key_common_event_index_attempt_serial",
    "attempt_serial_strictly_increases_for_every_fresh_or_retry_proposal": True,
    "attempt_ids_may_not_be_reused_after_retry_rollback_or_restart": True,
    "TDG6_ledger_hash_required_fields": ["last_accepted_time_hex", "accepted_macro_step_count", "cumulative_temporal_retry_count", "current_macro_step_temporal_retry_count", "last_accepted_macro_step_temporal_retry_count", "accumulated_debit_vector_hex", "serialized_temporal_rejections"],
    "TDG6_ledger_sha256_rule": "sha256_of_canonical_complete_TDG6_ledger_mapping",
    "journal_record_envelope_required_fields": ["record_kind", "journal_parent_sha256", "payload", "record_sha256"],
    "journal_record_kinds": ["TDG6_REJECTION", "CURSOR_TRANSITION", "ACCEPTED_FINE_COMMIT", "TERMINAL_INVALID", "COMMON_EVENT_ROLLBACK"],
    "journal_chain_root_parent_sha256": "0000000000000000000000000000000000000000000000000000000000000000",
    "journal_record_sha256_rule": "sha256_of_canonical_record_without_record_sha256",
    "cursor_journal_tip_semantics": "hash_of_latest_durable_journal_record_preceding_cursor_transition",
    "cursor_transition_record_sha256_becomes_checkpoint_journal_tip": True,
    "cursor_chain_root_parent_sha256": "0000000000000000000000000000000000000000000000000000000000000000",
    "cursor_chain_sha256_rule": "sha256_of_canonical_cursor_without_cursor_chain_sha256",
    "cursor_chain_parent_must_equal_immediately_preceding_same_member_cursor_chain_sha256": True,
    "cursor_generation_strictly_increases_on_every_durable_transition_including_rollback": True,
    "allowed_nonterminal_transitions": ["FRESH_READY_to_RETRY_PENDING_on_durable_nonexhausted_rejection", "RETRY_PENDING_to_RETRY_PENDING_on_durable_nonexhausted_rejection", "FRESH_READY_or_RETRY_PENDING_to_FRESH_READY_on_atomic_accepted_fine_commit", "any_nonterminal_mode_to_bound_snapshot_on_explicit_common_event_rollback"],
    "FRESH_READY_requires_retry_successor_payload_json_null": True,
    "FRESH_READY_requires_TDG6_current_retry_count_zero": True,
    "RETRY_PENDING_requires_exact_retry_successor_payload": True,
    "RETRY_PENDING_requires_TDG6_current_retry_count_positive": True,
    "retry_successor_required_fields": ["predecessor_plan", "successor_plan", "TDG6_rejection_evidence", "complete_rejection_prefix", "durable_rejection_record_sha256", "predecessor_journal_tip_sha256", "predecessor_attempt_id", "successor_attempt_id", "member_key", "method", "point_count", "accepted_state_sha256", "accepted_boundary_time", "previous_step_index", "previous_transaction_serial", "retry_count", "half_cap", "minimum_macro_step"],
    "retry_successor_semantic_relations": ["successor_attempt_serial_equals_predecessor_attempt_serial_plus_one", "retry_count_equals_predecessor_current_retry_count_plus_one", "successor_cap_is_exact_binary64_half_of_predecessor_cap", "accepted_member_method_grid_state_time_step_and_serial_unchanged", "TDG6_rejection_evidence_equals_complete_rejection_prefix_tail", "durable_rejection_record_payload_equals_TDG6_rejection_evidence", "RETRY_PENDING_cursor_journal_tip_equals_durable_rejection_record_sha256"],
    "retry_transition_record_required_fields": ["predecessor_cursor_chain_sha256", "durable_rejection_record_sha256", "successor_cursor_chain_sha256", "successor_attempt_id"],
    "TDG6_post_rejection_typed_outcomes": ["TDG6TemporalRetryRequired", "TDG6TemporalRetryExhausted"],
    "TDG6_typed_outcome_delivery": "typed_exception_caught_by_campaign_owner_before_escape",
    "campaign_owner_must_complete_retry_or_terminal_persistence_before_typed_outcome_escapes": True,
    "nonexhausted_retry_transition_order": ["durable_TDG6_rejection_journal_record", "typed_TDG6_retry_required_raised", "derive_exact_successor", "append_and_fsync_RETRY_PENDING_cursor_transition", "adopt_RETRY_PENDING_cursor", "invoke_retry_entrypoint"],
    "exhaustion_transition_order": ["durable_TDG6_rejection_journal_record", "typed_TDG6_retry_exhaustion_raised", "append_and_fsync_terminal_invalid_nonconverged_record", "checkpoint_terminal_state"],
    "terminal_exhaustion_required_fields": ["protocol_artifact_id", "campaign_id", "committed_common_event_index", "campaign_generation", "attempt_serial", "attempt_id", "member_key", "method", "point_count", "accepted_boundary_time", "accepted_state_sha256", "previous_step_index", "previous_transaction_serial", "TDG6_ledger_sha256", "durable_rejection_record_sha256", "predecessor_cursor_chain_sha256", "journal_tip_sha256", "exhaustion_reason", "classification", "terminal_lock"],
    "terminal_exhaustion_classification_exact": "invalid_implementation_or_nonconverged_run",
    "terminal_exhaustion_requires_terminal_lock_true": True,
    "terminal_checkpoint_required_fields": ["protocol_artifact_id", "campaign_id", "campaign_generation", "committed_common_event_index", "active_event_target_time", "disposition", "terminal_member_key", "terminal_record_sha256", "cursor_set_sha256", "ledger_set_sha256", "accepted_state_set_sha256", "journal_tip_sha256", "terminal_lock", "checkpoint_sha256"],
    "terminal_lock_forbids_fresh_and_retry_entrypoints_until_separately_authorized_new_campaign_id": True,
    "retry_must_be_durable_before_execution": True,
    "TDG6_retry_exhaustion_must_not_create_RETRY_PENDING": True,
    "TDG6_retry_exhaustion_is_invalid_nonconverged_not_physical": True,
    "fresh_entrypoint_forbids_RETRY_PENDING": True,
    "retry_entrypoint_requires_RETRY_PENDING": True,
    "accepted_fine_commit_must_advance_state_time_and_reset_current_TDG6_retry_count_before_FRESH_READY": True,
    "accepted_fine_state_TDG6_ledger_and_FRESH_READY_cursor_must_publish_as_one_atomic_checkpoint_generation": True,
    "RETRY_PENDING_may_not_reset_to_FRESH_READY_without_accepted_fine_commit_or_explicit_common_event_rollback": True,
    "cursor_set_member_record_fields": ["member_key", "cursor_chain_sha256", "TDG6_ledger_sha256", "accepted_state_sha256", "accepted_boundary_time", "cursor_generation", "committed_common_event_index"],
    "cursor_set_sha256_rule": "sha256_of_canonical_list_in_member_keys_in_canonical_order",
    "ledger_set_sha256_rule": "sha256_of_canonical_member_key_to_complete_TDG6_ledger_sha256_mapping",
    "accepted_state_set_sha256_rule": "sha256_of_canonical_member_key_to_accepted_state_sha256_mapping",
    "campaign_checkpoint_required_fields": ["protocol_artifact_id", "campaign_id", "campaign_generation", "committed_common_event_index", "active_event_target_time", "disposition", "member_keys", "cursor_set_sha256", "ledger_set_sha256", "accepted_state_set_sha256", "journal_tip_sha256", "terminal_lock", "checkpoint_sha256"],
    "campaign_generation_strictly_increases_on_every_atomic_checkpoint_publication": True,
    "checkpoint_must_bind_complete_six_member_cursor_ledger_state_sets_and_journal_tip": True,
    "mixed_or_incomplete_member_event_cursor_set_is_invalid": True,
    "durable_journal_append_requires_one_canonical_record_flush_and_fsync_before_return": True,
    "atomic_checkpoint_requires_file_fsync_atomic_replace_and_parent_directory_fsync": True,
    "common_event_snapshot_required_fields": ["campaign_id", "campaign_generation", "committed_common_event_index", "event_time", "member_keys", "cursor_set_sha256", "ledger_set_sha256", "accepted_state_set_sha256", "journal_tip_sha256", "checkpoint_sha256"],
    "common_event_abort_must_restore_complete_six_member_snapshot_and_append_explicit_rollback": True,
    "common_event_rollback_is_one_campaign_record_over_complete_six_member_snapshot": True,
    "rollback_must_append_new_generation_binding_restored_snapshot_not_reuse_old_cursor_bytes": True,
    "rollback_record_required_fields": ["aborted_attempt_id", "aborted_campaign_generation", "aborted_cursor_set_sha256", "aborted_ledger_set_sha256", "aborted_accepted_state_set_sha256", "restored_snapshot_campaign_generation", "restored_snapshot_checkpoint_sha256", "restored_snapshot_cursor_set_sha256", "restored_snapshot_ledger_set_sha256", "restored_snapshot_accepted_state_set_sha256", "new_campaign_generation", "new_cursor_set_sha256", "new_ledger_set_sha256", "new_accepted_state_set_sha256", "new_journal_parent_sha256", "new_journal_record_sha256"],
    "rollback_restored_snapshot_must_be_verified_ancestor_of_aborted_checkpoint": True,
    "rollback_must_increment_every_member_cursor_generation_and_campaign_generation": True,
    "rollback_must_issue_new_attempt_ids_and_may_not_reuse_aborted_attempt_id": True,
    "crash_before_new_durable_record_restores_exact_last_durable_cursor": True,
    "crash_after_durable_rejection_before_cursor_or_terminal_persistence_may_not_recompute_or_classify_fresh": True,
    "unmatched_durable_rejection_recovery_requires_exact_bound_evidence_and_typed_branch_reconstruction_or_invalid_stop": True,
    "mismatched_durable_rejection_or_successor_is_invalid_not_fresh": True,
    "crash_after_cursor_persistence_restores_exact_RETRY_PENDING": True,
    "recovery_authority_order": ["validate_highest_complete_atomic_checkpoint_generation", "validate_complete_six_member_cursor_ledger_state_sets", "scan_only_contiguous_canonical_journal_suffix_rooted_at_checkpoint_tip", "reconstruct_only_exact_bound_retry_required_or_retry_exhausted_branch", "materialize_one_new_atomic_checkpoint_generation_or_stop_invalid"],
    "checkpoint_controls_state_activation_except_exact_reconstructible_post_rejection_branch": True,
    "syntactically_present_uncheckpointed_record_is_not_otherwise_semantically_active": True,
    "torn_noncontiguous_forged_or_mixed_generation_journal_checkpoint_state_is_invalid": True,
    "old_or_new_checkpoint_after_replace_crash_must_be_selected_only_by_complete_hash_chain_and_generation_validation": True,
    "future_runtime_must_fault_inject_every_journal_cursor_terminal_checkpoint_and_six_member_rollback_boundary": True,
    "runtime_implementation_is_not_part_of_this_protocol_freeze": True,
}

PROTO15_INHERITANCE = {
    "PROTO14_all_physical_equations_source_backend_reference_map_CFL_dissipation_constraints_spatial_health_boundary_scale_null_and_terminal_rules_inherited_unchanged": True,
    "PROTO14_method_owned_ladders_and_all_numerical_thresholds_inherited_unchanged": True,
    "PROTO14_TDG6_channel_order_threshold_and_minimum_step_inherited_unchanged": True,
    "PROTO14_invalid_terminal_result_remains_historical_and_unreclassified": True,
    "TDG7_IMP3_stage_safe_lattice_repair_is_qualified_but_not_a_campaign_cursor": True,
    "successor_runtime_owner": "FGC-1-HLT13-MON13",
}

PROTO15_CLAIMS = {
    "PROTO15_frozen": True,
    "PROTO15_campaign_owned_cursor_contract_frozen": True,
    "PROTO15_six_restart_payloads_hash_bound": True,
    "PROTO14_invalid_runtime_history_preserved": True,
    "PROTO15_successor_runtime_implemented": False,
    "PROTO15_pretrajectory_runtime_authorized": False,
    "PROTO15_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO15_trajectory_read": False,
    "PROTO15_output_namespace_created": False,
    "PROTO15_calibration_completed": False,
    "GR0_case_eligible": False,
    "classical_spherical_diagnostic_authorized": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "DEF1_execution_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}


def _mapping(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


def _exact(name: str, value: Any, expected: Mapping[str, Any]) -> Mapping[str, Any]:
    answer = _mapping(name, value)
    if dict(answer) != dict(expected):
        raise ValueError(f"{name} differs")
    return answer


def _semantic_hash(value: object) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8")).hexdigest()


def validate_sf1_protocol_v15(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete fail-closed, freeze-only PROTO15 contract."""

    protocol = _mapping("protocol", protocol)
    if set(protocol) != {"schema_version", "artifact_id", "protocol_version", "project_version", "frozen", "units", "purpose", "amendment", "restart_inputs", "lifecycle", "inheritance", "claims"}:
        raise ValueError("PROTO15 root keys differ")
    identity = {key: protocol[key] for key in ("schema_version", "artifact_id", "protocol_version", "project_version", "frozen", "units", "purpose")}
    if identity != {"schema_version": 1, "artifact_id": SF1_PROTOCOL_V15_ARTIFACT_ID, "protocol_version": 15, "project_version": "0.11.0", "frozen": True, "units": "c=hbar=1; code length L0=4", "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification"}:
        raise ValueError("PROTO15 identity, version, units, purpose, or freeze differs")
    amendment = _exact("PROTO15 amendment", protocol["amendment"], PROTO15_AMENDMENT)
    restart_inputs = _exact("PROTO15 restart inputs", protocol["restart_inputs"], PROTO15_RESTART_INPUTS)
    lifecycle = _exact("PROTO15 lifecycle", protocol["lifecycle"], PROTO15_LIFECYCLE)
    inheritance = _exact("PROTO15 inheritance", protocol["inheritance"], PROTO15_INHERITANCE)
    claims = _exact("PROTO15 claims", protocol["claims"], PROTO15_CLAIMS)
    if any(value is not False for key, value in claims.items() if key not in {"PROTO15_frozen", "PROTO15_campaign_owned_cursor_contract_frozen", "PROTO15_six_restart_payloads_hash_bound", "PROTO14_invalid_runtime_history_preserved"}):
        raise ValueError("PROTO15 authorization and physical claims must remain fail-closed")
    contract = {"amendment": dict(amendment), "restart_inputs": dict(restart_inputs), "lifecycle": dict(lifecycle), "inheritance": dict(inheritance)}
    return {
        "artifact_id": SF1_PROTOCOL_V15_ARTIFACT_ID,
        "protocol_version": 15,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "campaign_owned_cursor_contract_defined": True,
        "historical_PROTO14_result_reclassified": False,
        "production_trajectory_authorized": False,
        "contract_sha256": _semantic_hash(contract),
        "claims": dict(claims),
    }


__all__ = ["PROTO15_AMENDMENT", "PROTO15_CLAIMS", "PROTO15_INHERITANCE", "PROTO15_LIFECYCLE", "PROTO15_RESTART_INPUTS", "SF1_PROTOCOL_V15_ARTIFACT_ID", "validate_sf1_protocol_v15"]
