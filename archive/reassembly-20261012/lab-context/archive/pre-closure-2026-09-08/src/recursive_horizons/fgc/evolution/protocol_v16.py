"""Pure validation for PROTO16's prospective common-event commit contract."""

from __future__ import annotations

from typing import Any, Mapping

SF1_PROTOCOL_V16_ARTIFACT_ID = "FGC-2-SF1-PROTO16"
MEMBERS = [
    "RK4-2049",
    "RK4-4097",
    "RK4-8193",
    "SSPRK3-4097",
    "SSPRK3-8193",
    "SSPRK3-16385",
]
RECEIPT_PAYLOAD_FIELDS = [
    "protocol_artifact_id",
    "campaign_id",
    "previous_campaign_generation",
    "previous_checkpoint_sha256",
    "previous_committed_common_event_index",
    "completed_common_event_index",
    "completed_event_time",
    "next_event_target_time",
    "member_keys",
    "prior_cursor_set_sha256",
    "prior_ledger_set_sha256",
    "prior_accepted_state_set_sha256",
]
GENESIS_SPEC_FIELDS = [
    "protocol_artifact_id",
    "protocol_config_sha256",
    "protocol_freeze_result_sha256",
    "hlt13_result_sha256",
    "run_plan_sha256",
    "runtime_module_sha256",
    "adapter_module_sha256",
    "runner_sha256",
    "numerical_environment",
    "campaign_id",
    "namespace",
    "journal_root_parent_sha256",
    "common_event_index",
    "active_event_target_time",
    "member_descriptors",
    "member_descriptor_set_sha256",
    "physical_bundle_location",
    "physical_bundle_sha256",
    "genesis_checkpoint_sha256",
    "cursor_set_sha256",
    "ledger_set_sha256",
    "accepted_state_set_sha256",
]
MEMBER_DESCRIPTOR_FIELDS = [
    "member_key",
    "method",
    "point_count",
    "source_checkpoint",
    "state_sha256",
    "restart_payload_sha256",
    "state_shape",
    "state_dtype",
    "state_layout",
    "accepted_boundary_time",
    "step_index",
    "transaction_serial",
    "initial_TDG6_ledger_sha256",
    "physical_bundle_member_sha256",
]
TRUE_CLAIMS = {
    "PROTO16_frozen",
    "PROTO16_common_event_commit_contract_frozen",
    "PROTO16_external_trusted_genesis_contract_frozen",
    "PROTO15_runtime_implementation_reused_but_PROTO16_transition_unimplemented",
}
FALSE_CLAIMS = {
    "PROTO16_successor_runtime_implemented",
    "PROTO16_pretrajectory_runtime_authorized",
    "PROTO16_fresh_GR0_dynamic_calibration_authorized",
    "PROTO16_trajectory_read",
    "PROTO16_output_namespace_created",
    "PROTO16_calibration_completed",
    "GR0_case_eligible",
    "classical_spherical_diagnostic_authorized",
    "SGBL_execution_authorized",
    "FGCQR_holdout_execution_authorized",
    "DEF1_execution_authorized",
    "retained_EFT_evolution_authorized",
    "physical_transition_claim_authorized",
    "FGCQR_mechanism_rejected",
    "general_gradient_route_rejected",
    "singularity_resolution_derived",
    "child_domain_or_topology_derived",
    "dark_sector_mechanism_derived",
    "varying_locally_measured_c_derived",
}


def validate_sf1_protocol_v16(protocol: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(protocol, Mapping):
        raise TypeError("PROTO16 must be a mapping")
    value = dict(protocol)
    expected_root = {
        "schema_version",
        "artifact_id",
        "protocol_version",
        "project_version",
        "frozen",
        "units",
        "purpose",
        "amendment",
        "restart_inputs",
        "event_commit",
        "trusted_genesis",
        "namespace",
        "claims",
    }
    if set(value) != expected_root:
        raise ValueError("PROTO16 root keys differ")
    if (
        value["schema_version"],
        value["artifact_id"],
        value["protocol_version"],
        value["project_version"],
        value["frozen"],
        value["units"],
    ) != (
        1,
        SF1_PROTOCOL_V16_ARTIFACT_ID,
        16,
        "0.11.0",
        True,
        "c=hbar=1; code length L0=4",
    ):
        raise ValueError("PROTO16 identity differs")
    if value["purpose"] != (
        "outcome-neutral_constraint-controlled_spherical_mechanism_falsification"
    ):
        raise ValueError("PROTO16 purpose differs")
    amendment = value["amendment"]
    expected_amendment = {
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO15",
        "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v15.toml",
        "predecessor_protocol_sha256": (
            "f1f1849d85b6cbe1eb61267392a63ca11389f551df127033ddf905abb38fdc3b"
        ),
        "predecessor_freeze_artifact_id": "FGC-1-PRO15-FRZ1",
        "predecessor_freeze_config": "configs/fgc/fgc-1-pro15-frz1.toml",
        "predecessor_freeze_config_sha256": (
            "1a45af49cca92701a79cf103fd6795ffebf6353af372bff779cd4aae90f71676"
        ),
        "predecessor_freeze_result": "results/fgc-1-pro15-frz1.json",
        "predecessor_freeze_result_sha256": (
            "99d24425cb5d40da2e248605265b83bcdc82a2126cfa4b5bce4239c6ef30635c"
        ),
        "hlt13_artifact_id": "FGC-1-HLT13-MON13",
        "hlt13_result": "results/fgc-1-hlt13-mon13.json",
        "hlt13_result_sha256": (
            "de55f8d55a1db0043945d3553c7d2070703809b82c1f7141bf923c036b9df5b7"
        ),
        "evidence_checkpoint_commit": ("fddc7455621e805f833bc7a89bcf6712f27c6204"),
        "revision_classification": (
            "premise_only_common_event_commit_and_trusted_genesis_contract_freeze"
        ),
        "all_physical_numerical_TDG6_TDG7_and_non_event_PROTO15_rules_inherited_unchanged": True,
        "historical_raw_campaign_bundle_must_not_be_read_by_this_freeze": True,
        "SGBL_evolution_outcomes_inspected_before_revision": False,
        "FGCQR_evolution_outcomes_inspected_before_revision": False,
        "mechanism_question_answered_before_revision": False,
    }
    if amendment != expected_amendment:
        raise ValueError("PROTO16 amendment differs")
    if value["restart_inputs"] != {
        "restart_coordinate_time": "23/16",
        "branch": "GR-0",
        "amplitude": "3",
        "member_keys_in_canonical_order": MEMBERS,
        "six_restart_descriptors_inherited_by_exact_PROTO15_hash": True,
        "no_state_reinitialization_reprojection_interpolation_or_refitting_is_permitted": True,
    }:
        raise ValueError("PROTO16 restart descriptors differ")
    event = value["event_commit"]
    required_event = {
        "journal_record_kind",
        "required_member_mode",
        "all_six_members_required",
        "accepted_boundary_time_must_equal_active_event_target",
        "current_temporal_retry_count_must_equal_zero",
        "completed_common_event_index_is_previous_plus_one",
        "completed_event_time_is_prior_active_event_target",
        "receipt_payload_required_fields",
        "self_successor_cursor_set_and_checkpoint_hashes_forbidden_in_receipt_payload",
        "receipt_sha256_is_outer_journal_record_sha256",
        "receipt_journal_parent_is_prior_checkpoint_journal_tip",
        "durable_canonical_journal_record_before_checkpoint",
        "checkpoint_must_follow_fsynced_commit_record",
        "new_checkpoint_binds_derived_successor_cursor_set_after_receipt_hash_exists",
        "new_checkpoint_keeps_state_and_ledger_set_hashes_unchanged",
        "new_checkpoint_increments_committed_event_index",
        "frozen_common_event_cadence",
        "new_target_is_exact_prior_target_plus_cadence",
        "all_cursor_generations_and_attempt_ids_are_new",
        "new_cursor_generation_and_attempt_serial_are_prior_plus_one",
        "new_attempt_id_uses_completed_common_event_index",
        "new_cursor_parent_is_prior_cursor_hash",
        "new_cursor_journal_tip_is_commit_record_hash",
        "only_lone_durable_common_event_commit_suffix_is_exactly_reconstructible",
        "any_other_uncheckpointed_suffix_is_invalid",
        "future_runtime_fault_injects_common_event_commit_record_and_checkpoint_windows",
    }
    excluded = {
        "journal_record_kind",
        "required_member_mode",
        "receipt_payload_required_fields",
        "frozen_common_event_cadence",
    }
    if (
        set(event) != required_event
        or event["journal_record_kind"] != "COMMON_EVENT_COMMIT"
        or event["required_member_mode"] != "FRESH_READY"
        or event["receipt_payload_required_fields"] != RECEIPT_PAYLOAD_FIELDS
        or event["frozen_common_event_cadence"] != "1/16"
        or any(event[key] is not True for key in required_event - excluded)
    ):
        raise ValueError("PROTO16 common-event contract differs")
    genesis = value["trusted_genesis"]
    expected_genesis = {
        "authority_artifact_id",
        "authority_must_be_committed_before_namespace_creation",
        "authority_binds_external_immutable_manifest_or_commit",
        "genesis_spec_location",
        "genesis_spec_required_fields",
        "member_descriptor_required_fields",
        "member_descriptors_are_complete_canonical_six_member_collection",
        "member_descriptor_set_sha256_is_checksum_not_substitute",
        "genesis_checkpoint_is_derivable_without_namespace_generated_value",
        "genesis_spec_and_checkpoint_must_not_embed_genesis_spec_sha256",
        "launch_authority_tuple_required_fields",
        "authorization_result_must_not_self_reference_its_own_commit_or_result_hash",
        "runner_must_prove_authorization_result_blob_is_tracked_at_launch_authorization_commit",
        "wholesale_self_consistent_store_replacement_is_not_authenticated_without_authority",
        "trusted_Git_history_is_external_root_assumption",
    }
    if (
        set(genesis) != expected_genesis
        or genesis["authority_artifact_id"] != "FGC-1-HLT14-MON14"
        or genesis["genesis_spec_location"] != "embedded_in_future_HLT14_result"
        or genesis["genesis_spec_required_fields"] != GENESIS_SPEC_FIELDS
        or genesis["member_descriptor_required_fields"] != MEMBER_DESCRIPTOR_FIELDS
        or genesis["launch_authority_tuple_required_fields"]
        != [
            "authorization_commit",
            "authorization_result_path",
            "authorization_result_sha256",
            "genesis_spec_sha256",
        ]
    ):
        raise ValueError("PROTO16 trusted-genesis schema differs")
    for key in expected_genesis - {
        "authority_artifact_id",
        "genesis_spec_location",
        "genesis_spec_required_fields",
        "member_descriptor_required_fields",
        "launch_authority_tuple_required_fields",
    }:
        if genesis[key] is not True:
            raise ValueError("PROTO16 trusted-genesis invariant differs")
    if value["namespace"] != {
        "calibration_output_root": "runs/fgc-2-sf1/proto16/calibration",
        "holdout_output_root": "runs/fgc-2-sf1/proto16/holdout",
        "both_roots_must_be_absent": True,
        "freeze_must_create_neither_root": True,
    }:
        raise ValueError("PROTO16 namespace differs")
    claims = value["claims"]
    if (
        set(claims) != TRUE_CLAIMS | FALSE_CLAIMS
        or any(claims[key] is not True for key in TRUE_CLAIMS)
        or any(claims[key] is not False for key in FALSE_CLAIMS)
    ):
        raise ValueError("PROTO16 claims differ")
    return value


__all__ = [
    "FALSE_CLAIMS",
    "GENESIS_SPEC_FIELDS",
    "MEMBER_DESCRIPTOR_FIELDS",
    "MEMBERS",
    "RECEIPT_PAYLOAD_FIELDS",
    "SF1_PROTOCOL_V16_ARTIFACT_ID",
    "TRUE_CLAIMS",
    "validate_sf1_protocol_v16",
]
