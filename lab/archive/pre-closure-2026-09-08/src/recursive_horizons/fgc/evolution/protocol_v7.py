"""Pure validation for the outcome-neutral PROTO7 transaction overlay.

PROTO7 consumes immutable PROTO6 and CAL3 evidence.  It replaces stage-name
retry ownership with one transaction boundary: an accepted-state source miss
is terminal; a source-only miss anywhere in a wholly unaccepted proposal may
rollback and retry; any non-source miss vetoes retry.  This module reads no run
output and performs no evolution.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V7_ARTIFACT_ID = "FGC-2-SF1-PROTO7"

PROTO7_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO6",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v6.toml",
    "predecessor_protocol_sha256": "e50093914f7b84a88bd593361f8b3eb43942cf9a8502f637fc38d8e61dd82448",
    "diagnosis_artifact_id": "FGC-1-CAL3-PREF5",
    "diagnosis_config": "configs/fgc/fgc-1-cal3-pref5.toml",
    "diagnosis_config_sha256": "cec13d7e50d2546d423e7ffea6690226d71fbe5174073e78ffd2ac37f9e44275",
    "diagnosis_result": "results/fgc-1-cal3-pref5.json",
    "diagnosis_result_sha256": "669ee5d9a17a8d3e0b213b06fcb3ff426e0908b331fe38429b5704bbbf49080f",
    "diagnosis_checkpoint_commit": "27c63e7aca1439e75e18b576c6da414e458280ff",
    "revision_classification": "premise_only_post_calibration_general_unaccepted_proposal_stop_ownership_revision",
    "revision_reasons": [
        "PROTO6_classifies_retryability_by_integrator_stage_label",
        "CAL3_last_accepted_state_passes_the_unchanged_raw_source_gate",
        "CAL3_terminal_failure_is_source_only_inside_a_wholly_unaccepted_proposal",
        "CAL3_one_further_halving_passes_the_unchanged_complete_transaction",
    ],
    "revision_changes": [
        "accepted_state_source_precheck_remains_terminal",
        "every_source_only_failure_in_a_wholly_unaccepted_proposal_becomes_retryable",
        "any_non_source_failure_vetoes_retry",
        "complete_terminal_and_retry_exhaustion_observability",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO6_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_initial_data_and_holdout_cases_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "raw_source_residual_and_kinetic_condition_thresholds_unchanged": True,
    "constraint_spectral_trapped_boundary_and_candidate_health_stops_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["5/2", "3"],
    "PROTO1_through_PROTO6_preserved": True,
    "PROTO6_GR0_calibration_trajectory_inspected_before_revision": True,
    "first_nonzero_common_event_completed_before_revision": False,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO7_TRANSACTION_STOP_OWNERSHIP = {
    "accepted_state_source_gate_evaluated_before_each_proposal": True,
    "accepted_state_source_gate_failure_is_terminal": True,
    "accepted_state_source_gate_failure_reason": "newton_residual_limit",
    "proposal_is_unaccepted_until_complete_transaction_commit": True,
    "unaccepted_proposal_includes_all_integrator_stages_and_candidate_endpoint": True,
    "every_source_only_failure_in_unaccepted_proposal_is_retryable": True,
    "any_non_source_failure_in_unaccepted_proposal_vetoes_retry": True,
    "non_source_failure_priority_is_inherited_unchanged": True,
    "retry_does_not_accept_or_advance_the_failed_proposal": True,
    "retry_restarts_from_the_bitwise_last_accepted_state": True,
    "retry_factor": "1/2",
    "maximum_retry_count": 32,
    "minimum_step_size": "1/1073741824",
    "raw_source_residual_tolerance_is_unchanged": True,
    "source_refinement_algorithm_is_unchanged": True,
    "kinetic_condition_gate_is_unchanged_and_terminal": True,
    "first_scientific_failure_remains_immutable": True,
    "complete_proposal_must_be_previewed_before_retry_classification": True,
    "every_failed_source_call_stage_time_serial_and_residual_must_be_serialized": True,
    "retry_record_must_bind_member_step_retry_count_accepted_state_hash_and_rollback": True,
    "terminal_and_retry_exhaustion_records_must_retain_complete_failure_evidence": True,
}

PROTO7_PARTITION = {
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro7-hld1.toml",
}

PROTO7_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto7/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto7/holdout",
    "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    "PROTO6_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO7_outcomes": True,
}

PROTO7_INHERITANCE = {
    "PROTO6_native_semidiscrete_initialization_inherited_unchanged": True,
    "PROTO6_method_owned_constraint_admission_inherited_unchanged": True,
    "PROTO6_roundoff_zero_classification_inherited_unchanged": True,
    "PROTO6_stop_applicability_partition_inherited_except_transaction_source_ownership": True,
    "PROTO6_weighted_spectral_admission_inherited_unchanged": True,
    "PROTO6_calibration_selection_and_trapped_sign_rules_inherited_unchanged": True,
    "PROTO6_retry_factor_budget_and_minimum_step_inherited_unchanged": True,
    "HLT4_runtime_evidence_remains_valid_for_unchanged_rules_only": True,
    "CAL3_replay_is_design_evidence_not_calibration_evidence": True,
    "successor_runtime_owner": "FGC-1-HLT5-MON5",
}

PROTO7_CLAIMS = {
    "PROTO7_successor_runtime_compositor_implemented": False,
    "PROTO7_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO7_resolved_holdout_manifest_authorized": False,
    "classical_spherical_diagnostic_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
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
    mapping = _mapping(name, value)
    if dict(mapping) != dict(expected):
        raise ValueError(f"{name} differs")
    return mapping


def _semantic_hash(value: object) -> str:
    return sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def validate_sf1_protocol_v7(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete, fail-closed PROTO7 overlay without I/O."""

    protocol = _mapping("protocol", protocol)
    expected_root = {
        "schema_version",
        "artifact_id",
        "protocol_version",
        "project_version",
        "frozen",
        "units",
        "purpose",
        "amendment",
        "replacement",
        "inheritance",
        "claims",
    }
    if set(protocol) != expected_root:
        raise ValueError("PROTO7 root keys differ")
    identity = {
        key: protocol[key]
        for key in (
            "schema_version",
            "artifact_id",
            "protocol_version",
            "project_version",
            "frozen",
            "units",
            "purpose",
        )
    }
    if identity != {
        "schema_version": 1,
        "artifact_id": SF1_PROTOCOL_V7_ARTIFACT_ID,
        "protocol_version": 7,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO7 identity, version, units, purpose, or frozen state differs")

    amendment = _exact("PROTO7 amendment", protocol["amendment"], PROTO7_AMENDMENT)
    replacement = _mapping("PROTO7 replacement", protocol["replacement"])
    if set(replacement) != {"transaction_stop_ownership", "partition", "provenance"}:
        raise ValueError("PROTO7 replacement keys differ")
    ownership = _exact(
        "PROTO7 transaction-stop-ownership replacement",
        replacement["transaction_stop_ownership"],
        PROTO7_TRANSACTION_STOP_OWNERSHIP,
    )
    partition = _exact("PROTO7 partition replacement", replacement["partition"], PROTO7_PARTITION)
    provenance = _exact("PROTO7 provenance replacement", replacement["provenance"], PROTO7_PROVENANCE)
    inheritance = _exact("PROTO7 inheritance", protocol["inheritance"], PROTO7_INHERITANCE)
    claims = _exact("PROTO7 claims", protocol["claims"], PROTO7_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO7 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "transaction_stop_ownership": ownership,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V7_ARTIFACT_ID,
        "protocol_version": 7,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO6",
        "diagnosis_artifact_id": "FGC-1-CAL3-PREF5",
        "diagnosis_checkpoint_commit": amendment["diagnosis_checkpoint_commit"],
        "eligible_calibration_amplitudes": ["5/2", "3"],
        "accepted_state_source_failure_terminal": True,
        "all_source_only_unaccepted_proposal_failures_retryable": True,
        "candidate_endpoint_in_unaccepted_proposal": True,
        "non_source_failure_vetoes_retry": True,
        "complete_terminal_observability_required": True,
        "raw_source_residual_tolerance_unchanged": True,
        "all_non_source_stops_unchanged": True,
        "retry_factor": "1/2",
        "maximum_retry_count": 32,
        "minimum_step_size": "1/1073741824",
        "PROTO6_calibration_trajectory_inspected_before_revision": True,
        "first_nonzero_common_event_completed_before_revision": False,
        "FGCQR_outcomes_inspected_before_revision": False,
        "SGBL_outcomes_inspected_before_revision": False,
        "mechanism_question_answered_before_revision": False,
        "resolved_holdout_manifest": partition["resolved_holdout_manifest"],
        "calibration_output_root": provenance["calibration_output_root"],
        "holdout_output_root": provenance["holdout_output_root"],
        "resolved_holdout_manifest_present": False,
        "semantic_holdout_contract_sha256": _semantic_hash(overlay_contract),
        "claims": dict(claims),
    }


__all__ = [
    "PROTO7_AMENDMENT",
    "PROTO7_CLAIMS",
    "PROTO7_INHERITANCE",
    "PROTO7_PARTITION",
    "PROTO7_PROVENANCE",
    "PROTO7_TRANSACTION_STOP_OWNERSHIP",
    "SF1_PROTOCOL_V7_ARTIFACT_ID",
    "validate_sf1_protocol_v7",
]
