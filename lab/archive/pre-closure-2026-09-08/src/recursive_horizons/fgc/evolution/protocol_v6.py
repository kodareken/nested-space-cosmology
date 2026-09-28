"""Pure validation for the outcome-neutral PROTO6 stop-ownership overlay.

PROTO6 consumes PROTO5 and the post-calibration CAL2 diagnosis as immutable
evidence.  It changes only who owns a source residual miss: an accepted-state
miss remains a terminal premise failure, while a miss confined to an
unaccepted internal Runge--Kutta trial may halve and retry from the last
accepted state.  This module reads no run output and performs no evolution.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V6_ARTIFACT_ID = "FGC-2-SF1-PROTO6"

PROTO6_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO5",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v5.toml",
    "predecessor_protocol_sha256": "516061e19ac14857eda26076fa18d9163db2fcc7c9d31aa6aebefd698377c9ab",
    "diagnosis_artifact_id": "FGC-1-CAL2-PREF4",
    "diagnosis_config": "configs/fgc/fgc-1-cal2-pref4.toml",
    "diagnosis_config_sha256": "d596a45899baecbb5ea1f9d5d26409782376327ce6ddb2fc961e03626abc5a55",
    "diagnosis_result": "results/fgc-1-cal2-pref4.json",
    "diagnosis_result_sha256": "0e56f3b91205d68a6f81f4dad56f8ecce6d93ce0184a5899a1bace6d42aca11b",
    "diagnosis_checkpoint_commit": "47f4c83bae2df26c51c3544eed8c06507a0cd9eb",
    "revision_classification": "premise_only_post_calibration_trial_stage_stop_ownership_revision",
    "revision_reasons": [
        "PROTO5_treats_an_unaccepted_internal_stage_source_miss_as_terminal",
        "CAL2_last_accepted_state_passes_the_unchanged_raw_source_gate",
        "CAL2_first_inherited_half_step_retry_passes_the_unchanged_raw_source_gate",
    ],
    "revision_changes": [
        "accepted_state_source_precheck_remains_terminal",
        "internal_unaccepted_stage_source_miss_becomes_retryable",
        "inherited_half_step_retry_budget_and_minimum_step",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO5_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_initial_data_and_holdout_cases_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "raw_source_residual_and_kinetic_condition_thresholds_unchanged": True,
    "constraint_spectral_trapped_boundary_and_candidate_health_stops_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["5/2", "3"],
    "PROTO1_through_PROTO5_preserved": True,
    "PROTO5_GR0_calibration_trajectory_inspected_before_revision": True,
    "first_nonzero_common_event_completed_before_revision": False,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO6_SOURCE_STOP_OWNERSHIP = {
    "accepted_state_source_gate_evaluated_before_each_proposal": True,
    "accepted_state_source_gate_failure_is_terminal": True,
    "accepted_state_source_gate_failure_reason": "newton_residual_limit",
    "internal_unaccepted_stage_source_gate_failure_is_retryable": True,
    "retry_does_not_accept_or_advance_the_failed_trial": True,
    "retry_restarts_from_the_bitwise_last_accepted_state": True,
    "retry_factor": "1/2",
    "maximum_retry_count": 32,
    "minimum_step_size": "1/1073741824",
    "raw_source_residual_tolerance_is_unchanged": True,
    "source_refinement_algorithm_is_unchanged": True,
    "kinetic_condition_gate_is_unchanged_and_terminal": True,
    "non_source_stage_failures_retain_PROTO5_ownership": True,
    "first_scientific_failure_remains_immutable": True,
    "every_rejected_trial_residual_stage_and_step_size_must_be_serialized": True,
}

PROTO6_PARTITION = {
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro6-hld1.toml",
}

PROTO6_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto6/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto6/holdout",
    "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    "PROTO5_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO6_outcomes": True,
}

PROTO6_INHERITANCE = {
    "PROTO5_native_semidiscrete_initialization_inherited_unchanged": True,
    "PROTO5_method_owned_constraint_admission_inherited_unchanged": True,
    "PROTO5_roundoff_zero_classification_inherited_unchanged": True,
    "PROTO5_stop_applicability_partition_inherited_except_source_trial_ownership": True,
    "PROTO5_weighted_spectral_admission_inherited_unchanged": True,
    "PROTO5_calibration_selection_and_trapped_sign_rules_inherited_unchanged": True,
    "HLT3_runtime_evidence_remains_valid_for_unchanged_rules_only": True,
    "CAL2_replay_is_design_evidence_not_calibration_evidence": True,
    "successor_runtime_owner": "FGC-1-HLT4-MON4",
}

PROTO6_CLAIMS = {
    "PROTO6_successor_runtime_compositor_implemented": False,
    "PROTO6_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO6_resolved_holdout_manifest_authorized": False,
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


def validate_sf1_protocol_v6(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete, fail-closed PROTO6 overlay without I/O."""

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
        raise ValueError("PROTO6 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V6_ARTIFACT_ID,
        "protocol_version": 6,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO6 identity, version, units, purpose, or frozen state differs")

    amendment = _exact("PROTO6 amendment", protocol["amendment"], PROTO6_AMENDMENT)
    replacement = _mapping("PROTO6 replacement", protocol["replacement"])
    if set(replacement) != {"source_stop_ownership", "partition", "provenance"}:
        raise ValueError("PROTO6 replacement keys differ")
    source_ownership = _exact(
        "PROTO6 source-stop-ownership replacement",
        replacement["source_stop_ownership"],
        PROTO6_SOURCE_STOP_OWNERSHIP,
    )
    partition = _exact(
        "PROTO6 partition replacement",
        replacement["partition"],
        PROTO6_PARTITION,
    )
    provenance = _exact(
        "PROTO6 provenance replacement",
        replacement["provenance"],
        PROTO6_PROVENANCE,
    )
    inheritance = _exact(
        "PROTO6 inheritance", protocol["inheritance"], PROTO6_INHERITANCE
    )
    claims = _exact("PROTO6 claims", protocol["claims"], PROTO6_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO6 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "source_stop_ownership": source_ownership,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V6_ARTIFACT_ID,
        "protocol_version": 6,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO5",
        "diagnosis_artifact_id": "FGC-1-CAL2-PREF4",
        "diagnosis_checkpoint_commit": amendment["diagnosis_checkpoint_commit"],
        "eligible_calibration_amplitudes": ["5/2", "3"],
        "accepted_state_source_failure_terminal": True,
        "internal_unaccepted_stage_source_failure_retryable": True,
        "raw_source_residual_tolerance_unchanged": True,
        "all_non_source_stops_unchanged": True,
        "retry_factor": "1/2",
        "maximum_retry_count": 32,
        "minimum_step_size": "1/1073741824",
        "PROTO5_calibration_trajectory_inspected_before_revision": True,
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
    "PROTO6_AMENDMENT",
    "PROTO6_CLAIMS",
    "PROTO6_INHERITANCE",
    "PROTO6_PARTITION",
    "PROTO6_PROVENANCE",
    "PROTO6_SOURCE_STOP_OWNERSHIP",
    "SF1_PROTOCOL_V6_ARTIFACT_ID",
    "validate_sf1_protocol_v6",
]
