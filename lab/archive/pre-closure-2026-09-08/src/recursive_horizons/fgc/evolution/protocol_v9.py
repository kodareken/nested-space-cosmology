"""Pure validation for the outcome-neutral PROTO9 resolution overlay.

PROTO9 consumes immutable PROTO8 and CAL5 evidence.  It changes only the
three-grid calibration ladder and versioned output namespace.  Every action,
initial physical input, method, equation, solver, tolerance, monitor,
constraint/spectral threshold, boundary rule, observable, and claim boundary
is inherited unchanged.  This module reads no run output and evolves no state.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V9_ARTIFACT_ID = "FGC-2-SF1-PROTO9"

PROTO9_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO8",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v8.toml",
    "predecessor_protocol_sha256": "5de88456ed041f8a4293150f49b4caea057560dd796c8ea2b9a8110ecf7770ab",
    "diagnosis_artifact_id": "FGC-1-CAL5-PREF7",
    "diagnosis_config": "configs/fgc/fgc-1-cal5-pref7.toml",
    "diagnosis_config_sha256": "02ab7bbe58d42661cbf40017acf4e7d935a82d47f3a4f5cf2ea52a03a56a1732",
    "diagnosis_result": "results/fgc-1-cal5-pref7.json",
    "diagnosis_result_sha256": "66404f44dae76021eb98268e6c48e72f17058423837177c410445d20c178b6fa",
    "diagnosis_checkpoint_commit": "a2c929e33ec8f2ee5464b086a57db8ed7f9277f0",
    "revision_classification": "premise_only_post_calibration_resolution_ladder_revision",
    "revision_reasons": [
        "PROTO8_campaign_terminated_normally_with_no_eligible_GR0_case",
        "PROTO8_common_event_ownership_repair_passed_on_evolved_amplitude_5over2_data",
        "CAL5_all_owned_constraint_norms_contract_monotonically",
        "CAL5_existing_2049_states_clear_the_unchanged_prospective_coarsest_guards",
        "CAL5_existing_4097_states_clear_all_unchanged_absolute_spectral_budgets",
        "CAL5_conditional_8193_projections_are_planning_only_not_admission",
        "CAL5_new_8193_constraint_and_tail_data_are_required",
        "CAL5_amplitude_3_phi_derivative_tail_failure_remains_public_and_unresolved",
    ],
    "revision_changes": [
        "calibration_resolution_ladder_from_1025_2049_4097_to_2049_4097_8193",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO8_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_initial_data_and_holdout_cases_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "source_solver_raw_residual_kinetic_transaction_and_retry_rules_unchanged": True,
    "methods_CFL_dissipation_and_final_time_unchanged": True,
    "constraint_ownership_magnitude_and_order_rules_unchanged": True,
    "spectral_ownership_absolute_and_nested_tail_rules_unchanged": True,
    "trapped_boundary_candidate_health_and_scale_stops_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["5/2", "3"],
    "PROTO1_through_PROTO8_preserved": True,
    "PROTO8_GR0_calibration_trajectory_inspected_before_revision": True,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO9_RESOLUTION = {
    "old_point_counts": [1025, 2049, 4097],
    "point_counts": [2049, 4097, 8193],
    "entire_nested_ladder_shifted_one_level": True,
    "coarsest_grid_remains_public_convergence_witness": True,
    "medium_and_finest_grids_own_unchanged_absolute_spectral_budgets": True,
    "all_adjacent_tail_ratios_remain_required": True,
    "coarsest_and_finest_constraint_magnitude_guards_unchanged": True,
    "minimum_finest_pair_constraint_order": "3/2",
    "conditional_Richardson_projection_is_not_admission": True,
    "new_8193_data_must_be_observed_before_any_pass": True,
    "no_old_state_may_be_relabelled_as_a_PROTO9_outcome": True,
}

PROTO9_PARTITION = {
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro9-hld1.toml",
}

PROTO9_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto9/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto9/holdout",
    "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    "PROTO8_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO9_outcomes": True,
}

PROTO9_INHERITANCE = {
    "PROTO8_common_event_constraint_ownership_inherited_unchanged": True,
    "PROTO8_common_event_spectral_ownership_inherited_unchanged": True,
    "PROTO7_general_transaction_stop_ownership_inherited_unchanged": True,
    "PROTO7_source_solver_raw_threshold_and_retry_budget_inherited_unchanged": True,
    "PROTO7_native_semidiscrete_initialization_inherited_unchanged": True,
    "PROTO8_stop_applicability_partition_inherited_except_resolution_ladder": True,
    "PROTO8_calibration_selection_and_trapped_sign_rules_inherited_unchanged": True,
    "PROTO8_boundary_health_scale_affine_and_outcome_rules_inherited_unchanged": True,
    "HLT6_runtime_evidence_remains_valid_for_unchanged_common_event_rules_only": True,
    "new_8193_inputs_and_runner_require_a_separate_HLT7_certificate": True,
    "CAL5_is_design_evidence_not_calibration_or_mechanism_evidence": True,
    "successor_runtime_owner": "FGC-1-HLT7-MON7",
}

PROTO9_CLAIMS = {
    "PROTO9_successor_runtime_implemented": False,
    "PROTO9_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO9_resolved_holdout_manifest_authorized": False,
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


def validate_sf1_protocol_v9(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete, fail-closed PROTO9 overlay without I/O."""

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
        raise ValueError("PROTO9 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V9_ARTIFACT_ID,
        "protocol_version": 9,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO9 identity, version, units, purpose, or frozen state differs")

    amendment = _exact("PROTO9 amendment", protocol["amendment"], PROTO9_AMENDMENT)
    replacement = _mapping("PROTO9 replacement", protocol["replacement"])
    if set(replacement) != {"resolution_ladder", "partition", "provenance"}:
        raise ValueError("PROTO9 replacement keys differ")
    resolution = _exact(
        "PROTO9 resolution replacement",
        replacement["resolution_ladder"],
        PROTO9_RESOLUTION,
    )
    partition = _exact(
        "PROTO9 partition replacement", replacement["partition"], PROTO9_PARTITION
    )
    provenance = _exact(
        "PROTO9 provenance replacement", replacement["provenance"], PROTO9_PROVENANCE
    )
    inheritance = _exact(
        "PROTO9 inheritance", protocol["inheritance"], PROTO9_INHERITANCE
    )
    claims = _exact("PROTO9 claims", protocol["claims"], PROTO9_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO9 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "resolution_ladder": resolution,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V9_ARTIFACT_ID,
        "protocol_version": 9,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO8",
        "diagnosis_artifact_id": "FGC-1-CAL5-PREF7",
        "diagnosis_checkpoint_commit": amendment["diagnosis_checkpoint_commit"],
        "eligible_calibration_amplitudes": ["5/2", "3"],
        "old_point_counts": resolution["old_point_counts"],
        "point_counts": resolution["point_counts"],
        "only_resolution_ladder_changes": True,
        "constraint_rules_unchanged": True,
        "spectral_rules_unchanged": True,
        "physical_and_numerical_thresholds_unchanged": True,
        "new_8193_data_required": True,
        "conditional_projection_is_admission": False,
        "PROTO8_calibration_trajectory_inspected_before_revision": True,
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
    "PROTO9_AMENDMENT",
    "PROTO9_CLAIMS",
    "PROTO9_INHERITANCE",
    "PROTO9_PARTITION",
    "PROTO9_PROVENANCE",
    "PROTO9_RESOLUTION",
    "SF1_PROTOCOL_V9_ARTIFACT_ID",
    "validate_sf1_protocol_v9",
]
