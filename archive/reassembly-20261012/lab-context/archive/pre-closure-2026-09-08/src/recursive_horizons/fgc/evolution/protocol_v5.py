"""Pure validation for the outcome-neutral PROTO5 semidiscrete overlay.

PROTO5 is intentionally implemented outside the hash-bound RUN1 protocol
module.  It consumes PROTO4 as immutable predecessor evidence and changes only
the numerical premises diagnosed by CAL1-PREF3.  This module reads no run
output and performs no evolution.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V5_ARTIFACT_ID = "FGC-2-SF1-PROTO5"

PROTO5_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO4",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v4.toml",
    "predecessor_protocol_sha256": "02a25337e040b97fdd115bc7da4aa87fd09755366599bcade064644ffb3cb87a",
    "diagnosis_artifact_id": "FGC-1-CAL1-PREF3",
    "diagnosis_config": "configs/fgc/fgc-1-cal1-pref3.toml",
    "diagnosis_config_sha256": "514861c679272cd182da9674e7640734a427c0349157c03f0899cbe439b475ff",
    "diagnosis_result": "results/fgc-1-cal1-pref3.json",
    "diagnosis_result_sha256": "3ea43baf0fb9830b348fc49a79af5c535d8f7089aa03a957a22f6927528a4395",
    "diagnosis_checkpoint_commit": "0781d7d218e9c129858b134c3181fad07ac692b0",
    "revision_classification": "premise_based_pre_trajectory_semidiscrete_contract_revision",
    "revision_reasons": [
        "PROTO4_omits_continuum_to_semidiscrete_q_initialization",
        "PROTO4_uses_one_constraint_magnitude_budget_for_different_method_orders",
        "PROTO4_exact_zero_only_convergence_logic_misclassifies_binary64_roundoff",
    ],
    "revision_changes": [
        "native_SBP_q_initialization_before_first_stage",
        "method_owned_raw_constraint_magnitude_guards",
        "public_roundoff_zero_enclosure_for_convergence_classification_only",
        "monotone_three_grid_refinement_with_finest_pair_order",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO4_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_and_holdout_cases_unchanged": True,
    "physical_u_and_p_initial_data_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["5/2", "3"],
    "PROTO1_through_PROTO4_preserved": True,
    "initial_common_event_diagnostics_inspected_before_revision": True,
    "fresh_GR0_calibration_trajectories_inspected_before_revision": False,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
}

PROTO5_INITIAL_STATE = {
    "physical_u_source": "FGC-1-ID2-ALL1_complete_grid_state",
    "physical_p_source": "FGC-1-ID2-ALL1_complete_grid_state",
    "physical_u_and_p_must_be_bitwise_preserved": True,
    "auxiliary_q_initialization": "q_equals_native_SBP_D_h_u_before_first_stage",
    "projection_applies_only_to_the_kinematic_auxiliary_q": True,
    "projection_is_a_discrete_reduction_constraint_initialization_not_a_physical_field_change": True,
    "primary_native_spatial_operator": "D4_2_SBP",
    "comparator_native_spatial_operator": "D2_1_SBP",
    "projection_must_be_hash_bound_at_every_resolution": True,
}

PROTO5_CONSTRAINT_ADMISSION = {
    "fields": [
        "Hamiltonian",
        "radial_momentum",
        "gauge_t",
        "gauge_r",
        "six_kinematic_reduction_fields",
    ],
    "normalization": "max_absolute_unredefined_projection_term_sum_or_M_Pl_squared_over_L0_squared_or_1e-30_componentwise",
    "common_event_alignment": "same_coordinate_time_for_GR0_calibration_and_same_affine_label_and_origin_for_DEF1",
    "raw_normalized_residuals_decide_all_magnitude_guards": True,
    "roundoff_zero_enclosure_formula": "4096_times_binary64_machine_epsilon",
    "roundoff_operation_budget": 4096,
    "roundoff_enclosure_use": "convergence_zero_classification_only_and_public_nonnegative_error_component",
    "roundoff_enclosure_must_not_be_subtracted_from_any_observable_or_residual": True,
    "primary_method": "RK4",
    "primary_spatial_order": 4,
    "primary_coarsest_normalized_guard_max": "1/50",
    "primary_finest_normalized_guard_max": "1/1000",
    "comparator_method": "SSPRK3",
    "comparator_spatial_order": 2,
    "comparator_coarsest_normalized_guard_max": "1/10",
    "comparator_finest_normalized_guard_max": "1/200",
    "minimum_nested_resolutions": 3,
    "all_three_resolutions_must_refine_monotonically": True,
    "convergence_order_pair": "finest_adjacent_pair_after_monotone_three_grid_refinement",
    "minimum_finest_adjacent_pair_order": "3/2",
    "coarse_to_middle_order_is_serialized_diagnostic_not_admission_boolean": True,
    "Richardson_error_enters_conservative_observable_error_sum": True,
    "constraint_error_must_not_be_subtracted_from_observable_margin": True,
    "single_resolution_smallness_cannot_establish_admission": True,
    "propagation_theorem_not_inferred_from_numerical_smallness": True,
}

PROTO5_PARTITION = {
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro5-hld1.toml",
}

PROTO5_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto5/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto5/holdout",
    "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    "earlier_proto4_outputs_if_any_are_development_only_and_forbidden_as_PROTO5_evidence": True,
}

PROTO5_INHERITANCE = {
    "PROTO4_stop_applicability_partition_inherited_unchanged": True,
    "PROTO4_weighted_spectral_admission_inherited_unchanged": True,
    "PROTO4_calibration_selection_and_trapped_sign_rules_inherited_unchanged": True,
    "HLT2_formula_evidence_remains_valid_for_unchanged_rules": True,
    "ID2_static_candidate_ledger_remains_valid": True,
    "CAL1_prospective_repair_is_design_evidence_not_calibration_evidence": True,
    "successor_runtime_owner": "FGC-1-HLT3-MON3",
}

PROTO5_CLAIMS = {
    "PROTO5_successor_runtime_compositor_implemented": False,
    "PROTO5_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO5_resolved_holdout_manifest_authorized": False,
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


def validate_sf1_protocol_v5(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete, fail-closed PROTO5 overlay without I/O."""

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
        raise ValueError("PROTO5 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V5_ARTIFACT_ID,
        "protocol_version": 5,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO5 identity, version, units, purpose, or frozen state differs")

    amendment = _exact("PROTO5 amendment", protocol["amendment"], PROTO5_AMENDMENT)
    replacement = _mapping("PROTO5 replacement", protocol["replacement"])
    if set(replacement) != {"initial_state", "constraint_admission", "partition", "provenance"}:
        raise ValueError("PROTO5 replacement keys differ")
    initial_state = _exact(
        "PROTO5 initial-state replacement",
        replacement["initial_state"],
        PROTO5_INITIAL_STATE,
    )
    constraints = _exact(
        "PROTO5 constraint-admission replacement",
        replacement["constraint_admission"],
        PROTO5_CONSTRAINT_ADMISSION,
    )
    partition = _exact(
        "PROTO5 partition replacement",
        replacement["partition"],
        PROTO5_PARTITION,
    )
    provenance = _exact(
        "PROTO5 provenance replacement",
        replacement["provenance"],
        PROTO5_PROVENANCE,
    )
    inheritance = _exact(
        "PROTO5 inheritance",
        protocol["inheritance"],
        PROTO5_INHERITANCE,
    )
    claims = _exact("PROTO5 claims", protocol["claims"], PROTO5_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO5 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "initial_state": initial_state,
        "constraint_admission": constraints,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V5_ARTIFACT_ID,
        "protocol_version": 5,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO4",
        "diagnosis_artifact_id": "FGC-1-CAL1-PREF3",
        "diagnosis_checkpoint_commit": amendment["diagnosis_checkpoint_commit"],
        "eligible_calibration_amplitudes": ["5/2", "3"],
        "native_semidiscrete_q_initialization_validated": True,
        "physical_u_and_p_unchanged": True,
        "method_owned_raw_constraint_guards_validated": True,
        "roundoff_zero_classification_only_validated": True,
        "monotone_three_grid_finest_pair_order_validated": True,
        "fresh_calibration_trajectories_inspected_before_revision": False,
        "FGCQR_outcomes_inspected_before_revision": False,
        "SGBL_outcomes_inspected_before_revision": False,
        "resolved_holdout_manifest": partition["resolved_holdout_manifest"],
        "calibration_output_root": provenance["calibration_output_root"],
        "holdout_output_root": provenance["holdout_output_root"],
        "resolved_holdout_manifest_present": False,
        "semantic_holdout_contract_sha256": _semantic_hash(overlay_contract),
        "claims": dict(claims),
    }
