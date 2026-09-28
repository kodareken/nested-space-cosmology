"""Pure validation for the outcome-neutral PROTO13 calibration overlay.

PROTO13 consumes the immutable PROTO12 contract and the prospectively frozen
RSP2/PREF14 result.  It changes only the selected GR-0 calibration amplitude,
the method-owned resolution ladders, the restart source at ``t=23/16``, and
the versioned output namespace.  It reads no trajectory and grants no runtime,
candidate, EFT, transition, or physical authorization.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V13_ARTIFACT_ID = "FGC-2-SF1-PROTO13"

PROTO13_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO12",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v12.toml",
    "predecessor_protocol_sha256": "2c33ef4a703e918683c886eb216bfc77af5a48d477c1e1556e0b5176d44fbe7f",
    "resolution_evidence_artifact_id": "FGC-1-RSP2-PREF14",
    "resolution_evidence_config": "configs/fgc/fgc-1-rsp2-pref14.toml",
    "resolution_evidence_config_sha256": "d9d9a6ebae797f5f59646b7ae7d998740516c775038e2e2ff356577b35373b42",
    "resolution_evidence_result": "results/fgc-1-rsp2-pref14.json",
    "resolution_evidence_result_sha256": "dd09cbfff2f3f22729adc6580fd577ab183a68960fedaefcc4e3a90d0bc1b0c1",
    "calibration_diagnosis_artifact_id": "FGC-1-CAL9-PREF13",
    "calibration_diagnosis_result": "results/fgc-1-cal9-pref13.json",
    "calibration_diagnosis_result_sha256": "69a2aeff732081d1d1be2e6fd49d020335e2fa13f737d5e21484d080fb1f5fd2",
    "evidence_checkpoint_commit": "ef5c78eb5c981883ea47dacb7a968bb660eff4dd",
    "revision_classification": "premise_only_post_RSP2_method_owned_ladder_and_restart_revision",
    "revision_reasons": [
        "CAL9_amplitude_three_missed_only_the_SSPRK3_radial_momentum_order_gate_on_4097_to_8193",
        "RSP2_prospectively_added_exactly_one_16385_SSPRK3_member",
        "RSP2_independently_reconstructed_8193_to_16385_order_1_9817436463819333",
        "RSP2_complete_constraint_admission_passed_with_minimum_finite_order_1_8428133958883794",
        "RSP2_had_zero_source_and_CFL_retries_and_no_typed_stop",
        "all_evidence_is_GR0_numerical_premise_evidence_only",
    ],
    "revision_changes": [
        "eligible_calibration_amplitudes_from_5_over_2_and_3_to_3_only",
        "one_common_resolution_ladder_to_method_owned_resolution_ladders",
        "fresh_time_zero_states_to_hash_bound_restart_states_at_t_23_over_16",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO12_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_and_initial_physical_data_unchanged": True,
    "physical_equations_source_backend_reference_map_and_stop_rules_unchanged": True,
    "methods_CFL_dissipation_retry_limits_final_time_and_common_event_schedule_unchanged": True,
    "constraint_ownership_magnitude_order_and_roundoff_rules_unchanged": True,
    "absolute_spectral_budgets_and_direct_ratio_ceiling_unchanged": True,
    "trapped_sign_null_orientation_boundary_health_and_scale_rules_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["3"],
    "primary_method_point_counts_in_frozen_order": [2049, 4097, 8193],
    "comparator_method_point_counts_in_frozen_order": [4097, 8193, 16385],
    "restart_coordinate_time": "23/16",
    "final_coordinate_time": "32",
    "minimum_consecutive_trapped_common_events": 8,
    "trapped_sign_margin_over_combined_error_factor": 4,
    "PROTO12_and_RSP2_GR0_trajectories_inspected_before_revision": True,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO13_METHOD_OWNED_LADDERS = {
    "primary_method": "RK4",
    "primary_point_counts": [2049, 4097, 8193],
    "primary_adjacent_pairs": [[2049, 4097], [4097, 8193]],
    "comparator_method": "SSPRK3",
    "comparator_point_counts": [4097, 8193, 16385],
    "comparator_adjacent_pairs": [[4097, 8193], [8193, 16385]],
    "each_method_must_independently_pass_its_complete_constraint_and_spatial_error_budget": True,
    "no_grid_from_one_method_may_substitute_for_a_missing_grid_from_the_other": True,
    "cross_method_comparison_may_not_assume_equal_finest_grid_spacing": True,
    "minimum_constraint_finest_pair_order": "3/2",
    "maximum_nested_tail_ratio": "1/4",
    "RSP2_preserved_coarse_pair_miss_remains_public": True,
    "RSP2_finer_pair_pass_is_design_evidence_not_fresh_calibration": True,
}

PROTO13_RESTART = {
    "branch": "GR-0",
    "amplitude": "3",
    "restart_coordinate_time": "23/16",
    "final_coordinate_time": "32",
    "primary_members_source": "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
    "comparator_4097_and_8193_source": "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
    "comparator_16385_source": "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
    "all_six_state_tracer_history_and_runtime_monitor_payloads_must_restore": True,
    "all_six_members_must_have_identical_restart_coordinate_time": True,
    "accepted_stage_counts_and_causal_debits_continue_from_the_bound_checkpoints": True,
    "source_and_CFL_retry_provenance_continue_without_reclassification": True,
    "no_state_reinitialization_reprojection_interpolation_or_refitting_is_permitted": True,
    "first_new_common_event_coordinate_time": "3/2",
    "common_output_interval": "1/16",
    "checkpoint_interval": "1/4",
}

PROTO13_CROSS_METHOD = {
    "comparison": "Richardson_or_interval_enclosed_observables_at_common_physical_radii",
    "primary_finest_spacing_may_differ_from_comparator_finest_spacing": True,
    "raw_finest_grid_values_may_not_be_compared_as_if_collocated": True,
    "each_method_must_supply_its_own_finest_pair_Richardson_enclosure": True,
    "cross_method_agreement_uses_overlap_of_extrapolated_or_interval_enclosed_observables": True,
    "fixed_future_null_orientation_inherited": True,
    "trapped_sign_margin_over_combined_error_factor": 4,
    "minimum_consecutive_qualified_trapped_common_events": 8,
    "constraint_and_spectral_admissions_remain_vetoes_not_Raychaudhuri_error_conversions": True,
}

PROTO13_PARTITION = {
    "development_branch": "GR-0",
    "calibration_amplitude": "3",
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro13-hld1.toml",
    "resolved_holdout_manifest_present": False,
}

PROTO13_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto13/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto13/holdout",
    "both_output_roots_must_be_absent_before_their_first_authorized_run": True,
    "PROTO12_and_RSP2_outputs_are_immutable_restart_evidence_not_PROTO13_outcomes": True,
}

PROTO13_INHERITANCE = {
    "PROTO12_unredefined_GR0_REF1_equations_and_SRC4_source_backend_inherited_unchanged": True,
    "PROTO12_reference_balanced_derivative_constraint_and_reduction_maps_inherited_unchanged": True,
    "PROTO12_CFL_dissipation_retries_health_boundary_scale_and_transaction_rules_inherited_unchanged": True,
    "PROTO12_raw_plus_pairwise_guarded_spatial_classifier_inherited_per_method": True,
    "PROTO12_complete_constraint_classifier_inherited_per_method": True,
    "PROTO12_affine_null_orientation_trapped_sign_and_error_factor_inherited_unchanged": True,
    "RSP2_raw_norm_enclosure_order_and_nonclaim_ledgers_remain_public": True,
    "fresh_restart_restoration_and_time_23_over_16_admission_required_before_continuation": True,
    "successor_runtime_owner": "FGC-1-HLT11-MON11",
}

PROTO13_CLAIMS = {
    "PROTO13_successor_runtime_implemented": False,
    "PROTO13_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO13_resolved_holdout_manifest_authorized": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "GR0_case_eligible": False,
    "classical_spherical_diagnostic_authorized": False,
    "SGBL_execution_authorized": False,
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
    answer = _mapping(name, value)
    if dict(answer) != dict(expected):
        raise ValueError(f"{name} differs")
    return answer


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


def validate_sf1_protocol_v13(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete fail-closed PROTO13 overlay without I/O."""

    protocol = _mapping("protocol", protocol)
    if set(protocol) != {
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
    }:
        raise ValueError("PROTO13 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V13_ARTIFACT_ID,
        "protocol_version": 13,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO13 identity, version, units, purpose, or freeze differs")
    amendment = _exact("PROTO13 amendment", protocol["amendment"], PROTO13_AMENDMENT)
    replacement = _mapping("PROTO13 replacement", protocol["replacement"])
    if set(replacement) != {
        "method_owned_ladders",
        "restart",
        "cross_method_observables",
        "partition",
        "provenance",
    }:
        raise ValueError("PROTO13 replacement keys differ")
    ladders = _exact(
        "PROTO13 method-owned ladders",
        replacement["method_owned_ladders"],
        PROTO13_METHOD_OWNED_LADDERS,
    )
    restart = _exact("PROTO13 restart", replacement["restart"], PROTO13_RESTART)
    cross_method = _exact(
        "PROTO13 cross-method observables",
        replacement["cross_method_observables"],
        PROTO13_CROSS_METHOD,
    )
    partition = _exact(
        "PROTO13 partition", replacement["partition"], PROTO13_PARTITION
    )
    provenance = _exact(
        "PROTO13 provenance", replacement["provenance"], PROTO13_PROVENANCE
    )
    inheritance = _exact(
        "PROTO13 inheritance", protocol["inheritance"], PROTO13_INHERITANCE
    )
    claims = _exact("PROTO13 claims", protocol["claims"], PROTO13_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO13 claims must remain fail-closed")
    contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "resolution_evidence_result_sha256": amendment[
            "resolution_evidence_result_sha256"
        ],
        "method_owned_ladders": ladders,
        "restart": restart,
        "cross_method_observables": cross_method,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V13_ARTIFACT_ID,
        "protocol_version": 13,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "eligible_calibration_amplitudes": ["3"],
        "primary_point_counts": [2049, 4097, 8193],
        "comparator_point_counts": [4097, 8193, 16385],
        "restart_coordinate_time": "23/16",
        "final_coordinate_time": "32",
        "FGCQR_outcomes_inspected_before_revision": False,
        "SGBL_outcomes_inspected_before_revision": False,
        "mechanism_question_answered_before_revision": False,
        "contract_sha256": _semantic_hash(contract),
        "claims": dict(claims),
    }


__all__ = [
    "PROTO13_AMENDMENT",
    "PROTO13_CLAIMS",
    "PROTO13_CROSS_METHOD",
    "PROTO13_INHERITANCE",
    "PROTO13_METHOD_OWNED_LADDERS",
    "PROTO13_PARTITION",
    "PROTO13_PROVENANCE",
    "PROTO13_RESTART",
    "SF1_PROTOCOL_V13_ARTIFACT_ID",
    "validate_sf1_protocol_v13",
]
