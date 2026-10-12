"""Pure validation for the outcome-neutral PROTO12 calibration overlay.

PROTO12 consumes immutable PROTO11, CAL8/PREF10, RSP1/PREF12, and SRC4/VEC1
evidence.  It changes only two numerical premises: the unchanged GR-0 REF1
source is evaluated through SRC4's independently certified tensor-contracted
backend, and a failed raw nested-tail ratio on either adjacent grid pair may
be labelled diagnostically saturated under the complete pairwise guard set.

The raw ratios, absolute budgets, equations, physical inputs, grid ladder,
methods, thresholds, stops, observables, and claims remain public and fixed.
This module reads no run output and evolves no state.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V12_ARTIFACT_ID = "FGC-2-SF1-PROTO12"

PROTO12_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO11",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v11.toml",
    "predecessor_protocol_sha256": "31d21bf6257d6afe033ea8a0d48b59360187a1d9baf2bc95c2fcb086e39ce440",
    "calibration_diagnosis_artifact_id": "FGC-1-CAL8-PREF10",
    "calibration_diagnosis_config": "configs/fgc/fgc-1-cal8-pref10.toml",
    "calibration_diagnosis_config_sha256": "a7646fb7ef2b229c36981a2f4a87e4428de1bdc5cbe324baedc965097c674cb3",
    "calibration_diagnosis_result": "results/fgc-1-cal8-pref10.json",
    "calibration_diagnosis_result_sha256": "10c6898e1ec6e7c07793180c4a4f8fa26a1b2f660115a0b2bc54497ff8bbb89d",
    "resolution_evidence_artifact_id": "FGC-1-RSP1-PREF12",
    "resolution_evidence_config": "configs/fgc/fgc-1-rsp1-pref12.toml",
    "resolution_evidence_config_sha256": "10248e01e83714a3e0b83d213748be8d5c6e1b234133b7bbe35a2809515e5bd1",
    "resolution_evidence_result": "results/fgc-1-rsp1-pref12.json",
    "resolution_evidence_result_sha256": "0bb184a519230d1b221d20fb88df70841800575b886df0b0642837deb26e4551",
    "source_instrument_artifact_id": "FGC-1-SRC4-VEC1",
    "source_instrument_config": "configs/fgc/fgc-1-src4-vec1.toml",
    "source_instrument_config_sha256": "39c1c11c703c58cf884e4f0f240db0dd1eaa085494e6d2cb14761fbcc48681d2",
    "source_instrument_result": "results/fgc-1-src4-vec1.json",
    "source_instrument_result_sha256": "fb510b0232ac956a73321b63a31d1f52a7732132ff736d41f6d08a4da125aa7a",
    "evidence_checkpoint_commit": "ffcc445354dc2bb2ee09afc3296f35ac8e473f46",
    "revision_classification": (
        "premise_only_post_calibration_source_instrument_and_pairwise_spectral_revision"
    ),
    "revision_reasons": [
        "PROTO11_completed_without_an_eligible_GR0_case",
        "CAL8_separated_a_binary64_source_wall_from_an_independent_direct_coarse_spectral_veto",
        "PREF11_proved_the_captured_source_wall_was_floating_cancellation_not_a_singular_exact_system",
        "SRC3_cleared_the_unchanged_full_grid_source_gate",
        "RSP1_all_six_amplitude_three_members_reached_the_independent_high_ladder_endpoint_without_source_retry",
        "RSP1_both_methods_cleared_the_direct_phi_derivative_target_on_both_adjacent_pairs",
        "RSP1_generic_all_field_raw_admission_remained_false_and_public",
        "RSP1_pairwise_recomputation_found_every_failed_raw_ratio_below_both_measured_map_scales_with_absolute_budget_and_profile_guards",
        "SRC4_preserved_the_SRC3_equations_branch_and_strict_gates_under_independent_exact_and_differential_controls",
        "all_evidence_is_GR0_numerical_premise_evidence_only",
    ],
    "revision_changes": [
        "GR0_REF1_source_backend_from_SRC3_explicit_index_loops_to_certified_SRC4_tensor_contractions",
        "nested_tail_admission_from_coarse_direct_plus_finest_direct_or_saturated_to_each_pair_direct_or_saturated",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO11_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_initial_data_and_resolution_ladder_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "raw_source_residual_kinetic_transaction_and_retry_thresholds_unchanged": True,
    "methods_CFL_dissipation_final_time_and_common_event_schedule_unchanged": True,
    "constraint_ownership_magnitude_and_order_rules_unchanged": True,
    "absolute_spectral_budgets_and_direct_ratio_ceiling_unchanged": True,
    "trapped_boundary_candidate_health_and_scale_stops_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["5/2", "3"],
    "calibration_point_counts_in_frozen_order": [2049, 4097, 8193],
    "RSP1_high_ladder_is_qualification_evidence_not_the_PROTO12_calibration_ladder": True,
    "PROTO1_through_PROTO11_preserved": True,
    "PROTO11_GR0_calibration_trajectory_inspected_before_revision": True,
    "RSP1_GR0_resolution_trajectory_inspected_before_revision": True,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO12_SOURCE_EVALUATOR = {
    "branch": "GR-0",
    "continuum_equations": "unredefined_ACT1_VAR1_GR0_REF1",
    "certifying_artifact": "FGC-1-SRC4-VEC1",
    "implementation": "src/recursive_horizons/fgc/evolution/src4_vectorized_reference_source.py",
    "backend": "NumPy_tensor_contractions_over_only_four_dimensional_spacetime_indices",
    "SRC3_explicit_index_backend_remains_the_independent_differential_oracle": True,
    "reference_state_derivative_map_remains_PROTO11": True,
    "complete_raw_residual_is_evaluated_before_any_exact_reference_identity_branch": True,
    "raw_source_residual_maximum": "1/1000000000000",
    "kinetic_condition_number_maximum": "10000000000",
    "maximum_complete_refinement_iterations": 16,
    "point_batch_size_maximum": 2048,
    "fixed_point_batching_may_not_change_any_pointwise_affine_system": True,
    "equation_reference_geometry_threshold_and_root_branch_changes_forbidden": True,
    "higher_precision_fallback_epsilon_floor_and_fitted_coefficients_forbidden": True,
    "wall_clock_speed_is_not_a_scientific_gate": True,
}

PROTO12_SPECTRAL_INTERPRETATION = {
    "calibration_point_counts": [2049, 4097, 8193],
    "maximum_nested_tail_ratio": "1/4",
    "adjacent_pairs": [[2049, 4097], [4097, 8193]],
    "each_field_and_derivative_pair_rule": (
        "direct_ratio_below_one_quarter_or_conditioning_guarded_diagnostic_saturation"
    ),
    "direct_ratio_remains_primary_on_every_pair": True,
    "every_raw_power_fraction_and_ratio_must_be_serialized": True,
    "medium_and_fine_individual_absolute_budgets_remain_globally_required": True,
    "a_saturated_pair_requires_both_pair_member_absolute_budgets": True,
    "a_saturated_pair_requires_both_pair_member_top_band_erasures_within_their_own_measured_map_scales": True,
    "all_three_round_trip_residuals_must_contract_strictly_or_be_exact_zero": True,
    "whole_profile_infinity_and_RMS_differences_must_contract_strictly_or_be_exact_zero": True,
    "all_map_residuals_tail_amplitudes_and_erasure_witnesses_must_be_serialized": True,
    "tail_erasure_is_an_orthogonal_FFT_diagnostic_input_witness": True,
    "round_trip_interpolation_diagnostic_is_not_a_continuum_error_bound": True,
    "diagnostic_saturation_is_not_physical_resolution_or_mechanism_evidence": True,
    "RSP1_generic_all_field_raw_failure_must_remain_public": True,
    "no_epsilon_floor_new_tolerance_fitted_asymptote_or_hidden_ratio_is_permitted": True,
}

PROTO12_PARTITION = {
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro12-hld1.toml",
}

PROTO12_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto12/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto12/holdout",
    "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    "PROTO11_and_RSP1_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO12_outcomes": True,
}

PROTO12_INHERITANCE = {
    "PROTO11_reference_state_differentiation_map_inherited_unchanged": True,
    "PROTO11_resolution_ladder_amplitudes_methods_and_schedule_inherited_unchanged": True,
    "PROTO11_constraint_physical_transaction_boundary_and_stop_rules_inherited_unchanged": True,
    "RSP1_high_ladder_target_pass_remains_scoped_design_evidence_only": True,
    "RSP1_generic_raw_all_field_failure_remains_public": True,
    "CAL8_RSP1_and_SRC4_are_numerical_premise_evidence_not_mechanism_evidence": True,
    "fresh_states_must_be_rebuilt_and_source_checked_under_SRC4": True,
    "fresh_evolved_states_must_reassess_every_source_constraint_pairwise_spectral_health_and_boundary_gate": True,
    "successor_runtime_owner": "FGC-1-HLT10-MON10",
}

PROTO12_CLAIMS = {
    "PROTO12_successor_runtime_implemented": False,
    "PROTO12_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO12_resolved_holdout_manifest_authorized": False,
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


def validate_sf1_protocol_v12(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete fail-closed PROTO12 overlay without I/O."""

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
        raise ValueError("PROTO12 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V12_ARTIFACT_ID,
        "protocol_version": 12,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO12 identity, version, units, purpose, or freeze differs")

    amendment = _exact("PROTO12 amendment", protocol["amendment"], PROTO12_AMENDMENT)
    replacement = _mapping("PROTO12 replacement", protocol["replacement"])
    if set(replacement) != {
        "source_evaluator",
        "spatial_spectral_interpretation",
        "partition",
        "provenance",
    }:
        raise ValueError("PROTO12 replacement keys differ")
    source = _exact(
        "PROTO12 source evaluator",
        replacement["source_evaluator"],
        PROTO12_SOURCE_EVALUATOR,
    )
    spectral = _exact(
        "PROTO12 spectral interpretation",
        replacement["spatial_spectral_interpretation"],
        PROTO12_SPECTRAL_INTERPRETATION,
    )
    partition = _exact(
        "PROTO12 partition", replacement["partition"], PROTO12_PARTITION
    )
    provenance = _exact(
        "PROTO12 provenance", replacement["provenance"], PROTO12_PROVENANCE
    )
    inheritance = _exact(
        "PROTO12 inheritance", protocol["inheritance"], PROTO12_INHERITANCE
    )
    claims = _exact("PROTO12 claims", protocol["claims"], PROTO12_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO12 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "source_evaluator": source,
        "spatial_spectral_interpretation": spectral,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V12_ARTIFACT_ID,
        "protocol_version": 12,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO11",
        "evidence_checkpoint_commit": amendment["evidence_checkpoint_commit"],
        "evidence_artifact_ids": [
            amendment["calibration_diagnosis_artifact_id"],
            amendment["resolution_evidence_artifact_id"],
            amendment["source_instrument_artifact_id"],
        ],
        "eligible_calibration_amplitudes": ["5/2", "3"],
        "point_counts": [2049, 4097, 8193],
        "RSP1_qualification_point_counts": [4097, 8193, 16385],
        "resolution_ladder_unchanged": True,
        "source_backend": source["backend"],
        "source_equations_thresholds_and_branch_unchanged": True,
        "pairwise_spectral_rule": spectral[
            "each_field_and_derivative_pair_rule"
        ],
        "direct_ratio_ceiling": spectral["maximum_nested_tail_ratio"],
        "raw_ratios_remain_public": True,
        "round_trip_is_a_continuum_error_bound": False,
        "diagnostic_saturation_is_physical_resolution_evidence": False,
        "new_numeric_tolerance_added": False,
        "physical_and_numerical_thresholds_unchanged": True,
        "PROTO11_GR0_calibration_inspected_before_revision": True,
        "RSP1_GR0_resolution_trajectory_inspected_before_revision": True,
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
    "PROTO12_AMENDMENT",
    "PROTO12_CLAIMS",
    "PROTO12_INHERITANCE",
    "PROTO12_PARTITION",
    "PROTO12_PROVENANCE",
    "PROTO12_SOURCE_EVALUATOR",
    "PROTO12_SPECTRAL_INTERPRETATION",
    "SF1_PROTOCOL_V12_ARTIFACT_ID",
    "validate_sf1_protocol_v12",
]
