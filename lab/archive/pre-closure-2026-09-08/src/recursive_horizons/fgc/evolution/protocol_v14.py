"""Pure validation for the outcome-neutral PROTO14 calibration overlay.

PROTO14 preserves the immutable PROTO13 GR-0 calibration problem while
replacing only the retired sampled temporal admission with the independently
bound and synthetically qualified TDG6 three-level compositor.  The historical
PROTO13/CAL10 result remains unchanged.  This module performs no I/O, advances
no state, and grants no runtime, candidate, EFT, transition, or physical
authorization.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V14_ARTIFACT_ID = "FGC-2-SF1-PROTO14"

PROTO14_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO13",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v13.toml",
    "predecessor_protocol_sha256": "fc9398bad6ea4448f659188641e6da6bd7c2e295a443fe73bf73eb0308480692",
    "predecessor_freeze_artifact_id": "FGC-1-PRO13-FRZ1",
    "predecessor_freeze_result": "results/fgc-1-pro13-frz1.json",
    "predecessor_freeze_result_sha256": "a6fd1ce1bf429fd39ff2d7975de009991199a0f41c73e44d46e57b88b35d74b4",
    "historical_terminal_artifact_id": "FGC-1-CAL10-PREF15",
    "historical_terminal_result": "results/fgc-1-cal10-pref15.json",
    "historical_terminal_result_sha256": "1ebbdc7222baee8134a6d466589770ee445aee594f10630fc09aa6264cb4ac18",
    "replacement_instrument_artifact_id": "FGC-1-TDG6-IMP2",
    "replacement_instrument_config": "configs/fgc/fgc-1-tdg6-imp2.toml",
    "replacement_instrument_config_sha256": "89803b43e9243558aed33473bd9ebd010142b02302f1a5b0b380e8c4f4d7d8c0",
    "replacement_instrument_result": "results/fgc-1-tdg6-imp2.json",
    "replacement_instrument_result_sha256": "6a2cac24f65e067aed7ff7de299f9c1f0f3c3319e9aee692763b8df1db4fd5b5",
    "evidence_checkpoint_commit": "27badf9ad3c0bc3aa1c8a4153a42a26bdaf4c00f",
    "revision_classification": "premise_only_replacement_temporal_admission_rerun_revision",
    "revision_reasons": [
        "PROTO13_reached_the_first_64_sample_temporal_gate_with_spatial_and_constraint_admissions_passing",
        "TDG4_proved_that_the_sampled_rule_cannot_bound_the_declared_continuum_quantity_on_the_unrestricted_smooth_class",
        "TDG5_and_TDG6_defined_bound_and_implemented_a_finite_dimensional_same_grid_replacement",
        "TDG6_IMP2_synthetically_qualified_complete_state_all_of_admission_atomic_commit_retry_and_checkpoint_semantics",
        "historical_PROTO13_and_CAL10_results_remain_valid_and_unreclassified",
        "all_evidence_remains_GR0_numerical_premise_evidence_only",
    ],
    "revision_changes": [
        "sampled_64_history_temporal_rule_from_admission_to_nonveto_diagnostic",
        "every_accepted_post_restart_macro_step_to_TDG6_one_two_four_same_grid_admission",
        "add_componentwise_outward_temporal_debit_and_durable_temporal_retry_ledger",
        "version_output_namespace_from_proto13_to_proto14",
    ],
    "all_unlisted_PROTO13_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_and_initial_physical_data_unchanged": True,
    "physical_equations_source_backend_reference_map_and_non_temporal_stops_unchanged": True,
    "method_owned_ladders_constraints_spatial_health_boundary_and_scale_rules_unchanged": True,
    "trapped_sign_null_orientation_common_event_schedule_and_final_time_unchanged": True,
    "PROTO13_GR0_trajectory_inspected_before_revision": True,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO14_METHOD_OWNED_LADDERS = {
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
    "maximum_nested_spatial_tail_ratio": "1/4",
}

PROTO14_RESTART = {
    "branch": "GR-0",
    "amplitude": "3",
    "restart_coordinate_time": "23/16",
    "final_coordinate_time": "32",
    "primary_members_source": "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
    "comparator_4097_and_8193_source": "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
    "comparator_16385_source": "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
    "all_six_state_tracer_history_and_runtime_monitor_payloads_must_restore": True,
    "all_six_members_must_have_identical_restart_coordinate_time": True,
    "accepted_stage_counts_and_causal_debits_continue_from_bound_checkpoints": True,
    "source_and_CFL_retry_provenance_continue_without_reclassification": True,
    "no_state_reinitialization_reprojection_interpolation_or_refitting_is_permitted": True,
    "PROTO13_terminal_checkpoint_may_not_resume": True,
    "restart_arrays_are_fixed_discrete_calibration_inputs": True,
    "pre_restart_temporal_error_bound_proved": False,
    "restart_may_seed_candidate_or_physical_claim": False,
    "TDG6_ledger_initialized_at_restart": True,
    "TDG6_accumulated_debit_initial_values": [0.0] * 18,
    "TDG6_accepted_macro_step_count_initial": 0,
    "TDG6_cumulative_temporal_retry_count_initial": 0,
    "TDG6_serialized_temporal_rejection_count_initial": 0,
    "first_new_common_event_coordinate_time": "3/2",
    "common_output_interval": "1/16",
    "checkpoint_interval": "1/4",
}

PROTO14_TEMPORAL_ADMISSION = {
    "instrument_artifact_id": "FGC-1-TDG6-IMP2",
    "runtime_module": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
    "runtime_module_sha256": "039de0bd008a2536bea6f7187afc3856bfc841d14790cbae5ca2fe27f7eb8e2d",
    "methods": ["RK4", "SSPRK3"],
    "level_step_counts": [1, 2, 4],
    "channel_order": [
        f"{block}:{field}"
        for block in ("u", "p", "q")
        for field in ("alpha", "v", "lambda", "R", "phi", "chi")
    ],
    "continuous_extension": "complete_state_piecewise_cubic_Hermite",
    "minimum_observed_order": "3/2",
    "exact_outward_order_test": "8*U12^2<=L01^2",
    "passing_classes": [
        "exact_zero",
        "enclosure_dominated_debit_only",
        "resolved_order_pass",
    ],
    "admission_is_all_of_all_18_channels": True,
    "every_accepted_post_restart_macro_step_must_pass": True,
    "only_four_quarter_fine_path_may_commit": True,
    "outer_and_medium_paths_are_evidence_only": True,
    "finest_pair_debit_accumulation": "componentwise_nonnegative_outward_sum",
    "debit_cancellation_permitted": False,
    "absolute_state_tolerance_used": False,
    "physical_signal_used_for_normalization": False,
    "temporal_retry_factor": "1/2",
    "maximum_temporal_retries_per_macro_step": 32,
    "minimum_macro_step": "1/1073741824",
    "temporal_rejection_must_be_durable_before_retry": True,
    "source_only_retry_remains_PROTO7_owned": True,
    "non_source_shadow_stop_remains_terminal_at_its_owner": True,
    "temporal_retry_exhaustion_is_invalid_not_physical": True,
    "checkpoint_must_round_trip_debit_counters_and_complete_rejections": True,
    "sampled_64_history_rule_retained_as_nonveto_diagnostic": True,
    "sampled_64_history_rule_may_admit_or_reject_PROTO14": False,
    "piecewise_cubic_is_exact_PDE_history": False,
    "trajectory_asymptotic_regime_proved": False,
    "global_PDE_error_bound_proved": False,
}

PROTO14_COMMON_EVENT = {
    "common_output_interval": "1/16",
    "minimum_consecutive_qualified_trapped_common_events": 8,
    "trapped_sign_margin_over_combined_spatial_error_factor": 4,
    "fixed_future_null_orientation_inherited": True,
    "primary_and_comparator_common_admissions_are_hard_vetoes": True,
    "each_method_supplies_its_own_finest_pair_Richardson_interval": True,
    "cross_method_Richardson_intervals_must_overlap": True,
    "raw_unequal_finest_grid_collocation_forbidden": True,
    "all_prior_macro_steps_must_have_TDG6_admission": True,
    "accumulated_temporal_debit_must_be_serialized_at_every_common_event": True,
    "accumulated_temporal_debit_is_not_yet_a_Raychaudhuri_error_map": True,
    "eligible_case_is_calibration_input_to_later_premise_gates_only": True,
}

PROTO14_TERMINAL = {
    "allowed_campaign_classifications": [
        "GR0_calibration_completed_selected_amplitude",
        "calibration_failed_no_eligible_GR0_case",
        "invalid_implementation_or_nonconverged_run",
    ],
    "selected_requires_eight_consecutive_fully_qualified_events": True,
    "final_time_without_qualified_interval_is_no_eligible_GR0_case": True,
    "constraint_spatial_health_boundary_or_scale_stop_cannot_be_relabelled_as_physics": True,
    "TDG6_retry_exhaustion_classification": "invalid_implementation_or_nonconverged_run",
    "unexpected_error_classification": "invalid_implementation_or_nonconverged_run",
    "candidate_branch_may_open_from_PROTO14_result_alone": False,
    "physical_question_may_be_answered_by_PROTO14": False,
}

PROTO14_PARTITION = {
    "development_branch": "GR-0",
    "calibration_amplitude": "3",
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro14-hld1.toml",
    "resolved_holdout_manifest_present": False,
}

PROTO14_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto14/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto14/holdout",
    "campaign_manifest_name": "manifest.json",
    "append_only_event_log_name": "events.jsonl",
    "checkpoint_name": "latest-checkpoint.npz",
    "result_name": "campaign-result.json",
    "both_output_roots_must_be_absent_before_their_first_authorized_run": True,
    "PROTO13_terminal_output_is_immutable_history_not_a_PROTO14_restart": True,
    "temporal_rejections_must_enter_append_only_log_before_retry": True,
}

PROTO14_INHERITANCE = {
    "PROTO13_unredefined_GR0_REF1_equations_and_SRC4_source_backend_inherited_unchanged": True,
    "PROTO13_reference_balanced_derivative_constraint_and_reduction_maps_inherited_unchanged": True,
    "PROTO13_CFL_dissipation_source_retries_health_boundary_scale_and_transaction_rules_inherited_unchanged": True,
    "PROTO13_method_owned_spatial_and_constraint_classifiers_inherited_unchanged": True,
    "PROTO13_affine_null_orientation_trapped_sign_and_spatial_error_factor_inherited_unchanged": True,
    "PROTO13_sampled_temporal_rule_retained_only_as_diagnostic": True,
    "CAL10_terminal_result_remains_historical_and_unreclassified": True,
    "TDG6_IMP2_runtime_contract_inherited_exactly": True,
    "fresh_restart_restoration_and_time_23_over_16_admission_required_before_continuation": True,
    "successor_runtime_owner": "FGC-1-HLT12-MON12",
}

PROTO14_CLAIMS = {
    "PROTO14_successor_runtime_implemented": False,
    "PROTO14_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO14_resolved_holdout_manifest_authorized": False,
    "fresh_GR0_dynamic_calibration_completed": False,
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
    return sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def validate_sf1_protocol_v14(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete fail-closed PROTO14 overlay without I/O."""

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
        raise ValueError("PROTO14 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V14_ARTIFACT_ID,
        "protocol_version": 14,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO14 identity, version, units, purpose, or freeze differs")
    amendment = _exact("PROTO14 amendment", protocol["amendment"], PROTO14_AMENDMENT)
    replacement = _mapping("PROTO14 replacement", protocol["replacement"])
    expected_replacement = {
        "method_owned_ladders": PROTO14_METHOD_OWNED_LADDERS,
        "restart": PROTO14_RESTART,
        "temporal_admission": PROTO14_TEMPORAL_ADMISSION,
        "common_event": PROTO14_COMMON_EVENT,
        "terminal_classification": PROTO14_TERMINAL,
        "partition": PROTO14_PARTITION,
        "provenance": PROTO14_PROVENANCE,
    }
    if set(replacement) != set(expected_replacement):
        raise ValueError("PROTO14 replacement keys differ")
    validated_replacement = {
        key: dict(_exact(f"PROTO14 {key}", replacement[key], expected))
        for key, expected in expected_replacement.items()
    }
    inheritance = _exact(
        "PROTO14 inheritance", protocol["inheritance"], PROTO14_INHERITANCE
    )
    claims = _exact("PROTO14 claims", protocol["claims"], PROTO14_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO14 claims must remain fail-closed")
    contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "historical_terminal_result_sha256": amendment[
            "historical_terminal_result_sha256"
        ],
        "replacement_instrument_result_sha256": amendment[
            "replacement_instrument_result_sha256"
        ],
        "replacement": validated_replacement,
        "inheritance": dict(inheritance),
    }
    return {
        "artifact_id": SF1_PROTOCOL_V14_ARTIFACT_ID,
        "protocol_version": 14,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "replacement_temporal_admission_defined": True,
        "historical_PROTO13_result_reclassified": False,
        "restart_coordinate_time": "23/16",
        "final_coordinate_time": "32",
        "production_trajectory_authorized": False,
        "FGCQR_outcomes_inspected_before_revision": False,
        "SGBL_outcomes_inspected_before_revision": False,
        "mechanism_question_answered_before_revision": False,
        "contract_sha256": _semantic_hash(contract),
        "claims": dict(claims),
    }


__all__ = [
    "PROTO14_AMENDMENT",
    "PROTO14_CLAIMS",
    "PROTO14_COMMON_EVENT",
    "PROTO14_INHERITANCE",
    "PROTO14_METHOD_OWNED_LADDERS",
    "PROTO14_PARTITION",
    "PROTO14_PROVENANCE",
    "PROTO14_RESTART",
    "PROTO14_TEMPORAL_ADMISSION",
    "PROTO14_TERMINAL",
    "SF1_PROTOCOL_V14_ARTIFACT_ID",
    "validate_sf1_protocol_v14",
]
