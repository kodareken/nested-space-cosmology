"""Pure validation for the outcome-neutral PROTO8 admission overlay.

PROTO8 consumes immutable PROTO7 and CAL4 evidence.  It changes only the
ownership of two common-event diagnostics: reduction constraints are judged on
the rows owned by the evolution operator, with accumulated roundoff affecting
order classification only; and the coarsest spectrum remains a public
convergence witness while the medium/fine pair owns absolute-budget admission.
This module reads no run output and performs no evolution.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V8_ARTIFACT_ID = "FGC-2-SF1-PROTO8"

PROTO8_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO7",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v7.toml",
    "predecessor_protocol_sha256": "ea0563b205bf5c6c90d8f728afc56a13f624953fffc67c6d6b9a97c03c35f3fb",
    "diagnosis_artifact_id": "FGC-1-CAL4-PREF6",
    "diagnosis_config": "configs/fgc/fgc-1-cal4-pref6.toml",
    "diagnosis_config_sha256": "687322014f16a5c48f6b3947f4cfd18a0f60f603f3c4a3abda842b19e3c8f1e2",
    "diagnosis_result": "results/fgc-1-cal4-pref6.json",
    "diagnosis_result_sha256": "e03d647e9d068b35be79ba817a3ab6c36b93f29c29de5e395a2106030d958815",
    "diagnosis_checkpoint_commit": "219ba9493899029bc51bd313e3c5f689f3809c2c",
    "revision_classification": "premise_only_post_calibration_common_event_diagnostic_ownership_revision",
    "revision_reasons": [
        "PROTO7_constraint_convergence_includes_projector_owned_outer_rows",
        "CAL4_localizes_all_old_reduction_order_failures_to_the_fixed_row_stencil_interface",
        "CAL4_evolution_owned_constraints_pass_for_amplitude_5over2_under_unchanged_guards_and_order",
        "PROTO7_requires_the_coarsest_convergence_witness_to_pass_resolved_grid_absolute_spectral_budgets",
        "CAL4_amplitude_5over2_has_a_resolved_medium_fine_pair_and_all_nested_tails_contract",
        "CAL4_amplitude_3_remains_inadmissible_under_independent_unchanged_guards",
    ],
    "revision_changes": [
        "evolution_owned_constraint_diagnostic_domain",
        "accumulated_stage_count_roundoff_enclosure_for_order_classification_only",
        "coarsest_spectral_grid_as_public_convergence_witness",
        "medium_and_finest_absolute_spectral_admission_with_all_adjacent_nested_contraction",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO7_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_initial_data_and_holdout_cases_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "source_solver_raw_residual_kinetic_and_retry_rules_unchanged": True,
    "constraint_magnitude_and_order_thresholds_unchanged": True,
    "spectral_absolute_and_nested_tail_thresholds_unchanged": True,
    "trapped_boundary_candidate_health_and_scale_stops_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["5/2", "3"],
    "PROTO1_through_PROTO7_preserved": True,
    "PROTO7_GR0_calibration_trajectory_inspected_before_revision": True,
    "first_nonzero_common_event_completed_before_revision": True,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO8_CONSTRAINT_OWNERSHIP = {
    "applies_only_to_committed_synchronized_common_events": True,
    "projector_fixed_outer_rows": 4,
    "owned_domain_excludes_fixed_rows_and_every_derivative_stencil_that_touches_them": True,
    "RK4_derivative_stencil_reach_rows": 3,
    "SSPRK3_derivative_stencil_reach_rows": 1,
    "full_domain_and_owned_domain_raw_residuals_must_be_serialized": True,
    "owned_domain_raw_norms_decide_unchanged_magnitude_guards": True,
    "coarsest_and_finest_constraint_magnitude_guards_unchanged": True,
    "minimum_finest_pair_constraint_order": "3/2",
    "roundoff_operation_budget_per_diagnostic_evaluation": 4096,
    "accumulated_roundoff_formula": "4096*binary64_epsilon*(1+accepted_stage_count)",
    "accumulated_roundoff_may_change_only_convergence_order_classification": True,
    "raw_magnitude_may_not_be_subtracted_from_physical_or_observable_margins": True,
    "interior_reduction_defect_above_the_enclosure_must_fail": True,
    "projector_owned_outer_defects_remain_public": True,
    "common_event_time_alignment_and_monotone_refinement_required": True,
}

PROTO8_SPECTRAL_OWNERSHIP = {
    "coarsest_grid_role": "public_convergence_witness",
    "coarsest_absolute_budget_failure_alone_does_not_veto": True,
    "medium_and_finest_grids_each_must_pass_every_unchanged_absolute_budget": True,
    "all_adjacent_field_and_derivative_tail_ratios_must_pass": True,
    "unresolved_top_fraction": "1/8",
    "maximum_top_band_field_power_fraction": "1/1048576",
    "maximum_top_band_derivative_power_fraction": "1/1024",
    "maximum_rms_scale_over_cutoff": "1/4",
    "maximum_nested_tail_ratio": "1/4",
    "noncontracting_tail_remains_terminal": True,
    "all_per_grid_budgets_and_tail_ratios_remain_public": True,
}

PROTO8_PARTITION = {
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro8-hld1.toml",
}

PROTO8_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto8/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto8/holdout",
    "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    "PROTO7_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO8_outcomes": True,
}

PROTO8_INHERITANCE = {
    "PROTO7_general_transaction_stop_ownership_inherited_unchanged": True,
    "PROTO7_source_solver_raw_threshold_and_retry_budget_inherited_unchanged": True,
    "PROTO7_native_semidiscrete_initialization_inherited_unchanged": True,
    "PROTO7_stop_applicability_partition_inherited_except_common_event_diagnostic_ownership": True,
    "PROTO7_calibration_selection_and_trapped_sign_rules_inherited_unchanged": True,
    "PROTO7_boundary_health_scale_affine_and_outcome_rules_inherited_unchanged": True,
    "HLT5_runtime_evidence_remains_valid_for_unchanged_rules_only": True,
    "CAL4_replay_is_design_evidence_not_calibration_or_mechanism_evidence": True,
    "successor_runtime_owner": "FGC-1-HLT6-MON6",
}

PROTO8_CLAIMS = {
    "PROTO8_successor_runtime_compositor_implemented": False,
    "PROTO8_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO8_resolved_holdout_manifest_authorized": False,
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


def validate_sf1_protocol_v8(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete, fail-closed PROTO8 overlay without I/O."""

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
        raise ValueError("PROTO8 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V8_ARTIFACT_ID,
        "protocol_version": 8,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO8 identity, version, units, purpose, or frozen state differs")

    amendment = _exact("PROTO8 amendment", protocol["amendment"], PROTO8_AMENDMENT)
    replacement = _mapping("PROTO8 replacement", protocol["replacement"])
    if set(replacement) != {
        "common_event_constraint_ownership",
        "common_event_spectral_ownership",
        "partition",
        "provenance",
    }:
        raise ValueError("PROTO8 replacement keys differ")
    constraint = _exact(
        "PROTO8 constraint-ownership replacement",
        replacement["common_event_constraint_ownership"],
        PROTO8_CONSTRAINT_OWNERSHIP,
    )
    spectrum = _exact(
        "PROTO8 spectral-ownership replacement",
        replacement["common_event_spectral_ownership"],
        PROTO8_SPECTRAL_OWNERSHIP,
    )
    partition = _exact(
        "PROTO8 partition replacement", replacement["partition"], PROTO8_PARTITION
    )
    provenance = _exact(
        "PROTO8 provenance replacement",
        replacement["provenance"],
        PROTO8_PROVENANCE,
    )
    inheritance = _exact(
        "PROTO8 inheritance", protocol["inheritance"], PROTO8_INHERITANCE
    )
    claims = _exact("PROTO8 claims", protocol["claims"], PROTO8_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO8 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "common_event_constraint_ownership": constraint,
        "common_event_spectral_ownership": spectrum,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V8_ARTIFACT_ID,
        "protocol_version": 8,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO7",
        "diagnosis_artifact_id": "FGC-1-CAL4-PREF6",
        "diagnosis_checkpoint_commit": amendment["diagnosis_checkpoint_commit"],
        "eligible_calibration_amplitudes": ["5/2", "3"],
        "projector_fixed_outer_rows": constraint["projector_fixed_outer_rows"],
        "RK4_excluded_outer_rows": (
            constraint["projector_fixed_outer_rows"]
            + constraint["RK4_derivative_stencil_reach_rows"]
        ),
        "SSPRK3_excluded_outer_rows": (
            constraint["projector_fixed_outer_rows"]
            + constraint["SSPRK3_derivative_stencil_reach_rows"]
        ),
        "accumulated_roundoff_changes_order_only": True,
        "raw_constraint_magnitude_guards_unchanged": True,
        "coarsest_spectrum_is_public_convergence_witness": True,
        "finest_pair_absolute_spectral_budgets_required": True,
        "all_adjacent_nested_tail_ratios_required": True,
        "spectral_thresholds_unchanged": True,
        "all_physical_and_numerical_thresholds_unchanged": True,
        "PROTO7_calibration_trajectory_inspected_before_revision": True,
        "first_nonzero_common_event_completed_before_revision": True,
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
    "PROTO8_AMENDMENT",
    "PROTO8_CLAIMS",
    "PROTO8_CONSTRAINT_OWNERSHIP",
    "PROTO8_INHERITANCE",
    "PROTO8_PARTITION",
    "PROTO8_PROVENANCE",
    "PROTO8_SPECTRAL_OWNERSHIP",
    "SF1_PROTOCOL_V8_ARTIFACT_ID",
    "validate_sf1_protocol_v8",
]
