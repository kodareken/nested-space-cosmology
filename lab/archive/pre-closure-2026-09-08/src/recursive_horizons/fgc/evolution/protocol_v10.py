"""Pure validation for the outcome-neutral PROTO10 spectral overlay.

PROTO10 consumes immutable PROTO9 and CAL6/PREF8 evidence.  It preserves the
complete resolution ladder, direct ratio ceiling, absolute budgets, action,
initial data, evolution method, solver, transaction, stop, observable, and
claim boundary.  Its only scientific amendment is a conditioning-guarded
interpretation of a failed finest-pair nested-tail ratio.  This module reads no
run output and evolves no state.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V10_ARTIFACT_ID = "FGC-2-SF1-PROTO10"

PROTO10_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO9",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v9.toml",
    "predecessor_protocol_sha256": "81ce8426b006518d58accae50e43fcbfe942d013dfd010a716ff265b57afb48e",
    "diagnosis_artifact_id": "FGC-1-CAL6-PREF8",
    "diagnosis_config": "configs/fgc/fgc-1-cal6-pref8.toml",
    "diagnosis_config_sha256": "940238188ac8c76527165be6d48fcd1ef0e9b84565a9fe7a7e53b12a8e8559ea",
    "diagnosis_result": "results/fgc-1-cal6-pref8.json",
    "diagnosis_result_sha256": "3d4624cca60363c9daf6f65f9393d959d7a7fafbcf6eac7f8dc9c2f237c24f5e",
    "diagnosis_checkpoint_commit": "dd3f5adb22f08aa7950b74cef2595aaca2e8cd6e",
    "revision_classification": (
        "premise_only_post_preflight_finest_pair_spectral_interpretation_revision"
    ),
    "revision_reasons": [
        "PROTO9_preflight_constructed_all_twelve_direct_states_without_opening_a_trajectory",
        "HLT7_all_raw_source_constraint_and_absolute_spectrum_gates_passed",
        "HLT7_all_four_t0_admissions_failed_only_the_finest_pair_nested_tail_ratio",
        "CAL6_preserved_all_sixteen_raw_finest_pair_failures",
        "CAL6_every_coarse_pair_ratio_passed_directly",
        "CAL6_every_medium_and_fine_absolute_spectrum_budget_passed_unchanged",
        "CAL6_all_failed_finest_pair_ratios_are_conditioned_below_the_measured_map_non_idempotence",
        "CAL6_complete_profiles_and_round_trip_residuals_contract_strictly",
        "CAL6_guarded_saturation_is_design_evidence_only",
    ],
    "revision_changes": [
        "finest_pair_nested_tail_admission_from_direct_only_to_direct_or_conditioning_guarded_saturation",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO9_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_initial_data_and_resolution_ladder_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "source_solver_raw_residual_kinetic_transaction_and_retry_rules_unchanged": True,
    "methods_CFL_dissipation_and_final_time_unchanged": True,
    "constraint_ownership_magnitude_and_order_rules_unchanged": True,
    "absolute_spectral_budgets_and_direct_ratio_ceiling_unchanged": True,
    "coarse_to_medium_direct_ratio_veto_unchanged": True,
    "trapped_boundary_candidate_health_and_scale_stops_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["5/2", "3"],
    "PROTO1_through_PROTO9_preserved": True,
    "PROTO9_GR0_trajectory_inspected_before_revision": False,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO10_SPECTRAL_INTERPRETATION = {
    "point_counts": [2049, 4097, 8193],
    "maximum_nested_tail_ratio": "1/4",
    "coarse_to_medium_field_and_derivative_ratios_must_pass_directly": True,
    "medium_and_fine_individual_absolute_budgets_must_pass": True,
    "finest_pair_rule": (
        "direct_ratio_below_one_quarter_or_conditioning_guarded_diagnostic_saturation"
    ),
    "finest_pair_direct_ratio_remains_primary": True,
    "failed_finest_pair_ratio_may_only_use_guarded_saturation": True,
    "medium_and_fine_top_band_erasure_must_fit_inside_each_measured_round_trip_scale": True,
    "all_three_round_trip_residuals_must_contract_strictly_or_be_exact_zero": True,
    "whole_profile_infinity_and_RMS_differences_must_contract_strictly_or_be_exact_zero": True,
    "all_raw_power_fractions_ratios_residuals_and_erasure_witnesses_must_be_serialized": True,
    "tail_erasure_is_an_orthogonal_FFT_diagnostic_input_witness": True,
    "round_trip_interpolation_diagnostic_is_not_a_continuum_error_bound": True,
    "diagnostic_saturation_is_not_physical_resolution_or_mechanism_evidence": True,
    "no_epsilon_floor_or_new_numeric_tolerance_is_permitted": True,
}

PROTO10_PARTITION = {
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro10-hld1.toml",
}

PROTO10_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto10/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto10/holdout",
    "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    "PROTO9_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO10_outcomes": True,
}

PROTO10_INHERITANCE = {
    "PROTO9_resolution_ladder_inherited_unchanged": True,
    "PROTO9_constraint_ownership_magnitude_and_order_rules_inherited_unchanged": True,
    "PROTO9_absolute_spectral_budgets_and_direct_ratio_ceiling_inherited_unchanged": True,
    "PROTO9_coarse_to_medium_direct_ratio_veto_inherited_unchanged": True,
    "PROTO9_physical_numerical_transaction_and_stop_rules_inherited_unchanged": True,
    "HLT7_direct_state_hashes_and_raw_admissions_remain_valid_as_diagnostic_inputs_only": True,
    "CAL6_is_design_evidence_not_calibration_or_mechanism_evidence": True,
    "fresh_evolved_states_must_reassess_every_saturation_guard": True,
    "successor_runtime_owner": "FGC-1-HLT8-MON8",
}

PROTO10_CLAIMS = {
    "PROTO10_successor_runtime_implemented": False,
    "PROTO10_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO10_resolved_holdout_manifest_authorized": False,
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


def validate_sf1_protocol_v10(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete, fail-closed PROTO10 overlay without I/O."""

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
        raise ValueError("PROTO10 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V10_ARTIFACT_ID,
        "protocol_version": 10,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO10 identity, version, units, purpose, or frozen state differs")

    amendment = _exact("PROTO10 amendment", protocol["amendment"], PROTO10_AMENDMENT)
    replacement = _mapping("PROTO10 replacement", protocol["replacement"])
    if set(replacement) != {"spectral_interpretation", "partition", "provenance"}:
        raise ValueError("PROTO10 replacement keys differ")
    spectral = _exact(
        "PROTO10 spectral interpretation",
        replacement["spectral_interpretation"],
        PROTO10_SPECTRAL_INTERPRETATION,
    )
    partition = _exact(
        "PROTO10 partition", replacement["partition"], PROTO10_PARTITION
    )
    provenance = _exact(
        "PROTO10 provenance", replacement["provenance"], PROTO10_PROVENANCE
    )
    inheritance = _exact(
        "PROTO10 inheritance", protocol["inheritance"], PROTO10_INHERITANCE
    )
    claims = _exact("PROTO10 claims", protocol["claims"], PROTO10_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO10 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "spectral_interpretation": spectral,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V10_ARTIFACT_ID,
        "protocol_version": 10,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO9",
        "diagnosis_artifact_id": "FGC-1-CAL6-PREF8",
        "diagnosis_checkpoint_commit": amendment["diagnosis_checkpoint_commit"],
        "eligible_calibration_amplitudes": ["5/2", "3"],
        "point_counts": list(spectral["point_counts"]),
        "resolution_ladder_unchanged": True,
        "only_finest_pair_spectral_interpretation_changes": True,
        "maximum_nested_tail_ratio": spectral["maximum_nested_tail_ratio"],
        "coarse_to_medium_direct_ratio_veto": True,
        "medium_and_fine_absolute_budgets_unchanged": True,
        "finest_pair_rule": spectral["finest_pair_rule"],
        "finest_pair_direct_ratio_remains_primary": True,
        "guarded_saturation_requires_tail_erasure_map_witness": True,
        "guarded_saturation_requires_round_trip_contraction": True,
        "guarded_saturation_requires_complete_profile_contraction": True,
        "round_trip_is_a_continuum_error_bound": False,
        "diagnostic_saturation_is_physical_resolution_evidence": False,
        "new_numeric_tolerance_added": False,
        "physical_and_numerical_thresholds_unchanged": True,
        "PROTO9_GR0_trajectory_inspected_before_revision": False,
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
    "PROTO10_AMENDMENT",
    "PROTO10_CLAIMS",
    "PROTO10_INHERITANCE",
    "PROTO10_PARTITION",
    "PROTO10_PROVENANCE",
    "PROTO10_SPECTRAL_INTERPRETATION",
    "SF1_PROTOCOL_V10_ARTIFACT_ID",
    "validate_sf1_protocol_v10",
]
