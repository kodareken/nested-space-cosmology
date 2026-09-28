"""Pure validation for the outcome-neutral PROTO11 numerical-map overlay.

PROTO11 consumes immutable PROTO10 and CAL7/PREF9 evidence. It changes only
the floating semidiscrete representation of derivatives around the exact
Minkowski spherical ADM reference. It does not change a continuum equation,
field value, raw source gate, solver, method, retry, physical input, spectrum,
observable, or claim. This module reads no run output and evolves no state.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Mapping


SF1_PROTOCOL_V11_ARTIFACT_ID = "FGC-2-SF1-PROTO11"

PROTO11_AMENDMENT = {
    "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO10",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v10.toml",
    "predecessor_protocol_sha256": "f06ea93e82a4cdda2a0474623bc0db9ac679e4044657afb32be0d5342fe437fa",
    "diagnosis_artifact_id": "FGC-1-CAL7-PREF9",
    "diagnosis_config": "configs/fgc/fgc-1-cal7-pref9.toml",
    "diagnosis_config_sha256": "f7d20efe795a0dd33264ef9cca4ae7977461fa5750eea2a01098b42a0b7c4544",
    "diagnosis_result": "results/fgc-1-cal7-pref9.json",
    "diagnosis_result_sha256": "aab4ea62d569e224bf0bff0f7c3368b2d8dd07007735c1d1993ea3f9ef5b5d2a",
    "diagnosis_checkpoint_commit": "3e34c5e1b0f4193b1178f09725be6e72e462c119",
    "revision_classification": (
        "premise_only_post_calibration_well_balanced_reference_map_revision"
    ),
    "revision_reasons": [
        "PROTO10_completed_without_an_eligible_GR0_case",
        "both_amplitudes_stopped_before_the_first_evolved_common_event",
        "both_stops_were_source_only_on_the_same_RK4_8193_member",
        "CAL7_localized_the_initial_floor_to_the_first_positive_alpha_row",
        "the_raw_floor_scaled_as_binary64_epsilon_over_r_squared",
        "the_floor_was_identical_across_matter_amplitudes",
        "the_well_balanced_t0_control_passed_the_unchanged_raw_source_gate_on_all_twelve_inputs",
        "the_control_is_numerical_premise_evidence_only",
    ],
    "revision_changes": [
        "semidiscrete_reference_state_differentiation_map",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_PROTO10_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_initial_data_and_resolution_ladder_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "source_solver_raw_residual_kinetic_transaction_and_retry_rules_unchanged": True,
    "methods_CFL_dissipation_and_final_time_unchanged": True,
    "constraint_magnitude_order_and_spectral_rules_unchanged": True,
    "direct_one_quarter_ratio_ceiling_and_coarse_veto_unchanged": True,
    "trapped_boundary_candidate_health_and_scale_stops_unchanged": True,
    "eligible_calibration_amplitudes_in_frozen_order": ["5/2", "3"],
    "PROTO1_through_PROTO10_preserved": True,
    "PROTO10_GR0_calibration_trajectory_inspected_before_revision": True,
    "PROTO10_candidate_or_holdout_outcome_observed_before_revision": False,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
    "mechanism_question_answered_before_revision": False,
}

PROTO11_REFERENCE_MAP = {
    "field_order": ["alpha", "shift", "lambda", "R", "phi", "chi"],
    "reference_u": ["1", "0", "1", "r", "0", "0"],
    "reference_q": ["0", "0", "0", "1", "0", "0"],
    "reference_state": "exact_Minkowski_spherical_ADM",
    "reference_is_fixed_and_has_no_fitted_or_evolved_parameters": True,
    "initial_q_rule": "D_h(u-u_ref)+q_ref",
    "source_q_r_rule": "D_h(q-q_ref)",
    "reduction_constraint_rule": "q-(D_h(u-u_ref)+q_ref)",
    "common_event_constraints_use_the_same_rules": True,
    "p_r_rule_is_unchanged": "D_h(p)",
    "du_dt_rule_is_unchanged": "p",
    "dq_dt_rule_is_unchanged": "D_h(p)",
    "q_remains_an_independent_evolved_field": True,
    "interior_q_is_not_reprojected_after_initialization_or_Runge_Kutta_stages": True,
    "reference_subtraction_is_not_added_to_the_continuum_equations": True,
    "reference_subtraction_is_not_damping_or_constraint_cleaning": True,
    "physical_u_p_and_observables_are_not_reference_subtracted": True,
    "raw_source_residual_is_computed_from_the_complete_unredefined_equations": True,
    "raw_source_residual_maximum": "1/1000000000000",
    "no_epsilon_floor_or_new_numeric_tolerance_is_permitted": True,
    "exact_Minkowski_must_be_bitwise_preserved_on_every_frozen_grid": True,
    "nontrivial_generic_REF1_and_reduction_defects_must_remain_visible": True,
    "all_rebuilt_projected_state_hashes_must_be_serialized": True,
}

PROTO11_PARTITION = {
    "resolved_holdout_manifest": "configs/fgc/fgc-1-pro11-hld1.toml",
}

PROTO11_PROVENANCE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto11/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto11/holdout",
    "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    "PROTO10_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO11_outcomes": True,
}

PROTO11_INHERITANCE = {
    "PROTO10_resolution_ladder_inherited_unchanged": True,
    "PROTO10_guarded_spectral_interpretation_inherited_unchanged": True,
    "PROTO10_constraint_magnitude_and_order_rules_inherited_unchanged": True,
    "PROTO10_physical_numerical_transaction_and_stop_rules_inherited_unchanged": True,
    "HLT8_raw_inputs_and_admissions_remain_diagnostic_history_only": True,
    "CAL7_is_numerical_premise_evidence_not_calibration_or_mechanism_evidence": True,
    "fresh_states_must_be_rebuilt_under_the_PROTO11_map": True,
    "fresh_evolved_states_must_reassess_every_source_constraint_spectral_and_health_gate": True,
    "successor_runtime_owner": "FGC-1-HLT9-MON9",
}

PROTO11_CLAIMS = {
    "PROTO11_successor_runtime_implemented": False,
    "PROTO11_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO11_resolved_holdout_manifest_authorized": False,
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


def validate_sf1_protocol_v11(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete, fail-closed PROTO11 overlay without I/O."""

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
        raise ValueError("PROTO11 root keys differ")
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
        "artifact_id": SF1_PROTOCOL_V11_ARTIFACT_ID,
        "protocol_version": 11,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO11 identity, version, units, purpose, or frozen state differs")

    amendment = _exact("PROTO11 amendment", protocol["amendment"], PROTO11_AMENDMENT)
    replacement = _mapping("PROTO11 replacement", protocol["replacement"])
    if set(replacement) != {"reference_state_differentiation", "partition", "provenance"}:
        raise ValueError("PROTO11 replacement keys differ")
    reference_map = _exact(
        "PROTO11 reference-state differentiation",
        replacement["reference_state_differentiation"],
        PROTO11_REFERENCE_MAP,
    )
    partition = _exact(
        "PROTO11 partition", replacement["partition"], PROTO11_PARTITION
    )
    provenance = _exact(
        "PROTO11 provenance", replacement["provenance"], PROTO11_PROVENANCE
    )
    inheritance = _exact(
        "PROTO11 inheritance", protocol["inheritance"], PROTO11_INHERITANCE
    )
    claims = _exact("PROTO11 claims", protocol["claims"], PROTO11_CLAIMS)
    if any(value is not False for value in claims.values()):
        raise ValueError("PROTO11 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "reference_state_differentiation": reference_map,
        "partition": partition,
        "provenance": provenance,
        "inheritance": inheritance,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V11_ARTIFACT_ID,
        "protocol_version": 11,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "predecessor_protocol_artifact_id": "FGC-2-SF1-PROTO10",
        "diagnosis_artifact_id": "FGC-1-CAL7-PREF9",
        "diagnosis_checkpoint_commit": amendment["diagnosis_checkpoint_commit"],
        "eligible_calibration_amplitudes": ["5/2", "3"],
        "point_counts": [2049, 4097, 8193],
        "resolution_ladder_unchanged": True,
        "only_reference_state_differentiation_changes": True,
        "reference_u": list(reference_map["reference_u"]),
        "reference_q": list(reference_map["reference_q"]),
        "initial_q_rule": reference_map["initial_q_rule"],
        "source_q_r_rule": reference_map["source_q_r_rule"],
        "reduction_constraint_rule": reference_map["reduction_constraint_rule"],
        "interior_q_reprojection_permitted": False,
        "raw_source_residual_maximum": reference_map[
            "raw_source_residual_maximum"
        ],
        "new_numeric_tolerance_added": False,
        "physical_and_numerical_thresholds_unchanged": True,
        "PROTO10_GR0_calibration_inspected_before_revision": True,
        "candidate_or_holdout_outcome_observed_before_revision": False,
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
    "PROTO11_AMENDMENT",
    "PROTO11_CLAIMS",
    "PROTO11_INHERITANCE",
    "PROTO11_PARTITION",
    "PROTO11_PROVENANCE",
    "PROTO11_REFERENCE_MAP",
    "SF1_PROTOCOL_V11_ARTIFACT_ID",
    "validate_sf1_protocol_v11",
]
