"""Fail-closed authorization for one classical spherical FGC diagnostic.

RUN1 is deliberately weaker than the retained-EFT authorization in EFT1 and
stronger than a permission to run an arbitrary numerical experiment.  It
composes hash-bound analytic predecessors, one frozen outcome-neutral study
protocol, and future continuum/numerical certificates into an all-of decision.
Missing evidence is false.  No configuration boolean can stand in for a
separate result artifact.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from typing import Any, Mapping, Sequence


RUN1_SYM1_ARTIFACT_ID = "FGC-1-RUN1-SYM1"
SF1_PROTOCOL_V1_ARTIFACT_ID = "FGC-2-SF1-PROTO1"
SF1_PROTOCOL_V2_ARTIFACT_ID = "FGC-2-SF1-PROTO2"
SF1_PROTOCOL_ARTIFACT_ID = "FGC-2-SF1-PROTO3"
SF1_PROTOCOL_V4_ARTIFACT_ID = "FGC-2-SF1-PROTO4"

SF1_PROTOCOL_V2_AMENDMENT = {
    "predecessor_protocol_artifact_id": SF1_PROTOCOL_V1_ARTIFACT_ID,
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol.toml",
    "preflight_artifact_id": "FGC-1-ID0-PREF1",
    "preflight_config": "configs/fgc/fgc-1-id0-pref1.toml",
    "preflight_result": "results/fgc-1-id0-pref1.json",
    "revision_classification": "premise_based_pre_holdout_revision",
    "revision_reason": "PROTO1_omitted_Pi_phi_and_original_r_chi_amplitudes_missed_initial_compactness_floor",
    "minimal_completion_used_for_preflight": "Pi_phi=0",
    "original_protocol_preserved": True,
    "FGCQR_outcomes_inspected_before_revision": False,
    "SGBL_outcomes_inspected_before_revision": False,
    "repair_scan_source": "GR0_constraints_only",
    "repair_scan_candidate_amplitudes": ["2", "5/2", "3", "7/2", "4", "9/2", "5"],
    "preflight_compactness_absolute_upper_bound": "3927454576397784264080597/276701161105643274240000000",
    "protocol_initial_compactness_minimum": "1/10",
}

SF1_PROTOCOL_V3_AMENDMENT = {
    "predecessor_protocol_artifact_id": SF1_PROTOCOL_V2_ARTIFACT_ID,
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v2.toml",
    "diagnosis_artifact_id": "FGC-1-ID1-FAM1",
    "diagnosis_result": "results/fgc-1-id1-fam1.json",
    "diagnosis_checkpoint_commit": "fae2c5b76ba53680f9e70188f79784b160826734",
    "diagnosis_result_sha256": "aeda6543854ca54456dc58648ffc6824be297eb3930f0fbe7c8b8beae5e8075f",
    "predecessor_protocol_sha256": "091d6d7dc8109d21d9297ea6eb37a7296e0ea45fc5e05269bba1b488159f09ba",
    "revision_classification": "premise_based_pre_holdout_revision",
    "revision_reason": "PROTO2_nine_eighths_width_case_missed_declared_center_buffer_by_1/16_L0",
    "original_minimum_center_buffer_over_L0": "5/2",
    "widest_declared_half_width_over_L0": "9/16",
    "widest_case_exact_center_buffer_over_L0": "39/16",
    "amended_minimum_center_buffer_over_L0": "39/16",
    "widest_case_coarsest_grid_intervals": 78,
    "analytic_center_requirement": "strictly_positive_exact_vacuum_neighborhood",
    "static_all_case_initial_premises_required_before_GR0_dynamic_eligibility": True,
    "PROTO1_and_PROTO2_preserved": True,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
}

SF1_PROTOCOL_V4_AMENDMENT = {
    "predecessor_protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "predecessor_protocol_sha256": "57b28eb9abe80c81f9e47ae6bea6ce2d69e103542df500e97c6c93805a8ac13b",
    "diagnosis_artifact_id": "FGC-1-CAL0-PREF2",
    "diagnosis_config": "configs/fgc/fgc-1-cal0-pref2.toml",
    "diagnosis_result": "results/fgc-1-cal0-pref2.json",
    "diagnosis_result_sha256": "d952e737d22d63b33e8d87e366538925bd28acbf82a3062563de49867a229138",
    "diagnosis_checkpoint_commit": "b5aa4934b58afd505b22e008202285c475942892",
    "revision_classification": "premise_based_pre_holdout_numerical_contract_revision",
    "revision_reasons": [
        "PROTO3_last_bin_alias_rule_rejects_declared_compact_pulse_at_t0",
        "PROTO3_conflates_GR0_calibration_with_FGCQR_target_health",
        "PROTO3_constraint_admission_lacks_common_event_error_budget",
    ],
    "revision_changes": [
        "branch_specific_stop_applicability",
        "weighted_spectral_tail_and_RMS_admission",
        "nested_common_event_constraint_admission",
        "versioned_manifest_and_output_namespace",
    ],
    "all_unlisted_contract_fields_inherited_by_exact_predecessor_hash": True,
    "action_couplings_profiles_amplitude_order_and_holdout_cases_unchanged": True,
    "physical_equations_null_observable_outcome_and_robustness_rules_unchanged": True,
    "PROTO1_PROTO2_and_PROTO3_preserved": True,
    "GR0_exploratory_outcomes_inspected_before_revision": True,
    "GR0_exploratory_runs_are_development_only_and_cannot_enter_calibration_result": True,
    "FGCQR_evolution_outcomes_inspected_before_revision": False,
    "SGBL_evolution_outcomes_inspected_before_revision": False,
}

PROTO4_LEGACY_STOP_IDS = (
    "nonpositive_lapse",
    "nonpositive_radial_metric",
    "nonpositive_areal_radius_away_from_center",
    "hat_cone_not_Lorentzian",
    "newton_residual_limit",
    "newton_iteration_limit",
    "newton_residual_not_monotonic",
    "kinetic_condition_limit",
    "branch_displacement_limit",
    "acceleration_limit",
    "all_covector_companion_deformation_limit",
    "action_energy_deformation_limit",
    "kinetic_deformation_limit",
    "effective_planck_ratio_limit",
    "characteristic_imaginary_part_limit",
    "eigenframe_condition_limit",
    "symmetrizer_coercivity_limit",
    "normalized_constraint_limit",
    "constraint_convergence_order_limit",
    "phi_cutoff_proxy_limit",
    "proper_frequency_cutoff_proxy_limit",
    "proper_wavenumber_cutoff_proxy_limit",
    "temporal_alias_band_occupied",
    "radial_alias_band_occupied",
    "curvature_cutoff_proxy_limit",
    "boundary_causal_buffer",
)

PROTO4_UNIVERSAL_STOP_IDS = (
    "nonpositive_lapse",
    "nonpositive_radial_metric",
    "nonpositive_areal_radius_away_from_center",
    "hat_cone_not_Lorentzian",
    "newton_residual_limit",
    "newton_iteration_limit",
    "newton_residual_not_monotonic",
    "kinetic_condition_limit",
    "boundary_causal_buffer",
)

PROTO4_CANDIDATE_BRANCH_STOP_IDS = (
    "branch_displacement_limit",
    "acceleration_limit",
    "all_covector_companion_deformation_limit",
    "action_energy_deformation_limit",
    "kinetic_deformation_limit",
    "effective_planck_ratio_limit",
    "characteristic_imaginary_part_limit",
    "eigenframe_condition_limit",
    "symmetrizer_coercivity_limit",
    "phi_cutoff_proxy_limit",
    "curvature_cutoff_proxy_limit",
)

PROTO4_REPLACED_STOP_IDS = (
    "normalized_constraint_limit",
    "constraint_convergence_order_limit",
    "proper_frequency_cutoff_proxy_limit",
    "proper_wavenumber_cutoff_proxy_limit",
    "temporal_alias_band_occupied",
    "radial_alias_band_occupied",
)

FUTURE_EVIDENCE_REQUIREMENTS: dict[str, tuple[str, str]] = {
    "protocol_holdout": (
        "FGC-1-PRO3-HLD1",
        "resolved_holdout_manifest_frozen_before_FGCQR_outcomes",
    ),
    "run_domain": (
        "FGC-1-DOM4-RUN1",
        "nonzero_classical_spherical_run_envelope_passed",
    ),
    "multidirectional_health": (
        "FGC-1-HYP2-MD1",
        "quantitative_all_covector_weak_coupling_health_envelope_passed",
    ),
    "nonlinear_source": (
        "FGC-1-SRC1-NL1",
        "nonlinear_REF1_source_and_branch_solver_verified",
    ),
    "physical_constraints": (
        "FGC-1-CON4-PHY1",
        "physical_gauge_reduction_constraint_system_closed",
    ),
    "regular_center": (
        "FGC-1-CTR1-REG1",
        "regular_center_formulation_verified",
    ),
    "initial_data": (
        "FGC-1-ID1-FAM1",
        "nonzero_width_finite_mass_constraint_compatible_family_constructed",
    ),
    "boundary_control": (
        "FGC-1-BND2-CP1",
        "spherical_boundary_or_domain_of_dependence_control_passed",
    ),
    "health_monitor": (
        "FGC-1-HLT1-MON1",
        "classical_health_and_typed_stop_monitoring_verified",
    ),
    "numerical_validation": (
        "FGC-1-NUM1-VAL1",
        "independent_solver_validation_and_measurement_contract_passed",
    ),
}

FUTURE_EVIDENCE_CLASSIFICATIONS = {
    "protocol_holdout": "resolved_outcome_neutral_holdout_manifest_certificate",
    "run_domain": "classical_spherical_run_domain_certificate",
    "multidirectional_health": "classical_all_covector_early_kill_certificate",
    "nonlinear_source": "nonlinear_unredefined_source_and_branch_solver_certificate",
    "physical_constraints": "physical_gauge_reduction_constraint_closure_certificate",
    "regular_center": "regular_center_formulation_certificate",
    "initial_data": "finite_mass_constraint_compatible_data_family_certificate",
    "boundary_control": "spherical_domain_of_dependence_and_boundary_certificate",
    "health_monitor": "classical_health_and_typed_stop_monitor_certificate",
    "numerical_validation": "independent_solver_and_affine_measurement_validation_certificate",
}

PROTOCOL_CLAIMS = {
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}

PROTO4_CLAIMS = {
    "classical_spherical_diagnostic_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    **PROTOCOL_CLAIMS,
}

FUTURE_COMMON_NONCLAIMS = dict(PROTOCOL_CLAIMS)
FUTURE_SCOPE_BINDING_KEYS = {
    "run1_artifact_id",
    "run1_config_sha256",
    "protocol_artifact_id",
    "protocol_config_sha256",
    "protocol_semantic_holdout_contract_sha256",
    "action_artifact_id",
    "action_config_sha256",
    "action_result_sha256",
    "variation_artifact_id",
    "variation_config_sha256",
    "variation_result_sha256",
    "target_branch",
    "physical_equations",
    "declared_run_envelope_sha256",
}

FUTURE_PAYLOAD_TRUE_FIELDS: dict[str, tuple[str, ...]] = {
    "run_domain": (
        "nonzero_parameter_volume",
        "FGCQR_target_branch_only",
        "all_protocol_cases_contained",
        "branch_continuity_quantified",
        "classical_scope_only",
    ),
    "multidirectional_health": (
        "all_nonzero_spatial_covectors_covered",
        "strong_hyperbolicity_margin_positive",
        "independent_chi_sector_included",
        "classical_principal_health_only_not_EFT_validity",
    ),
    "nonlinear_source": (
        "unredefined_ACT1_VAR1_equations_used",
        "analytic_residual_and_acceleration_Jacobian_used",
        "exact_fraction_and_floating_controls_agree",
        "finite_difference_Jacobian_crosscheck_passed",
        "branch_continuity_and_typed_failure_verified",
    ),
    "physical_constraints": (
        "Hamiltonian_and_momentum_projections_derived",
        "metric_defined_gauge_constraints_closed",
        "kinematic_reduction_constraints_closed",
        "subsidiary_system_closed",
        "canonical_normalized_monitors_defined",
        "smallness_on_one_run_not_used_as_proof",
    ),
    "regular_center": (
        "R_equals_r_times_A_and_v_equals_r_times_V",
        "all_removable_singular_limits_derived_analytically",
        "negative_power_guard_edge_fails_closed",
        "parity_and_elementary_flatness_enforced",
        "finite_curvature_series_controls_passed",
        "first_grid_point_convergence_passed",
    ),
    "initial_data": (
        "nonzero_width_family",
        "finite_Misner_Sharp_mass",
        "independent_chi_pulse",
        "explicit_nonzero_phi_seed",
        "physical_and_gauge_constraints_solved",
        "exact_outer_vacuum_buffer",
        "continuation_failure_reasons_typed",
    ),
    "boundary_control": (
        "evolving_physical_and_auxiliary_cones_included",
        "SBP_stencil_reach_included",
        "measured_affine_interval_outside_boundary_domain_of_dependence",
        "outer_boundary_move_check_passed",
        "incoming_constraint_characteristics_controlled",
    ),
    "health_monitor": (
        "canonical_dimensionless_norms_executable",
        "all_protocol_thresholds_implemented_strictly",
        "every_accepted_stage_monitored",
        "nonmonotone_internal_stage_abscissae_supported",
        "strict_global_stage_transaction_order_enforced",
        "stop_occurs_before_threshold_crossing",
        "first_failed_premise_is_preserved",
        "scale_stops_not_promoted_to_EFT_validity",
    ),
    "numerical_validation": (
        "Minkowski_preservation_passed",
        "regular_center_manufactured_solution_passed",
        "GR0_scalar_control_passed",
        "exact_vs_floating_backend_passed",
        "constraint_convergence_passed",
        "boundary_isolation_passed",
        "all_typed_stops_injected",
        "checkpoint_restart_equivalence_passed",
        "direct_and_Raychaudhuri_routes_agree",
        "primary_and_comparator_methods_passed",
    ),
}


def _mapping(name: str, value: object) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _boolean(name: str, value: object, expected: bool | None = None) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be boolean")
    if expected is not None and value is not expected:
        raise ValueError(f"{name} must be {str(expected).lower()}")
    return value


def _integer(name: str, value: object, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _fraction(
    name: str,
    value: object,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> Fraction:
    if not isinstance(value, str):
        raise ValueError(f"{name} must be an exact rational string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be an exact rational string") from exc
    if str(result) != value:
        raise ValueError(f"{name} must be a canonical rational string")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _strings(name: str, value: object, *, nonempty: bool = True) -> list[str]:
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise ValueError(f"{name} must be a list of strings")
    if nonempty and not value:
        raise ValueError(f"{name} must not be empty")
    if len(set(value)) != len(value):
        raise ValueError(f"{name} must not contain duplicates")
    return value


def _rational_strings(name: str, value: object) -> list[Fraction]:
    strings = _strings(name, value)
    return [_fraction(f"{name}[{index}]", item) for index, item in enumerate(strings)]


def _semantic_hash(value: Mapping[str, Any]) -> str:
    source = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return sha256(source.encode("utf-8")).hexdigest()


def _sha256_string(name: str, value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _nonempty_string(name: str, value: object) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _validate_sf1_protocol_v1(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the frozen PROTO1 schema and derive its holdout contract hash."""

    protocol = _mapping("protocol", protocol)
    _keys(
        "protocol",
        protocol,
        {
            "schema_version",
            "artifact_id",
            "protocol_version",
            "project_version",
            "frozen",
            "units",
            "purpose",
            "action",
            "dimensionless_parameters",
            "initial_data",
            "calibration",
            "gauge",
            "numerics",
            "health_stops",
            "analysis",
            "norms_and_estimators",
            "partition",
            "provenance",
            "robustness",
            "negative_scope",
            "classification",
            "revision_policy",
            "claims",
        },
    )
    if {
        "schema_version": protocol["schema_version"],
        "artifact_id": protocol["artifact_id"],
        "protocol_version": protocol["protocol_version"],
        "project_version": protocol["project_version"],
        "frozen": protocol["frozen"],
        "units": protocol["units"],
        "purpose": protocol["purpose"],
    } != {
        "schema_version": 1,
        "artifact_id": SF1_PROTOCOL_V1_ARTIFACT_ID,
        "protocol_version": 1,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO1 identity, version, units, purpose, or frozen state differs")

    action = _mapping("protocol.action", protocol["action"])
    _keys(
        "protocol.action",
        action,
        {
            "branch_order",
            "physics_holdout_branch",
            "physical_equations",
            "gauge_extension_role",
            "action_config",
            "variation_config",
            "eft_cutoff_config",
        },
    )
    if action != {
        "branch_order": ["GR-0", "SGB-L", "FGC-QR"],
        "physics_holdout_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "gauge_extension_role": "formulation_only_and_zero_on_the_metric_defined_gauge_surface",
        "action_config": "configs/fgc/fgc-1-action-gate.toml",
        "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
        "eft_cutoff_config": "configs/fgc/fgc-1-eft0-led1.toml",
    }:
        raise ValueError("PROTO1 action contract differs")

    parameters = _mapping("protocol.dimensionless_parameters", protocol["dimensionless_parameters"])
    _keys(
        "protocol.dimensionless_parameters",
        parameters,
        {
            "cutoff_Lambda",
            "M_Pl_over_Lambda",
            "mu_over_Lambda",
            "beta",
            "g4",
            "eta_times_Lambda_squared",
            "alpha_times_Lambda_SGB_L",
            "mu_times_L0",
            "eta_over_L0_squared",
            "fixture_values_are_not_physical_matching",
            "cutoff_is_external_and_not_inferred_from_fixture",
        },
    )
    exact = {
        name: _fraction(f"protocol.dimensionless_parameters.{name}", parameters[name])
        for name in (
            "cutoff_Lambda",
            "M_Pl_over_Lambda",
            "mu_over_Lambda",
            "beta",
            "g4",
            "eta_times_Lambda_squared",
            "alpha_times_Lambda_SGB_L",
            "mu_times_L0",
            "eta_over_L0_squared",
        )
    }
    expected_exact = {
        "cutoff_Lambda": Fraction(16),
        "M_Pl_over_Lambda": Fraction(1, 8),
        "mu_over_Lambda": Fraction(3, 16),
        "beta": Fraction(-1, 4),
        "g4": Fraction(1, 2),
        "eta_times_Lambda_squared": Fraction(128),
        "alpha_times_Lambda_SGB_L": Fraction(-4),
        "mu_times_L0": Fraction(12),
        "eta_over_L0_squared": Fraction(1, 32),
    }
    if exact != expected_exact:
        raise ValueError("PROTO1 dimensionless parameters differ from ACT1/EFT0")
    _boolean(
        "protocol.dimensionless_parameters.fixture_values_are_not_physical_matching",
        parameters["fixture_values_are_not_physical_matching"],
        True,
    )
    _boolean(
        "protocol.dimensionless_parameters.cutoff_is_external_and_not_inferred_from_fixture",
        parameters["cutoff_is_external_and_not_inferred_from_fixture"],
        True,
    )

    initial_data = _mapping("protocol.initial_data", protocol["initial_data"])
    _keys("protocol.initial_data", initial_data, {"chi", "phi_seed"})
    chi = _mapping("protocol.initial_data.chi", initial_data["chi"])
    _keys(
        "protocol.initial_data.chi",
        chi,
        {
            "profile",
            "unit_bump_formula",
            "profile_formula",
            "center_value_rule",
            "amplitude_candidates",
            "amplitude_selection_rule",
            "center_over_L0",
            "half_width_over_L0",
            "momentum_relation",
            "momentum_scale",
            "compactness_definition",
            "initial_compactness_min",
            "initial_compactness_max",
            "minimum_center_buffer_over_L0",
            "minimum_outer_vacuum_buffer_over_L0",
            "require_no_initial_trapped_sphere",
            "require_GR0_trapped_sphere_before_final_time",
            "same_physical_data_across_branches",
            "constraint_solve_may_change_geometry_but_not_declared_chi_profile",
        },
    )
    amplitudes = _rational_strings("protocol.initial_data.chi.amplitude_candidates", chi["amplitude_candidates"])
    if any(value <= 0 for value in amplitudes) or amplitudes != sorted(amplitudes):
        raise ValueError("PROTO1 chi amplitudes must be strictly positive and ordered")
    if len(set(amplitudes)) != len(amplitudes):
        raise ValueError("PROTO1 chi amplitudes must be distinct")
    chi_center = _fraction("protocol.initial_data.chi.center_over_L0", chi["center_over_L0"], positive=True)
    chi_width = _fraction("protocol.initial_data.chi.half_width_over_L0", chi["half_width_over_L0"], positive=True)
    center_buffer = _fraction(
        "protocol.initial_data.chi.minimum_center_buffer_over_L0",
        chi["minimum_center_buffer_over_L0"],
        positive=True,
    )
    if chi_center - chi_width < center_buffer:
        raise ValueError("PROTO1 chi support violates its center buffer")
    compactness_min = _fraction(
        "protocol.initial_data.chi.initial_compactness_min",
        chi["initial_compactness_min"],
        positive=True,
    )
    compactness_max = _fraction(
        "protocol.initial_data.chi.initial_compactness_max",
        chi["initial_compactness_max"],
        positive=True,
    )
    if not compactness_min < compactness_max < 1:
        raise ValueError("PROTO1 initial compactness interval must lie strictly inside (0,1)")
    _fraction("protocol.initial_data.chi.momentum_scale", chi["momentum_scale"], positive=True)
    _fraction(
        "protocol.initial_data.chi.minimum_outer_vacuum_buffer_over_L0",
        chi["minimum_outer_vacuum_buffer_over_L0"],
        positive=True,
    )
    for key in ("require_no_initial_trapped_sphere", "require_GR0_trapped_sphere_before_final_time"):
        _boolean(f"protocol.initial_data.chi.{key}", chi[key], True)
    if chi["profile"] != "compact_C_infinity_bump_of_r_times_chi":
        raise ValueError("PROTO1 chi profile differs")
    if chi["unit_bump_formula"] != (
        "B(x)=exp(1-1/(1-x^2))_for_abs(x)<1_and_0_otherwise"
    ) or chi["profile_formula"] != (
        "r_times_chi(r)=amplitude*B((r-center)/half_width)"
    ) or chi["center_value_rule"] != (
        "chi(0)=0_because_the_declared_support_has_a_strict_center_buffer"
    ):
        raise ValueError("PROTO1 chi bump formula differs")
    if chi["amplitude_selection_rule"] != (
        "smallest_GR0_candidate_with_no_initial_trapped_sphere_and_a_converged_trapped_sphere_before_final_time"
    ):
        raise ValueError("PROTO1 chi calibration rule differs")
    if chi["momentum_relation"] != (
        "partial_t(r_chi)=partial_r(r_chi)_for_the_declared_future_ingoing_orientation"
    ):
        raise ValueError("PROTO1 ingoing momentum relation differs")
    if chi["compactness_definition"] != "max_r(2*m_Misner_Sharp/R)":
        raise ValueError("PROTO1 compactness definition differs")
    if chi["same_physical_data_across_branches"] != (
        "same_areal_radius_chi_profile_and_future_unit_normal_derivative_before_each_branch_constraint_solve"
    ):
        raise ValueError("PROTO1 cross-branch physical-data map differs")
    _boolean(
        "protocol.initial_data.chi.constraint_solve_may_change_geometry_but_not_declared_chi_profile",
        chi["constraint_solve_may_change_geometry_but_not_declared_chi_profile"],
        True,
    )

    phi_seed = _mapping("protocol.initial_data.phi_seed", initial_data["phi_seed"])
    _keys(
        "protocol.initial_data.phi_seed",
        phi_seed,
        {
            "profile",
            "unit_bump_formula",
            "profile_formula",
            "amplitude",
            "center_over_L0",
            "half_width_over_L0",
            "nonzero_seed_required",
            "roundoff_activation_forbidden",
        },
    )
    phi_amplitude = _fraction(
        "protocol.initial_data.phi_seed.amplitude", phi_seed["amplitude"], positive=True
    )
    phi_center = _fraction(
        "protocol.initial_data.phi_seed.center_over_L0", phi_seed["center_over_L0"], positive=True
    )
    phi_width = _fraction(
        "protocol.initial_data.phi_seed.half_width_over_L0", phi_seed["half_width_over_L0"], positive=True
    )
    if phi_seed["profile"] != "compact_C_infinity_bump" or phi_center - phi_width <= 0:
        raise ValueError("PROTO1 phi seed profile or support differs")
    if phi_seed["unit_bump_formula"] != (
        "B(x)=exp(1-1/(1-x^2))_for_abs(x)<1_and_0_otherwise"
    ) or phi_seed["profile_formula"] != (
        "phi(r)=amplitude*B((r-center)/half_width)"
    ):
        raise ValueError("PROTO1 phi seed formula differs")
    for key in ("nonzero_seed_required", "roundoff_activation_forbidden"):
        _boolean(f"protocol.initial_data.phi_seed.{key}", phi_seed[key], True)

    calibration_contract = _mapping("protocol.calibration", protocol["calibration"])
    _keys(
        "protocol.calibration",
        calibration_contract,
        {
            "branch",
            "candidate_order",
            "uses_only_GR0_outcomes",
            "primary_method",
            "selected_case_requires_comparator_confirmation",
            "required_nested_resolutions",
            "trapped_detection",
            "minimum_consecutive_trapped_samples",
            "all_RUN1_health_and_constraint_stops_apply",
            "tie_break",
            "no_eligible_outcome",
            "no_eligible_outcome_forbids_FGCQR_execution",
            "selected_amplitude_is_not_a_fit_to_FGCQR_outcomes",
            "manifest_must_serialize_every_candidate_config_result_hash_and_eligibility_reason",
        },
    )
    expected_calibration_strings = {
        "branch": "GR-0",
        "candidate_order": "strictly_ascending_as_declared_in_initial_data.chi.amplitude_candidates",
        "primary_method": "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure_plus_RK4",
        "trapped_detection": "theta_plus_and_theta_minus_each_strictly_negative_by_more_than_four_times_their_combined_error_bound_at_eight_consecutive_common_events",
        "tie_break": "first_eligible_amplitude_in_declared_order",
        "no_eligible_outcome": "calibration_failed_no_eligible_GR0_case",
    }
    if any(
        calibration_contract[key] != value
        for key, value in expected_calibration_strings.items()
    ):
        raise ValueError("PROTO1 calibration algorithm differs")
    if _integer(
        "protocol.calibration.required_nested_resolutions",
        calibration_contract["required_nested_resolutions"],
        minimum=3,
    ) != 3 or _integer(
        "protocol.calibration.minimum_consecutive_trapped_samples",
        calibration_contract["minimum_consecutive_trapped_samples"],
        minimum=2,
    ) != 8:
        raise ValueError("PROTO1 calibration resolution or persistence rule differs")
    for key in (
        "uses_only_GR0_outcomes",
        "selected_case_requires_comparator_confirmation",
        "all_RUN1_health_and_constraint_stops_apply",
        "no_eligible_outcome_forbids_FGCQR_execution",
        "selected_amplitude_is_not_a_fit_to_FGCQR_outcomes",
        "manifest_must_serialize_every_candidate_config_result_hash_and_eligibility_reason",
    ):
        _boolean(f"protocol.calibration.{key}", calibration_contract[key], True)

    gauge = _mapping("protocol.gauge", protocol["gauge"])
    _keys(
        "protocol.gauge",
        gauge,
        {
            "formulation",
            "tilde_normal_factor",
            "hat_normal_factor",
            "normal_factor_order",
            "reference_policy",
            "physical_metric_defines_null_observables",
            "auxiliary_cones_are_not_physical_light_cones",
        },
    )
    if gauge["formulation"] != "metric_defined_reference_modified_harmonic":
        raise ValueError("PROTO1 gauge formulation differs")
    if (
        _integer("protocol.gauge.tilde_normal_factor", gauge["tilde_normal_factor"], minimum=2),
        _integer("protocol.gauge.hat_normal_factor", gauge["hat_normal_factor"], minimum=2),
        gauge["normal_factor_order"],
    ) != (4, 9, "1<tilde<hat"):
        raise ValueError("PROTO1 auxiliary cone choice differs")
    if gauge["reference_policy"] != (
        "derive_a_regular_center_reference_or_rederive_every_affected_REF1_descendant"
    ):
        raise ValueError("PROTO1 regular-center reference policy differs")
    for key in ("physical_metric_defines_null_observables", "auxiliary_cones_are_not_physical_light_cones"):
        _boolean(f"protocol.gauge.{key}", gauge[key], True)

    numerics = _mapping("protocol.numerics", protocol["numerics"])
    _keys(
        "protocol.numerics",
        numerics,
        {
            "radial_domain_min_over_L0",
            "radial_domain_max_over_L0",
            "measurement_radius_max_over_L0",
            "final_time_over_L0",
            "resolutions",
            "cfl",
            "ko_dissipation",
            "checkpoint_interval_over_L0",
            "output_interval_over_L0",
            "primary_spatial_method",
            "primary_time_method",
            "comparator_spatial_method",
            "comparator_time_method",
            "no_excision_in_measurement_region",
            "monitor_every_accepted_stage",
            "numpy_is_only_planned_numerical_dependency",
        },
    )
    radial_min = _fraction(
        "protocol.numerics.radial_domain_min_over_L0",
        numerics["radial_domain_min_over_L0"],
        nonnegative=True,
    )
    radial_max = _fraction(
        "protocol.numerics.radial_domain_max_over_L0",
        numerics["radial_domain_max_over_L0"],
        positive=True,
    )
    measurement_max = _fraction(
        "protocol.numerics.measurement_radius_max_over_L0",
        numerics["measurement_radius_max_over_L0"],
        positive=True,
    )
    final_time = _fraction(
        "protocol.numerics.final_time_over_L0", numerics["final_time_over_L0"], positive=True
    )
    if radial_min != 0 or not measurement_max < radial_max or final_time <= 0:
        raise ValueError("PROTO1 numerical domain differs or is empty")
    resolutions = numerics["resolutions"]
    if not isinstance(resolutions, list) or any(
        isinstance(value, bool) or not isinstance(value, int) for value in resolutions
    ):
        raise ValueError("PROTO1 resolutions must be integer grid sizes")
    if len(resolutions) != 3 or any(
        resolutions[index + 1] != 2 * resolutions[index] - 1
        for index in range(len(resolutions) - 1)
    ):
        raise ValueError("PROTO1 resolutions must be three exactly nested grids")
    cfl = _fraction("protocol.numerics.cfl", numerics["cfl"], positive=True)
    if cfl >= 1:
        raise ValueError("PROTO1 CFL must be strictly below one")
    for key in ("ko_dissipation", "checkpoint_interval_over_L0", "output_interval_over_L0"):
        _fraction(f"protocol.numerics.{key}", numerics[key], positive=True)
    expected_methods = {
        "primary_spatial_method": "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure",
        "primary_time_method": "RK4",
        "comparator_spatial_method": "second_order_diagonal_norm_SBP",
        "comparator_time_method": "SSPRK3",
    }
    if any(numerics[key] != value for key, value in expected_methods.items()):
        raise ValueError("PROTO1 numerical method contract differs")
    for key in (
        "no_excision_in_measurement_region",
        "monitor_every_accepted_stage",
        "numpy_is_only_planned_numerical_dependency",
    ):
        _boolean(f"protocol.numerics.{key}", numerics[key], True)

    health = _mapping("protocol.health_stops", protocol["health_stops"])
    health_rationals = {
        "newton_residual_infinity_max",
        "kinetic_condition_number_max",
        "normalized_branch_displacement_max",
        "effective_planck_ratio_min",
        "characteristic_imaginary_part_max",
        "eigenframe_condition_number_max",
        "symmetrizer_coercivity_min",
        "normalized_constraint_infinity_max",
        "minimum_observed_constraint_convergence_order",
        "phi_over_Lambda_max",
        "proper_frequency_over_Lambda_max",
        "proper_wavenumber_over_Lambda_max",
        "curvature_scale_over_Lambda_squared_max",
        "minimum_boundary_causal_buffer_over_L0",
    }
    _keys(
        "protocol.health_stops",
        health,
        health_rationals
        | {
            "newton_iteration_max",
            "newton_residual_must_decrease_monotonically",
            "boundary_isolation_uses_evolving_all_cone_speed_envelope",
            "boundary_isolation_includes_SBP_stencil_reach",
            "outer_boundary_move_check_required",
            "stop_before_threshold_crossing",
            "first_failed_premise_determines_stop",
            "declared_scale_stops_are_not_a_complete_EFT_remainder_bound",
        },
    )
    for key in health_rationals:
        _fraction(f"protocol.health_stops.{key}", health[key], positive=True)
    _integer("protocol.health_stops.newton_iteration_max", health["newton_iteration_max"], minimum=1)
    for key in (
        "newton_residual_must_decrease_monotonically",
        "boundary_isolation_uses_evolving_all_cone_speed_envelope",
        "boundary_isolation_includes_SBP_stencil_reach",
        "outer_boundary_move_check_required",
        "stop_before_threshold_crossing",
        "first_failed_premise_determines_stop",
        "declared_scale_stops_are_not_a_complete_EFT_remainder_bound",
    ):
        _boolean(f"protocol.health_stops.{key}", health[key], True)

    analysis = _mapping("protocol.analysis", protocol["analysis"])
    _keys(
        "protocol.analysis",
        analysis,
        {
            "future_null_orientation",
            "null_generator",
            "affine_launch_slice",
            "affine_generator_label",
            "affine_initial_normalization",
            "affine_origin",
            "affine_transport_equation",
            "affine_normalization_residual_max",
            "cross_resolution_trajectory_alignment",
            "expansion_definition",
            "trapped_definition",
            "raychaudhuri_rhs",
            "activation_observable",
            "minimum_resolved_activation_factor",
            "activation_requires_growth_over_seed_and_matched_controls",
            "matched_GR0_and_SGBL_controls_required",
            "error_components",
            "error_combination",
            "direct_Raychaudhuri_agreement_residual_max",
            "positive_margin_over_error_factor",
            "minimum_consecutive_positive_samples",
            "minimum_affine_interval_over_L0",
            "minimum_nested_resolutions",
            "minimum_independent_methods",
            "minimum_observed_solution_convergence_order",
            "direct_and_Raychaudhuri_routes_must_agree",
            "Ricci_term_alone_is_not_defocusing",
        },
    )
    if analysis["future_null_orientation"] != (
        "outgoing_k_plus_has_positive_areal_derivative_in_the_Minkowski_control"
    ):
        raise ValueError("PROTO1 null orientation differs")
    if analysis["null_generator"] != "affinely_normalized_physical_metric_radial_null_generator":
        raise ValueError("PROTO1 affine generator differs")
    expected_affine = {
        "affine_launch_slice": "initial_constraint_solved_slice_t=0",
        "affine_generator_label": "initial_areal_radius_R0",
        "affine_initial_normalization": "minus_g_ab_k_plus^a_n_slice^b=1_at_the_launch_slice",
        "affine_origin": "lambda_plus=0_at_the_launch_slice",
        "affine_transport_equation": "k_plus^b*nabla_b*k_plus^a=0",
        "cross_resolution_trajectory_alignment": "same_initial_areal_radius_R0_and_same_affine_origin",
    }
    if any(analysis[key] != value for key, value in expected_affine.items()):
        raise ValueError("PROTO1 affine normalization or alignment differs")
    _fraction(
        "protocol.analysis.affine_normalization_residual_max",
        analysis["affine_normalization_residual_max"],
        positive=True,
    )
    if analysis["expansion_definition"] != "theta_plus=(2/R)*dR/dlambda_plus":
        raise ValueError("PROTO1 expansion definition differs")
    if analysis["trapped_definition"] != "theta_plus<0_and_theta_minus<0_with_one_fixed_orientation":
        raise ValueError("PROTO1 trapped-sphere definition differs")
    if analysis["raychaudhuri_rhs"] != "-theta_squared/2-shear_squared-R_ab_k^a_k^b":
        raise ValueError("PROTO1 Raychaudhuri definition differs")
    if analysis["activation_observable"] != (
        "max_measurement_domain_abs_phi_divided_by_initial_seed_peak"
    ):
        raise ValueError("PROTO1 regulator activation observable differs")
    if _fraction(
        "protocol.analysis.minimum_resolved_activation_factor",
        analysis["minimum_resolved_activation_factor"],
        positive=True,
    ) < 2:
        raise ValueError("PROTO1 regulator activation factor is too weak")
    for key in (
        "activation_requires_growth_over_seed_and_matched_controls",
        "matched_GR0_and_SGBL_controls_required",
    ):
        _boolean(f"protocol.analysis.{key}", analysis[key], True)
    if _strings("protocol.analysis.error_components", analysis["error_components"]) != [
        "spatial_temporal_Richardson",
        "constraint",
        "initial_data",
        "nonlinear_source",
        "gauge",
        "affine",
        "null",
        "trajectory_alignment",
        "interpolation",
        "extraction_window",
        "boundary",
    ]:
        raise ValueError("PROTO1 error budget differs")
    if analysis["error_combination"] != "conservative_sum_of_componentwise_upper_bounds":
        raise ValueError("PROTO1 error combination must remain conservative")
    _fraction(
        "protocol.analysis.direct_Raychaudhuri_agreement_residual_max",
        analysis["direct_Raychaudhuri_agreement_residual_max"],
        positive=True,
    )
    if _integer(
        "protocol.analysis.positive_margin_over_error_factor",
        analysis["positive_margin_over_error_factor"],
        minimum=1,
    ) != 4:
        raise ValueError("PROTO1 positive margin factor must remain four")
    _integer(
        "protocol.analysis.minimum_consecutive_positive_samples",
        analysis["minimum_consecutive_positive_samples"],
        minimum=2,
    )
    _fraction(
        "protocol.analysis.minimum_affine_interval_over_L0",
        analysis["minimum_affine_interval_over_L0"],
        positive=True,
    )
    if _integer(
        "protocol.analysis.minimum_nested_resolutions",
        analysis["minimum_nested_resolutions"],
        minimum=1,
    ) < 3:
        raise ValueError("PROTO1 requires at least three nested resolutions")
    if _integer(
        "protocol.analysis.minimum_independent_methods",
        analysis["minimum_independent_methods"],
        minimum=1,
    ) < 2:
        raise ValueError("PROTO1 requires at least two numerical methods")
    _fraction(
        "protocol.analysis.minimum_observed_solution_convergence_order",
        analysis["minimum_observed_solution_convergence_order"],
        positive=True,
    )
    for key in ("direct_and_Raychaudhuri_routes_must_agree", "Ricci_term_alone_is_not_defocusing"):
        _boolean(f"protocol.analysis.{key}", analysis[key], True)

    norms = _mapping("protocol.norms_and_estimators", protocol["norms_and_estimators"])
    _keys(
        "protocol.norms_and_estimators",
        norms,
        {
            "state_scaling",
            "matrix_norm",
            "branch_displacement_reference",
            "constraint_normalization",
            "zero_denominator_policy",
            "curvature_scale",
            "proper_frequency",
            "proper_wavenumber",
            "window_and_alias_policy",
            "definition_owner",
            "thresholds_are_classical_stop_proxies_not_complete_EFT_control",
        },
    )
    expected_norms = {
        "state_scaling": "dimensionless_U_uses_L0_and_canonical_mass_dimensions_for_u_p_q",
        "matrix_norm": "induced_infinity_norm_after_declared_dimensionless_row_and_column_scaling",
        "branch_displacement_reference": "continuation_root_from_the_preceding_accepted_Runge_Kutta_stage",
        "constraint_normalization": "componentwise_max_of_declared_physical_source_scale_and_1e-30_before_infinity_norm",
        "zero_denominator_policy": "zero_numerator_and_denominator_maps_to_zero_otherwise_fail_closed",
        "curvature_scale": "maximum_measurement_domain_fourth_root_of_abs_Kretschmann_and_square_root_of_abs_Ricci_eigenvalues",
        "proper_frequency": "maximum_over_declared_fields_of_windowed_local_orthonormal_normal_derivative_spectral_support",
        "proper_wavenumber": "maximum_over_declared_fields_of_windowed_proper_radial_derivative_spectral_support",
        "window_and_alias_policy": "compact_vacuum_buffer_window_with_top_one_eighth_of_resolved_spectrum_treated_as_unresolved_and_stopped",
        "definition_owner": "FGC-1-HLT1-MON1_must_serialize_executable_formulas_before_holdout",
    }
    if any(norms[key] != value for key, value in expected_norms.items()):
        raise ValueError("PROTO1 norm or invariant-estimator contract differs")
    _boolean(
        "protocol.norms_and_estimators.thresholds_are_classical_stop_proxies_not_complete_EFT_control",
        norms["thresholds_are_classical_stop_proxies_not_complete_EFT_control"],
        True,
    )

    partition = _mapping("protocol.partition", protocol["partition"])
    _keys(
        "protocol.partition",
        partition,
        {
            "development_control_ids",
            "calibration_case_ids",
            "held_out_case_ids",
            "held_out_phi_seed_factors",
            "held_out_chi_width_factors",
            "resolved_holdout_manifest",
            "resolved_holdout_manifest_required_before_authorization",
            "held_out_outcomes_inspected_before_freeze",
        },
    )
    development = _strings("protocol.partition.development_control_ids", partition["development_control_ids"])
    calibration_cases = _strings(
        "protocol.partition.calibration_case_ids", partition["calibration_case_ids"]
    )
    holdout = _strings("protocol.partition.held_out_case_ids", partition["held_out_case_ids"])
    if (
        set(development) & set(calibration_cases)
        or set(development) & set(holdout)
        or set(calibration_cases) & set(holdout)
    ):
        raise ValueError("PROTO1 development, calibration, and holdout partitions must be disjoint")
    for key in ("held_out_phi_seed_factors", "held_out_chi_width_factors"):
        factors = _rational_strings(f"protocol.partition.{key}", partition[key])
        if any(value <= 0 for value in factors) or Fraction(1) not in factors:
            raise ValueError(f"PROTO1 {key} must contain positive factors including one")
    if partition["resolved_holdout_manifest"] != "configs/fgc/fgc-1-pro2-hld1.toml":
        raise ValueError("PROTO1 resolved holdout manifest path differs")
    _boolean(
        "protocol.partition.resolved_holdout_manifest_required_before_authorization",
        partition["resolved_holdout_manifest_required_before_authorization"],
        True,
    )
    _boolean(
        "protocol.partition.held_out_outcomes_inspected_before_freeze",
        partition["held_out_outcomes_inspected_before_freeze"],
        False,
    )

    provenance = _mapping("protocol.provenance", protocol["provenance"])
    _keys(
        "protocol.provenance",
        provenance,
        {
            "machine_claim_scope",
            "human_private_observation_not_machine_verifiable",
            "holdout_output_root",
            "manifest_requires_preexisting_output_inventory_to_be_empty",
            "runner_refuses_existing_output_or_overwrite",
            "manifest_requires_ordered_expanded_case_config_hashes",
            "manifest_requires_append_only_event_log_hash",
            "every_holdout_result_embeds_manifest_result_hash",
            "all_scope_bindings_precede_FGCQR_execution",
        },
    )
    if provenance["machine_claim_scope"] != (
        "repository_level_sequencing_not_private_human_blindness"
    ) or provenance["holdout_output_root"] != "runs/fgc-2-sf1/holdout":
        raise ValueError("PROTO1 provenance scope or output root differs")
    for key in (
        "human_private_observation_not_machine_verifiable",
        "manifest_requires_preexisting_output_inventory_to_be_empty",
        "runner_refuses_existing_output_or_overwrite",
        "manifest_requires_ordered_expanded_case_config_hashes",
        "manifest_requires_append_only_event_log_hash",
        "every_holdout_result_embeds_manifest_result_hash",
        "all_scope_bindings_precede_FGCQR_execution",
    ):
        _boolean(f"protocol.provenance.{key}", provenance[key], True)

    robustness = _mapping("protocol.robustness", protocol["robustness"])
    _keys(
        "protocol.robustness",
        robustness,
        {
            "physical_box_design",
            "finite_samples_alone_do_not_prove_an_open_neighborhood",
            "open_neighborhood_requires_interval_or_continuity_certificate",
            "all_declared_physical_box_axes_have_nonzero_half_width",
            "ROB1_separates_physical_box_from_numerical_gauge_diagnostics",
            "gauge_cone_pairs",
            "cfl_values",
            "ko_dissipation_values",
            "outer_boundary_radius_over_L0_values",
            "relative_half_widths",
        },
    )
    if robustness["physical_box_design"] != (
        "closed_rational_relative_box_around_the_resolved_central_case"
    ):
        raise ValueError("PROTO1 ROB1 physical box differs")
    for key in (
        "finite_samples_alone_do_not_prove_an_open_neighborhood",
        "open_neighborhood_requires_interval_or_continuity_certificate",
        "all_declared_physical_box_axes_have_nonzero_half_width",
        "ROB1_separates_physical_box_from_numerical_gauge_diagnostics",
    ):
        _boolean(f"protocol.robustness.{key}", robustness[key], True)
    if _strings("protocol.robustness.gauge_cone_pairs", robustness["gauge_cone_pairs"]) != [
        "3:7",
        "4:9",
        "5:11",
    ]:
        raise ValueError("PROTO1 gauge robustness cases differ")
    expected_robustness_values = {
        "cfl_values": [Fraction(1, 16), Fraction(3, 32), Fraction(1, 8)],
        "ko_dissipation_values": [Fraction(1, 128), Fraction(1, 64), Fraction(3, 128)],
        "outer_boundary_radius_over_L0_values": [Fraction(24), Fraction(32), Fraction(40)],
    }
    for key, expected in expected_robustness_values.items():
        if _rational_strings(f"protocol.robustness.{key}", robustness[key]) != expected:
            raise ValueError(f"PROTO1 {key} robustness cases differ")
    half_widths = _mapping(
        "protocol.robustness.relative_half_widths", robustness["relative_half_widths"]
    )
    _keys(
        "protocol.robustness.relative_half_widths",
        half_widths,
        {
            "mu",
            "beta",
            "g4",
            "eta",
            "resolved_chi_amplitude",
            "chi_half_width",
            "phi_seed_amplitude",
        },
    )
    if any(
        _fraction(
            f"protocol.robustness.relative_half_widths.{key}", value, positive=True
        ) != Fraction(1, 32)
        for key, value in half_widths.items()
    ):
        raise ValueError("PROTO1 physical robustness half-widths differ")

    negative = _mapping("protocol.negative_scope", protocol["negative_scope"])
    _keys(
        "protocol.negative_scope",
        negative,
        {
            "completed_no_defocusing_scope",
            "branch_rejection_scope",
            "single_stop_or_solver_failure_cannot_reject_branch",
            "general_gradient_phase_transition_variable_cone_or_nested_domain_rejection_forbidden",
            "obstruction_requires_invariant_or_canonical_dimensionless_quantity",
            "obstruction_location_and_threshold_must_converge_under_refinement",
            "obstruction_requires_both_numerical_methods",
            "obstruction_requires_gauge_domain_and_boundary_perturbation",
            "obstruction_requires_analytic_inequality_or_independent_solver_formulation",
            "obstruction_certificate_required_for_branch_rejected_label",
        },
    )
    if negative["completed_no_defocusing_scope"] != (
        "only_the_frozen_finite_holdout_cases_time_window_and_measurement_domain"
    ) or negative["branch_rejection_scope"] != (
        "explicit_parameter_box_times_data_family_times_formulation_named_by_the_obstruction_certificate"
    ):
        raise ValueError("PROTO1 negative-result scope differs")
    for key, value in negative.items():
        if key not in {"completed_no_defocusing_scope", "branch_rejection_scope"}:
            _boolean(f"protocol.negative_scope.{key}", value, True)

    classification = _mapping("protocol.classification", protocol["classification"])
    _keys(
        "protocol.classification",
        classification,
        {
            "positive_label",
            "negative_labels",
            "invalid_label",
            "typed_stop_labels",
            "calibration_failure_label",
            "measurement_stop_label",
            "single_stop_is_only_an_obstruction_candidate",
            "scientific_negative_requires_family_completion_or_robust_converged_obstruction",
        },
    )
    if classification["positive_label"] != "completed_defocusing_candidate":
        raise ValueError("PROTO1 positive outcome label differs")
    if _strings("protocol.classification.negative_labels", classification["negative_labels"]) != [
        "completed_no_defocusing",
        "branch_rejected_by_converged_health_or_constraint_obstruction",
    ]:
        raise ValueError("PROTO1 negative outcome labels differ")
    if classification["invalid_label"] != "invalid_implementation_or_nonconverged_run":
        raise ValueError("PROTO1 invalid outcome label differs")
    if _strings("protocol.classification.typed_stop_labels", classification["typed_stop_labels"]) != [
        "stopped_branch_loss",
        "stopped_hyperbolicity_loss",
        "stopped_constraint_loss",
        "stopped_scale_control",
        "stopped_boundary_contamination",
        "stopped_affine_or_measurement_contract_loss",
    ]:
        raise ValueError("PROTO1 typed stop labels differ")
    if classification["calibration_failure_label"] != (
        "calibration_failed_no_eligible_GR0_case"
    ) or classification["measurement_stop_label"] != (
        "stopped_affine_or_measurement_contract_loss"
    ):
        raise ValueError("PROTO1 calibration or measurement stop label differs")
    for key in (
        "single_stop_is_only_an_obstruction_candidate",
        "scientific_negative_requires_family_completion_or_robust_converged_obstruction",
    ):
        _boolean(f"protocol.classification.{key}", classification[key], True)

    revision = _mapping("protocol.revision_policy", protocol["revision_policy"])
    _keys(
        "protocol.revision_policy",
        revision,
        {
            "outcome_based_parameter_changes_forbidden",
            "premise_based_revision_requires_new_protocol_version",
            "post_unblinding_revision_reclassifies_affected_cases_as_exploratory",
            "implementation_fault_rerun_requires_identical_configuration_hash",
            "all_resolved_holdout_hashes_must_precede_FGCQR_execution",
        },
    )
    for key, value in revision.items():
        _boolean(f"protocol.revision_policy.{key}", value, True)

    claims = _mapping("protocol.claims", protocol["claims"])
    if dict(claims) != PROTOCOL_CLAIMS:
        raise ValueError("PROTO1 claims must remain fail-closed")

    holdout_contract = {
        "action": action,
        "dimensionless_parameters": parameters,
        "initial_data": initial_data,
        "calibration": calibration_contract,
        "gauge": gauge,
        "numerics": numerics,
        "health_stops": health,
        "analysis": analysis,
        "norms_and_estimators": norms,
        "partition": partition,
        "provenance": provenance,
        "robustness": robustness,
        "negative_scope": negative,
        "classification": classification,
        "revision_policy": revision,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V1_ARTIFACT_ID,
        "protocol_version": 1,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "held_out_outcomes_inspected_before_freeze": False,
        "explicit_nonzero_phi_seed": str(phi_amplitude),
        "roundoff_activation_forbidden": True,
        "calibration_case_count": len(calibration_cases),
        "held_out_case_count": len(holdout),
        "amplitude_candidates": [str(value) for value in amplitudes],
        "calibration_case_ids": calibration_cases,
        "held_out_case_ids": holdout,
        "resolved_holdout_manifest_required": True,
        "resolved_holdout_manifest_present": False,
        "semantic_holdout_contract_sha256": _semantic_hash(holdout_contract),
    }


def _validate_sf1_protocol_v2(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the premise-corrected PROTO2 contract.

    PROTO2 changes only the pre-holdout initial-data premise and the versioned
    output namespace.  Every other physics, numerical, health, affine,
    robustness, outcome, and nonclaim rule is validated by reducing an exact
    copy to the already sealed PROTO1 schema.  The amendment itself remains in
    the PROTO2 semantic hash, so a premise or provenance change invalidates all
    downstream scope bindings.
    """

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
        "action",
        "dimensionless_parameters",
        "initial_data",
        "calibration",
        "gauge",
        "numerics",
        "health_stops",
        "analysis",
        "norms_and_estimators",
        "partition",
        "provenance",
        "robustness",
        "negative_scope",
        "classification",
        "revision_policy",
        "claims",
    }
    _keys("protocol", protocol, expected_root)
    if {
        "schema_version": protocol["schema_version"],
        "artifact_id": protocol["artifact_id"],
        "protocol_version": protocol["protocol_version"],
        "project_version": protocol["project_version"],
        "frozen": protocol["frozen"],
        "units": protocol["units"],
        "purpose": protocol["purpose"],
    } != {
        "schema_version": 1,
        "artifact_id": SF1_PROTOCOL_V2_ARTIFACT_ID,
        "protocol_version": 2,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError("PROTO2 identity, version, units, purpose, or frozen state differs")

    amendment = _mapping("protocol.amendment", protocol["amendment"])
    expected_amendment = SF1_PROTOCOL_V2_AMENDMENT
    _keys("protocol.amendment", amendment, set(expected_amendment))
    if dict(amendment) != expected_amendment:
        raise ValueError("PROTO2 pre-holdout amendment or provenance differs")

    initial_data = _mapping("protocol.initial_data", protocol["initial_data"])
    chi = _mapping("protocol.initial_data.chi", initial_data.get("chi"))
    amplitudes = _rational_strings(
        "protocol.initial_data.chi.amplitude_candidates",
        chi.get("amplitude_candidates"),
    )
    expected_amplitudes = [
        Fraction(2),
        Fraction(5, 2),
        Fraction(3),
        Fraction(7, 2),
        Fraction(4),
        Fraction(9, 2),
        Fraction(5),
    ]
    if amplitudes != expected_amplitudes:
        raise ValueError("PROTO2 chi amplitudes differ from the pre-holdout GR0 repair scan")
    phi_seed = _mapping("protocol.initial_data.phi_seed", initial_data.get("phi_seed"))
    expected_phi_keys = {
        "profile",
        "unit_bump_formula",
        "profile_formula",
        "amplitude",
        "center_over_L0",
        "half_width_over_L0",
        "unit_normal_momentum",
        "unit_normal_momentum_frame",
        "nonzero_seed_required",
        "roundoff_activation_forbidden",
    }
    _keys("protocol.initial_data.phi_seed", phi_seed, expected_phi_keys)
    if _fraction(
        "protocol.initial_data.phi_seed.unit_normal_momentum",
        phi_seed["unit_normal_momentum"],
        nonnegative=True,
    ) != 0 or phi_seed["unit_normal_momentum_frame"] != "future_slice_unit_normal":
        raise ValueError("PROTO2 must freeze Pi_phi=0 in the future slice unit-normal frame")
    provenance = _mapping("protocol.provenance", protocol["provenance"])
    if provenance.get("holdout_output_root") != "runs/fgc-2-sf1/proto2/holdout":
        raise ValueError("PROTO2 holdout output namespace differs")

    normalized = deepcopy(dict(protocol))
    normalized.pop("amendment")
    normalized["artifact_id"] = SF1_PROTOCOL_V1_ARTIFACT_ID
    normalized["protocol_version"] = 1
    normalized["initial_data"]["chi"]["amplitude_candidates"] = [
        "1/64",
        "3/128",
        "1/32",
        "3/64",
        "1/16",
        "3/32",
        "1/8",
    ]
    normalized["initial_data"]["phi_seed"].pop("unit_normal_momentum")
    normalized["initial_data"]["phi_seed"].pop("unit_normal_momentum_frame")
    normalized["provenance"]["holdout_output_root"] = "runs/fgc-2-sf1/holdout"
    inherited = _validate_sf1_protocol_v1(normalized)

    holdout_contract = {
        key: protocol[key]
        for key in (
            "amendment",
            "action",
            "dimensionless_parameters",
            "initial_data",
            "calibration",
            "gauge",
            "numerics",
            "health_stops",
            "analysis",
            "norms_and_estimators",
            "partition",
            "provenance",
            "robustness",
            "negative_scope",
            "classification",
            "revision_policy",
        )
    }
    return {
        **inherited,
        "artifact_id": SF1_PROTOCOL_V2_ARTIFACT_ID,
        "protocol_version": 2,
        "amendment_validated": True,
        "predecessor_protocol_artifact_id": SF1_PROTOCOL_V1_ARTIFACT_ID,
        "preflight_artifact_id": "FGC-1-ID0-PREF1",
        "premise_revision_only": True,
        "FGCQR_outcomes_inspected_before_revision": False,
        "SGBL_outcomes_inspected_before_revision": False,
        "explicit_phi_unit_normal_momentum": "0",
        "amplitude_candidates": [str(value) for value in amplitudes],
        "semantic_holdout_contract_sha256": _semantic_hash(holdout_contract),
    }


def _validate_sf1_protocol_v3(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the prospective, premise-only PROTO3 centre-buffer repair.

    The amendment changes no field equation, pulse, coupling, numerical method,
    stop, outcome rule, or robustness box. It replaces a round-number centre
    buffer that contradicted PROTO2's widest already-declared support with the
    exact positive vacuum neighbourhood that support actually leaves. The RUN1
    reproducer separately reads and verifies the immutable PROTO2 and ID1 blobs
    from the amendment commit; this pure validator locks their identifiers and
    hashes into the semantic contract and reduces every inherited field to the
    sealed PROTO2 validator.
    """

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
        "action",
        "dimensionless_parameters",
        "initial_data",
        "calibration",
        "gauge",
        "numerics",
        "health_stops",
        "analysis",
        "norms_and_estimators",
        "partition",
        "provenance",
        "robustness",
        "negative_scope",
        "classification",
        "revision_policy",
        "claims",
    }
    _keys("protocol", protocol, expected_root)
    identity = {
        "schema_version": protocol["schema_version"],
        "artifact_id": protocol["artifact_id"],
        "protocol_version": protocol["protocol_version"],
        "project_version": protocol["project_version"],
        "frozen": protocol["frozen"],
        "units": protocol["units"],
        "purpose": protocol["purpose"],
    }
    if identity != {
        "schema_version": 1,
        "artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "protocol_version": 3,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError(
            "PROTO3 identity, version, units, purpose, or frozen state differs"
        )

    amendment = _mapping("protocol.amendment", protocol["amendment"])
    _keys("protocol.amendment", amendment, set(SF1_PROTOCOL_V3_AMENDMENT))
    if dict(amendment) != SF1_PROTOCOL_V3_AMENDMENT:
        raise ValueError(
            "PROTO3 pre-holdout amendment or immutable diagnosis provenance differs"
        )

    initial_data = _mapping("protocol.initial_data", protocol["initial_data"])
    chi = _mapping("protocol.initial_data.chi", initial_data.get("chi"))
    partition = _mapping("protocol.partition", protocol["partition"])
    numerics = _mapping("protocol.numerics", protocol["numerics"])
    provenance = _mapping("protocol.provenance", protocol["provenance"])

    center = _fraction(
        "protocol.initial_data.chi.center_over_L0",
        chi.get("center_over_L0"),
        positive=True,
    )
    base_half_width = _fraction(
        "protocol.initial_data.chi.half_width_over_L0",
        chi.get("half_width_over_L0"),
        positive=True,
    )
    width_factors = _rational_strings(
        "protocol.partition.held_out_chi_width_factors",
        partition.get("held_out_chi_width_factors"),
    )
    phi_seed_factors = _rational_strings(
        "protocol.partition.held_out_phi_seed_factors",
        partition.get("held_out_phi_seed_factors"),
    )
    held_out_case_ids = _strings(
        "protocol.partition.held_out_case_ids",
        partition.get("held_out_case_ids"),
    )
    if width_factors != [Fraction(1), Fraction(7, 8), Fraction(9, 8)]:
        raise ValueError("PROTO3 held-out chi-width factors differ")
    if phi_seed_factors != [Fraction(1), Fraction(1, 2), Fraction(2)]:
        raise ValueError("PROTO3 held-out phi-seed factors differ")
    if held_out_case_ids != [
        "FGCQR-CENTRAL",
        "FGCQR-SEED-HALF",
        "FGCQR-SEED-DOUBLE",
        "FGCQR-WIDTH-SEVEN-EIGHTHS",
        "FGCQR-WIDTH-NINE-EIGHTHS",
    ]:
        raise ValueError("PROTO3 held-out case partition differs")
    widest_half_width = base_half_width * max(width_factors)
    exact_center_buffer = center - widest_half_width
    amended_floor = _fraction(
        "protocol.initial_data.chi.minimum_center_buffer_over_L0",
        chi.get("minimum_center_buffer_over_L0"),
        positive=True,
    )
    original_floor = _fraction(
        "protocol.amendment.original_minimum_center_buffer_over_L0",
        amendment["original_minimum_center_buffer_over_L0"],
        positive=True,
    )
    if (
        widest_half_width != Fraction(9, 16)
        or exact_center_buffer != Fraction(39, 16)
        or amended_floor != exact_center_buffer
        or original_floor - exact_center_buffer != Fraction(1, 16)
    ):
        raise ValueError("PROTO3 exact centre-buffer repair arithmetic differs")

    resolutions = numerics.get("resolutions")
    if not isinstance(resolutions, list) or not resolutions or any(
        not isinstance(item, int) or item < 2 for item in resolutions
    ):
        raise ValueError("PROTO3 numerical resolutions must be integer grids")
    radial_min = _fraction(
        "protocol.numerics.radial_domain_min_over_L0",
        numerics.get("radial_domain_min_over_L0"),
        nonnegative=True,
    )
    radial_max = _fraction(
        "protocol.numerics.radial_domain_max_over_L0",
        numerics.get("radial_domain_max_over_L0"),
        positive=True,
    )
    coarsest_spacing = (radial_max - radial_min) / (min(resolutions) - 1)
    coarsest_intervals = exact_center_buffer / coarsest_spacing
    if (
        coarsest_intervals != 78
        or amendment["widest_case_coarsest_grid_intervals"] != 78
    ):
        raise ValueError("PROTO3 widest-case coarsest-grid centre buffer differs")
    if (
        amendment["analytic_center_requirement"]
        != "strictly_positive_exact_vacuum_neighborhood"
    ):
        raise ValueError("PROTO3 analytic regular-centre requirement differs")
    _boolean(
        "protocol.amendment.static_all_case_initial_premises_required_before_GR0_dynamic_eligibility",
        amendment[
            "static_all_case_initial_premises_required_before_GR0_dynamic_eligibility"
        ],
        True,
    )
    if (
        partition.get("resolved_holdout_manifest")
        != "configs/fgc/fgc-1-pro3-hld1.toml"
    ):
        raise ValueError("PROTO3 resolved holdout manifest namespace differs")
    if provenance.get("holdout_output_root") != "runs/fgc-2-sf1/proto3/holdout":
        raise ValueError("PROTO3 holdout output namespace differs")

    normalized = deepcopy(dict(protocol))
    normalized["artifact_id"] = SF1_PROTOCOL_V2_ARTIFACT_ID
    normalized["protocol_version"] = 2
    normalized["amendment"] = deepcopy(SF1_PROTOCOL_V2_AMENDMENT)
    normalized["initial_data"]["chi"]["minimum_center_buffer_over_L0"] = "5/2"
    normalized["partition"]["resolved_holdout_manifest"] = (
        "configs/fgc/fgc-1-pro2-hld1.toml"
    )
    normalized["provenance"]["holdout_output_root"] = (
        "runs/fgc-2-sf1/proto2/holdout"
    )
    inherited = _validate_sf1_protocol_v2(normalized)

    holdout_contract = {
        key: protocol[key]
        for key in (
            "amendment",
            "action",
            "dimensionless_parameters",
            "initial_data",
            "calibration",
            "gauge",
            "numerics",
            "health_stops",
            "analysis",
            "norms_and_estimators",
            "partition",
            "provenance",
            "robustness",
            "negative_scope",
            "classification",
            "revision_policy",
        )
    }
    return {
        **inherited,
        "artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "protocol_version": 3,
        "amendment_validated": True,
        "predecessor_protocol_artifact_id": SF1_PROTOCOL_V2_ARTIFACT_ID,
        "diagnosis_artifact_id": "FGC-1-ID1-FAM1",
        "premise_revision_only": True,
        "PROTO2_amendment_validated": True,
        "PROTO1_and_PROTO2_preserved": True,
        "FGCQR_outcomes_inspected_before_revision": False,
        "SGBL_outcomes_inspected_before_revision": False,
        "exact_widest_center_buffer_over_L0": "39/16",
        "widest_case_coarsest_grid_intervals": 78,
        "static_all_case_initial_premises_required_before_GR0_dynamic_eligibility": True,
        "semantic_holdout_contract_sha256": _semantic_hash(holdout_contract),
    }


def _validate_sf1_protocol_v4(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the premise-only PROTO4 numerical-contract overlay.

    PROTO4 deliberately does not duplicate PROTO3's 270-line physical study
    contract.  It inherits every unlisted field by the exact predecessor-file
    digest and serializes only the changes required by CAL0.  Repository and
    reproduction checks independently read that predecessor and the immutable
    diagnosis commit; this pure validator locks the overlay itself.
    """

    protocol = _mapping("protocol", protocol)
    _keys(
        "protocol",
        protocol,
        {
            "schema_version",
            "artifact_id",
            "protocol_version",
            "project_version",
            "frozen",
            "units",
            "purpose",
            "amendment",
            "replacement",
            "stop_applicability",
            "numerical_admission",
            "claims",
        },
    )
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
        "artifact_id": SF1_PROTOCOL_V4_ARTIFACT_ID,
        "protocol_version": 4,
        "project_version": "0.11.0",
        "frozen": True,
        "units": "c=hbar=1; code length L0=4",
        "purpose": "outcome-neutral_constraint-controlled_spherical_mechanism_falsification",
    }:
        raise ValueError(
            "PROTO4 identity, version, units, purpose, or frozen state differs"
        )

    amendment = _mapping("protocol.amendment", protocol["amendment"])
    _keys("protocol.amendment", amendment, set(SF1_PROTOCOL_V4_AMENDMENT))
    if dict(amendment) != SF1_PROTOCOL_V4_AMENDMENT:
        raise ValueError("PROTO4 amendment or immutable CAL0 provenance differs")

    replacement = _mapping("protocol.replacement", protocol["replacement"])
    _keys(
        "protocol.replacement",
        replacement,
        {"calibration", "norms_and_estimators", "partition", "provenance"},
    )
    calibration = _mapping(
        "protocol.replacement.calibration", replacement["calibration"]
    )
    expected_calibration = {
        "all_RUN1_health_and_constraint_stops_apply": False,
        "stop_applicability_table_required": True,
        "three_nested_common_event_admission_required": True,
        "selected_case_requires_fresh_primary_and_comparator_runs": True,
        "development_runs_forbidden_as_calibration_evidence": True,
    }
    _keys(
        "protocol.replacement.calibration",
        calibration,
        set(expected_calibration),
    )
    if dict(calibration) != expected_calibration:
        raise ValueError("PROTO4 calibration replacement differs")
    norms = _mapping(
        "protocol.replacement.norms_and_estimators",
        replacement["norms_and_estimators"],
    )
    expected_norms = {
        "window_and_alias_policy": "compact_window_with_RMS_scale_top_band_power_budget_and_nested_convergence_last_bin_is_diagnostic_only",
        "constraint_normalization": "componentwise_max_of_absolute_unredefined_projection_term_sum_M_Pl_squared_over_L0_squared_and_1e-30",
        "definition_owner": "FGC-1-HLT2-MON2_must_serialize_executable_formulas_before_calibration",
    }
    _keys("protocol.replacement.norms_and_estimators", norms, set(expected_norms))
    if dict(norms) != expected_norms:
        raise ValueError("PROTO4 norm/estimator replacement differs")
    partition = _mapping(
        "protocol.replacement.partition", replacement["partition"]
    )
    if dict(partition) != {
        "resolved_holdout_manifest": "configs/fgc/fgc-1-pro4-hld1.toml"
    }:
        raise ValueError("PROTO4 resolved holdout manifest namespace differs")
    provenance = _mapping(
        "protocol.replacement.provenance", replacement["provenance"]
    )
    if dict(provenance) != {
        "calibration_output_root": "runs/fgc-2-sf1/proto4/calibration",
        "holdout_output_root": "runs/fgc-2-sf1/proto4/holdout",
        "both_output_roots_must_be_empty_before_their_first_authorized_run": True,
    }:
        raise ValueError("PROTO4 output namespaces or overwrite rule differ")

    applicability = _mapping(
        "protocol.stop_applicability", protocol["stop_applicability"]
    )
    expected_applicability_keys = {
        "legacy_monitor_artifact_id",
        "successor_monitor_artifact_id",
        "universal_runtime_stop_ids",
        "candidate_branch_only_stop_ids",
        "replaced_by_PROTO4_convergent_admission_ids",
        "all_legacy_stop_ids_partitioned_exactly_once",
        "universal_stops_apply_to_GR0_SGBL_and_FGCQR",
        "candidate_branch_stops_require_branch_owned_definitions_for_SGBL_and_FGCQR",
        "candidate_branch_stops_forbidden_as_GR0_calibration_rejection_reasons",
        "replaced_stops_remain_historical_controls_not_active_PROTO4_decisions",
    }
    _keys(
        "protocol.stop_applicability",
        applicability,
        expected_applicability_keys,
    )
    universal = _strings(
        "protocol.stop_applicability.universal_runtime_stop_ids",
        applicability["universal_runtime_stop_ids"],
    )
    candidate_only = _strings(
        "protocol.stop_applicability.candidate_branch_only_stop_ids",
        applicability["candidate_branch_only_stop_ids"],
    )
    replaced = _strings(
        "protocol.stop_applicability.replaced_by_PROTO4_convergent_admission_ids",
        applicability["replaced_by_PROTO4_convergent_admission_ids"],
    )
    if (
        applicability["legacy_monitor_artifact_id"] != "FGC-1-HLT1-MON1"
        or applicability["successor_monitor_artifact_id"] != "FGC-1-HLT2-MON2"
        or universal != list(PROTO4_UNIVERSAL_STOP_IDS)
        or candidate_only != list(PROTO4_CANDIDATE_BRANCH_STOP_IDS)
        or replaced != list(PROTO4_REPLACED_STOP_IDS)
        or len(set(universal + candidate_only + replaced)) != len(PROTO4_LEGACY_STOP_IDS)
        or set(universal + candidate_only + replaced) != set(PROTO4_LEGACY_STOP_IDS)
    ):
        raise ValueError("PROTO4 stop-applicability partition differs")
    for key in (
        "all_legacy_stop_ids_partitioned_exactly_once",
        "universal_stops_apply_to_GR0_SGBL_and_FGCQR",
        "candidate_branch_stops_require_branch_owned_definitions_for_SGBL_and_FGCQR",
        "candidate_branch_stops_forbidden_as_GR0_calibration_rejection_reasons",
        "replaced_stops_remain_historical_controls_not_active_PROTO4_decisions",
    ):
        _boolean(f"protocol.stop_applicability.{key}", applicability[key], True)

    admission = _mapping(
        "protocol.numerical_admission", protocol["numerical_admission"]
    )
    _keys(
        "protocol.numerical_admission",
        admission,
        {"spectrum", "constraints", "calibration"},
    )
    spectrum = _mapping(
        "protocol.numerical_admission.spectrum", admission["spectrum"]
    )
    expected_spectrum_keys = {
        "representation",
        "spatial_window",
        "temporal_window",
        "relative_amplitude_floor_for_last_bin_diagnostic",
        "absolute_amplitude_floor_for_last_bin_diagnostic",
        "unresolved_top_fraction",
        "maximum_top_band_field_power_fraction",
        "maximum_top_band_derivative_weighted_power_fraction",
        "maximum_nested_refinement_tail_ratio",
        "maximum_RMS_proper_frequency_over_Lambda",
        "maximum_RMS_proper_wavenumber_over_Lambda",
        "minimum_temporal_samples",
        "last_bin_alias_boolean_is_diagnostic_only",
        "zero_signal_has_zero_tail_and_RMS_scale",
        "every_nonzero_denominator_must_be_finite_and_positive",
    }
    _keys("protocol.numerical_admission.spectrum", spectrum, expected_spectrum_keys)
    if (
        spectrum["representation"]
        != "dimensionless_Minkowski_reference_deviations_alpha_minus_1_shift_over_r_lambda_minus_1_R_over_r_minus_1_phi_over_Lambda_chi_over_Lambda_with_analytic_center_limits"
        or spectrum["spatial_window"]
        != "C_infinity_compact_vacuum_buffer_window"
        or spectrum["temporal_window"] != "causal_past_only_proper_time_window"
        or _fraction(
            "protocol.numerical_admission.spectrum.relative_amplitude_floor_for_last_bin_diagnostic",
            spectrum["relative_amplitude_floor_for_last_bin_diagnostic"],
            positive=True,
        )
        != Fraction(1, 1099511627776)
        or _fraction(
            "protocol.numerical_admission.spectrum.absolute_amplitude_floor_for_last_bin_diagnostic",
            spectrum["absolute_amplitude_floor_for_last_bin_diagnostic"],
            positive=True,
        )
        != Fraction(1, 10**30)
        or _fraction(
            "protocol.numerical_admission.spectrum.unresolved_top_fraction",
            spectrum["unresolved_top_fraction"],
            positive=True,
        )
        != Fraction(1, 8)
        or _fraction(
            "protocol.numerical_admission.spectrum.maximum_top_band_field_power_fraction",
            spectrum["maximum_top_band_field_power_fraction"],
            positive=True,
        )
        != Fraction(1, 1048576)
        or _fraction(
            "protocol.numerical_admission.spectrum.maximum_top_band_derivative_weighted_power_fraction",
            spectrum["maximum_top_band_derivative_weighted_power_fraction"],
            positive=True,
        )
        != Fraction(1, 1024)
        or _fraction(
            "protocol.numerical_admission.spectrum.maximum_nested_refinement_tail_ratio",
            spectrum["maximum_nested_refinement_tail_ratio"],
            positive=True,
        )
        != Fraction(1, 4)
        or _fraction(
            "protocol.numerical_admission.spectrum.maximum_RMS_proper_frequency_over_Lambda",
            spectrum["maximum_RMS_proper_frequency_over_Lambda"],
            positive=True,
        )
        != Fraction(1, 4)
        or _fraction(
            "protocol.numerical_admission.spectrum.maximum_RMS_proper_wavenumber_over_Lambda",
            spectrum["maximum_RMS_proper_wavenumber_over_Lambda"],
            positive=True,
        )
        != Fraction(1, 4)
        or _integer(
            "protocol.numerical_admission.spectrum.minimum_temporal_samples",
            spectrum["minimum_temporal_samples"],
            minimum=16,
        )
        != 64
    ):
        raise ValueError("PROTO4 spectral admission differs")
    for key in (
        "last_bin_alias_boolean_is_diagnostic_only",
        "zero_signal_has_zero_tail_and_RMS_scale",
        "every_nonzero_denominator_must_be_finite_and_positive",
    ):
        _boolean(f"protocol.numerical_admission.spectrum.{key}", spectrum[key], True)

    constraints = _mapping(
        "protocol.numerical_admission.constraints", admission["constraints"]
    )
    expected_constraint_keys = {
        "fields",
        "normalization",
        "common_event_alignment",
        "coarsest_normalized_guard_max",
        "finest_normalized_guard_max",
        "minimum_observed_convergence_order",
        "minimum_nested_resolutions",
        "Richardson_error_enters_conservative_observable_error_sum",
        "constraint_error_must_not_be_subtracted_from_observable_margin",
        "single_resolution_smallness_cannot_establish_admission",
        "propagation_theorem_not_inferred_from_numerical_smallness",
    }
    _keys(
        "protocol.numerical_admission.constraints",
        constraints,
        expected_constraint_keys,
    )
    if (
        _strings(
            "protocol.numerical_admission.constraints.fields",
            constraints["fields"],
        )
        != [
            "Hamiltonian",
            "radial_momentum",
            "gauge_t",
            "gauge_r",
            "six_kinematic_reduction_fields",
        ]
        or constraints["normalization"]
        != "max_absolute_unredefined_projection_term_sum_or_M_Pl_squared_over_L0_squared_or_1e-30_componentwise"
        or constraints["common_event_alignment"]
        != "same_coordinate_time_for_GR0_calibration_and_same_affine_label_and_origin_for_DEF1"
        or _fraction(
            "protocol.numerical_admission.constraints.coarsest_normalized_guard_max",
            constraints["coarsest_normalized_guard_max"],
            positive=True,
        )
        != Fraction(1, 50)
        or _fraction(
            "protocol.numerical_admission.constraints.finest_normalized_guard_max",
            constraints["finest_normalized_guard_max"],
            positive=True,
        )
        != Fraction(1, 1000)
        or _fraction(
            "protocol.numerical_admission.constraints.minimum_observed_convergence_order",
            constraints["minimum_observed_convergence_order"],
            positive=True,
        )
        != Fraction(3, 2)
        or _integer(
            "protocol.numerical_admission.constraints.minimum_nested_resolutions",
            constraints["minimum_nested_resolutions"],
            minimum=3,
        )
        != 3
    ):
        raise ValueError("PROTO4 constraint admission differs")
    for key in (
        "Richardson_error_enters_conservative_observable_error_sum",
        "constraint_error_must_not_be_subtracted_from_observable_margin",
        "single_resolution_smallness_cannot_establish_admission",
        "propagation_theorem_not_inferred_from_numerical_smallness",
    ):
        _boolean(
            f"protocol.numerical_admission.constraints.{key}",
            constraints[key],
            True,
        )

    calibration_admission = _mapping(
        "protocol.numerical_admission.calibration", admission["calibration"]
    )
    expected_calibration_admission = {
        "static_initial_admission_required_for_every_declared_candidate": True,
        "dynamic_candidates_run_in_strict_declared_amplitude_order": True,
        "first_eligible_candidate_stops_later_calibration_runs": True,
        "minimum_consecutive_trapped_common_events": 8,
        "trapped_sign_margin_over_total_error_factor": 4,
        "primary_and_comparator_agreement_required_for_selected_candidate": True,
        "no_candidate_outcome_can_change_this_protocol": True,
    }
    _keys(
        "protocol.numerical_admission.calibration",
        calibration_admission,
        set(expected_calibration_admission),
    )
    if dict(calibration_admission) != expected_calibration_admission:
        raise ValueError("PROTO4 calibration admission differs")

    claims = _mapping("protocol.claims", protocol["claims"])
    if dict(claims) != PROTO4_CLAIMS:
        raise ValueError("PROTO4 claims must remain fail-closed")

    overlay_contract = {
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "replacement": replacement,
        "stop_applicability": applicability,
        "numerical_admission": admission,
    }
    return {
        "artifact_id": SF1_PROTOCOL_V4_ARTIFACT_ID,
        "protocol_version": 4,
        "frozen": True,
        "outcome_neutral_contract_validated": True,
        "premise_revision_only": True,
        "overlay_schema_and_hash_binding_validated": True,
        "predecessor_protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "diagnosis_artifact_id": "FGC-1-CAL0-PREF2",
        "diagnosis_checkpoint_commit": amendment["diagnosis_checkpoint_commit"],
        "GR0_exploratory_outcomes_disclosed": True,
        "development_runs_forbidden_as_calibration_evidence": True,
        "FGCQR_outcomes_inspected_before_revision": False,
        "SGBL_outcomes_inspected_before_revision": False,
        "amplitude_candidates": ["2", "5/2", "3", "7/2", "4", "9/2", "5"],
        "held_out_case_ids": [
            "FGCQR-CENTRAL",
            "FGCQR-SEED-HALF",
            "FGCQR-SEED-DOUBLE",
            "FGCQR-WIDTH-SEVEN-EIGHTHS",
            "FGCQR-WIDTH-NINE-EIGHTHS",
        ],
        "legacy_stop_count": len(PROTO4_LEGACY_STOP_IDS),
        "branch_specific_stop_partition_validated": True,
        "weighted_spectral_admission_validated": True,
        "common_event_constraint_admission_validated": True,
        "resolved_holdout_manifest": partition["resolved_holdout_manifest"],
        "calibration_output_root": provenance["calibration_output_root"],
        "holdout_output_root": provenance["holdout_output_root"],
        "resolved_holdout_manifest_present": False,
        "semantic_holdout_contract_sha256": _semantic_hash(overlay_contract),
        "claims": claims,
    }


def validate_sf1_protocol(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Validate immutable PROTO1--PROTO3 history or frozen PROTO4 overlay."""

    candidate = _mapping("protocol", protocol)
    artifact_id = candidate.get("artifact_id")
    if artifact_id == SF1_PROTOCOL_V1_ARTIFACT_ID:
        return _validate_sf1_protocol_v1(candidate)
    if artifact_id == SF1_PROTOCOL_V2_ARTIFACT_ID:
        return _validate_sf1_protocol_v2(candidate)
    if artifact_id == SF1_PROTOCOL_ARTIFACT_ID:
        return _validate_sf1_protocol_v3(candidate)
    if artifact_id == SF1_PROTOCOL_V4_ARTIFACT_ID:
        return _validate_sf1_protocol_v4(candidate)
    raise ValueError("unknown SF1 protocol artifact")


def _artifact(name: str, value: Mapping[str, Any], expected_id: str) -> None:
    if value.get("artifact_id") != expected_id:
        raise ValueError(f"{name} artifact identity differs")
    nonclaims = _mapping(f"{name}.nonclaims", value.get("nonclaims"))
    if not nonclaims or any(item is not False for item in nonclaims.values()):
        raise ValueError(f"{name} contains an absent or promoted nonclaim ledger")


def _gate(name: str, value: Mapping[str, Any], key: str, expected: bool) -> None:
    status = _mapping(f"{name}.gate_status", value.get("gate_status"))
    if status.get(key) is not expected:
        raise ValueError(f"{name} gate {key} must be {str(expected).lower()}")


def _future_evidence_status(
    name: str,
    value: Mapping[str, Any] | None,
    protocol_certificate: Mapping[str, Any],
) -> dict[str, Any]:
    expected_id, expected_gate = FUTURE_EVIDENCE_REQUIREMENTS[name]
    if value is None:
        return {
            "artifact_id": expected_id,
            "required_gate": expected_gate,
            "status": "absent",
            "passed": False,
        }
    value = _mapping(f"future_evidence.{name}", value)
    _keys(
        f"future_evidence.{name}",
        value,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "classification",
            "generated_by",
            "derivation_document",
            "derivation_document_sha256",
            "source_config_sha256",
            "predecessor_sha256",
            "implementation_sha256",
            "scope_bindings",
            "certificate_contract",
            "artifact_payload",
            "gate_status",
            "nonclaims",
        },
    )
    if value["schema_version"] != 1 or value["project_version"] != "0.11.0":
        raise ValueError(f"future_evidence.{name} schema or project version differs")
    if value["classification"] != FUTURE_EVIDENCE_CLASSIFICATIONS[name]:
        raise ValueError(f"future_evidence.{name} classification differs")
    _nonempty_string(f"future_evidence.{name}.generated_by", value["generated_by"])
    _nonempty_string(
        f"future_evidence.{name}.derivation_document", value["derivation_document"]
    )
    _sha256_string(
        f"future_evidence.{name}.derivation_document_sha256",
        value["derivation_document_sha256"],
    )
    for ledger_name in (
        "source_config_sha256",
        "predecessor_sha256",
        "implementation_sha256",
    ):
        ledger = _mapping(
            f"future_evidence.{name}.{ledger_name}", value[ledger_name]
        )
        if not ledger:
            raise ValueError(f"future_evidence.{name}.{ledger_name} must not be empty")
        for path, digest in ledger.items():
            _nonempty_string(f"future_evidence.{name}.{ledger_name}.path", path)
            _sha256_string(
                f"future_evidence.{name}.{ledger_name}.{path}", digest
            )

    certificate_contract = _mapping(
        f"future_evidence.{name}.certificate_contract", value["certificate_contract"]
    )
    _keys(
        f"future_evidence.{name}.certificate_contract",
        certificate_contract,
        {
            "artifact_specific_payload_validated",
            "canonical_reproduction_passed",
            "scope_bindings_verified",
            "outcome_data_not_used_to_select_certificate_contract",
        },
    )
    for key, item in certificate_contract.items():
        _boolean(f"future_evidence.{name}.certificate_contract.{key}", item, True)

    scope = _mapping(f"future_evidence.{name}.scope_bindings", value["scope_bindings"])
    _keys(f"future_evidence.{name}.scope_bindings", scope, FUTURE_SCOPE_BINDING_KEYS)
    expected_scope_identity = {
        "run1_artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "protocol_semantic_holdout_contract_sha256": protocol_certificate[
            "semantic_holdout_contract_sha256"
        ],
        "action_artifact_id": "FGC-1-ACT1",
        "variation_artifact_id": "FGC-1-VAR1",
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
    }
    if any(scope[key] != expected for key, expected in expected_scope_identity.items()):
        raise ValueError(f"future_evidence.{name} belongs to a different RUN1 scope")
    for key in FUTURE_SCOPE_BINDING_KEYS - set(expected_scope_identity):
        _sha256_string(f"future_evidence.{name}.scope_bindings.{key}", scope[key])

    _artifact(f"future_evidence.{name}", value, expected_id)
    nonclaims = _mapping(f"future_evidence.{name}.nonclaims", value["nonclaims"])
    for claim, expected in FUTURE_COMMON_NONCLAIMS.items():
        if nonclaims.get(claim) is not expected:
            raise ValueError(
                f"future_evidence.{name} omits or promotes common nonclaim {claim}"
            )
    payload = _mapping(f"future_evidence.{name}.artifact_payload", value["artifact_payload"])
    if not payload:
        raise ValueError(f"future_evidence.{name}.artifact_payload must not be empty")
    if name == "protocol_holdout":
        _keys(
            "future_evidence.protocol_holdout.artifact_payload",
            payload,
            {
                "selected_GR0_amplitude",
                "calibration_candidate_ledger",
                "static_all_case_initial_premise_ledger",
                "ordered_held_out_case_ids",
                "held_out_case_config_sha256",
                "repository_sequence_receipt",
            },
        )
        selected = payload["selected_GR0_amplitude"]
        if selected not in protocol_certificate["amplitude_candidates"]:
            raise ValueError("PRO3-HLD1 selected GR-0 amplitude is not a protocol candidate")
        candidate_ledger = payload["calibration_candidate_ledger"]
        if not isinstance(candidate_ledger, list) or len(candidate_ledger) != len(
            protocol_certificate["amplitude_candidates"]
        ):
            raise ValueError("PRO3-HLD1 calibration ledger does not cover every candidate")
        eligible_amplitudes: list[str] = []
        for index, item in enumerate(candidate_ledger):
            item = _mapping(f"PRO3-HLD1.calibration_candidate_ledger[{index}]", item)
            _keys(
                f"PRO3-HLD1.calibration_candidate_ledger[{index}]",
                item,
                {
                    "amplitude",
                    "config_sha256",
                    "result_sha256",
                    "eligible",
                    "eligibility_reason",
                },
            )
            if item["amplitude"] != protocol_certificate["amplitude_candidates"][index]:
                raise ValueError("PRO3-HLD1 calibration candidate order differs")
            _sha256_string("PRO3-HLD1 calibration config hash", item["config_sha256"])
            _sha256_string("PRO3-HLD1 calibration result hash", item["result_sha256"])
            eligible = _boolean("PRO3-HLD1 calibration eligibility", item["eligible"])
            _nonempty_string("PRO3-HLD1 eligibility reason", item["eligibility_reason"])
            if eligible:
                eligible_amplitudes.append(item["amplitude"])
        if not eligible_amplitudes or selected != eligible_amplitudes[0]:
            raise ValueError("PRO3-HLD1 does not select the first eligible GR-0 amplitude")
        expected_static_cases = [
            ("FGCQR-CENTRAL", "1", "1", "5/2"),
            ("FGCQR-SEED-HALF", "1", "1/2", "5/2"),
            ("FGCQR-SEED-DOUBLE", "1", "2", "5/2"),
            ("FGCQR-WIDTH-SEVEN-EIGHTHS", "7/8", "1", "41/16"),
            ("FGCQR-WIDTH-NINE-EIGHTHS", "9/8", "1", "39/16"),
        ]
        static_ledger = payload["static_all_case_initial_premise_ledger"]
        if not isinstance(static_ledger, list) or len(static_ledger) != len(
            expected_static_cases
        ):
            raise ValueError("PRO3-HLD1 static initial-premise ledger is incomplete")
        for index, (item, expected_case) in enumerate(
            zip(static_ledger, expected_static_cases, strict=True)
        ):
            item = _mapping(
                f"PRO3-HLD1.static_all_case_initial_premise_ledger[{index}]",
                item,
            )
            _keys(
                f"PRO3-HLD1.static_all_case_initial_premise_ledger[{index}]",
                item,
                {
                    "case_id",
                    "chi_width_factor",
                    "phi_seed_factor",
                    "exact_center_buffer_over_L0",
                    "config_sha256",
                    "initial_data_result_sha256",
                    "constraint_compatible",
                    "regular_center",
                    "finite_mass",
                    "initial_compactness_in_protocol_window",
                    "no_initial_trapped_sphere",
                    "static_premises_passed",
                    "FGCQR_evolution_outcome_inspected",
                },
            )
            identity = (
                item["case_id"],
                item["chi_width_factor"],
                item["phi_seed_factor"],
                item["exact_center_buffer_over_L0"],
            )
            if identity != expected_case:
                raise ValueError("PRO3-HLD1 static initial-premise case differs")
            _sha256_string(
                "PRO3-HLD1 static case config hash", item["config_sha256"]
            )
            _sha256_string(
                "PRO3-HLD1 static initial-data result hash",
                item["initial_data_result_sha256"],
            )
            for key in (
                "constraint_compatible",
                "regular_center",
                "finite_mass",
                "initial_compactness_in_protocol_window",
                "no_initial_trapped_sphere",
                "static_premises_passed",
            ):
                _boolean(
                    f"PRO3-HLD1 static initial premise {key}", item[key], True
                )
            _boolean(
                "PRO3-HLD1 static FGCQR outcome access",
                item["FGCQR_evolution_outcome_inspected"],
                False,
            )
        if payload["ordered_held_out_case_ids"] != protocol_certificate["held_out_case_ids"]:
            raise ValueError("PRO3-HLD1 held-out case order differs from active PROTO3")
        case_hashes = _mapping(
            "PRO3-HLD1.held_out_case_config_sha256",
            payload["held_out_case_config_sha256"],
        )
        if set(case_hashes) != set(protocol_certificate["held_out_case_ids"]):
            raise ValueError("PRO3-HLD1 held-out case hashes are incomplete")
        for case, digest in case_hashes.items():
            _sha256_string(f"PRO3-HLD1 held-out case {case}", digest)
        receipt = _mapping(
            "PRO3-HLD1.repository_sequence_receipt",
            payload["repository_sequence_receipt"],
        )
        _keys(
            "PRO3-HLD1.repository_sequence_receipt",
            receipt,
            {
                "preexisting_FGCQR_output_paths",
                "runner_refuses_overwrite",
                "append_only_event_log_sha256",
                "sealed_before_FGCQR_execution",
                "human_private_observation_not_machine_verifiable",
            },
        )
        if receipt["preexisting_FGCQR_output_paths"] != []:
            raise ValueError("PRO3-HLD1 found pre-existing FGC-QR output")
        for key in (
            "runner_refuses_overwrite",
            "sealed_before_FGCQR_execution",
            "human_private_observation_not_machine_verifiable",
        ):
            _boolean(f"PRO3-HLD1.repository_sequence_receipt.{key}", receipt[key], True)
        _sha256_string(
            "PRO3-HLD1.repository_sequence_receipt.append_only_event_log_sha256",
            receipt["append_only_event_log_sha256"],
        )
    else:
        true_fields = FUTURE_PAYLOAD_TRUE_FIELDS[name]
        _keys(
            f"future_evidence.{name}.artifact_payload",
            payload,
            {
                "run_envelope_semantic_sha256",
                "quantitative_evidence",
                *true_fields,
            },
        )
        if payload["run_envelope_semantic_sha256"] != scope[
            "declared_run_envelope_sha256"
        ]:
            raise ValueError(
                f"future_evidence.{name} payload and scope bind different run envelopes"
            )
        quantitative = _mapping(
            f"future_evidence.{name}.artifact_payload.quantitative_evidence",
            payload["quantitative_evidence"],
        )
        if not quantitative:
            raise ValueError(
                f"future_evidence.{name} quantitative evidence must not be empty"
            )
        for key in true_fields:
            _boolean(
                f"future_evidence.{name}.artifact_payload.{key}", payload[key], True
            )
    _gate(f"future_evidence.{name}", value, expected_gate, True)
    return {
        "artifact_id": expected_id,
        "required_gate": expected_gate,
        "status": "scope_bound_hash_bound_gate_passed",
        "declared_run_envelope_sha256": scope["declared_run_envelope_sha256"],
        "passed": True,
    }


def _predicate(
    identifier: str,
    passed: bool,
    evidence: Sequence[str],
    failure_consequence: str,
) -> dict[str, Any]:
    return {
        "id": identifier,
        "passed": passed,
        "evidence": list(evidence),
        "failure_consequence": failure_consequence,
    }


def scoped_classical_spherical_run_audit(
    *,
    action: Mapping[str, Any],
    variation: Mapping[str, Any],
    uhyp1: Mapping[str, Any],
    con2: Mapping[str, Any],
    con3: Mapping[str, Any],
    bnd1: Mapping[str, Any],
    def0: Mapping[str, Any],
    eft1: Mapping[str, Any],
    protocol: Mapping[str, Any],
    future_evidence: Mapping[str, Mapping[str, Any] | None],
) -> dict[str, Any]:
    """Compose inherited evidence and future certificates into RUN1."""

    action = _mapping("ACT1", action)
    variation = _mapping("VAR1", variation)
    uhyp1 = _mapping("DOM3-UHYP1", uhyp1)
    con2 = _mapping("CON2-MPROP1", con2)
    con3 = _mapping("CON3-CAU1", con3)
    bnd1 = _mapping("BND1-MD1", bnd1)
    def0 = _mapping("DEF0-OBS1", def0)
    eft1 = _mapping("EFT1-OPEN1", eft1)
    protocol_certificate = validate_sf1_protocol(protocol)

    _artifact("ACT1", action, "FGC-1-ACT1")
    _artifact("VAR1", variation, "FGC-1-VAR1")
    _artifact("DOM3-UHYP1", uhyp1, "FGC-1-HYP1-DOM3-UHYP1")
    _artifact("CON2-MPROP1", con2, "FGC-1-HYP1-CON2-MPROP1")
    _artifact("CON3-CAU1", con3, "FGC-1-HYP1-CON3-CAU1")
    _artifact("BND1-MD1", bnd1, "FGC-1-HYP1-BND1-MD1")
    _artifact("DEF0-OBS1", def0, "FGC-1-DEF0-OBS1")
    _artifact("EFT1-OPEN1", eft1, "FGC-1-EFT1-OPEN1")

    _gate("ACT1", action, "declared_action_and_scalar_algebra_certificate_passed", True)
    _gate("VAR1", variation, "var1_covariant_metric_equation_derived_and_cross_checked", True)
    _gate(
        "DOM3-UHYP1",
        uhyp1,
        "compact_nonflat_radial_strong_hyperbolicity_on_implicit_branch_graph_passed",
        True,
    )
    _gate("DOM3-UHYP1", uhyp1, "multidirectional_strong_hyperbolicity_proven", False)
    _gate("CON2-MPROP1", con2, "local_differential_identity", True)
    _gate("CON2-MPROP1", con2, "complete_kinematic_1plus1_reduction_subsidiary", True)
    _gate(
        "CON3-CAU1",
        con3,
        "conditional_boundary_free_gauge_Cauchy_uniqueness_theorem_derived",
        True,
    )
    _gate("CON3-CAU1", con3, "constraint_preserving_ACT1_IBVP_proven", False)
    _gate(
        "BND1-MD1",
        bnd1,
        "uniform_frozen_radial_main_system_boundary_dissipation_passed",
        True,
    )
    _gate("BND1-MD1", bnd1, "constraint_preserving_ACT1_IBVP_proven", False)
    _gate("DEF0-OBS1", def0, "exact_metric_null_observable_contract_and_controls_passed", True)
    _gate("DEF0-OBS1", def0, "metric_null_affine_defocusing_DEF1_derived", False)
    _gate("EFT1-OPEN1", eft1, "audit_definition_and_fail_closed_composition_passed", True)
    _gate("EFT1-OPEN1", eft1, "retained_EFT_open_run_envelope_passed", False)
    _gate("EFT1-OPEN1", eft1, "evolution_authorized", False)

    fixtures = action.get("fixtures")
    if not isinstance(fixtures, list) or [item.get("model_id") for item in fixtures] != [
        "GR-0",
        "SGB-L",
        "FGC-QR",
    ]:
        raise ValueError("ACT1 frozen branch order differs")
    for fixture in fixtures:
        action_gate = _mapping(f"ACT1.{fixture['model_id']}.action_gate", fixture.get("action_gate"))
        if any(value is not True for value in action_gate.values()):
            raise ValueError(f"ACT1 {fixture['model_id']} action gate is incomplete")

    eft_audit = _mapping("EFT1-OPEN1.open_run_authorization_audit", eft1.get("open_run_authorization_audit"))
    if (
        eft_audit.get("passed_predicate_count") != 3
        or eft_audit.get("required_predicate_count") != 12
        or eft_audit.get("retained_EFT_open_run_envelope_passed") is not False
        or eft_audit.get("evolution_authorized") is not False
    ):
        raise ValueError("EFT1 negative 3/12 authorization state differs")

    future_evidence = _mapping("future_evidence", future_evidence)
    if set(future_evidence) != set(FUTURE_EVIDENCE_REQUIREMENTS):
        raise ValueError("future evidence slots differ")
    future = {
        name: _future_evidence_status(
            name, future_evidence[name], protocol_certificate
        )
        for name in FUTURE_EVIDENCE_REQUIREMENTS
    }
    envelope_hashes = {
        item["declared_run_envelope_sha256"]
        for item in future.values()
        if item["passed"]
    }
    if len(envelope_hashes) > 1:
        raise ValueError("future evidence certificates bind different run envelopes")
    protocol_certificate = dict(protocol_certificate)
    protocol_certificate["resolved_holdout_manifest_present"] = future[
        "protocol_holdout"
    ]["passed"]

    predicates = [
        _predicate(
            "frozen_outcome_neutral_protocol_and_resolved_holdout_hash",
            future["protocol_holdout"]["passed"],
            [
                "PROTO3 schema, immutable amendment provenance, and semantic holdout contract hash validate",
                f"PRO3-HLD1 status={future['protocol_holdout']['status']}",
            ],
            "the FGC-QR cases are not fully resolved and sealed before outcome access",
        ),
        _predicate(
            "classical_spherical_domain_branch_and_multidirectional_early_kill_envelope",
            future["run_domain"]["passed"] and future["multidirectional_health"]["passed"],
            [
                "DOM3 supplies a compact radial branch graph only",
                f"DOM4-RUN1 status={future['run_domain']['status']}",
                f"HYP2-MD1 status={future['multidirectional_health']['status']}",
            ],
            "the intended run can leave the healthy branch or fail in an untested covector direction",
        ),
        _predicate(
            "nonlinear_unredefined_REF1_source_and_acceleration_solver_verified",
            future["nonlinear_source"]["passed"],
            [
                "ACT1/VAR1 equations are frozen and REF1 is the gauge formulation",
                f"SRC1-NL1 status={future['nonlinear_source']['status']}",
            ],
            "a numerical trajectory would not be a verified solution of the declared equations",
        ),
        _predicate(
            "physical_gauge_and_reduction_constraint_system_closed",
            future["physical_constraints"]["passed"],
            [
                "CON2 supplies local identities and CON3 only conditional boundary-free uniqueness",
                f"CON4-PHY1 status={future['physical_constraints']['status']}",
            ],
            "constraint drift could masquerade as the proposed mechanism",
        ),
        _predicate(
            "regular_center_finite_mass_constraint_compatible_initial_data_family",
            future["regular_center"]["passed"] and future["initial_data"]["passed"],
            [
                "COMP1 is one activated two-jet at r=4, not initial data",
                f"CTR1-REG1 status={future['regular_center']['status']}",
                f"ID1-FAM1 status={future['initial_data']['status']}",
            ],
            "the collapse problem lacks a regular centre or a nonzero-width admissible data family",
        ),
        _predicate(
            "spherical_boundary_or_domain_of_dependence_control",
            future["boundary_control"]["passed"],
            [
                "BND1 supplies frozen annular main-system dissipation only",
                f"BND2-CP1 status={future['boundary_control']['status']}",
            ],
            "incoming or reflected boundary error could contaminate the measured interval",
        ),
        _predicate(
            "classical_health_and_typed_fail_closed_monitoring",
            future["health_monitor"]["passed"],
            [
                "PROTO3 preserves the frozen branch, constraint, characteristic, invariant, scale, and boundary stops",
                f"HLT1-MON1 status={future['health_monitor']['status']}",
            ],
            "a run could cross a failed premise and continue as if it remained interpretable",
        ),
        _predicate(
            "independent_solver_validation_and_affine_measurement_contract",
            future["numerical_validation"]["passed"],
            [
                "DEF0 supplies the affine physical-null observable but no ACT1 trajectory",
                f"NUM1-VAL1 status={future['numerical_validation']['status']}",
            ],
            "implementation or measurement error could be mistaken for a scientific outcome",
        ),
    ]
    missing = [item["id"] for item in predicates if item["passed"] is not True]
    classical_authorized = not missing
    retained_eft_authorized = False
    physical_transition_authorized = False

    if retained_eft_authorized and not classical_authorized:
        raise ValueError("retained-EFT authorization cannot bypass scoped prerequisites")
    if physical_transition_authorized and not (
        retained_eft_authorized and classical_authorized
    ):
        raise ValueError("physical-transition authorization violates the implication contract")

    physical_predicates = [
        "retained_EFT_evolution_authorized",
        "COL1_converged_constraint_controlled_evolution",
        "finite_trapped_interval_with_fixed_orientation",
        "affine_null_and_complete_Raychaudhuri_route_agreement",
        "strict_positive_complete_Raychaudhuri_margin",
        "combined_error_budget_below_one_quarter_margin",
        "three_resolution_convergence",
        "independent_method_agreement",
        "no_health_constraint_branch_scale_or_boundary_stop",
        "ROB1_nonzero_open_neighborhood_persistence",
        "finite_invariant_transition_layer_derived",
        "regular_trapped_to_expanding_causal_continuation_derived",
        "exterior_junction_and_quasilocal_flux_closure_derived",
        "generalized_entropy_ledger_closed",
        "upper_closure_or_finite_handoff_derived",
    ]

    return {
        "artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "classification": "fail_closed_scoped_classical_spherical_diagnostic_authorization_contract_not_EFT_validity_evolution_or_physical_transition_certificate",
        "authorization_logic": "all_eight_scoped_predicates_must_pass; one_false_predicate_stops_FGCQR_holdout_execution",
        "protocol_validation": protocol_certificate,
        "validated_predecessor_evidence": {
            "ACT1_VAR1_unredefined_equations_frozen": True,
            "DOM3_compact_radial_branch_graph": True,
            "CON2_metric_gauge_and_kinematic_identities": True,
            "CON3_conditional_boundary_free_gauge_uniqueness": True,
            "BND1_frozen_radial_main_system_dissipation": True,
            "DEF0_affine_physical_null_observable_contract": True,
            "EFT1_negative_3_of_12_state_preserved": True,
        },
        "future_evidence_status": future,
        "authorization_audits": {
            "classical_spherical_diagnostic": {
                "logic": "all_required_predicates_must_pass",
                "required_predicates": predicates,
                "passed_predicate_count": len(predicates) - len(missing),
                "required_predicate_count": len(predicates),
                "missing_predicate_ids": missing,
                "authorized": classical_authorized,
            },
            "retained_EFT_evolution": {
                "logic": "delegated_without_relaxation_to_FGC-1-EFT1-OPEN1",
                "EFT1_passed_predicate_count": 3,
                "EFT1_required_predicate_count": 12,
                "RUN1_version_locked_to_current_EFT1_negative_state": True,
                "authorized": retained_eft_authorized,
            },
            "physical_transition_claim": {
                "logic": "reserved_false_in_RUN1; a_future_artifact_requires_all_local_global_flux_entropy_handoff_predicates_plus_both_prior_authorizations",
                "required_predicate_ids": physical_predicates,
                "missing_predicate_ids": physical_predicates,
                "RUN1_can_never_set_this_gate_true": True,
                "authorized": physical_transition_authorized,
            },
        },
        "gate_status": {
            "authorization_contract_and_protocol_validation_passed": True,
            "classical_spherical_diagnostic_authorized": classical_authorized,
            "retained_EFT_evolution_authorized": retained_eft_authorized,
            "physical_transition_claim_authorized": physical_transition_authorized,
            "FGCQR_holdout_execution_authorized": classical_authorized,
        },
        "claim_scope": {
            "allowed_only_if_scoped_gate_later_becomes_true": [
                "bounded classical spherical diagnostic of the frozen FGC-QR equations",
                "outcome-neutral numerical mechanism test",
                "typed fail-closed stop classification",
            ],
            "never_implied_by_RUN1": [
                "retained EFT validity",
                "UV completion",
                "physical transition",
                "singularity resolution",
                "child domain",
                "dark-sector mechanism",
                "varying locally measured c",
                "observational applicability",
            ],
        },
        "stop_mask": {
            "execute_FGCQR_holdout": classical_authorized,
            "report_scoped_classical_mechanism_result": False,
            "report_retained_EFT_evolution": False,
            "report_physical_transition": False,
            "report_singularity_resolution": False,
        },
    }
