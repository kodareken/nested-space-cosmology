#!/usr/bin/env python3
"""Regenerate the FGC-1-HLT1-MON1 transactional monitor certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass, replace
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite, pi
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt1-mon1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt1-mon1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt1-mon1.md"

from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    CONSTRAINT_NORMALIZATION_FLOOR,
    FIELD_MASS_DIMENSIONS,
    FIELD_ORDER,
    STOP_PRIORITY,
    ClassicalHealthMonitor,
    ClassicalHealthStop,
    HealthThresholds,
    StageHealthSnapshot,
    compact_vacuum_buffer_window,
    curvature_scale_and_ratio,
    dimensionless_state_components,
    induced_scaled_infinity_norm,
    normalized_branch_displacement,
    normalized_componentwise_infinity,
    windowed_spectral_support,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_COMMON_NONCLAIMS,
    FUTURE_EVIDENCE_CLASSIFICATIONS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-HLT1-MON1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "classical_health_and_typed_stop_monitoring_verified"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "run1_config": "configs/fgc/fgc-1-run1-sym1.toml",
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "nonlinear_source_config": "configs/fgc/fgc-1-src1-nl1.toml",
    "nonlinear_source_result": "results/fgc-1-src1-nl1.json",
    "physical_constraints_config": "configs/fgc/fgc-1-con4-phy1.toml",
    "physical_constraints_result": "results/fgc-1-con4-phy1.json",
    "run_domain_config": "configs/fgc/fgc-1-dom4-run1.toml",
    "run_domain_result": "results/fgc-1-dom4-run1.json",
    "multidirectional_health_config": "configs/fgc/fgc-1-hyp2-md1.toml",
    "multidirectional_health_result": "results/fgc-1-hyp2-md1.json",
    "boundary_control_config": "configs/fgc/fgc-1-bnd2-cp1.toml",
    "boundary_control_result": "results/fgc-1-bnd2-cp1.json",
}
RESULT_REQUIREMENTS = {
    "action_result": ("FGC-1-ACT1", "declared_action_and_scalar_algebra_certificate_passed"),
    "variation_result": ("FGC-1-VAR1", "var1_covariant_metric_equation_derived_and_cross_checked"),
    "nonlinear_source_result": ("FGC-1-SRC1-NL1", "nonlinear_REF1_source_and_branch_solver_verified"),
    "physical_constraints_result": ("FGC-1-CON4-PHY1", "physical_gauge_reduction_constraint_system_closed"),
    "run_domain_result": ("FGC-1-DOM4-RUN1", "nonzero_classical_spherical_run_envelope_passed"),
    "multidirectional_health_result": ("FGC-1-HYP2-MD1", "quantitative_all_covector_weak_coupling_health_envelope_passed"),
    "boundary_control_result": ("FGC-1-BND2-CP1", "spherical_boundary_or_domain_of_dependence_control_passed"),
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/evolution/boundary_domain.py",
        "src/recursive_horizons/fgc/evolution/health_monitor.py",
        "scripts/reproduce_fgc_hlt1_mon1.py",
    )
)


def _pairs(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
        answer[key] = value
    return answer


def _mapping(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a table")
    return value


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(
            f"{name} keys differ: missing={sorted(expected - set(value))}, "
            f"extra={sorted(set(value) - expected)}"
        )


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _fraction(name: str, value: Any, *, positive: bool = False, nonnegative: bool = False) -> Fraction:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a canonical rational string")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a canonical rational string") from exc
    if value != _fraction_text(answer):
        raise ValueError(f"{name} must use canonical rational encoding")
    if positive and answer <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and answer < 0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"{name} must name {expected}")
    path = (REPOSITORY / expected).resolve()
    if not path.is_file() or path.relative_to(REPOSITORY).as_posix() != expected:
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return path


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _load_json(path: Path, name: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _serial(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_serial(item) for item in value]
    if isinstance(value, np.generic):
        return _serial(value.item())
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("nonfinite value cannot enter canonical evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported evidence value {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(_serial(value), indent=2, sort_keys=True, allow_nan=False) + "\n"


def _threshold_values(raw: Mapping[str, Any]) -> tuple[dict[str, Fraction | int], HealthThresholds]:
    table = _mapping("thresholds", raw["thresholds"])
    expected = {
        "newton_residual_infinity_max", "newton_iteration_max",
        "kinetic_condition_number_max", "normalized_branch_displacement_max",
        "acceleration_infinity_max", "companion_deformation_operator_2_max",
        "action_energy_deformation_operator_2_max", "kinetic_deformation_operator_2_max",
        "effective_planck_ratio_min", "characteristic_imaginary_part_max",
        "eigenframe_condition_number_max", "symmetrizer_coercivity_min",
        "normalized_constraint_infinity_max", "minimum_observed_constraint_convergence_order",
        "phi_over_Lambda_max", "proper_frequency_over_Lambda_max",
        "proper_wavenumber_over_Lambda_max", "curvature_scale_squared_over_Lambda_squared_max",
        "boundary_margin_over_required_buffer_min",
    }
    _keys("thresholds", table, expected)
    values: dict[str, Fraction | int] = {
        key: (table[key] if key == "newton_iteration_max" else _fraction(f"thresholds.{key}", table[key], positive=key != "boundary_margin_over_required_buffer_min", nonnegative=key == "boundary_margin_over_required_buffer_min"))
        for key in table
    }
    expected_values = {
        "newton_residual_infinity_max": Q(1, 10**12),
        "newton_iteration_max": 16,
        "kinetic_condition_number_max": Q(10**10),
        "normalized_branch_displacement_max": Q(1, 16),
        "acceleration_infinity_max": Q(4),
        "companion_deformation_operator_2_max": Q(1, 4096),
        "action_energy_deformation_operator_2_max": Q(1, 8192),
        "kinetic_deformation_operator_2_max": Q(1, 64),
        "effective_planck_ratio_min": Q(1, 2),
        "characteristic_imaginary_part_max": Q(1, 10**10),
        "eigenframe_condition_number_max": Q(10**8),
        "symmetrizer_coercivity_min": Q(1, 10**10),
        "normalized_constraint_infinity_max": Q(1, 10**7),
        "minimum_observed_constraint_convergence_order": Q(3, 2),
        "phi_over_Lambda_max": Q(1, 4),
        "proper_frequency_over_Lambda_max": Q(1, 4),
        "proper_wavenumber_over_Lambda_max": Q(1, 4),
        "curvature_scale_squared_over_Lambda_squared_max": Q(1, 4),
        "boundary_margin_over_required_buffer_min": Q(0),
    }
    if values != expected_values:
        raise ValueError("HLT1 thresholds differ")
    thresholds = HealthThresholds(
        newton_residual_max=float(values["newton_residual_infinity_max"]),
        newton_iteration_max=int(values["newton_iteration_max"]),
        kinetic_condition_max=float(values["kinetic_condition_number_max"]),
        branch_displacement_max=float(values["normalized_branch_displacement_max"]),
        acceleration_max=float(values["acceleration_infinity_max"]),
        companion_deformation_max=float(values["companion_deformation_operator_2_max"]),
        action_energy_deformation_max=float(values["action_energy_deformation_operator_2_max"]),
        kinetic_deformation_max=float(values["kinetic_deformation_operator_2_max"]),
        effective_planck_ratio_min=float(values["effective_planck_ratio_min"]),
        characteristic_imaginary_max=float(values["characteristic_imaginary_part_max"]),
        eigenframe_condition_max=float(values["eigenframe_condition_number_max"]),
        symmetrizer_coercivity_min=float(values["symmetrizer_coercivity_min"]),
        normalized_constraint_max=float(values["normalized_constraint_infinity_max"]),
        constraint_convergence_order_min=float(values["minimum_observed_constraint_convergence_order"]),
        phi_over_cutoff_max=float(values["phi_over_Lambda_max"]),
        proper_frequency_over_cutoff_max=float(values["proper_frequency_over_Lambda_max"]),
        proper_wavenumber_over_cutoff_max=float(values["proper_wavenumber_over_Lambda_max"]),
        curvature_scale_squared_over_cutoff_squared_max=float(values["curvature_scale_squared_over_Lambda_squared_max"]),
    )
    return values, thresholds


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    _keys("config", raw, {"schema_version", "artifact_id", "project_version", "metric_signature", "riemann_convention", *EXPECTED_PATHS, "scope", "canonical_norms", "spectral_estimators", "curvature_estimator", "thresholds", "strict_comparison", "transaction", "proof_contract", "claims"})
    if raw["schema_version"] != 1 or raw["artifact_id"] != ARTIFACT_ID or raw["project_version"] != PROJECT_VERSION or raw["metric_signature"] != "-+++" or raw["riemann_convention"] != "plus_partial_mu_gamma_nu":
        raise ValueError("HLT1 identity or convention differs")
    paths = {name: _path(name, raw[name], expected) for name, expected in EXPECTED_PATHS.items()}
    scope = _mapping("scope", raw["scope"])
    if scope != {"target_branch": "FGC-QR", "physical_equations": "unredefined_ACT1_VAR1", "monitor_role": "pre_acceptance_classical_health_and_scale_stop_transaction", "accepted_trajectory_evaluated": False, "holdout_execution_performed": False, "retained_EFT_remainder_control_supplied": False}:
        raise ValueError("HLT1 scope differs")
    norms = _mapping("canonical_norms", raw["canonical_norms"])
    _keys("canonical_norms", norms, {"field_order", "field_mass_dimensions", "state_formula", "matrix_formula", "branch_formula", "constraint_formula", "constraint_floor", "zero_denominator_policy"})
    if tuple(norms["field_order"]) != FIELD_ORDER or tuple(norms["field_mass_dimensions"]) != tuple(int(item) for item in FIELD_MASS_DIMENSIONS) or _fraction("canonical_norms.constraint_floor", norms["constraint_floor"], positive=True) != Q(1, 10**30):
        raise ValueError("HLT1 canonical norm basis differs")
    spectra = _mapping("spectral_estimators", raw["spectral_estimators"])
    _keys("spectral_estimators", spectra, {"minimum_samples", "proper_coordinate_requirement", "spatial_window", "temporal_window", "relative_amplitude_floor", "absolute_amplitude_floor", "unresolved_top_fraction", "support_frequency", "angular_conversion", "alias_rule"})
    spectral_values = {
        "relative": _fraction("spectral relative floor", spectra["relative_amplitude_floor"], positive=True),
        "absolute": _fraction("spectral absolute floor", spectra["absolute_amplitude_floor"], positive=True),
        "top": _fraction("spectral unresolved fraction", spectra["unresolved_top_fraction"], positive=True),
    }
    if spectra["minimum_samples"] != 16 or spectral_values != {"relative": Q(1, 2**40), "absolute": Q(1, 10**30), "top": Q(1, 8)}:
        raise ValueError("HLT1 spectral estimator differs")
    curvature = _mapping("curvature_estimator", raw["curvature_estimator"])
    if curvature != {"linear_scale": "kappa=max_abs_Kretschmann_pow_one_quarter_and_sqrt_abs_Ricci_eigenvalue", "dimensionless_monitor": "kappa_squared_over_Lambda_squared", "protocol_key_interpretation": "curvature_scale_over_Lambda_squared_means_square_of_linear_scale_ratio"}:
        raise ValueError("HLT1 curvature estimator differs")
    values, thresholds = _threshold_values(raw)
    strict = _mapping("strict_comparison", raw["strict_comparison"])
    if strict != {"continuous_maximum_rule": "value_strictly_less_than_limit", "continuous_minimum_rule": "value_strictly_greater_than_limit", "integer_iteration_rule": "value_less_than_or_equal_to_limit", "boolean_premise_rule": "required_boolean_must_be_true", "alias_rule": "occupied_must_be_false"}:
        raise ValueError("HLT1 strict comparison differs")
    transaction = _mapping("transaction", raw["transaction"])
    _keys(
        "transaction",
        transaction,
        {
            "monitor_every_proposed_accepted_Runge_Kutta_stage",
            "stage_index_is_a_global_transaction_serial",
            "stage_index_must_increase_strictly",
            "physical_stage_abscissae_may_be_nonmonotone",
            "evaluate_all_premises_before_state_commit",
            "failed_trial_never_advances_last_accepted_state",
            "first_failed_premise_uses_frozen_priority",
            "first_failed_premise_is_immutable_after_stop",
            "further_acceptance_after_stop_forbidden",
        },
    )
    if not transaction or any(value is not True for value in transaction.values()):
        raise ValueError("every HLT1 transaction flag must be true")
    proof = _mapping("proof_contract", raw["proof_contract"])
    _keys(
        "proof_contract",
        proof,
        {
            "all_PROTO3_health_thresholds_must_be_consumed",
            "SRC1_branch_thresholds_must_be_consumed",
            "HYP2_sufficient_ball_thresholds_must_be_consumed",
            "BND2_causal_margin_must_be_consumed",
            "canonical_norms_must_have_direct_controls",
            "proper_spectral_estimators_must_have_low_and_alias_controls",
            "every_stop_reason_must_be_injected_independently",
            "simultaneous_failures_must_preserve_frozen_priority",
            "equal_continuous_threshold_must_fail",
            "canonical_hash_bound_result_required",
            "no_holdout_execution",
        },
    )
    if not proof or any(value is not True for value in proof.values()):
        raise ValueError("every HLT1 proof flag must be true")
    claims = _mapping("claims", raw["claims"])
    if claims != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("HLT1 claims must remain fail-closed")

    with paths["protocol_config"].open("rb") as handle:
        protocol_raw = tomllib.load(handle)
    protocol = validate_sf1_protocol(protocol_raw)
    protocol_health = protocol_raw["health_stops"]
    protocol_matches = {
        "newton_residual_infinity_max": "newton_residual_infinity_max",
        "newton_iteration_max": "newton_iteration_max",
        "kinetic_condition_number_max": "kinetic_condition_number_max",
        "normalized_branch_displacement_max": "normalized_branch_displacement_max",
        "effective_planck_ratio_min": "effective_planck_ratio_min",
        "characteristic_imaginary_part_max": "characteristic_imaginary_part_max",
        "eigenframe_condition_number_max": "eigenframe_condition_number_max",
        "symmetrizer_coercivity_min": "symmetrizer_coercivity_min",
        "normalized_constraint_infinity_max": "normalized_constraint_infinity_max",
        "minimum_observed_constraint_convergence_order": "minimum_observed_constraint_convergence_order",
        "phi_over_Lambda_max": "phi_over_Lambda_max",
        "proper_frequency_over_Lambda_max": "proper_frequency_over_Lambda_max",
        "proper_wavenumber_over_Lambda_max": "proper_wavenumber_over_Lambda_max",
    }
    for local, source in protocol_matches.items():
        expected = protocol_health[source]
        observed = raw["thresholds"][local]
        if observed != expected:
            raise ValueError(f"HLT1 threshold {local} differs from PROTO3")
    if raw["thresholds"]["curvature_scale_squared_over_Lambda_squared_max"] != protocol_health["curvature_scale_over_Lambda_squared_max"]:
        raise ValueError("HLT1 curvature threshold differs from PROTO3")

    results: dict[str, dict[str, Any]] = {}
    for name, (artifact_id, gate) in RESULT_REQUIREMENTS.items():
        result = _load_json(paths[name], name)
        if result.get("artifact_id") != artifact_id or result.get("gate_status", {}).get(gate) is not True:
            raise ValueError(f"{name} identity or required gate differs")
        results[name] = result
    return {"raw": raw, "paths": paths, "protocol": protocol, "results": results, "threshold_values": values, "thresholds": thresholds, "spectral_values": spectral_values}


def _scope_bindings(loaded: Mapping[str, Any]) -> dict[str, Any]:
    paths = loaded["paths"]
    protocol = loaded["protocol"]
    return {
        "run1_artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "run1_config_sha256": _sha(paths["run1_config"]),
        "protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "protocol_config_sha256": _sha(paths["protocol_config"]),
        "protocol_semantic_holdout_contract_sha256": protocol["semantic_holdout_contract_sha256"],
        "action_artifact_id": "FGC-1-ACT1",
        "action_config_sha256": _sha(paths["action_config"]),
        "action_result_sha256": _sha(paths["action_result"]),
        "variation_artifact_id": "FGC-1-VAR1",
        "variation_config_sha256": _sha(paths["variation_config"]),
        "variation_result_sha256": _sha(paths["variation_result"]),
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "declared_run_envelope_sha256": protocol["semantic_holdout_contract_sha256"],
    }


def _baseline_snapshot(**changes: Any) -> StageHealthSnapshot:
    base = StageHealthSnapshot(
        time=0.125,
        stage_index=1,
        minimum_lapse=1.0,
        minimum_radial_metric=1.0,
        minimum_areal_radius_away_from_center=0.01,
        minimum_hat_lorentzian_margin=0.5,
        newton_residual_infinity=0.0,
        newton_iterations=2,
        newton_residual_monotonic=True,
        kinetic_condition_infinity=100.0,
        normalized_branch_displacement=0.001,
        acceleration_infinity=0.5,
        companion_deformation_operator_2=1.0 / 8192.0,
        action_energy_deformation_operator_2=1.0 / 16384.0,
        kinetic_deformation_operator_2=1.0 / 128.0,
        effective_planck_ratio=1.0,
        characteristic_imaginary_part=0.0,
        eigenframe_condition_number=100.0,
        symmetrizer_coercivity=0.5,
        normalized_constraint_infinity=0.0,
        observed_constraint_convergence_order=2.0,
        phi_over_cutoff=0.01,
        proper_frequency_over_cutoff=0.05,
        proper_wavenumber_over_cutoff=0.05,
        curvature_scale_squared_over_cutoff_squared=0.05,
        temporal_alias_band_occupied=False,
        radial_alias_band_occupied=False,
        boundary_margin_over_required_buffer=1.0,
    )
    return replace(base, **changes)


def _stop_injections(thresholds: HealthThresholds) -> dict[str, Any]:
    mutations = {
        "nonpositive_lapse": {"minimum_lapse": 0.0},
        "nonpositive_radial_metric": {"minimum_radial_metric": 0.0},
        "nonpositive_areal_radius_away_from_center": {"minimum_areal_radius_away_from_center": 0.0},
        "hat_cone_not_Lorentzian": {"minimum_hat_lorentzian_margin": 0.0},
        "newton_residual_limit": {"newton_residual_infinity": thresholds.newton_residual_max},
        "newton_iteration_limit": {"newton_iterations": thresholds.newton_iteration_max + 1},
        "newton_residual_not_monotonic": {"newton_residual_monotonic": False},
        "kinetic_condition_limit": {"kinetic_condition_infinity": thresholds.kinetic_condition_max},
        "branch_displacement_limit": {"normalized_branch_displacement": thresholds.branch_displacement_max},
        "acceleration_limit": {"acceleration_infinity": thresholds.acceleration_max},
        "all_covector_companion_deformation_limit": {"companion_deformation_operator_2": thresholds.companion_deformation_max},
        "action_energy_deformation_limit": {"action_energy_deformation_operator_2": thresholds.action_energy_deformation_max},
        "kinetic_deformation_limit": {"kinetic_deformation_operator_2": thresholds.kinetic_deformation_max},
        "effective_planck_ratio_limit": {"effective_planck_ratio": thresholds.effective_planck_ratio_min},
        "characteristic_imaginary_part_limit": {"characteristic_imaginary_part": thresholds.characteristic_imaginary_max},
        "eigenframe_condition_limit": {"eigenframe_condition_number": thresholds.eigenframe_condition_max},
        "symmetrizer_coercivity_limit": {"symmetrizer_coercivity": thresholds.symmetrizer_coercivity_min},
        "normalized_constraint_limit": {"normalized_constraint_infinity": thresholds.normalized_constraint_max},
        "constraint_convergence_order_limit": {"observed_constraint_convergence_order": thresholds.constraint_convergence_order_min},
        "phi_cutoff_proxy_limit": {"phi_over_cutoff": thresholds.phi_over_cutoff_max},
        "proper_frequency_cutoff_proxy_limit": {"proper_frequency_over_cutoff": thresholds.proper_frequency_over_cutoff_max},
        "proper_wavenumber_cutoff_proxy_limit": {"proper_wavenumber_over_cutoff": thresholds.proper_wavenumber_over_cutoff_max},
        "temporal_alias_band_occupied": {"temporal_alias_band_occupied": True},
        "radial_alias_band_occupied": {"radial_alias_band_occupied": True},
        "curvature_cutoff_proxy_limit": {"curvature_scale_squared_over_cutoff_squared": thresholds.curvature_scale_squared_over_cutoff_squared_max},
        "boundary_causal_buffer": {"boundary_margin_over_required_buffer": 0.0},
    }
    if tuple(mutations) != STOP_PRIORITY:
        raise ValueError("HLT1 injected stop order differs from implementation")
    records = {}
    for reason, mutation in mutations.items():
        monitor = ClassicalHealthMonitor(thresholds)
        try:
            monitor.assess_and_accept(_baseline_snapshot(**mutation))
        except ClassicalHealthStop as stop:
            if stop.reason != reason:
                raise ValueError(f"HLT1 injection {reason} stopped as {stop.reason}")
            records[reason] = {
                "typed_reason": stop.reason,
                "all_failures": stop.failures,
                "accepted_stage_count": stop.state.accepted_stage_count,
                "last_accepted_time": stop.state.last_accepted_time,
                "passed": True,
            }
        else:
            raise ValueError(f"HLT1 injection {reason} was accepted")

    transaction = ClassicalHealthMonitor(thresholds)
    accepted = transaction.assess_and_accept(_baseline_snapshot(time=0.125, stage_index=1))
    try:
        transaction.assess_and_accept(_baseline_snapshot(time=0.25, stage_index=2, minimum_lapse=0.0, newton_residual_infinity=thresholds.newton_residual_max))
    except ClassicalHealthStop as first:
        first_state = first.state
    else:
        raise ValueError("HLT1 simultaneous injected failure was accepted")
    try:
        transaction.assess_and_accept(_baseline_snapshot(time=0.375, stage_index=3))
    except ClassicalHealthStop as repeated:
        repeated_state = repeated.state
    else:
        raise ValueError("HLT1 accepted a stage after its first stop")

    ssprk_transaction = ClassicalHealthMonitor(thresholds)
    ssprk_states = tuple(
        ssprk_transaction.assess_and_accept(
            _baseline_snapshot(time=stage_time, stage_index=serial)
        )
        for serial, stage_time in enumerate((0.0, 1.0, 0.5, 1.0), start=1)
    )
    before_duplicate = ssprk_transaction.state
    try:
        ssprk_transaction.assess_and_accept(
            _baseline_snapshot(time=2.0, stage_index=4)
        )
    except ValueError as error:
        duplicate_rejected = "strictly increasing global" in str(error)
    else:
        duplicate_rejected = False
    if not duplicate_rejected or ssprk_transaction.state != before_duplicate:
        raise ValueError("HLT1 global transaction ordering control failed")
    return {
        "records": records,
        "injected_reason_order": STOP_PRIORITY,
        "every_reason_independently_stopped": len(records) == len(STOP_PRIORITY),
        "transaction_control": {
            "accepted_before_failure": accepted.accepted_stage_count,
            "simultaneous_first_reason": first_state.first_failed_premise,
            "last_accepted_stage_after_failure": first_state.last_accepted_stage_index,
            "last_accepted_time_after_failure": first_state.last_accepted_time,
            "first_failure_immutable_on_later_attempt": repeated_state == first_state,
            "further_acceptance_after_stop_forbidden": True,
            "SSPRK3_nonmonotone_stage_abscissae_supported": len(ssprk_states) == 4,
            "strict_global_stage_transaction_order_enforced": duplicate_rejected,
            "ordering_error_does_not_mutate_last_accepted_state": ssprk_transaction.state == before_duplicate,
        },
    }


def _numerical_formula_controls(loaded: Mapping[str, Any]) -> dict[str, Any]:
    u = np.asarray((2.0, 3.0, 4.0, 8.0, 5.0, 6.0))
    p = np.ones(6)
    q = 2.0 * np.ones(6)
    du, dp, dq = dimensionless_state_components(u, p, q, length_unit=4.0)
    matrix_norm = induced_scaled_infinity_norm(np.asarray(((1.0, -2.0), (3.0, 4.0))), row_scales=(2.0, 0.5), column_scales=(1.0, 2.0))
    constraints = normalized_componentwise_infinity((0.0, 2.0e-8), (0.0, 2.0))
    branch = normalized_branch_displacement((1.0, 3.0), (0.0, 1.0), component_scales=(2.0, 4.0))
    kappa, curvature_ratio = curvature_scale_and_ratio((256.0,), (4.0,), cutoff=16.0)

    sample_count = loaded["raw"]["spectral_estimators"]["minimum_samples"] * 8
    proper = np.arange(sample_count, dtype=np.float64) / sample_count
    low = np.sin(2.0 * pi * 2.0 * proper)
    high = np.sin(2.0 * pi * 60.0 * proper)
    spectral_kwargs = {
        "relative_amplitude_floor": float(loaded["spectral_values"]["relative"]),
        "absolute_amplitude_floor": float(loaded["spectral_values"]["absolute"]),
        "unresolved_top_fraction": float(loaded["spectral_values"]["top"]),
    }
    low_support = windowed_spectral_support(low, proper, **spectral_kwargs)
    high_support = windowed_spectral_support(high, proper, **spectral_kwargs)
    window = compact_vacuum_buffer_window(proper, support_minimum=0.0, support_maximum=(sample_count - 1) / sample_count, taper_width=0.125)
    if low_support.support_bin != 2 or low_support.alias_band_occupied or high_support.support_bin != 60 or not high_support.alias_band_occupied:
        raise ValueError("HLT1 spectral support controls failed")
    return {
        "dimensionless_state": {"u": du, "p": dp, "q": dq, "field_order": FIELD_ORDER, "mass_dimensions": FIELD_MASS_DIMENSIONS},
        "induced_scaled_infinity_norm": matrix_norm,
        "normalized_constraint_infinity": constraints,
        "normalized_branch_displacement": branch,
        "curvature": {"linear_scale": kappa, "scale_squared_over_cutoff_squared": curvature_ratio},
        "spectral": {"sample_count": sample_count, "low_mode": low_support, "injected_alias_mode": high_support, "compact_window_exact_zero_at_both_edges": bool(window[0] == 0.0 and window[-1] == 0.0), "compact_window_has_unit_plateau": bool(np.max(window) == 1.0)},
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    loaded = load_config(config_path)
    paths = loaded["paths"]
    scope = _scope_bindings(loaded)
    formula_controls = _numerical_formula_controls(loaded)
    injections = _stop_injections(loaded["thresholds"])
    predecessor_names = tuple(name for name in EXPECTED_PATHS if name not in {"protocol_config", "run1_config"})
    quantitative = {
        "canonical_formulas": {
            "state": "u_i*L0^d_i; p_i,q_i*L0^(d_i+1)",
            "matrix": "max_row_sum(abs(diag(row_scale)*A*diag(column_scale)^-1))",
            "branch": "max_i(abs(a_i-a_previous_accepted_i)/scale_i)",
            "constraint": "max_i(abs(C_i)/max(abs(source_i),1e-30))",
            "curvature": "(max(abs(K)^1/4,sqrt(abs(Ricci_eigenvalue)))/Lambda)^2",
            "spectral": "largest_windowed_rfft_bin_above_max(absolute_floor,relative_floor*peak)",
        },
        "formula_controls": formula_controls,
        "thresholds": loaded["threshold_values"],
        "strict_comparison_contract": {
            "continuous_upper_equality_stops": True,
            "continuous_lower_equality_stops": True,
            "Newton_iteration_count_equal_to_16_is_allowed": True,
            "alias_occupation_stops": True,
        },
        "typed_stop_injections": injections,
        "source_predecessor_margins": {
            "HYP2_witness_extrema": loaded["results"]["multidirectional_health_result"]["artifact_payload"]["quantitative_evidence"]["witness_extrema"],
            "BND2_smallest_reference_margin": loaded["results"]["boundary_control_result"]["artifact_payload"]["quantitative_evidence"]["outer_boundary_move_control"]["smallest_reference_margin_over_required_buffer"],
            "CON4_absolute_threshold_was_deferred_to_HLT1_NUM1": loaded["results"]["physical_constraints_result"]["artifact_payload"]["quantitative_evidence"]["monitor_contract"]["HLT1_and_NUM1_must_add_absolute_scale_and_thresholds"],
        },
        "epistemic_boundary": {
            "monitor_formulas_verified_without_an_accepted_FGCQR_trajectory": True,
            "future_solver_must_call_transaction_before_every_stage_commit": True,
            "spectral_window_is_not_an_exact_continuum_bandlimit": True,
            "scale_thresholds_are_classical_stop_proxies_not_EFT_remainder_control": True,
            "time_evolution_performed": False,
            "FGCQR_holdout_outcome_inspected": False,
        },
    }
    return _serial({
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": FUTURE_EVIDENCE_CLASSIFICATIONS["health_monitor"],
        "generated_by": "scripts/reproduce_fgc_hlt1_mon1.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "predecessor_sha256": {_rel(paths[name]): _sha(paths[name]) for name in predecessor_names},
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": scope,
        "certificate_contract": {"artifact_specific_payload_validated": True, "canonical_reproduction_passed": True, "scope_bindings_verified": True, "outcome_data_not_used_to_select_certificate_contract": True},
        "artifact_payload": {
            "run_envelope_semantic_sha256": scope["declared_run_envelope_sha256"],
            "canonical_dimensionless_norms_executable": True,
            "all_protocol_thresholds_implemented_strictly": True,
            "every_accepted_stage_monitored": True,
            "nonmonotone_internal_stage_abscissae_supported": True,
            "strict_global_stage_transaction_order_enforced": True,
            "stop_occurs_before_threshold_crossing": True,
            "first_failed_premise_is_preserved": True,
            "scale_stops_not_promoted_to_EFT_validity": True,
            "quantitative_evidence": quantitative,
        },
        "gate_status": {REQUIRED_GATE: True},
        "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
    })


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    value = _load_json(path, "HLT1 result")
    if path.read_text(encoding="utf-8") != _canonical(value):
        raise ValueError("HLT1 result must use canonical sorted JSON")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    arguments.output.write_text(_canonical(record(arguments.config)), encoding="utf-8")


if __name__ == "__main__":
    main()
