#!/usr/bin/env python3
"""Regenerate the outcome-blind FGC-1-NUM1-VAL1 validation certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite, log
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-num1-val1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-num1-val1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-num1-val1.md"

from scripts.reproduce_fgc_hlt1_mon1 import (  # noqa: E402
    load_config as load_hlt1_config,
)
from scripts.reproduce_fgc_src1_nl1 import (  # noqa: E402
    load_config as load_src1_config,
)
from recursive_horizons.fgc.action import ActionParameters, ModelID  # noqa: E402
from recursive_horizons.fgc.evolution.affine_measurement import (  # noqa: E402
    analyze_affine_null,
    polynomial_flrw_history,
)
from recursive_horizons.fgc.evolution.nonlinear_source import (  # noqa: E402
    solve_accelerations,
)
from recursive_horizons.fgc.evolution.numerical_controls import (  # noqa: E402
    WaveControl,
    boundary_isolation_control,
    checkpoint_restart_control,
    convergence_control,
    minkowski_preservation_control,
    runtime_transaction_control,
    sbp_and_dissipation_controls,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
    adm_accelerations_from_base_metric,
    adm_lower_jets_from_base_state,
    regular_center_acceleration_limit,
    solve_grid_accelerations,
)
from recursive_horizons.fgc.evolution.health_monitor import STOP_PRIORITY  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.modified_harmonic_modes import (  # noqa: E402
    exact_comp1_acceleration_root,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_COMMON_NONCLAIMS,
    FUTURE_EVIDENCE_CLASSIFICATIONS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    validate_sf1_protocol,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-NUM1-VAL1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "independent_solver_validation_and_measurement_contract_passed"
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
    "regular_center_config": "configs/fgc/fgc-1-ctr1-reg1.toml",
    "regular_center_result": "results/fgc-1-ctr1-reg1.json",
    "initial_data_config": "configs/fgc/fgc-1-id1-fam1.toml",
    "initial_data_result": "results/fgc-1-id1-fam1.json",
    "run_domain_config": "configs/fgc/fgc-1-dom4-run1.toml",
    "run_domain_result": "results/fgc-1-dom4-run1.json",
    "multidirectional_health_config": "configs/fgc/fgc-1-hyp2-md1.toml",
    "multidirectional_health_result": "results/fgc-1-hyp2-md1.json",
    "boundary_control_config": "configs/fgc/fgc-1-bnd2-cp1.toml",
    "boundary_control_result": "results/fgc-1-bnd2-cp1.json",
    "health_monitor_config": "configs/fgc/fgc-1-hlt1-mon1.toml",
    "health_monitor_result": "results/fgc-1-hlt1-mon1.json",
    "observable_config": "configs/fgc/fgc-1-def0-obs1.toml",
    "observable_result": "results/fgc-1-def0-obs1.json",
}
RESULT_REQUIREMENTS = {
    "action_result": ("FGC-1-ACT1", "declared_action_and_scalar_algebra_certificate_passed"),
    "variation_result": ("FGC-1-VAR1", "var1_covariant_metric_equation_derived_and_cross_checked"),
    "nonlinear_source_result": ("FGC-1-SRC1-NL1", "nonlinear_REF1_source_and_branch_solver_verified"),
    "physical_constraints_result": ("FGC-1-CON4-PHY1", "physical_gauge_reduction_constraint_system_closed"),
    "regular_center_result": ("FGC-1-CTR1-REG1", "regular_center_formulation_verified"),
    "initial_data_result": ("FGC-1-ID1-FAM1", "nonzero_width_finite_mass_constraint_compatible_family_constructed"),
    "run_domain_result": ("FGC-1-DOM4-RUN1", "nonzero_classical_spherical_run_envelope_passed"),
    "multidirectional_health_result": ("FGC-1-HYP2-MD1", "quantitative_all_covector_weak_coupling_health_envelope_passed"),
    "boundary_control_result": ("FGC-1-BND2-CP1", "spherical_boundary_or_domain_of_dependence_control_passed"),
    "health_monitor_result": ("FGC-1-HLT1-MON1", "classical_health_and_typed_stop_monitoring_verified"),
    "observable_result": ("FGC-1-DEF0-OBS1", "exact_metric_null_observable_contract_and_controls_passed"),
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/evolution/numerical_engine.py",
        "src/recursive_horizons/fgc/evolution/numerical_controls.py",
        "src/recursive_horizons/fgc/evolution/runtime_transaction.py",
        "src/recursive_horizons/fgc/evolution/affine_measurement.py",
        "src/recursive_horizons/fgc/evolution/vectorized_source.py",
        "src/recursive_horizons/fgc/evolution/nonlinear_source.py",
        "src/recursive_horizons/fgc/modified_harmonic_constraints.py",
        "src/recursive_horizons/fgc/modified_harmonic_modes.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "scripts/reproduce_fgc_hlt1_mon1.py",
        "scripts/reproduce_fgc_src1_nl1.py",
        "scripts/reproduce_fgc_num1_val1.py",
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


def _fraction(
    name: str,
    value: Any,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> Fraction:
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
    try:
        relative = path.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} escapes the repository") from exc
    if relative != expected or not path.is_file():
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
    if isinstance(value, (tuple, list)):
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


def _parse_thresholds(raw: Mapping[str, Any]) -> dict[str, Fraction]:
    table = _mapping("thresholds", raw["thresholds"])
    expected = {
        "SBP_identity_residual_max",
        "Minkowski_binary64_drift_max",
        "primary_minimum_solution_order",
        "comparator_minimum_solution_order",
        "minimum_reduction_constraint_order",
        "exact_scalar_batch_root_absolute_error_max",
        "batch_REF1_residual_infinity_max",
        "center_limit_minimum_observed_order",
        "boundary_measurement_difference_max",
        "affine_normalization_residual_max",
        "affine_geodesic_residual_max",
        "null_residual_max",
        "direct_Raychaudhuri_agreement_residual_max",
        "cross_method_theta_difference_max",
    }
    _keys("thresholds", table, expected)
    values = {
        key: _fraction(
            f"thresholds.{key}",
            value,
            positive=key != "Minkowski_binary64_drift_max",
            nonnegative=key == "Minkowski_binary64_drift_max",
        )
        for key, value in table.items()
    }
    required = {
        "SBP_identity_residual_max": Q(1, 10**14),
        "Minkowski_binary64_drift_max": Q(0),
        "primary_minimum_solution_order": Q(5, 2),
        "comparator_minimum_solution_order": Q(19, 10),
        "minimum_reduction_constraint_order": Q(3, 2),
        "exact_scalar_batch_root_absolute_error_max": Q(1, 10**12),
        "batch_REF1_residual_infinity_max": Q(1, 10**12),
        "center_limit_minimum_observed_order": Q(5),
        "boundary_measurement_difference_max": Q(1, 10**13),
        "affine_normalization_residual_max": Q(1, 10**10),
        "affine_geodesic_residual_max": Q(1, 10**10),
        "null_residual_max": Q(1, 10**12),
        "direct_Raychaudhuri_agreement_residual_max": Q(1, 10**8),
        "cross_method_theta_difference_max": Q(1, 10**8),
    }
    if values != required:
        raise ValueError("NUM1 thresholds differ")
    return values


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    _keys(
        "config",
        raw,
        {
            "schema_version", "artifact_id", "project_version", "metric_signature",
            "riemann_convention", *EXPECTED_PATHS, "scope", "methods", "controls",
            "thresholds", "proof_contract", "claims",
        },
    )
    if raw["schema_version"] != 1 or raw["artifact_id"] != ARTIFACT_ID:
        raise ValueError("NUM1 identity differs")
    if raw["project_version"] != PROJECT_VERSION or raw["metric_signature"] != "-+++" or raw["riemann_convention"] != "plus_partial_mu_gamma_nu":
        raise ValueError("NUM1 version or convention differs")
    paths = {name: _path(name, raw[name], expected) for name, expected in EXPECTED_PATHS.items()}
    with paths["protocol_config"].open("rb") as handle:
        protocol_raw = tomllib.load(handle)
    protocol = validate_sf1_protocol(protocol_raw)

    scope = _mapping("scope", raw["scope"])
    expected_scope = {
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "validation_role": "outcome_blind_numerical_engine_source_adapter_and_affine_measurement_controls",
        "FGCQR_trajectory_evaluated": False,
        "collapse_evolution_performed": False,
        "holdout_outcome_inspected": False,
        "retained_EFT_remainder_control_supplied": False,
    }
    if scope != expected_scope:
        raise ValueError("NUM1 scope differs")

    methods = _mapping("methods", raw["methods"])
    if methods != {
        "primary": PRIMARY_METHOD,
        "comparator": COMPARATOR_METHOD,
        "array_backend": "numpy_binary64",
        "only_new_numerical_runtime_dependency": "numpy==2.5.1",
    } or np.__version__ != "2.5.1":
        raise ValueError("NUM1 methods or pinned NumPy runtime differ")
    controls = _mapping("controls", raw["controls"])
    _keys(
        "controls", controls,
        {
            "convergence_point_counts", "convergence_final_time", "wave_domain_minimum",
            "wave_domain_maximum", "Minkowski_point_count", "Minkowski_step_count",
            "checkpoint_total_steps", "boundary_aligned_spacing", "boundary_outer_radii",
            "boundary_measurement_radius", "affine_history_time_count",
            "affine_history_radius_count", "affine_final_time", "affine_maximum_radius",
            "affine_launch_radius", "affine_step_size", "center_limit_spacings",
            "typed_stop_count",
        },
    )
    if controls["convergence_point_counts"] != [33, 65, 129] or controls["Minkowski_point_count"] != 65 or controls["Minkowski_step_count"] != 16 or controls["checkpoint_total_steps"] != 16 or controls["affine_history_time_count"] != 65 or controls["affine_history_radius_count"] != 129 or controls["typed_stop_count"] != len(STOP_PRIORITY):
        raise ValueError("NUM1 integer controls differ")
    fractions = {
        key: _fraction(
            f"controls.{key}", controls[key],
            positive=key != "wave_domain_minimum",
            nonnegative=key == "wave_domain_minimum",
        )
        for key in (
            "convergence_final_time", "wave_domain_minimum", "wave_domain_maximum",
            "boundary_aligned_spacing", "boundary_measurement_radius", "affine_final_time",
            "affine_maximum_radius", "affine_launch_radius", "affine_step_size",
        )
    }
    outer_radii = tuple(_fraction("controls.boundary_outer_radii", value, positive=True) for value in controls["boundary_outer_radii"])
    center_spacings = tuple(_fraction("controls.center_limit_spacings", value, positive=True) for value in controls["center_limit_spacings"])
    if fractions != {
        "convergence_final_time": Q(1, 5), "wave_domain_minimum": Q(0),
        "wave_domain_maximum": Q(3141592653589793, 10**15),
        "boundary_aligned_spacing": Q(1, 32), "boundary_measurement_radius": Q(5, 2),
        "affine_final_time": Q(1, 4), "affine_maximum_radius": Q(4),
        "affine_launch_radius": Q(1), "affine_step_size": Q(1, 2048),
    } or outer_radii != (Q(6), Q(8)) or center_spacings != (Q(1, 8), Q(1, 16), Q(1, 32)):
        raise ValueError("NUM1 rational controls differ")

    proof = _mapping("proof_contract", raw["proof_contract"])
    expected_proof = {
        "both_frozen_methods_required", "candidate_endpoint_evaluated_separately",
        "regular_center_parity_and_limit_controls_required", "same_unredefined_REF1_evaluator_required",
        "exact_scalar_and_batched_ADM_source_routes_required",
        "all_HLT1_typed_stops_and_atomic_runtime_composition_required",
        "checkpoint_restart_must_be_bitwise", "outer_boundary_isolation_must_use_aligned_spacing",
        "physical_metric_null_and_affine_residuals_required",
        "direct_and_complete_Raychaudhuri_routes_required", "canonical_hash_bound_result_required",
        "no_FGCQR_holdout_execution",
    }
    _keys("proof_contract", proof, expected_proof)
    if not all(value is True for value in proof.values()):
        raise ValueError("NUM1 proof contract must be all true")
    if raw["claims"] != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("NUM1 nonclaim ledger differs")

    protocol_methods = protocol_raw["numerics"]
    if (
        protocol_methods["primary_spatial_method"] + "_plus_" + protocol_methods["primary_time_method"] != PRIMARY_METHOD
        or protocol_methods["comparator_spatial_method"] + "_plus_" + protocol_methods["comparator_time_method"] != COMPARATOR_METHOD
    ):
        raise ValueError("NUM1 methods differ from PROTO3")
    protocol_analysis = protocol_raw["analysis"]
    if _fraction("PROTO3 affine residual", protocol_analysis["affine_normalization_residual_max"], positive=True) != _parse_thresholds(raw)["affine_normalization_residual_max"] or _fraction("PROTO3 route residual", protocol_analysis["direct_Raychaudhuri_agreement_residual_max"], positive=True) != _parse_thresholds(raw)["direct_Raychaudhuri_agreement_residual_max"]:
        raise ValueError("NUM1 affine thresholds differ from PROTO3")

    results = {}
    for name, (artifact_id, gate) in RESULT_REQUIREMENTS.items():
        result = _load_json(paths[name], name)
        if result.get("artifact_id") != artifact_id or result.get("gate_status", {}).get(gate) is not True:
            raise ValueError(f"{name} identity or required gate differs")
        results[name] = result
    return {
        "raw": raw,
        "paths": paths,
        "protocol": protocol,
        "protocol_raw": protocol_raw,
        "scope": scope,
        "controls": controls,
        "control_fractions": fractions,
        "center_spacings": center_spacings,
        "thresholds": _parse_thresholds(raw),
        "results": results,
    }


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


def _run_summary(run: Any) -> dict[str, Any]:
    return {
        "point_count": run.point_count,
        "spacing": run.spacing,
        "step_size": run.step_size,
        "step_count": run.step_count,
        "solution_l2_error": run.solution_l2_error,
        "solution_infinity_error": run.solution_infinity_error,
        "reduction_constraint_l2": run.reduction_constraint_l2,
        "reduction_constraint_infinity": run.reduction_constraint_infinity,
        "relative_energy_drift": run.relative_energy_drift,
        "endpoint_sha256": run.endpoint_sha256,
    }


def _convergence_summary(control: WaveControl, method: str, loaded: Mapping[str, Any]) -> dict[str, Any]:
    result = convergence_control(
        control,
        method=method,
        point_counts=loaded["controls"]["convergence_point_counts"],
        final_time=float(loaded["control_fractions"]["convergence_final_time"]),
    )
    return {
        "control_id": result["control_id"],
        "method": result["method"],
        "runs": [_run_summary(run) for run in result["runs"]],
        "solution_orders": result["solution_orders"],
        "constraint_orders": result["constraint_orders"],
        "minimum_solution_order": result["minimum_solution_order"],
        "minimum_constraint_order": result["minimum_constraint_order"],
    }


def _center_limit_control(loaded: Mapping[str, Any]) -> dict[str, Any]:
    exact = np.arange(1.0, 7.0)
    exact[ADM_CENTER_PARITIES == -1] = 0.0
    records = []
    errors = []
    for spacing in loaded["center_spacings"]:
        h = float(spacing)
        radii = h * np.arange(1.0, 5.0)
        values = np.empty((4, 6))
        for field in range(6):
            values[:, field] = (
                field + 1.0
                + (field + 2.0) * radii**2
                + (field + 3.0) * radii**4
                + (field + 4.0) * radii**6
                + (field + 5.0) * radii**8
            )
        limit = regular_center_acceleration_limit(values)
        error = float(np.max(np.abs(limit.acceleration - exact), initial=0.0))
        errors.append(error)
        records.append(
            {
                "spacing": spacing,
                "four_point_error_infinity": error,
                "three_vs_four_estimator_infinity": limit.estimator_infinity,
                "odd_accelerations_exactly_zero": bool(
                    np.array_equal(limit.acceleration[ADM_CENTER_PARITIES == -1], np.zeros(2))
                ),
            }
        )
    orders = tuple(log(coarse / fine, 2.0) for coarse, fine in zip(errors, errors[1:]))
    return {
        "analytic_family": "even_polynomial_through_r_pow_8_with_exact_center_values",
        "records": records,
        "observed_orders": orders,
        "minimum_observed_order": min(orders),
    }


def _source_backend_control(loaded: Mapping[str, Any]) -> dict[str, Any]:
    src = load_src1_config(loaded["paths"]["nonlinear_source_config"])
    qift = src["qift"]
    parameters = qift.flat_fixture["action_parameters"]
    action = ActionParameters(
        model_id=ModelID.FGC_QR,
        planck_mass=float(parameters["planck_mass"]),
        scalar_mass=float(parameters["scalar_mass"]),
        quartic_coupling=float(parameters["quartic_coupling"]),
        pulse_width=1.0,
        ricci_coupling=float(parameters["ricci_coupling"]),
        quadratic_gb_coupling=float(parameters["quadratic_gb_coupling"]),
    )
    flat = state_from_generalized_adm_pg_fixture(qift.flat_fixture)
    activated = activated_compatible_state(qift.flat_fixture)["state"]
    reference = flat_spherical_annulus_reference(
        radial_domain_minimum=qift.radial_domain_minimum
    )
    arguments = {
        "reference": reference,
        "coordinate_radius": qift.coordinate_radius,
        "tilde_normal_factor": qift.tilde_normal_factor,
        "hat_normal_factor": qift.hat_normal_factor,
    }
    u, p, q, p_r, q_r = adm_lower_jets_from_base_state(activated)
    radius = float(qift.coordinate_radius)
    batch = solve_grid_accelerations(
        u[None, :], p[None, :], q[None, :], p_r[None, :], q_r[None, :],
        np.asarray((radius,)), action=action,
    )
    scalar = solve_accelerations(
        activated, branch_center=flat, config=src["solver"], **arguments
    )
    exact_base = exact_comp1_acceleration_root(
        activated,
        **arguments,
        qift_acceleration_half_width=qift.acceleration_half_width,
    )["root_acceleration"]
    scalar_adm = adm_accelerations_from_base_metric(u, p, scalar.accelerations)
    exact_adm = adm_accelerations_from_base_metric(
        u, p, tuple(float(value) for value in exact_base)
    )
    scalar_error = float(np.max(np.abs(batch.accelerations[0] - scalar_adm), initial=0.0))
    exact_error = float(np.max(np.abs(batch.accelerations[0] - exact_adm), initial=0.0))

    radii = np.linspace(0.5, 4.0, 12)
    flat_u = np.zeros((radii.size, 6))
    flat_u[:, 0] = 1.0
    flat_u[:, 2] = 1.0
    flat_u[:, 3] = radii
    flat_p = np.zeros_like(flat_u)
    flat_q = np.zeros_like(flat_u)
    flat_q[:, 3] = 1.0
    flat_grid = solve_grid_accelerations(
        flat_u, flat_p, flat_q, np.zeros_like(flat_u), np.zeros_like(flat_u),
        radii, action=action,
    )
    return {
        "activated_fixture_id": "COMP1-activated-compatible-point",
        "scalar_vs_batch_ADM_root_absolute_error": scalar_error,
        "exact_vs_batch_ADM_root_absolute_error": exact_error,
        "batch_REF1_residual_infinity": batch.residual_infinity,
        "batch_kinetic_condition_infinity": batch.condition_infinity_maximum,
        "batch_branch_displacement_infinity": batch.branch_displacement_infinity,
        "batch_affine_in_accelerations_verified": batch.affine_in_accelerations_verified,
        "flat_annulus_point_count": radii.size,
        "flat_zero_root_residual_infinity": flat_grid.residual_infinity,
        "flat_zero_root_acceleration_infinity": float(np.max(np.abs(flat_grid.accelerations), initial=0.0)),
        "flat_zero_root_iterations": flat_grid.iterations,
    }


def _affine_controls(loaded: Mapping[str, Any]) -> dict[str, Any]:
    controls = loaded["controls"]
    fractions = loaded["control_fractions"]
    history = polynomial_flrw_history(
        time_count=controls["affine_history_time_count"],
        radius_count=controls["affine_history_radius_count"],
        final_time=float(fractions["affine_final_time"]),
        maximum_radius=float(fractions["affine_maximum_radius"]),
    )
    records = []
    raw = []
    for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
        result = analyze_affine_null(
            history,
            method=method,
            launch_time=0.0,
            launch_radius=float(fractions["affine_launch_radius"]),
            final_time=float(fractions["affine_final_time"]),
            step_size=float(fractions["affine_step_size"]),
        )
        raw.append(result)
        records.append(
            {
                "method": method,
                "sample_count": result.theta.size,
                "initial_normalization_residual": result.initial_normalization_residual,
                "maximum_null_residual": result.maximum_null_residual,
                "maximum_affine_geodesic_residual": result.maximum_affine_residual,
                "maximum_direct_Raychaudhuri_disagreement": result.maximum_route_disagreement,
                "trajectory_sha256": array_content_sha256(
                    result.coordinate_times, result.affine_parameters,
                    result.coordinate_radii, result.k_t, result.k_r, result.theta,
                ),
                "measurement_sha256": array_content_sha256(
                    result.direct_dtheta_dlambda, result.complete_raychaudhuri_rhs,
                    result.ricci_null, result.null_residual,
                    result.affine_residual_t, result.affine_residual_r,
                ),
                "radial_screen_shear_squared": 0.0,
                "hypersurface_twist_squared": 0.0,
            }
        )
    cross = float(np.max(np.abs(raw[0].theta - raw[1].theta), initial=0.0))
    return {
        "geometry": "exact_flat_FLRW_with_scale_factor_a_of_t_equals_1_plus_t",
        "records": records,
        "cross_method_theta_difference_infinity": cross,
    }


def _quantitative_certificate(loaded: Mapping[str, Any]) -> dict[str, Any]:
    thresholds = {key: float(value) for key, value in loaded["thresholds"].items()}
    methods = (PRIMARY_METHOD, COMPARATOR_METHOD)
    manufactured = WaveControl("MMS-CENTER", 1, 1.0, True)
    gr0 = WaveControl("GR0-SCALAR", -1, 1.0, False)
    method_records = []
    for method in methods:
        minkowski = minkowski_preservation_control(
            method=method,
            point_count=loaded["controls"]["Minkowski_point_count"],
            step_count=loaded["controls"]["Minkowski_step_count"],
        )
        mms = _convergence_summary(manufactured, method, loaded)
        scalar = _convergence_summary(gr0, method, loaded)
        restart = checkpoint_restart_control(
            method=method,
            total_steps=loaded["controls"]["checkpoint_total_steps"],
        )
        boundary = boundary_isolation_control(method=method)
        solution_minimum = (
            thresholds["primary_minimum_solution_order"]
            if method == PRIMARY_METHOD
            else thresholds["comparator_minimum_solution_order"]
        )
        passed = (
            minkowski["bitwise_preserved"]
            and minkowski["maximum_binary64_drift"] <= thresholds["Minkowski_binary64_drift_max"]
            and mms["minimum_solution_order"] >= solution_minimum
            and scalar["minimum_solution_order"] >= solution_minimum
            and mms["minimum_constraint_order"] > thresholds["minimum_reduction_constraint_order"]
            and scalar["minimum_constraint_order"] > thresholds["minimum_reduction_constraint_order"]
            and restart["canonical_roundtrip_bitwise"]
            and restart["continuous_and_restart_endpoint_bitwise_equal"]
            and boundary["maximum_measurement_difference"] <= thresholds["boundary_measurement_difference_max"]
        )
        if not passed:
            raise ValueError(f"NUM1 method controls failed for {method}")
        method_records.append(
            {
                "method": method,
                "Minkowski": minkowski,
                "regular_center_manufactured_solution": mms,
                "GR0_scalar": scalar,
                "checkpoint_restart": restart,
                "boundary_isolation": boundary,
                "all_method_controls_passed": True,
            }
        )

    sbp = sbp_and_dissipation_controls()
    if any(
        record["sbp_identity_residual"] > thresholds["SBP_identity_residual_max"]
        or not record["dissipation_nonpositive"]
        for record in sbp["records"]
    ):
        raise ValueError("NUM1 SBP or dissipation control failed")
    center = _center_limit_control(loaded)
    if center["minimum_observed_order"] < thresholds["center_limit_minimum_observed_order"] or not all(record["odd_accelerations_exactly_zero"] for record in center["records"]):
        raise ValueError("NUM1 centre acceleration limit control failed")
    source = _source_backend_control(loaded)
    if max(source["scalar_vs_batch_ADM_root_absolute_error"], source["exact_vs_batch_ADM_root_absolute_error"]) > thresholds["exact_scalar_batch_root_absolute_error_max"] or max(source["batch_REF1_residual_infinity"], source["flat_zero_root_residual_infinity"]) > thresholds["batch_REF1_residual_infinity_max"] or source["flat_zero_root_acceleration_infinity"] > thresholds["batch_REF1_residual_infinity_max"] or not source["batch_affine_in_accelerations_verified"]:
        raise ValueError("NUM1 exact/scalar/batched source control failed")

    hlt1 = loaded["results"]["health_monitor_result"]
    injections = hlt1["artifact_payload"]["quantitative_evidence"]["typed_stop_injections"]
    if (
        len(injections["records"]) != loaded["controls"]["typed_stop_count"]
        or set(injections["records"]) != set(STOP_PRIORITY)
        or not injections["every_reason_independently_stopped"]
    ):
        raise ValueError("NUM1 did not consume all HLT1 typed stops")
    runtime = runtime_transaction_control(load_hlt1_config(loaded["paths"]["health_monitor_config"])["thresholds"])
    if not runtime["all_transactions_passed"]:
        raise ValueError("NUM1 runtime transaction control failed")

    affine = _affine_controls(loaded)
    for record in affine["records"]:
        if abs(record["initial_normalization_residual"]) > thresholds["affine_normalization_residual_max"] or record["maximum_null_residual"] > thresholds["null_residual_max"] or record["maximum_affine_geodesic_residual"] > thresholds["affine_geodesic_residual_max"] or record["maximum_direct_Raychaudhuri_disagreement"] > thresholds["direct_Raychaudhuri_agreement_residual_max"]:
            raise ValueError("NUM1 affine measurement control failed")
    if affine["cross_method_theta_difference_infinity"] > thresholds["cross_method_theta_difference_max"]:
        raise ValueError("NUM1 affine comparator disagreement exceeded threshold")

    return {
        "runtime": {
            "python_implementation": sys.implementation.name,
            "python_version": ".".join(str(value) for value in sys.version_info[:3]),
            "numpy_version": np.__version__,
            "floating_type": "IEEE_754_binary64",
        },
        "thresholds": loaded["thresholds"],
        "SBP_and_dissipation": sbp,
        "method_controls": method_records,
        "regular_center_acceleration_limit": center,
        "exact_scalar_and_batched_source": source,
        "typed_stop_count": len(injections["records"]),
        "typed_stop_priority_order": STOP_PRIORITY,
        "runtime_transaction": runtime,
        "affine_measurement": affine,
        "epistemic_boundary": {
            "controls_are_independent_of_FGCQR_holdout_outcomes": True,
            "FGCQR_trajectory_evaluated": False,
            "collapse_evolution_performed": False,
            "trapped_interval_measured": False,
            "regulator_activation_measured": False,
            "defocusing_margin_measured": False,
            "retained_EFT_validity_inferred": False,
        },
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    loaded = load_config(config_path)
    paths = loaded["paths"]
    scope = _scope_bindings(loaded)
    quantitative = _quantitative_certificate(loaded)
    predecessor_names = tuple(
        name for name in EXPECTED_PATHS if name not in {"protocol_config", "run1_config"}
    )
    return _serial(
        {
            "schema_version": 1,
            "artifact_id": ARTIFACT_ID,
            "project_version": PROJECT_VERSION,
            "classification": FUTURE_EVIDENCE_CLASSIFICATIONS["numerical_validation"],
            "generated_by": "scripts/reproduce_fgc_num1_val1.py",
            "derivation_document": _rel(OWNER_DOCUMENT),
            "derivation_document_sha256": _sha(OWNER_DOCUMENT),
            "source_config_sha256": {_rel(config_path): _sha(config_path)},
            "predecessor_sha256": {
                _rel(paths[name]): _sha(paths[name]) for name in predecessor_names
            },
            "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
            "scope_bindings": scope,
            "certificate_contract": {
                "artifact_specific_payload_validated": True,
                "canonical_reproduction_passed": True,
                "scope_bindings_verified": True,
                "outcome_data_not_used_to_select_certificate_contract": True,
            },
            "artifact_payload": {
                "run_envelope_semantic_sha256": scope["declared_run_envelope_sha256"],
                "Minkowski_preservation_passed": True,
                "regular_center_manufactured_solution_passed": True,
                "GR0_scalar_control_passed": True,
                "exact_vs_floating_backend_passed": True,
                "constraint_convergence_passed": True,
                "boundary_isolation_passed": True,
                "all_typed_stops_injected": True,
                "checkpoint_restart_equivalence_passed": True,
                "direct_and_Raychaudhuri_routes_agree": True,
                "primary_and_comparator_methods_passed": True,
                "quantitative_evidence": quantitative,
            },
            "gate_status": {REQUIRED_GATE: True},
            "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
        }
    )


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    value = _load_json(path, "NUM1 result")
    if path.read_text(encoding="utf-8") != _canonical(value):
        raise ValueError("NUM1 result must use canonical sorted JSON")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    arguments.output.write_text(_canonical(record(arguments.config)), encoding="utf-8")


if __name__ == "__main__":
    main()
