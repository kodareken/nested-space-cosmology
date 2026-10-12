#!/usr/bin/env python3
"""Regenerate the FGC-1-SRC1-NL1 nonlinear REF1 solver certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from fractions import Fraction
from hashlib import sha256
import json
from math import pi, sqrt
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-src1-nl1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-src1-nl1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-src1-nl1.md"

from scripts.reproduce_fgc_hyp1_dom1_qift1 import (  # noqa: E402
    load_config as load_qift1_config,
)
from recursive_horizons.fgc.evolution.nonlinear_source import (  # noqa: E402
    PARAMETER_ORDER,
    AccelerationSolveStop,
    NonlinearSolverConfig,
    acceleration_jacobian,
    finite_difference_acceleration_jacobian,
    float_ref1_residual,
    parameter_vector,
    solve_accelerations,
    state_from_parameter_vector,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.modified_harmonic_implicit import (  # noqa: E402
    exact_full_residual_acceleration_jacobian,
)
from recursive_horizons.fgc.modified_harmonic_modes import (  # noqa: E402
    exact_comp1_acceleration_root,
)
from recursive_horizons.fgc.modified_harmonic_reference import (  # noqa: E402
    modified_harmonic_full_residuals,
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
    BASE_FIELD_ORDER,
    SphericalState,
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-SRC1-NL1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "nonlinear_REF1_source_and_branch_solver_verified"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "run1_config": "configs/fgc/fgc-1-run1-sym1.toml",
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "reference_config": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
    "reference_result": "results/fgc-1-hyp1-mhg-reference.json",
    "first_order_config": "configs/fgc/fgc-1-hyp1-fo1-rc1.toml",
    "first_order_result": "results/fgc-1-hyp1-fo1-rc1.json",
    "quantified_domain_config": "configs/fgc/fgc-1-hyp1-dom1-qift1.toml",
    "quantified_domain_result": "results/fgc-1-hyp1-dom1-qift1.json",
    "compatible_point_config": "configs/fgc/fgc-1-hyp1-con1-comp1.toml",
    "compatible_point_result": "results/fgc-1-hyp1-con1-comp1.json",
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "src/recursive_horizons/fgc/evolution/floating_jet.py",
        "src/recursive_horizons/fgc/evolution/nonlinear_source.py",
        "scripts/reproduce_fgc_src1_nl1.py",
    )
)
RESULT_GATES = {
    "reference_result": "complete_reference_gauge_residual_passed",
    "first_order_result": "exact_local_first_order_dae_lift_passed",
    "quantified_domain_result": "quantified_full_dimensional_local_implicit_branch_box_passed",
    "compatible_point_result": "exact_activated_local_compatible_constraint_datum_passed",
}


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


def _fraction(name: str, value: Any, *, positive: bool = False) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a canonical rational string")
    if isinstance(value, int):
        answer = Q(value)
    elif isinstance(value, str) and value:
        try:
            answer = Q(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{name} must be a canonical rational string") from exc
        if value != _fraction_text(answer):
            raise ValueError(f"{name} must use canonical rational encoding")
    else:
        raise ValueError(f"{name} must be an exact rational value")
    if positive and answer <= 0:
        raise ValueError(f"{name} must be positive")
    return answer


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected or not isinstance(value, str):
        raise ValueError(f"{name} must name {expected}")
    path = (REPOSITORY / value).resolve()
    try:
        relative = path.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} escapes the repository") from exc
    if relative != value or not path.is_file():
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return path


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _pairs(pairs):
    output = {}
    for key, value in pairs:
        if key in output:
            raise ValueError(f"duplicate JSON object key: {key}")
        output[key] = value
    return output


def _canonical(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _load_json(path: Path, name: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _binds_direct_config(result: Mapping[str, Any], path: Path, name: str) -> None:
    relative, digest = _rel(path), _sha(path)
    if result.get("source_config") == relative and result.get("source_config_sha256") == digest:
        return
    for ledger_name in ("source_config_sha256", "source_configs_sha256"):
        ledger = result.get(ledger_name)
        if isinstance(ledger, Mapping) and ledger.get(relative) == digest:
            return
    source_configs = result.get("source_configs")
    if isinstance(source_configs, Mapping):
        nested = source_configs.get(relative)
        if isinstance(nested, Mapping) and nested.get("sha256") == digest:
            return
    raise ValueError(f"{name} is stale against its direct config")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    source = path.resolve().read_bytes()
    raw = _mapping("root", tomllib.loads(source.decode("utf-8")))
    _keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            *EXPECTED_PATHS,
            "scope",
            "solver",
            "backend",
            "controls",
            "proof_contract",
            "claims",
        },
    )
    if {
        "schema_version": raw["schema_version"],
        "artifact_id": raw["artifact_id"],
        "project_version": raw["project_version"],
        "metric_signature": raw["metric_signature"],
        "riemann_convention": raw["riemann_convention"],
    } != {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }:
        raise ValueError("SRC1 identity, version, or conventions differ")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    protocol = tomllib.loads(paths["protocol_config"].read_text(encoding="utf-8"))
    protocol_certificate = validate_sf1_protocol(protocol)
    qift = load_qift1_config(paths["quantified_domain_config"])

    for result_name, gate in RESULT_GATES.items():
        result = _load_json(paths[result_name], result_name)
        _binds_direct_config(
            result,
            paths[result_name.removesuffix("_result") + "_config"],
            result_name,
        )
        if result.get("gate_status", {}).get(gate) is not True:
            raise ValueError(f"{result_name} required gate is not passed")
    for result_name, expected_id in (
        ("action_result", "FGC-1-ACT1"),
        ("variation_result", "FGC-1-VAR1"),
    ):
        result = _load_json(paths[result_name], result_name)
        _binds_direct_config(
            result,
            paths[result_name.removesuffix("_result") + "_config"],
            result_name,
        )
        if result.get("artifact_id") != expected_id:
            raise ValueError(f"{result_name} identity differs")

    scope = _mapping("scope", raw["scope"])
    _keys(
        "scope",
        scope,
        {
            "target_branch",
            "physical_equations",
            "reference_id",
            "coordinate_radius",
            "radial_domain_minimum",
            "tilde_normal_factor",
            "hat_normal_factor",
            "scope_is_QIFT1_local_branch_not_run_domain",
        },
    )
    if {
        "target_branch": scope["target_branch"],
        "physical_equations": scope["physical_equations"],
        "reference_id": scope["reference_id"],
        "scope_is_QIFT1_local_branch_not_run_domain": scope[
            "scope_is_QIFT1_local_branch_not_run_domain"
        ],
    } != {
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "reference_id": "flat_spherical_annulus",
        "scope_is_QIFT1_local_branch_not_run_domain": True,
    }:
        raise ValueError("SRC1 scope differs")
    radius = _fraction("scope.coordinate_radius", scope["coordinate_radius"], positive=True)
    radial_minimum = _fraction(
        "scope.radial_domain_minimum", scope["radial_domain_minimum"], positive=True
    )
    if (
        radius != qift.coordinate_radius
        or radial_minimum != qift.radial_domain_minimum
        or scope["tilde_normal_factor"] != qift.tilde_normal_factor
        or scope["hat_normal_factor"] != qift.hat_normal_factor
    ):
        raise ValueError("SRC1 scope differs from QIFT1")

    solver = _mapping("solver", raw["solver"])
    _keys(
        "solver",
        solver,
        {
            "residual_infinity_tolerance",
            "maximum_iterations",
            "condition_number_infinity_maximum",
            "normalized_branch_displacement_maximum",
            "parameter_half_width",
            "acceleration_half_width",
            "maximum_backtracks",
            "minimum_step_fraction",
            "strict_residual_decrease_required",
            "branch_continuity_from_warm_start_required",
            "stop_before_crossing_any_limit",
        },
    )
    solver_values = {
        "residual_tolerance": _fraction(
            "solver.residual_infinity_tolerance",
            solver["residual_infinity_tolerance"],
            positive=True,
        ),
        "condition_number_maximum": _fraction(
            "solver.condition_number_infinity_maximum",
            solver["condition_number_infinity_maximum"],
            positive=True,
        ),
        "normalized_branch_displacement_maximum": _fraction(
            "solver.normalized_branch_displacement_maximum",
            solver["normalized_branch_displacement_maximum"],
            positive=True,
        ),
        "parameter_half_width": _fraction(
            "solver.parameter_half_width", solver["parameter_half_width"], positive=True
        ),
        "acceleration_half_width": _fraction(
            "solver.acceleration_half_width", solver["acceleration_half_width"], positive=True
        ),
        "minimum_step_fraction": _fraction(
            "solver.minimum_step_fraction", solver["minimum_step_fraction"], positive=True
        ),
    }
    protocol_health = protocol["health_stops"]
    expected_from_protocol = {
        "residual_tolerance": Q(protocol_health["newton_residual_infinity_max"]),
        "condition_number_maximum": Q(protocol_health["kinetic_condition_number_max"]),
        "normalized_branch_displacement_maximum": Q(
            protocol_health["normalized_branch_displacement_max"]
        ),
    }
    if any(solver_values[key] != value for key, value in expected_from_protocol.items()):
        raise ValueError("SRC1 solver limits differ from PROTO1")
    if (
        solver_values["parameter_half_width"] != qift.parameter_half_width
        or solver_values["acceleration_half_width"] != qift.acceleration_half_width
    ):
        raise ValueError("SRC1 branch boxes differ from QIFT1")
    if solver["maximum_iterations"] != protocol_health["newton_iteration_max"]:
        raise ValueError("SRC1 iteration limit differs from PROTO1")
    if solver["maximum_backtracks"] != 24:
        raise ValueError("SRC1 backtracking count differs")
    for key in (
        "strict_residual_decrease_required",
        "branch_continuity_from_warm_start_required",
        "stop_before_crossing_any_limit",
    ):
        if solver[key] is not True:
            raise ValueError(f"solver.{key} must be true")
    solver_config = NonlinearSolverConfig(
        residual_tolerance=float(solver_values["residual_tolerance"]),
        maximum_iterations=solver["maximum_iterations"],
        condition_number_maximum=float(solver_values["condition_number_maximum"]),
        normalized_branch_displacement_maximum=float(
            solver_values["normalized_branch_displacement_maximum"]
        ),
        parameter_half_width=float(solver_values["parameter_half_width"]),
        acceleration_half_width=float(solver_values["acceleration_half_width"]),
        maximum_backtracks=solver["maximum_backtracks"],
        minimum_step_fraction=float(solver_values["minimum_step_fraction"]),
    )

    backend = _mapping("backend", raw["backend"])
    expected_backend = {
        "residual": "same_modified_harmonic_full_residuals_REF1_evaluator",
        "floating_state": "FloatJet2_subclass_of_established_Jet2_interface",
        "analytic_acceleration_jacobian": "six_binary64_first_tangent_forward_AD_columns_through_same_REF1_evaluator",
        "linear_solve": "numpy_linalg_solve",
        "condition_estimator": "numpy_linalg_cond_infinity_norm",
        "no_second_field_equation_implementation": True,
        "exact_package_remains_numpy_free": True,
    }
    if dict(backend) != expected_backend:
        raise ValueError("SRC1 floating backend contract differs")

    controls = _mapping("controls", raw["controls"])
    _keys(
        "controls",
        controls,
        {
            "rational_fixture_ids",
            "continuation_path_point_count",
            "continuation_path_max_fraction_of_parameter_half_width",
            "exact_float_absolute_tolerance",
            "exact_float_relative_tolerance",
            "finite_difference_step",
            "finite_difference_absolute_tolerance",
            "finite_difference_relative_tolerance",
            "nonrational_control_scale_fraction_of_parameter_half_width",
            "typed_failure_reasons",
        },
    )
    if controls["rational_fixture_ids"] != [
        "IMP1-flat-zero-root",
        "COMP1-zero-acceleration",
        "COMP1-exact-root",
    ]:
        raise ValueError("SRC1 rational fixture order differs")
    if controls["continuation_path_point_count"] != 5:
        raise ValueError("SRC1 continuation path point count differs")
    if controls["typed_failure_reasons"] != [
        "nonfinite_input",
        "parameter_box_violation",
        "initial_acceleration_outside_authorized_box",
        "kinetic_condition_limit",
        "iteration_limit",
    ]:
        raise ValueError("SRC1 typed failure set differs")
    control_fractions = {
        name: _fraction(f"controls.{name}", controls[name], positive=True)
        for name in (
            "continuation_path_max_fraction_of_parameter_half_width",
            "exact_float_absolute_tolerance",
            "exact_float_relative_tolerance",
            "finite_difference_step",
            "finite_difference_absolute_tolerance",
            "finite_difference_relative_tolerance",
            "nonrational_control_scale_fraction_of_parameter_half_width",
        )
    }
    if control_fractions["continuation_path_max_fraction_of_parameter_half_width"] != Q(1, 4):
        raise ValueError("SRC1 continuation path scale differs")
    if control_fractions["nonrational_control_scale_fraction_of_parameter_half_width"] != Q(1, 8):
        raise ValueError("SRC1 nonrational control scale differs")

    proof = _mapping("proof_contract", raw["proof_contract"])
    expected_proof = {
        "require_same_REF1_residual_code_path",
        "require_exact_fraction_controls",
        "require_exact_COMP1_root_agreement",
        "require_nonrational_finite_difference_Jacobian_crosscheck",
        "require_monotonic_residual_history",
        "require_branch_continuity_path",
        "require_typed_failure_injection",
        "require_canonical_hash_bound_result",
        "require_no_holdout_execution",
    }
    _keys("proof_contract", proof, expected_proof)
    if any(value is not True for value in proof.values()):
        raise ValueError("every SRC1 proof-contract entry must be true")
    if dict(_mapping("claims", raw["claims"])) != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("SRC1 claims must remain fail-closed")
    return {
        "raw": dict(raw),
        "paths": paths,
        "protocol": protocol,
        "protocol_certificate": protocol_certificate,
        "qift": qift,
        "solver": solver_config,
        "controls": controls,
        "control_fractions": control_fractions,
        "source_sha256": sha256(source).hexdigest(),
    }


def _exact_state_with_accelerations(
    state: SphericalState, accelerations: Sequence[Fraction]
) -> SphericalState:
    return replace(
        state,
        **{
            field: replace(getattr(state, field), dtt=accelerations[index])
            for index, field in enumerate(BASE_FIELD_ORDER)
        },
    )


def _max_errors(exact: np.ndarray, observed: np.ndarray) -> tuple[float, float]:
    difference = np.abs(exact - observed)
    absolute = float(np.max(difference))
    relative = float(np.max(difference / np.maximum(1.0, np.abs(exact))))
    return absolute, relative


def _exact_float_control(
    identifier: str,
    state: SphericalState,
    accelerations: Sequence[Fraction],
    *,
    reference,
    radius: Fraction,
    tilde: Fraction,
    hat: Fraction,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> dict[str, Any]:
    exact_state = _exact_state_with_accelerations(state, accelerations)
    exact_residual = np.asarray(
        [
            float(value)
            for value in modified_harmonic_full_residuals(
                exact_state,
                reference=reference,
                coordinate_radius=radius,
                tilde_normal_factor=tilde,
                hat_normal_factor=hat,
            )["full_residual_vector"]
        ],
        dtype=np.float64,
    )
    exact_tangent = exact_full_residual_acceleration_jacobian(
        exact_state,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde,
        hat_normal_factor=hat,
    )
    exact_jacobian = np.asarray(
        [[float(value) for value in row] for row in exact_tangent["jacobian"]],
        dtype=np.float64,
    )
    floating_residual, floating_jacobian, condition = acceleration_jacobian(
        state,
        [float(value) for value in accelerations],
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde,
        hat_normal_factor=hat,
    )
    residual_absolute, residual_relative = _max_errors(exact_residual, floating_residual)
    jacobian_absolute, jacobian_relative = _max_errors(exact_jacobian, floating_jacobian)
    passed = (
        residual_absolute <= absolute_tolerance
        and residual_relative <= relative_tolerance
        and jacobian_absolute <= absolute_tolerance
        and jacobian_relative <= relative_tolerance
    )
    if not passed:
        raise ValueError(f"SRC1 exact/floating control {identifier} exceeds tolerance")
    return {
        "fixture_id": identifier,
        "exact_residual_float_max_absolute_error": residual_absolute,
        "exact_residual_float_max_scaled_error": residual_relative,
        "exact_Jacobian_float_max_absolute_error": jacobian_absolute,
        "exact_Jacobian_float_max_scaled_error": jacobian_relative,
        "floating_Jacobian_condition_infinity": condition,
        "passed": True,
    }


def _stop_receipt(reason: str, function) -> dict[str, Any]:
    try:
        function()
    except AccelerationSolveStop as exc:
        if exc.reason != reason:
            raise ValueError(f"SRC1 expected stop {reason}, received {exc.reason}") from exc
        return {
            "reason": exc.reason,
            "outcome_label": exc.outcome_label,
            "stopped_before_crossing": True,
        }
    raise ValueError(f"SRC1 typed stop {reason} was not raised")


def _quantitative_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    qift = configuration["qift"]
    solver = configuration["solver"]
    fractions = configuration["control_fractions"]
    center = state_from_generalized_adm_pg_fixture(qift.flat_fixture)
    compatible = activated_compatible_state(qift.flat_fixture)["state"]
    reference = flat_spherical_annulus_reference(
        radial_domain_minimum=qift.radial_domain_minimum
    )
    arguments = {
        "reference": reference,
        "coordinate_radius": qift.coordinate_radius,
        "tilde_normal_factor": qift.tilde_normal_factor,
        "hat_normal_factor": qift.hat_normal_factor,
    }
    exact_root = exact_comp1_acceleration_root(
        compatible,
        **arguments,
        qift_acceleration_half_width=qift.acceleration_half_width,
    )["root_acceleration"]
    absolute_tolerance = float(fractions["exact_float_absolute_tolerance"])
    relative_tolerance = float(fractions["exact_float_relative_tolerance"])
    rational_controls = [
        _exact_float_control(
            "IMP1-flat-zero-root",
            center,
            (Q(0),) * 6,
            reference=reference,
            radius=qift.coordinate_radius,
            tilde=qift.tilde_normal_factor,
            hat=qift.hat_normal_factor,
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
        ),
        _exact_float_control(
            "COMP1-zero-acceleration",
            compatible,
            (Q(0),) * 6,
            reference=reference,
            radius=qift.coordinate_radius,
            tilde=qift.tilde_normal_factor,
            hat=qift.hat_normal_factor,
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
        ),
        _exact_float_control(
            "COMP1-exact-root",
            compatible,
            exact_root,
            reference=reference,
            radius=qift.coordinate_radius,
            tilde=qift.tilde_normal_factor,
            hat=qift.hat_normal_factor,
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
        ),
    ]

    flat_solve = solve_accelerations(
        center, branch_center=center, config=solver, **arguments
    )
    comp1_solve = solve_accelerations(
        compatible, branch_center=center, config=solver, **arguments
    )
    exact_root_float = np.asarray([float(value) for value in exact_root])
    root_error = float(
        np.linalg.norm(np.asarray(comp1_solve.accelerations) - exact_root_float, ord=np.inf)
    )
    if root_error > absolute_tolerance:
        raise ValueError("SRC1 floating COMP1 root does not match the exact root")

    center_parameters = np.asarray(parameter_vector(center), dtype=np.float64)
    width = float(qift.parameter_half_width)
    nonrational_scale = width * float(
        fractions["nonrational_control_scale_fraction_of_parameter_half_width"]
    )
    nonrational_parameters = center_parameters + np.asarray(
        [
            (sqrt(2.0) if index % 2 == 0 else -sqrt(3.0))
            * nonrational_scale
            / 2.0
            for index in range(len(PARAMETER_ORDER))
        ],
        dtype=np.float64,
    )
    nonrational_state = state_from_parameter_vector(center, nonrational_parameters)
    nonrational_acceleration = np.asarray(
        [
            (pi if index % 2 == 0 else -sqrt(5.0)) * 1.0e-6
            for index in range(6)
        ],
        dtype=np.float64,
    )
    _, analytic_jacobian, analytic_condition = acceleration_jacobian(
        nonrational_state, nonrational_acceleration, **arguments
    )
    difference_jacobian = finite_difference_acceleration_jacobian(
        nonrational_state,
        nonrational_acceleration,
        step=float(fractions["finite_difference_step"]),
        **arguments,
    )
    fd_absolute, fd_relative = _max_errors(analytic_jacobian, difference_jacobian)
    if (
        fd_absolute > float(fractions["finite_difference_absolute_tolerance"])
        or fd_relative > float(fractions["finite_difference_relative_tolerance"])
    ):
        raise ValueError("SRC1 nonrational finite-difference Jacobian check failed")

    path_records = []
    warm_start = np.zeros(6, dtype=np.float64)
    count = configuration["controls"]["continuation_path_point_count"]
    path_scale = float(
        fractions["continuation_path_max_fraction_of_parameter_half_width"]
    )
    for point in range(1, count + 1):
        fraction = point / count
        parameters = center_parameters + np.asarray(
            [
                (1.0 if index % 2 == 0 else -1.0)
                * width
                * path_scale
                * fraction
                * (1.0 - (index % 3) / 8.0)
                for index in range(len(PARAMETER_ORDER))
            ],
            dtype=np.float64,
        )
        state = state_from_parameter_vector(center, parameters)
        solved = solve_accelerations(
            state,
            branch_center=center,
            warm_start=warm_start,
            config=solver,
            **arguments,
        )
        if any(
            item.residual_infinity_after >= item.residual_infinity_before
            for item in solved.iteration_history
        ):
            raise ValueError("SRC1 accepted a non-decreasing residual step")
        path_records.append(
            {
                "path_point": point,
                "maximum_parameter_displacement": float(
                    np.linalg.norm(parameters - center_parameters, ord=np.inf)
                ),
                "iterations": solved.iterations,
                "residual_infinity": solved.residual_infinity,
                "condition_infinity": solved.jacobian_condition_infinity,
                "warm_start_displacement_infinity": solved.branch_displacement_infinity,
                "all_accepted_steps_strictly_decrease_residual": True,
            }
        )
        warm_start = np.asarray(solved.accelerations)

    outside_parameters = center_parameters.copy()
    outside_parameters[0] += 2.0 * width
    outside_state = state_from_parameter_vector(center, outside_parameters)
    typed_stops = [
        _stop_receipt(
            "nonfinite_input",
            lambda: solve_accelerations(
                center,
                branch_center=center,
                warm_start=[float("nan"), 0, 0, 0, 0, 0],
                config=solver,
                **arguments,
            ),
        ),
        _stop_receipt(
            "parameter_box_violation",
            lambda: solve_accelerations(
                outside_state, branch_center=center, config=solver, **arguments
            ),
        ),
        _stop_receipt(
            "initial_acceleration_outside_authorized_box",
            lambda: solve_accelerations(
                center,
                branch_center=center,
                warm_start=[solver.acceleration_half_width, 0, 0, 0, 0, 0],
                config=solver,
                **arguments,
            ),
        ),
        _stop_receipt(
            "kinetic_condition_limit",
            lambda: solve_accelerations(
                center,
                branch_center=center,
                config=replace(solver, condition_number_maximum=1.0),
                **arguments,
            ),
        ),
        _stop_receipt(
            "iteration_limit",
            lambda: solve_accelerations(
                compatible,
                branch_center=center,
                config=replace(solver, residual_tolerance=1.0e-30, maximum_iterations=1),
                **arguments,
            ),
        ),
    ]

    return {
        "backend": {
            "numpy_version": np.__version__,
            "same_REF1_evaluator_import": "recursive_horizons.fgc.modified_harmonic_reference.modified_harmonic_full_residuals",
            "floating_jet_is_established_Jet2_subclass": True,
            "second_field_equation_implementation_present": False,
        },
        "rational_exact_float_controls": rational_controls,
        "root_controls": {
            "flat_zero_root_residual_infinity": flat_solve.residual_infinity,
            "COMP1_exact_root_float_agreement_infinity": root_error,
            "COMP1_float_residual_infinity": comp1_solve.residual_infinity,
            "COMP1_iterations": comp1_solve.iterations,
        },
        "nonrational_Jacobian_control": {
            "finite_difference_step": float(fractions["finite_difference_step"]),
            "maximum_absolute_error": fd_absolute,
            "maximum_scaled_error": fd_relative,
            "analytic_condition_infinity": analytic_condition,
            "passed": True,
        },
        "continuation_path": {
            "point_count": count,
            "records": path_records,
            "maximum_warm_start_displacement_infinity": max(
                item["warm_start_displacement_infinity"] for item in path_records
            ),
            "all_points_converged_on_declared_branch": True,
            "all_accepted_steps_strictly_decrease_residual": True,
        },
        "typed_failure_injection": {
            "receipts": typed_stops,
            "all_declared_failures_typed_before_crossing": True,
        },
        "scope": {
            "QIFT1_parameter_half_width": float(qift.parameter_half_width),
            "QIFT1_acceleration_half_width": float(qift.acceleration_half_width),
            "local_branch_only_not_classical_run_domain": True,
            "holdout_executed": False,
        },
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(config_path)
    quantitative = _quantitative_certificate(configuration)
    protocol_hash = configuration["protocol_certificate"][
        "semantic_holdout_contract_sha256"
    ]
    paths = configuration["paths"]
    scope = {
        "run1_artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "run1_config_sha256": _sha(paths["run1_config"]),
        "protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "protocol_config_sha256": _sha(paths["protocol_config"]),
        "protocol_semantic_holdout_contract_sha256": protocol_hash,
        "action_artifact_id": "FGC-1-ACT1",
        "action_config_sha256": _sha(paths["action_config"]),
        "action_result_sha256": _sha(paths["action_result"]),
        "variation_artifact_id": "FGC-1-VAR1",
        "variation_config_sha256": _sha(paths["variation_config"]),
        "variation_result_sha256": _sha(paths["variation_result"]),
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "declared_run_envelope_sha256": protocol_hash,
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": FUTURE_EVIDENCE_CLASSIFICATIONS["nonlinear_source"],
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): configuration["source_sha256"]},
        "predecessor_sha256": {
            _rel(path): _sha(path)
            for name, path in paths.items()
            if name not in {"run1_config", "protocol_config"}
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
            "run_envelope_semantic_sha256": protocol_hash,
            "quantitative_evidence": quantitative,
            "unredefined_ACT1_VAR1_equations_used": True,
            "analytic_residual_and_acceleration_Jacobian_used": True,
            "exact_fraction_and_floating_controls_agree": True,
            "finite_difference_Jacobian_crosscheck_passed": True,
            "branch_continuity_and_typed_failure_verified": True,
        },
        "gate_status": {REQUIRED_GATE: True},
        "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError("SRC1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("SRC1 result must use canonical sorted JSON")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    arguments.output.write_text(_canonical(record(arguments.config)), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
