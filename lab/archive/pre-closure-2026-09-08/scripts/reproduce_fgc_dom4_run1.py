#!/usr/bin/env python3
"""Regenerate the FGC-1-DOM4-RUN1 classical run-domain certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-dom4-run1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-dom4-run1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-dom4-run1.md"

from recursive_horizons.fgc.evolution.initial_state_bridge import (  # noqa: E402
    solve_initial_second_jet,
)
from recursive_horizons.fgc.evolution.nonlinear_source import (  # noqa: E402
    NonlinearSolverConfig,
)
from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    InitialDataParameters,
    solve_initial_data,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_COMMON_NONCLAIMS,
    FUTURE_EVIDENCE_CLASSIFICATIONS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-DOM4-RUN1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "nonzero_classical_spherical_run_envelope_passed"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "run1_config": "configs/fgc/fgc-1-run1-sym1.toml",
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "radial_domain_config": "configs/fgc/fgc-1-hyp1-dom3-uhyp1.toml",
    "radial_domain_result": "results/fgc-1-hyp1-dom3-uhyp1.json",
    "nonlinear_source_config": "configs/fgc/fgc-1-src1-nl1.toml",
    "nonlinear_source_result": "results/fgc-1-src1-nl1.json",
    "physical_constraints_config": "configs/fgc/fgc-1-con4-phy1.toml",
    "physical_constraints_result": "results/fgc-1-con4-phy1.json",
    "regular_center_config": "configs/fgc/fgc-1-ctr1-reg1.toml",
    "regular_center_result": "results/fgc-1-ctr1-reg1.json",
    "initial_data_config": "configs/fgc/fgc-1-id1-fam1.toml",
    "initial_data_result": "results/fgc-1-id1-fam1.json",
}
RESULT_REQUIREMENTS = {
    "action_result": (
        "FGC-1-ACT1",
        "declared_action_and_scalar_algebra_certificate_passed",
    ),
    "variation_result": (
        "FGC-1-VAR1",
        "var1_covariant_metric_equation_derived_and_cross_checked",
    ),
    "radial_domain_result": (
        "FGC-1-HYP1-DOM3-UHYP1",
        "compact_nonflat_radial_strong_hyperbolicity_on_implicit_branch_graph_passed",
    ),
    "nonlinear_source_result": (
        "FGC-1-SRC1-NL1",
        "nonlinear_REF1_source_and_branch_solver_verified",
    ),
    "physical_constraints_result": (
        "FGC-1-CON4-PHY1",
        "physical_gauge_reduction_constraint_system_closed",
    ),
    "regular_center_result": (
        "FGC-1-CTR1-REG1",
        "regular_center_formulation_verified",
    ),
    "initial_data_result": (
        "FGC-1-ID1-FAM1",
        "nonzero_width_finite_mass_constraint_compatible_family_constructed",
    ),
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/initial_data_family.py",
        "src/recursive_horizons/fgc/evolution/floating_jet.py",
        "src/recursive_horizons/fgc/evolution/nonlinear_source.py",
        "src/recursive_horizons/fgc/evolution/initial_state_bridge.py",
        "scripts/reproduce_fgc_dom4_run1.py",
    )
)


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


def _fraction_list(
    name: str,
    value: Any,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> tuple[Fraction, ...]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{name} must be a nonempty list")
    answer = tuple(
        _fraction(
            f"{name}[{index}]",
            item,
            positive=positive,
            nonnegative=nonnegative,
        )
        for index, item in enumerate(value)
    )
    if len(set(answer)) != len(answer):
        raise ValueError(f"{name} must contain unique values")
    return answer


def _path(name: str, value: Any, expected: str) -> Path:
    if not isinstance(value, str) or value != expected:
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
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
        answer[key] = value
    return answer


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serial(item) for item in value]
    if isinstance(value, np.generic):
        return _serial(value.item())
    if isinstance(value, float) and not isfinite(value):
        raise ValueError("nonfinite value cannot enter canonical evidence")
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"unsupported evidence value {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(_serial(value), indent=2, sort_keys=True, allow_nan=False) + "\n"


def _load_json(path: Path, name: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _binds_direct_config(result: Mapping[str, Any], config: Path, name: str) -> None:
    relative = _rel(config)
    digest = _sha(config)
    if (
        result.get("source_config") == relative
        and result.get("source_config_sha256") == digest
    ):
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


def _all_true_table(name: str, value: Any, expected: set[str]) -> dict[str, bool]:
    table = dict(_mapping(name, value))
    _keys(name, table, expected)
    if any(item is not True for item in table.values()):
        raise ValueError(f"every {name} entry must be true")
    return table


def _semantic_definition(
    raw: Mapping[str, Any], protocol_certificate: Mapping[str, Any]
) -> dict[str, Any]:
    """Return the domain definition independently hashed inside DOM4/HYP2."""

    return {
        "schema_version": 1,
        "protocol_semantic_holdout_contract_sha256": protocol_certificate[
            "semantic_holdout_contract_sha256"
        ],
        "scope": raw["scope"],
        "parameter_envelope": raw["parameter_envelope"],
        "spacetime_envelope": raw["spacetime_envelope"],
        "formulation": raw["formulation"],
        "source_branch": raw["source_branch"],
    }


def domain_definition_sha256(
    raw: Mapping[str, Any], protocol_certificate: Mapping[str, Any]
) -> str:
    return sha256(
        json.dumps(
            _semantic_definition(raw, protocol_certificate),
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


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
            "parameter_envelope",
            "spacetime_envelope",
            "formulation",
            "source_branch",
            "witnesses",
            "proof_contract",
            "claims",
        },
    )
    identity = {
        key: raw[key]
        for key in (
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
        )
    }
    if identity != {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }:
        raise ValueError("DOM4 identity, version, or conventions differ")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    protocol = tomllib.loads(paths["protocol_config"].read_text(encoding="utf-8"))
    protocol_certificate = validate_sf1_protocol(protocol)
    if protocol_certificate["artifact_id"] != SF1_PROTOCOL_ARTIFACT_ID:
        raise ValueError("DOM4 requires active PROTO3")
    for result_name, (artifact_id, gate) in RESULT_REQUIREMENTS.items():
        result = _load_json(paths[result_name], result_name)
        _binds_direct_config(
            result,
            paths[result_name.removesuffix("_result") + "_config"],
            result_name,
        )
        if result.get("artifact_id") != artifact_id:
            raise ValueError(f"{result_name} identity differs")
        if result.get("gate_status", {}).get(gate) is not True:
            raise ValueError(f"{result_name} required gate is not passed")

    scope = dict(_mapping("scope", raw["scope"]))
    if scope != {
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "formulation": "metric_defined_reference_modified_harmonic",
        "state": "u_p_q_with_six_REF1_accelerations",
        "initial_slice": "unit_lapse_zero_shift_maximal_polar_areal",
        "run_envelope_is_classical_not_retained_EFT": True,
        "trajectory_containment_is_deferred_to_HLT1": True,
        "holdout_execution_performed": False,
    }:
        raise ValueError("DOM4 scope differs")

    parameters = dict(_mapping("parameter_envelope", raw["parameter_envelope"]))
    _keys(
        "parameter_envelope",
        parameters,
        {
            "chi_amplitude_minimum",
            "chi_amplitude_maximum",
            "chi_amplitude_candidates",
            "chi_half_width_minimum",
            "chi_half_width_central",
            "chi_half_width_maximum",
            "phi_seed_base_amplitude",
            "phi_seed_factor_minimum",
            "phi_seed_factor_central",
            "phi_seed_factor_maximum",
            "dimensionless_coordinate_volume",
            "physical_phi_amplitude_coordinate_volume",
            "initial_compactness_minimum",
            "initial_compactness_maximum",
            "parameter_container_does_not_waive_per_case_static_admission",
            "resolved_selected_amplitude_is_deferred_to_PRO3_HLD1",
        },
    )
    candidates = _fraction_list(
        "parameter_envelope.chi_amplitude_candidates",
        parameters["chi_amplitude_candidates"],
        positive=True,
    )
    numeric_parameters = {
        key: _fraction(f"parameter_envelope.{key}", parameters[key], positive=True)
        for key in (
            "chi_amplitude_minimum",
            "chi_amplitude_maximum",
            "chi_half_width_minimum",
            "chi_half_width_central",
            "chi_half_width_maximum",
            "phi_seed_base_amplitude",
            "phi_seed_factor_minimum",
            "phi_seed_factor_central",
            "phi_seed_factor_maximum",
            "dimensionless_coordinate_volume",
            "physical_phi_amplitude_coordinate_volume",
            "initial_compactness_minimum",
            "initial_compactness_maximum",
        )
    }
    expected_candidates = tuple(Q(value) for value in protocol["initial_data"]["chi"]["amplitude_candidates"])
    if candidates != expected_candidates or candidates != (
        Q(2), Q(5, 2), Q(3), Q(7, 2), Q(4), Q(9, 2), Q(5)
    ):
        raise ValueError("DOM4 amplitude container differs from PROTO3")
    expected_numeric = {
        "chi_amplitude_minimum": Q(2),
        "chi_amplitude_maximum": Q(5),
        "chi_half_width_minimum": Q(7, 4),
        "chi_half_width_central": Q(2),
        "chi_half_width_maximum": Q(9, 4),
        "phi_seed_base_amplitude": Q(1, 131072),
        "phi_seed_factor_minimum": Q(1, 2),
        "phi_seed_factor_central": Q(1),
        "phi_seed_factor_maximum": Q(2),
        "dimensionless_coordinate_volume": Q(9, 4),
        "physical_phi_amplitude_coordinate_volume": Q(9, 524288),
        "initial_compactness_minimum": Q(1, 10),
        "initial_compactness_maximum": Q(3, 4),
    }
    if numeric_parameters != expected_numeric:
        raise ValueError("DOM4 parameter envelope differs")
    for key in (
        "parameter_container_does_not_waive_per_case_static_admission",
        "resolved_selected_amplitude_is_deferred_to_PRO3_HLD1",
    ):
        if parameters[key] is not True:
            raise ValueError(f"parameter_envelope.{key} must be true")

    spacetime = dict(_mapping("spacetime_envelope", raw["spacetime_envelope"]))
    _keys(
        "spacetime_envelope",
        spacetime,
        {
            "radial_domain_minimum",
            "radial_domain_maximum",
            "measurement_radius_maximum",
            "final_time",
            "minimum_boundary_causal_buffer",
            "length_unit_L0",
        },
    )
    spacetime_values = {
        key: _fraction(
            f"spacetime_envelope.{key}",
            spacetime[key],
            positive=key != "radial_domain_minimum",
            nonnegative=key == "radial_domain_minimum",
        )
        for key in spacetime
    }
    if spacetime_values != {
        "radial_domain_minimum": Q(0),
        "radial_domain_maximum": Q(128),
        "measurement_radius_maximum": Q(24),
        "final_time": Q(32),
        "minimum_boundary_causal_buffer": Q(16),
        "length_unit_L0": Q(4),
    }:
        raise ValueError("DOM4 spacetime envelope differs from PROTO3")

    formulation = dict(_mapping("formulation", raw["formulation"]))
    if formulation != {
        "tilde_normal_factor": 4,
        "hat_normal_factor": 9,
        "normal_factor_order": "1<tilde<hat",
        "independent_chi_sector": "canonical_scalar_direct_sum",
        "reference": "flat_spherical_reference_with_exact_regular_center_descendant",
    }:
        raise ValueError("DOM4 formulation differs")

    source_branch = dict(_mapping("source_branch", raw["source_branch"]))
    _keys(
        "source_branch",
        source_branch,
        {
            "residual_infinity_maximum",
            "maximum_iterations",
            "kinetic_condition_infinity_maximum",
            "branch_displacement_infinity_maximum",
            "acceleration_infinity_maximum",
            "parameter_half_width_for_same_state_solve",
            "maximum_backtracks",
            "minimum_step_fraction",
            "strict_residual_decrease_required",
            "stop_before_crossing_any_limit",
        },
    )
    source_values = {
        key: _fraction(f"source_branch.{key}", source_branch[key], positive=True)
        for key in (
            "residual_infinity_maximum",
            "kinetic_condition_infinity_maximum",
            "branch_displacement_infinity_maximum",
            "acceleration_infinity_maximum",
            "parameter_half_width_for_same_state_solve",
            "minimum_step_fraction",
        )
    }
    health = protocol["health_stops"]
    if source_values != {
        "residual_infinity_maximum": Q(health["newton_residual_infinity_max"]),
        "kinetic_condition_infinity_maximum": Q(health["kinetic_condition_number_max"]),
        "branch_displacement_infinity_maximum": Q(health["normalized_branch_displacement_max"]),
        "acceleration_infinity_maximum": Q(4),
        "parameter_half_width_for_same_state_solve": Q(1),
        "minimum_step_fraction": Q(1, 16777216),
    }:
        raise ValueError("DOM4 source/branch limits differ")
    if (
        source_branch["maximum_iterations"] != health["newton_iteration_max"]
        or source_branch["maximum_backtracks"] != 24
        or source_branch["strict_residual_decrease_required"] is not True
        or source_branch["stop_before_crossing_any_limit"] is not True
    ):
        raise ValueError("DOM4 source/branch discrete contract differs")
    solver = NonlinearSolverConfig(
        residual_tolerance=float(source_values["residual_infinity_maximum"]),
        maximum_iterations=source_branch["maximum_iterations"],
        condition_number_maximum=float(source_values["kinetic_condition_infinity_maximum"]),
        normalized_branch_displacement_maximum=float(source_values["branch_displacement_infinity_maximum"]),
        parameter_half_width=float(source_values["parameter_half_width_for_same_state_solve"]),
        acceleration_half_width=float(source_values["acceleration_infinity_maximum"]),
        maximum_backtracks=source_branch["maximum_backtracks"],
        minimum_step_fraction=float(source_values["minimum_step_fraction"]),
    )

    witnesses = dict(_mapping("witnesses", raw["witnesses"]))
    _keys(
        "witnesses",
        witnesses,
        {
            "constraint_step_count",
            "branch_seed_factors",
            "branch_witness_chi_amplitude",
            "branch_witness_chi_half_width",
            "stress_chi_amplitude",
            "stress_chi_half_width",
            "stress_phi_seed_factor",
            "support_points",
            "strict_witness_plus_continuity_establishes_only_local_open_state_neighborhoods",
        },
    )
    witness_values = {
        key: _fraction(f"witnesses.{key}", witnesses[key], positive=True)
        for key in (
            "branch_witness_chi_amplitude",
            "branch_witness_chi_half_width",
            "stress_chi_amplitude",
            "stress_chi_half_width",
            "stress_phi_seed_factor",
        )
    }
    seed_factors = _fraction_list(
        "witnesses.branch_seed_factors",
        witnesses["branch_seed_factors"],
        nonnegative=True,
    )
    if (
        witnesses["constraint_step_count"] != 512
        or seed_factors != (Q(0), Q(1), Q(2))
        or witness_values
        != {
            "branch_witness_chi_amplitude": Q(2),
            "branch_witness_chi_half_width": Q(2),
            "stress_chi_amplitude": Q(5),
            "stress_chi_half_width": Q(9, 4),
            "stress_phi_seed_factor": Q(2),
        }
        or witnesses["support_points"] != ["profile_center", "peak_compactness"]
        or witnesses[
            "strict_witness_plus_continuity_establishes_only_local_open_state_neighborhoods"
        ]
        is not True
    ):
        raise ValueError("DOM4 witness contract differs")

    _all_true_table(
        "proof_contract",
        raw["proof_contract"],
        {
            "all_protocol_parameter_coordinates_must_lie_inside_container",
            "container_membership_must_not_substitute_for_static_case_admission",
            "strict_nonsingular_REF1_witnesses_required",
            "seed_continuation_must_preserve_one_acceleration_branch",
            "stress_initial_slice_must_be_regular_finite_mass_and_untrapped",
            "radial_DOM3_is_consumed_but_not_promoted_to_all_covectors",
            "HYP2_must_independently_close_all_covectors",
            "HLT1_must_monitor_trajectory_containment_at_every_accepted_stage",
            "canonical_hash_bound_result_required",
            "no_evolution_or_holdout_outcome_access",
        },
    )
    if dict(_mapping("claims", raw["claims"])) != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("DOM4 claims must remain fail-closed")
    return {
        "raw": dict(raw),
        "paths": paths,
        "protocol": protocol,
        "protocol_certificate": protocol_certificate,
        "parameters": numeric_parameters,
        "candidates": candidates,
        "spacetime": spacetime_values,
        "source": source_values,
        "solver": solver,
        "seed_factors": seed_factors,
        "witness_values": witness_values,
        "constraint_step_count": witnesses["constraint_step_count"],
        "domain_definition_sha256": domain_definition_sha256(raw, protocol_certificate),
        "source_sha256": sha256(source).hexdigest(),
    }


def _point_index(solution, target: float) -> int:
    return min(
        range(2, len(solution.points) - 2),
        key=lambda index: abs(solution.points[index].radius - target),
    )


def _solve_summary(jet) -> dict[str, Any]:
    result = jet.acceleration_solve
    return {
        "coordinate_radius": jet.coordinate_radius,
        "accelerations": result.accelerations,
        "residual_infinity": result.residual_infinity,
        "iterations": result.iterations,
        "kinetic_condition_infinity": result.jacobian_condition_infinity,
        "branch_displacement_infinity": result.branch_displacement_infinity,
        "acceleration_infinity": max(abs(value) for value in result.accelerations),
        "monotonic_residual_history": all(
            item.residual_infinity_after < item.residual_infinity_before
            for item in result.iteration_history
        ),
        "branch_continuity_preserved": result.branch_continuity_preserved,
        "spatial_second_derivative_estimators": {
            "lambda_rr": jet.spatial_diagnostics.lambda_rr_estimator,
            "k_rr": jet.spatial_diagnostics.k_rr_estimator,
        },
    }


def _parameters(
    *, amplitude: Fraction, width: Fraction, seed_factor: Fraction
) -> InitialDataParameters:
    return InitialDataParameters(
        chi_amplitude=float(amplitude),
        chi_half_width=float(width),
        phi_amplitude=float(Q(1, 131072) * seed_factor),
    )


def _quantitative_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    solver = configuration["solver"]
    steps = configuration["constraint_step_count"]
    witness = configuration["witness_values"]
    seed_records = []
    warm_start: Sequence[float] | None = None
    for factor in configuration["seed_factors"]:
        parameters = _parameters(
            amplitude=witness["branch_witness_chi_amplitude"],
            width=witness["branch_witness_chi_half_width"],
            seed_factor=factor,
        )
        solution = solve_initial_data(parameters, step_count=steps, method="RK4")
        index = _point_index(solution, parameters.center)
        jet = solve_initial_second_jet(
            solution,
            index,
            warm_start=warm_start,
            solver_config=solver,
        )
        warm_start = jet.acceleration_solve.accelerations
        seed_records.append(
            {
                "phi_seed_factor": factor,
                "initial_peak_compactness": solution.peak_compactness,
                "initial_regular_center": solution.regular_center,
                "initial_finite_mass": solution.finite_mass,
                "initial_untrapped": solution.no_initial_trapped_sphere,
                **_solve_summary(jet),
            }
        )

    stress_parameters = _parameters(
        amplitude=witness["stress_chi_amplitude"],
        width=witness["stress_chi_half_width"],
        seed_factor=witness["stress_phi_seed_factor"],
    )
    stress = solve_initial_data(stress_parameters, step_count=steps, method="RK4")
    stress_records = []
    for identifier, target in (
        ("profile_center", stress_parameters.center),
        ("peak_compactness", stress.peak_radius),
    ):
        jet = solve_initial_second_jet(
            stress,
            _point_index(stress, target),
            solver_config=solver,
        )
        stress_records.append({"support_point": identifier, **_solve_summary(jet)})

    source_records = seed_records + stress_records
    residual_limit = float(configuration["source"]["residual_infinity_maximum"])
    condition_limit = float(
        configuration["source"]["kinetic_condition_infinity_maximum"]
    )
    acceleration_limit = float(
        configuration["source"]["acceleration_infinity_maximum"]
    )
    displacement_limit = float(
        configuration["source"]["branch_displacement_infinity_maximum"]
    )
    strict_source_margins = all(
        record["residual_infinity"] < residual_limit
        and record["kinetic_condition_infinity"] < condition_limit
        and record["acceleration_infinity"] < acceleration_limit
        and record["branch_displacement_infinity"] <= displacement_limit
        and record["monotonic_residual_history"]
        and record["branch_continuity_preserved"]
        for record in source_records
    )
    if not strict_source_margins:
        raise ValueError("DOM4 witness missed a strict REF1 source/branch margin")
    if not (
        stress.regular_center
        and stress.finite_mass
        and stress.no_initial_trapped_sphere
        and stress.compactness_inside_protocol_window
    ):
        raise ValueError("DOM4 stress initial slice is not protocol-admissible")

    protocol = configuration["protocol"]
    held_width_factors = tuple(
        Q(value) for value in protocol["partition"]["held_out_chi_width_factors"]
    )
    held_seed_factors = tuple(
        Q(value) for value in protocol["partition"]["held_out_phi_seed_factors"]
    )
    p = configuration["parameters"]
    all_coordinates_contained = (
        all(p["chi_amplitude_minimum"] <= value <= p["chi_amplitude_maximum"] for value in configuration["candidates"])
        and all(
            p["chi_half_width_minimum"]
            <= p["chi_half_width_central"] * value
            <= p["chi_half_width_maximum"]
            for value in held_width_factors
        )
        and all(
            p["phi_seed_factor_minimum"] <= value <= p["phi_seed_factor_maximum"]
            for value in held_seed_factors
        )
    )
    if not all_coordinates_contained:
        raise ValueError("DOM4 does not contain every PROTO3 parameter coordinate")

    return {
        "domain_definition_sha256": configuration["domain_definition_sha256"],
        "parameter_container": {
            "chi_amplitude": [p["chi_amplitude_minimum"], p["chi_amplitude_maximum"]],
            "chi_half_width": [p["chi_half_width_minimum"], p["chi_half_width_maximum"]],
            "phi_seed_factor": [p["phi_seed_factor_minimum"], p["phi_seed_factor_maximum"]],
            "dimensionless_coordinate_volume": p["dimensionless_coordinate_volume"],
            "physical_phi_amplitude_coordinate_volume": p[
                "physical_phi_amplitude_coordinate_volume"
            ],
            "all_PROTO3_candidate_and_holdout_coordinates_contained": True,
            "membership_is_not_static_case_admission": True,
        },
        "spacetime_container": configuration["spacetime"],
        "seed_branch_continuation": {
            "records": seed_records,
            "strict_nonsingular_branch_margins": True,
            "smooth_REF1_coefficients_plus_invertible_acceleration_Jacobian_give_local_open_neighborhoods": True,
        },
        "stress_initial_slice": {
            "chi_amplitude": witness["stress_chi_amplitude"],
            "chi_half_width": witness["stress_chi_half_width"],
            "phi_seed_factor": witness["stress_phi_seed_factor"],
            "peak_compactness": stress.peak_compactness,
            "peak_radius": stress.peak_radius,
            "outer_mass": stress.outer_mass,
            "regular_center": stress.regular_center,
            "finite_mass": stress.finite_mass,
            "no_initial_trapped_sphere": stress.no_initial_trapped_sphere,
            "inside_protocol_compactness_window": stress.compactness_inside_protocol_window,
            "support_records": stress_records,
        },
        "strict_source_and_branch_margins_passed": True,
        "epistemic_boundary": {
            "parameter_container_is_not_an_all_case_initial_data_certificate": True,
            "local_open_state_neighborhoods_are_not_a_global_run_tube": True,
            "multidirectional_health_is_deferred_to_HYP2": True,
            "trajectory_containment_is_deferred_to_HLT1": True,
            "time_evolution_performed": False,
            "FGCQR_holdout_outcome_inspected": False,
        },
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(config_path)
    quantitative = _quantitative_certificate(configuration)
    paths = configuration["paths"]
    protocol_hash = configuration["protocol_certificate"][
        "semantic_holdout_contract_sha256"
    ]
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
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": FUTURE_EVIDENCE_CLASSIFICATIONS["run_domain"],
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
            "nonzero_parameter_volume": True,
            "FGCQR_target_branch_only": True,
            "all_protocol_cases_contained": True,
            "branch_continuity_quantified": True,
            "classical_scope_only": True,
        },
        "gate_status": {REQUIRED_GATE: True},
        "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
    }
    return _serial(result)


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError("DOM4 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("DOM4 result must use canonical sorted JSON")
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
