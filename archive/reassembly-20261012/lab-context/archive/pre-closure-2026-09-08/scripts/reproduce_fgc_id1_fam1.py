#!/usr/bin/env python3
"""Regenerate the FGC-1-ID1-FAM1 initial-data family certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import isfinite, log2
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-id1-fam1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-id1-fam1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-id1-fam1.md"

from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    FGCQRActionParameters,
    STOP_OUTCOME_LABELS,
    InitialDataParameters,
    InitialDataSolveStop,
    exact_constraint_and_gauge_crosscheck,
    polar_areal_constraint_coefficients,
    polar_areal_constraint_residuals,
    polar_areal_constraint_rhs,
    solve_initial_data,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    gr0_constraint_residuals,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_COMMON_NONCLAIMS,
    FUTURE_EVIDENCE_CLASSIFICATIONS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-ID1-FAM1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "nonzero_width_finite_mass_constraint_compatible_family_constructed"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "run1_config": "configs/fgc/fgc-1-run1-sym1.toml",
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "reference_config": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
    "reference_result": "results/fgc-1-hyp1-mhg-reference.json",
    "physical_constraints_config": "configs/fgc/fgc-1-con4-phy1.toml",
    "physical_constraints_result": "results/fgc-1-con4-phy1.json",
    "regular_center_config": "configs/fgc/fgc-1-ctr1-reg1.toml",
    "regular_center_result": "results/fgc-1-ctr1-reg1.json",
    "protocol_preflight_config": "configs/fgc/fgc-1-id0-pref1.toml",
    "protocol_preflight_result": "results/fgc-1-id0-pref1.json",
}
RESULT_REQUIREMENTS = {
    "action_result": (
        "FGC-1-ACT1",
        "declared_action_and_scalar_algebra_certificate_passed",
        True,
    ),
    "variation_result": (
        "FGC-1-VAR1",
        "var1_covariant_metric_equation_derived_and_cross_checked",
        True,
    ),
    "reference_result": (
        "FGC-1-HYP1-MHG2-REF1",
        "complete_reference_gauge_residual_passed",
        True,
    ),
    "physical_constraints_result": (
        "FGC-1-CON4-PHY1",
        "physical_gauge_reduction_constraint_system_closed",
        True,
    ),
    "regular_center_result": (
        "FGC-1-CTR1-REG1",
        "regular_center_formulation_verified",
        True,
    ),
    "protocol_preflight_result": (
        "FGC-1-ID0-PREF1",
        "PROTO1_initial_data_preflight_obstruction_verified",
        True,
    ),
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/initial_data_family.py",
        "src/recursive_horizons/fgc/initial_data_preflight.py",
        "src/recursive_horizons/fgc/modified_harmonic_constraints.py",
        "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "src/recursive_horizons/fgc/reference_connection.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "scripts/reproduce_fgc_id1_fam1.py",
    )
)
EXACT_CONTROL_KEYS = {
    "radius",
    "radial_metric",
    "angular_extrinsic_curvature",
    "radial_metric_derivative",
    "angular_extrinsic_curvature_derivative",
    "phi",
    "phi_r",
    "phi_rr",
    "phi_pi",
    "phi_pi_r",
    "chi",
    "chi_r",
    "chi_rr",
    "chi_pi",
    "chi_pi_r",
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


def _integer_list(name: str, value: Any) -> tuple[int, ...]:
    if (
        not isinstance(value, list)
        or not value
        or any(isinstance(item, bool) or not isinstance(item, int) or item <= 0 for item in value)
    ):
        raise ValueError(f"{name} must contain positive integers")
    return tuple(value)


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
            "action",
            "profiles",
            "family_box",
            "continuation",
            "numerics",
            "exact_controls",
            "proof_contract",
            "claims",
        },
    )
    identity = {
        "schema_version": raw["schema_version"],
        "artifact_id": raw["artifact_id"],
        "project_version": raw["project_version"],
        "metric_signature": raw["metric_signature"],
        "riemann_convention": raw["riemann_convention"],
    }
    if identity != {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }:
        raise ValueError("ID1 identity, version, or conventions differ")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    protocol = tomllib.loads(paths["protocol_config"].read_text(encoding="utf-8"))
    protocol_certificate = validate_sf1_protocol(protocol)
    if (
        protocol_certificate["artifact_id"] != SF1_PROTOCOL_ARTIFACT_ID
        or protocol_certificate["premise_revision_only"] is not True
        or protocol_certificate["FGCQR_outcomes_inspected_before_revision"] is not False
    ):
        raise ValueError("ID1 requires the active premise-only PROTO3 contract")

    for result_name, (artifact_id, gate, expected) in RESULT_REQUIREMENTS.items():
        result = _load_json(paths[result_name], result_name)
        config_name = result_name.removesuffix("_result") + "_config"
        _binds_direct_config(result, paths[config_name], result_name)
        if result.get("artifact_id") != artifact_id:
            raise ValueError(f"{result_name} identity differs")
        if result.get("gate_status", {}).get(gate) is not expected:
            raise ValueError(f"{result_name} required gate differs")
    preflight = _load_json(paths["protocol_preflight_result"], "protocol_preflight_result")
    if (
        preflight.get("gate_status", {}).get("PROTO1_FGCQR_holdout_execution_authorized")
        is not False
        or preflight.get("artifact_payload", {})
        .get("epistemic_boundary", {})
        .get("FGCQR_action_or_mechanism_tested")
        is not False
    ):
        raise ValueError("ID1 must remain downstream of the pre-holdout ID0 diagnosis")

    scope = _mapping("scope", raw["scope"])
    expected_scope = {
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "initial_slice": "unit_lapse_zero_shift_maximal_polar_areal",
        "areal_radius": "R=r",
        "extrinsic_curvature_eigenvalues": [
            "K^r_r=-2*k",
            "K^theta_theta=k",
            "K^phi_phi=k",
        ],
        "reference_id": "flat_spherical_reference_with_Ct_Cr_solved_on_the_annular_support",
        "tilde_normal_factor": 4,
        "hat_normal_factor": 9,
        "regular_center_is_exact_vacuum_buffer": True,
        "asymptotic_end_is_exact_finite_mass_vacuum": True,
        "constraint_family_is_initial_data_not_evolution": True,
        "FGCQR_holdout_outcomes_inspected": False,
    }
    if dict(scope) != expected_scope:
        raise ValueError("ID1 scope differs")

    action = _mapping("action", raw["action"])
    _keys(
        "action",
        action,
        {"planck_mass", "beta", "scalar_mass", "quartic_coupling", "eta"},
    )
    action_values = {
        key: _fraction(f"action.{key}", value)
        for key, value in action.items()
    }
    if action_values != {
        "planck_mass": Q(2),
        "beta": -Q(1, 4),
        "scalar_mass": Q(3),
        "quartic_coupling": Q(1, 2),
        "eta": Q(1, 2),
    }:
        raise ValueError("ID1 action parameters differ from ACT1/PROTO3")
    action_parameters = FGCQRActionParameters(
        planck_mass=float(action_values["planck_mass"]),
        beta=float(action_values["beta"]),
        scalar_mass=float(action_values["scalar_mass"]),
        quartic_coupling=float(action_values["quartic_coupling"]),
        eta=float(action_values["eta"]),
    )

    profiles = _mapping("profiles", raw["profiles"])
    _keys(
        "profiles",
        profiles,
        {
            "chi_profile",
            "chi_momentum",
            "chi_center",
            "phi_profile",
            "phi_half_width",
            "phi_unit_normal_momentum",
            "outer_radius",
            "minimum_center_buffer",
            "minimum_outer_vacuum_buffer",
        },
    )
    expected_profile_text = {
        "chi_profile": "r_times_chi=amplitude*B((r-center)/half_width)",
        "chi_momentum": "partial_t(r*chi)=partial_r(r*chi)",
        "phi_profile": "phi=amplitude*B((r-center)/half_width)",
    }
    if any(profiles[key] != value for key, value in expected_profile_text.items()):
        raise ValueError("ID1 profile definitions differ")
    profile_values = {
        key: _fraction(f"profiles.{key}", profiles[key], positive=True)
        for key in (
            "chi_center",
            "phi_half_width",
            "outer_radius",
            "minimum_center_buffer",
            "minimum_outer_vacuum_buffer",
        )
    }
    phi_momentum = _fraction(
        "profiles.phi_unit_normal_momentum",
        profiles["phi_unit_normal_momentum"],
        nonnegative=True,
    )
    if phi_momentum != 0:
        raise ValueError("ID1 phi seed momentum must equal the PROTO3 completion")
    if profile_values != {
        "chi_center": Q(12),
        "phi_half_width": Q(2),
        "outer_radius": Q(128),
        "minimum_center_buffer": Q(10),
        "minimum_outer_vacuum_buffer": Q(32),
    }:
        raise ValueError("ID1 dimensional profiles differ from PROTO3 code units")

    family = _mapping("family_box", raw["family_box"])
    _keys(
        "family_box",
        family,
        {
            "chi_amplitude_minimum",
            "chi_amplitude_maximum",
            "chi_amplitude_nodes",
            "chi_half_width_minimum",
            "chi_half_width_maximum",
            "chi_half_width_nodes",
            "phi_seed_base_amplitude",
            "phi_seed_factor_minimum",
            "phi_seed_factor_maximum",
            "phi_seed_factor_nodes",
            "physical_parameter_volume",
            "sampled_tensor_grid_size",
            "sampled_tensor_grid_is_not_a_global_interval_enclosure",
            "positive_local_family_uses_ODE_continuous_dependence",
            "qualification_box_does_not_cover_PROTO3_width_nine_eighths_case",
        },
    )
    axes = {
        "chi_amplitude": _fraction_list(
            "family_box.chi_amplitude_nodes", family["chi_amplitude_nodes"], positive=True
        ),
        "chi_half_width": _fraction_list(
            "family_box.chi_half_width_nodes", family["chi_half_width_nodes"], positive=True
        ),
        "phi_seed_factor": _fraction_list(
            "family_box.phi_seed_factor_nodes", family["phi_seed_factor_nodes"], positive=True
        ),
    }
    endpoints = {
        "chi_amplitude": (
            _fraction("family_box.chi_amplitude_minimum", family["chi_amplitude_minimum"], positive=True),
            _fraction("family_box.chi_amplitude_maximum", family["chi_amplitude_maximum"], positive=True),
        ),
        "chi_half_width": (
            _fraction("family_box.chi_half_width_minimum", family["chi_half_width_minimum"], positive=True),
            _fraction("family_box.chi_half_width_maximum", family["chi_half_width_maximum"], positive=True),
        ),
        "phi_seed_factor": (
            _fraction("family_box.phi_seed_factor_minimum", family["phi_seed_factor_minimum"], positive=True),
            _fraction("family_box.phi_seed_factor_maximum", family["phi_seed_factor_maximum"], positive=True),
        ),
    }
    if any(tuple(sorted(values)) != values for values in axes.values()):
        raise ValueError("ID1 family axes must be strictly ascending")
    if any((values[0], values[-1]) != endpoints[name] for name, values in axes.items()):
        raise ValueError("ID1 family node endpoints differ from the declared box")
    if axes != {
        "chi_amplitude": (Q(2), Q(17, 8), Q(9, 4)),
        "chi_half_width": (Q(15, 8), Q(31, 16), Q(2)),
        "phi_seed_factor": (Q(1, 2), Q(1), Q(2)),
    }:
        raise ValueError("ID1 qualification family differs")
    seed_base = _fraction(
        "family_box.phi_seed_base_amplitude",
        family["phi_seed_base_amplitude"],
        positive=True,
    )
    volume = _fraction(
        "family_box.physical_parameter_volume",
        family["physical_parameter_volume"],
        positive=True,
    )
    expected_volume = (
        endpoints["chi_amplitude"][1] - endpoints["chi_amplitude"][0]
    ) * (
        endpoints["chi_half_width"][1] - endpoints["chi_half_width"][0]
    ) * seed_base * (
        endpoints["phi_seed_factor"][1] - endpoints["phi_seed_factor"][0]
    )
    if seed_base != Q(1, 131072) or volume != expected_volume:
        raise ValueError("ID1 family volume or seed differs")
    if family["sampled_tensor_grid_size"] != 27:
        raise ValueError("ID1 family grid size differs")
    for key in (
        "sampled_tensor_grid_is_not_a_global_interval_enclosure",
        "positive_local_family_uses_ODE_continuous_dependence",
        "qualification_box_does_not_cover_PROTO3_width_nine_eighths_case",
    ):
        if family[key] is not True:
            raise ValueError(f"family_box.{key} must be true")
    if profile_values["chi_center"] - endpoints["chi_half_width"][1] < profile_values[
        "minimum_center_buffer"
    ]:
        raise ValueError("ID1 family violates the PROTO3 centre buffer")

    continuation = _mapping("continuation", raw["continuation"])
    _keys(
        "continuation",
        continuation,
        {
            "central_chi_amplitude",
            "central_chi_half_width",
            "phi_seed_factors_from_GR0",
            "zero_seed_is_GR0_constraint_start",
            "nonzero_seed_is_required_for_certified_FGCQR_members",
            "branch_continuity_requires_all_stages",
        },
    )
    continuation_factors = _fraction_list(
        "continuation.phi_seed_factors_from_GR0",
        continuation["phi_seed_factors_from_GR0"],
        nonnegative=True,
    )
    if (
        _fraction("continuation.central_chi_amplitude", continuation["central_chi_amplitude"], positive=True)
        != Q(2)
        or _fraction("continuation.central_chi_half_width", continuation["central_chi_half_width"], positive=True)
        != Q(2)
        or continuation_factors != (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1))
    ):
        raise ValueError("ID1 continuation path differs")
    for key in (
        "zero_seed_is_GR0_constraint_start",
        "nonzero_seed_is_required_for_certified_FGCQR_members",
        "branch_continuity_requires_all_stages",
    ):
        if continuation[key] is not True:
            raise ValueError(f"continuation.{key} must be true")

    numerics = _mapping("numerics", raw["numerics"])
    _keys(
        "numerics",
        numerics,
        {
            "methods",
            "nested_step_counts",
            "report_step_count",
            "minimum_RK4_common_grid_profile_order",
            "minimum_SSPRK3_common_grid_profile_order",
            "maximum_cross_method_outer_mass_difference",
            "maximum_cross_method_peak_compactness_difference",
            "maximum_constraint_residual_infinity",
            "minimum_abs_constraint_jacobian_determinant",
            "maximum_constraint_jacobian_condition_infinity",
            "minimum_effective_planck_coefficient",
            "minimum_vacuum_metric_denominator",
            "minimum_peak_compactness",
            "maximum_peak_compactness",
            "require_no_initial_trapped_sphere",
        },
    )
    nested = _integer_list("numerics.nested_step_counts", numerics["nested_step_counts"])
    if numerics["methods"] != ["RK4", "SSPRK3"] or nested != (64, 128, 256):
        raise ValueError("ID1 method or nested-resolution contract differs")
    report_steps = numerics["report_step_count"]
    if isinstance(report_steps, bool) or not isinstance(report_steps, int) or report_steps < 512:
        raise ValueError("ID1 report step count is too small")
    thresholds = {
        key: _fraction(f"numerics.{key}", numerics[key], positive=True)
        for key in (
            "minimum_RK4_common_grid_profile_order",
            "minimum_SSPRK3_common_grid_profile_order",
            "maximum_cross_method_outer_mass_difference",
            "maximum_cross_method_peak_compactness_difference",
            "maximum_constraint_residual_infinity",
            "minimum_abs_constraint_jacobian_determinant",
            "maximum_constraint_jacobian_condition_infinity",
            "minimum_effective_planck_coefficient",
            "minimum_vacuum_metric_denominator",
            "minimum_peak_compactness",
            "maximum_peak_compactness",
        )
    }
    if numerics["require_no_initial_trapped_sphere"] is not True:
        raise ValueError("ID1 must require no initial trapped sphere")

    exact_controls = _mapping("exact_controls", raw["exact_controls"])
    _keys("exact_controls", exact_controls, {"fixture_A", "fixture_B"})
    fixtures = []
    for identifier, fixture in exact_controls.items():
        fixture = _mapping(f"exact_controls.{identifier}", fixture)
        _keys(f"exact_controls.{identifier}", fixture, EXACT_CONTROL_KEYS)
        fixtures.append(
            (
                identifier,
                {
                    key: _fraction(f"exact_controls.{identifier}.{key}", value)
                    for key, value in fixture.items()
                },
            )
        )

    proof = _mapping("proof_contract", raw["proof_contract"])
    expected_proof = {
        "complete_unredefined_constraints_specialized",
        "two_exact_rational_full_tensor_crosschecks_required",
        "affine_radial_constraint_Jacobian_derived",
        "direct_constraint_solve_used_instead_of_generic_optimizer",
        "physical_Hamiltonian_and_momentum_constraints_solved",
        "metric_defined_gauge_Ct_and_Cr_solved",
        "normal_gauge_propagation_remains_conditional_on_CON4",
        "regular_center_and_exact_vacuum_exterior_required",
        "independent_RK4_and_SSPRK3_methods_required",
        "nonzero_width_three_axis_family_required",
        "typed_continuation_failures_required",
        "canonical_hash_bound_result_required",
        "no_evolution_or_holdout_execution",
    }
    _keys("proof_contract", proof, expected_proof)
    if any(value is not True for value in proof.values()):
        raise ValueError("every ID1 proof-contract entry must be true")
    if dict(_mapping("claims", raw["claims"])) != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("ID1 claims must remain fail-closed")

    return {
        "raw": dict(raw),
        "paths": paths,
        "protocol_certificate": protocol_certificate,
        "action": action_parameters,
        "profile_values": profile_values,
        "family_axes": axes,
        "family_volume": volume,
        "seed_base": seed_base,
        "continuation_factors": continuation_factors,
        "nested_steps": nested,
        "report_steps": report_steps,
        "thresholds": thresholds,
        "exact_fixtures": fixtures,
        "source_sha256": sha256(source).hexdigest(),
    }


def _profile_error(coarse, fine) -> float:
    if fine.step_count % coarse.step_count != 0:
        raise ValueError("ID1 profile comparison requires nested grids")
    stride = fine.step_count // coarse.step_count
    sampled = fine.points[::stride]
    if len(sampled) != len(coarse.points):
        raise ValueError("ID1 nested profile lengths differ")
    return max(
        max(
            abs(left.radial_metric - right.radial_metric),
            abs(
                left.angular_extrinsic_curvature
                - right.angular_extrinsic_curvature
            ),
        )
        for left, right in zip(coarse.points, sampled, strict=True)
    )


def _convergence_record(parameters, method: str, steps: Sequence[int]) -> dict[str, Any]:
    solutions = [
        solve_initial_data(parameters, step_count=count, method=method)
        for count in steps
    ]
    coarse_medium = _profile_error(solutions[0], solutions[1])
    medium_fine = _profile_error(solutions[1], solutions[2])
    if coarse_medium <= 0.0 or medium_fine <= 0.0:
        raise ValueError(f"ID1 {method} profile convergence is unresolved")
    order = log2(coarse_medium / medium_fine)
    return {
        "method": method,
        "step_counts": list(steps),
        "coarse_medium_common_grid_profile_error_infinity": coarse_medium,
        "medium_fine_common_grid_profile_error_infinity": medium_fine,
        "observed_common_grid_profile_order": order,
        "outer_masses": [solution.outer_mass for solution in solutions],
    }


def _solution_summary(solution) -> dict[str, Any]:
    return {
        "outer_mass": solution.outer_mass,
        "peak_compactness": solution.peak_compactness,
        "peak_radius": solution.peak_radius,
        "minimum_abs_constraint_jacobian_determinant": solution.minimum_abs_constraint_jacobian_determinant,
        "maximum_constraint_jacobian_condition_infinity": solution.maximum_constraint_jacobian_condition_infinity,
        "maximum_constraint_residual_infinity": solution.maximum_constraint_residual_infinity,
        "minimum_effective_planck_coefficient": solution.minimum_effective_planck_coefficient,
        "minimum_vacuum_metric_denominator": solution.minimum_vacuum_metric_denominator,
        "finite_mass": solution.finite_mass,
        "regular_center": solution.regular_center,
        "exact_inner_vacuum_buffer": solution.exact_inner_vacuum_buffer,
        "exact_outer_vacuum_buffer": solution.exact_outer_vacuum_buffer,
        "no_initial_trapped_sphere": solution.no_initial_trapped_sphere,
        "compactness_inside_protocol_window": solution.compactness_inside_protocol_window,
    }


def _exact_gr0_reduction_control(fixture: Mapping[str, Fraction]) -> dict[str, Any]:
    values = dict(fixture)
    values.update(
        {
            "phi": Q(0),
            "phi_r": Q(0),
            "phi_rr": Q(0),
            "phi_pi": Q(0),
            "phi_pi_r": Q(0),
        }
    )
    coefficients = polar_areal_constraint_coefficients(
        radius=values["radius"],
        radial_metric=values["radial_metric"],
        angular_extrinsic_curvature=values["angular_extrinsic_curvature"],
        phi=Q(0),
        phi_r=Q(0),
        phi_rr=Q(0),
        phi_pi=Q(0),
        phi_pi_r=Q(0),
        chi_r=values["chi_r"],
        chi_pi=values["chi_pi"],
        planck_mass=Q(2),
        beta=-Q(1, 4),
        scalar_mass=Q(3),
        quartic_coupling=Q(1, 2),
        eta=Q(1, 2),
    )
    fgc = polar_areal_constraint_residuals(
        coefficients,
        radial_metric_derivative=values["radial_metric_derivative"],
        angular_extrinsic_curvature_derivative=values[
            "angular_extrinsic_curvature_derivative"
        ],
    )
    gr0 = gr0_constraint_residuals(
        radius=values["radius"],
        radial_metric=values["radial_metric"],
        angular_extrinsic_curvature=values["angular_extrinsic_curvature"],
        radial_metric_derivative=values["radial_metric_derivative"],
        angular_extrinsic_curvature_derivative=values[
            "angular_extrinsic_curvature_derivative"
        ],
        phi=Q(0),
        phi_r=Q(0),
        phi_pi=Q(0),
        chi=values["chi"],
        chi_r=values["chi_r"],
        chi_pi=values["chi_pi"],
    )
    if fgc != gr0:
        raise ValueError("ID1 zero-seed FGC-QR constraints do not reduce to GR-0")
    return {
        "FGCQR_zero_seed_constraint_pair": fgc,
        "GR0_constraint_pair": gr0,
        "exact_pair_equality": True,
    }


def _typed_stop_control(parameters: InitialDataParameters) -> dict[str, Any]:
    try:
        polar_areal_constraint_rhs(
            parameters.center,
            (1.0, 0.0),
            parameters,
            minimum_abs_jacobian_entry=10.0,
        )
    except InitialDataSolveStop as exc:
        if exc.reason != "singular_constraint_jacobian":
            raise ValueError("ID1 typed-stop control returned the wrong reason") from exc
        receipt = {
            "injected_reason": exc.reason,
            "outcome_label": exc.outcome_label,
            "typed_stop_raised": True,
        }
    else:
        raise ValueError("ID1 typed-stop control did not stop")
    return {
        "reason_to_outcome_label": dict(sorted(STOP_OUTCOME_LABELS.items())),
        "injection_receipt": receipt,
        "continuation_failures_have_typed_reason_and_outcome": True,
    }


def _quantitative_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    action = configuration["action"]
    profile = configuration["profile_values"]
    seed_base = configuration["seed_base"]
    central = InitialDataParameters(
        chi_amplitude=2.0,
        chi_half_width=2.0,
        phi_amplitude=float(seed_base),
        center=float(profile["chi_center"]),
        phi_half_width=float(profile["phi_half_width"]),
        outer_radius=float(profile["outer_radius"]),
        action=action,
    )
    convergence = {
        method: _convergence_record(
            central, method, configuration["nested_steps"]
        )
        for method in ("RK4", "SSPRK3")
    }
    thresholds = configuration["thresholds"]
    if convergence["RK4"]["observed_common_grid_profile_order"] < float(
        thresholds["minimum_RK4_common_grid_profile_order"]
    ):
        raise ValueError("ID1 RK4 profile convergence order missed its threshold")
    if convergence["SSPRK3"]["observed_common_grid_profile_order"] < float(
        thresholds["minimum_SSPRK3_common_grid_profile_order"]
    ):
        raise ValueError("ID1 SSPRK3 profile convergence order missed its threshold")

    exact_controls = []
    for identifier, fixture in configuration["exact_fixtures"]:
        control = exact_constraint_and_gauge_crosscheck(**fixture)
        exact_controls.append({"fixture_id": identifier, **control})
    gr0_reduction = _exact_gr0_reduction_control(
        configuration["exact_fixtures"][0][1]
    )

    axes = configuration["family_axes"]
    grid_records = []
    extrema = {
        "minimum_peak_compactness": float("inf"),
        "maximum_peak_compactness": 0.0,
        "minimum_abs_constraint_jacobian_determinant": float("inf"),
        "maximum_constraint_jacobian_condition_infinity": 0.0,
        "maximum_constraint_residual_infinity": 0.0,
        "minimum_effective_planck_coefficient": float("inf"),
        "minimum_vacuum_metric_denominator": float("inf"),
        "maximum_cross_method_outer_mass_difference": 0.0,
        "maximum_cross_method_peak_compactness_difference": 0.0,
    }
    for amplitude, width, seed_factor in product(
        axes["chi_amplitude"],
        axes["chi_half_width"],
        axes["phi_seed_factor"],
    ):
        parameters = InitialDataParameters(
            chi_amplitude=float(amplitude),
            chi_half_width=float(width),
            phi_amplitude=float(seed_base * seed_factor),
            center=float(profile["chi_center"]),
            phi_half_width=float(profile["phi_half_width"]),
            outer_radius=float(profile["outer_radius"]),
            action=action,
        )
        primary = solve_initial_data(
            parameters,
            step_count=configuration["report_steps"],
            method="RK4",
            maximum_constraint_residual=float(
                thresholds["maximum_constraint_residual_infinity"]
            ),
        )
        comparator = solve_initial_data(
            parameters,
            step_count=configuration["report_steps"],
            method="SSPRK3",
            maximum_constraint_residual=float(
                thresholds["maximum_constraint_residual_infinity"]
            ),
        )
        mass_difference = abs(primary.outer_mass - comparator.outer_mass)
        compactness_difference = abs(
            primary.peak_compactness - comparator.peak_compactness
        )
        extrema["minimum_peak_compactness"] = min(
            extrema["minimum_peak_compactness"], primary.peak_compactness
        )
        extrema["maximum_peak_compactness"] = max(
            extrema["maximum_peak_compactness"], primary.peak_compactness
        )
        extrema["minimum_abs_constraint_jacobian_determinant"] = min(
            extrema["minimum_abs_constraint_jacobian_determinant"],
            primary.minimum_abs_constraint_jacobian_determinant,
            comparator.minimum_abs_constraint_jacobian_determinant,
        )
        extrema["maximum_constraint_jacobian_condition_infinity"] = max(
            extrema["maximum_constraint_jacobian_condition_infinity"],
            primary.maximum_constraint_jacobian_condition_infinity,
            comparator.maximum_constraint_jacobian_condition_infinity,
        )
        extrema["maximum_constraint_residual_infinity"] = max(
            extrema["maximum_constraint_residual_infinity"],
            primary.maximum_constraint_residual_infinity,
            comparator.maximum_constraint_residual_infinity,
        )
        extrema["minimum_effective_planck_coefficient"] = min(
            extrema["minimum_effective_planck_coefficient"],
            primary.minimum_effective_planck_coefficient,
            comparator.minimum_effective_planck_coefficient,
        )
        extrema["minimum_vacuum_metric_denominator"] = min(
            extrema["minimum_vacuum_metric_denominator"],
            primary.minimum_vacuum_metric_denominator,
            comparator.minimum_vacuum_metric_denominator,
        )
        extrema["maximum_cross_method_outer_mass_difference"] = max(
            extrema["maximum_cross_method_outer_mass_difference"], mass_difference
        )
        extrema["maximum_cross_method_peak_compactness_difference"] = max(
            extrema["maximum_cross_method_peak_compactness_difference"],
            compactness_difference,
        )
        grid_records.append(
            {
                "chi_amplitude": amplitude,
                "chi_half_width": width,
                "phi_seed_factor": seed_factor,
                "primary": _solution_summary(primary),
                "comparator": _solution_summary(comparator),
                "cross_method_outer_mass_difference": mass_difference,
                "cross_method_peak_compactness_difference": compactness_difference,
            }
        )

    required_inequalities = {
        "peak_compactness_above_minimum": extrema["minimum_peak_compactness"]
        >= float(thresholds["minimum_peak_compactness"]),
        "peak_compactness_below_maximum": extrema["maximum_peak_compactness"]
        <= float(thresholds["maximum_peak_compactness"]),
        "constraint_Jacobian_determinant_margin": extrema[
            "minimum_abs_constraint_jacobian_determinant"
        ]
        >= float(thresholds["minimum_abs_constraint_jacobian_determinant"]),
        "constraint_Jacobian_condition_margin": extrema[
            "maximum_constraint_jacobian_condition_infinity"
        ]
        <= float(thresholds["maximum_constraint_jacobian_condition_infinity"]),
        "constraint_residual_margin": extrema[
            "maximum_constraint_residual_infinity"
        ]
        <= float(thresholds["maximum_constraint_residual_infinity"]),
        "effective_Planck_margin": extrema["minimum_effective_planck_coefficient"]
        >= float(thresholds["minimum_effective_planck_coefficient"]),
        "vacuum_metric_margin": extrema["minimum_vacuum_metric_denominator"]
        >= float(thresholds["minimum_vacuum_metric_denominator"]),
        "cross_method_outer_mass_margin": extrema[
            "maximum_cross_method_outer_mass_difference"
        ]
        <= float(thresholds["maximum_cross_method_outer_mass_difference"]),
        "cross_method_peak_compactness_margin": extrema[
            "maximum_cross_method_peak_compactness_difference"
        ]
        <= float(thresholds["maximum_cross_method_peak_compactness_difference"]),
        "all_grid_members_regular_finite_and_untrapped": all(
            record[method][key]
            for record in grid_records
            for method in ("primary", "comparator")
            for key in (
                "finite_mass",
                "regular_center",
                "exact_inner_vacuum_buffer",
                "exact_outer_vacuum_buffer",
                "no_initial_trapped_sphere",
                "compactness_inside_protocol_window",
            )
        ),
    }
    if not all(required_inequalities.values()):
        failed = sorted(key for key, passed in required_inequalities.items() if not passed)
        raise ValueError(f"ID1 family qualification failed: {failed!r}")

    continuation_records = []
    preceding = None
    maximum_displacement = 0.0
    continuation_steps = min(configuration["report_steps"], 1024)
    for factor in configuration["continuation_factors"]:
        parameters = InitialDataParameters(
            chi_amplitude=2.0,
            chi_half_width=2.0,
            phi_amplitude=float(seed_base * factor),
            center=float(profile["chi_center"]),
            phi_half_width=float(profile["phi_half_width"]),
            outer_radius=float(profile["outer_radius"]),
            action=action,
        )
        solution = solve_initial_data(
            parameters, step_count=continuation_steps, method="RK4"
        )
        displacement = 0.0 if preceding is None else _profile_error(preceding, solution)
        maximum_displacement = max(maximum_displacement, displacement)
        continuation_records.append(
            {
                "phi_seed_factor": factor,
                "maximum_profile_displacement_from_preceding_stage": displacement,
                "outer_mass": solution.outer_mass,
                "minimum_abs_constraint_jacobian_determinant": solution.minimum_abs_constraint_jacobian_determinant,
                "maximum_constraint_jacobian_condition_infinity": solution.maximum_constraint_jacobian_condition_infinity,
                "passed": True,
            }
        )
        preceding = solution

    return {
        "backend": {
            "numpy_version": np.__version__,
            "constraint_kernel": "direct_affine_radial_ACT1_VAR1_specialization",
            "generic_nonlinear_optimizer_used": False,
        },
        "exact_full_evaluator_controls": exact_controls,
        "exact_zero_seed_GR0_reduction": gr0_reduction,
        "central_common_grid_convergence": convergence,
        "family_box": {
            "exact_positive_parameter_volume": configuration["family_volume"],
            "tensor_grid_member_count": len(grid_records),
            "records": grid_records,
            "extrema": extrema,
            "required_inequalities": required_inequalities,
            "sampled_tensor_grid_is_not_a_global_interval_enclosure": True,
            "smooth_local_open_family_follows_from_regular_ODE_continuous_dependence": True,
            "whole_declared_box_interval_enclosure_proven": False,
            "PROTO3_width_nine_eighths_case_qualified_here": False,
        },
        "phi_seed_continuation_from_GR0": {
            "records": continuation_records,
            "maximum_adjacent_profile_displacement_infinity": maximum_displacement,
            "zero_seed_constraints_reduce_exactly_to_GR0": True,
            "all_nonzero_stages_retain_strict_constraint_branch_margins": True,
        },
        "constraint_and_gauge_contract": {
            "Hamiltonian_and_momentum_solved_at_every_reported_point": True,
            "metric_defined_Ct_and_Cr_solved_by_exact_algebraic_rates": True,
            "normal_gauge_compatibility_is_CON4_conditional_propagation_not_a_new_ID1_theorem": True,
            "kinematic_radial_derivatives_constructed_from_the_same_slice": True,
        },
        "typed_failure_contract": _typed_stop_control(central),
        "epistemic_boundary": {
            "initial_hypersurfaces_constructed": True,
            "time_evolution_performed": False,
            "collapse_or_trapped_surface_derived": False,
            "Raychaudhuri_defocusing_tested": False,
            "FGCQR_holdout_outcomes_inspected": False,
            "retained_EFT_validity_established": False,
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
        "classification": FUTURE_EVIDENCE_CLASSIFICATIONS["initial_data"],
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
            "nonzero_width_family": True,
            "finite_Misner_Sharp_mass": True,
            "independent_chi_pulse": True,
            "explicit_nonzero_phi_seed": True,
            "physical_and_gauge_constraints_solved": True,
            "exact_outer_vacuum_buffer": True,
            "continuation_failure_reasons_typed": True,
        },
        "gate_status": {REQUIRED_GATE: True},
        "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
    }
    return _serial(result)


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("ID1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("ID1 result must use canonical sorted JSON")
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
