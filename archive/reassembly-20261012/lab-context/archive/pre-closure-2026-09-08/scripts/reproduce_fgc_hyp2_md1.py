#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP2-MD1 all-covector health certificate."""

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
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hyp2-md1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hyp2-md1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hyp2-md1.md"

from scripts.reproduce_fgc_dom4_run1 import (  # noqa: E402
    load_config as load_dom4_config,
)
from scripts.reproduce_fgc_hyp1_dom1_qift1 import (  # noqa: E402
    load_config as load_qift1_config,
)
from recursive_horizons.fgc.evolution.covariant_principal_health import (  # noqa: E402
    background_from_spherical_state,
    second_order_coefficient_matrices,
)
from recursive_horizons.fgc.evolution.initial_state_bridge import (  # noqa: E402
    solve_initial_second_jet,
)
from recursive_horizons.fgc.evolution.multidirectional_health import (  # noqa: E402
    CLUSTER_SPECS,
    WeakCouplingThresholds,
    canonical_health_monitor_values,
    esf_reference_cluster_certificate,
    weak_coupling_health_certificate,
)
from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    InitialDataParameters,
    solve_initial_data,
)
from recursive_horizons.fgc.modified_harmonic import (  # noqa: E402
    modified_harmonic_symbol,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_COMMON_NONCLAIMS,
    FUTURE_EVIDENCE_CLASSIFICATIONS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-HYP2-MD1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "quantitative_all_covector_weak_coupling_health_envelope_passed"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "run1_config": "configs/fgc/fgc-1-run1-sym1.toml",
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "modified_harmonic_config": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
    "modified_harmonic_result": "results/fgc-1-hyp1-modified-harmonic.json",
    "radial_domain_config": "configs/fgc/fgc-1-hyp1-dom3-uhyp1.toml",
    "radial_domain_result": "results/fgc-1-hyp1-dom3-uhyp1.json",
    "run_domain_config": "configs/fgc/fgc-1-dom4-run1.toml",
    "run_domain_result": "results/fgc-1-dom4-run1.json",
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
    "modified_harmonic_result": (
        "FGC-1-HYP1-MHG1",
        "exact_pointwise_modified_harmonic_principal_gate_passed",
    ),
    "radial_domain_result": (
        "FGC-1-HYP1-DOM3-UHYP1",
        "compact_nonflat_radial_strong_hyperbolicity_on_implicit_branch_graph_passed",
    ),
    "run_domain_result": (
        "FGC-1-DOM4-RUN1",
        "nonzero_classical_spherical_run_envelope_passed",
    ),
    "initial_data_result": (
        "FGC-1-ID1-FAM1",
        "nonzero_width_finite_mass_constraint_compatible_family_constructed",
    ),
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/modified_harmonic.py",
        "src/recursive_horizons/fgc/evolution/covariant_principal_health.py",
        "src/recursive_horizons/fgc/evolution/initial_state_bridge.py",
        "src/recursive_horizons/fgc/evolution/multidirectional_health.py",
        "scripts/reproduce_fgc_hyp2_md1.py",
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
            "auxiliary_cones",
            "weak_coupling_thresholds",
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
        raise ValueError("HYP2 identity, version, or conventions differ")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    protocol = tomllib.loads(paths["protocol_config"].read_text(encoding="utf-8"))
    protocol_certificate = validate_sf1_protocol(protocol)
    if protocol_certificate["artifact_id"] != SF1_PROTOCOL_ARTIFACT_ID:
        raise ValueError("HYP2 requires active PROTO3")
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
    dom4 = load_dom4_config(paths["run_domain_config"])
    dom4_result = _load_json(paths["run_domain_result"], "run_domain_result")
    dom4_quantitative = dom4_result.get("artifact_payload", {}).get(
        "quantitative_evidence", {}
    )
    if (
        dom4_quantitative.get("domain_definition_sha256")
        != dom4["domain_definition_sha256"]
        or dom4_result.get("scope_bindings", {}).get(
            "protocol_semantic_holdout_contract_sha256"
        )
        != protocol_certificate["semantic_holdout_contract_sha256"]
    ):
        raise ValueError("HYP2 and DOM4 bind different run-domain definitions")

    scope = dict(_mapping("scope", raw["scope"]))
    if scope != {
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "principal_system": "full_four_dimensional_metric_phi_chi_modified_harmonic_symbol",
        "metric_field_components": 10,
        "scalar_field_components": 2,
        "first_order_characteristic_dimension": 24,
        "all_nonzero_spatial_covectors_by_homogeneity": True,
        "classical_principal_health_only": True,
        "trajectory_containment_is_deferred_to_HLT1": True,
        "retained_EFT_remainder_control_supplied": False,
        "holdout_execution_performed": False,
    }:
        raise ValueError("HYP2 scope differs")

    cones = dict(_mapping("auxiliary_cones", raw["auxiliary_cones"]))
    _keys(
        "auxiliary_cones",
        cones,
        {
            "tilde_normal_factor",
            "hat_normal_factor",
            "cluster_centers",
            "cluster_multiplicities",
            "physical_cluster_radius",
            "auxiliary_cluster_radius",
            "physical_resolvent_singular_lower_bound",
            "auxiliary_resolvent_singular_lower_bound",
            "reference_contour_samples",
        },
    )
    centers = tuple(_fraction(f"auxiliary_cones.cluster_centers[{i}]", value) for i, value in enumerate(cones["cluster_centers"]))
    if (
        cones["tilde_normal_factor"] != 4
        or cones["hat_normal_factor"] != 9
        or centers != (-Q(1), -Q(1, 2), -Q(1, 3), Q(1, 3), Q(1, 2), Q(1))
        or cones["cluster_multiplicities"] != [4, 4, 4, 4, 4, 4]
        or _fraction("auxiliary_cones.physical_cluster_radius", cones["physical_cluster_radius"], positive=True) != Q(1, 5)
        or _fraction("auxiliary_cones.auxiliary_cluster_radius", cones["auxiliary_cluster_radius"], positive=True) != Q(1, 20)
        or _fraction("auxiliary_cones.physical_resolvent_singular_lower_bound", cones["physical_resolvent_singular_lower_bound"], positive=True) != Q(1, 32)
        or _fraction("auxiliary_cones.auxiliary_resolvent_singular_lower_bound", cones["auxiliary_resolvent_singular_lower_bound"], positive=True) != Q(1, 1024)
        or cones["reference_contour_samples"] != 512
    ):
        raise ValueError("HYP2 auxiliary-cone contract differs")
    expected_clusters = tuple(
        (float(center), float(radius), multiplicity, float(lower))
        for _name, center, radius, multiplicity, lower in CLUSTER_SPECS
    )
    configured_clusters = tuple(
        (
            float(center),
            float(Q(1, 5) if abs(center) == 1 else Q(1, 20)),
            4,
            float(Q(1, 32) if abs(center) == 1 else Q(1, 1024)),
        )
        for center in centers
    )
    if configured_clusters != expected_clusters:
        raise ValueError("HYP2 config and implementation cluster contracts differ")

    thresholds_raw = dict(
        _mapping("weak_coupling_thresholds", raw["weak_coupling_thresholds"])
    )
    _keys(
        "weak_coupling_thresholds",
        thresholds_raw,
        {
            "companion_deformation_operator_2_maximum",
            "action_energy_deformation_operator_2_maximum",
            "kinetic_deformation_operator_2_maximum",
            "effective_planck_squared_minimum",
            "action_hessian_symmetry_defect_maximum",
            "pure_gauge_identity_coefficient_bound_maximum",
            "physical_Riesz_projector_displacement_maximum",
            "physical_energy_coercivity_minimum",
        },
    )
    threshold_values = {
        key: _fraction(
            f"weak_coupling_thresholds.{key}",
            thresholds_raw[key],
            positive=key != "physical_energy_coercivity_minimum",
            nonnegative=key == "physical_energy_coercivity_minimum",
        )
        for key in thresholds_raw
    }
    expected_thresholds = {
        "companion_deformation_operator_2_maximum": Q(1, 4096),
        "action_energy_deformation_operator_2_maximum": Q(1, 8192),
        "kinetic_deformation_operator_2_maximum": Q(1, 64),
        "effective_planck_squared_minimum": Q(2),
        "action_hessian_symmetry_defect_maximum": Q(1, 10000000000),
        "pure_gauge_identity_coefficient_bound_maximum": Q(1, 10000000000),
        "physical_Riesz_projector_displacement_maximum": Q(1),
        "physical_energy_coercivity_minimum": Q(0),
    }
    if threshold_values != expected_thresholds:
        raise ValueError("HYP2 weak-coupling thresholds differ")
    thresholds = WeakCouplingThresholds(
        companion_deformation_maximum=float(
            threshold_values["companion_deformation_operator_2_maximum"]
        ),
        action_energy_deformation_maximum=float(
            threshold_values["action_energy_deformation_operator_2_maximum"]
        ),
        kinetic_deformation_maximum=float(
            threshold_values["kinetic_deformation_operator_2_maximum"]
        ),
        effective_planck_squared_minimum=float(
            threshold_values["effective_planck_squared_minimum"]
        ),
        action_hessian_symmetry_defect_maximum=float(
            threshold_values["action_hessian_symmetry_defect_maximum"]
        ),
        gauge_identity_defect_maximum=float(
            threshold_values["pure_gauge_identity_coefficient_bound_maximum"]
        ),
    )

    witnesses = dict(_mapping("witnesses", raw["witnesses"]))
    _keys(
        "witnesses",
        witnesses,
        {
            "constraint_step_count",
            "case_ids",
            "central_chi_amplitude",
            "central_chi_half_width",
            "central_phi_seed_factor",
            "double_phi_seed_factor",
            "stress_chi_amplitude",
            "stress_chi_half_width",
            "stress_phi_seed_factor",
            "stress_points",
            "witnesses_demonstrate_nonempty_envelope_not_trajectory_containment",
        },
    )
    witness_values = {
        key: _fraction(f"witnesses.{key}", witnesses[key], positive=True)
        for key in (
            "central_chi_amplitude",
            "central_chi_half_width",
            "central_phi_seed_factor",
            "double_phi_seed_factor",
            "stress_chi_amplitude",
            "stress_chi_half_width",
            "stress_phi_seed_factor",
        )
    }
    if (
        witnesses["constraint_step_count"] != 512
        or witnesses["case_ids"]
        != ["CENTRAL-SEED", "CENTRAL-SEED-DOUBLE", "STRESS-CENTER", "STRESS-PEAK"]
        or witness_values
        != {
            "central_chi_amplitude": Q(2),
            "central_chi_half_width": Q(2),
            "central_phi_seed_factor": Q(1),
            "double_phi_seed_factor": Q(2),
            "stress_chi_amplitude": Q(5),
            "stress_chi_half_width": Q(9, 4),
            "stress_phi_seed_factor": Q(2),
        }
        or witnesses["stress_points"] != ["profile_center", "peak_compactness"]
        or witnesses[
            "witnesses_demonstrate_nonempty_envelope_not_trajectory_containment"
        ]
        is not True
    ):
        raise ValueError("HYP2 witness contract differs")

    _all_true_table(
        "proof_contract",
        raw["proof_contract"],
        {
            "VAR1_action_Hessian_row_normalization_must_be_used",
            "Frobenius_metric_basis_must_make_spatial_rotations_orthogonal",
            "complete_covector_coefficient_tensor_must_be_extracted",
            "directional_scan_cannot_substitute_for_coefficient_norm_bound",
            "six_Riesz_contours_must_partition_all_characteristics",
            "physical_action_energy_form_must_remain_coercive",
            "pure_gauge_identity_must_be_checked_coefficient_by_coefficient",
            "independent_chi_sector_must_be_retained",
            "exact_spherical_radial_symbol_crosscheck_required",
            "every_declared_witness_must_pass_strictly",
            "HLT1_must_stop_before_exit_from_the_health_ball",
            "canonical_hash_bound_result_required",
            "no_evolution_or_holdout_outcome_access",
        },
    )
    if dict(_mapping("claims", raw["claims"])) != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("HYP2 claims must remain fail-closed")
    return {
        "raw": dict(raw),
        "paths": paths,
        "protocol_certificate": protocol_certificate,
        "dom4": dom4,
        "thresholds": thresholds,
        "threshold_values": threshold_values,
        "witness_values": witness_values,
        "case_ids": tuple(witnesses["case_ids"]),
        "constraint_step_count": witnesses["constraint_step_count"],
        "reference_contour_samples": cones["reference_contour_samples"],
        "source_sha256": sha256(source).hexdigest(),
    }


def _point_index(solution, target: float) -> int:
    return min(
        range(2, len(solution.points) - 2),
        key=lambda index: abs(solution.points[index].radius - target),
    )


def _parameters(
    *, amplitude: Fraction, width: Fraction, seed_factor: Fraction
) -> InitialDataParameters:
    return InitialDataParameters(
        chi_amplitude=float(amplitude),
        chi_half_width=float(width),
        phi_amplitude=float(Q(1, 131072) * seed_factor),
    )


def _radial_symbol_crosscheck() -> dict[str, Any]:
    exact_state = activated_compatible_state(load_qift1_config().flat_fixture)[
        "state"
    ]
    background = background_from_spherical_state(exact_state)
    a_matrix, b_matrix, c_matrix = second_order_coefficient_matrices(
        background, (1.0, 0.0, 0.0)
    )
    embedding = np.zeros((12, 6))
    for column, row in enumerate((0, 1, 4)):
        embedding[row, column] = 1.0
    embedding[7, 3] = embedding[9, 3] = 1.0
    embedding[10, 4] = embedding[11, 5] = 1.0
    reduced_a = embedding.T @ a_matrix @ embedding
    reduced_b = embedding.T @ b_matrix @ embedding
    reduced_c = embedding.T @ c_matrix @ embedding
    covariant = np.block(
        [
            [np.zeros((6, 6)), np.eye(6)],
            [
                -np.linalg.solve(reduced_a, reduced_c),
                -np.linalg.solve(reduced_a, reduced_b),
            ],
        ]
    )
    covariant_roots = np.linalg.eigvals(covariant)

    exact_symbol = modified_harmonic_symbol(
        exact_state, tilde_normal_factor=4, hat_normal_factor=9
    )
    exact_a = np.asarray(
        [[float(entry[2] if len(entry) > 2 else 0) for entry in row] for row in exact_symbol]
    )
    exact_b = np.asarray(
        [[float(entry[1] if len(entry) > 1 else 0) for entry in row] for row in exact_symbol]
    )
    exact_c = np.asarray(
        [[float(entry[0] if entry else 0) for entry in row] for row in exact_symbol]
    )
    exact_companion = np.block(
        [
            [np.zeros((6, 6)), np.eye(6)],
            [-np.linalg.solve(exact_a, exact_c), -np.linalg.solve(exact_a, exact_b)],
        ]
    )
    coordinate_speeds = np.linalg.eigvals(exact_companion)
    radial_metric = float(exact_state.h_rr.value) ** 0.5
    expected = -radial_metric * coordinate_speeds
    ordered_covariant = np.sort_complex(covariant_roots)
    ordered_expected = np.sort_complex(expected)
    maximum = float(np.max(np.abs(ordered_covariant - ordered_expected)))
    if maximum > 2.0e-13:
        raise ValueError("HYP2 covariant symbol misses the independent exact radial symbol")
    return {
        "fixture_id": "COMP1-activated-compatible-radial-restriction",
        "root_count": len(covariant_roots),
        "maximum_absolute_root_difference": maximum,
        "tolerance": 2.0e-13,
        "independent_exact_spherical_symbol_agreement": True,
    }


def _health_record(
    identifier: str,
    solution,
    point_index: int,
    thresholds: WeakCouplingThresholds,
) -> dict[str, Any]:
    jet = solve_initial_second_jet(solution, point_index)
    background = background_from_spherical_state(jet.accelerated_state)
    certificate = weak_coupling_health_certificate(
        background, thresholds=thresholds
    )
    if not certificate[
        "quantitative_all_covector_weak_coupling_health_envelope_passed"
    ]:
        raise ValueError(f"HYP2 witness {identifier} misses the health envelope")
    return {
        "case_id": identifier,
        "coordinate_radius": jet.coordinate_radius,
        "initial_peak_compactness": solution.peak_compactness,
        "initial_outer_mass": solution.outer_mass,
        "source_residual_infinity": jet.acceleration_solve.residual_infinity,
        "source_kinetic_condition_infinity": jet.acceleration_solve.jacobian_condition_infinity,
        "background": certificate["background"],
        "deformation_bounds": certificate["deformation_bounds"],
        "cluster_resolvent_singular_margins": certificate[
            "cluster_resolvent_singular_margins"
        ],
        "physical_projector_displacement_upper": certificate[
            "physical_projector_displacement_upper"
        ],
        "transported_reference_energy_coercivity_lower": certificate[
            "transported_reference_energy_coercivity_lower"
        ],
        "physical_energy_transport_error_upper": certificate[
            "physical_energy_transport_error_upper"
        ],
        "physical_energy_coercivity_lower": certificate[
            "physical_energy_coercivity_lower"
        ],
        "monitor_values": canonical_health_monitor_values(certificate),
        "predicates": certificate["predicates"],
        "passed": True,
    }


def _quantitative_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    witness = configuration["witness_values"]
    steps = configuration["constraint_step_count"]
    thresholds = configuration["thresholds"]
    central = solve_initial_data(
        _parameters(
            amplitude=witness["central_chi_amplitude"],
            width=witness["central_chi_half_width"],
            seed_factor=witness["central_phi_seed_factor"],
        ),
        step_count=steps,
    )
    double = solve_initial_data(
        _parameters(
            amplitude=witness["central_chi_amplitude"],
            width=witness["central_chi_half_width"],
            seed_factor=witness["double_phi_seed_factor"],
        ),
        step_count=steps,
    )
    stress = solve_initial_data(
        _parameters(
            amplitude=witness["stress_chi_amplitude"],
            width=witness["stress_chi_half_width"],
            seed_factor=witness["stress_phi_seed_factor"],
        ),
        step_count=steps,
    )
    records = [
        _health_record(
            "CENTRAL-SEED",
            central,
            _point_index(central, central.parameters.center),
            thresholds,
        ),
        _health_record(
            "CENTRAL-SEED-DOUBLE",
            double,
            _point_index(double, double.parameters.center),
            thresholds,
        ),
        _health_record(
            "STRESS-CENTER",
            stress,
            _point_index(stress, stress.parameters.center),
            thresholds,
        ),
        _health_record(
            "STRESS-PEAK",
            stress,
            _point_index(stress, stress.peak_radius),
            thresholds,
        ),
    ]
    if tuple(item["case_id"] for item in records) != configuration["case_ids"]:
        raise ValueError("HYP2 witness order differs")
    values = configuration["threshold_values"]
    strict = all(
        record["deformation_bounds"]["companion_operator_2"]
        < float(values["companion_deformation_operator_2_maximum"])
        and record["deformation_bounds"]["action_energy_operator_2"]
        < float(values["action_energy_deformation_operator_2_maximum"])
        and record["deformation_bounds"]["kinetic_operator_2"]
        < float(values["kinetic_deformation_operator_2_maximum"])
        and record["background"]["effective_planck_squared"]
        > float(values["effective_planck_squared_minimum"])
        and record["physical_projector_displacement_upper"]
        < float(values["physical_Riesz_projector_displacement_maximum"])
        and record["physical_energy_coercivity_lower"]
        > float(values["physical_energy_coercivity_minimum"])
        and all(record["predicates"].values())
        for record in records
    )
    if not strict:
        raise ValueError("HYP2 witness misses a strict all-covector margin")

    reference = esf_reference_cluster_certificate(
        configuration["reference_contour_samples"]
    )
    extrema = {
        "maximum_companion_deformation_operator_2": max(
            item["deformation_bounds"]["companion_operator_2"] for item in records
        ),
        "maximum_action_energy_deformation_operator_2": max(
            item["deformation_bounds"]["action_energy_operator_2"] for item in records
        ),
        "maximum_kinetic_deformation_operator_2": max(
            item["deformation_bounds"]["kinetic_operator_2"] for item in records
        ),
        "maximum_action_hessian_symmetry_defect": max(
            item["deformation_bounds"]["action_hessian_symmetry_defect"] for item in records
        ),
        "maximum_pure_gauge_identity_coefficient_bound": max(
            item["deformation_bounds"]["gauge_identity_defect"] for item in records
        ),
        "maximum_physical_projector_displacement_upper": max(
            item["physical_projector_displacement_upper"] for item in records
        ),
        "minimum_physical_energy_coercivity_lower": min(
            item["physical_energy_coercivity_lower"] for item in records
        ),
        "minimum_cluster_resolvent_singular_margin": min(
            min(item["cluster_resolvent_singular_margins"].values()) for item in records
        ),
    }
    return {
        "domain_definition_sha256": configuration["dom4"][
            "domain_definition_sha256"
        ],
        "theorem_specialization": {
            "source": "Kovacs_and_Reall_modified_harmonic_weak_coupling_continuity_argument",
            "full_field_content": "ten_metric_components_plus_phi_plus_independent_chi",
            "first_order_dimension": 24,
            "all_direction_mechanism": "Frobenius_orthogonal_rotation_covariance_plus_complete_covector_coefficient_tensor_norm",
            "spectral_separation": "six_Riesz_contours_with_Neumann_resolvent_margin",
            "physical_reality_and_diagonalizability": "positive_ungauged_action_energy_form_on_transported_physical_Riesz_subspaces",
            "sufficient_binary64_certificate_not_universal_Horndeski_epsilon": True,
        },
        "reference_cluster_certificate": reference,
        "independent_exact_radial_crosscheck": _radial_symbol_crosscheck(),
        "witness_records": records,
        "witness_extrema": extrema,
        "strict_all_covector_margins_passed": True,
        "open_health_ball": {
            "strict_inequalities_and_continuous_coefficient_map_give_nonzero_open_state_neighborhoods": True,
            "every_nonzero_spatial_covector_covered_by_homogeneity": True,
            "finite_directional_sampling_used_for_candidate_backgrounds": False,
            "reference_contours_use_binary64_SVD_with_Lipschitz_chord_and_roundoff_allowance": True,
        },
        "epistemic_boundary": {
            "witnesses_do_not_enclose_the_complete_initial_hypersurfaces": True,
            "trajectory_containment_is_deferred_to_HLT1": True,
            "retained_EFT_remainder_or_frequency_control_supplied": False,
            "nonspherical_evolution_or_robustness_proven": False,
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
        "classification": FUTURE_EVIDENCE_CLASSIFICATIONS[
            "multidirectional_health"
        ],
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
            "all_nonzero_spatial_covectors_covered": True,
            "strong_hyperbolicity_margin_positive": True,
            "independent_chi_sector_included": True,
            "classical_principal_health_only_not_EFT_validity": True,
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
        raise ValueError("HYP2 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("HYP2 result must use canonical sorted JSON")
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
