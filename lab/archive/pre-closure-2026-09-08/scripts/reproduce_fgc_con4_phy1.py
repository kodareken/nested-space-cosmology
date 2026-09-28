#!/usr/bin/env python3
"""Regenerate the FGC-1-CON4-PHY1 constraint-closure certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-con4-phy1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-con4-phy1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-con4-phy1.md"

from recursive_horizons.fgc.constraint_system import (  # noqa: E402
    acceleration_independence_control,
    conditional_constraint_closure_statement,
    constraint_acceleration_structure_certificate,
    constraint_monitor_contract,
    gauge_constraint_characteristics,
    gauge_constraint_monitor_bundle,
    physical_constraint_term_certificate,
    physical_to_normal_gauge_map,
    reduction_constraint_monitor_bundle,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.modified_harmonic_reduction_subsidiary import (  # noqa: E402
    ReductionDifferentialFieldJet,
    ReductionDifferentialState,
)
from recursive_horizons.fgc.modified_harmonic_reference import (  # noqa: E402
    modified_harmonic_gauge_constraint,
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
    Jet2,
    SphericalState,
)
from scripts.reproduce_fgc_hyp1_fo1_rc1 import (  # noqa: E402
    load_config as load_fo1_config,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-CON4-PHY1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "physical_gauge_reduction_constraint_system_closed"
PREDECESSORS = {
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "compatible_point_config": "configs/fgc/fgc-1-hyp1-con1-comp1.toml",
    "compatible_point_result": "results/fgc-1-hyp1-con1-comp1.json",
    "metric_subsidiary_config": "configs/fgc/fgc-1-hyp1-con2-mprop1.toml",
    "metric_subsidiary_result": "results/fgc-1-hyp1-con2-mprop1.json",
    "gauge_cauchy_config": "configs/fgc/fgc-1-hyp1-con3-cau1.toml",
    "gauge_cauchy_result": "results/fgc-1-hyp1-con3-cau1.json",
}
EXPECTED_RESULTS = {
    "action_result": ("FGC-1-ACT1", ("declared_action_and_scalar_algebra_certificate_passed",)),
    "variation_result": ("FGC-1-VAR1", ("var1_covariant_metric_equation_derived_and_cross_checked",)),
    "compatible_point_result": (
        "FGC-1-HYP1-CON1-COMP1",
        ("exact_activated_local_compatible_constraint_datum_passed",),
    ),
    "metric_subsidiary_result": (
        "FGC-1-HYP1-CON2-MPROP1",
        ("local_differential_identity", "complete_kinematic_1plus1_reduction_subsidiary"),
    ),
    "gauge_cauchy_result": (
        "FGC-1-HYP1-CON3-CAU1",
        ("conditional_boundary_free_gauge_Cauchy_uniqueness_theorem_derived",),
    ),
}
EXPECTED_STATE_IDS = ("COMP1-compatible", "nonflat-shifted-A", "nonflat-shifted-B")
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/constraint_system.py",
        "src/recursive_horizons/fgc/modified_harmonic_constraints.py",
        "src/recursive_horizons/fgc/modified_harmonic_metric_propagation.py",
        "src/recursive_horizons/fgc/modified_harmonic_reduction_subsidiary.py",
        "src/recursive_horizons/fgc/modified_harmonic_cauchy.py",
        "scripts/reproduce_fgc_con4_phy1.py",
    )
)


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _path(name: str, value: Any, expected: str) -> Path:
    if not isinstance(value, str) or value != expected:
        raise ValueError(f"{name} must name {expected}")
    path = (REPOSITORY / value).resolve()
    try:
        relative = path.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside repository") from exc
    if relative != value or not path.is_file():
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return path


def _text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _fraction(name: str, value: Any) -> Fraction:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a canonical rational string")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a canonical rational string") from exc
    if _text(answer) != value:
        raise ValueError(f"{name} must be a canonical rational string")
    return answer


def _pairs(pairs):
    output = {}
    for key, value in pairs:
        if key in output:
            raise ValueError(f"duplicate JSON object key: {key}")
        output[key] = value
    return output


def _load_json(path: Path, name: str) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _text(value)
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"unsupported exact record value {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(_serial(value), indent=2, sort_keys=True, allow_nan=False) + "\n"


def _binds(record: Mapping[str, Any], config: Path, name: str) -> None:
    direct = record.get("source_config_sha256")
    plural = record.get("source_configs_sha256")
    expected = _sha(config)
    if isinstance(direct, str):
        passed = direct == expected
    elif isinstance(direct, Mapping):
        passed = direct.get(_rel(config)) == expected
    elif isinstance(plural, Mapping):
        passed = plural.get(_rel(config)) == expected
    else:
        passed = False
    if not passed:
        raise ValueError(f"{name} result is stale against its source config")


def _validate_predecessor(
    name: str,
    record: Mapping[str, Any],
    *,
    artifact_id: str,
    gates: tuple[str, ...],
) -> None:
    if record.get("artifact_id") != artifact_id or record.get("project_version") != PROJECT_VERSION:
        raise ValueError(f"{name} predecessor identity differs")
    status = record.get("gate_status")
    if not isinstance(status, Mapping) or any(status.get(gate) is not True for gate in gates):
        raise ValueError(f"{name} predecessor gate is not passed")
    nonclaims = record.get("nonclaims")
    if isinstance(nonclaims, Mapping) and any(value is not False for value in nonclaims.values()):
        raise ValueError(f"{name} predecessor promoted a nonclaim")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    source = path.resolve().read_bytes()
    raw = tomllib.loads(source.decode("utf-8"))
    _keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "protocol_config",
            "run1_config",
            *PREDECESSORS,
            "scope",
            "controls",
            "monitor_contract",
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
        raise ValueError("CON4-PHY1 identity, version, or conventions differ")

    protocol_path = _path(
        "protocol_config",
        raw["protocol_config"],
        "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    )
    run1_path = _path(
        "run1_config",
        raw["run1_config"],
        "configs/fgc/fgc-1-run1-sym1.toml",
    )
    paths = {name: _path(name, raw[name], expected) for name, expected in PREDECESSORS.items()}
    protocol_raw = tomllib.loads(protocol_path.read_text(encoding="utf-8"))
    protocol = validate_sf1_protocol(protocol_raw)

    scope = raw["scope"]
    _keys(
        "scope",
        scope,
        {
            "target_branch",
            "physical_equations",
            "reference_id",
            "radial_domain_minimum",
            "tilde_normal_factor",
            "hat_normal_factor",
            "boundary_free_domain_of_dependence_only",
            "smooth_solution_is_a_conditional_premise",
        },
    )
    if (
        scope["target_branch"] != "FGC-QR"
        or scope["physical_equations"] != "unredefined_ACT1_VAR1"
        or scope["reference_id"] != "flat_spherical_annulus"
        or _fraction("scope.radial_domain_minimum", scope["radial_domain_minimum"]) != Q(1, 2)
        or scope["tilde_normal_factor"] != 4
        or scope["hat_normal_factor"] != 9
        or scope["boundary_free_domain_of_dependence_only"] is not True
        or scope["smooth_solution_is_a_conditional_premise"] is not True
    ):
        raise ValueError("CON4-PHY1 scope differs")

    controls = raw["controls"]
    _keys(
        "controls",
        controls,
        {
            "exact_state_ids",
            "exact_state_count",
            "acceleration_directions_per_state",
            "expected_double_dual_H_support_count",
            "expected_double_dual_M_support_count",
            "expected_active_spherical_gauge_components",
            "expected_hat_root_signs_at_each_control",
            "on_shell_and_off_shell_reduction_monitor_controls",
        },
    )
    if (
        tuple(controls["exact_state_ids"]) != EXPECTED_STATE_IDS
        or controls["exact_state_count"] != 3
        or controls["acceleration_directions_per_state"] != 6
        or controls["expected_double_dual_H_support_count"] != 36
        or controls["expected_double_dual_M_support_count"] != 36
        or controls["expected_active_spherical_gauge_components"] != 2
        or controls["expected_hat_root_signs_at_each_control"] != [-1, 1]
        or controls["on_shell_and_off_shell_reduction_monitor_controls"] is not True
    ):
        raise ValueError("CON4-PHY1 exact control contract differs")

    monitor = raw["monitor_contract"]
    _keys(
        "monitor_contract",
        monitor,
        {
            "formula",
            "zero_scale_policy",
            "raw_value_and_denominator_recorded",
            "common_rescaling_invariant",
            "absolute_error_scale_supplied_here",
            "HLT1_and_NUM1_must_add_absolute_thresholds",
            "smallness_on_one_run_is_not_proof",
        },
    )
    expected_monitor = {
        "formula": "abs(sum_signed_terms)/sum_abs_signed_terms",
        "zero_scale_policy": "all_terms_exactly_zero_returns_zero_with_explicit_flag",
        "raw_value_and_denominator_recorded": True,
        "common_rescaling_invariant": True,
        "absolute_error_scale_supplied_here": False,
        "HLT1_and_NUM1_must_add_absolute_thresholds": True,
        "smallness_on_one_run_is_not_proof": True,
    }
    if monitor != expected_monitor:
        raise ValueError("CON4-PHY1 monitor contract differs")
    if not isinstance(raw["proof_contract"], Mapping) or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("every CON4-PHY1 proof-contract premise must remain true")
    if raw["claims"] != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("CON4-PHY1 claims must remain fail-closed")

    records: dict[str, dict[str, Any]] = {}
    for result_name, (artifact_id, gates) in EXPECTED_RESULTS.items():
        record = _load_json(paths[result_name], result_name)
        _validate_predecessor(result_name, record, artifact_id=artifact_id, gates=gates)
        config_name = result_name.removesuffix("_result") + "_config"
        _binds(record, paths[config_name], result_name)
        records[result_name] = record

    return {
        "raw": raw,
        "source_sha256": sha256(source).hexdigest(),
        "protocol_path": protocol_path,
        "run1_path": run1_path,
        "paths": paths,
        "records": records,
        "protocol": protocol,
    }


def _control_states() -> tuple[tuple[str, SphericalState, Fraction], ...]:
    fixture = load_fo1_config().flat_fixture
    compatible = activated_compatible_state(fixture)["state"]
    state_a = SphericalState(
        h_tt=Jet2(-Q(3, 2), Q(1, 13), -Q(1, 17), Q(1, 19), Q(1, 23), -Q(1, 29)),
        h_tr=Jet2(Q(1, 5), -Q(1, 31), Q(1, 37), Q(1, 41), -Q(1, 43), Q(1, 47)),
        h_rr=Jet2(Q(4, 3), Q(1, 53), -Q(1, 59), Q(1, 61), Q(1, 67), -Q(1, 71)),
        areal_radius=Jet2(Q(5, 2), -Q(1, 73), Q(1, 79), Q(1, 83), -Q(1, 89), Q(1, 97)),
        phi=Jet2(Q(1, 7), Q(1, 101), -Q(1, 103), Q(1, 107), Q(1, 109), -Q(1, 113)),
        chi=Jet2(-Q(1, 11), -Q(1, 127), Q(1, 131), Q(1, 137), -Q(1, 139), Q(1, 149)),
        planck_mass=Q(1), beta=Q(1, 3), mu=Q(1), g4=Q(1), eta=Q(1, 5), branch="FGC-QR",
    )
    state_b = SphericalState(
        h_tt=Jet2(-Q(5, 4), -Q(1, 17), Q(1, 19), -Q(1, 23), Q(1, 29), Q(1, 31)),
        h_tr=Jet2(-Q(1, 6), Q(1, 37), -Q(1, 41), Q(1, 43), Q(1, 47), -Q(1, 53)),
        h_rr=Jet2(Q(3, 2), -Q(1, 59), Q(1, 61), Q(1, 67), -Q(1, 71), Q(1, 73)),
        areal_radius=Jet2(Q(7, 3), Q(1, 79), -Q(1, 83), -Q(1, 89), Q(1, 97), Q(1, 101)),
        phi=Jet2(-Q(1, 9), Q(1, 103), Q(1, 107), -Q(1, 109), Q(1, 113), Q(1, 127)),
        chi=Jet2(Q(2, 13), -Q(1, 131), -Q(1, 137), Q(1, 139), Q(1, 149), -Q(1, 151)),
        planck_mass=Q(1), beta=-Q(1, 5), mu=Q(1), g4=Q(3, 2), eta=Q(1, 7), branch="FGC-QR",
    )
    return (
        (EXPECTED_STATE_IDS[0], compatible, Q(4)),
        (EXPECTED_STATE_IDS[1], state_a, Q(5, 2)),
        (EXPECTED_STATE_IDS[2], state_b, Q(3)),
    )


def _monitor_summary(bundle: Mapping[str, Any]) -> dict[str, Any]:
    monitors = [item["monitor"] for item in bundle["constraint_C"]]
    monitors.extend(item["monitor"] for item in bundle["covariant_derivative_nabla_C"])
    return {
        "monitor_count": len(monitors),
        "maximum_relative_cancellation": max(item["relative_cancellation"] for item in monitors),
        "nonzero_raw_value_count": sum(item["raw_value"] != 0 for item in monitors),
        "zero_scale_convention_count": sum(item["zero_scale_convention_applied"] for item in monitors),
        "all_values_in_closed_unit_interval": all(Q(0) <= item["relative_cancellation"] <= Q(1) for item in monitors),
        "all_components_recomposed_exactly": bundle["all_components_recomposed_exactly"],
    }


def _reduction_state(*, on_shell: bool) -> ReductionDifferentialState:
    jet = (
        ReductionDifferentialFieldJet(Q(2), Q(2), Q(3), Q(3), Q(5), Q(5), Q(5))
        if on_shell
        else ReductionDifferentialFieldJet(Q(2), Q(1), Q(3), Q(1), Q(5), Q(2), Q(7))
    )
    return ReductionDifferentialState(**{field: jet for field in ("h_tt", "h_tr", "h_rr", "areal_radius", "phi", "chi")})


def _reduction_summary(*, on_shell: bool) -> dict[str, Any]:
    bundle = reduction_constraint_monitor_bundle(_reduction_state(on_shell=on_shell))
    records = bundle["records"]
    return {
        "on_kinematic_shell": on_shell,
        "D_values": tuple(item["D_equals_partial_t_u_minus_p"]["raw_value"] for item in records),
        "C_values": tuple(item["C_equals_q_minus_partial_r_u"]["raw_value"] for item in records),
        "K_values": tuple(item["K_equals_partial_t_q_minus_partial_r_p"]["raw_value"] for item in records),
        "all_subsidiary_identity_residuals_zero": all(
            item["subsidiary_identity_d_t_C_plus_d_r_D_minus_K"]["raw_value"] == 0
            for item in records
        ),
        "zero_radial_principal_speed": bundle["subsidiary_radial_principal_speed"] == 0,
        "no_incoming_reduction_fields": bundle["incoming_reduction_constraint_fields_at_either_boundary"] == (),
    }


def _quantitative_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    states = _control_states()
    structure = constraint_acceleration_structure_certificate()
    controls = configuration["raw"]["controls"]
    if (
        structure["Gauss_Bonnet_sector"]["Hamiltonian_support_count"]
        != controls["expected_double_dual_H_support_count"]
        or structure["Gauss_Bonnet_sector"]["momentum_support_count_for_one_spatial_direction"]
        != controls["expected_double_dual_M_support_count"]
    ):
        raise ValueError("double-dual support count differs from the frozen control")
    acceleration = acceleration_independence_control(tuple(state for _, state, _ in states))
    if acceleration["state_count"] != controls["exact_state_count"]:
        raise ValueError("acceleration-control state count differs")

    state_records: list[dict[str, Any]] = []
    reference = flat_spherical_annulus_reference(radial_domain_minimum=Q(1, 2))
    for state_id, state, radius in states:
        gauge = modified_harmonic_gauge_constraint(
            state,
            reference=reference,
            coordinate_radius=radius,
            tilde_normal_factor=Q(4),
        )
        physical = physical_constraint_term_certificate(state)
        normal_map = physical_to_normal_gauge_map(state, gauge, hat_normal_factor=Q(9))
        characteristics = gauge_constraint_characteristics(state, hat_normal_factor=Q(9))
        if list(characteristics["root_signs"]) != controls["expected_hat_root_signs_at_each_control"]:
            raise ValueError(f"{state_id} hat root sign split differs")
        state_records.append(
            {
                "state_id": state_id,
                "coordinate_radius": radius,
                "shift": physical["normalization"]["shift_s"],
                "lapse_squared": physical["normalization"]["lapse_squared"],
                "h_rr": physical["normalization"]["radial_spatial_metric_h_rr"],
                "physical_constraints": {
                    "term_order": physical["term_order"],
                    "H": physical["Hamiltonian"],
                    "M": physical["momentum"],
                    "H_terms": physical["Hamiltonian_terms"],
                    "M_terms": physical["momentum_terms"],
                    "H_relative_cancellation": physical["Hamiltonian_monitor"]["relative_cancellation"],
                    "M_relative_cancellation": physical["momentum_monitor"]["relative_cancellation"],
                    "exact_recomposition": physical["Hamiltonian_exact_recomposition"] and physical["momentum_exact_recomposition"],
                },
                "normal_gauge_map": {
                    "matrix": normal_map["computed_matrix"],
                    "analytic_matrix_equal": normal_map["direct_and_analytic_matrices_equal"],
                    "determinant": normal_map["determinant"],
                    "determinant_strictly_negative": normal_map["determinant_strictly_negative"],
                    "H_M_zero_iff_normal_nabla_C_zero": normal_map["normal_derivative_equivalence"]["H_and_M_zero_iff_normal_nabla_C_zero_on_REF1_shell"],
                },
                "gauge_characteristics": {
                    "hat_null_polynomial": characteristics["hat_null_polynomial"],
                    "root_isolations": characteristics["root_isolations"],
                    "root_signs": characteristics["root_signs"],
                    "outer_incoming_component_count": characteristics["outer_incoming_component_count_at_control"],
                    "inner_incoming_component_count": characteristics["inner_incoming_component_count_at_control"],
                    "constraint_preserving_boundary_map_derived": characteristics["constraint_preserving_boundary_map_to_main_fields_derived"],
                },
                "gauge_monitor_summary": _monitor_summary(gauge_constraint_monitor_bundle(gauge)),
            }
        )

    con2 = configuration["records"]["metric_subsidiary_result"]
    con3 = configuration["records"]["gauge_cauchy_result"]
    closure = conditional_constraint_closure_statement()
    monitor = constraint_monitor_contract()
    reduction_on = _reduction_summary(on_shell=True)
    reduction_off = _reduction_summary(on_shell=False)
    if not reduction_on["all_subsidiary_identity_residuals_zero"] or not reduction_off["all_subsidiary_identity_residuals_zero"]:
        raise ValueError("reduction monitor controls lost the off-shell identity")
    if any(reduction_on[name] != (Q(0),) * 6 for name in ("D_values", "C_values", "K_values")):
        raise ValueError("on-shell reduction control is not on shell")
    if not any(value != 0 for name in ("D_values", "C_values", "K_values") for value in reduction_off[name]):
        raise ValueError("off-shell reduction control did not exercise a nonzero residual")

    return {
        "acceleration_structure": structure,
        "exact_acceleration_controls": acceleration,
        "state_controls": tuple(state_records),
        "predecessor_composition": {
            "CON2_artifact_id": con2["artifact_id"],
            "CON2_local_metric_derived_identity_passed": con2["gate_status"]["local_differential_identity"],
            "CON2_complete_kinematic_reduction_subsidiary_passed": con2["gate_status"]["complete_kinematic_1plus1_reduction_subsidiary"],
            "CON3_artifact_id": con3["artifact_id"],
            "CON3_conditional_boundary_free_gauge_uniqueness_passed": con3["gate_status"]["conditional_boundary_free_gauge_Cauchy_uniqueness_theorem_derived"],
            "conditional_closure": closure,
        },
        "reduction_monitor_controls": {
            "on_shell": reduction_on,
            "off_shell": reduction_off,
        },
        "monitor_contract": monitor,
        "scope": {
            "boundary_free_domain_of_dependence_only": True,
            "smooth_full_REF1_plus_two_scalar_solution_is_conditional": True,
            "constraint_preserving_boundary_map_deferred_to_BND2": True,
            "compatible_nonzero_width_slice_deferred_to_ID1": True,
            "holdout_executed": False,
        },
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(config_path)
    quantitative = _quantitative_certificate(configuration)
    paths = configuration["paths"]
    protocol_hash = configuration["protocol"]["semantic_holdout_contract_sha256"]
    scope = {
        "run1_artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "run1_config_sha256": _sha(configuration["run1_path"]),
        "protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "protocol_config_sha256": _sha(configuration["protocol_path"]),
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
    return _serial({
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": FUTURE_EVIDENCE_CLASSIFICATIONS["physical_constraints"],
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): configuration["source_sha256"]},
        "predecessor_sha256": {
            _rel(path): _sha(path)
            for path in paths.values()
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
            "Hamiltonian_and_momentum_projections_derived": True,
            "metric_defined_gauge_constraints_closed": True,
            "kinematic_reduction_constraints_closed": True,
            "subsidiary_system_closed": True,
            "canonical_normalized_monitors_defined": True,
            "smallness_on_one_run_not_used_as_proof": True,
        },
        "gate_status": {REQUIRED_GATE: True},
        "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
    })


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError("CON4-PHY1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("CON4-PHY1 result must use canonical sorted JSON")
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
