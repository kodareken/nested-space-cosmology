#!/usr/bin/env python3
"""Regenerate the FGC-1-CTR1-REG1 regular-centre certificate."""

from __future__ import annotations

import argparse
from copy import deepcopy
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
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-ctr1-reg1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-ctr1-reg1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-ctr1-reg1.md"

from recursive_horizons.fgc.regular_center import (  # noqa: E402
    CERTIFIED_MAX_POWER,
    CERTIFIED_MIN_POWER,
    REGULAR_EQUATION_ORDER,
    REGULAR_PRIMITIVE_ORDER,
    SERIES_MAX_POWER,
    SERIES_MIN_POWER,
    TIME_LAYER_ORDER,
    analytic_defect_ledger,
    annular_regularized_evaluation,
    first_grid_point_convergence_certificate,
    regular_center_series_certificate,
    validate_regular_profile,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_COMMON_NONCLAIMS,
    FUTURE_EVIDENCE_CLASSIFICATIONS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-CTR1-REG1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "regular_center_formulation_verified"
EXPECTED_PROOF_CONTRACT = {
    "R_equals_r_times_A_and_v_equals_r_times_V",
    "reference_one_over_r_and_one_over_r_squared_differences_rederived",
    "complete_unredefined_ACT1_VAR1_residual_consumed",
    "complete_affected_REF1_extension_rederived",
    "all_removable_singular_limits_checked_in_exact_Laurent_arithmetic",
    "negative_power_guard_edge_fails_closed",
    "center_Taylor_system_serialized",
    "Minkowski_series_control_exact",
    "smooth_nontrivial_series_control_exact",
    "parity_and_elementary_flatness_mutations_fail_closed",
    "finite_curvature_series_required",
    "first_grid_point_convergence_required",
    "canonical_hash_bound_result_required",
    "no_initial_data_family_inferred",
    "no_holdout_execution",
}
PREDECESSORS = {
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "ref1_config": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
    "ref1_result": "results/fgc-1-hyp1-mhg-reference.json",
    "src1_config": "configs/fgc/fgc-1-src1-nl1.toml",
    "src1_result": "results/fgc-1-src1-nl1.json",
    "con4_config": "configs/fgc/fgc-1-con4-phy1.toml",
    "con4_result": "results/fgc-1-con4-phy1.json",
}
EXPECTED_RESULTS = {
    "action_result": ("FGC-1-ACT1", ("declared_action_and_scalar_algebra_certificate_passed",)),
    "variation_result": ("FGC-1-VAR1", ("var1_covariant_metric_equation_derived_and_cross_checked",)),
    "ref1_result": (
        "FGC-1-HYP1-MHG2-REF1",
        (
            "complete_reference_gauge_residual_passed",
            "exact_gauge_surface_equivalence_controls_passed",
            "mhg1_principal_gauge_block_regression_passed",
        ),
    ),
    "src1_result": ("FGC-1-SRC1-NL1", ("nonlinear_REF1_source_and_branch_solver_verified",)),
    "con4_result": ("FGC-1-CON4-PHY1", ("physical_gauge_reduction_constraint_system_closed",)),
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/regular_center.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/reference_connection.py",
        "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "src/recursive_horizons/fgc/scoped_run_authorization.py",
        "scripts/reproduce_fgc_ctr1_reg1.py",
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
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
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


def _profile_from_control(name: str, value: Mapping[str, Any]) -> dict[str, Any]:
    _keys(name, value, {"fields"})
    fields = value["fields"]
    if not isinstance(fields, Mapping) or set(fields) != set(REGULAR_PRIMITIVE_ORDER):
        raise ValueError(f"{name}.fields differ")
    output_fields: dict[str, dict[str, dict[int, Fraction]]] = {}
    expected_length = 1 if name == "minkowski" else 3
    for field in REGULAR_PRIMITIVE_ORDER:
        layers = fields[field]
        if not isinstance(layers, Mapping) or set(layers) != set(TIME_LAYER_ORDER):
            raise ValueError(f"{name}.{field} layers differ")
        output_fields[field] = {}
        for layer in TIME_LAYER_ORDER:
            raw = layers[layer]
            if not isinstance(raw, list) or len(raw) != expected_length:
                raise ValueError(f"{name}.{field}.{layer} coefficient count differs")
            coefficients = {
                2 * index: _fraction(f"{name}.{field}.{layer}[{index}]", item)
                for index, item in enumerate(raw)
            }
            output_fields[field][layer] = coefficients
    profile = {
        "profile_id": name,
        "fields": output_fields,
        "action": {
            "branch": "FGC-QR",
            "planck_mass": Q(2),
            "beta": -Q(1, 4),
            "mu": Q(3),
            "g4": Q(1, 2),
            "alpha": Q(0),
            "eta": Q(1, 2),
        },
    }
    validate_regular_profile(profile)
    return profile


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
            "series",
            "convergence",
            "negative_controls",
            "proof_contract",
            "controls",
            "claims",
        },
    )
    if {
        key: raw[key]
        for key in ("schema_version", "artifact_id", "project_version", "metric_signature", "riemann_convention")
    } != {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }:
        raise ValueError("CTR1-REG1 identity, version, or conventions differ")
    protocol_path = _path("protocol_config", raw["protocol_config"], "configs/fgc/fgc-2-sf1-protocol-v3.toml")
    run1_path = _path("run1_config", raw["run1_config"], "configs/fgc/fgc-1-run1-sym1.toml")
    paths = {name: _path(name, raw[name], expected) for name, expected in PREDECESSORS.items()}
    protocol = validate_sf1_protocol(tomllib.loads(protocol_path.read_text(encoding="utf-8")))

    scope = raw["scope"]
    _keys(
        "scope",
        scope,
        {
            "target_branch",
            "physical_equations",
            "gauge_formulation",
            "regular_variables",
            "shift_formula",
            "areal_radius_formula",
            "primitive_parity",
            "radial_derivative_parity",
            "elementary_flatness",
            "tilde_normal_factor",
            "hat_normal_factor",
            "affected_REF1_descendants_rederived",
            "existing_annular_reference_not_evaluated_at_r_equals_zero",
        },
    )
    if scope != {
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "gauge_formulation": "metric_defined_reference_modified_harmonic",
        "regular_variables": ["alpha", "V", "lambda", "A", "phi", "chi"],
        "shift_formula": "v=r*V",
        "areal_radius_formula": "R=r*A",
        "primitive_parity": "all_six_regular_variables_even",
        "radial_derivative_parity": "all_six_radial_derivatives_odd",
        "elementary_flatness": "A(t,0)=lambda(t,0)_through_two_time_derivatives",
        "tilde_normal_factor": 4,
        "hat_normal_factor": 9,
        "affected_REF1_descendants_rederived": True,
        "existing_annular_reference_not_evaluated_at_r_equals_zero": True,
    }:
        raise ValueError("CTR1-REG1 scope differs")

    series = raw["series"]
    _keys(
        "series",
        series,
        {
            "computed_minimum_power",
            "computed_maximum_power",
            "certified_minimum_power",
            "certified_maximum_power",
            "regular_equation_order",
            "curvature_controls",
            "require_no_negative_certified_powers",
            "require_even_regularized_equations",
            "require_even_regularized_gauge_components",
        },
    )
    if series != {
        "computed_minimum_power": SERIES_MIN_POWER,
        "computed_maximum_power": SERIES_MAX_POWER,
        "certified_minimum_power": CERTIFIED_MIN_POWER,
        "certified_maximum_power": CERTIFIED_MAX_POWER,
        "regular_equation_order": list(REGULAR_EQUATION_ORDER),
        "curvature_controls": ["Ricci_scalar", "Ricci_squared", "Kretschmann", "Gauss_Bonnet"],
        "require_no_negative_certified_powers": True,
        "require_even_regularized_equations": True,
        "require_even_regularized_gauge_components": True,
    }:
        raise ValueError("CTR1-REG1 series contract differs")

    convergence = raw["convergence"]
    _keys(
        "convergence",
        convergence,
        {"radii", "minimum_coarse_to_fine_error_ratio", "interpretation", "not_a_PDE_discretization_convergence_claim"},
    )
    radii = tuple(_fraction(f"convergence.radii[{index}]", value) for index, value in enumerate(convergence["radii"]))
    if (
        radii != (Q(1, 16), Q(1, 32), Q(1, 64))
        or _fraction("convergence.minimum_coarse_to_fine_error_ratio", convergence["minimum_coarse_to_fine_error_ratio"]) != Q(15, 4)
        or convergence["interpretation"] != "exact_first_positive_grid_point_approach_to_formal_center_limit"
        or convergence["not_a_PDE_discretization_convergence_claim"] is not True
    ):
        raise ValueError("CTR1-REG1 convergence contract differs")

    negative = raw["negative_controls"]
    _keys(
        "negative_controls",
        negative,
        {
            "conical_A0_replacement",
            "odd_phi_r_coefficient",
            "elementary_flatness_dt_A0_replacement",
            "require_nonzero_Ricci_r_minus_2_defect",
            "require_nonzero_box_phi_r_minus_1_defect",
            "all_invalid_profiles_fail_before_certificate",
        },
    )
    negative_values = {
        "conical_A0_replacement": _fraction("negative_controls.conical_A0_replacement", negative["conical_A0_replacement"]),
        "odd_phi_r_coefficient": _fraction("negative_controls.odd_phi_r_coefficient", negative["odd_phi_r_coefficient"]),
        "elementary_flatness_dt_A0_replacement": _fraction("negative_controls.elementary_flatness_dt_A0_replacement", negative["elementary_flatness_dt_A0_replacement"]),
    }
    if negative_values != {
        "conical_A0_replacement": Q(11, 10),
        "odd_phi_r_coefficient": Q(1, 17),
        "elementary_flatness_dt_A0_replacement": Q(1, 39),
    } or any(negative[key] is not True for key in (
        "require_nonzero_Ricci_r_minus_2_defect",
        "require_nonzero_box_phi_r_minus_1_defect",
        "all_invalid_profiles_fail_before_certificate",
    )):
        raise ValueError("CTR1-REG1 negative-control contract differs")
    proof_contract = raw["proof_contract"]
    if not isinstance(proof_contract, Mapping):
        raise ValueError("CTR1-REG1 proof contract must be a mapping")
    _keys("proof_contract", proof_contract, EXPECTED_PROOF_CONTRACT)
    if any(value is not True for value in proof_contract.values()):
        raise ValueError("every CTR1-REG1 proof-contract premise must remain true")
    if raw["claims"] != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("CTR1-REG1 claims must remain fail-closed")

    controls = raw["controls"]
    if not isinstance(controls, Mapping) or set(controls) != {"minkowski", "nontrivial"}:
        raise ValueError("CTR1-REG1 controls differ")
    profiles = {
        name: _profile_from_control(name, controls[name])
        for name in ("minkowski", "nontrivial")
    }

    records: dict[str, dict[str, Any]] = {}
    for result_name, (artifact_id, gates) in EXPECTED_RESULTS.items():
        record = _load_json(paths[result_name], result_name)
        _validate_predecessor(result_name, record, artifact_id=artifact_id, gates=gates)
        _binds(record, paths[result_name.removesuffix("_result") + "_config"], result_name)
        records[result_name] = record
    return {
        "raw": raw,
        "source_sha256": sha256(source).hexdigest(),
        "protocol_path": protocol_path,
        "run1_path": run1_path,
        "paths": paths,
        "records": records,
        "protocol": protocol,
        "profiles": profiles,
        "radii": radii,
        "negative_values": negative_values,
    }


def _expected_failure(profile: Mapping[str, Any], pattern: str) -> str:
    try:
        validate_regular_profile(profile)
    except (TypeError, ValueError) as exc:
        message = str(exc)
        if pattern not in message:
            raise ValueError(f"negative control failed for the wrong reason: {message}") from exc
        return message
    raise ValueError("invalid regular-centre profile was accepted")


def _quantitative_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    minkowski_profile = configuration["profiles"]["minkowski"]
    nontrivial_profile = configuration["profiles"]["nontrivial"]
    minkowski = regular_center_series_certificate(minkowski_profile)
    nontrivial = regular_center_series_certificate(nontrivial_profile)
    if any(minkowski["center_limits"]) or any(minkowski["gauge_center_limits"]):
        raise ValueError("Minkowski centre residual or gauge limit is nonzero")
    if any(minkowski["curvature_series"][name] for name in minkowski["curvature_series"]):
        raise ValueError("Minkowski centre curvature is nonzero")
    minkowski_annular = annular_regularized_evaluation(
        minkowski_profile,
        coordinate_radius=configuration["radii"][0],
        reference_radial_minimum=configuration["radii"][-1] / 4,
    )
    if any(minkowski_annular["regular_equations"]) or any(minkowski_annular["regular_gauge"]):
        raise ValueError("unchanged annular REF1 evaluator lost the exact Minkowski control")
    if not all(nontrivial["checks"].values()):
        raise ValueError("nontrivial regular-centre series check failed")
    if not any(nontrivial["curvature_series"]["Kretschmann"].values()):
        raise ValueError("nontrivial centre fixture did not exercise curvature")
    if nontrivial["center_limits"][2] != nontrivial["center_limits"][3]:
        raise ValueError("regular centre lost radial/angular metric isotropy")
    convergence = first_grid_point_convergence_certificate(
        nontrivial_profile,
        radii=configuration["radii"],
        formal_certificate=nontrivial,
    )

    conical = deepcopy(nontrivial_profile)
    conical["fields"]["A"]["value"][0] = configuration["negative_values"]["conical_A0_replacement"]
    conical_ledger = analytic_defect_ledger(conical)
    if conical_ledger["Ricci_scalar_r_minus_2_from_conical_defect"] == 0:
        raise ValueError("conical negative control did not produce its analytic curvature pole")
    conical_failure = _expected_failure(conical, "elementary flatness")

    odd_scalar = deepcopy(nontrivial_profile)
    odd_scalar["fields"]["phi"]["value"][1] = configuration["negative_values"]["odd_phi_r_coefficient"]
    odd_ledger = analytic_defect_ledger(odd_scalar)
    if odd_ledger["box_phi_r_minus_1_from_odd_scalar_term"] == 0:
        raise ValueError("odd-scalar negative control did not produce its analytic wave pole")
    odd_failure = _expected_failure(odd_scalar, "even centre parity")

    time_defect = deepcopy(nontrivial_profile)
    time_defect["fields"]["A"]["dt"][0] = configuration["negative_values"]["elementary_flatness_dt_A0_replacement"]
    time_failure = _expected_failure(time_defect, "elementary flatness")

    return {
        "Minkowski_control": {
            "formal_certificate": minkowski,
            "annular_radius": minkowski_annular["coordinate_radius"],
            "annular_regular_equations": minkowski_annular["regular_equations"],
            "annular_regular_gauge": minkowski_annular["regular_gauge"],
            "exact_zero_residual_gauge_and_curvature": True,
        },
        "smooth_nontrivial_control": {
            "formal_certificate": nontrivial,
            "center_radial_angular_metric_isotropy_exact": True,
            "nonzero_finite_Kretschmann_control": True,
            "first_grid_point_convergence": convergence,
        },
        "negative_controls": {
            "conical_defect": {
                "analytic_ledger": conical_ledger,
                "failure": conical_failure,
                "rejected_before_certificate": True,
            },
            "odd_scalar_parity": {
                "analytic_ledger": odd_ledger,
                "failure": odd_failure,
                "rejected_before_certificate": True,
            },
            "time_derivative_elementary_flatness": {
                "failure": time_failure,
                "rejected_before_certificate": True,
            },
        },
        "scope": {
            "affected_REF1_descendants_rederived": True,
            "annular_reference_evaluated_at_zero": False,
            "initial_data_family_constructed": False,
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
    return _serial(
        {
            "schema_version": 1,
            "artifact_id": ARTIFACT_ID,
            "project_version": PROJECT_VERSION,
            "classification": FUTURE_EVIDENCE_CLASSIFICATIONS["regular_center"],
            "generated_by": _rel(Path(__file__)),
            "derivation_document": _rel(OWNER_DOCUMENT),
            "derivation_document_sha256": _sha(OWNER_DOCUMENT),
            "source_config_sha256": {_rel(config_path): configuration["source_sha256"]},
            "predecessor_sha256": {_rel(path): _sha(path) for path in paths.values()},
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
                "R_equals_r_times_A_and_v_equals_r_times_V": True,
                "all_removable_singular_limits_derived_analytically": True,
                "negative_power_guard_edge_fails_closed": True,
                "parity_and_elementary_flatness_enforced": True,
                "finite_curvature_series_controls_passed": True,
                "first_grid_point_convergence_passed": True,
            },
            "gate_status": {REQUIRED_GATE: True},
            "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
        }
    )


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError("CTR1-REG1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("CTR1-REG1 result must use canonical sorted JSON")
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
