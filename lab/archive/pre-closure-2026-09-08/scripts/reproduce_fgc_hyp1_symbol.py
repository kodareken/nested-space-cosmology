#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP1-SYM1 exact physical-symbol certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "src"))

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-symbol.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-hyp1-symbol.json"
ACTION_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-action-gate.toml"
VARIATION_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-metric-variation.toml"
REDUCTION_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-reduction.toml"
OWNER_DOCUMENT = REPOSITORY / "docs" / "fgc-hyp1-symbol.md"
IMPLEMENTATION_FILES = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "exact_linear_algebra.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_symbol.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_reduction.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "__init__.py",
    REPOSITORY / "scripts" / "reproduce_fgc_action.py",
    REPOSITORY / "scripts" / "reproduce_fgc_metric_variation.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_reduction.py",
    Path(__file__).resolve(),
)

from scripts.reproduce_fgc_action import load_config as load_action_config
from scripts.reproduce_fgc_metric_variation import load_config as load_variation_config
from scripts.reproduce_fgc_hyp1_reduction import load_config as load_reduction_config
from recursive_horizons.fgc.spherical_symbol import spherical_symbol_certificate


@dataclass(frozen=True, slots=True)
class HYP1SymbolConfig:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    riemann_convention: str
    action_config: Path
    action_config_relative: str
    variation_config: Path
    variation_config_relative: str
    reduction_config: Path
    reduction_config_relative: str
    configuration: dict[str, Any]
    source_sha256: str


def _mapping(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a table")
    return value


def _exact_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    actual = set(value)
    if actual != expected:
        raise ValueError(
            f"{name} keys differ: missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )


def _repository_path(name: str, value: Any) -> tuple[Path, str]:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a non-empty repository-relative path")
    relative = Path(value)
    if relative.is_absolute():
        raise ValueError(f"{name} must be repository relative")
    resolved = (REPOSITORY / relative).resolve()
    try:
        canonical = resolved.relative_to(REPOSITORY)
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside the repository") from exc
    if canonical != relative:
        raise ValueError(f"{name} must be canonical and traversal free")
    if not resolved.is_file():
        raise ValueError(f"{name} does not exist")
    return resolved, canonical.as_posix()


def _string_list(name: str, value: Any, expected: list[str]) -> None:
    if not isinstance(value, list) or value != expected:
        raise ValueError(f"{name} is not the frozen ordered string list")


def load_config(path: Path = DEFAULT_CONFIG) -> HYP1SymbolConfig:
    """Load the frozen SYM1 configuration with no implicit defaults."""

    path = path.resolve()
    source = path.read_bytes()
    raw = _mapping("root", tomllib.loads(source.decode("utf-8")))
    _exact_keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "action_config",
            "variation_config",
            "reduction_config",
            "symbol",
            "adm_preflight",
            "proof_contract",
            "fixture_expectations",
        },
    )
    if raw["schema_version"] != 1:
        raise ValueError("schema_version must equal 1")
    if raw["artifact_id"] != "FGC-1-HYP1-SYM1":
        raise ValueError("artifact_id must equal FGC-1-HYP1-SYM1")
    if raw["project_version"] != "0.11.0":
        raise ValueError("project_version must equal 0.11.0")
    if raw["metric_signature"] != "-+++":
        raise ValueError("metric_signature must equal -+++")
    if raw["riemann_convention"] != "plus_partial_mu_gamma_nu":
        raise ValueError("riemann_convention must equal plus_partial_mu_gamma_nu")

    action_path, action_relative = _repository_path("action_config", raw["action_config"])
    variation_path, variation_relative = _repository_path(
        "variation_config", raw["variation_config"]
    )
    reduction_path, reduction_relative = _repository_path(
        "reduction_config", raw["reduction_config"]
    )
    if action_relative != "configs/fgc/fgc-1-action-gate.toml":
        raise ValueError("action_config must name the frozen ACT1 configuration")
    if variation_relative != "configs/fgc/fgc-1-metric-variation.toml":
        raise ValueError("variation_config must name the frozen VAR1 configuration")
    if reduction_relative != "configs/fgc/fgc-1-hyp1-reduction.toml":
        raise ValueError("reduction_config must name the frozen RED1 configuration")

    action = load_action_config(action_path)
    variation = load_variation_config(variation_path)
    reduction = load_reduction_config(reduction_path)
    for name, version in (
        ("ACT1", action.project_version),
        ("VAR1", variation.project_version),
        ("RED1", reduction.project_version),
    ):
        if version != raw["project_version"]:
            raise ValueError(f"{name} and SYM1 project versions differ")
    for name, signature in (
        ("ACT1", action.metric_signature),
        ("VAR1", variation.metric_signature),
        ("RED1", reduction.metric_signature),
    ):
        if signature != raw["metric_signature"]:
            raise ValueError(f"{name} and SYM1 metric signatures differ")
    if variation.riemann_convention != raw["riemann_convention"]:
        raise ValueError("VAR1 and SYM1 Riemann conventions differ")
    if reduction.riemann_convention != raw["riemann_convention"]:
        raise ValueError("RED1 and SYM1 Riemann conventions differ")

    symbol = _mapping("symbol", raw["symbol"])
    _exact_keys(
        "symbol",
        symbol,
        {
            "covector",
            "base_field_order",
            "equation_order",
            "right_gauge_generator_order",
            "left_bianchi_generator_order",
            "quotient_variable_order",
            "quotient_rr_equation_order",
            "quotient_tt_equation_order",
            "metric_elimination_variable_order",
            "scalar_order",
            "metric_patch_singularity",
            "alternate_metric_patch",
            "physical_metric",
            "equation_modification",
        },
    )
    scalar_expected = {
        "covector": "xi_A=(-c,1)",
        "metric_patch_singularity": "c=-shift",
        "alternate_metric_patch": "metric_tt_row_covers_c_equals_minus_shift",
        "physical_metric": "g_ab_matter_coupled_metric",
        "equation_modification": "none_unredefined_covariant_equations",
    }
    for key, expected in scalar_expected.items():
        if symbol[key] != expected:
            raise ValueError(f"symbol.{key} is not frozen")
    _string_list(
        "symbol.base_field_order",
        symbol["base_field_order"],
        ["h_tt", "h_tr", "h_rr", "areal_radius", "phi", "chi"],
    )
    _string_list(
        "symbol.equation_order",
        symbol["equation_order"],
        [
            "metric_tt",
            "metric_tr",
            "metric_rr",
            "metric_theta_theta",
            "scalar_phi",
            "scalar_chi",
        ],
    )
    _string_list(
        "symbol.right_gauge_generator_order",
        symbol["right_gauge_generator_order"],
        ["X_t", "X_r"],
    )
    _string_list(
        "symbol.left_bianchi_generator_order",
        symbol["left_bianchi_generator_order"],
        ["b=t", "b=r"],
    )
    _string_list(
        "symbol.quotient_variable_order",
        symbol["quotient_variable_order"],
        ["h_tt_representative", "areal_radius", "phi", "chi"],
    )
    _string_list(
        "symbol.quotient_rr_equation_order",
        symbol["quotient_rr_equation_order"],
        ["metric_rr", "metric_theta_theta", "scalar_phi", "scalar_chi"],
    )
    _string_list(
        "symbol.quotient_tt_equation_order",
        symbol["quotient_tt_equation_order"],
        ["metric_tt", "metric_theta_theta", "scalar_phi", "scalar_chi"],
    )
    _string_list(
        "symbol.metric_elimination_variable_order",
        symbol["metric_elimination_variable_order"],
        ["h_tt_representative", "areal_radius"],
    )
    _string_list("symbol.scalar_order", symbol["scalar_order"], ["phi", "chi"])

    adm = _mapping("adm_preflight", raw["adm_preflight"])
    expected_adm = {
        "chart": "generalized_radial_adm_dynamic_lambda_dynamic_areal_radius",
        "normal_derivative": "partial_perp=partial_t-shift*partial_r",
        "extrinsic_curvature_sign": "K_ij=-one_half_L_n_gamma_ij",
        "normal_variable_order": [
            "A_perp",
            "A_r",
            "L_r",
            "B_perp",
            "B_r",
            "S_r",
            "P_phi_perp",
            "P_phi_r",
            "Q_phi_r",
            "P_chi_perp",
            "P_chi_r",
            "Q_chi_r",
        ],
        "candidate_constraint_order": ["E_nn", "E_ns"],
        "evolution_projection_order": ["E_ss", "E_theta_theta", "E_phi", "E_chi"],
        "kinetic_variable_order": [
            "A_perp",
            "B_perp",
            "P_phi_perp",
            "P_chi_perp",
        ],
        "gauge_source_second_jet_order": [
            "alpha.dtt",
            "alpha.dtr",
            "alpha.drr",
            "shift.dtt",
            "shift.dtr",
            "shift.drr",
        ],
        "naive_comparator": "prescribed_lapse_shift_unconstrained_dynamic_R_first_order_reduction",
    }
    if dict(adm) != expected_adm:
        raise ValueError("adm_preflight is not the frozen SYM1 preflight")

    proof = _mapping("proof_contract", raw["proof_contract"])
    expected_proof_keys = {
        "rational_encoding",
        "fixture_ids",
        "branch_order",
        "require_exact_right_gauge_null_identity",
        "require_exact_left_bianchi_null_identity",
        "require_metric_block_patch_factor",
        "require_metric_block_invertible_at_physical_roots",
        "require_chi_metric_null_factor",
        "require_quadratic_regulator_factor",
        "require_positive_pointwise_discriminants",
        "require_positive_exact_scalar_companion_symmetrizers",
        "require_candidate_constraints_free_of_normal_accelerations",
        "require_candidate_constraints_free_of_gauge_source_second_jets",
        "require_nonzero_pointwise_metric_phi_kinetic_determinant",
        "require_naive_gr_comparator_rejected",
    }
    _exact_keys("proof_contract", proof, expected_proof_keys)
    if proof["rational_encoding"] != "canonical_fraction_string_numerator_over_denominator":
        raise ValueError("proof_contract.rational_encoding is not frozen")
    fixture_ids = [
        "GR0_regular",
        "SGBL_constant_f",
        "FGCQR_schwarzschild_exterior",
        "FGCQR_activated_generic",
    ]
    branches = ["GR-0", "SGB-L", "FGC-QR", "FGC-QR"]
    _string_list("proof_contract.fixture_ids", proof["fixture_ids"], fixture_ids)
    _string_list("proof_contract.branch_order", proof["branch_order"], branches)
    for key in expected_proof_keys - {"rational_encoding", "fixture_ids", "branch_order"}:
        if proof[key] is not True:
            raise ValueError(f"proof_contract.{key} must be true")

    expectations = raw["fixture_expectations"]
    if not isinstance(expectations, list) or len(expectations) != 4:
        raise ValueError("fixture_expectations must contain four frozen entries")
    expected_relations = [
        "proportional_to_metric_null",
        "proportional_to_metric_null",
        "proportional_to_metric_null",
        "not_proportional_to_metric_null",
    ]
    expected_distinct = [False, False, False, True]
    for index, expectation in enumerate(expectations):
        item = _mapping(f"fixture_expectations[{index}]", expectation)
        _exact_keys(
            f"fixture_expectations[{index}]",
            item,
            {
                "fixture_id",
                "branch",
                "regulator_factor_relation",
                "regulator_cone_distinct_from_metric",
            },
        )
        if item["fixture_id"] != fixture_ids[index] or item["branch"] != branches[index]:
            raise ValueError("fixture expectation identity or branch order differs")
        if item["regulator_factor_relation"] != expected_relations[index]:
            raise ValueError("fixture regulator-factor expectation differs")
        if item["regulator_cone_distinct_from_metric"] is not expected_distinct[index]:
            raise ValueError("fixture cone-distinction expectation differs")

    reduction_ids = [fixture["fixture_id"] for fixture in reduction.configuration["fixtures"]]
    reduction_branches = [fixture["model_id"] for fixture in reduction.configuration["fixtures"]]
    if reduction_ids != fixture_ids or reduction_branches != branches:
        raise ValueError("RED1 fixtures do not match the frozen SYM1 proof contract")

    return HYP1SymbolConfig(
        schema_version=raw["schema_version"],
        artifact_id=raw["artifact_id"],
        project_version=raw["project_version"],
        metric_signature=raw["metric_signature"],
        riemann_convention=raw["riemann_convention"],
        action_config=action_path,
        action_config_relative=action_relative,
        variation_config=variation_path,
        variation_config_relative=variation_relative,
        reduction_config=reduction_path,
        reduction_config_relative=reduction_relative,
        configuration=dict(raw),
        source_sha256=sha256(source).hexdigest(),
    )


def _fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _serialize_exact(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if isinstance(value, Mapping):
        return {str(key): _serialize_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serialize_exact(item) for item in value]
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise ValueError(f"certificate contains unsupported value {type(value).__name__}")


def _relative(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _object_sha256(value: Any) -> str:
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return sha256(canonical.encode("utf-8")).hexdigest()


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    reduction = load_reduction_config(config.reduction_config)
    certificate = spherical_symbol_certificate(
        {"fixtures": reduction.configuration["fixtures"]}
    )
    if not certificate["all_declared_exact_checks_pass"]:
        raise ValueError("SYM1 exact certificate failed a required proof check")

    expectations = {
        item["fixture_id"]: item
        for item in config.configuration["fixture_expectations"]
    }
    for item in certificate["fixture_records"]:
        expected = expectations[item["fixture_id"]]
        distinct = not item["regulator_cone_proportional_to_metric"]
        if distinct is not expected["regulator_cone_distinct_from_metric"]:
            raise ValueError(f"{item['fixture_id']} cone relation differs from config")
        expected_relation = expected["regulator_factor_relation"]
        actual_relation = (
            "proportional_to_metric_null"
            if item["regulator_cone_proportional_to_metric"]
            else "not_proportional_to_metric_null"
        )
        if actual_relation != expected_relation:
            raise ValueError(f"{item['fixture_id']} factor relation differs from config")

    serialized_certificate = _serialize_exact(certificate)
    source_paths = (
        config.action_config,
        config.variation_config,
        config.reduction_config,
        config_path.resolve(),
    )
    output = {
        "schema_version": config.schema_version,
        "artifact_id": config.artifact_id,
        "artifact_label": "FGC-1 HYP1 exact covariant physical-symbol and ADM formulation preflight",
        "project_version": config.project_version,
        "classification": certificate["classification"],
        "development_status": "pointwise_physical_scalar_quotient_complete_full_hyp1_open",
        "generated_by": _relative(Path(__file__)),
        "source_configs": {
            "action": config.action_config_relative,
            "variation": config.variation_config_relative,
            "reduction": config.reduction_config_relative,
            "symbol": _relative(config_path),
        },
        "source_config_sha256": {
            _relative(path): _file_sha256(path) for path in source_paths
        },
        "derivation_document": _relative(OWNER_DOCUMENT),
        "derivation_document_sha256": _file_sha256(OWNER_DOCUMENT),
        "implementation_sha256": {
            _relative(path): _file_sha256(path) for path in IMPLEMENTATION_FILES
        },
        "formulation": _serialize_exact(config.configuration["symbol"]),
        "adm_preflight": _serialize_exact(config.configuration["adm_preflight"]),
        "proof_contract": _serialize_exact(config.configuration["proof_contract"]),
        "input_fixture_ids": [
            fixture["fixture_id"] for fixture in reduction.configuration["fixtures"]
        ],
        "input_branch_order": [
            fixture["model_id"] for fixture in reduction.configuration["fixtures"]
        ],
        "symbol_certificate": serialized_certificate,
        "symbol_certificate_sha256": _object_sha256(serialized_certificate),
        "gate_status": {
            "exact_covariant_physical_scalar_quotient_passed": True,
            "pointwise_scalar_companion_symmetrizers_passed": True,
            "adm_principal_constraint_kinetic_preflight_passed": True,
            "naive_fixed_gauge_comparator_rejected": True,
            "full_fgc1_hyp1_health_gate_passed": False,
            "evolution_authorized": False,
            "publication_local_defocusing_gate_passed": False,
        },
        "nonclaims": serialized_certificate["nonclaims"],
    }
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    payload = record(arguments.config)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
