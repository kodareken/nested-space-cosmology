#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP1-MHG1 exact modified-harmonic certificate."""

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

DEFAULT_CONFIG = (
    REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-modified-harmonic.toml"
)
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-hyp1-modified-harmonic.json"
OWNER_DOCUMENT = REPOSITORY / "docs" / "fgc-hyp1-modified-harmonic.md"
IMPLEMENTATION_FILES = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "exact_linear_algebra.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_reduction.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_symbol.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "modified_harmonic.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "__init__.py",
    REPOSITORY / "scripts" / "reproduce_fgc_action.py",
    REPOSITORY / "scripts" / "reproduce_fgc_metric_variation.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_reduction.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_symbol.py",
    Path(__file__).resolve(),
)

from scripts.reproduce_fgc_action import load_config as load_action_config
from scripts.reproduce_fgc_metric_variation import load_config as load_variation_config
from scripts.reproduce_fgc_hyp1_reduction import load_config as load_reduction_config
from scripts.reproduce_fgc_hyp1_symbol import load_config as load_symbol_config
from recursive_horizons.fgc.modified_harmonic import modified_harmonic_certificate


@dataclass(frozen=True, slots=True)
class HYP1ModifiedHarmonicConfig:
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
    symbol_config: Path
    symbol_config_relative: str
    tilde_normal_factor: Fraction
    hat_normal_factor: Fraction
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


def _positive_integer(name: str, value: Any) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return Fraction(value)


def load_config(
    path: Path = DEFAULT_CONFIG,
) -> HYP1ModifiedHarmonicConfig:
    """Load the frozen MHG1 configuration with no implicit defaults."""

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
            "symbol_config",
            "formulation",
            "proof_contract",
            "open_gates",
        },
    )
    expected_scalars = {
        "schema_version": 1,
        "artifact_id": "FGC-1-HYP1-MHG1",
        "project_version": "0.11.0",
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }
    for name, expected in expected_scalars.items():
        if raw[name] != expected:
            raise ValueError(f"{name} must equal {expected}")

    path_fields = {
        "action_config": "configs/fgc/fgc-1-action-gate.toml",
        "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction_config": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol_config": "configs/fgc/fgc-1-hyp1-symbol.toml",
    }
    resolved_paths: dict[str, Path] = {}
    relative_paths: dict[str, str] = {}
    for name, expected in path_fields.items():
        resolved, relative = _repository_path(name, raw[name])
        if relative != expected:
            raise ValueError(f"{name} must name the frozen predecessor configuration")
        resolved_paths[name] = resolved
        relative_paths[name] = relative

    predecessors = (
        ("ACT1", load_action_config(resolved_paths["action_config"])),
        ("VAR1", load_variation_config(resolved_paths["variation_config"])),
        ("RED1", load_reduction_config(resolved_paths["reduction_config"])),
        ("SYM1", load_symbol_config(resolved_paths["symbol_config"])),
    )
    for name, predecessor in predecessors:
        if predecessor.project_version != raw["project_version"]:
            raise ValueError(f"{name} and MHG1 project versions differ")
        if predecessor.metric_signature != raw["metric_signature"]:
            raise ValueError(f"{name} and MHG1 metric signatures differ")
    for name, predecessor in predecessors[1:]:
        if predecessor.riemann_convention != raw["riemann_convention"]:
            raise ValueError(f"{name} and MHG1 Riemann conventions differ")

    formulation = _mapping("formulation", raw["formulation"])
    _exact_keys(
        "formulation",
        formulation,
        {
            "physical_equations",
            "gauge_condition",
            "gauge_fixing",
            "auxiliary_metric_family",
            "tilde_normal_factor",
            "hat_normal_factor",
            "normal_factor_order",
            "covector",
            "field_order",
            "equation_order",
            "first_order_state_order",
            "gauge_constraint_order",
            "equation_modification_on_gauge_surface",
        },
    )
    frozen_formulation = {
        "physical_equations": "unredefined_ACT1_VAR1",
        "gauge_condition": "H^mu=-tilde_g^rho_sigma*Gamma^mu_rho_sigma=0_up_to_declared_reference_connection",
        "gauge_fixing": "F*hat_P_alpha^beta_mu_nu*partial_beta_H^alpha",
        "auxiliary_metric_family": "g_aux^mu_nu=g^mu_nu-(q-1)*n^mu*n^nu",
        "normal_factor_order": "1<tilde<hat",
        "covector": "xi_A=(-c,1)",
        "equation_modification_on_gauge_surface": "none",
    }
    for name, expected in frozen_formulation.items():
        if formulation[name] != expected:
            raise ValueError(f"formulation.{name} is not frozen")
    tilde_factor = _positive_integer(
        "formulation.tilde_normal_factor", formulation["tilde_normal_factor"]
    )
    hat_factor = _positive_integer(
        "formulation.hat_normal_factor", formulation["hat_normal_factor"]
    )
    if not 1 < tilde_factor < hat_factor:
        raise ValueError("formulation requires 1 < tilde factor < hat factor")
    _string_list(
        "formulation.field_order",
        formulation["field_order"],
        ["h_tt", "h_tr", "h_rr", "areal_radius", "phi", "chi"],
    )
    _string_list(
        "formulation.equation_order",
        formulation["equation_order"],
        [
            "metric_tt_mhg",
            "metric_tr_mhg",
            "metric_rr_mhg",
            "metric_theta_theta_mhg",
            "scalar_phi",
            "scalar_chi",
        ],
    )
    _string_list(
        "formulation.first_order_state_order",
        formulation["first_order_state_order"],
        [
            "dt_h_tt",
            "dt_h_tr",
            "dt_h_rr",
            "dt_areal_radius",
            "dt_phi",
            "dt_chi",
            "dr_h_tt",
            "dr_h_tr",
            "dr_h_rr",
            "dr_areal_radius",
            "dr_phi",
            "dr_chi",
        ],
    )
    _string_list(
        "formulation.gauge_constraint_order",
        formulation["gauge_constraint_order"],
        ["H^t", "H^r"],
    )

    proof = _mapping("proof_contract", raw["proof_contract"])
    proof_keys = {
        "fixture_ids",
        "require_positive_effective_planck_coefficient",
        "require_nonzero_gauge_fixing_symbol",
        "require_exact_complete_determinant_factorization",
        "require_tilde_multiplicity",
        "require_hat_multiplicity",
        "require_auxiliary_cone_nesting_and_separation",
        "require_hat_cone_gauge_constraint_principal_propagation",
        "require_invertible_coordinate_time_kinetic_matrix",
        "require_exact_standard_first_order_reduction",
        "require_pointwise_complete_real_first_order_eigenbasis",
        "require_direct_x_zero_regularity_control",
        "require_sym1_physical_factor_regression",
    }
    _exact_keys("proof_contract", proof, proof_keys)
    fixture_ids = [
        "GR0_regular",
        "SGBL_constant_f",
        "FGCQR_schwarzschild_exterior",
        "FGCQR_activated_generic",
    ]
    _string_list("proof_contract.fixture_ids", proof["fixture_ids"], fixture_ids)
    for name in proof_keys - {
        "fixture_ids",
        "require_tilde_multiplicity",
        "require_hat_multiplicity",
    }:
        if proof[name] is not True:
            raise ValueError(f"proof_contract.{name} must be true")
    if proof["require_tilde_multiplicity"] != 2:
        raise ValueError("proof_contract.require_tilde_multiplicity must equal 2")
    if proof["require_hat_multiplicity"] != 2:
        raise ValueError("proof_contract.require_hat_multiplicity must equal 2")

    open_gates = _mapping("open_gates", raw["open_gates"])
    expected_open = {
        "complete_lower_order_first_order_sources_derived",
        "gauge_constraint_lower_order_coefficients_serialized",
        "reduction_constraint_propagation_proven",
        "uniform_open_domain_symmetrizer_proven",
        "open_retained_eft_domain_proven",
        "boundary_or_regular_center_system_derived",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
        "singularity_resolution_derived",
    }
    _exact_keys("open_gates", open_gates, expected_open)
    if any(value is not False for value in open_gates.values()):
        raise ValueError("every MHG1 open gate must remain false")

    reduction = predecessors[2][1]
    actual_fixture_ids = [
        fixture["fixture_id"] for fixture in reduction.configuration["fixtures"]
    ]
    if actual_fixture_ids != fixture_ids:
        raise ValueError("RED1 fixture order differs from MHG1")

    return HYP1ModifiedHarmonicConfig(
        schema_version=raw["schema_version"],
        artifact_id=raw["artifact_id"],
        project_version=raw["project_version"],
        metric_signature=raw["metric_signature"],
        riemann_convention=raw["riemann_convention"],
        action_config=resolved_paths["action_config"],
        action_config_relative=relative_paths["action_config"],
        variation_config=resolved_paths["variation_config"],
        variation_config_relative=relative_paths["variation_config"],
        reduction_config=resolved_paths["reduction_config"],
        reduction_config_relative=relative_paths["reduction_config"],
        symbol_config=resolved_paths["symbol_config"],
        symbol_config_relative=relative_paths["symbol_config"],
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
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
    certificate = modified_harmonic_certificate(
        {
            "fixtures": reduction.configuration["fixtures"],
            "tilde_normal_factor": config.tilde_normal_factor,
            "hat_normal_factor": config.hat_normal_factor,
        }
    )
    if not certificate["all_declared_exact_checks_pass"]:
        raise ValueError("MHG1 exact certificate failed a required proof check")
    serialized_certificate = _serialize_exact(certificate)
    source_paths = (
        config.action_config,
        config.variation_config,
        config.reduction_config,
        config.symbol_config,
        config_path.resolve(),
    )
    output = {
        "schema_version": config.schema_version,
        "artifact_id": config.artifact_id,
        "artifact_label": "FGC-1 HYP1 exact modified-harmonic spherical principal and gauge-constraint preflight",
        "project_version": config.project_version,
        "classification": certificate["classification"],
        "development_status": "complete_pointwise_mhg_principal_preflight_full_open_domain_hyp1_open",
        "generated_by": _relative(Path(__file__)),
        "source_configs": {
            "action": config.action_config_relative,
            "variation": config.variation_config_relative,
            "reduction": config.reduction_config_relative,
            "symbol": config.symbol_config_relative,
            "modified_harmonic": _relative(config_path),
        },
        "source_config_sha256": {
            _relative(path): _file_sha256(path) for path in source_paths
        },
        "derivation_document": _relative(OWNER_DOCUMENT),
        "derivation_document_sha256": _file_sha256(OWNER_DOCUMENT),
        "implementation_sha256": {
            _relative(path): _file_sha256(path) for path in IMPLEMENTATION_FILES
        },
        "formulation": _serialize_exact(config.configuration["formulation"]),
        "proof_contract": _serialize_exact(config.configuration["proof_contract"]),
        "input_fixture_ids": [
            fixture["fixture_id"] for fixture in reduction.configuration["fixtures"]
        ],
        "input_branch_order": [
            fixture["model_id"] for fixture in reduction.configuration["fixtures"]
        ],
        "modified_harmonic_certificate": serialized_certificate,
        "modified_harmonic_certificate_sha256": _object_sha256(
            serialized_certificate
        ),
        "gate_status": {
            "exact_pointwise_modified_harmonic_principal_gate_passed": True,
            "gauge_constraint_principal_propagation_gate_passed": True,
            "pointwise_complete_real_first_order_eigenbasis_passed": True,
            "direct_x_zero_regular_coefficient_gate_passed": True,
            "complete_lower_order_first_order_formulation_passed": False,
            "uniform_open_domain_hyperbolicity_gate_passed": False,
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
