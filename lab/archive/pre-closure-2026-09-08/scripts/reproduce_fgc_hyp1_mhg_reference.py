#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP1-MHG2-REF1 reference-gauge certificate."""

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

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-mhg-reference.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-hyp1-mhg-reference.json"
OWNER_DOCUMENT = REPOSITORY / "docs" / "fgc-hyp1-mhg-reference.md"
IMPLEMENTATION_FILES = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "exact_linear_algebra.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_reduction.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_symbol.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "modified_harmonic.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "reference_connection.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "modified_harmonic_reference.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "__init__.py",
    REPOSITORY / "scripts" / "reproduce_fgc_action.py",
    REPOSITORY / "scripts" / "reproduce_fgc_metric_variation.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_reduction.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_symbol.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_modified_harmonic.py",
    Path(__file__).resolve(),
)

from scripts.reproduce_fgc_action import load_config as load_action_config
from scripts.reproduce_fgc_metric_variation import load_config as load_variation_config
from scripts.reproduce_fgc_hyp1_reduction import load_config as load_reduction_config
from scripts.reproduce_fgc_hyp1_symbol import load_config as load_symbol_config
from scripts.reproduce_fgc_hyp1_modified_harmonic import (
    load_config as load_modified_harmonic_config,
)
from recursive_horizons.fgc.modified_harmonic_reference import (
    modified_harmonic_reference_certificate,
    required_ref1_nonclaims,
)
from recursive_horizons.fgc.reference_connection import (
    REFERENCE_COORDINATE_ORDER,
    SPHERICAL_FLAT_REFERENCE_ID,
    flat_spherical_annulus_reference,
)


@dataclass(frozen=True, slots=True)
class HYP1MHGReferenceConfig:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    riemann_convention: str
    source_paths: dict[str, Path]
    source_relatives: dict[str, str]
    radial_domain_minimum: Fraction
    tilde_normal_factor: Fraction
    hat_normal_factor: Fraction
    coordinate_radii: dict[str, Fraction]
    fixture_ids: tuple[str, ...]
    gauge_surface_fixture_ids: tuple[str, ...]
    off_gauge_fixture_ids: tuple[str, ...]
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


def _fraction(name: str, value: Any, *, positive: bool = False) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a canonical rational string")
    if isinstance(value, int):
        result = Fraction(value)
    elif isinstance(value, str) and value:
        try:
            result = Fraction(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{name} must be a canonical rational string") from exc
        canonical = (
            str(result.numerator)
            if result.denominator == 1
            else f"{result.numerator}/{result.denominator}"
        )
        if value != canonical:
            raise ValueError(f"{name} must use canonical rational encoding")
    else:
        raise ValueError(f"{name} must be an exact rational value")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _positive_integer(name: str, value: Any) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return Fraction(value)


def _string_list(name: str, value: Any, expected: list[str]) -> None:
    if not isinstance(value, list) or value != expected:
        raise ValueError(f"{name} is not the frozen ordered string list")


def _repository_path(name: str, value: Any) -> tuple[Path, str]:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a repository-relative path")
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


def load_config(path: Path = DEFAULT_CONFIG) -> HYP1MHGReferenceConfig:
    """Load the frozen REF1 configuration and reconcile every predecessor."""

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
            "modified_harmonic_config",
            "reference",
            "formulation",
            "coordinate_radii",
            "proof_contract",
            "open_gates",
        },
    )
    expected_scalars = {
        "schema_version": 1,
        "artifact_id": "FGC-1-HYP1-MHG2-REF1",
        "project_version": "0.11.0",
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }
    for name, expected in expected_scalars.items():
        if raw[name] != expected:
            raise ValueError(f"{name} must equal {expected}")

    expected_paths = {
        "action_config": "configs/fgc/fgc-1-action-gate.toml",
        "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction_config": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol_config": "configs/fgc/fgc-1-hyp1-symbol.toml",
        "modified_harmonic_config": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
    }
    source_paths: dict[str, Path] = {}
    source_relatives: dict[str, str] = {}
    for name, expected in expected_paths.items():
        resolved, relative = _repository_path(name, raw[name])
        if relative != expected:
            raise ValueError(f"{name} must name the frozen predecessor")
        source_paths[name] = resolved
        source_relatives[name] = relative

    predecessors = (
        ("ACT1", load_action_config(source_paths["action_config"])),
        ("VAR1", load_variation_config(source_paths["variation_config"])),
        ("RED1", load_reduction_config(source_paths["reduction_config"])),
        ("SYM1", load_symbol_config(source_paths["symbol_config"])),
        (
            "MHG1",
            load_modified_harmonic_config(source_paths["modified_harmonic_config"]),
        ),
    )
    for name, predecessor in predecessors:
        if predecessor.project_version != raw["project_version"]:
            raise ValueError(f"{name} and REF1 project versions differ")
        if predecessor.metric_signature != raw["metric_signature"]:
            raise ValueError(f"{name} and REF1 metric signatures differ")
    for name, predecessor in predecessors[1:]:
        if predecessor.riemann_convention != raw["riemann_convention"]:
            raise ValueError(f"{name} and REF1 Riemann conventions differ")

    reference = _mapping("reference", raw["reference"])
    _exact_keys(
        "reference",
        reference,
        {
            "reference_id",
            "coordinate_order",
            "angular_evaluation",
            "radial_domain_minimum",
            "prescribed_gauge_source",
            "center_included",
        },
    )
    if reference["reference_id"] != SPHERICAL_FLAT_REFERENCE_ID:
        raise ValueError("reference.reference_id is not frozen")
    _string_list(
        "reference.coordinate_order",
        reference["coordinate_order"],
        list(REFERENCE_COORDINATE_ORDER),
    )
    if reference["angular_evaluation"] != "theta=pi/2":
        raise ValueError("reference.angular_evaluation is not frozen")
    if reference["prescribed_gauge_source"] != "zero":
        raise ValueError("reference.prescribed_gauge_source is not frozen")
    if reference["center_included"] is not False:
        raise ValueError("REF1 must not claim that its spherical annulus includes the center")
    radial_minimum = _fraction(
        "reference.radial_domain_minimum",
        reference["radial_domain_minimum"],
        positive=True,
    )

    formulation = _mapping("formulation", raw["formulation"])
    _exact_keys(
        "formulation",
        formulation,
        {
            "physical_equations",
            "gauge_constraint",
            "gauge_fixing",
            "auxiliary_metric_family",
            "tilde_normal_factor",
            "hat_normal_factor",
            "normal_factor_order",
            "field_order",
            "equation_order",
            "gauge_constraint_order",
            "scalar_equations",
            "equation_modification_on_gauge_surface",
        },
    )
    frozen_formulation = {
        "physical_equations": "unredefined_ACT1_VAR1",
        "gauge_constraint": "C^mu=-tilde_g^rho_sigma*(Gamma^mu_rho_sigma-bar_Gamma^mu_rho_sigma)",
        "gauge_fixing": "F*hat_P_alpha^beta_mu_nu*nabla_beta_C^alpha",
        "auxiliary_metric_family": "g_aux^mu_nu=g^mu_nu-(q-1)*n^mu*n^nu",
        "normal_factor_order": "1<tilde<hat",
        "scalar_equations": "unmodified",
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
    mhg = predecessors[4][1]
    if (
        tilde_factor != mhg.tilde_normal_factor
        or hat_factor != mhg.hat_normal_factor
    ):
        raise ValueError("REF1 auxiliary factors differ from MHG1")
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
        "formulation.gauge_constraint_order",
        formulation["gauge_constraint_order"],
        ["C^t", "C^r", "C^theta", "C^phi"],
    )

    proof = _mapping("proof_contract", raw["proof_contract"])
    proof_keys = {
        "fixture_ids",
        "gauge_surface_fixture_ids",
        "off_gauge_fixture_ids",
        "require_positive_reference_annulus",
        "require_positive_effective_planck_coefficient",
        "require_exact_gauge_surface_extension_zero",
        "require_nonzero_off_gauge_extension",
        "require_exact_full_residual_decomposition",
        "require_unmodified_scalar_equations",
        "require_symmetric_spherical_extension",
        "require_exact_mhg1_principal_regression",
        "require_schwarzschild_unredefined_solution_control",
        "require_spherical_contravariant_gauge_components",
    }
    _exact_keys("proof_contract", proof, proof_keys)
    fixture_ids = (
        "GR0_regular",
        "SGBL_constant_f",
        "FGCQR_schwarzschild_exterior",
        "FGCQR_activated_generic",
    )
    gauge_surface_ids = ("GR0_regular", "SGBL_constant_f")
    off_gauge_ids = (
        "FGCQR_schwarzschild_exterior",
        "FGCQR_activated_generic",
    )
    _string_list("proof_contract.fixture_ids", proof["fixture_ids"], list(fixture_ids))
    _string_list(
        "proof_contract.gauge_surface_fixture_ids",
        proof["gauge_surface_fixture_ids"],
        list(gauge_surface_ids),
    )
    _string_list(
        "proof_contract.off_gauge_fixture_ids",
        proof["off_gauge_fixture_ids"],
        list(off_gauge_ids),
    )
    for name in proof_keys - {
        "fixture_ids",
        "gauge_surface_fixture_ids",
        "off_gauge_fixture_ids",
    }:
        if proof[name] is not True:
            raise ValueError(f"proof_contract.{name} must be true")

    radii = _mapping("coordinate_radii", raw["coordinate_radii"])
    _exact_keys("coordinate_radii", radii, set(fixture_ids))
    coordinate_radii = {
        fixture_id: _fraction(
            f"coordinate_radii.{fixture_id}", radii[fixture_id], positive=True
        )
        for fixture_id in fixture_ids
    }
    if any(value <= radial_minimum for value in coordinate_radii.values()):
        raise ValueError("every coordinate radius must lie strictly inside the annulus")

    open_gates = _mapping("open_gates", raw["open_gates"])
    _exact_keys("open_gates", open_gates, set(required_ref1_nonclaims()))
    if any(value is not False for value in open_gates.values()):
        raise ValueError("every REF1 open gate must remain false")

    reduction = predecessors[2][1]
    actual_fixture_ids = tuple(
        fixture["fixture_id"] for fixture in reduction.configuration["fixtures"]
    )
    if actual_fixture_ids != fixture_ids:
        raise ValueError("RED1 fixture order differs from REF1")

    return HYP1MHGReferenceConfig(
        schema_version=raw["schema_version"],
        artifact_id=raw["artifact_id"],
        project_version=raw["project_version"],
        metric_signature=raw["metric_signature"],
        riemann_convention=raw["riemann_convention"],
        source_paths=source_paths,
        source_relatives=source_relatives,
        radial_domain_minimum=radial_minimum,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
        coordinate_radii=coordinate_radii,
        fixture_ids=fixture_ids,
        gauge_surface_fixture_ids=gauge_surface_ids,
        off_gauge_fixture_ids=off_gauge_ids,
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
    reduction = load_reduction_config(config.source_paths["reduction_config"])
    reference = flat_spherical_annulus_reference(
        radial_domain_minimum=config.radial_domain_minimum
    )
    certificate = modified_harmonic_reference_certificate(
        {
            "fixtures": reduction.configuration["fixtures"],
            "coordinate_radii": config.coordinate_radii,
            "reference": reference,
            "tilde_normal_factor": config.tilde_normal_factor,
            "hat_normal_factor": config.hat_normal_factor,
            "gauge_surface_fixture_ids": config.gauge_surface_fixture_ids,
            "off_gauge_fixture_ids": config.off_gauge_fixture_ids,
        }
    )
    if not certificate["all_declared_exact_checks_pass"]:
        raise ValueError("REF1 exact certificate failed a required proof check")
    serialized_certificate = _serialize_exact(certificate)
    source_paths = tuple(config.source_paths.values()) + (config_path.resolve(),)
    return {
        "schema_version": config.schema_version,
        "artifact_id": config.artifact_id,
        "artifact_label": "FGC-1 HYP1 exact full reference-connection modified-harmonic residual gate",
        "project_version": config.project_version,
        "classification": certificate["classification"],
        "development_status": "complete_reference_gauge_residual_implicit_propagation_and_domain_gates_open",
        "generated_by": _relative(Path(__file__)),
        "source_configs": {
            "action": config.source_relatives["action_config"],
            "variation": config.source_relatives["variation_config"],
            "reduction": config.source_relatives["reduction_config"],
            "symbol": config.source_relatives["symbol_config"],
            "modified_harmonic": config.source_relatives[
                "modified_harmonic_config"
            ],
            "reference": _relative(config_path),
        },
        "source_config_sha256": {
            _relative(path): _file_sha256(path) for path in source_paths
        },
        "derivation_document": _relative(OWNER_DOCUMENT),
        "derivation_document_sha256": _file_sha256(OWNER_DOCUMENT),
        "implementation_sha256": {
            _relative(path): _file_sha256(path) for path in IMPLEMENTATION_FILES
        },
        "reference": _serialize_exact(config.configuration["reference"]),
        "formulation": _serialize_exact(config.configuration["formulation"]),
        "proof_contract": _serialize_exact(config.configuration["proof_contract"]),
        "input_fixture_ids": list(config.fixture_ids),
        "input_branch_order": [
            fixture["model_id"] for fixture in reduction.configuration["fixtures"]
        ],
        "reference_gauge_certificate": serialized_certificate,
        "reference_gauge_certificate_sha256": _object_sha256(
            serialized_certificate
        ),
        "gate_status": {
            "complete_reference_gauge_residual_passed": True,
            "exact_gauge_surface_equivalence_controls_passed": True,
            "mhg1_principal_gauge_block_regression_passed": True,
            "implicit_second_time_derivative_branch_passed": False,
            "complete_gauge_propagation_passed": False,
            "complete_first_order_reduction_passed": False,
            "uniform_open_domain_hyperbolicity_gate_passed": False,
            "full_fgc1_hyp1_health_gate_passed": False,
            "evolution_authorized": False,
            "publication_local_defocusing_gate_passed": False,
        },
        "nonclaims": serialized_certificate["nonclaims"],
    }


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
