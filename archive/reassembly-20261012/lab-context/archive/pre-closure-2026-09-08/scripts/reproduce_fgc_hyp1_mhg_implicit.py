#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP1-MHG3-IMP1 local implicit-branch certificate."""

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

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-mhg-implicit.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-hyp1-mhg-implicit.json"
OWNER_DOCUMENT = REPOSITORY / "docs" / "fgc-hyp1-mhg-implicit.md"
IMPLEMENTATION_FILES = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "exact_linear_algebra.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "exact_tangent.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_reduction.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "spherical_symbol.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "modified_harmonic.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "reference_connection.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "modified_harmonic_reference.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "modified_harmonic_implicit.py",
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "__init__.py",
    REPOSITORY / "scripts" / "reproduce_fgc_action.py",
    REPOSITORY / "scripts" / "reproduce_fgc_metric_variation.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_reduction.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_symbol.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_modified_harmonic.py",
    REPOSITORY / "scripts" / "reproduce_fgc_hyp1_mhg_reference.py",
    Path(__file__).resolve(),
)

from scripts.reproduce_fgc_action import load_config as load_action_config
from scripts.reproduce_fgc_metric_variation import load_config as load_variation_config
from scripts.reproduce_fgc_hyp1_reduction import load_config as load_reduction_config
from scripts.reproduce_fgc_hyp1_symbol import load_config as load_symbol_config
from scripts.reproduce_fgc_hyp1_modified_harmonic import (
    load_config as load_modified_harmonic_config,
)
from scripts.reproduce_fgc_hyp1_mhg_reference import (
    load_config as load_reference_config,
)
from recursive_horizons.fgc.action import ModelID
from recursive_horizons.fgc.modified_harmonic_implicit import (
    EXPECTED_FGCQR_PARAMETERS,
    IMPLICIT_ACCELERATION_ORDER,
    IMPLICIT_EQUATION_ORDER,
    modified_harmonic_implicit_certificate,
    required_imp1_nonclaims,
)
from recursive_horizons.fgc.reference_connection import (
    SPHERICAL_FLAT_REFERENCE_ID,
    flat_spherical_annulus_reference,
)


@dataclass(frozen=True, slots=True)
class HYP1MHGImplicitConfig:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    riemann_convention: str
    source_paths: dict[str, Path]
    source_relatives: dict[str, str]
    radial_domain_minimum: Fraction
    coordinate_radius: Fraction
    tilde_normal_factor: Fraction
    hat_normal_factor: Fraction
    fixture: dict[str, Any]
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


def _string_list(name: str, value: Any, expected: tuple[str, ...]) -> None:
    if not isinstance(value, list) or value != list(expected):
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


def _parse_fixture(raw: Mapping[str, Any]) -> dict[str, Any]:
    _exact_keys(
        "solution_fixture",
        raw,
        {"fixture_id", "model_id", "purpose", "action_parameters", "state"},
    )
    frozen_strings = {
        "fixture_id": "FGCQR_flat_reference_vacuum",
        "model_id": "FGC-QR",
        "purpose": "exact_flat_reference_gauge_vacuum_ift_root",
    }
    for name, expected in frozen_strings.items():
        if raw[name] != expected:
            raise ValueError(f"solution_fixture.{name} is not frozen")

    parameters = _mapping("solution_fixture.action_parameters", raw["action_parameters"])
    parameter_keys = {
        "planck_mass",
        "scalar_mass",
        "quartic_coupling",
        "ricci_coupling",
        "linear_gb_coupling",
        "quadratic_gb_coupling",
    }
    _exact_keys("solution_fixture.action_parameters", parameters, parameter_keys)
    parsed_parameters = {
        name: _fraction(f"solution_fixture.action_parameters.{name}", parameters[name])
        for name in parameter_keys
    }
    expected_parameters = {
        "planck_mass": EXPECTED_FGCQR_PARAMETERS["planck_mass"],
        "scalar_mass": EXPECTED_FGCQR_PARAMETERS["mu"],
        "quartic_coupling": EXPECTED_FGCQR_PARAMETERS["g4"],
        "ricci_coupling": EXPECTED_FGCQR_PARAMETERS["beta"],
        "linear_gb_coupling": EXPECTED_FGCQR_PARAMETERS["alpha"],
        "quadratic_gb_coupling": EXPECTED_FGCQR_PARAMETERS["eta"],
    }
    if parsed_parameters != expected_parameters:
        raise ValueError("solution_fixture action parameters are not the frozen FGC-QR values")

    state = _mapping("solution_fixture.state", raw["state"])
    state_keys = {"alpha", "shift", "lambda", "areal_radius", "phi", "chi"}
    _exact_keys("solution_fixture.state", state, state_keys)
    parsed_state: dict[str, dict[str, Fraction]] = {}
    for field in ("alpha", "shift", "lambda", "phi", "chi"):
        table = _mapping(f"solution_fixture.state.{field}", state[field])
        _exact_keys(f"solution_fixture.state.{field}", table, {"value"})
        parsed_state[field] = {
            "value": _fraction(f"solution_fixture.state.{field}.value", table["value"])
        }
    radius = _mapping("solution_fixture.state.areal_radius", state["areal_radius"])
    _exact_keys("solution_fixture.state.areal_radius", radius, {"value", "dr"})
    parsed_state["areal_radius"] = {
        "value": _fraction("solution_fixture.state.areal_radius.value", radius["value"]),
        "dr": _fraction("solution_fixture.state.areal_radius.dr", radius["dr"]),
    }
    expected_state = {
        "alpha": {"value": Fraction(1)},
        "shift": {"value": Fraction(0)},
        "lambda": {"value": Fraction(1)},
        "areal_radius": {"value": Fraction(4), "dr": Fraction(1)},
        "phi": {"value": Fraction(0)},
        "chi": {"value": Fraction(0)},
    }
    if parsed_state != expected_state:
        raise ValueError("solution_fixture.state is not the frozen flat vacuum root")
    return {
        **frozen_strings,
        "action_parameters": parsed_parameters,
        "state": parsed_state,
    }


def load_config(path: Path = DEFAULT_CONFIG) -> HYP1MHGImplicitConfig:
    """Load IMP1 and reconcile its exact root against all predecessor gates."""

    path = path.resolve()
    source = path.read_bytes()
    raw = _mapping("root", tomllib.loads(source.decode("utf-8")))
    path_keys = {
        "action_config",
        "variation_config",
        "reduction_config",
        "symbol_config",
        "modified_harmonic_config",
        "reference_config",
    }
    _exact_keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            *path_keys,
            "reference",
            "formulation",
            "solution_fixture",
            "proof_contract",
            "open_gates",
        },
    )
    expected_scalars = {
        "schema_version": 1,
        "artifact_id": "FGC-1-HYP1-MHG3-IMP1",
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
        "reference_config": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
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
        ("MHG1", load_modified_harmonic_config(source_paths["modified_harmonic_config"])),
        ("REF1", load_reference_config(source_paths["reference_config"])),
    )
    for name, predecessor in predecessors:
        if predecessor.project_version != raw["project_version"]:
            raise ValueError(f"{name} and IMP1 project versions differ")
        if predecessor.metric_signature != raw["metric_signature"]:
            raise ValueError(f"{name} and IMP1 metric signatures differ")
    for name, predecessor in predecessors[1:]:
        if predecessor.riemann_convention != raw["riemann_convention"]:
            raise ValueError(f"{name} and IMP1 Riemann conventions differ")

    reference = _mapping("reference", raw["reference"])
    _exact_keys(
        "reference",
        reference,
        {
            "reference_id",
            "radial_domain_minimum",
            "coordinate_radius",
            "center_included",
        },
    )
    if reference["reference_id"] != SPHERICAL_FLAT_REFERENCE_ID:
        raise ValueError("reference.reference_id is not frozen")
    if reference["center_included"] is not False:
        raise ValueError("IMP1 must not claim that its annulus includes the center")
    radial_minimum = _fraction(
        "reference.radial_domain_minimum",
        reference["radial_domain_minimum"],
        positive=True,
    )
    coordinate_radius = _fraction(
        "reference.coordinate_radius", reference["coordinate_radius"], positive=True
    )
    if coordinate_radius <= radial_minimum:
        raise ValueError("reference.coordinate_radius must lie inside the annulus")

    formulation = _mapping("formulation", raw["formulation"])
    _exact_keys(
        "formulation",
        formulation,
        {
            "physical_equations",
            "gauge_equations",
            "acceleration_variables",
            "acceleration_order",
            "equation_order",
            "differentiation_method",
            "implicit_theorem",
            "tilde_normal_factor",
            "hat_normal_factor",
        },
    )
    frozen_formulation = {
        "physical_equations": "unredefined_ACT1_VAR1",
        "gauge_equations": "REF1_full_reference_connection_modified_harmonic",
        "acceleration_variables": "base_metric_and_scalar_coordinate_time_second_derivatives",
        "differentiation_method": "exact_first_tangent_forward_ad",
        "implicit_theorem": "finite_dimensional_real_implicit_function_theorem",
    }
    for name, expected in frozen_formulation.items():
        if formulation[name] != expected:
            raise ValueError(f"formulation.{name} is not frozen")
    _string_list(
        "formulation.acceleration_order",
        formulation["acceleration_order"],
        IMPLICIT_ACCELERATION_ORDER,
    )
    _string_list(
        "formulation.equation_order",
        formulation["equation_order"],
        IMPLICIT_EQUATION_ORDER,
    )
    tilde_factor = _positive_integer(
        "formulation.tilde_normal_factor", formulation["tilde_normal_factor"]
    )
    hat_factor = _positive_integer(
        "formulation.hat_normal_factor", formulation["hat_normal_factor"]
    )
    if not 1 < tilde_factor < hat_factor:
        raise ValueError("formulation requires 1 < tilde factor < hat factor")

    reference_predecessor = predecessors[5][1]
    if (
        radial_minimum != reference_predecessor.radial_domain_minimum
        or tilde_factor != reference_predecessor.tilde_normal_factor
        or hat_factor != reference_predecessor.hat_normal_factor
    ):
        raise ValueError("IMP1 reference or auxiliary factors differ from REF1")

    fixture = _parse_fixture(_mapping("solution_fixture", raw["solution_fixture"]))
    if fixture["state"]["areal_radius"]["value"] != coordinate_radius:
        raise ValueError("IMP1 fixture radius differs from reference coordinate radius")

    action = predecessors[0][1]
    action_model = next(
        model for model in action.models if model.model_id is ModelID.FGC_QR
    )
    action_values = {
        "planck_mass": Fraction(str(action_model.planck_mass)),
        "scalar_mass": Fraction(str(action_model.scalar_mass)),
        "quartic_coupling": Fraction(str(action_model.quartic_coupling)),
        "ricci_coupling": Fraction(str(action_model.ricci_coupling)),
        "linear_gb_coupling": Fraction(str(action_model.linear_gb_coupling)),
        "quadratic_gb_coupling": Fraction(str(action_model.quadratic_gb_coupling)),
    }
    if fixture["action_parameters"] != action_values:
        raise ValueError("IMP1 FGC-QR parameters differ from ACT1")

    proof = _mapping("proof_contract", raw["proof_contract"])
    proof_keys = {
        "require_exact_fgcqr_action_identity",
        "require_exact_flat_reference_gauge_root",
        "require_exact_complete_ref1_residual_zero",
        "require_exact_first_tangent_primal_regression",
        "require_exact_mhg1_coordinate_time_kinetic_regression",
        "require_exact_nonzero_acceleration_jacobian_determinant",
        "require_exact_full_rank_acceleration_jacobian",
        "require_exact_two_sided_inverse_identities",
        "require_all_smooth_domain_denominators_regular_at_root",
        "require_local_implicit_function_theorem_hypotheses",
    }
    _exact_keys("proof_contract", proof, proof_keys)
    if any(value is not True for value in proof.values()):
        raise ValueError("every IMP1 proof contract entry must be true")

    open_gates = _mapping("open_gates", raw["open_gates"])
    _exact_keys("open_gates", open_gates, set(required_imp1_nonclaims()))
    if any(value is not False for value in open_gates.values()):
        raise ValueError("every IMP1 open gate must remain false")

    return HYP1MHGImplicitConfig(
        schema_version=raw["schema_version"],
        artifact_id=raw["artifact_id"],
        project_version=raw["project_version"],
        metric_signature=raw["metric_signature"],
        riemann_convention=raw["riemann_convention"],
        source_paths=source_paths,
        source_relatives=source_relatives,
        radial_domain_minimum=radial_minimum,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
        fixture=fixture,
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
    reference = flat_spherical_annulus_reference(
        radial_domain_minimum=config.radial_domain_minimum
    )
    certificate = modified_harmonic_implicit_certificate(
        {
            "fixture": config.fixture,
            "reference": reference,
            "coordinate_radius": config.coordinate_radius,
            "tilde_normal_factor": config.tilde_normal_factor,
            "hat_normal_factor": config.hat_normal_factor,
        }
    )
    if (
        not certificate["all_declared_exact_checks_pass"]
        or not certificate[
            "local_smooth_implicit_coordinate_time_acceleration_branch_proven"
        ]
    ):
        raise ValueError("IMP1 exact certificate failed a required proof check")
    serialized = _serialize_exact(certificate)
    source_paths = tuple(config.source_paths.values()) + (config_path.resolve(),)
    return {
        "schema_version": config.schema_version,
        "artifact_id": config.artifact_id,
        "artifact_label": "FGC-1 HYP1 exact local implicit coordinate-time acceleration branch gate",
        "project_version": config.project_version,
        "classification": certificate["classification"],
        "development_status": "flat_fgcqr_local_implicit_acceleration_branch_complete_propagation_domain_and_evolution_gates_open",
        "generated_by": _relative(Path(__file__)),
        "source_configs": {
            "action": config.source_relatives["action_config"],
            "variation": config.source_relatives["variation_config"],
            "reduction": config.source_relatives["reduction_config"],
            "symbol": config.source_relatives["symbol_config"],
            "modified_harmonic": config.source_relatives[
                "modified_harmonic_config"
            ],
            "reference": config.source_relatives["reference_config"],
            "implicit": _relative(config_path),
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
        "input_fixture_id": config.fixture["fixture_id"],
        "input_model_id": config.fixture["model_id"],
        "implicit_acceleration_certificate": serialized,
        "implicit_acceleration_certificate_sha256": _object_sha256(serialized),
        "gate_status": {
            "flat_fgcqr_local_implicit_acceleration_branch_passed": True,
            "quantified_implicit_branch_neighborhood_passed": False,
            "complete_gauge_propagation_passed": False,
            "complete_first_order_reduction_passed": False,
            "uniform_open_domain_hyperbolicity_gate_passed": False,
            "full_fgc1_hyp1_health_gate_passed": False,
            "evolution_authorized": False,
            "publication_local_defocusing_gate_passed": False,
        },
        "nonclaims": serialized["nonclaims"],
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
