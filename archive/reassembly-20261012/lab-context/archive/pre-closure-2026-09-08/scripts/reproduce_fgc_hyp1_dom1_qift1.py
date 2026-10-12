#!/usr/bin/env python3
"""Regenerate the FGC-1-HYP1-DOM1-QIFT1 interval branch certificate."""

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

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-dom1-qift1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-hyp1-dom1-qift1.json"
OWNER_DOCUMENT = REPOSITORY / "docs" / "fgc-hyp1-dom1-qift1.md"

from scripts.reproduce_fgc_hyp1_fo1_rc1 import (
    IMPLEMENTATION_FILES as FO1_IMPLEMENTATION_FILES,
    load_config as load_fo1_config,
)
from recursive_horizons.fgc.exact_interval import Interval
from recursive_horizons.fgc.modified_harmonic_quantified_domain import (
    QIFT1_ACCELERATION_ORDER,
    QIFT1_PARAMETER_ORDER,
    quantified_implicit_branch_certificate,
    required_qift1_nonclaims,
)
from recursive_horizons.fgc.reference_connection import (
    SPHERICAL_FLAT_REFERENCE_ID,
    flat_spherical_annulus_reference,
)


Q = Fraction
IMPLEMENTATION_FILES = tuple(
    dict.fromkeys(
        FO1_IMPLEMENTATION_FILES
        + (
            REPOSITORY / "src" / "recursive_horizons" / "fgc" / "exact_interval.py",
            REPOSITORY / "src" / "recursive_horizons" / "fgc" / "interval_tangent.py",
            REPOSITORY
            / "src"
            / "recursive_horizons"
            / "fgc"
            / "modified_harmonic_quantified_domain.py",
            Path(__file__).resolve(),
        )
    )
)


@dataclass(frozen=True, slots=True)
class HYP1DOM1QIFT1Config:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    riemann_convention: str
    first_order_path: Path
    first_order_relative: str
    predecessor_source_paths: dict[str, Path]
    predecessor_source_relatives: dict[str, str]
    radial_domain_minimum: Fraction
    coordinate_radius: Fraction
    tilde_normal_factor: Fraction
    hat_normal_factor: Fraction
    parameter_half_width: Fraction
    acceleration_half_width: Fraction
    flat_fixture: dict[str, Any]
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
        result = Q(value)
    elif isinstance(value, str) and value:
        try:
            result = Q(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{name} must be a canonical rational string") from exc
        canonical = _fraction_text(result)
        if value != canonical:
            raise ValueError(f"{name} must use canonical rational encoding")
    else:
        raise ValueError(f"{name} must be an exact rational value")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


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
    if canonical != relative or not resolved.is_file():
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return resolved, canonical.as_posix()


def _string_list(name: str, value: Any, expected: tuple[str, ...]) -> None:
    if not isinstance(value, list) or value != list(expected):
        raise ValueError(f"{name} is not the frozen ordered string list")


def load_config(path: Path = DEFAULT_CONFIG) -> HYP1DOM1QIFT1Config:
    """Load QIFT1 and reconcile it through the complete FO1 predecessor chain."""

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
            "first_order_config",
            "reference",
            "formulation",
            "quantified_parameter_box",
            "acceleration_box",
            "proof_contract",
            "open_gates",
        },
    )
    expected_scalars = {
        "schema_version": 1,
        "artifact_id": "FGC-1-HYP1-DOM1-QIFT1",
        "project_version": "0.11.0",
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }
    for name, expected in expected_scalars.items():
        if raw[name] != expected:
            raise ValueError(f"{name} must equal {expected}")

    first_order_path, first_order_relative = _repository_path(
        "first_order_config", raw["first_order_config"]
    )
    if first_order_relative != "configs/fgc/fgc-1-hyp1-fo1-rc1.toml":
        raise ValueError("first_order_config must name the frozen FO1-RC1 predecessor")
    predecessor = load_fo1_config(first_order_path)
    if (
        predecessor.project_version != raw["project_version"]
        or predecessor.metric_signature != raw["metric_signature"]
        or predecessor.riemann_convention != raw["riemann_convention"]
    ):
        raise ValueError("QIFT1 and FO1-RC1 conventions differ")

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
        raise ValueError("QIFT1 reference ID is not frozen")
    if reference["center_included"] is not False:
        raise ValueError("QIFT1 annulus must exclude the regular center")
    radial_minimum = _fraction(
        "reference.radial_domain_minimum",
        reference["radial_domain_minimum"],
        positive=True,
    )
    coordinate_radius = _fraction(
        "reference.coordinate_radius", reference["coordinate_radius"], positive=True
    )
    if radial_minimum != predecessor.radial_domain_minimum:
        raise ValueError("QIFT1 and FO1-RC1 reference annuli differ")
    if coordinate_radius != predecessor.coordinate_radii[
        "FGCQR_flat_reference_vacuum"
    ]:
        raise ValueError("QIFT1 coordinate root differs from FO1-RC1")

    formulation = _mapping("formulation", raw["formulation"])
    _exact_keys(
        "formulation",
        formulation,
        {
            "residual",
            "parameter_definition",
            "acceleration_definition",
            "parameter_order",
            "acceleration_order",
            "norm",
            "theorem_route",
            "interval_arithmetic",
            "differentiation_method",
            "tilde_normal_factor",
            "hat_normal_factor",
        },
    )
    expected_formulation = {
        "residual": "complete_REF1_residual_R(a;z)=0",
        "parameter_definition": "z=(u,p,q,p_r,q_r)_in_frozen_base_field_order",
        "acceleration_definition": "a=p_t_in_frozen_base_field_order",
        "norm": "unscaled_componentwise_infinity_norm_in_frozen_field_order",
        "theorem_route": "uniform_banach_contraction_with_krawczyk_box_inclusion",
        "interval_arithmetic": "closed_exact_rational_endpoints_no_floating_point",
        "differentiation_method": "exact_rational_interval_first_tangent_forward_ad",
    }
    for name, expected in expected_formulation.items():
        if formulation[name] != expected:
            raise ValueError(f"formulation.{name} is not frozen")
    _string_list("formulation.parameter_order", formulation["parameter_order"], QIFT1_PARAMETER_ORDER)
    _string_list(
        "formulation.acceleration_order",
        formulation["acceleration_order"],
        QIFT1_ACCELERATION_ORDER,
    )
    tilde_factor = _fraction(
        "formulation.tilde_normal_factor",
        formulation["tilde_normal_factor"],
        positive=True,
    )
    hat_factor = _fraction(
        "formulation.hat_normal_factor",
        formulation["hat_normal_factor"],
        positive=True,
    )
    if (
        tilde_factor != predecessor.tilde_normal_factor
        or hat_factor != predecessor.hat_normal_factor
        or not 1 < tilde_factor < hat_factor
    ):
        raise ValueError("QIFT1 auxiliary factors differ from FO1-RC1")

    parameter_box = _mapping(
        "quantified_parameter_box", raw["quantified_parameter_box"]
    )
    _exact_keys(
        "quantified_parameter_box",
        parameter_box,
        {
            "dimension",
            "uniform_half_width",
            "all_axes_have_nonzero_width",
            "coordinate_radius_is_fixed_not_an_axis",
        },
    )
    if parameter_box["dimension"] != len(QIFT1_PARAMETER_ORDER):
        raise ValueError("QIFT1 parameter dimension must be thirty")
    if parameter_box["all_axes_have_nonzero_width"] is not True:
        raise ValueError("QIFT1 every parameter axis must have nonzero width")
    if parameter_box["coordinate_radius_is_fixed_not_an_axis"] is not True:
        raise ValueError("QIFT1 coordinate radius must remain a fixed parameter")
    parameter_half_width = _fraction(
        "quantified_parameter_box.uniform_half_width",
        parameter_box["uniform_half_width"],
        positive=True,
    )

    acceleration_box = _mapping("acceleration_box", raw["acceleration_box"])
    _exact_keys(
        "acceleration_box",
        acceleration_box,
        {
            "dimension",
            "uniform_half_width",
            "center_is_flat_root_zero_acceleration",
        },
    )
    if acceleration_box["dimension"] != len(QIFT1_ACCELERATION_ORDER):
        raise ValueError("QIFT1 acceleration dimension must be six")
    if acceleration_box["center_is_flat_root_zero_acceleration"] is not True:
        raise ValueError("QIFT1 acceleration box must be centered on the IMP1 root")
    acceleration_half_width = _fraction(
        "acceleration_box.uniform_half_width",
        acceleration_box["uniform_half_width"],
        positive=True,
    )

    proof = _mapping("proof_contract", raw["proof_contract"])
    proof_keys = {
        "require_exact_imp1_flat_root_and_inverse",
        "require_full_thirty_dimensional_parameter_box",
        "require_every_parameter_axis_nonzero_width",
        "require_exact_interval_ref1_residual_enclosure",
        "require_exact_interval_acceleration_jacobian_enclosure",
        "require_center_jacobian_containment",
        "require_strict_regular_chart_margins",
        "require_contraction_infinity_norm_strictly_below_one",
        "require_krawczyk_image_strictly_inside_acceleration_box",
        "require_unique_acceleration_root_for_every_parameter_point",
        "require_no_floating_point_or_sampled_proof",
    }
    _exact_keys("proof_contract", proof, proof_keys)
    if any(value is not True for value in proof.values()):
        raise ValueError("every QIFT1 proof contract entry must be true")

    open_gates = _mapping("open_gates", raw["open_gates"])
    _exact_keys("open_gates", open_gates, set(required_qift1_nonclaims()))
    if any(value is not False for value in open_gates.values()):
        raise ValueError("every QIFT1 open gate must remain false")

    return HYP1DOM1QIFT1Config(
        schema_version=raw["schema_version"],
        artifact_id=raw["artifact_id"],
        project_version=raw["project_version"],
        metric_signature=raw["metric_signature"],
        riemann_convention=raw["riemann_convention"],
        first_order_path=first_order_path,
        first_order_relative=first_order_relative,
        predecessor_source_paths=predecessor.source_paths,
        predecessor_source_relatives=predecessor.source_relatives,
        radial_domain_minimum=radial_minimum,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
        parameter_half_width=parameter_half_width,
        acceleration_half_width=acceleration_half_width,
        flat_fixture=predecessor.flat_fixture,
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
    if isinstance(value, Interval):
        return {
            "lower": _fraction_text(value.lower),
            "upper": _fraction_text(value.upper),
        }
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


def _reject_duplicate_json_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def _canonical_json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    """Load one generated result only if its JSON bytes are canonical and unique-keyed."""

    try:
        source = path.read_text(encoding="utf-8")
        payload = json.loads(source, object_pairs_hook=_reject_duplicate_json_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("QIFT1 result must be valid unique-key JSON") from exc
    if not isinstance(payload, dict):
        raise ValueError("QIFT1 result must be a JSON object")
    if source != _canonical_json_text(payload):
        raise ValueError("QIFT1 result must use canonical sorted indented JSON encoding")
    return payload


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    certificate = quantified_implicit_branch_certificate(
        {
            "flat_fixture": config.flat_fixture,
            "reference": flat_spherical_annulus_reference(
                radial_domain_minimum=config.radial_domain_minimum
            ),
            "coordinate_radius": config.coordinate_radius,
            "tilde_normal_factor": config.tilde_normal_factor,
            "hat_normal_factor": config.hat_normal_factor,
            "parameter_half_width": config.parameter_half_width,
            "acceleration_half_width": config.acceleration_half_width,
        }
    )
    if not (
        certificate["all_declared_exact_checks_pass"]
        and certificate[
            "quantified_full_dimensional_local_implicit_branch_box_derived"
        ]
    ):
        raise ValueError("QIFT1 exact interval certificate failed")
    serialized = _serialize_exact(certificate)
    predecessor_paths = tuple(config.predecessor_source_paths.values())
    source_paths = predecessor_paths + (config.first_order_path, config_path.resolve())
    source_configs = {
        "action": config.predecessor_source_relatives["action_config"],
        "variation": config.predecessor_source_relatives["variation_config"],
        "reduction": config.predecessor_source_relatives["reduction_config"],
        "symbol": config.predecessor_source_relatives["symbol_config"],
        "modified_harmonic": config.predecessor_source_relatives[
            "modified_harmonic_config"
        ],
        "reference": config.predecessor_source_relatives["reference_config"],
        "implicit": config.predecessor_source_relatives["implicit_config"],
        "propagation": config.predecessor_source_relatives["propagation_config"],
        "first_order": config.first_order_relative,
        "quantified_domain": _relative(config_path),
    }
    return {
        "schema_version": config.schema_version,
        "artifact_id": config.artifact_id,
        "artifact_label": "FGC-1 HYP1 full-dimensional exact rational-interval quantified implicit-branch gate",
        "project_version": config.project_version,
        "classification": certificate["classification"],
        "development_status": "full_dimensional_local_implicit_branch_box_complete_constraint_hyperbolicity_eft_and_evolution_gates_open",
        "generated_by": _relative(Path(__file__)),
        "source_configs": source_configs,
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
        "quantified_parameter_box": _serialize_exact(
            config.configuration["quantified_parameter_box"]
        ),
        "acceleration_box": _serialize_exact(
            config.configuration["acceleration_box"]
        ),
        "proof_contract": _serialize_exact(config.configuration["proof_contract"]),
        "input_fixture_ids": [config.flat_fixture["fixture_id"]],
        "input_model_ids": [config.flat_fixture["model_id"]],
        "quantified_branch_certificate": serialized,
        "quantified_branch_certificate_sha256": _object_sha256(serialized),
        "gate_status": {
            "quantified_full_dimensional_local_implicit_branch_box_passed": True,
            "complete_metric_gauge_physical_constraint_propagation_passed": False,
            "uniform_radial_strong_hyperbolicity_box_passed": False,
            "retained_eft_validity_envelope_passed": False,
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
    arguments.output.write_text(_canonical_json_text(payload), encoding="utf-8")
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
