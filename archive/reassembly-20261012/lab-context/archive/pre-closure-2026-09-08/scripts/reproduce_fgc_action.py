#!/usr/bin/env python3
"""Regenerate the FGC-1-ACT1 covariant-action/background certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from math import isfinite
from numbers import Real
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc import (  # noqa: E402
    ActionParameters,
    ModelID,
    action_certificate,
    linear_effective_mass_squared,
    required_nonclaims,
    schwarzschild_activation_radius,
    schwarzschild_gauss_bonnet,
)

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-action-gate.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-action-gate.json"
ACTION_SOURCE = REPOSITORY / "src" / "recursive_horizons" / "fgc" / "action.py"


@dataclass(frozen=True, slots=True)
class TournamentConfig:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    ricci_scalar: float
    schwarzschild_radius: float
    sample_radius_over_activation: float
    models: tuple[ActionParameters, ...]
    source_sha256: str


def _mapping(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a table")
    return value


def _exact_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    actual = set(value)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(f"{name} keys differ: missing={missing}, extra={extra}")


def _finite(name: str, value: Any, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result) or (positive and result <= 0.0):
        qualifier = "positive finite" if positive else "finite"
        raise ValueError(f"{name} must be a {qualifier} real number")
    return result


def load_config(path: Path = DEFAULT_CONFIG) -> TournamentConfig:
    """Load and strictly validate the frozen TOML tournament fixture."""

    source = path.resolve().read_bytes()
    raw = _mapping("root", tomllib.loads(source.decode("utf-8")))
    _exact_keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "background",
            "models",
        },
    )
    if raw["schema_version"] != 1:
        raise ValueError("schema_version must equal 1")
    if raw["artifact_id"] != "FGC-1-ACT1":
        raise ValueError("artifact_id must equal FGC-1-ACT1")
    if not isinstance(raw["project_version"], str) or not raw["project_version"]:
        raise ValueError("project_version must be a non-empty string")
    if raw["metric_signature"] != "-+++":
        raise ValueError("metric_signature must equal -+++")

    background = _mapping("background", raw["background"])
    _exact_keys(
        "background",
        background,
        {"ricci_scalar", "schwarzschild_radius", "sample_radius_over_activation"},
    )
    ricci = _finite("background.ricci_scalar", background["ricci_scalar"])
    radius_s = _finite(
        "background.schwarzschild_radius",
        background["schwarzschild_radius"],
        positive=True,
    )
    sample_ratio = _finite(
        "background.sample_radius_over_activation",
        background["sample_radius_over_activation"],
        positive=True,
    )

    model_tables = _mapping("models", raw["models"])
    expected_models = {model.value for model in ModelID}
    _exact_keys("models", model_tables, expected_models)
    parameter_keys = {
        "planck_mass",
        "scalar_mass",
        "quartic_coupling",
        "pulse_width",
        "scalar_field",
        "ricci_coupling",
        "linear_gb_coupling",
        "quadratic_gb_coupling",
    }
    models: list[ActionParameters] = []
    for model_id in ModelID:
        table = _mapping(f"models.{model_id.value}", model_tables[model_id.value])
        _exact_keys(f"models.{model_id.value}", table, parameter_keys)
        models.append(ActionParameters(model_id=model_id, **table))

    return TournamentConfig(
        schema_version=1,
        artifact_id="FGC-1-ACT1",
        project_version=raw["project_version"],
        metric_signature="-+++",
        ricci_scalar=ricci,
        schwarzschild_radius=radius_s,
        sample_radius_over_activation=sample_ratio,
        models=tuple(models),
        source_sha256=sha256(source).hexdigest(),
    )


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, object]:
    """Build the combined tournament record from one frozen configuration."""

    config_path = config_path.resolve()
    config = load_config(config_path)
    by_id = {parameters.model_id: parameters for parameters in config.models}
    fgc = by_id[ModelID.FGC_QR]
    activation_radius = schwarzschild_activation_radius(
        fgc, config.schwarzschild_radius
    )
    sample_radius = config.sample_radius_over_activation * activation_radius
    sample_invariant = schwarzschild_gauss_bonnet(
        config.schwarzschild_radius, sample_radius
    )
    certificates = [
        action_certificate(
            parameters,
            ricci_scalar=config.ricci_scalar,
            schwarzschild_radius=config.schwarzschild_radius,
            areal_radius=sample_radius,
        )
        for parameters in config.models
    ]

    activation_invariant = schwarzschild_gauss_bonnet(
        config.schwarzschild_radius, activation_radius
    )
    activation_mass = linear_effective_mass_squared(
        fgc,
        ricci_scalar=0.0,
        gauss_bonnet=activation_invariant,
        phi=0.0,
    )
    relative_residual = abs(activation_mass) / (fgc.scalar_mass**2)
    certificate_by_id = {entry["model_id"]: entry for entry in certificates}
    gr = certificate_by_id["GR-0"]
    linear = certificate_by_id["SGB-L"]
    qr = certificate_by_id["FGC-QR"]
    nonclaims = required_nonclaims()

    try:
        source_config = str(config_path.relative_to(REPOSITORY))
    except ValueError:
        source_config = str(config_path)
    return {
        "schema_version": config.schema_version,
        "project_version": config.project_version,
        "artifact_id": config.artifact_id,
        "artifact": "fgc1_covariant_action_and_linear_activation_certificate",
        "classification": "action_and_background_algebra_not_principal_symbol_or_collapse_result",
        "development_status": "pre_publication_gate",
        "source_config": source_config,
        "source_config_sha256": config.source_sha256,
        "implementation_sources_sha256": {
            "src/recursive_horizons/fgc/action.py": sha256(
                ACTION_SOURCE.read_bytes()
            ).hexdigest(),
            "scripts/reproduce_fgc_action.py": sha256(
                Path(__file__).resolve().read_bytes()
            ).hexdigest(),
        },
        "matter_regulator_assignment": {
            "chi": "independent_canonical_collapsing_matter_control",
            "phi": "curvature_sensitive_regulator_candidate",
            "physical_metric": "g_ab_matter_coupled_metric",
        },
        "fixtures": certificates,
        "schwarzschild_control": {
            "ricci_scalar": config.ricci_scalar,
            "schwarzschild_radius": config.schwarzschild_radius,
            "activation_radius": activation_radius,
            "sample_radius_over_activation": config.sample_radius_over_activation,
            "sample_areal_radius": sample_radius,
            "sample_gauss_bonnet": sample_invariant,
            "activation_gauss_bonnet": activation_invariant,
            "activation_mass_squared_relative_residual": relative_residual,
        },
        "verified_algebraic_checks": {
            "frozen_branch_order": [entry["model_id"] for entry in certificates]
            == [model.value for model in ModelID],
            "gr0_exact_action_at_zero_regulator": gr["low_gradient_reference_gate"][
                "exact_gr0_action"
            ],
            "sgb_l_is_curvature_sourced_at_zero_regulator": linear["background"][
                "scalar_equation_algebraic_residual"
            ]
            != 0.0,
            "fgc_qr_zero_regulator_solves_declared_scalar_equation": qr[
                "low_gradient_reference_gate"
            ]["zero_regulator_solves_declared_scalar_equation"],
            "fgc_qr_inside_sample_is_linearly_tachyonic": qr["background"][
                "linear_tachyonic_at_reference_point"
            ],
            "schwarzschild_activation_zero_within_relative_tolerance_1e_12": relative_residual
            <= 1.0e-12,
            "all_required_nonclaims_false": all(value is False for value in nonclaims.values()),
        },
        "gate_status": {
            "declared_action_and_scalar_algebra_certificate_passed": True,
            "full_fgc1_health_gate_passed": False,
            "publication_local_defocusing_gate_passed": False,
        },
        "nonclaims": nonclaims,
        "scope": "Declares and checks the exact GR-0, SGB-L, and FGC-QR action fixtures, scalar equation, dimensions, and Schwarzschild linear activation scale. It does not derive the expanded metric source, spherical principal symbol, collapse, affine defocusing, singularity resolution, child spacetime, dark sector, variable c, particle spectrum, or a theory of everything.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        record(arguments.config), indent=2, sort_keys=True, allow_nan=False
    ) + "\n"
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
