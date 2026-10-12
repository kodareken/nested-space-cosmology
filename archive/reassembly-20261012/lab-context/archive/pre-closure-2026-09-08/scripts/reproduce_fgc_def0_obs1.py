#!/usr/bin/env python3
"""Regenerate the exact FGC-1-DEF0-OBS1 geometry controls."""

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
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-def0-obs1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-def0-obs1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-def0-obs1.md"

from scripts.reproduce_fgc_hyp1_reduction import load_config as load_red1_config  # noqa: E402
from recursive_horizons.fgc.metric_null_observable import (  # noqa: E402
    SphericalADMGeometryJet,
    metric_null_raychaudhuri_point_certificate,
)
from recursive_horizons.fgc.spherical_reduction import Jet2  # noqa: E402


Q = Fraction
ARTIFACT_ID = "FGC-1-DEF0-OBS1"
PROJECT_VERSION = "0.11.0"
RED1_CONFIG = "configs/fgc/fgc-1-hyp1-reduction.toml"
RED1_RESULT = "results/fgc-1-hyp1-reduction.json"
NONCLAIM_NAMES = (
    "HYP1_complete",
    "constraint_preserving_ACT1_IBVP_proven",
    "COL1_evolution_solution_consumed",
    "finite_trapped_interval_derived",
    "resolution_independent_positive_Raychaudhuri_margin_derived",
    "metric_null_affine_defocusing_DEF1_derived",
    "ROB1_open_parameter_neighborhood_derived",
    "non_spherical_shear_robustness_derived",
    "finite_invariant_transition_surface_derived",
    "trapped_to_anti_trapped_transition_derived",
    "singularity_resolution_derived",
    "physical_wall_derived",
)
NONCLAIMS = {name: False for name in NONCLAIM_NAMES}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/metric_null_observable.py",
        "scripts/reproduce_fgc_def0_obs1.py",
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
    answer = (REPOSITORY / value).resolve()
    if answer.relative_to(REPOSITORY).as_posix() != value or not answer.is_file():
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return answer


def _text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _fraction(name: str, value: Any) -> Fraction:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a canonical rational string")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a canonical rational string") from exc
    if value != _text(answer):
        raise ValueError(f"{name} must be a canonical rational string")
    return answer


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _text(value)
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported exact result value {type(value).__name__}")


def _pairs(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
        answer[key] = value
    return answer


def _canonical(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _load_red1_result(path: Path) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_pairs)
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("RED1 predecessor must use canonical unique-key JSON")
    return value


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
            "physical_metric_reduction_config",
            "physical_metric_reduction_result",
            "controls",
            "proof_contract",
            "open_gates",
        },
    )
    if {key: raw[key] for key in ("schema_version", "artifact_id", "project_version", "metric_signature", "riemann_convention")} != {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }:
        raise ValueError("DEF0-OBS1 identity, version, or conventions differ")
    paths = {
        "physical_metric_reduction_config": _path("physical_metric_reduction_config", raw["physical_metric_reduction_config"], RED1_CONFIG),
        "physical_metric_reduction_result": _path("physical_metric_reduction_result", raw["physical_metric_reduction_result"], RED1_RESULT),
    }
    load_red1_config(paths["physical_metric_reduction_config"])
    controls = raw["controls"]
    _keys("controls", controls, {"minkowski_spherical", "contracting_flat_de_sitter", "synthetic_negative_Rkk"})
    _keys("controls.minkowski_spherical", controls["minkowski_spherical"], {"areal_radius", "branch"})
    _keys("controls.contracting_flat_de_sitter", controls["contracting_flat_de_sitter"], {"Hubble_rate", "comoving_radius", "branch"})
    _keys("controls.synthetic_negative_Rkk", controls["synthetic_negative_Rkk"], {"areal_radius", "areal_radius_dt", "areal_radius_dr", "areal_radius_dtt", "branch"})
    if controls["minkowski_spherical"]["branch"] != "outgoing" or _fraction("Minkowski radius", controls["minkowski_spherical"]["areal_radius"]) != 4:
        raise ValueError("Minkowski control differs")
    if controls["contracting_flat_de_sitter"]["branch"] != "outgoing" or _fraction("de Sitter H", controls["contracting_flat_de_sitter"]["Hubble_rate"]) != -1 or _fraction("de Sitter radius", controls["contracting_flat_de_sitter"]["comoving_radius"]) != 2:
        raise ValueError("contracting de Sitter control differs")
    synthetic = controls["synthetic_negative_Rkk"]
    if synthetic["branch"] != "outgoing" or tuple(_fraction(f"synthetic.{name}", synthetic[name]) for name in ("areal_radius", "areal_radius_dt", "areal_radius_dr", "areal_radius_dtt")) != (Q(4), Q(-1), Q(1), Q(1)):
        raise ValueError("synthetic negative-Rkk control differs")
    if any(value is not True for value in raw["proof_contract"].values()):
        raise ValueError("every DEF0-OBS1 proof-contract gate must remain true")
    if raw["open_gates"] != NONCLAIMS:
        raise ValueError("DEF0-OBS1 open gate contract differs")
    return {"raw": raw, "paths": paths, "source_sha256": sha256(source).hexdigest()}


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    predecessor = _load_red1_result(config["paths"]["physical_metric_reduction_result"])
    if predecessor.get("artifact_id") != "FGC-1-HYP1-RED1" or predecessor.get("project_version") != PROJECT_VERSION:
        raise ValueError("RED1 predecessor identity differs")
    if predecessor.get("gate_status", {}).get("spherical_covariant_equations_evaluated_at_local_jets") is not True:
        raise ValueError("RED1 physical metric/curvature gate is not passed")
    source_ledger = predecessor.get(
        "source_config_sha256", predecessor.get("source_configs_sha256")
    )
    config_relative = _rel(config["paths"]["physical_metric_reduction_config"])
    if not isinstance(source_ledger, Mapping) or source_ledger.get(config_relative) != _sha(config["paths"]["physical_metric_reduction_config"]):
        raise ValueError("RED1 predecessor is stale against its direct config")

    minkowski = metric_null_raychaudhuri_point_certificate(
        SphericalADMGeometryJet(Jet2.constant(1), Jet2.constant(0), Jet2.constant(1), Jet2(4, dr=1)),
        branch="outgoing",
        affine_scale=Jet2.constant(1),
    )
    contracting = metric_null_raychaudhuri_point_certificate(
        SphericalADMGeometryJet(Jet2.constant(1), Jet2.constant(0), Jet2(1, dt=-1, dtt=1), Jet2(2, dt=-2, dr=1, dtt=2, dtr=-1)),
        branch="outgoing",
        affine_scale=Jet2(1, dt=1, dtt=1),
    )
    synthetic = metric_null_raychaudhuri_point_certificate(
        SphericalADMGeometryJet(Jet2.constant(1), Jet2.constant(0), Jet2.constant(1), Jet2(4, dt=-1, dr=1, dtt=1)),
        branch="outgoing",
        affine_scale=Jet2.constant(1),
    )
    controls = {
        "minkowski_spherical": minkowski,
        "contracting_flat_de_sitter": contracting,
        "synthetic_negative_Rkk": synthetic,
    }
    if any(control["nonclaims"] != NONCLAIMS for control in controls.values()):
        raise ValueError("DEF0-OBS1 core and artifact nonclaim contracts differ")
    if minkowski["raychaudhuri"]["complete_rhs"] != Q(-1, 8) or contracting["null_frame"]["round_sphere_classification"] != "trapped" or synthetic["raychaudhuri"]["complete_rhs"] != Q(1, 2):
        raise ValueError("DEF0-OBS1 exact control values differ")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "classification": "exact_physical_metric_affine_radial_null_Raychaudhuri_observable_controls_not_COL1_DEF1_or_transition",
        "project_version": PROJECT_VERSION,
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): config["source_sha256"]},
        "predecessor_sha256": {_rel(path): _sha(path) for path in config["paths"].values()},
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "frozen_configuration": _serial(config["raw"]),
        "controls": _serial(controls),
        "gate_status": {
            "exact_metric_null_observable_contract_and_controls_passed": True,
            "metric_null_affine_defocusing_DEF1_derived": False,
            "evolution_authorized": False,
        },
        "nonclaims": NONCLAIMS,
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    try:
        source = path.read_text(encoding="utf-8")
        value = json.loads(source, object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("DEF0-OBS1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("DEF0-OBS1 result must use canonical sorted JSON")
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
