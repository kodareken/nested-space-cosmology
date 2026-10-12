#!/usr/bin/env python3
"""Regenerate the conditional FGC-1-HYP1-CON3-CAU1 record."""

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
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hyp1-con3-cau1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hyp1-con3-cau1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hyp1-con3-cau1.md"

from scripts.reproduce_fgc_hyp1_con1_comp1 import (  # noqa: E402
    load_canonical_result as load_comp1_result,
    load_config as load_comp1_config,
)
from scripts.reproduce_fgc_hyp1_con2_mprop1 import (  # noqa: E402
    load_canonical_result as load_con2_result,
    load_config as load_con2_config,
)
from recursive_horizons.fgc.modified_harmonic_cauchy import (  # noqa: E402
    conditional_gauge_cauchy_certificate,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-HYP1-CON3-CAU1"
PROJECT_VERSION = "0.11.0"
PREDECESSORS = {
    "metric_subsidiary_config": "configs/fgc/fgc-1-hyp1-con2-mprop1.toml",
    "metric_subsidiary_result": "results/fgc-1-hyp1-con2-mprop1.json",
    "compatible_data_config": "configs/fgc/fgc-1-hyp1-con1-comp1.toml",
    "compatible_data_result": "results/fgc-1-hyp1-con1-comp1.json",
}
NONCLAIM_NAMES = (
    "smooth_nonlinear_REF1_solution_constructed",
    "nonzero_width_compatible_initial_hypersurface_constructed",
    "physical_initial_constraint_manifold_solved",
    "uniform_normal_gauge_map_on_UHYP1_graph_proven",
    "metric_derived_incoming_gauge_boundary_map_to_main_fields_derived",
    "physical_constraint_propagation_proven",
    "constraint_preserving_ACT1_IBVP_proven",
    "quasilinear_local_existence_proven",
    "evolution_authorized",
    "collapse_solution_derived",
    "finite_invariant_transition_surface_derived",
    "singularity_resolution_derived",
)
NONCLAIMS = {name: False for name in NONCLAIM_NAMES}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/modified_harmonic_cauchy.py",
        "src/recursive_horizons/fgc/modified_harmonic_metric_propagation.py",
        "src/recursive_horizons/fgc/modified_harmonic_constraints.py",
        "scripts/reproduce_fgc_hyp1_con3_cau1.py",
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
    try:
        relative = answer.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside repository") from exc
    if relative != value or not answer.is_file():
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
    raise TypeError(f"unsupported exact record value {type(value).__name__}")


def _pairs(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
        answer[key] = value
    return answer


def _canonical(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _binds(result: Mapping[str, Any], config: Path, name: str) -> None:
    ledger = result.get("source_config_sha256")
    if not isinstance(ledger, Mapping) or ledger.get(_rel(config)) != _sha(config):
        raise ValueError(f"{name} result is stale against its direct source config")


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
            *PREDECESSORS,
            "reference",
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
        raise ValueError("CON3-CAU1 identity, version, or conventions differ")
    paths = {name: _path(name, raw[name], expected) for name, expected in PREDECESSORS.items()}
    comp1 = load_comp1_config(paths["compatible_data_config"])
    load_con2_config(paths["metric_subsidiary_config"])
    _keys(
        "reference",
        raw["reference"],
        {"reference_id", "radial_domain_minimum", "coordinate_radius", "tilde_normal_factor", "hat_normal_factor"},
    )
    reference = raw["reference"]
    if (
        reference["reference_id"] != "flat_spherical_annulus"
        or _fraction("reference.radial_domain_minimum", reference["radial_domain_minimum"]) != Q(1, 2)
        or _fraction("reference.coordinate_radius", reference["coordinate_radius"]) != Q(4)
        or _fraction("reference.tilde_normal_factor", reference["tilde_normal_factor"]) != Q(4)
        or _fraction("reference.hat_normal_factor", reference["hat_normal_factor"]) != Q(9)
        or comp1.coordinate_radius != Q(4)
    ):
        raise ValueError("CON3-CAU1 reference contract differs")
    if any(value is not True for value in raw["proof_contract"].values()):
        raise ValueError("every CON3-CAU1 proof-contract gate must remain true")
    if raw["open_gates"] != NONCLAIMS:
        raise ValueError("CON3-CAU1 open gate contract differs")
    return {"raw": raw, "paths": paths, "comp1": comp1, "source_sha256": sha256(source).hexdigest()}


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    con2 = load_con2_result(config["paths"]["metric_subsidiary_result"])
    comp1 = load_comp1_result(config["paths"]["compatible_data_result"])
    if con2.get("artifact_id") != "FGC-1-HYP1-CON2-MPROP1" or con2.get("project_version") != PROJECT_VERSION:
        raise ValueError("CON2-MPROP1 result identity differs")
    if comp1.get("artifact_id") != "FGC-1-HYP1-CON1-COMP1" or comp1.get("project_version") != PROJECT_VERSION:
        raise ValueError("CON1-COMP1 result identity differs")
    if con2.get("gate_status", {}).get("local_differential_identity") is not True:
        raise ValueError("CON2-MPROP1 local differential identity is not passed")
    if comp1.get("gate_status", {}).get("exact_activated_local_compatible_constraint_datum_passed") is not True:
        raise ValueError("CON1-COMP1 compatible-data gate is not passed")
    if any(value is not False for value in con2.get("nonclaims", {}).values()) or any(value is not False for value in comp1.get("nonclaims", {}).values()):
        raise ValueError("a predecessor nonclaim was promoted")
    _binds(con2, config["paths"]["metric_subsidiary_config"], "CON2-MPROP1")
    _binds(comp1, config["paths"]["compatible_data_config"], "CON1-COMP1")
    state = activated_compatible_state(config["comp1"].qift1_flat_fixture)["state"]
    metric = con2["controls"]["activated_comp1_deterministic_third_jet"]["metric_derived_gauge_subsidiary"]
    reference_config = config["raw"]["reference"]
    certificate = conditional_gauge_cauchy_certificate(
        state,
        reference=flat_spherical_annulus_reference(radial_domain_minimum=Q(1, 2)),
        coordinate_radius=_fraction("coordinate_radius", reference_config["coordinate_radius"]),
        tilde_normal_factor=_fraction("tilde_normal_factor", reference_config["tilde_normal_factor"]),
        hat_normal_factor=_fraction("hat_normal_factor", reference_config["hat_normal_factor"]),
        metric_subsidiary_record=metric,
        comp1_constraint_record=comp1,
    )
    if certificate["nonclaims"] != NONCLAIMS:
        raise ValueError("CON3-CAU1 core and artifact nonclaim contracts differ")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "classification": "exact_COMP1_compatible_point_witness_and_conditional_boundary_free_gauge_Cauchy_uniqueness_not_solution_hypersurface_or_IBVP",
        "project_version": PROJECT_VERSION,
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): config["source_sha256"]},
        "predecessor_sha256": {_rel(path): _sha(path) for path in config["paths"].values()},
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "frozen_configuration": _serial(config["raw"]),
        "conditional_gauge_Cauchy_certificate": _serial(certificate),
        "gate_status": {
            "conditional_boundary_free_gauge_Cauchy_uniqueness_theorem_derived": True,
            "exact_COMP1_local_compatible_Cauchy_data_witness_passed": True,
            "constraint_preserving_ACT1_IBVP_proven": False,
            "evolution_authorized": False,
        },
        "nonclaims": NONCLAIMS,
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    try:
        source = path.read_text(encoding="utf-8")
        value = json.loads(source, object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("CON3-CAU1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("CON3-CAU1 result must use canonical sorted JSON")
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
