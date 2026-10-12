#!/usr/bin/env python3
"""Regenerate the fail-closed FGC-1-EFT1-OPEN1 audit record."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-eft1-open1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-eft1-open1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-eft1-open1.md"

from recursive_horizons.fgc.eft_run_envelope import (  # noqa: E402
    EFT1_OPEN1_ARTIFACT_ID,
    retained_eft_open_run_audit,
)


PROJECT_VERSION = "0.11.0"
PREDECESSORS = {
    "eft0_config": "configs/fgc/fgc-1-eft0-led1.toml",
    "eft0_result": "results/fgc-1-eft0-led1.json",
    "uhyp1_config": "configs/fgc/fgc-1-hyp1-dom3-uhyp1.toml",
    "uhyp1_result": "results/fgc-1-hyp1-dom3-uhyp1.json",
    "bnd1_config": "configs/fgc/fgc-1-hyp1-bnd1-md1.toml",
    "bnd1_result": "results/fgc-1-hyp1-bnd1-md1.json",
    "con3_config": "configs/fgc/fgc-1-hyp1-con3-cau1.toml",
    "con3_result": "results/fgc-1-hyp1-con3-cau1.json",
}
CLAIMS = {
    "complete_Wilsonian_EFT_defined": False,
    "UV_completion_derived": False,
    "radiative_stability_proven": False,
    "retained_EFT_open_run_envelope_proven": False,
    "constraint_complete_IBVP_proven": False,
    "initial_data_constructed": False,
    "evolution_authorized": False,
    "collapse_solution_derived": False,
    "affine_null_defocusing_derived": False,
    "finite_invariant_transition_surface_derived": False,
    "singularity_resolution_derived": False,
}
AUTHORIZATION_CONTRACT = {
    "all_required_predicates_must_pass": True,
    "one_false_predicate_stops_the_run": True,
    "local_component_bounds_are_not_frequency_bounds": True,
    "radial_hyperbolicity_is_not_multidirectional_hyperbolicity": True,
    "frozen_main_boundary_dissipation_is_not_constraint_preserving_IBVP": True,
    "conditional_Cauchy_uniqueness_is_not_solution_existence": True,
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/eft_run_envelope.py",
        "scripts/reproduce_fgc_eft1_open1.py",
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


def _pairs(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
        answer[key] = value
    return answer


def _canonical(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _load_unique_json(path: Path, name: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _binds(result: Mapping[str, Any], config: Path, name: str) -> None:
    ledger = result.get("source_config_sha256")
    if isinstance(ledger, Mapping):
        value = ledger.get(_rel(config))
    else:
        source_configs = result.get("source_configs")
        value = None
        if isinstance(source_configs, Mapping):
            nested = source_configs.get(_rel(config))
            if isinstance(nested, Mapping):
                value = nested.get("sha256")
    if value != _sha(config):
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
            *PREDECESSORS,
            "authorization_contract",
            "claims",
        },
    )
    if {
        "schema_version": raw["schema_version"],
        "artifact_id": raw["artifact_id"],
        "project_version": raw["project_version"],
    } != {
        "schema_version": 1,
        "artifact_id": EFT1_OPEN1_ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
    }:
        raise ValueError("EFT1-OPEN1 identity or version differs")
    paths = {name: _path(name, raw[name], expected) for name, expected in PREDECESSORS.items()}
    if raw["authorization_contract"] != AUTHORIZATION_CONTRACT:
        raise ValueError("EFT1-OPEN1 authorization contract differs")
    if raw["claims"] != CLAIMS:
        raise ValueError("EFT1-OPEN1 claims must remain fail-closed")
    return {
        "raw": raw,
        "paths": paths,
        "source_sha256": sha256(source).hexdigest(),
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    records = {
        name: _load_unique_json(config["paths"][f"{name}_result"], name)
        for name in ("eft0", "uhyp1", "bnd1", "con3")
    }
    for name in records:
        _binds(records[name], config["paths"][f"{name}_config"], name)
    audit = retained_eft_open_run_audit(**records)
    if audit["nonclaims"] != CLAIMS:
        raise ValueError("EFT1-OPEN1 core and artifact nonclaim contracts differ")
    if audit["retained_EFT_open_run_envelope_passed"] is not False:
        raise ValueError("EFT1-OPEN1 unexpectedly authorized a run")
    return {
        "schema_version": 1,
        "artifact_id": EFT1_OPEN1_ARTIFACT_ID,
        "classification": audit["classification"],
        "project_version": PROJECT_VERSION,
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): config["source_sha256"]},
        "predecessor_sha256": {
            _rel(path): _sha(path) for path in config["paths"].values()
        },
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "frozen_configuration": config["raw"],
        "open_run_authorization_audit": audit,
        "gate_status": {
            "audit_definition_and_fail_closed_composition_passed": True,
            "retained_EFT_open_run_envelope_passed": False,
            "initial_data_authorized": False,
            "evolution_authorized": False,
            "collapse_authorized": False,
            "DEF1_authorized": False,
        },
        "nonclaims": CLAIMS,
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    try:
        source = path.read_text(encoding="utf-8")
        value = json.loads(source, object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("EFT1-OPEN1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("EFT1-OPEN1 result must use canonical sorted JSON")
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
