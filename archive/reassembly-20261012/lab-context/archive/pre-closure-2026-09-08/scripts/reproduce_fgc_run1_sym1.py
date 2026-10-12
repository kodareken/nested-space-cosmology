#!/usr/bin/env python3
"""Regenerate the fail-closed FGC-1-RUN1-SYM1 authorization record."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-run1-sym1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-run1-sym1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-run1-sym1.md"

from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_EVIDENCE_REQUIREMENTS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    SF1_PROTOCOL_V2_ARTIFACT_ID,
    scoped_classical_spherical_run_audit,
    validate_sf1_protocol,
)


PROJECT_VERSION = "0.11.0"
PROTOCOL_PATH = "configs/fgc/fgc-2-sf1-protocol-v3.toml"
PREDECESSORS = {
    "protocol_preflight_config": "configs/fgc/fgc-1-id0-pref1.toml",
    "protocol_preflight_result": "results/fgc-1-id0-pref1.json",
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "uhyp1_config": "configs/fgc/fgc-1-hyp1-dom3-uhyp1.toml",
    "uhyp1_result": "results/fgc-1-hyp1-dom3-uhyp1.json",
    "con2_config": "configs/fgc/fgc-1-hyp1-con2-mprop1.toml",
    "con2_result": "results/fgc-1-hyp1-con2-mprop1.json",
    "con3_config": "configs/fgc/fgc-1-hyp1-con3-cau1.toml",
    "con3_result": "results/fgc-1-hyp1-con3-cau1.json",
    "bnd1_config": "configs/fgc/fgc-1-hyp1-bnd1-md1.toml",
    "bnd1_result": "results/fgc-1-hyp1-bnd1-md1.json",
    "def0_config": "configs/fgc/fgc-1-def0-obs1.toml",
    "def0_result": "results/fgc-1-def0-obs1.json",
    "eft1_config": "configs/fgc/fgc-1-eft1-open1.toml",
    "eft1_result": "results/fgc-1-eft1-open1.json",
}
FUTURE_PATHS = {
    "protocol_holdout": (
        "configs/fgc/fgc-1-pro3-hld1.toml",
        "results/fgc-1-pro3-hld1.json",
    ),
    "run_domain": (
        "configs/fgc/fgc-1-dom4-run1.toml",
        "results/fgc-1-dom4-run1.json",
    ),
    "multidirectional_health": (
        "configs/fgc/fgc-1-hyp2-md1.toml",
        "results/fgc-1-hyp2-md1.json",
    ),
    "nonlinear_source": (
        "configs/fgc/fgc-1-src1-nl1.toml",
        "results/fgc-1-src1-nl1.json",
    ),
    "physical_constraints": (
        "configs/fgc/fgc-1-con4-phy1.toml",
        "results/fgc-1-con4-phy1.json",
    ),
    "regular_center": (
        "configs/fgc/fgc-1-ctr1-reg1.toml",
        "results/fgc-1-ctr1-reg1.json",
    ),
    "initial_data": (
        "configs/fgc/fgc-1-id1-fam1.toml",
        "results/fgc-1-id1-fam1.json",
    ),
    "boundary_control": (
        "configs/fgc/fgc-1-bnd2-cp1.toml",
        "results/fgc-1-bnd2-cp1.json",
    ),
    "health_monitor": (
        "configs/fgc/fgc-1-hlt1-mon1.toml",
        "results/fgc-1-hlt1-mon1.json",
    ),
    "numerical_validation": (
        "configs/fgc/fgc-1-num1-val1.toml",
        "results/fgc-1-num1-val1.json",
    ),
}
AUTHORIZATION_CONTRACT = {
    "all_scoped_predicates_must_pass": True,
    "one_false_scoped_predicate_stops_holdout_execution": True,
    "scoped_classical_diagnostic_is_not_retained_EFT_evolution": True,
    "scoped_classical_diagnostic_is_not_a_physical_transition_claim": True,
    "physical_transition_implies_retained_EFT_and_scoped_authorization": True,
    "no_reverse_implication_is_assumed": True,
    "missing_evidence_is_false_not_inferred": True,
    "bare_config_booleans_cannot_substitute_for_hashed_evidence": True,
}
CLAIMS = {
    "complete_Wilsonian_EFT_defined": False,
    "UV_completion_derived": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "collapse_solution_derived": False,
    "affine_null_defocusing_derived": False,
    "finite_invariant_transition_surface_derived": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/scoped_run_authorization.py",
        "scripts/reproduce_fgc_run1_sym1.py",
    )
)


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _canonical_path(name: str, value: Any, expected: str, *, must_exist: bool) -> Path:
    if not isinstance(value, str) or value != expected:
        raise ValueError(f"{name} must name {expected}")
    path = (REPOSITORY / value).resolve()
    try:
        relative = path.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside repository") from exc
    if relative != value or (must_exist and not path.is_file()):
        suffix = "existing" if must_exist else "declared"
        raise ValueError(f"{name} must be canonical, traversal free, and {suffix}")
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


def _validate_repository_hash_ledger(
    result: Mapping[str, Any], ledger_name: str, name: str
) -> None:
    ledger = result.get(ledger_name)
    if not isinstance(ledger, Mapping) or not ledger:
        raise ValueError(f"{name} {ledger_name} ledger must be nonempty")
    for relative, digest in ledger.items():
        if not isinstance(relative, str) or not isinstance(digest, str):
            raise ValueError(f"{name} {ledger_name} contains a non-string entry")
        path = (REPOSITORY / relative).resolve()
        try:
            canonical = path.relative_to(REPOSITORY).as_posix()
        except ValueError as exc:
            raise ValueError(f"{name} {ledger_name} escapes repository") from exc
        if canonical != relative or not path.is_file() or _sha(path) != digest:
            raise ValueError(f"{name} {ledger_name} is stale or noncanonical for {relative}")


def _validate_future_result_integrity(
    result: Mapping[str, Any], result_path: Path, name: str
) -> None:
    source = result_path.read_text(encoding="utf-8")
    if source != _canonical(result):
        raise ValueError(f"future_evidence.{name} result must use canonical sorted JSON")
    for ledger_name in (
        "source_config_sha256",
        "predecessor_sha256",
        "implementation_sha256",
    ):
        _validate_repository_hash_ledger(result, ledger_name, f"future_evidence.{name}")
    document = result.get("derivation_document")
    document_digest = result.get("derivation_document_sha256")
    if not isinstance(document, str) or not isinstance(document_digest, str):
        raise ValueError(f"future_evidence.{name} derivation binding is absent")
    document_path = (REPOSITORY / document).resolve()
    try:
        canonical_document = document_path.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"future_evidence.{name} derivation document escapes repository") from exc
    if (
        canonical_document != document
        or not document_path.is_file()
        or _sha(document_path) != document_digest
    ):
        raise ValueError(f"future_evidence.{name} derivation binding is stale")
    generated_by = result.get("generated_by")
    implementation = result.get("implementation_sha256")
    if (
        not isinstance(generated_by, str)
        or not isinstance(implementation, Mapping)
        or generated_by not in implementation
    ):
        raise ValueError(f"future_evidence.{name} generator is not implementation-bound")


def _git_blob(commit: str, relative: str, name: str) -> bytes:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{relative}"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise ValueError(f"cannot read immutable {name} blob from {commit}")
    return completed.stdout


def _validate_protocol_diagnosis(
    protocol: Mapping[str, Any], protocol_certificate: Mapping[str, Any]
) -> dict[str, Any]:
    amendment = protocol.get("amendment")
    if not isinstance(amendment, Mapping):
        raise ValueError("PROTO3 amendment is absent")
    commit = amendment.get("diagnosis_checkpoint_commit")
    predecessor_path = amendment.get("predecessor_protocol_config")
    diagnosis_path = amendment.get("diagnosis_result")
    if not all(isinstance(item, str) and item for item in (commit, predecessor_path, diagnosis_path)):
        raise ValueError("PROTO3 immutable diagnosis coordinates are absent")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if ancestor.returncode != 0:
        raise ValueError("PROTO3 diagnosis checkpoint is not an ancestor of HEAD")

    predecessor_blob = _git_blob(commit, predecessor_path, "PROTO2 protocol")
    if sha256(predecessor_blob).hexdigest() != amendment.get(
        "predecessor_protocol_sha256"
    ):
        raise ValueError("PROTO3 immutable PROTO2 blob hash differs")
    try:
        predecessor = tomllib.loads(predecessor_blob.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ValueError("PROTO3 immutable PROTO2 blob is not valid TOML") from exc
    predecessor_certificate = validate_sf1_protocol(predecessor)
    if (
        predecessor_certificate.get("artifact_id") != SF1_PROTOCOL_V2_ARTIFACT_ID
        or predecessor_certificate.get("protocol_version") != 2
    ):
        raise ValueError("PROTO3 immutable predecessor is not sealed PROTO2")

    diagnosis_blob = _git_blob(commit, diagnosis_path, "ID1 diagnosis")
    if sha256(diagnosis_blob).hexdigest() != amendment.get("diagnosis_result_sha256"):
        raise ValueError("PROTO3 immutable ID1 diagnosis blob hash differs")
    try:
        diagnosis = json.loads(
            diagnosis_blob.decode("utf-8"), object_pairs_hook=_pairs
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("PROTO3 immutable ID1 diagnosis is not unique-key JSON") from exc
    if not isinstance(diagnosis, dict) or diagnosis_blob.decode("utf-8") != _canonical(
        diagnosis
    ):
        raise ValueError("PROTO3 immutable ID1 diagnosis is not canonical JSON")

    status = diagnosis.get("gate_status")
    scope = diagnosis.get("scope_bindings")
    payload = diagnosis.get("artifact_payload")
    quantitative = payload.get("quantitative_evidence") if isinstance(payload, Mapping) else None
    family_box = quantitative.get("family_box") if isinstance(quantitative, Mapping) else None
    epistemic = quantitative.get("epistemic_boundary") if isinstance(quantitative, Mapping) else None
    if (
        diagnosis.get("artifact_id") != "FGC-1-ID1-FAM1"
        or diagnosis.get("classification")
        != "finite_mass_constraint_compatible_data_family_certificate"
        or not isinstance(status, Mapping)
        or status.get(
            "nonzero_width_finite_mass_constraint_compatible_family_constructed"
        )
        is not True
        or not isinstance(scope, Mapping)
        or scope.get("protocol_artifact_id") != SF1_PROTOCOL_V2_ARTIFACT_ID
        or scope.get("protocol_config_sha256")
        != amendment.get("predecessor_protocol_sha256")
        or not isinstance(family_box, Mapping)
        or family_box.get("PROTO2_width_nine_eighths_case_qualified_here")
        is not False
        or not isinstance(epistemic, Mapping)
        or epistemic.get("FGCQR_holdout_outcomes_inspected") is not False
        or epistemic.get("time_evolution_performed") is not False
        or epistemic.get("Raychaudhuri_defocusing_tested") is not False
        or protocol_certificate.get("diagnosis_artifact_id")
        != "FGC-1-ID1-FAM1"
        or protocol_certificate.get("PROTO1_and_PROTO2_preserved") is not True
        or protocol_certificate.get("premise_revision_only") is not True
    ):
        raise ValueError("PROTO3 is not bound to the immutable premise-only ID1 diagnosis")
    return {
        "artifact_id": "FGC-1-ID1-FAM1",
        "historical_checkpoint_commit": commit,
        "immutable_PROTO2_config_sha256": amendment["predecessor_protocol_sha256"],
        "immutable_ID1_result_sha256": amendment["diagnosis_result_sha256"],
        "PROTO2_width_nine_eighths_case_unqualified_at_checkpoint": True,
        "PROTO3_revision_is_premise_only": True,
        "FGCQR_evolution_outcomes_inspected_before_revision": False,
        "SGBL_evolution_outcomes_inspected_before_revision": False,
    }


def _validate_protocol_preflight(
    result: Mapping[str, Any],
    result_path: Path,
    config_path: Path,
    protocol: Mapping[str, Any],
    protocol_certificate: Mapping[str, Any],
) -> dict[str, Any]:
    if result_path.read_text(encoding="utf-8") != _canonical(result):
        raise ValueError("protocol preflight result must use canonical sorted JSON")
    for ledger_name in (
        "source_config_sha256",
        "predecessor_sha256",
        "implementation_sha256",
    ):
        _validate_repository_hash_ledger(result, ledger_name, "protocol_preflight")
    _binds(result, config_path, "protocol_preflight")
    if (
        result.get("artifact_id") != "FGC-1-ID0-PREF1"
        or result.get("project_version") != PROJECT_VERSION
        or result.get("classification")
        != "pre_holdout_GR0_protocol_obstruction_not_an_FGCQR_mechanism_result"
    ):
        raise ValueError("protocol preflight identity or classification differs")
    status = result.get("gate_status")
    scope = result.get("scope_bindings")
    history = result.get("historical_preflight_checkpoint")
    certificate = result.get("certificate_contract")
    if (
        not isinstance(status, Mapping)
        or status.get("PROTO1_initial_data_preflight_obstruction_verified") is not True
        or status.get("PROTO1_FGCQR_holdout_execution_authorized") is not False
        or not isinstance(scope, Mapping)
        or scope.get("target_protocol_artifact_id") != "FGC-2-SF1-PROTO1"
        or scope.get("FGCQR_holdout_outcomes_inspected") is not False
        or not isinstance(history, Mapping)
        or history.get("original_ID0_canonical_ledgers_verified_at_commit") is not True
        or history.get("original_RUN1_two_of_eight_stop_verified_at_commit") is not True
        or not isinstance(certificate, Mapping)
        or certificate.get("historical_pre_PROTO2_checkpoint_verified") is not True
        or protocol_certificate.get("preflight_artifact_id") != "FGC-1-ID0-PREF1"
        or protocol_certificate.get("PROTO2_amendment_validated") is not True
        or protocol_certificate.get("PROTO1_and_PROTO2_preserved") is not True
        or protocol_certificate.get("premise_revision_only") is not True
    ):
        raise ValueError(
            "PROTO3 does not preserve the verified PROTO2/ID0 premise-repair lineage"
        )
    return {
        "artifact_id": "FGC-1-ID0-PREF1",
        "historical_checkpoint_commit": history["commit"],
        "PROTO1_initial_data_preflight_obstruction_verified": True,
        "PROTO2_revision_is_premise_only": True,
        "PROTO3_preserves_PROTO2_ID0_lineage": True,
        "FGCQR_outcomes_inspected_before_revision": False,
        "SGBL_outcomes_inspected_before_revision": False,
    }


def _expected_scope_bindings(
    config_path: Path,
    config: Mapping[str, Any],
    protocol_certificate: Mapping[str, Any],
) -> dict[str, str]:
    paths = config["predecessor_paths"]
    return {
        "run1_artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "run1_config_sha256": _sha(config_path),
        "protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "protocol_config_sha256": _sha(config["protocol_path"]),
        "protocol_semantic_holdout_contract_sha256": protocol_certificate[
            "semantic_holdout_contract_sha256"
        ],
        "action_artifact_id": "FGC-1-ACT1",
        "action_config_sha256": _sha(paths["action_config"]),
        "action_result_sha256": _sha(paths["action_result"]),
        "variation_artifact_id": "FGC-1-VAR1",
        "variation_config_sha256": _sha(paths["variation_config"]),
        "variation_result_sha256": _sha(paths["variation_result"]),
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
    }


def _binds(result: Mapping[str, Any], config: Path, name: str) -> None:
    relative = _rel(config)
    digest = _sha(config)
    if (
        result.get("source_config") == relative
        and result.get("source_config_sha256") == digest
    ):
        return
    for ledger_name in ("source_config_sha256", "source_configs_sha256"):
        ledger = result.get(ledger_name)
        if isinstance(ledger, Mapping) and ledger.get(relative) == digest:
            return
    source_configs = result.get("source_configs")
    if isinstance(source_configs, Mapping):
        nested = source_configs.get(relative)
        if isinstance(nested, Mapping) and nested.get("sha256") == digest:
            return
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
            "protocol_config",
            *PREDECESSORS,
            "authorization_contract",
            "future_evidence",
            "claims",
        },
    )
    if {
        "schema_version": raw["schema_version"],
        "artifact_id": raw["artifact_id"],
        "project_version": raw["project_version"],
    } != {
        "schema_version": 1,
        "artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
    }:
        raise ValueError("RUN1-SYM1 identity or version differs")
    protocol_path = _canonical_path(
        "protocol_config", raw["protocol_config"], PROTOCOL_PATH, must_exist=True
    )
    predecessor_paths = {
        name: _canonical_path(name, raw[name], expected, must_exist=True)
        for name, expected in PREDECESSORS.items()
    }
    if raw["authorization_contract"] != AUTHORIZATION_CONTRACT:
        raise ValueError("RUN1-SYM1 authorization contract differs")
    if raw["claims"] != CLAIMS:
        raise ValueError("RUN1-SYM1 claims must remain fail-closed")

    future = raw["future_evidence"]
    if not isinstance(future, Mapping) or set(future) != set(FUTURE_EVIDENCE_REQUIREMENTS):
        raise ValueError("RUN1-SYM1 future evidence slots differ")
    future_paths: dict[str, dict[str, Path]] = {}
    for name, (expected_id, expected_gate) in FUTURE_EVIDENCE_REQUIREMENTS.items():
        specification = future[name]
        if not isinstance(specification, Mapping):
            raise ValueError(f"future_evidence.{name} must be a table")
        _keys(
            f"future_evidence.{name}",
            specification,
            {"artifact_id", "config", "result", "required_gate"},
        )
        expected_config, expected_result = FUTURE_PATHS[name]
        if (
            specification["artifact_id"] != expected_id
            or specification["required_gate"] != expected_gate
        ):
            raise ValueError(f"future_evidence.{name} identity or gate differs")
        future_paths[name] = {
            "config": _canonical_path(
                f"future_evidence.{name}.config",
                specification["config"],
                expected_config,
                must_exist=False,
            ),
            "result": _canonical_path(
                f"future_evidence.{name}.result",
                specification["result"],
                expected_result,
                must_exist=False,
            ),
        }
    return {
        "raw": raw,
        "protocol_path": protocol_path,
        "predecessor_paths": predecessor_paths,
        "future_paths": future_paths,
        "source_sha256": sha256(source).hexdigest(),
    }


def _load_future_records(
    paths: Mapping[str, Mapping[str, Path]],
    expected_scope: Mapping[str, str] | None = None,
) -> tuple[dict[str, dict[str, Any] | None], dict[str, str]]:
    records: dict[str, dict[str, Any] | None] = {}
    hashes: dict[str, str] = {}
    for name in FUTURE_EVIDENCE_REQUIREMENTS:
        config_path = paths[name]["config"]
        result_path = paths[name]["result"]
        config_exists = config_path.is_file()
        result_exists = result_path.is_file()
        if config_exists is not result_exists:
            raise ValueError(
                f"future_evidence.{name} config and result must either both exist or both be absent"
            )
        if not config_exists:
            records[name] = None
            continue
        if expected_scope is None:
            raise ValueError("future evidence scope hashes are required for present records")
        record = _load_unique_json(result_path, f"future_evidence.{name}")
        _validate_future_result_integrity(record, result_path, name)
        _binds(record, config_path, f"future_evidence.{name}")
        scope = record.get("scope_bindings")
        if not isinstance(scope, Mapping) or any(
            scope.get(key) != value for key, value in expected_scope.items()
        ):
            raise ValueError(f"future_evidence.{name} scope hash binding differs")
        records[name] = record
        hashes[_rel(config_path)] = _sha(config_path)
        hashes[_rel(result_path)] = _sha(result_path)
    return records, hashes


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    records = {
        name: _load_unique_json(
            config["predecessor_paths"][f"{name}_result"], name
        )
        for name in ("action", "variation", "uhyp1", "con2", "con3", "bnd1", "def0", "eft1")
    }
    for name in records:
        _binds(records[name], config["predecessor_paths"][f"{name}_config"], name)
    protocol = tomllib.loads(config["protocol_path"].read_text(encoding="utf-8"))
    protocol_certificate = validate_sf1_protocol(protocol)
    diagnosis_summary = _validate_protocol_diagnosis(protocol, protocol_certificate)
    preflight_result_path = config["predecessor_paths"]["protocol_preflight_result"]
    preflight = _load_unique_json(preflight_result_path, "protocol_preflight")
    preflight_summary = _validate_protocol_preflight(
        preflight,
        preflight_result_path,
        config["predecessor_paths"]["protocol_preflight_config"],
        protocol,
        protocol_certificate,
    )
    expected_scope = _expected_scope_bindings(
        config_path, config, protocol_certificate
    )
    future_records, future_hashes = _load_future_records(
        config["future_paths"], expected_scope
    )
    audit = scoped_classical_spherical_run_audit(
        **records,
        protocol=protocol,
        future_evidence=future_records,
    )
    if audit["gate_status"]["retained_EFT_evolution_authorized"] is not False:
        raise ValueError("RUN1-SYM1 cannot relax the negative EFT1 authorization")
    if audit["gate_status"]["physical_transition_claim_authorized"] is not False:
        raise ValueError("RUN1-SYM1 cannot authorize a physical transition")
    predecessor_hashes = {
        _rel(path): _sha(path) for path in config["predecessor_paths"].values()
    }
    predecessor_hashes.update(future_hashes)
    return {
        "schema_version": 1,
        "artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "classification": audit["classification"],
        "project_version": PROJECT_VERSION,
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): config["source_sha256"],
            _rel(config["protocol_path"]): _sha(config["protocol_path"]),
        },
        "predecessor_sha256": predecessor_hashes,
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "frozen_configuration": config["raw"],
        "protocol_amendment_diagnosis": diagnosis_summary,
        "protocol_amendment_preflight": preflight_summary,
        "scoped_run_authorization_audit": audit,
        "gate_status": audit["gate_status"],
        "nonclaims": CLAIMS,
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    try:
        source = path.read_text(encoding="utf-8")
        value = json.loads(source, object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("RUN1-SYM1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("RUN1-SYM1 result must use canonical sorted JSON")
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
