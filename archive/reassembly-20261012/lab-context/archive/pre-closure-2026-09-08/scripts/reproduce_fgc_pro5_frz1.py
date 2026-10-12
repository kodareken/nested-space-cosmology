#!/usr/bin/env python3
"""Regenerate the immutable, outcome-neutral PROTO5 freeze certificate."""

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
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro5-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro5-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro5-frz1.md"

from recursive_horizons.fgc.evolution.protocol_v5 import (  # noqa: E402
    PROTO5_AMENDMENT,
    PROTO5_CLAIMS,
    SF1_PROTOCOL_V5_ARTIFACT_ID,
    validate_sf1_protocol_v5,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    SF1_PROTOCOL_V4_ARTIFACT_ID,
    validate_sf1_protocol,
)


ARTIFACT_ID = "FGC-1-PRO5-FRZ1"
PROJECT_VERSION = "0.11.0"
EXPECTED_SCOPE = {
    "target_protocol": SF1_PROTOCOL_V5_ARTIFACT_ID,
    "predecessor_protocol": SF1_PROTOCOL_V4_ARTIFACT_ID,
    "diagnosis_artifact": "FGC-1-CAL1-PREF3",
    "calibration_branch": "GR-0",
    "holdout_branch": "FGC-QR",
    "freeze_role": "premise_only_pre_trajectory_semidiscrete_protocol_and_lineage_certificate",
    "initial_common_event_diagnostics_disclosed": True,
    "fresh_GR0_calibration_trajectory_inspected": False,
    "FGCQR_evolution_outcomes_inspected": False,
    "SGBL_evolution_outcomes_inspected": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": "0781d7d218e9c129858b134c3181fad07ac692b0",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v4.toml",
    "predecessor_protocol_sha256": "02a25337e040b97fdd115bc7da4aa87fd09755366599bcade064644ffb3cb87a",
    "diagnosis_config": "configs/fgc/fgc-1-cal1-pref3.toml",
    "diagnosis_config_sha256": "514861c679272cd182da9674e7640734a427c0349157c03f0899cbe439b475ff",
    "diagnosis_result": "results/fgc-1-cal1-pref3.json",
    "diagnosis_result_sha256": "3ea43baf0fb9830b348fc49a79af5c535d8f7089aa03a957a22f6927528a4395",
}
EXPECTED_PROOF = {
    "PROTO5_overlay_must_validate_exactly",
    "checkpoint_must_be_ancestor_of_HEAD",
    "predecessor_and_diagnosis_must_be_read_from_checkpoint",
    "immutable_blob_hashes_must_match_PROTO5_amendment",
    "immutable_PROTO4_must_validate_unchanged",
    "immutable_CAL1_must_be_canonical_and_premise_only",
    "physical_experiment_and_outcome_rules_must_be_inherited",
    "only_q_may_be_projected_and_u_p_must_be_bitwise_preserved",
    "raw_residuals_must_decide_magnitude_guards",
    "roundoff_enclosure_must_be_public_nonnegative_and_zero_classification_only",
    "method_owned_guards_and_finest_pair_rule_must_be_frozen",
    "new_output_namespaces_must_be_declared_but_not_accessed",
    "canonical_hash_bound_result_required",
}
EXPECTED_CLAIMS = dict(PROTO5_CLAIMS)
IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v5.py",
    REPOSITORY / "src/recursive_horizons/fgc/scoped_run_authorization.py",
    REPOSITORY / "scripts/reproduce_fgc_pro5_frz1.py",
)


def _canonical(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def _unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def _repo_path(name: str, value: object) -> Path:
    if not isinstance(value, str) or not value:
        raise TypeError(f"{name} must be a repository-relative path")
    path = (REPOSITORY / value).resolve()
    try:
        relative = path.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} escapes the repository") from exc
    if relative != value or not path.is_file():
        raise ValueError(f"{name} is absent or noncanonical")
    return path


def _load_toml_bytes(source: bytes, name: str) -> dict[str, Any]:
    try:
        value = tomllib.loads(source.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ValueError(f"{name} is not valid UTF-8 TOML") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} root must be a table")
    return value


def _load_json_bytes(source: bytes, name: str) -> dict[str, Any]:
    try:
        value = json.loads(source.decode("utf-8"), object_pairs_hook=_unique_pairs)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{name} is not valid unique-key JSON") from exc
    if not isinstance(value, dict) or source.decode("utf-8") != _canonical(value):
        raise ValueError(f"{name} is not canonical sorted JSON")
    return value


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json_bytes(path.read_bytes(), "PRO5-FRZ1 result")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    source = path.read_bytes()
    raw = _load_toml_bytes(source, "PRO5-FRZ1 config")
    expected_root = {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        "protocol_config",
        "scope",
        "immutable_lineage",
        "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("PRO5-FRZ1 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("PRO5-FRZ1 identity or convention differs")
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("PRO5-FRZ1 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("PRO5-FRZ1 immutable lineage differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("PRO5-FRZ1 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("PRO5-FRZ1 claims must remain fail-closed")

    protocol_path = _repo_path("protocol_config", raw["protocol_config"])
    protocol_source = protocol_path.read_bytes()
    protocol_raw = _load_toml_bytes(protocol_source, "PROTO5")
    protocol_certificate = validate_sf1_protocol_v5(protocol_raw)
    if protocol_raw["amendment"] != PROTO5_AMENDMENT:
        raise ValueError("PROTO5 amendment differs from its validator")
    amendment_lineage = {
        "checkpoint_commit": protocol_raw["amendment"]["diagnosis_checkpoint_commit"],
        "predecessor_protocol_config": protocol_raw["amendment"]["predecessor_protocol_config"],
        "predecessor_protocol_sha256": protocol_raw["amendment"]["predecessor_protocol_sha256"],
        "diagnosis_config": protocol_raw["amendment"]["diagnosis_config"],
        "diagnosis_config_sha256": protocol_raw["amendment"]["diagnosis_config_sha256"],
        "diagnosis_result": protocol_raw["amendment"]["diagnosis_result"],
        "diagnosis_result_sha256": protocol_raw["amendment"]["diagnosis_result_sha256"],
    }
    if amendment_lineage != EXPECTED_LINEAGE:
        raise ValueError("PROTO5 amendment and PRO5-FRZ1 lineage differ")
    return {
        "raw": raw,
        "source_sha256": sha256(source).hexdigest(),
        "protocol_path": protocol_path,
        "protocol_source": protocol_source,
        "protocol_raw": protocol_raw,
        "protocol_certificate": protocol_certificate,
    }


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


def _immutable_lineage(configuration: Mapping[str, Any]) -> dict[str, Any]:
    lineage = configuration["raw"]["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if ancestor.returncode != 0:
        raise ValueError("PROTO5 diagnosis checkpoint is not an ancestor of HEAD")

    predecessor_blob = _git_blob(
        commit, lineage["predecessor_protocol_config"], "PROTO4 protocol"
    )
    if sha256(predecessor_blob).hexdigest() != lineage["predecessor_protocol_sha256"]:
        raise ValueError("immutable PROTO4 blob hash differs")
    predecessor = validate_sf1_protocol(
        _load_toml_bytes(predecessor_blob, "immutable PROTO4")
    )
    if (
        predecessor["artifact_id"] != SF1_PROTOCOL_V4_ARTIFACT_ID
        or predecessor["protocol_version"] != 4
        or predecessor["frozen"] is not True
        or predecessor["outcome_neutral_contract_validated"] is not True
        or not isinstance(predecessor.get("claims"), Mapping)
        or any(value is not False for value in predecessor["claims"].values())
    ):
        raise ValueError("immutable PROTO4 is not frozen and fail-closed")

    diagnosis_config_blob = _git_blob(
        commit, lineage["diagnosis_config"], "CAL1 config"
    )
    if sha256(diagnosis_config_blob).hexdigest() != lineage["diagnosis_config_sha256"]:
        raise ValueError("immutable CAL1 config hash differs")
    diagnosis_blob = _git_blob(
        commit, lineage["diagnosis_result"], "CAL1 diagnosis"
    )
    if sha256(diagnosis_blob).hexdigest() != lineage["diagnosis_result_sha256"]:
        raise ValueError("immutable CAL1 result hash differs")
    diagnosis = _load_json_bytes(diagnosis_blob, "immutable CAL1")
    gates = diagnosis.get("gate_status", {})
    payload = diagnosis.get("artifact_payload", {})
    boundary = payload.get("epistemic_boundary", {})
    candidates = payload.get("candidate_compositions", [])
    prospective = [
        method.get("projected_state_under_prospective_repair", {}).get(
            "admission_passed"
        )
        for candidate in candidates
        if isinstance(candidate, Mapping)
        for method in candidate.get("methods", [])
        if isinstance(method, Mapping)
    ]
    if (
        diagnosis.get("artifact_id") != "FGC-1-CAL1-PREF3"
        or diagnosis.get("classification")
        != "pre_trajectory_semidiscrete_test_contract_obstruction_not_mechanism_result"
        or gates.get("PROTO4_semidiscrete_common_event_contract_obstructed") is not True
        or gates.get("PROTO5_premise_revision_required") is not True
        or gates.get("PROTO4_fresh_GR0_dynamic_calibration_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or boundary.get("trajectory_read") is not False
        or boundary.get("FGCQR_outcome_read") is not False
        or boundary.get("mechanism_question_answered") is not False
        or [candidate.get("amplitude") for candidate in candidates] != ["5/2", "3"]
        or prospective != [True, True, True, True]
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable CAL1 is not the premise-only diagnosis")
    config_hash_ledger = diagnosis.get("source_config_sha256", {})
    if config_hash_ledger.get(lineage["diagnosis_config"]) != lineage[
        "diagnosis_config_sha256"
    ]:
        raise ValueError("immutable CAL1 does not bind its diagnosis config")

    for relative, expected, name in (
        (
            lineage["predecessor_protocol_config"],
            lineage["predecessor_protocol_sha256"],
            "current PROTO4",
        ),
        (
            lineage["diagnosis_config"],
            lineage["diagnosis_config_sha256"],
            "current CAL1 config",
        ),
        (
            lineage["diagnosis_result"],
            lineage["diagnosis_result_sha256"],
            "current CAL1 result",
        ),
    ):
        if sha256((REPOSITORY / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f"{name} differs from its immutable checkpoint blob")

    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "predecessor_protocol": {
            "artifact_id": predecessor["artifact_id"],
            "protocol_version": predecessor["protocol_version"],
            "blob_sha256": lineage["predecessor_protocol_sha256"],
            "semantic_holdout_contract_sha256": predecessor[
                "semantic_holdout_contract_sha256"
            ],
            "claims_all_false": all(
                value is False for value in predecessor["claims"].values()
            ),
        },
        "diagnosis": {
            "artifact_id": diagnosis["artifact_id"],
            "classification": diagnosis["classification"],
            "config_blob_sha256": lineage["diagnosis_config_sha256"],
            "result_blob_sha256": lineage["diagnosis_result_sha256"],
            "PROTO4_obstruction_verified": True,
            "premise_revision_required": True,
            "eligible_amplitudes": ["5/2", "3"],
            "prospective_t0_compositions_passed": 4,
            "trajectory_read": False,
            "FGCQR_outcome_read": False,
        },
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(path)
    lineage = _immutable_lineage(configuration)
    protocol = configuration["protocol_certificate"]
    protocol_relative = configuration["protocol_path"].relative_to(REPOSITORY).as_posix()
    raw = configuration["raw"]
    gate_status = {
        "PROTO5_outcome_neutral_protocol_frozen": True,
        "PROTO5_immutable_lineage_verified": True,
        "PROTO5_native_semidiscrete_initialization_contract_frozen": True,
        "PROTO5_method_owned_constraint_admission_frozen": True,
        "PROTO5_roundoff_zero_classification_contract_frozen": True,
        **raw["claims"],
    }
    nonclaims = {
        "fresh_GR0_calibration_completed": False,
        "SGBL_comparison_completed": False,
        "FGCQR_holdout_evolved": False,
        "trapped_sphere_observed": False,
        "regulator_activation_observed": False,
        "positive_Raychaudhuri_margin_observed": False,
        "FGCQR_mechanism_rejected": False,
        "general_gradient_route_rejected": False,
        "retained_EFT_validity": False,
        "physical_transition": False,
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "immutable_outcome_neutral_PROTO5_semidiscrete_protocol_freeze_not_run_authorization",
        "generated_by": "scripts/reproduce_fgc_pro5_frz1.py",
        "scope_bindings": dict(raw["scope"]),
        "source_config_sha256": {
            path.relative_to(REPOSITORY).as_posix(): configuration["source_sha256"],
            protocol_relative: sha256(configuration["protocol_source"]).hexdigest(),
        },
        "predecessor_sha256": {
            key: value
            for key, value in (
                (
                    raw["immutable_lineage"]["predecessor_protocol_config"],
                    raw["immutable_lineage"]["predecessor_protocol_sha256"],
                ),
                (
                    raw["immutable_lineage"]["diagnosis_config"],
                    raw["immutable_lineage"]["diagnosis_config_sha256"],
                ),
                (
                    raw["immutable_lineage"]["diagnosis_result"],
                    raw["immutable_lineage"]["diagnosis_result_sha256"],
                ),
            )
        },
        "implementation_sha256": {
            item.relative_to(REPOSITORY).as_posix(): sha256(item.read_bytes()).hexdigest()
            for item in IMPLEMENTATION
        },
        "derivation_document": OWNER_DOCUMENT.relative_to(REPOSITORY).as_posix(),
        "derivation_document_sha256": sha256(OWNER_DOCUMENT.read_bytes()).hexdigest(),
        "certificate_contract": sorted(raw["proof_contract"]),
        "gate_status": gate_status,
        "artifact_payload": {
            "protocol_certificate": protocol,
            "immutable_lineage": lineage,
            "frozen_repair": {
                "only_q_is_projected": True,
                "u_and_p_must_be_bitwise_preserved": True,
                "raw_norms_decide_magnitude_guards": True,
                "roundoff_operation_budget": 4096,
                "roundoff_enclosure_is_zero_classification_only": True,
                "primary_guards": {"coarsest": "1/50", "finest": "1/1000"},
                "comparator_guards": {"coarsest": "1/10", "finest": "1/200"},
                "minimum_finest_pair_order": "3/2",
                "all_three_resolutions_must_refine_monotonically": True,
                "eligible_amplitudes_in_order": ["5/2", "3"],
            },
            "epistemic_boundary": {
                "initial_common_event_diagnostics_were_seen": True,
                "fresh_calibration_trajectory_read": False,
                "FGCQR_outcome_read": False,
                "SGBL_outcome_read": False,
                "mechanism_question_answered": False,
                "physical_experiment_changed": False,
            },
            "decision": "PROTO5_is_frozen_but_HLT3_fresh_calibration_and_every_candidate_outcome_remain_closed",
        },
        "nonclaims": nonclaims,
    }


def write(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_canonical(payload), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = record(args.config.resolve())
    write(args.output.resolve(), payload)
    print(
        f"wrote {args.output} (PROTO5 frozen=true; fresh calibration authorized=false; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
