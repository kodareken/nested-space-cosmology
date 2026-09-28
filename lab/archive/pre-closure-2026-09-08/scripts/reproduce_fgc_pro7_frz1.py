#!/usr/bin/env python3
"""Regenerate the immutable, outcome-neutral PROTO7 freeze certificate."""

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
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro7-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro7-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro7-frz1.md"

from recursive_horizons.fgc.evolution.protocol_v6 import (  # noqa: E402
    SF1_PROTOCOL_V6_ARTIFACT_ID,
    validate_sf1_protocol_v6,
)
from recursive_horizons.fgc.evolution.protocol_v7 import (  # noqa: E402
    PROTO7_AMENDMENT,
    PROTO7_CLAIMS,
    SF1_PROTOCOL_V7_ARTIFACT_ID,
    validate_sf1_protocol_v7,
)


ARTIFACT_ID = "FGC-1-PRO7-FRZ1"
PROJECT_VERSION = "0.11.0"
EXPECTED_SCOPE = {
    "target_protocol": SF1_PROTOCOL_V7_ARTIFACT_ID,
    "predecessor_protocol": SF1_PROTOCOL_V6_ARTIFACT_ID,
    "diagnosis_artifact": "FGC-1-CAL3-PREF5",
    "calibration_branch": "GR-0",
    "holdout_branch": "FGC-QR",
    "freeze_role": "premise_only_post_calibration_general_unaccepted_proposal_ownership_protocol_and_lineage_certificate",
    "PROTO6_GR0_calibration_trajectory_disclosed": True,
    "first_nonzero_common_event_completed": False,
    "FGCQR_evolution_outcomes_inspected": False,
    "SGBL_evolution_outcomes_inspected": False,
    "mechanism_question_answered": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": "27c63e7aca1439e75e18b576c6da414e458280ff",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v6.toml",
    "predecessor_protocol_sha256": "e50093914f7b84a88bd593361f8b3eb43942cf9a8502f637fc38d8e61dd82448",
    "diagnosis_config": "configs/fgc/fgc-1-cal3-pref5.toml",
    "diagnosis_config_sha256": "cec13d7e50d2546d423e7ffea6690226d71fbe5174073e78ffd2ac37f9e44275",
    "diagnosis_result": "results/fgc-1-cal3-pref5.json",
    "diagnosis_result_sha256": "669ee5d9a17a8d3e0b213b06fcb3ff426e0908b331fe38429b5704bbbf49080f",
}
EXPECTED_PROOF = {
    "PROTO7_overlay_must_validate_exactly",
    "checkpoint_must_be_ancestor_of_HEAD",
    "predecessor_and_diagnosis_must_be_read_from_checkpoint",
    "immutable_blob_hashes_must_match_PROTO7_amendment",
    "immutable_PROTO6_must_validate_unchanged",
    "immutable_CAL3_must_be_canonical_and_nonmechanism",
    "accepted_state_source_failure_must_remain_terminal",
    "every_source_only_unaccepted_proposal_failure_may_retry",
    "candidate_endpoint_must_remain_unaccepted_until_commit",
    "any_non_source_failure_must_veto_retry",
    "raw_source_threshold_solver_and_retry_budget_must_remain_unchanged",
    "all_non_source_physical_and_numerical_gates_must_remain_unchanged",
    "complete_retry_terminal_and_exhaustion_evidence_must_be_required",
    "new_output_namespaces_must_be_declared_but_not_accessed",
    "canonical_hash_bound_result_required",
}
IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v6.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v7.py",
    Path(__file__).resolve(),
)


def _canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n"


def _unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


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


def _repo_file(name: str, value: object) -> Path:
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


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json_bytes(path.read_bytes(), "PRO7-FRZ1 result")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    source = path.read_bytes()
    raw = _load_toml_bytes(source, "PRO7-FRZ1 config")
    expected_root = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", "protocol_config", "scope", "immutable_lineage",
        "proof_contract", "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("PRO7-FRZ1 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("PRO7-FRZ1 identity or convention differs")
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("PRO7-FRZ1 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("PRO7-FRZ1 immutable lineage differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("PRO7-FRZ1 proof contract differs")
    if raw["claims"] != PROTO7_CLAIMS:
        raise ValueError("PRO7-FRZ1 claims must remain fail-closed")

    protocol_path = _repo_file("protocol_config", raw["protocol_config"])
    protocol_source = protocol_path.read_bytes()
    protocol_raw = _load_toml_bytes(protocol_source, "PROTO7")
    protocol_certificate = validate_sf1_protocol_v7(protocol_raw)
    if protocol_raw["amendment"] != PROTO7_AMENDMENT:
        raise ValueError("PROTO7 amendment differs from its validator")
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
        raise ValueError("PROTO7 amendment and PRO7-FRZ1 lineage differ")
    return {
        "raw": raw,
        "source_sha256": sha256(source).hexdigest(),
        "protocol_path": protocol_path,
        "protocol_source": protocol_source,
        "protocol_certificate": protocol_certificate,
    }


def _immutable_lineage(configuration: Mapping[str, Any]) -> dict[str, Any]:
    lineage = configuration["raw"]["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("PROTO7 diagnosis checkpoint is not an ancestor of HEAD")

    predecessor_blob = _git_blob(commit, lineage["predecessor_protocol_config"], "PROTO6 protocol")
    if sha256(predecessor_blob).hexdigest() != lineage["predecessor_protocol_sha256"]:
        raise ValueError("immutable PROTO6 blob hash differs")
    predecessor = validate_sf1_protocol_v6(_load_toml_bytes(predecessor_blob, "immutable PROTO6"))
    if (
        predecessor["artifact_id"] != SF1_PROTOCOL_V6_ARTIFACT_ID
        or predecessor["protocol_version"] != 6
        or predecessor["frozen"] is not True
        or predecessor["outcome_neutral_contract_validated"] is not True
        or any(value is not False for value in predecessor["claims"].values())
    ):
        raise ValueError("immutable PROTO6 is not frozen and fail-closed")

    diagnosis_config_blob = _git_blob(commit, lineage["diagnosis_config"], "CAL3 config")
    diagnosis_blob = _git_blob(commit, lineage["diagnosis_result"], "CAL3 diagnosis")
    if sha256(diagnosis_config_blob).hexdigest() != lineage["diagnosis_config_sha256"]:
        raise ValueError("immutable CAL3 config hash differs")
    if sha256(diagnosis_blob).hexdigest() != lineage["diagnosis_result_sha256"]:
        raise ValueError("immutable CAL3 result hash differs")
    diagnosis = _load_json_bytes(diagnosis_blob, "immutable CAL3")
    gates = diagnosis.get("gate_status", {})
    payload = diagnosis.get("artifact_payload", {})
    replay = payload.get("deterministic_terminal_replay", {})
    smaller = replay.get("one_further_halving_control", {})
    decision = payload.get("decision", {})
    boundary = payload.get("epistemic_boundary", {})
    if (
        diagnosis.get("artifact_id") != "FGC-1-CAL3-PREF5"
        or diagnosis.get("classification")
        != "post_PROTO6_GR0_unaccepted_candidate_endpoint_source_stop_ownership_diagnosis"
        or gates.get("PROTO6_unaccepted_candidate_endpoint_stop_ownership_contract_obstructed") is not True
        or gates.get("PROTO7_general_transaction_revision_required") is not True
        or gates.get("fresh_GR0_dynamic_calibration_completed") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or replay.get("last_accepted_state_passes_raw_gate") is not True
        or replay.get("all_terminal_failures_are_source_only") is not True
        or replay.get("failed_terminal_stages") != ["rk4_k4", "candidate_endpoint"]
        or smaller.get("accepted") is not True
        or smaller.get("retry_returned") is not False
        or decision.get("new_protocol_required") != SF1_PROTOCOL_V7_ARTIFACT_ID
        or boundary.get("first_nonzero_common_event_completed") is not False
        or boundary.get("mechanism_question_answered") is not False
        or boundary.get("candidate_or_general_gradient_route_rejected") is not False
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable CAL3 is not the bounded nonmechanism diagnosis")
    if diagnosis.get("source_config_sha256", {}).get(lineage["diagnosis_config"]) != lineage["diagnosis_config_sha256"]:
        raise ValueError("immutable CAL3 does not bind its diagnosis config")

    for relative, expected, name in (
        (lineage["predecessor_protocol_config"], lineage["predecessor_protocol_sha256"], "current PROTO6"),
        (lineage["diagnosis_config"], lineage["diagnosis_config_sha256"], "current CAL3 config"),
        (lineage["diagnosis_result"], lineage["diagnosis_result_sha256"], "current CAL3 result"),
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
            "semantic_holdout_contract_sha256": predecessor["semantic_holdout_contract_sha256"],
            "claims_all_false": all(value is False for value in predecessor["claims"].values()),
        },
        "diagnosis": {
            "artifact_id": diagnosis["artifact_id"],
            "classification": diagnosis["classification"],
            "config_blob_sha256": lineage["diagnosis_config_sha256"],
            "result_blob_sha256": lineage["diagnosis_result_sha256"],
            "PROTO6_transaction_ownership_obstructed": True,
            "PROTO7_general_transaction_revision_required": True,
            "accepted_state_source_gate_passed": True,
            "terminal_failures_source_only": True,
            "failed_terminal_stages": ["rk4_k4", "candidate_endpoint"],
            "one_further_halving_complete_transaction_passed": True,
            "first_nonzero_common_event_completed": False,
            "FGCQR_outcome_read": False,
            "mechanism_question_answered": False,
        },
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(path)
    lineage = _immutable_lineage(configuration)
    protocol = configuration["protocol_certificate"]
    raw = configuration["raw"]
    protocol_relative = configuration["protocol_path"].relative_to(REPOSITORY).as_posix()
    gate_status = {
        "PROTO7_outcome_neutral_protocol_frozen": True,
        "PROTO7_immutable_lineage_verified": True,
        "PROTO7_accepted_state_source_terminal_contract_frozen": True,
        "PROTO7_unaccepted_proposal_source_retry_contract_frozen": True,
        "PROTO7_non_source_retry_veto_contract_frozen": True,
        "PROTO7_complete_failure_observability_contract_frozen": True,
        "PROTO7_non_source_gates_inherited_unchanged": True,
        **raw["claims"],
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "immutable_outcome_neutral_PROTO7_general_transaction_stop_ownership_protocol_freeze_not_run_authorization",
        "generated_by": "scripts/reproduce_fgc_pro7_frz1.py",
        "scope_bindings": dict(raw["scope"]),
        "source_config_sha256": {
            path.relative_to(REPOSITORY).as_posix(): configuration["source_sha256"],
            protocol_relative: sha256(configuration["protocol_source"]).hexdigest(),
        },
        "predecessor_sha256": {
            raw["immutable_lineage"]["predecessor_protocol_config"]: raw["immutable_lineage"]["predecessor_protocol_sha256"],
            raw["immutable_lineage"]["diagnosis_config"]: raw["immutable_lineage"]["diagnosis_config_sha256"],
            raw["immutable_lineage"]["diagnosis_result"]: raw["immutable_lineage"]["diagnosis_result_sha256"],
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
                "accepted_state_source_failure_is_terminal": True,
                "all_source_only_unaccepted_proposal_failures_are_retryable": True,
                "candidate_endpoint_is_unaccepted_until_transaction_commit": True,
                "any_non_source_failure_vetoes_retry": True,
                "complete_proposal_is_previewed_before_classification": True,
                "retry_restarts_from_bitwise_last_accepted_state": True,
                "retry_factor": "1/2",
                "maximum_retry_count": 32,
                "minimum_step_size": "1/1073741824",
                "raw_source_residual_tolerance_is_unchanged": True,
                "source_refinement_algorithm_is_unchanged": True,
                "all_non_source_stops_are_unchanged": True,
                "terminal_and_retry_exhaustion_evidence_is_complete": True,
                "eligible_amplitudes_in_order": ["5/2", "3"],
            },
            "epistemic_boundary": {
                "PROTO6_GR0_calibration_trajectory_was_seen": True,
                "first_nonzero_common_event_completed": False,
                "PROTO7_trajectory_read": False,
                "FGCQR_outcome_read": False,
                "SGBL_outcome_read": False,
                "mechanism_question_answered": False,
                "physical_experiment_changed": False,
            },
            "decision": "PROTO7_is_frozen_but_HLT5_fresh_calibration_and_every_candidate_outcome_remain_closed",
        },
        "nonclaims": {
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
        },
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
        f"wrote {args.output} "
        "(PROTO7 frozen=true; fresh calibration authorized=false; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
