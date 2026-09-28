#!/usr/bin/env python3
"""Regenerate the FGC-1-PRO4-FRZ1 immutable PROTO4 freeze certificate."""

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
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro4-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro4-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro4-frz1.md"

from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    PROTO4_CLAIMS,
    SF1_PROTOCOL_ARTIFACT_ID,
    SF1_PROTOCOL_V4_ARTIFACT_ID,
    validate_sf1_protocol,
)


ARTIFACT_ID = "FGC-1-PRO4-FRZ1"
PROJECT_VERSION = "0.11.0"
EXPECTED_SCOPE = {
    "target_protocol": SF1_PROTOCOL_V4_ARTIFACT_ID,
    "predecessor_protocol": SF1_PROTOCOL_ARTIFACT_ID,
    "diagnosis_artifact": "FGC-1-CAL0-PREF2",
    "calibration_branch": "GR-0",
    "holdout_branch": "FGC-QR",
    "freeze_role": "premise_only_pre_holdout_protocol_and_lineage_certificate",
    "GR0_exploratory_outcomes_disclosed": True,
    "GR0_development_runs_may_enter_calibration_evidence": False,
    "FGCQR_evolution_outcomes_inspected": False,
    "SGBL_evolution_outcomes_inspected": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": "b5aa4934b58afd505b22e008202285c475942892",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "predecessor_protocol_sha256": "57b28eb9abe80c81f9e47ae6bea6ce2d69e103542df500e97c6c93805a8ac13b",
    "diagnosis_result": "results/fgc-1-cal0-pref2.json",
    "diagnosis_result_sha256": "d952e737d22d63b33e8d87e366538925bd28acbf82a3062563de49867a229138",
}
EXPECTED_PROOF = {
    "PROTO4_overlay_must_validate_exactly",
    "checkpoint_must_be_ancestor_of_HEAD",
    "predecessor_and_diagnosis_must_be_read_from_checkpoint",
    "immutable_blob_hashes_must_match_PROTO4_amendment",
    "immutable_PROTO3_must_validate_unchanged",
    "immutable_CAL0_must_be_canonical_and_premise_only",
    "action_couplings_profiles_amplitude_order_and_holdout_cases_must_be_inherited",
    "physical_equations_null_observable_outcomes_and_robustness_must_be_inherited",
    "every_legacy_stop_must_be_partitioned_exactly_once",
    "replacement_admission_must_be_convergence_and_error_budgeted",
    "new_output_namespaces_must_be_declared_but_not_accessed",
    "canonical_hash_bound_result_required",
}
EXPECTED_CLAIMS = {
    "PROTO4_resolved_holdout_manifest_authorized": False,
    "HLT2_successor_monitor_implemented": False,
    "classical_spherical_diagnostic_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}
IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/scoped_run_authorization.py",
    REPOSITORY / "scripts/reproduce_fgc_pro4_frz1.py",
)
INHERITED_PHYSICAL_SECTIONS = (
    "action",
    "dimensionless_parameters",
    "initial_data",
    "gauge",
    "numerics",
    "analysis",
    "robustness",
    "negative_scope",
    "classification",
    "revision_policy",
)


def _canonical(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def _semantic_hash(value: object) -> str:
    return sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _repo_path(name: str, value: object) -> Path:
    if not isinstance(value, str) or not value:
        raise TypeError(f"{name} must be a repository-relative path")
    candidate = (REPOSITORY / value).resolve()
    try:
        candidate.relative_to(REPOSITORY.resolve())
    except ValueError as exc:
        raise ValueError(f"{name} escapes the repository") from exc
    if not candidate.is_file():
        raise ValueError(f"{name} does not exist")
    return candidate


def _load_toml_bytes(source: bytes, name: str) -> dict[str, Any]:
    try:
        value = tomllib.loads(source.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ValueError(f"{name} is not valid UTF-8 TOML") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} root must be a table")
    return value


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_unique_pairs)
    except json.JSONDecodeError as exc:
        raise ValueError("PRO4-FRZ1 result is not valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("PRO4-FRZ1 result is not canonical sorted JSON")
    return value


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        raw = tomllib.loads(source)
    except tomllib.TOMLDecodeError as exc:
        raise ValueError("PRO4-FRZ1 config is not valid TOML") from exc
    _keys(
        "PRO4-FRZ1 root",
        raw,
        {
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
        },
    )
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("PRO4-FRZ1 identity or convention differs")
    if not isinstance(raw["scope"], Mapping) or dict(raw["scope"]) != EXPECTED_SCOPE:
        raise ValueError("PRO4-FRZ1 scope differs")
    if (
        not isinstance(raw["immutable_lineage"], Mapping)
        or dict(raw["immutable_lineage"]) != EXPECTED_LINEAGE
    ):
        raise ValueError("PRO4-FRZ1 immutable lineage differs")
    proof = raw["proof_contract"]
    if not isinstance(proof, Mapping):
        raise TypeError("proof_contract must be a table")
    _keys("proof_contract", proof, EXPECTED_PROOF)
    if any(value is not True for value in proof.values()):
        raise ValueError("every PRO4-FRZ1 proof-contract premise must remain true")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("PRO4-FRZ1 claims must remain fail-closed")

    protocol_path = _repo_path("protocol_config", raw["protocol_config"])
    protocol_source = protocol_path.read_bytes()
    protocol_raw = _load_toml_bytes(protocol_source, "PROTO4")
    protocol = validate_sf1_protocol(protocol_raw)
    if (
        protocol.get("artifact_id") != SF1_PROTOCOL_V4_ARTIFACT_ID
        or protocol.get("protocol_version") != 4
        or protocol.get("frozen") is not True
    ):
        raise ValueError("PRO4-FRZ1 did not load frozen PROTO4")
    amendment = protocol_raw.get("amendment")
    if not isinstance(amendment, Mapping):
        raise ValueError("PROTO4 amendment is absent")
    amendment_lineage = {
        "checkpoint_commit": amendment.get("diagnosis_checkpoint_commit"),
        "predecessor_protocol_config": amendment.get("predecessor_protocol_config"),
        "predecessor_protocol_sha256": amendment.get("predecessor_protocol_sha256"),
        "diagnosis_result": amendment.get("diagnosis_result"),
        "diagnosis_result_sha256": amendment.get("diagnosis_result_sha256"),
    }
    if amendment_lineage != EXPECTED_LINEAGE:
        raise ValueError("PROTO4 amendment and PRO4-FRZ1 lineage differ")
    if dict(protocol_raw.get("claims", {})) != PROTO4_CLAIMS:
        raise ValueError("PROTO4 protocol claims are not fail-closed")
    return {
        "raw": raw,
        "source_sha256": sha256(source.encode("utf-8")).hexdigest(),
        "protocol_path": protocol_path,
        "protocol_source": protocol_source,
        "protocol_raw": protocol_raw,
        "protocol_certificate": protocol,
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
        raise ValueError("PROTO4 diagnosis checkpoint is not an ancestor of HEAD")

    predecessor_blob = _git_blob(
        commit,
        lineage["predecessor_protocol_config"],
        "PROTO3 protocol",
    )
    predecessor_digest = sha256(predecessor_blob).hexdigest()
    if predecessor_digest != lineage["predecessor_protocol_sha256"]:
        raise ValueError("immutable PROTO3 blob hash differs")
    predecessor_raw = _load_toml_bytes(predecessor_blob, "immutable PROTO3")
    predecessor = validate_sf1_protocol(predecessor_raw)
    if (
        predecessor.get("artifact_id") != SF1_PROTOCOL_ARTIFACT_ID
        or predecessor.get("protocol_version") != 3
        or predecessor.get("frozen") is not True
    ):
        raise ValueError("immutable predecessor is not sealed PROTO3")

    diagnosis_blob = _git_blob(
        commit,
        lineage["diagnosis_result"],
        "CAL0 diagnosis",
    )
    diagnosis_digest = sha256(diagnosis_blob).hexdigest()
    if diagnosis_digest != lineage["diagnosis_result_sha256"]:
        raise ValueError("immutable CAL0 blob hash differs")
    try:
        diagnosis = json.loads(
            diagnosis_blob.decode("utf-8"),
            object_pairs_hook=_unique_pairs,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("immutable CAL0 is not valid unique-key JSON") from exc
    if not isinstance(diagnosis, dict) or diagnosis_blob.decode("utf-8") != _canonical(
        diagnosis
    ):
        raise ValueError("immutable CAL0 is not canonical sorted JSON")
    gates = diagnosis.get("gate_status")
    payload = diagnosis.get("artifact_payload")
    boundary = payload.get("epistemic_boundary") if isinstance(payload, Mapping) else None
    if (
        diagnosis.get("artifact_id") != "FGC-1-CAL0-PREF2"
        or diagnosis.get("classification")
        != "pre_holdout_protocol_numerical_contract_obstruction"
        or not isinstance(gates, Mapping)
        or gates.get("PROTO3_pre_holdout_numerical_contract_obstruction_verified")
        is not True
        or gates.get("PROTO3_resolved_holdout_manifest_authorized") is not False
        or gates.get("PROTO4_premise_revision_required") is not True
        or not isinstance(boundary, Mapping)
        or boundary.get("FGCQR_evolution_outcome_inspected") is not False
        or boundary.get("SGBL_evolution_outcome_inspected") is not False
        or boundary.get("FGCQR_action_or_mechanism_tested") is not False
        or boundary.get("general_gradient_mechanism_rejected") is not False
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable CAL0 is not the premise-only diagnosis")

    current_predecessor = REPOSITORY / lineage["predecessor_protocol_config"]
    if sha256(current_predecessor.read_bytes()).hexdigest() != predecessor_digest:
        raise ValueError("current PROTO3 path differs from its immutable checkpoint blob")

    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "predecessor_protocol": {
            "artifact_id": predecessor["artifact_id"],
            "protocol_version": predecessor["protocol_version"],
            "blob_sha256": predecessor_digest,
            "semantic_holdout_contract_sha256": predecessor[
                "semantic_holdout_contract_sha256"
            ],
        },
        "diagnosis": {
            "artifact_id": diagnosis["artifact_id"],
            "classification": diagnosis["classification"],
            "blob_sha256": diagnosis_digest,
            "PROTO3_obstruction_verified": True,
            "premise_revision_required": True,
            "FGCQR_outcome_inspected": False,
            "SGBL_outcome_inspected": False,
        },
        "predecessor_raw": predecessor_raw,
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(path)
    lineage = _immutable_lineage(configuration)
    predecessor_raw = lineage.pop("predecessor_raw")
    protocol = configuration["protocol_certificate"]
    protocol_raw = configuration["protocol_raw"]
    source_relative = str(path.resolve().relative_to(REPOSITORY.resolve()))
    protocol_relative = str(
        configuration["protocol_path"].resolve().relative_to(REPOSITORY.resolve())
    )

    inherited_hashes = {
        name: _semantic_hash(predecessor_raw[name])
        for name in INHERITED_PHYSICAL_SECTIONS
    }
    immutable_cases = {
        "amplitude_candidates": predecessor_raw["initial_data"]["chi"][
            "amplitude_candidates"
        ],
        "held_out_case_ids": predecessor_raw["partition"]["held_out_case_ids"],
        "held_out_phi_seed_factors": predecessor_raw["partition"]
        ["held_out_phi_seed_factors"],
        "held_out_chi_width_factors": predecessor_raw["partition"]
        ["held_out_chi_width_factors"],
    }
    if (
        immutable_cases["amplitude_candidates"] != protocol["amplitude_candidates"]
        or immutable_cases["held_out_case_ids"] != protocol["held_out_case_ids"]
    ):
        raise ValueError("PROTO4 declared cases differ from immutable PROTO3")

    applicability = protocol_raw["stop_applicability"]
    admission = protocol_raw["numerical_admission"]
    replacement = protocol_raw["replacement"]
    artifact_payload = {
        "protocol_certificate": protocol,
        "immutable_lineage": lineage,
        "inherited_physical_contract": {
            "section_semantic_sha256": inherited_hashes,
            "immutable_case_coordinates": immutable_cases,
            "action_couplings_profiles_and_equations_changed": False,
            "null_observable_outcome_rules_and_robustness_changed": False,
            "inheritance_is_by_exact_PROTO3_blob_hash": True,
        },
        "numerical_contract_revision": {
            "legacy_stop_count": protocol["legacy_stop_count"],
            "universal_runtime_stop_count": len(
                applicability["universal_runtime_stop_ids"]
            ),
            "candidate_branch_only_stop_count": len(
                applicability["candidate_branch_only_stop_ids"]
            ),
            "replaced_stop_count": len(
                applicability["replaced_by_PROTO4_convergent_admission_ids"]
            ),
            "partition_is_exact_and_disjoint": True,
            "spectral_admission_semantic_sha256": _semantic_hash(
                admission["spectrum"]
            ),
            "constraint_admission_semantic_sha256": _semantic_hash(
                admission["constraints"]
            ),
            "calibration_admission_semantic_sha256": _semantic_hash(
                admission["calibration"]
            ),
            "last_bin_boolean_is_diagnostic_only": admission["spectrum"]
            ["last_bin_alias_boolean_is_diagnostic_only"],
            "three_nested_common_event_admission_required": replacement[
                "calibration"
            ]["three_nested_common_event_admission_required"],
            "Richardson_error_enters_conservative_observable_error_sum": admission[
                "constraints"
            ]["Richardson_error_enters_conservative_observable_error_sum"],
        },
        "namespace_contract": {
            "resolved_holdout_manifest": protocol["resolved_holdout_manifest"],
            "calibration_output_root": protocol["calibration_output_root"],
            "holdout_output_root": protocol["holdout_output_root"],
            "output_namespace_contents_accessed_by_this_certificate": False,
            "emptiness_must_be_proved_by_future_manifest_before_first_run": True,
        },
        "epistemic_boundary": {
            "GR0_exploratory_outcomes_were_seen_before_freeze": True,
            "GR0_exploratory_runs_admissible_as_calibration_evidence": False,
            "fresh_primary_and_comparator_calibration_runs_required": True,
            "FGCQR_evolution_outcome_inspected": False,
            "SGBL_evolution_outcome_inspected": False,
            "HLT2_monitor_implemented": False,
            "resolved_holdout_manifest_constructed": False,
            "collapse_or_trapped_interval_derived": False,
            "regulator_activation_or_defocusing_tested": False,
            "FGCQR_action_or_mechanism_tested": False,
            "FGCQR_action_or_mechanism_rejected": False,
            "general_gradient_route_rejected": False,
        },
    }
    if (
        artifact_payload["numerical_contract_revision"][
            "universal_runtime_stop_count"
        ]
        + artifact_payload["numerical_contract_revision"][
            "candidate_branch_only_stop_count"
        ]
        + artifact_payload["numerical_contract_revision"]["replaced_stop_count"]
        != 26
    ):
        raise ValueError("PROTO4 stop partition does not cover 26 legacy stops")

    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "outcome_neutral_successor_protocol_freeze_certificate",
        "generated_by": "scripts/reproduce_fgc_pro4_frz1.py",
        "derivation_document": str(OWNER_DOCUMENT.relative_to(REPOSITORY)),
        "derivation_document_sha256": sha256(OWNER_DOCUMENT.read_bytes()).hexdigest(),
        "source_config_sha256": {
            source_relative: configuration["source_sha256"],
            protocol_relative: sha256(configuration["protocol_source"]).hexdigest(),
        },
        "predecessor_sha256": {
            EXPECTED_LINEAGE["predecessor_protocol_config"]: EXPECTED_LINEAGE[
                "predecessor_protocol_sha256"
            ],
        },
        "implementation_sha256": {
            str(item.relative_to(REPOSITORY)): sha256(item.read_bytes()).hexdigest()
            for item in IMPLEMENTATION
        },
        "certificate_contract": {
            "artifact_specific_payload_validated": True,
            "canonical_reproduction_passed": True,
            "scope_bindings_verified": True,
            "outcome_data_not_used_to_select_certificate_contract": True,
        },
        "artifact_payload": artifact_payload,
        "gate_status": {
            "PROTO4_outcome_neutral_protocol_frozen": True,
            "PROTO4_immutable_lineage_verified": True,
            "PROTO4_branch_specific_stop_partition_verified": True,
            "PROTO4_convergent_admission_contract_frozen": True,
            "PROTO4_resolved_holdout_manifest_authorized": False,
            "HLT2_successor_monitor_implemented": False,
            "classical_spherical_diagnostic_authorized": False,
            "FGCQR_holdout_execution_authorized": False,
            "retained_EFT_evolution_authorized": False,
            "physical_transition_claim_authorized": False,
        },
        "nonclaims": dict(EXPECTED_CLAIMS),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    payload = record(arguments.config)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(_canonical(payload), encoding="utf-8")
    print(
        "wrote "
        f"{arguments.output} "
        f"(PROTO4 frozen={payload['gate_status']['PROTO4_outcome_neutral_protocol_frozen']}, "
        "FGCQR holdout authorized="
        f"{payload['gate_status']['FGCQR_holdout_execution_authorized']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
