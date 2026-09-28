#!/usr/bin/env python3
"""Construct or verify CAL11's read-only PROTO14 terminal binder."""
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
sys.path[:0] = [str(REPOSITORY / "src")]
from recursive_horizons.fgc.evolution.cal11_proto14_campaign_diagnosis import (  # noqa: E402
    reconstruct_terminal, sha, source_pinned_static_replay,
)

ARTIFACT_ID = "FGC-1-CAL11-PREF22"
CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal11-pref22.toml"
OUTPUT = REPOSITORY / "results/fgc-1-cal11-pref22.json"
DOCUMENT = REPOSITORY / "docs/fgc-cal11-pref22.md"
IMPLEMENTATION = (
    Path(__file__).resolve(),
    REPOSITORY / "src/recursive_horizons/fgc/evolution/cal11_proto14_campaign_diagnosis.py",
)
EXPECTED_TRUE = {
    "PROTO14_campaign_terminal_result_bound",
    "PROTO14_terminal_checkpoint_recomputed",
    "PROTO14_terminal_checkpoint_owns_external_result",
    "PROTO14_terminal_event_identity_bound",
    "PROTO14_terminal_abort_rolled_back_to_event_23",
    "PROTO14_post_restart_accepted_macro_steps_zero",
    "PROTO14_invalid_runtime_or_nonconverged_terminal_observed",
}
EXPECTED_FALSE = {
    "PROTO14_runtime_fault_cause_fully_derived",
    "PROTO14_runtime_conformance_repaired",
    "PROTO14_terminal_checkpoint_may_resume",
    "historical_abort_is_scientific_GR0_stop",
    "historical_abort_is_temporal_admission_failure",
    "historical_abort_is_constraint_or_spatial_stop",
    "fresh_GR0_dynamic_calibration_completed",
    "GR0_case_eligible",
    "classical_spherical_diagnostic_authorized",
    "SGBL_execution_authorized",
    "FGCQR_holdout_execution_authorized",
    "FGCQR_mechanism_rejected",
    "DEF1_execution_authorized",
    "retained_EFT_evolution_authorized",
    "physical_transition_claim_authorized",
    "general_gradient_route_rejected",
    "singularity_resolution_derived",
    "child_domain_or_topology_derived",
    "dark_sector_mechanism_derived",
    "varying_locally_measured_c_derived",
}


def canonical(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("utf-8")


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args], cwd=REPOSITORY, check=check,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )


def _git_blob(commit: str, relative: str) -> bytes:
    return _git("show", f"{commit}:{relative}").stdout


def _load_canonical_bytes(payload: bytes, name: str) -> dict[str, Any]:
    value = json.loads(payload)
    if not isinstance(value, dict) or canonical(value) != payload:
        raise ValueError(f"{name} is not canonical JSON")
    return value


def validate_config_data(config: Mapping[str, Any]) -> None:
    expected_top = {
        "schema_version", "artifact_id", "project_version",
        "run_plan_config", "runtime_authorization_result",
        "protocol_freeze_result", "campaign_manifest", "campaign_event_log",
        "campaign_checkpoint", "campaign_result", "scope",
        "immutable_campaign", "terminal_identity", "authorization_lineage",
        "source_pinned_replay", "successor_boundary", "claims",
    }
    if set(config) != expected_top:
        raise ValueError("CAL11 configuration keys differ")
    expected_paths = {
        "run_plan_config": "configs/fgc/fgc-1-pro14-run1.toml",
        "runtime_authorization_result": "results/fgc-1-hlt12-mon12.json",
        "protocol_freeze_result": "results/fgc-1-pro14-frz1.json",
        "campaign_manifest": "runs/fgc-2-sf1/proto14/calibration/manifest.json",
        "campaign_event_log": "runs/fgc-2-sf1/proto14/calibration/events.jsonl",
        "campaign_checkpoint": "runs/fgc-2-sf1/proto14/calibration/latest-checkpoint.npz",
        "campaign_result": "runs/fgc-2-sf1/proto14/calibration/campaign-result.json",
    }
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != "0.11.0"
        or any(config.get(key) != value for key, value in expected_paths.items())
    ):
        raise ValueError("CAL11 top-level configuration contract differs")
    if config.get("scope") != {
        "target_protocol": "FGC-2-SF1-PROTO14",
        "role": "post_run_terminal_checkpoint_recomputed_invalid_runtime_result",
        "calibration_branch": "GR-0", "amplitude": "3",
        "GR0_numerical_trajectory_read": True,
        "SGBL_trajectory_read": False, "FGCQR_trajectory_read": False,
        "mechanism_question_answered": False,
        "physical_obstruction_inferred": False,
    }:
        raise ValueError("CAL11 scope differs")
    immutable = config.get("immutable_campaign", {})
    identity = config.get("terminal_identity", {})
    lineage = config.get("authorization_lineage", {})
    source_replay = config.get("source_pinned_replay", {})
    successor = config.get("successor_boundary", {})
    claims = config.get("claims", {})
    expected_immutable = {
        "authorization_commit": "7534a1662d8078a63c98025ecd464bb068fa012f",
        "campaign_id": "038c0128e76f48118b97d985384a1e7239edf7c97349cfc7552ee7d33e289e98",
        "manifest_sha256": "9396faf40b6b3681087c6336ffe1163048a1a45ef1de1f8acbd4025a1768e261",
        "event_log_sha256": "99e08c042820f582d4dc8c54beb00b96111aa4241d1977db5152ddca84fe0071",
        "checkpoint_sha256": "e7dc31a25a775a9377fbc45d1b239deb492bdc9ebbef85f94913513f885b015d",
        "campaign_result_sha256": "3d15623004a2044030b9039842ddbce04273949fbeb12b6b65cd61e3137bec69",
        "event_log_line_count": 2,
        "event_log_byte_count": 62410,
    }
    expected_lineage = {
        "run_plan_sha256": "76a8c07132b91cefc0cc466a6108557463841eca8ee1e4a081e718b9d9d472ba",
        "runtime_authorization_sha256": "7f4068a3116ff1838ce7f0165ae18e72d2eea68933b87f056b26d8149c935731",
        "protocol_freeze_result_sha256": "a0bac04e78fe9ed5af851694057e51a0de62e8682a2fbfb734e63c98537c305b",
        "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v14.toml",
        "protocol_config_sha256": "3043cbc23e7a78af764b6e74397512cb3f7bc78eace64485ab188963ba35f3a3",
        "runner": "scripts/run_fgc_gr0_calibration_v14.py",
        "runner_sha256": "4c13eaf86cf31f55e28b7f60661db78ad22e1bb5b7460adf7aeff455efdc8e45",
        "runner_id": "FGC-1-CAL11-RUN1-RUNNER",
        "runtime_module": "src/recursive_horizons/fgc/evolution/proto14_runtime.py",
        "runtime_module_sha256": "f5abcfcaf1d55c15fbb1b8825d348230f86c15f54dc50a2328533eb54b6ecca4",
        "tdg6_runtime_module": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
        "tdg6_runtime_module_sha256": "039de0bd008a2536bea6f7187afc3856bfc841d14790cbae5ca2fe27f7eb8e2d",
        "protocol_module": "src/recursive_horizons/fgc/evolution/protocol_v14.py",
        "protocol_module_sha256": "a4d766baf1cafa87f32f93594d4c771ac0ada3225007312546afc9325208d86c",
    }
    if (
        immutable != expected_immutable
        or identity != {
            "classification": "invalid_implementation_or_nonconverged_run",
            "reason": "invalid_implementation_or_nonconverged_run",
            "error_type": "ValueError",
            "error_message": "TDG6 subdivision is not bitwise uniform",
            "active_member": "RK4-2049", "target_common_event_index": 24,
            "last_fully_completed_common_event_index": 23,
            "last_fully_completed_coordinate_time": "23/16",
            "cross_member_state_rolled_back": True,
            "post_restart_accepted_macro_step_count_per_member": 0,
            "terminal_checkpoint_may_not_resume": True,
            "historical_abort_is_not_physical_classification": True,
        }
        or lineage != expected_lineage
        or source_replay != {
            "runtime_source": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
            "runtime_source_sha256": "039de0bd008a2536bea6f7187afc3856bfc841d14790cbae5ca2fe27f7eb8e2d",
            "required_function": "_subdivision_boundaries",
            "required_guard": "if not _same_binary64(right - left, step):",
            "required_exception": "raise ValueError(\"TDG6 subdivision is not bitwise uniform\")",
            "execution_scope": "source_pinned_static_pre_shadow_arithmetic_guard_only",
            "actual_arithmetic_invocation": False,
            "campaign_execution": False, "campaign_resume": False,
            "output_namespace_creation": False,
        }
        or successor != {
            "TDG7_subdivision_uniformity_diagnosis_design_authorized": True,
            "TDG7_runtime_repair_implemented": False,
            "PROTO14_runtime_mutation_authorized": False,
            "fresh_replacement_protocol_frozen": False,
            "fresh_GR0_calibration_authorized": False,
        }
        or set(claims) != EXPECTED_TRUE | EXPECTED_FALSE
        or any(claims.get(key) is not True for key in EXPECTED_TRUE)
        or any(claims.get(key) is not False for key in EXPECTED_FALSE)
    ):
        raise ValueError("CAL11 configuration contract differs")
    for section, names in (
        (immutable, ("manifest_sha256", "event_log_sha256", "checkpoint_sha256", "campaign_result_sha256")),
        (lineage, ("run_plan_sha256", "runtime_authorization_sha256", "protocol_freeze_result_sha256", "protocol_config_sha256", "runner_sha256", "runtime_module_sha256", "tdg6_runtime_module_sha256", "protocol_module_sha256")),
    ):
        if any(not isinstance(section.get(name), str) or len(section[name]) != 64 for name in names):
            raise ValueError("CAL11 hash contract differs")


def load_config(path: Path = CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _raw_paths(config: dict[str, object]) -> dict[str, Path]:
    return {key: REPOSITORY / str(config[key]) for key in (
        "campaign_manifest", "campaign_event_log", "campaign_checkpoint", "campaign_result"
    )}


def _raw_bundle_is_complete_or_absent(paths: Mapping[str, Path]) -> bool:
    present = [path.is_file() for path in paths.values()]
    if not any(present):
        return False
    if not all(present):
        raise ValueError("PROTO14 raw campaign bundle is partial")
    return True


def _authorization_lineage(
    config: Mapping[str, Any], manifest: Mapping[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, Mapping[str, Any]], bytes]:
    immutable = config["immutable_campaign"]
    lineage = config["authorization_lineage"]
    commit = immutable["authorization_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("CAL11 authorization commit is not an ancestor of HEAD")
    bindings = {
        config["run_plan_config"]: lineage["run_plan_sha256"],
        config["runtime_authorization_result"]: lineage["runtime_authorization_sha256"],
        config["protocol_freeze_result"]: lineage["protocol_freeze_result_sha256"],
        lineage["protocol_config"]: lineage["protocol_config_sha256"],
        lineage["runner"]: lineage["runner_sha256"],
        lineage["runtime_module"]: lineage["runtime_module_sha256"],
        lineage["tdg6_runtime_module"]: lineage["tdg6_runtime_module_sha256"],
        lineage["protocol_module"]: lineage["protocol_module_sha256"],
    }
    blob_hashes: dict[str, str] = {}
    blob_payloads: dict[str, bytes] = {}
    for relative, expected_hash in bindings.items():
        payload = _git_blob(commit, relative)
        digest = sha256(payload).hexdigest()
        if digest != expected_hash:
            raise ValueError(f"CAL11 authorization blob differs: {relative}")
        blob_hashes[relative] = digest
        blob_payloads[relative] = payload
    authorization = _load_canonical_bytes(
        blob_payloads[config["runtime_authorization_result"]],
        "authorization-era HLT12 result",
    )
    if authorization.get("artifact_id") != "FGC-1-HLT12-MON12":
        raise ValueError("CAL11 authorization artifact differs")
    members = authorization.get("artifact_payload", {}).get("restart_payloads")
    if not isinstance(members, list) or len(members) != 6:
        raise ValueError("CAL11 authorization restart payloads differ")
    expected_members = {item.get("key"): item for item in members if isinstance(item, dict)}
    if set(expected_members) != {
        "RK4-2049", "RK4-4097", "RK4-8193",
        "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
    }:
        raise ValueError("CAL11 authorization member set differs")
    if manifest is not None and (
        manifest.get("authorization_checkpoint_commit") != commit
        or manifest.get("authorization_sha256") != lineage["runtime_authorization_sha256"]
        or manifest.get("plan_sha256") != lineage["run_plan_sha256"]
    ):
        raise ValueError("CAL11 manifest authorization lineage differs")
    return ({
        "authorization_commit": commit,
        "authorization_commit_is_ancestor_of_HEAD": True,
        "authorization_blob_sha256": blob_hashes,
        "restart_member_count": 6,
    }, expected_members, blob_payloads[lineage["tdg6_runtime_module"]])


def build(config: dict[str, Any]) -> dict[str, object]:
    immutable = config["immutable_campaign"]
    paths = _raw_paths(config)
    raw_loaded = _raw_bundle_is_complete_or_absent(paths)
    manifest = None
    if raw_loaded:
        expected_hashes = {
            "campaign_manifest": "manifest_sha256", "campaign_event_log": "event_log_sha256",
            "campaign_checkpoint": "checkpoint_sha256", "campaign_result": "campaign_result_sha256",
        }
        for path_key, hash_key in expected_hashes.items():
            if sha(paths[path_key]) != immutable[hash_key]:
                raise ValueError(f"{path_key} hash differs")
        manifest = _load_canonical_bytes(paths["campaign_manifest"].read_bytes(), "PROTO14 manifest")
        _load_canonical_bytes(paths["campaign_result"].read_bytes(), "PROTO14 result")
    authorization, expected_members, historical_tdg6_source = _authorization_lineage(
        config, manifest,
    )
    member_summaries = {
        key: {
            "state_sha256": item["state_sha256"],
            "restart_payload_sha256": item["restart_payload_sha256"],
            "accepted_macro_step_count": 0,
            "temporal_retry_count": 0,
            "accumulated_debit_sha256": "a2f69ec892ac5ccef01ceb3ffe952383a3c96eae9c538e1263dee545c12e7af9",
        }
        for key, item in sorted(expected_members.items())
    }
    expected_raw = {
        "campaign_id": immutable["campaign_id"],
        "event_log_sha256": immutable["event_log_sha256"],
        "terminal_result_sha256": immutable["campaign_result_sha256"],
        "completed_common_event_index": 23,
        "aborted_target_common_event_index": 24,
        "all_members_zero_progress": True,
        "checkpoint_embeds_event_log_and_result": True,
        "member_summaries": member_summaries,
    }
    if raw_loaded:
        observed_raw = reconstruct_terminal(
            manifest_path=paths["campaign_manifest"], event_log_path=paths["campaign_event_log"],
            checkpoint_path=paths["campaign_checkpoint"], result_path=paths["campaign_result"],
            expected={
                **immutable,
                "terminal_identity": config["terminal_identity"],
                "authorization_lineage": config["authorization_lineage"],
            },
            expected_members=expected_members,
        )
        if observed_raw != expected_raw:
            raise ValueError("CAL11 raw reconstruction summary differs")
    source_replay = source_pinned_static_replay(
        source_payload=historical_tdg6_source,
        expected=config["source_pinned_replay"],
    )
    source_hashes = {path.relative_to(REPOSITORY).as_posix(): sha(path) for path in (CONFIG,)}
    implementation_hashes = {path.relative_to(REPOSITORY).as_posix(): sha(path) for path in IMPLEMENTATION}
    claims = config["claims"]
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "classification": "completed_PROTO14_invalid_runtime_terminal_result",
        "generated_by": "scripts/reproduce_fgc_cal11_pref22.py",
        "source_config_sha256": source_hashes,
        "implementation_sha256": implementation_hashes,
        "derivation_document": "docs/fgc-cal11-pref22.md",
        "derivation_document_sha256": sha(DOCUMENT),
        "artifact_payload": {
            "immutable_campaign": immutable,
            "authorization_lineage": authorization,
            "raw_bundle_presence_policy": {
                "all_or_none": True,
                "canonical_raw_bundle_was_bound": True,
                "portable_absence_skips_only_raw_reconstruction": True,
            },
            "terminal_checkpoint_replay": expected_raw,
            "source_pinned_static_replay": source_replay,
            "successor_boundary": config["successor_boundary"],
            "claim_boundary": claims,
        },
        "gate_status": claims,
        "nonclaims": {key: False for key in sorted(EXPECTED_FALSE)},
        "predecessor_sha256": {
            config["runtime_authorization_result"]: config["authorization_lineage"]["runtime_authorization_sha256"],
            config["protocol_freeze_result"]: config["authorization_lineage"]["protocol_freeze_result_sha256"],
            config["run_plan_config"]: config["authorization_lineage"]["run_plan_sha256"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    config = load_config()
    result = build(config)
    payload = canonical(result)
    output = args.output.resolve()
    if args.check:
        if not output.is_file() or output.read_bytes() != payload:
            raise ValueError("CAL11 canonical result differs")
    else:
        output.write_bytes(payload)
    print(json.dumps({
        "artifact_id": ARTIFACT_ID,
        "raw_bundle_loaded_for_this_verification": _raw_bundle_is_complete_or_absent(_raw_paths(config)),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
