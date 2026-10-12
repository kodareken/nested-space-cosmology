#!/usr/bin/env python3
"""Reproduce the PROTO15 campaign-owned lifecycle freeze without opening a run.

This certificate deliberately reads only compact, tracked lineage artifacts.
It neither reads the historical PROTO14 raw campaign nor opens a restart
payload, checkpoint, run directory, or numerical trajectory.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.protocol_v15 import (  # noqa: E402
    PROTO15_CLAIMS,
    PROTO15_LIFECYCLE,
    PROTO15_RESTART_INPUTS,
    validate_sf1_protocol_v15,
)


ARTIFACT_ID = "FGC-1-PRO15-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro15-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro15-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro15-frz1.md"
IMPLEMENTATION = (Path(__file__).resolve(), REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v15.py")

EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO15",
    "predecessor_protocol": "FGC-2-SF1-PROTO14",
    "historical_invalid_terminal": "FGC-1-CAL11-PREF22",
    "stage_safe_repair": "FGC-1-TDG7-IMP3",
    "freeze_role": "premise_only_campaign_owned_temporal_lifecycle_contract_certificate",
    "calibration_branch": "GR-0",
    "calibration_amplitude": "3",
    "holdout_branch": "FGC-QR",
    "PROTO15_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "mechanism_question_answered": False,
}
EXPECTED_LINEAGE = {
    "evidence_checkpoint_commit": "64316d8f80f8532524c306817239f823a73b14be",
    "proto14_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v14.toml",
    "proto14_protocol_sha256": "3043cbc23e7a78af764b6e74397512cb3f7bc78eace64485ab188963ba35f3a3",
    "proto14_freeze_result": "results/fgc-1-pro14-frz1.json",
    "proto14_freeze_result_sha256": "a0bac04e78fe9ed5af851694057e51a0de62e8682a2fbfb734e63c98537c305b",
    "hlt12_runtime_authorization_result": "results/fgc-1-hlt12-mon12.json",
    "hlt12_runtime_authorization_result_sha256": "7f4068a3116ff1838ce7f0165ae18e72d2eea68933b87f056b26d8149c935731",
    "cal11_terminal_result": "results/fgc-1-cal11-pref22.json",
    "cal11_terminal_result_sha256": "1339792198ec78ad83bc3e11a2805f09f961607ab6c1f6f06e79bb88ba6b58b8",
    "tdg7_repair_config": "configs/fgc/fgc-1-tdg7-imp3.toml",
    "tdg7_repair_config_sha256": "1c989c6e698e40a3237931508b558fe1c8528316e35e8c29c8fc0e2e320c85e8",
    "tdg7_repair_result": "results/fgc-1-tdg7-imp3.json",
    "tdg7_repair_result_sha256": "364223f0f0f93c446b9d5b6cf09d3aa3cbe5b48bd1c175faf96a34459fd5b779",
}
EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto15/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto15/holdout",
    "both_roots_must_be_absent": True,
    "freeze_must_create_neither_root": True,
}
EXPECTED_PROTO15_INPUTS = {
    "protocol_config_sha256": "f1f1849d85b6cbe1eb61267392a63ca11389f551df127033ddf905abb38fdc3b",
    "protocol_core_sha256": "c40936d2c8cb5eb7a8f36d8ba5c449b3ee797338ad01cc59e13df5474e2208c3",
}
EXPECTED_PROOF_KEYS = {
    "sealed_checkpoint_must_be_ancestor_of_HEAD",
    "historical_tracked_inputs_must_match_commit_and_worktree",
    "freeze_consumes_compact_tracked_artifacts_only",
    "historical_raw_runs_history_and_checkpoints_must_not_be_opened",
    "six_restart_records_must_match_PROTO14_and_HLT12_compact_lineage",
    "historical_CAL11_invalid_runtime_boundary_must_remain_public",
    "TDG7_stage_safe_repair_must_be_synthetic_only_and_cursor_incomplete",
    "PROTO15_cursor_lifecycle_and_crash_order_must_validate_exactly",
    "fresh_and_retry_cursor_invariants_must_be_fail_closed",
    "future_runtime_checkpoint_and_common_event_rollback_contract_must_remain_public",
    "canonical_cursor_journal_and_TDG6_ledger_hash_rules_must_remain_public",
    "retry_exhaustion_terminal_record_and_checkpoint_contract_must_remain_public",
    "rollback_must_create_new_generation_not_reuse_old_cursor_bytes",
    "complete_six_member_checkpoint_and_rollback_contract_must_remain_public",
    "retry_successor_semantic_linkage_and_terminal_lock_must_remain_public",
    "recovery_authority_and_fault_injection_contract_must_remain_public",
    "new_output_namespaces_must_be_absent_and_uncreated",
    "adversarial_mutations_must_fail_closed",
    "canonical_duplicate_key_safe_result_required",
}
COMPACT_MEMBER_KEYS = (
    "key", "source_checkpoint", "method", "point_count", "coordinate_time",
    "state_sha256", "restart_payload_sha256", "input_hash", "step_index",
    "transaction_serial", "accepted_stage_count", "source_retry_count",
    "CFL_retry_count", "event_sample_count", "tracer_count", "state_shape",
)


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n"


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_pairs)
    if not isinstance(value, dict):
        raise ValueError(f"{_rel(path)} must contain an object")
    return value


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(["git", *args], cwd=REPOSITORY, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def validate_config_value(
    value: Mapping[str, Any], *, verify_live_hashes: bool
) -> dict[str, Any]:
    """Apply one exact config contract to both file loads and mutations."""

    if not isinstance(value, Mapping):
        raise ValueError("PRO15-FRZ1 config must be a mapping")
    value = dict(value)
    expected_keys = {"schema_version", "artifact_id", "project_version", "protocol_config", "protocol_core", "proto15_inputs", "scope", "immutable_lineage", "namespace_precondition", "proof_contract", "claims"}
    if set(value) != expected_keys:
        raise ValueError("PRO15-FRZ1 config root differs")
    if (value["schema_version"], value["artifact_id"], value["project_version"], value["protocol_config"], value["protocol_core"]) != (1, ARTIFACT_ID, PROJECT_VERSION, "configs/fgc/fgc-2-sf1-protocol-v15.toml", "src/recursive_horizons/fgc/evolution/protocol_v15.py"):
        raise ValueError("PRO15-FRZ1 config identity differs")
    if value["proto15_inputs"] != EXPECTED_PROTO15_INPUTS or value["scope"] != EXPECTED_SCOPE or value["immutable_lineage"] != EXPECTED_LINEAGE or value["namespace_precondition"] != EXPECTED_NAMESPACE or value["claims"] != PROTO15_CLAIMS:
        raise ValueError("PRO15-FRZ1 frozen bindings differ")
    if verify_live_hashes and (
        _sha(REPOSITORY / value["protocol_config"])
        != value["proto15_inputs"]["protocol_config_sha256"]
        or _sha(REPOSITORY / value["protocol_core"])
        != value["proto15_inputs"]["protocol_core_sha256"]
    ):
        raise ValueError("PROTO15 protocol/core bytes differ from freeze")
    proof = value["proof_contract"]
    if set(proof) != EXPECTED_PROOF_KEYS or any(flag is not True for flag in proof.values()):
        raise ValueError("PRO15-FRZ1 proof contract differs")
    return value


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    return validate_config_value(value, verify_live_hashes=True)


def load_protocol(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    validate_sf1_protocol_v15(value)
    return value


def _immutable_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["evidence_checkpoint_commit"]
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("sealed PRO15 evidence checkpoint is not an ancestor of HEAD")
    bindings = (
        ("proto14_protocol_config", "proto14_protocol_sha256"),
        ("proto14_freeze_result", "proto14_freeze_result_sha256"),
        ("hlt12_runtime_authorization_result", "hlt12_runtime_authorization_result_sha256"),
        ("cal11_terminal_result", "cal11_terminal_result_sha256"),
        ("tdg7_repair_config", "tdg7_repair_config_sha256"),
        ("tdg7_repair_result", "tdg7_repair_result_sha256"),
    )
    records = []
    for path_key, hash_key in bindings:
        relative, expected = lineage[path_key], lineage[hash_key]
        commit_bytes = _git("show", f"{commit}:{relative}").stdout
        live = REPOSITORY / relative
        if sha256(commit_bytes).hexdigest() != expected or _sha(live) != expected:
            raise ValueError(f"PRO15 immutable lineage differs: {relative}")
        records.append({"path": relative, "sha256": expected, "matches_sealed_commit_and_worktree": True})
    return {"evidence_checkpoint_commit": commit, "checkpoint_is_ancestor_of_HEAD": True, "tracked_compact_records": records}


def _restart_members(protocol: Mapping[str, Any], config: Mapping[str, Any]) -> list[dict[str, Any]]:
    lineage = config["immutable_lineage"]
    pro14 = _load_json(REPOSITORY / lineage["proto14_freeze_result"])
    hlt12 = _load_json(REPOSITORY / lineage["hlt12_runtime_authorization_result"])
    expected = protocol["restart_inputs"]["members"]
    pro14_members = pro14.get("artifact_payload", {}).get("restart_members")
    hlt12_members = hlt12.get("artifact_payload", {}).get("restart_payloads")
    if pro14.get("artifact_id") != "FGC-1-PRO14-FRZ1" or hlt12.get("artifact_id") != "FGC-1-HLT12-MON12":
        raise ValueError("PRO15 compact restart lineage identity differs")
    if not isinstance(pro14_members, list) or not isinstance(hlt12_members, list) or len(pro14_members) != 6 or len(hlt12_members) != 6:
        raise ValueError("PRO15 compact restart record count differs")
    observed: list[dict[str, Any]] = []
    for expected_member, older, authorized in zip(expected, pro14_members, hlt12_members, strict=True):
        compact = {key: older.get(key) for key in COMPACT_MEMBER_KEYS}
        if compact != {key: authorized.get(key) for key in COMPACT_MEMBER_KEYS}:
            raise ValueError(f"PRO15 compact PRO14/HLT12 restart disagreement: {expected_member['key']}")
        for field in ("key", "method", "point_count", "state_sha256", "restart_payload_sha256", "step_index", "transaction_serial"):
            if compact[field] != expected_member[field]:
                raise ValueError(f"PRO15 restart input differs: {expected_member['key']} {field}")
        if compact["source_checkpoint"] != expected_member["source"]:
            raise ValueError(f"PRO15 restart input differs: {expected_member['key']} source")
        ledger = authorized.get("TDG6_ledger_initialization")
        expected_ledger = {"accepted_macro_step_count": 0, "accumulated_debit_channel_count": 18, "all_debits_exact_zero": True, "cumulative_temporal_retry_count": 0, "serialized_temporal_rejection_count": 0}
        if compact["coordinate_time"] != 23 / 16 or ledger != expected_ledger:
            raise ValueError(f"PRO15 restart cursor origin differs: {expected_member['key']}")
        observed.append({**compact, "TDG6_ledger_initialization": ledger, "cursor_origin_template": {"not_a_serialized_complete_cursor": True, "mode": "FRESH_READY", "retry_successor_payload_or_none": None, "current_retry_count": 0, "accepted_boundary_time_rational": "23/16", "accepted_boundary_time_binary64_hex": float(23 / 16).hex(), "accepted_state_sha256": compact["state_sha256"], "method": compact["method"], "point_count": compact["point_count"], "cursor_generation": 0, "cursor_chain_parent_sha256": protocol_lifecycle_root_hash(), "journal_tip_sha256": protocol_lifecycle_root_hash()}})
    return observed


def _historical_boundary(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    cal11 = _load_json(REPOSITORY / lineage["cal11_terminal_result"])
    tdg7 = _load_json(REPOSITORY / lineage["tdg7_repair_result"])
    cal = cal11.get("gate_status", {})
    tdg = tdg7.get("gate_status", {})
    if (cal11.get("artifact_id") != "FGC-1-CAL11-PREF22" or cal.get("PROTO14_invalid_runtime_or_nonconverged_terminal_observed") is not True or cal.get("historical_abort_is_scientific_GR0_stop") is not False or cal.get("PROTO14_terminal_checkpoint_may_resume") is not False):
        raise ValueError("PRO15 CAL11 invalid-runtime boundary differs")
    if (tdg7.get("artifact_id") != "FGC-1-TDG7-IMP3" or tdg.get("TDG7_runtime_repair_synthetic_qualification_passed") is not True or tdg.get("fresh_vs_retry_protocol_discriminator_implemented") is not False or tdg.get("future_protocol_cursor_required_before_calibration") is not True or tdg.get("fresh_GR0_calibration_authorized") is not False):
        raise ValueError("PRO15 TDG7 limitation boundary differs")
    return {"PROTO14_invalid_runtime_or_nonconverged_not_physics": True, "PROTO14_terminal_checkpoint_may_not_resume": True, "TDG7_stage_safe_repair_synthetically_qualified": True, "TDG7_stateless_adapter_is_not_campaign_cursor": True, "campaign_owned_cursor_required_before_calibration": True, "historical_raw_campaign_not_opened": True}


def _namespace_evidence(config: Mapping[str, Any]) -> dict[str, Any]:
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = config["namespace_precondition"][key]
        if (REPOSITORY / relative).exists():
            raise ValueError(f"PRO15 prospective namespace already exists: {relative}")
        records.append({"path": relative, "absent_before_freeze": True, "created_by_freeze": False})
    return {"both_new_namespaces_absent": True, "freeze_created_no_namespace": True, "records": records}


def _mutation_controls(protocol: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, bool]:
    protocol_attacks = {
        "nonexhausted_retry_order": lambda v: v["lifecycle"]["nonexhausted_retry_transition_order"].reverse(),
        "exhaustion_order": lambda v: v["lifecycle"]["exhaustion_transition_order"].reverse(),
        "exhaustion_cannot_retry": lambda v: v["lifecycle"].__setitem__("TDG6_retry_exhaustion_must_not_create_RETRY_PENDING", False),
        "terminal_exhaustion_fields": lambda v: v["lifecycle"]["terminal_exhaustion_required_fields"].remove("journal_tip_sha256"),
        "journal_envelope": lambda v: v["lifecycle"]["journal_record_envelope_required_fields"].remove("journal_parent_sha256"),
        "journal_tip_semantics": lambda v: v["lifecycle"].__setitem__("cursor_transition_record_sha256_becomes_checkpoint_journal_tip", False),
        "ledger_hash_rule": lambda v: v["lifecycle"].__setitem__("TDG6_ledger_sha256_rule", "unspecified"),
        "journal_durability": lambda v: v["lifecycle"].__setitem__("durable_journal_append_requires_one_canonical_record_flush_and_fsync_before_return", False),
        "checkpoint_directory_fsync": lambda v: v["lifecycle"].__setitem__("atomic_checkpoint_requires_file_fsync_atomic_replace_and_parent_directory_fsync", False),
        "accepted_commit_generation": lambda v: v["lifecycle"].__setitem__("accepted_fine_state_TDG6_ledger_and_FRESH_READY_cursor_must_publish_as_one_atomic_checkpoint_generation", False),
        "pre_record_crash_restore": lambda v: v["lifecycle"].__setitem__("crash_before_new_durable_record_restores_exact_last_durable_cursor", False),
        "member_set": lambda v: v["lifecycle"]["member_keys_in_canonical_order"].pop(),
        "retry_semantics": lambda v: v["lifecycle"]["retry_successor_semantic_relations"].remove("successor_cap_is_exact_binary64_half_of_predecessor_cap"),
        "terminal_lock": lambda v: v["lifecycle"].__setitem__("terminal_exhaustion_requires_terminal_lock_true", False),
        "common_event_set": lambda v: v["lifecycle"].__setitem__("common_event_rollback_is_one_campaign_record_over_complete_six_member_snapshot", False),
        "recovery_order": lambda v: v["lifecycle"]["recovery_authority_order"].reverse(),
        "attempt_semantics": lambda v: v["lifecycle"].__setitem__("attempt_serial_strictly_increases_for_every_fresh_or_retry_proposal", False),
        "typed_outcome_owner": lambda v: v["lifecycle"].__setitem__("campaign_owner_must_complete_retry_or_terminal_persistence_before_typed_outcome_escapes", False),
        "rollback_generation": lambda v: v["lifecycle"].__setitem__("rollback_must_append_new_generation_binding_restored_snapshot_not_reuse_old_cursor_bytes", False),
        "fresh_payload": lambda v: v["lifecycle"].__setitem__("FRESH_READY_requires_retry_successor_payload_json_null", False),
        "fresh_retry_count": lambda v: v["lifecycle"].__setitem__("FRESH_READY_requires_TDG6_current_retry_count_zero", False),
        "retry_count": lambda v: v["lifecycle"].__setitem__("RETRY_PENDING_requires_TDG6_current_retry_count_positive", False),
        "crash_recovery": lambda v: v["lifecycle"].__setitem__("unmatched_durable_rejection_recovery_requires_exact_bound_evidence_and_typed_branch_reconstruction_or_invalid_stop", False),
        "cursor_reset": lambda v: v["lifecycle"].__setitem__("RETRY_PENDING_may_not_reset_to_FRESH_READY_without_accepted_fine_commit_or_explicit_common_event_rollback", False),
        "claim": lambda v: v["claims"].__setitem__("FGCQR_holdout_execution_authorized", True),
    }
    controls: dict[str, bool] = {}
    for name, mutate in protocol_attacks.items():
        attacked = deepcopy(protocol); mutate(attacked)
        try:
            validate_sf1_protocol_v15(attacked)
        except (TypeError, ValueError):
            controls[name] = True
        else:
            controls[name] = False
    config_attacks = {
        "restart_payload": lambda v: v["immutable_lineage"].__setitem__("proto14_freeze_result_sha256", "0" * 64),
        "namespace": lambda v: v["namespace_precondition"].__setitem__("calibration_output_root", "runs/fgc-2-sf1/proto14/calibration"),
        "authorization": lambda v: v["claims"].__setitem__("PROTO15_fresh_GR0_dynamic_calibration_authorized", True),
        "extra_root": lambda v: v.__setitem__("unexpected", True),
        "removed_proof_key": lambda v: v["proof_contract"].pop("adversarial_mutations_must_fail_closed"),
        "false_proof_flag": lambda v: v["proof_contract"].__setitem__("adversarial_mutations_must_fail_closed", False),
        "protocol_core_path": lambda v: v.__setitem__("protocol_core", "src/recursive_horizons/fgc/evolution/protocol_v14.py"),
    }
    for name, mutate in config_attacks.items():
        attacked = deepcopy(config); mutate(attacked)
        try:
            load_config_from_value(attacked)
        except ValueError:
            controls[name] = True
        else:
            controls[name] = False
    if not all(controls.values()):
        raise RuntimeError("PRO15 adversarial mutation control failed")
    return controls


def load_config_from_value(value: Mapping[str, Any]) -> None:
    validate_config_value(value, verify_live_hashes=False)


def protocol_lifecycle_root_hash() -> str:
    return PROTO15_LIFECYCLE["cursor_chain_root_parent_sha256"]


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    protocol_path = REPOSITORY / config["protocol_config"]
    protocol = load_protocol(protocol_path)
    protocol_validation = validate_sf1_protocol_v15(protocol)
    lineage = _immutable_lineage(config)
    restart_members = _restart_members(protocol, config)
    boundary = _historical_boundary(config)
    namespaces = _namespace_evidence(config)
    mutations = _mutation_controls(protocol, config)
    gates = dict(config["claims"])
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "prospective_campaign_owned_temporal_lifecycle_protocol_freeze_without_runtime_or_trajectory",
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path), _rel(protocol_path): _sha(protocol_path)},
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": dict(config["scope"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "protocol_validation": protocol_validation,
            "historical_boundary": boundary,
            "restart_members": restart_members,
            "lifecycle_contract": dict(protocol["lifecycle"]),
            "namespace_precondition": namespaces,
            "mutation_controls": mutations,
            "decision": {"PROTO15_frozen": True, "campaign_owned_cursor_contract_frozen": True, "successor_runtime_owner": protocol["inheritance"]["successor_runtime_owner"], "successor_runtime_implemented": False, "calibration_authorized": False, "candidate_execution_authorized": False},
            "epistemic_boundary": {"freeze_consumes_compact_tracked_artifacts_only": True, "raw_historical_campaign_or_checkpoint_opened": False, "freeze_advances_no_trajectory": True, "fresh_GR0_case_selected": False, "SGBL_or_FGCQR_outcome_read": False, "TDG7_synthetic_repair_is_not_a_GR0_calibration": True, "cursor_runtime_and_crash_recovery_are_future_implementation_obligations": True},
            "future_runtime_obligations": {"canonical_record_and_binary64_hex_identity": True, "persist_complete_TDG6_ledger_and_cursor_hash_chain": True, "bind_every_retry_to_exact_rejection_hash_and_semantic_successor_relations": True, "append_and_fsync_retry_transition_after_exact_durable_TDG6_rejection": True, "append_and_fsync_terminal_invalid_record_after_TDG6_exhaustion": True, "restore_exact_last_durable_cursor_before_any_new_record": True, "refuse_fresh_classification_after_unmatched_durable_rejection": True, "publish_accepted_state_ledger_and_FRESH_READY_cursor_as_one_atomic_checkpoint_generation": True, "bind_complete_six_member_cursor_ledger_state_sets_and_journal_tip_into_atomic_fsynced_checkpoints": True, "append_new_campaign_generation_common_event_rollback_over_complete_six_member_snapshot": True, "recover_by_highest_complete_checkpoint_then_contiguous_journal_suffix_only": True, "fault_inject_every_durability_boundary": True, "separately_qualify_runtime_before_namespace_creation": True},
        },
        "gate_status": gates,
        "nonclaims": {name: value for name, value in gates.items() if value is False},
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_reject_duplicate_pairs)
    if not isinstance(value, dict) or value.get("artifact_id") != ARTIFACT_ID or source != _canonical(value):
        raise ValueError("PRO15 canonical result differs")
    return value


def verify_canonical(path: Path = DEFAULT_OUTPUT, config_path: Path = DEFAULT_CONFIG) -> None:
    if load_canonical_result(path) != record(config_path):
        raise ValueError("PRO15 result differs from fresh compact reproduction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    if arguments.verify:
        verify_canonical(arguments.output, arguments.config)
        print(f"verified {arguments.output}")
    else:
        _atomic_write(arguments.output, _canonical(record(arguments.config)).encode("utf-8"))
        print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
