#!/usr/bin/env python3
"""Reproduce HLT13's compact, synthetic-only PROTO15 runtime certificate.

The reproducer reads only tracked compact lineage and the runtime's deterministic
synthetic qualification summary.  It never opens a PROTO14 raw campaign bundle
or creates a PROTO15 output namespace.
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

from recursive_horizons.fgc.evolution import proto15_runtime  # noqa: E402


ARTIFACT_ID = "FGC-1-HLT13-MON13"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt13-mon13.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt13-mon13.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt13-mon13.md"
RUNTIME_PATH = REPOSITORY / "src/recursive_horizons/fgc/evolution/proto15_runtime.py"

COMPACT_LINEAGE = {
    "configs/fgc/fgc-1-pro15-frz1.toml": "1a45af49cca92701a79cf103fd6795ffebf6353af372bff779cd4aae90f71676",
    "configs/fgc/fgc-2-sf1-protocol-v15.toml": "f1f1849d85b6cbe1eb61267392a63ca11389f551df127033ddf905abb38fdc3b",
    "results/fgc-1-pro15-frz1.json": "99d24425cb5d40da2e248605265b83bcdc82a2126cfa4b5bce4239c6ef30635c",
    "scripts/reproduce_fgc_pro15_frz1.py": "0f06c51c4056edffd1a654834b9a9636c9fb34d46e1e15d9994ac0384eb2ca88",
    "src/recursive_horizons/fgc/evolution/protocol_v15.py": "c40936d2c8cb5eb7a8f36d8ba5c449b3ee797338ad01cc59e13df5474e2208c3",
    "configs/fgc/fgc-1-tdg7-imp3.toml": "1c989c6e698e40a3237931508b558fe1c8528316e35e8c29c8fc0e2e320c85e8",
    "results/fgc-1-tdg7-imp3.json": "364223f0f0f93c446b9d5b6cf09d3aa3cbe5b48bd1c175faf96a34459fd5b779",
}
TRUE_CLAIMS = {
    "PROTO15_frozen", "PROTO15_campaign_owned_cursor_contract_frozen",
    "PROTO15_six_restart_payloads_hash_bound", "PROTO14_invalid_runtime_history_preserved",
    "PROTO15_successor_runtime_implemented", "PROTO15_runtime_synthetic_qualification_passed",
}
FALSE_CLAIMS = {
    "PROTO15_pretrajectory_runtime_authorized", "PROTO15_fresh_GR0_dynamic_calibration_authorized",
    "PROTO15_trajectory_read", "PROTO15_output_namespace_created", "PROTO15_calibration_completed",
    "GR0_case_eligible", "classical_spherical_diagnostic_authorized", "SGBL_execution_authorized",
    "FGCQR_holdout_execution_authorized", "DEF1_execution_authorized",
    "retained_EFT_evolution_authorized", "physical_transition_claim_authorized",
    "FGCQR_mechanism_rejected", "general_gradient_route_rejected",
    "singularity_resolution_derived", "child_domain_or_topology_derived",
    "dark_sector_mechanism_derived", "varying_locally_measured_c_derived",
}
MEMBERS = ["RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"]
FAULT_HOOKS = ["before_write", "after_write", "before_flush", "after_flush", "before_fsync", "after_fsync", "before_replace", "after_replace", "before_dir_fsync", "after_dir_fsync"]
FAULT_SCENARIOS = ["state_object", "accepted_journal", "accepted_checkpoint", "rejection_journal", "cursor_transition", "terminal_record", "terminal_checkpoint", "rollback_record", "rollback_checkpoint"]
FAULT_PHASE_COVERAGE = sorted(FAULT_HOOKS)
FAULT_TARGETS = ["checkpoint", "journal_record", "state_object"]
ROLLBACK_JOURNAL_FIELD_LOCATIONS = {
    "new_journal_parent_sha256": "journal_envelope.journal_parent_sha256",
    "new_journal_record_sha256": "journal_envelope.record_sha256",
}
ACCEPTED_RECORD_FIELDS = [
    "member_key", "predecessor_cursor_chain_sha256",
    "successor_cursor_chain_sha256", "accepted_state_sha256",
    "TDG6_ledger_sha256", "executed_plan", "cursor",
]
EXECUTED_QUALIFICATION = {
    "qualification_case_count": 110,
    "nominal_case_count": 6,
    "mutation_case_count": 14,
    "fault_case_count": 90,
    "fault_scenario_count": 9,
    "fault_scenario_coverage": FAULT_SCENARIOS,
    "fault_phase_coverage": FAULT_PHASE_COVERAGE,
    "fault_target_coverage": FAULT_TARGETS,
    "ambiguous_replace_crash_case_count": 18,
    "fault_recovery_image_count": 108,
    "qualification_case_digest": "0acc88f1783f54f4c2b57fedd629f0a940945d1a01580054c7e35f099e8d8435",
    "all_qualification_cases_passed": True,
    "accepted_pending_retry_plan_bound": True,
    "rollback_journal_hashes_bound_by_record_envelope": True,
    "old_or_new_replace_crash_images_exercised": True,
}


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n"


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def _git_show_sha(commit: str, relative: str) -> str:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{relative}"], cwd=REPOSITORY,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if completed.returncode:
        raise ValueError(f"missing immutable HLT13 input: {relative}")
    return sha256(completed.stdout).hexdigest()


def _strict_mapping(name: str, actual: Mapping[str, Any], expected: Mapping[str, Any]) -> None:
    if dict(actual) != dict(expected):
        raise ValueError(f"HLT13 {name} differs")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    required = {
        "schema_version", "artifact_id", "project_version", "protocol_config",
        "protocol_freeze_config", "protocol_freeze_result", "runtime_module", "scope",
        "immutable_lineage", "runtime_binding", "synthetic_qualification", "namespace", "claims",
    }
    if set(config) != required:
        raise ValueError("HLT13 config root differs")
    if (config["schema_version"], config["artifact_id"], config["project_version"]) != (1, ARTIFACT_ID, PROJECT_VERSION):
        raise ValueError("HLT13 identity differs")
    if config["protocol_config"] != "configs/fgc/fgc-2-sf1-protocol-v15.toml" or config["protocol_freeze_config"] != "configs/fgc/fgc-1-pro15-frz1.toml" or config["protocol_freeze_result"] != "results/fgc-1-pro15-frz1.json" or config["runtime_module"] != _rel(RUNTIME_PATH):
        raise ValueError("HLT13 input paths differ")
    _strict_mapping("scope", config["scope"], {
        "target_protocol": "FGC-2-SF1-PROTO15", "runtime_owner": ARTIFACT_ID,
        "role": "synthetic_only_campaign_owned_durable_cursor_runtime_implementation_and_qualification",
        "PROTO15_trajectory_read": False, "historical_PROTO14_raw_bundle_read": False,
        "fresh_GR0_case_selected": False, "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False, "mechanism_question_answered": False,
    })
    lineage = config["immutable_lineage"]
    if lineage.get("checkpoint_commit") != "f3812d0b93e318890e973b4b21f681dbde8f437e" or lineage.get("checkpoint_must_be_ancestor_of_HEAD") is not True or lineage.get("tracked_compact_inputs_must_match_checkpoint_and_worktree") is not True:
        raise ValueError("HLT13 immutable checkpoint differs")
    expected_lineage_values = {
        "proto15_freeze_config_sha256": COMPACT_LINEAGE["configs/fgc/fgc-1-pro15-frz1.toml"],
        "proto15_protocol_config_sha256": COMPACT_LINEAGE["configs/fgc/fgc-2-sf1-protocol-v15.toml"],
        "proto15_freeze_result_sha256": COMPACT_LINEAGE["results/fgc-1-pro15-frz1.json"],
        "proto15_reproducer_sha256": COMPACT_LINEAGE["scripts/reproduce_fgc_pro15_frz1.py"],
        "proto15_protocol_core_sha256": COMPACT_LINEAGE["src/recursive_horizons/fgc/evolution/protocol_v15.py"],
        "tdg7_repair_config_sha256": COMPACT_LINEAGE["configs/fgc/fgc-1-tdg7-imp3.toml"],
        "tdg7_repair_result_sha256": COMPACT_LINEAGE["results/fgc-1-tdg7-imp3.json"],
    }
    if {key: lineage.get(key) for key in expected_lineage_values} != expected_lineage_values or set(lineage) != {"checkpoint_commit", "checkpoint_must_be_ancestor_of_HEAD", "tracked_compact_inputs_must_match_checkpoint_and_worktree", *expected_lineage_values}:
        raise ValueError("HLT13 compact lineage differs")
    binding = config["runtime_binding"]
    required_binding = {
        "runtime_module_sha256", "compact_synthetic_qualification_function", "content_addressed_state_objects",
        "immutable_one_record_journal_framing", "journal_record_envelope", "rational_and_lowercase_binary64_hex_time_identity",
        "persisted_attempt_and_generation_high_water_marks", "generational_checkpoint_selection_requires_complete_hash_chain",
        "same_generation_checkpoint_fork_is_invalid", "terminal_checkpoint_preserves_unchanged_nonterminal_cursor_set",
        "terminal_checkpoint_embeds_updated_terminal_member_ledger", "only_exact_post_rejection_suffix_is_reconstructible",
        "accepted_or_rollback_uncheckpointed_suffix_is_invalid", "every_accepted_fine_commit_publishes_complete_six_member_checkpoint",
        "campaign_rollback_issues_new_campaign_cursor_and_attempt_identities",
        "accepted_fine_commit_binds_exact_executed_tdg7_plan",
        "accepted_fine_commit_record_payload_fields",
        "rollback_protocol_journal_field_locations",
        "rollback_payload_omits_recursive_self_hash_by_design",
        "hash_integrity_requires_external_trusted_genesis",
        "wholesale_campaign_replacement_authenticated", "named_fault_hooks",
    }
    if set(binding) != required_binding or binding["compact_synthetic_qualification_function"] != "synthetic_qualification" or binding["journal_record_envelope"] != ["record_kind", "journal_parent_sha256", "payload", "record_sha256"]:
        raise ValueError("HLT13 runtime binding differs")
    if not isinstance(binding["runtime_module_sha256"], str) or len(binding["runtime_module_sha256"]) != 64:
        raise ValueError("HLT13 runtime hash differs")
    non_boolean_binding = {
        "runtime_module_sha256", "compact_synthetic_qualification_function",
        "journal_record_envelope", "accepted_fine_commit_record_payload_fields",
        "rollback_protocol_journal_field_locations", "named_fault_hooks",
        "wholesale_campaign_replacement_authenticated",
    }
    if any(binding[key] is not True for key in required_binding - non_boolean_binding):
        raise ValueError("HLT13 runtime invariant differs")
    if (binding["named_fault_hooks"] != FAULT_HOOKS
            or binding["accepted_fine_commit_record_payload_fields"] != ACCEPTED_RECORD_FIELDS
            or binding["rollback_protocol_journal_field_locations"] != ROLLBACK_JOURNAL_FIELD_LOCATIONS
            or binding["wholesale_campaign_replacement_authenticated"] is not False):
        raise ValueError("HLT13 fault hooks differ")
    synthetic = config["synthetic_qualification"]
    expected_synthetic = {
        "member_keys_in_canonical_order": MEMBERS, "no_historical_raw_state_or_checkpoint_may_be_opened": True,
        "deterministic_synthetic_state_objects_only": True,
        "fresh_retry_exhaustion_recovery_checkpoint_and_rollback_controls_required": True,
        "all_named_fault_hooks_must_be_exercised": True, "canonical_hash_bound_result_required": True,
        "stale_hash_partial_write_fork_and_cross_record_mutations_must_fail_closed": True,
        "wholesale_persistence_replacement_is_outside_integrity_scope": True,
        "external_trusted_genesis_is_required_for_authentication": True,
        **EXECUTED_QUALIFICATION,
    }
    if synthetic != expected_synthetic:
        raise ValueError("HLT13 synthetic qualification contract differs")
    if config["namespace"] != {"calibration_output_root": "runs/fgc-2-sf1/proto15/calibration", "holdout_output_root": "runs/fgc-2-sf1/proto15/holdout", "both_roots_must_remain_absent": True, "qualification_must_create_neither_root": True}:
        raise ValueError("HLT13 namespace boundary differs")
    claims = config["claims"]
    if set(claims) != TRUE_CLAIMS | FALSE_CLAIMS or any(claims[key] is not True for key in TRUE_CLAIMS) or any(claims[key] is not False for key in FALSE_CLAIMS):
        raise ValueError("HLT13 claim boundary differs")
    return config


def _verify_immutable_lineage(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    commit = config["immutable_lineage"]["checkpoint_commit"]
    ancestor = subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=REPOSITORY, check=False).returncode == 0
    if not ancestor:
        raise ValueError("HLT13 checkpoint is not an ancestor of HEAD")
    records = []
    for relative, digest in COMPACT_LINEAGE.items():
        live = REPOSITORY / relative
        if not live.is_file() or _sha(live) != digest or _git_show_sha(commit, relative) != digest:
            raise ValueError(f"HLT13 immutable compact lineage differs: {relative}")
        records.append({"path": relative, "sha256": digest, "matches_checkpoint_and_worktree": True})
    return records


def _runtime_summary(config: Mapping[str, Any]) -> dict[str, Any]:
    binding = config["runtime_binding"]
    if _sha(RUNTIME_PATH) != binding["runtime_module_sha256"]:
        raise ValueError("HLT13 runtime module hash differs")
    summary = getattr(proto15_runtime, binding["compact_synthetic_qualification_function"])()
    if not isinstance(summary, dict):
        raise ValueError("HLT13 synthetic runtime summary must be a mapping")
    required_true = {"synthetic_only", "content_addressed_state_objects", "immutable_one_record_journal_framing", "rational_and_hex_time_identity", "high_water_marks_persisted", "checkpoint_fork_rejected", "terminal_cursor_set_preserved", "terminal_ledger_embedded", "exact_rejection_suffix_recovery_only", "uncheckpointed_accepted_or_rollback_suffix_invalid", "six_member_fine_checkpoint", "accepted_pending_retry_plan_bound", "rollback_new_identities", "rollback_journal_hashes_bound_by_record_envelope", "all_fault_hooks_exercised", "old_or_new_replace_crash_images_exercised", "hash_integrity_requires_external_trusted_genesis", "all_qualification_cases_passed"}
    if summary.get("artifact_id") != ARTIFACT_ID or summary.get("six_member_contract") is not True or summary.get("six_member_order") != MEMBERS or any(summary.get(key) is not True for key in required_true):
        raise ValueError("HLT13 synthetic runtime qualification differs")
    if (summary.get("fault_hook_windows") != binding["named_fault_hooks"]
            or summary.get("rollback_protocol_journal_field_locations")
            != binding["rollback_protocol_journal_field_locations"]
            or summary.get("wholesale_campaign_replacement_authenticated") is not False
            or summary.get("physical_classification") is not False
            or any(summary.get(key) != value for key, value in EXECUTED_QUALIFICATION.items())):
        raise ValueError("HLT13 synthetic fault-hook qualification differs")
    if (summary["qualification_case_count"]
            != summary["nominal_case_count"] + summary["mutation_case_count"] + summary["fault_case_count"]
            or summary["fault_case_count"]
            != summary["fault_scenario_count"] * len(FAULT_HOOKS)
            or summary["ambiguous_replace_crash_case_count"]
            != summary["fault_scenario_count"] * 2
            or summary["fault_recovery_image_count"]
            != summary["fault_case_count"] + summary["ambiguous_replace_crash_case_count"]):
        raise ValueError("HLT13 executed qualification arithmetic differs")
    return summary


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    namespace_records = []
    for relative in (config["namespace"]["calibration_output_root"], config["namespace"]["holdout_output_root"]):
        if (REPOSITORY / relative).exists():
            raise ValueError(f"HLT13 namespace exists: {relative}")
        namespace_records.append({"path": relative, "absent_before_and_after_qualification": True, "created_by_qualification": False})
    lineage = _verify_immutable_lineage(config)
    summary = _runtime_summary(config)
    claims = dict(config["claims"])
    return {
        "schema_version": 1, "artifact_id": ARTIFACT_ID, "project_version": PROJECT_VERSION,
        "classification": "synthetic_only_campaign_owned_durable_cursor_runtime_implementation_and_qualification_without_trajectory_or_namespace",
        "generated_by": _rel(Path(__file__)), "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path), config["protocol_config"]: _sha(REPOSITORY / config["protocol_config"]), config["protocol_freeze_config"]: _sha(REPOSITORY / config["protocol_freeze_config"])},
        "implementation_sha256": {_rel(Path(__file__)): _sha(Path(__file__)), _rel(RUNTIME_PATH): _sha(RUNTIME_PATH)},
        "scope_bindings": dict(config["scope"]),
        "artifact_payload": {"immutable_lineage": {"evidence_checkpoint_commit": config["immutable_lineage"]["checkpoint_commit"], "checkpoint_is_ancestor_of_HEAD": True, "tracked_compact_records": lineage}, "runtime_binding": dict(config["runtime_binding"]), "synthetic_qualification_contract": dict(config["synthetic_qualification"]), "synthetic_qualification_summary": summary, "namespace_precondition": {"both_new_namespaces_absent": True, "qualification_created_no_namespace": True, "records": namespace_records}, "historical_boundary": {"historical_PROTO14_raw_bundle_opened": False, "PROTO14_invalid_runtime_history_preserved": True, "PROTO14_terminal_checkpoint_may_not_resume": True, "PROTO15_compact_lineage_only": True}, "decision": {"successor_runtime_owner": ARTIFACT_ID, "successor_runtime_implemented": True, "runtime_synthetic_qualification_passed": True, "pretrajectory_runtime_authorized": False, "calibration_authorized": False, "candidate_execution_authorized": False}},
        "gate_status": claims, "nonclaims": {key: False for key in FALSE_CLAIMS},
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_reject_duplicate_pairs)
    if not isinstance(value, dict) or value.get("artifact_id") != ARTIFACT_ID or source != _canonical(value):
        raise ValueError("HLT13 canonical result differs")
    return value


def verify_canonical(path: Path = DEFAULT_OUTPUT, config_path: Path = DEFAULT_CONFIG) -> None:
    if load_canonical_result(path) != record(config_path):
        raise ValueError("HLT13 result differs from fresh compact reproduction")


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload); handle.flush(); os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    if arguments.verify:
        verify_canonical(arguments.output, arguments.config); print(f"verified {arguments.output}")
    else:
        _atomic_write(arguments.output, _canonical(record(arguments.config)).encode("utf-8")); print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
