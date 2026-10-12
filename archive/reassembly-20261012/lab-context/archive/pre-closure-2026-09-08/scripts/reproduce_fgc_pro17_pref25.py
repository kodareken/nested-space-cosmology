#!/usr/bin/env python3
"""Reproduce the independent compact FGC-1-PRO17-PREF25 binder."""

from __future__ import annotations

import argparse
import ast
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


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from recursive_horizons.fgc.evolution import proto17_construction_theorem as theorem  # noqa: E402


ARTIFACT = "FGC-1-PRO17-PREF25"
SEALED_COMMIT = "f4d6e766aeb6351dc34ea482ddb3bd8eaac8d96d"
CONFIG = ROOT / "configs/fgc/fgc-1-pro17-pref25.toml"
OUTPUT = ROOT / "results/fgc-1-pro17-pref25.json"
DOC = ROOT / "docs/fgc-pro17-pref25.md"

SEALED = {
    "configs/fgc/fgc-2-sf1-protocol-v17.toml": "520466fb2f798159f7d0183938327eec0fa1ccf8b49c41101924f1ff01fcd484",
    "configs/fgc/fgc-1-pro17-frz1.toml": "6ff3cfad5eebb7824434e3b7f97da3ce7a8610312523e3461c8999532678d2ec",
    "results/fgc-1-pro17-frz1.json": "8e4b3f77015374cf825c141dff34d343eb3d5f4ccebfed5082bc283c1d414fea",
    "scripts/reproduce_fgc_pro17_frz1.py": "1481d8da60d2684a692fe9526236a0b82e6d63299e922545c9b9ee75b9c1a3cf",
    "src/recursive_horizons/fgc/evolution/protocol_v17.py": "2933d60d90607e0361d7deaec5395d8ad70a15ec41ab98ba2dd6aedbe9c5085c",
    "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py": "5fbd7f990d9f821344bbaae80840070386b3ae483fc278a85ae2a3ee5de8fa1d",
    "docs/fgc-pro17-frz1.md": "9b24df6c00860a5f00fe777d5efe32971b007493aef71cf8f27fc29a1b6efe20",
    "tests/test_fgc_protocol_v17.py": "e3e9c769f34768cd235a7c36b73f9a0df4c66abdf431fa5a04d43026e0b1c975",
    "tests/test_fgc_proto17_pure_construction.py": "76d628b352cfd35c362a04a686ef51f4685ddbd551945ff04b78e2471cdcf2ab",
    "tests/test_fgc_pro17_frz1_reproduction.py": "44dfee4557784d16e99e044f381fcad5f433b4d5f4a7aa2e1d45d41ac0638061",
}

GENESIS_FIELDS = (
    "genesis_algorithm_id", "canonical_json_encoding", "state_object_schema_id",
    "tdg6_ledger_schema_id", "cursor_schema_id", "checkpoint_schema_id",
    "protocol_artifact_id", "protocol_config_sha256", "protocol_freeze_result_sha256",
    "hlt13_result_sha256", "run_plan_sha256", "runtime_module_sha256",
    "adapter_module_sha256", "runner_sha256", "numerical_environment", "campaign_id",
    "namespace", "journal_root_parent_sha256", "common_event_index", "restart_coordinate_time",
    "active_event_target_time", "member_keys_in_canonical_order", "physical_input_manifest",
    "physical_input_manifest_sha256", "member_descriptors", "member_descriptor_set_sha256",
    "genesis_checkpoint_sha256", "cursor_set_sha256", "ledger_set_sha256",
    "accepted_state_set_sha256",
)
MEMBER_DESCRIPTOR_FIELDS = (
    "member_key", "method", "point_count", "source_checkpoint", "source_bundle_key",
    "accepted_boundary_time", "step_index", "transaction_serial", "state_object",
    "state_object_canonical_sha256", "initial_TDG6_ledger", "initial_TDG6_ledger_sha256",
)
STATE_FIELDS = (
    "schema_id", "member_key", "method", "point_count", "coordinate_time",
    "accepted_boundary_time", "step_index", "transaction_serial", "input_hash",
    "source_retry_count", "CFL_retry_count", "runtime_monitor_state", "causal_state",
    "tracer_state", "event_history_state", "physical_arrays", "physical_state_sha256",
    "restart_payload_sha256",
)
ARRAY_FIELDS = (
    "logical_name", "source_bundle_key", "npz_storage_key", "dtype", "layout",
    "shape", "little_endian_c_bytes_sha256",
)
LEDGER_FIELDS = (
    "last_accepted_time_hex", "accepted_macro_step_count",
    "cumulative_temporal_retry_count", "current_macro_step_temporal_retry_count",
    "last_accepted_macro_step_temporal_retry_count", "accumulated_debit_vector_hex",
    "serialized_temporal_rejections",
)
MEMBERS = ["RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"]
CONTROL_NAMES = [
    "array_storage_key", "foreign_protocol", "frozen_restart", "input_hash",
    "manifest_digest", "nonzero_generation_root_journal", "physical_state_digest",
    "receipt_cycle_injection", "receipt_sequence", "state_object_digest",
    "successor_checkpoint_hash", "successor_high_water", "successor_untouched_method",
    "zero_ledger",
]
FUTURE = {
    "raw_bundle_byte_and_hash_recomputation": {"required": True, "completed": False},
    "external_Git_authority_blob_and_live_import_validation": {"required": True, "completed": False},
    "semantic_historical_journal_payload_replay": {"required": True, "completed": False},
    "namespace_reuse_and_foreign_store_rejection": {"required": True, "completed": False},
}
TRUE = {
    "PROTO17_immutable_freeze_independently_bound",
    "PROTO17_independent_schema_and_construction_theorem_bound",
    "PROTO17_independent_generation_zero_abstract_fixture_bound",
    "PROTO17_independent_common_event_abstract_fixture_bound",
    "PROTO17_all_listed_synthetic_mutations_rejected",
    "PROTO17_future_HLT14_requirements_preserved_uncompleted",
    "HLT14_implementation_and_synthetic_qualification_authorized",
}


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canon(value: Any) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n"


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def _load_config() -> dict[str, Any]:
    with CONFIG.open("rb") as handle:
        value = tomllib.load(handle)
    if not isinstance(value, dict):
        raise ValueError("PREF25 configuration is not object-valued")
    _validate_config(value)
    return value


def _git_show(path: str) -> bytes:
    return subprocess.run(
        ["git", "show", f"{SEALED_COMMIT}:{path}"], cwd=ROOT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
    ).stdout


def _validate_commit(commit: object) -> str:
    if commit != SEALED_COMMIT:
        raise ValueError("PREF25 sealed commit differs")
    if subprocess.run(["git", "merge-base", "--is-ancestor", SEALED_COMMIT, "HEAD"], cwd=ROOT).returncode:
        raise ValueError("sealed PROTO17 commit is not an ancestor of HEAD")
    return SEALED_COMMIT


def _validate_sealed_blob(path: str, expected: object, historical: bytes, live: Path) -> None:
    if expected != SEALED[path] or sha256(historical).hexdigest() != expected or not live.is_file() or _sha(live) != expected:
        raise ValueError(f"sealed PROTO17 byte drift: {path}")


def _validate_config(value: Mapping[str, Any]) -> None:
    expected_scope = {
        "target_protocol": "FGC-2-SF1-PROTO17",
        "role": "post_freeze_independent_abstract_schema_genesis_and_common_event_construction_theorem_binder",
        "sealed_compact_PROTO17_blobs_read_from_immutable_commit": True,
        "raw_history_or_physical_bundle_opened": False,
        "campaign_checkpoint_loaded": False,
        "runtime_or_runner_implemented": False,
        "namespace_created": False,
        "trajectory_read": False,
        "synthetic_structural_fixtures_only": True,
        "mechanism_question_answered": False,
    }
    expected_lineage = {
        "checkpoint_commit": SEALED_COMMIT,
        "checkpoint_must_be_ancestor_of_HEAD": True,
        "all_bound_PROTO17_blobs_must_match_checkpoint_commit_and_worktree": True,
        "protocol_config_sha256": SEALED["configs/fgc/fgc-2-sf1-protocol-v17.toml"],
        "protocol_freeze_config_sha256": SEALED["configs/fgc/fgc-1-pro17-frz1.toml"],
        "protocol_freeze_result_sha256": SEALED["results/fgc-1-pro17-frz1.json"],
        "protocol_freeze_reproducer_sha256": SEALED["scripts/reproduce_fgc_pro17_frz1.py"],
        "protocol_schema_module_sha256": SEALED["src/recursive_horizons/fgc/evolution/protocol_v17.py"],
        "protocol_construction_module_sha256": SEALED["src/recursive_horizons/fgc/evolution/proto17_pure_construction.py"],
        "protocol_freeze_document_sha256": SEALED["docs/fgc-pro17-frz1.md"],
        "protocol_schema_test_sha256": SEALED["tests/test_fgc_protocol_v17.py"],
        "protocol_construction_test_sha256": SEALED["tests/test_fgc_proto17_pure_construction.py"],
        "protocol_freeze_test_sha256": SEALED["tests/test_fgc_pro17_frz1_reproduction.py"],
        "sealed_result_must_be_duplicate_safe_canonical_JSON": True,
        "sealed_result_namespace_absence_is_historical_not_reobserved": True,
        "trusted_Git_history_remains_an_explicit_external_assumption": True,
    }
    theorem_keys = {
        "stdlib_only", "must_not_import_protocol_v17_or_proto17_pure_construction",
        "must_parse_sealed_TOML_bytes_independently", "must_bind_all_five_exact_schema_field_tuples",
        "must_bind_canonical_six_member_order", "must_construct_complete_visible_synthetic_GenesisSpec",
        "must_construct_generation_zero_checkpoint_and_store_plan",
        "must_construct_structural_six_record_accepted_fine_prefix",
        "must_construct_predecessor_only_COMMON_EVENT_COMMIT",
        "must_construct_exact_successor_checkpoint", "must_recover_only_one_exact_lone_receipt_suffix",
        "must_match_sealed_genesis_receipt_and_successor_hashes",
        "historical_record_payload_semantics_are_not_replayed",
        "raw_physical_array_bytes_are_not_opened_or_recomputed",
    }
    qualified = {
        "sealed_control_count": 14, "sealed_control_names": CONTROL_NAMES,
        "independent_equivalent_control_corpus_required": True,
        "duplicate_safe_canonical_JSON_control_required": True,
        "sealed_blob_and_commit_mutation_controls_required": True,
        "schema_tuple_and_member_order_mutation_controls_required": True,
        "deferred_requirement_promotion_control_required": True,
        "runtime_claim_promotion_control_required": True,
    }
    authorization = {
        "authorized_object": "implementation_and_synthetic_qualification_of_FGC-1-HLT14-MON14_against_the_sealed_PROTO17_construction",
        "HLT14_implementation_and_synthetic_qualification_authorized": True,
        "HLT14_runtime_implemented": False,
        "pretrajectory_operation_authorized": False,
        "production_namespace_creation_authorized": False,
        "fresh_GR0_calibration_authorized": False,
        "candidate_execution_authorized": False,
        "physical_claim_authorized": False,
    }
    proof_keys = {
        "immutable_PROTO17_commit_and_all_bound_blobs_must_match",
        "independent_module_import_ast_must_exclude_PROTO17_implementation",
        "sealed_schema_bindings_and_qualification_must_be_rederived_not_trusted_by_name",
        "independent_fixture_must_match_all_three_sealed_construction_hashes",
        "all_listed_independent_and_binder_mutations_must_fail_closed",
        "four_future_HLT14_requirements_must_remain_required_and_uncompleted",
        "no_run_directory_or_raw_campaign_bundle_may_be_opened",
        "canonical_hash_bound_result_required",
    }
    expected_top = {
        "schema_version", "artifact_id", "project_version", "protocol_config", "protocol_freeze_config",
        "protocol_freeze_result", "protocol_freeze_reproducer", "protocol_schema_module",
        "protocol_construction_module", "protocol_freeze_document", "protocol_schema_test",
        "protocol_construction_test", "protocol_freeze_test", "binder_module", "scope",
        "immutable_lineage", "independent_theorem", "qualified_controls",
        "future_HLT14_runtime_requirements", "authorization", "proof_contract", "claims",
    }
    paths = {
        "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v17.toml",
        "protocol_freeze_config": "configs/fgc/fgc-1-pro17-frz1.toml",
        "protocol_freeze_result": "results/fgc-1-pro17-frz1.json",
        "protocol_freeze_reproducer": "scripts/reproduce_fgc_pro17_frz1.py",
        "protocol_schema_module": "src/recursive_horizons/fgc/evolution/protocol_v17.py",
        "protocol_construction_module": "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py",
        "protocol_freeze_document": "docs/fgc-pro17-frz1.md",
        "protocol_schema_test": "tests/test_fgc_protocol_v17.py",
        "protocol_construction_test": "tests/test_fgc_proto17_pure_construction.py",
        "protocol_freeze_test": "tests/test_fgc_pro17_frz1_reproduction.py",
        "binder_module": "src/recursive_horizons/fgc/evolution/proto17_construction_theorem.py",
    }
    if (
        set(value) != expected_top
        or (value.get("schema_version"), value.get("artifact_id"), value.get("project_version")) != (1, ARTIFACT, "0.11.0")
        or any(value.get(key) != path for key, path in paths.items())
        or value.get("scope") != expected_scope
        or value.get("immutable_lineage") != expected_lineage
        or set(value.get("independent_theorem", {})) != theorem_keys
        or any(item is not True for item in value.get("independent_theorem", {}).values())
        or value.get("qualified_controls") != qualified
        or value.get("future_HLT14_runtime_requirements") != FUTURE
        or value.get("authorization") != authorization
        or set(value.get("proof_contract", {})) != proof_keys
        or any(item is not True for item in value.get("proof_contract", {}).values())
    ):
        raise ValueError("PREF25 configuration contract differs")
    _claims(value.get("claims"))


def _sealed_bytes() -> dict[str, bytes]:
    _validate_commit(SEALED_COMMIT)
    answer: dict[str, bytes] = {}
    for path, expected in SEALED.items():
        historical = _git_show(path)
        live = ROOT / path
        _validate_sealed_blob(path, expected, historical, live)
        answer[path] = historical
    return answer


def _sealed_json(blob: bytes) -> dict[str, Any]:
    value = json.loads(blob.decode("utf-8"), object_pairs_hook=_reject_duplicates)
    if not isinstance(value, dict) or _canon(value).encode() != blob:
        raise ValueError("sealed FRZ1 result is not duplicate-safe canonical JSON")
    return value


def _schema(protocol: Mapping[str, Any]) -> dict[str, Any]:
    genesis = protocol.get("trusted_genesis", {})
    restart = protocol.get("restart_inputs", {})
    expected = {
        "genesis_spec_required_fields": list(GENESIS_FIELDS),
        "member_descriptor_required_fields": list(MEMBER_DESCRIPTOR_FIELDS),
        "state_object_required_fields": list(STATE_FIELDS),
        "physical_array_required_fields": list(ARRAY_FIELDS),
        "tdg6_ledger_required_fields": list(LEDGER_FIELDS),
        "member_keys_in_canonical_order": MEMBERS,
        "toml_bindings_match_core": True,
    }
    observed = {
        "genesis_spec_required_fields": genesis.get("genesis_spec_required_fields"),
        "member_descriptor_required_fields": genesis.get("member_descriptor_required_fields"),
        "state_object_required_fields": genesis.get("state_object_required_fields"),
        "physical_array_required_fields": genesis.get("physical_array_required_fields"),
        "tdg6_ledger_required_fields": genesis.get("tdg6_ledger_required_fields"),
        "member_keys_in_canonical_order": restart.get("member_keys_in_canonical_order"),
        "toml_bindings_match_core": True,
    }
    if observed != expected:
        raise ValueError("sealed PROTO17 schema bindings differ")
    return expected


def _no_oracle_import(theorem_source: bytes) -> bool:
    tree = ast.parse(theorem_source.decode("utf-8"))
    forbidden = {
        "recursive_horizons.fgc.evolution.protocol_v17",
        "recursive_horizons.fgc.evolution.proto17_pure_construction",
    }
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if (module in forbidden or module.split(".")[-1] in {"protocol_v17", "proto17_pure_construction"}
                    or any(name.name.split(".")[-1] in {"protocol_v17", "proto17_pure_construction"} for name in node.names)):
                return False
        if isinstance(node, ast.Import):
            if any(alias.name in forbidden or alias.name.endswith((".protocol_v17", ".proto17_pure_construction")) for alias in node.names):
                return False
    return True


def _future(value: object) -> dict[str, dict[str, bool]]:
    if value != FUTURE:
        raise ValueError("future HLT14 boundary differs")
    return deepcopy(FUTURE)


def _claims(value: object) -> dict[str, bool]:
    if not isinstance(value, Mapping) or set(value) != TRUE | {
        "HLT14_runtime_implemented", "PROTO17_successor_runtime_implemented",
        "PROTO17_pretrajectory_runtime_authorized", "PROTO17_output_namespace_created",
        "PROTO17_real_production_state_loaded", "PROTO17_raw_physical_array_bytes_recomputed",
        "PROTO17_historical_journal_payload_semantics_replayed",
        "PROTO17_filesystem_durability_or_restart_qualified",
        "PROTO17_external_Git_authority_validated_at_launch", "PROTO17_trajectory_read",
        "PROTO17_fresh_GR0_dynamic_calibration_authorized", "PROTO17_calibration_completed",
        "GR0_case_eligible", "classical_spherical_diagnostic_authorized",
        "SGBL_execution_authorized", "FGCQR_holdout_execution_authorized",
        "FGCQR_mechanism_rejected", "DEF1_execution_authorized",
        "retained_EFT_evolution_authorized", "physical_transition_claim_authorized",
        "general_gradient_route_rejected", "singularity_resolution_derived",
        "child_domain_or_topology_derived", "dark_sector_mechanism_derived",
        "varying_locally_measured_c_derived",
    }:
        raise ValueError("PREF25 claim key set differs")
    claims = dict(value)
    if any(claims[key] is not True for key in TRUE) or any(
        claims[key] is not False for key in claims if key not in TRUE
    ):
        raise ValueError("PREF25 claim boundary differs")
    return claims


def _binder_mutations(sealed_result: Mapping[str, Any], protocol: Mapping[str, Any]) -> dict[str, bool]:
    outcome: dict[str, bool] = {}
    try:
        _sealed_json(b'{"x":1,"x":2}')
    except ValueError:
        outcome["duplicate_JSON"] = True
    else:
        outcome["duplicate_JSON"] = False
    try:
        blob = _git_show("results/fgc-1-pro17-frz1.json")
        _validate_sealed_blob("results/fgc-1-pro17-frz1.json", "0" * 64, blob, ROOT / "results/fgc-1-pro17-frz1.json")
    except (ValueError, UnicodeDecodeError, json.JSONDecodeError):
        outcome["sealed_blob_or_hash"] = True
    else:
        outcome["sealed_blob_or_hash"] = False
    try:
        _validate_commit("0" * 40)
    except ValueError:
        outcome["sealed_commit"] = True
    else:
        outcome["sealed_commit"] = False
    mutated = deepcopy(protocol)
    mutated["trusted_genesis"]["state_object_required_fields"] = ["wrong"]
    try:
        _schema(mutated)
    except ValueError:
        outcome["schema_tuple"] = True
    else:
        outcome["schema_tuple"] = False
    mutated = deepcopy(protocol)
    mutated["restart_inputs"]["member_keys_in_canonical_order"] = list(reversed(MEMBERS))
    try:
        _schema(mutated)
    except ValueError:
        outcome["member_order"] = True
    else:
        outcome["member_order"] = False
    boundary = deepcopy(sealed_result["artifact_payload"]["mutation_contract"]["future_HLT14_runtime_requirements"])
    boundary["raw_bundle_byte_and_hash_recomputation"]["completed"] = True
    try:
        _future(boundary)
    except ValueError:
        outcome["deferred_completed"] = True
    else:
        outcome["deferred_completed"] = False
    source = (ROOT / _load_config()["binder_module"]).read_text(encoding="utf-8")
    forbidden = source + "\nimport recursive_horizons.fgc.evolution.protocol_v17 as decoy\n"
    outcome["forbidden_import"] = not _no_oracle_import(forbidden.encode("utf-8"))
    config = _load_config()
    claimed = deepcopy(config["claims"])
    claimed["HLT14_runtime_implemented"] = True
    try:
        _claims(claimed)
    except ValueError:
        outcome["runtime_claim_promotion"] = True
    else:
        outcome["runtime_claim_promotion"] = False
    return outcome


def _validate_sealed_frz1(
    sealed_result: Mapping[str, Any], freeze: Mapping[str, Any], schema: Mapping[str, Any], qualification: Mapping[str, Any],
) -> None:
    expected_controls = {name: True for name in CONTROL_NAMES}
    if (
        sealed_result.get("artifact_id") != "FGC-1-PRO17-FRZ1"
        or sealed_result.get("gate_status") != freeze.get("claims")
        or sealed_result.get("nonclaims") != {
            key: False for key, value in freeze.get("claims", {}).items() if value is False
        }
        or sealed_result.get("artifact_payload", {}).get("schema_bindings") != schema
        or sealed_result.get("artifact_payload", {}).get("synthetic_construction_qualification", {}).get("mutations_rejected") != expected_controls
        or sealed_result.get("artifact_payload", {}).get("mutation_contract", {}).get("qualified_pure_oracle_controls") != expected_controls
        or sealed_result.get("artifact_payload", {}).get("mutation_contract", {}).get("future_HLT14_runtime_requirements") != FUTURE
        or any(sealed_result.get("artifact_payload", {}).get("synthetic_construction_qualification", {}).get(key) != qualification.get(key) for key in (
            "genesis_checkpoint_sha256", "common_event_record_sha256", "successor_checkpoint_sha256",
        ))
    ):
        raise ValueError("sealed FRZ1 schema, controls, future boundary, or claim boundary differs")


def record() -> dict[str, Any]:
    config = _load_config()
    blobs = _sealed_bytes()
    protocol = tomllib.loads(blobs["configs/fgc/fgc-2-sf1-protocol-v17.toml"].decode("utf-8"))
    freeze = tomllib.loads(blobs["configs/fgc/fgc-1-pro17-frz1.toml"].decode("utf-8"))
    sealed_result = _sealed_json(blobs["results/fgc-1-pro17-frz1.json"])
    if not _no_oracle_import((ROOT / config["binder_module"]).read_bytes()):
        raise ValueError("independent theorem imports PROTO17 implementation")
    schema = _schema(protocol)
    qualification = theorem.qualification(
        blobs["configs/fgc/fgc-2-sf1-protocol-v17.toml"],
        blobs["configs/fgc/fgc-1-pro17-frz1.toml"],
    )
    _validate_sealed_frz1(sealed_result, freeze, schema, qualification)
    hash_keys = ("genesis_checkpoint_sha256", "common_event_record_sha256", "successor_checkpoint_sha256")
    sealed_qualification = sealed_result["artifact_payload"]["synthetic_construction_qualification"]
    if (
        not qualification.get("sealed_hashes_reproduced")
        or any(qualification.get(key) != sealed_qualification.get(key) for key in hash_keys)
        or qualification.get("qualified_mutation_names") != CONTROL_NAMES
        or qualification.get("qualified_mutation_count") != 14
        or qualification.get("all_mutations_rejected") is not True
        or sealed_qualification.get("qualified_mutation_names") != CONTROL_NAMES
        or sealed_qualification.get("qualified_mutation_count") != 14
        or sealed_qualification.get("all_qualified_mutations_rejected") is not True
    ):
        raise ValueError("independent theorem does not bind sealed construction controls")
    future = _future(sealed_result["artifact_payload"]["mutation_contract"]["future_HLT14_runtime_requirements"])
    binder_mutations = _binder_mutations(sealed_result, protocol)
    if not all(binder_mutations.values()):
        raise ValueError("binder mutation did not fail closed")
    claims = _claims(config["claims"])
    rows = [{"path": path, "sha256": digest, "matches_commit_and_worktree": True} for path, digest in SEALED.items()]
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT,
        "project_version": "0.11.0",
        "classification": "independent_PROTO17_schema_generation_zero_and_common_event_construction_binder_without_runtime_or_trajectory",
        "generated_by": "scripts/reproduce_fgc_pro17_pref25.py",
        "derivation_document": "docs/fgc-pro17-pref25.md",
        "derivation_document_sha256": _sha(DOC),
        "source_config_sha256": {"configs/fgc/fgc-1-pro17-pref25.toml": _sha(CONFIG)},
        "implementation_sha256": {
            "scripts/reproduce_fgc_pro17_pref25.py": _sha(Path(__file__)),
            "src/recursive_horizons/fgc/evolution/proto17_construction_theorem.py": _sha(ROOT / config["binder_module"]),
        },
        "scope_bindings": config["scope"],
        "artifact_payload": {
            "immutable_lineage": {"checkpoint_commit": SEALED_COMMIT, "records": rows},
            "sealed_proto17_contract": {"protocol": protocol, "freeze": freeze},
            "independent_schema_and_construction_audit": {
                "schema_bindings": schema,
                "theorem_imports_PROTO17_implementation": False,
                "independent_qualification": {key: qualification[key] for key in (
                    *hash_keys, "sealed_hashes_reproduced", "lone_suffix_reconstructs_exactly",
                    "store_plan_is_complete", "all_mutations_rejected", "qualified_mutation_names",
                    "qualified_mutation_count", "historical_record_payload_semantics_deferred",
                    "raw_bytes_deferred", "structural_only", "runtime_or_physical_claim",
                )},
            },
            "qualified_controls_audit": {
                "sealed_control_names": CONTROL_NAMES,
                "sealed_control_count": 14,
                "sealed_controls_all_true": True,
                "independent_equivalent_control_corpus_passed": True,
                "binder_specific_mutations": binder_mutations,
            },
            "future_HLT14_boundary": future,
            "decision": {
                "HLT14_implementation_and_synthetic_qualification_authorized": True,
                "HLT14_runtime_implemented": False,
                "pretrajectory_authorized": False,
                "production_namespace_creation_authorized": False,
                "calibration_authorized": False,
                "candidate_execution_authorized": False,
            },
        },
        "gate_status": claims,
        "nonclaims": {key: False for key in claims if key not in TRUE},
    }


def verify(path: Path = OUTPUT) -> None:
    source = path.read_text(encoding="utf-8")
    observed = json.loads(source, object_pairs_hook=_reject_duplicates)
    if source != _canon(observed) or observed != record():
        raise ValueError("PREF25 canonical result differs")


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.verify:
        verify(args.output)
        print(f"verified {args.output}")
    else:
        _atomic_write(args.output, _canon(record()).encode("utf-8"))
        print(f"wrote {args.output}")
