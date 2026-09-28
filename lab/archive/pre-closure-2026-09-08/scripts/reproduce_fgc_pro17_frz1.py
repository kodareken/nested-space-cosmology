#!/usr/bin/env python3
"""Reproduce the schema-only PROTO17 construction-contract freeze.

This command deliberately never opens an historical restart bundle, fabricates
a physical GenesisSpec, imports a future runtime, or creates an output root.
It freezes the inputs a future HLT14 authority must expose before it can do any
of those things.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import tempfile
import tomllib
from typing import Any, Mapping

import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from recursive_horizons.fgc.evolution.protocol_v17 import (  # noqa: E402
    ARRAY_DESCRIPTOR_FIELDS, CADENCE_RATIONAL, GENESIS_INPUT_FIELDS,
    LEDGER_FIELDS, MEMBER_DESCRIPTOR_FIELDS, MEMBER_KEYS, RECEIPT_FIELDS,
    STATE_OBJECT_FIELDS, validate_proto17_schema,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import (  # noqa: E402
    run_pure_construction_qualification,
)

CONFIG = ROOT / "configs/fgc/fgc-1-pro17-frz1.toml"
PROTOCOL = ROOT / "configs/fgc/fgc-2-sf1-protocol-v17.toml"
OUTPUT = ROOT / "results/fgc-1-pro17-frz1.json"
DOCUMENT = ROOT / "docs/fgc-pro17-frz1.md"
PROTOCOL_MODULE = ROOT / "src/recursive_horizons/fgc/evolution/protocol_v17.py"
ORACLE_MODULE = ROOT / "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py"

LINEAGE = {
    "results/fgc-1-hlt13-mon13.json": "de55f8d55a1db0043945d3553c7d2070703809b82c1f7141bf923c036b9df5b7",
    "configs/fgc/fgc-2-sf1-protocol-v16.toml": "3b0ab8663ed0568f03841167a2db096419b0978a71b534eefa9f9bb9a6787c75",
    "results/fgc-1-pro16-frz1.json": "f8e14deda1638bda6e98831ae1ea80c8d1fdffb23957e183acd8b57a154ce43f",
    "configs/fgc/fgc-1-pro16-pref24.toml": "c8d983513d123b9d04b8a0574bfea1ce2184f9c9cd1c853a156830cb60f9dc10",
    "results/fgc-1-pro16-pref24.json": "168f995ee03851d62154661a8575d4af0aa068af128be5874f8e333e4341995f",
    "src/recursive_horizons/fgc/evolution/proto15_runtime.py": "73f07ae26e90642bc6e5158d9a85bfa499d8082c13e082fc8f19e968d8e624bc",
}
TRUE = {
    "PROTO17_frozen",
    "PROTO17_exact_common_event_successor_construction_frozen",
    "PROTO17_external_trusted_genesis_construction_frozen",
    "PROTO16_under_specification_preserved",
}


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n"


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} is not object-valued")
    return value


def _git_show(commit: str, relative: str) -> bytes:
    return subprocess.run(
        ["git", "show", f"{commit}:{relative}"], cwd=ROOT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
    ).stdout


def _assert_namespace_absent(protocol: Mapping[str, Any]) -> list[dict[str, Any]]:
    namespace = protocol["namespace"]
    roots = [namespace["calibration_output_root"], namespace["holdout_output_root"]]
    records = []
    for relative in roots:
        path = ROOT / relative
        if path.exists():
            raise ValueError(f"PROTO17 freeze requires absent namespace: {relative}")
        records.append({"path": relative, "absent_before_freeze": True})
    return records


def validate_contract(config: Mapping[str, Any], protocol: Mapping[str, Any]) -> None:
    if config.get("artifact_id") != "FGC-1-PRO17-FRZ1" or protocol.get("artifact_id") != "FGC-2-SF1-PROTO17":
        raise ValueError("PROTO17 identity differs")
    if protocol.get("frozen") is not True or protocol.get("protocol_version") != 17:
        raise ValueError("PROTO17 freeze identity differs")
    restart = protocol.get("restart_inputs", {})
    if (
        restart.get("restart_coordinate_time") != "23/16"
        or restart.get("genesis_accepted_boundary_time") != "23/16"
        or restart.get("genesis_active_event_target_time") != "3/2"
    ):
        raise ValueError("PROTO17 genesis time identities differ")
    if config.get("claims") != protocol.get("claims"):
        raise ValueError("PROTO17 claims differ between protocol and freeze")
    claims = config["claims"]
    if any(claims.get(key) is not True for key in TRUE) or any(
        value is not False for key, value in claims.items() if key not in TRUE
    ):
        raise ValueError("PROTO17 claim boundary differs")
    genesis = protocol.get("trusted_genesis", {})
    schema_bindings = {
        "genesis_spec_required_fields": GENESIS_INPUT_FIELDS,
        "member_descriptor_required_fields": MEMBER_DESCRIPTOR_FIELDS,
        "state_object_required_fields": STATE_OBJECT_FIELDS,
        "physical_array_required_fields": ARRAY_DESCRIPTOR_FIELDS,
        "tdg6_ledger_required_fields": LEDGER_FIELDS,
    }
    if (
        not isinstance(genesis, Mapping)
        or any(tuple(genesis.get(name, ())) != expected for name, expected in schema_bindings.items())
        or tuple(restart.get("member_keys_in_canonical_order", ())) != MEMBER_KEYS
        or "genesis_spec_sha256" in genesis.get("genesis_spec_required_fields", [])
    ):
        raise ValueError("PROTO17 exact schema bindings differ")
    for key in (
        "all_state_objects_and_ledgers_are_visible_construction_inputs",
        "physical_array_raw_byte_hashes_and_semantic_descriptors_are_required_as_declared_external_inputs",
        "physical_array_raw_bytes_are_not_opened_or_recomputed_by_this_freeze",
        "physical_state_sha256_is_u_p_q_array_content_hash_not_cursor_content_address",
        "physical_state_and_restart_payload_hashes_are_not_recomputed_without_the_raw_bundle",
        "state_object_canonical_sha256_is_cursor_and_checkpoint_accepted_state_content_address",
        "physical_input_manifest_sha256_is_sha256_of_canonical_embedded_manifest_without_its_sibling_field",
        "generation_zero_checkpoint_is_byte_derivable_without_namespace_generated_value",
        "authorization_result_must_be_loaded_from_authorization_commit_not_worktree",
        "foreign_self_consistent_store_is_not_authenticated_without_external_authority",
    ):
        if genesis.get(key) is not True:
            raise ValueError(f"PROTO17 genesis rule differs: {key}")
    successor = protocol.get("common_event_successor", {})
    if not all(successor.get(key) is True for key in (
        "receipt_is_predecessor_only", "outer_record_sha256_must_exist_before_successor_construction",
        "successor_high_water_maps_equal_new_cursor_values",
        "only_one_lone_durable_common_event_suffix_is_reconstructible",
        "receipt_sequence_is_hash_contiguous_supplied_prior_journal_length_plus_one",
    )):
        raise ValueError("PROTO17 successor construction contract differs")
    proof = config.get("proof_contract", {})
    if not proof or any(value is not True for value in proof.values()):
        raise ValueError("PROTO17 proof contract differs")


def _schema_bindings(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Serialize the exact core/TOML vocabulary used by this freeze."""
    genesis = protocol["trusted_genesis"]
    return {
        "genesis_spec_required_fields": list(GENESIS_INPUT_FIELDS),
        "member_descriptor_required_fields": list(MEMBER_DESCRIPTOR_FIELDS),
        "state_object_required_fields": list(STATE_OBJECT_FIELDS),
        "physical_array_required_fields": list(ARRAY_DESCRIPTOR_FIELDS),
        "tdg6_ledger_required_fields": list(LEDGER_FIELDS),
        "member_keys_in_canonical_order": list(MEMBER_KEYS),
        "toml_bindings_match_core": all(
            tuple(genesis[name]) == expected
            for name, expected in (
                ("genesis_spec_required_fields", GENESIS_INPUT_FIELDS),
                ("member_descriptor_required_fields", MEMBER_DESCRIPTOR_FIELDS),
                ("state_object_required_fields", STATE_OBJECT_FIELDS),
                ("physical_array_required_fields", ARRAY_DESCRIPTOR_FIELDS),
                ("tdg6_ledger_required_fields", LEDGER_FIELDS),
            )
        ),
    }


def _synthetic_oracle_schema() -> dict[str, Any]:
    """The pure oracle's own small schema, separate from production authority."""
    return {
        "schema_version": 1,
        "artifact_id": "FGC-2-SF1-PROTO17",
        "protocol_version": 17,
        "frozen": True,
        "member_keys": list(MEMBER_KEYS),
        "cadence": CADENCE_RATIONAL,
        "genesis_input_fields": list(GENESIS_INPUT_FIELDS),
        "member_descriptor_fields": list(MEMBER_DESCRIPTOR_FIELDS),
        "receipt_fields": list(RECEIPT_FIELDS),
        "claims": {
            "PROTO17_pure_reference_oracle_implemented": True,
            "PROTO17_runtime_implemented": False,
            "PROTO17_namespace_authorized": False,
            "PROTO17_trajectory_read": False,
        },
    }


def record(*, require_absent_namespace: bool) -> dict[str, Any]:
    config = _load_toml(CONFIG)
    protocol = _load_toml(PROTOCOL)
    validate_contract(config, protocol)
    validate_proto17_schema(_synthetic_oracle_schema())
    qualification = run_pure_construction_qualification()
    if (
        qualification.get("side_effect_free_reference_oracle_only") is not True
        or qualification.get("lone_suffix_reconstructs_byte_identically") is not True
        or qualification.get("all_qualified_mutations_rejected") is not True
        or qualification.get("qualified_mutation_count") != 14
    ):
        raise ValueError("PROTO17 synthetic pure-construction qualification failed")
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=ROOT).returncode:
        raise ValueError("PROTO17 lineage checkpoint is not an ancestor of HEAD")
    records = []
    for relative, expected in LINEAGE.items():
        historical = _git_show(commit, relative)
        live = ROOT / relative
        if sha256(historical).hexdigest() != expected or _sha(live) != expected:
            raise ValueError(f"sealed/live lineage differs: {relative}")
        records.append({"path": relative, "sha256": expected, "matches_commit_and_worktree": True})
    temporal = (
        _assert_namespace_absent(protocol)
        if require_absent_namespace
        else [{"path": path, "absent_before_freeze": True} for path in (
            protocol["namespace"]["calibration_output_root"], protocol["namespace"]["holdout_output_root"],
        )]
    )
    return {
        "schema_version": 1,
        "artifact_id": "FGC-1-PRO17-FRZ1",
        "project_version": "0.11.0",
        "classification": "prospective_exact_successor_and_external_generation_zero_construction_freeze_without_runtime_or_trajectory",
        "generated_by": "scripts/reproduce_fgc_pro17_frz1.py",
        "derivation_document": "docs/fgc-pro17-frz1.md",
        "derivation_document_sha256": _sha(DOCUMENT),
        "source_config_sha256": {
            "configs/fgc/fgc-2-sf1-protocol-v17.toml": _sha(PROTOCOL),
            "configs/fgc/fgc-1-pro17-frz1.toml": _sha(CONFIG),
        },
        "implementation_sha256": {
            "scripts/reproduce_fgc_pro17_frz1.py": _sha(Path(__file__)),
            "src/recursive_horizons/fgc/evolution/protocol_v17.py": _sha(PROTOCOL_MODULE),
            "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py": _sha(ORACLE_MODULE),
        },
        "scope_bindings": config["scope"],
        "artifact_payload": {
            "immutable_lineage": {"checkpoint_commit": commit, "records": records},
            "schema_bindings": _schema_bindings(protocol),
            "proto17_contract": {
                "common_event_successor": protocol["common_event_successor"],
                "trusted_genesis": protocol["trusted_genesis"],
                "physical_state_sha256_is_not_persisted_state_object_content_address": True,
                "cursor_and_checkpoint_accepted_state_sha256_mean_state_object_canonical_sha256": True,
                "concrete_production_GenesisSpec_or_physical_bundle_supplied": False,
                "concrete_generation_zero_checkpoint_supplied": False,
            },
            "synthetic_construction_qualification": qualification,
            "temporal_namespace_precondition": {
                "both_new_namespaces_absent_at_freeze_boundary": True,
                "freeze_created_no_namespace": True,
                "records": temporal,
            },
            "mutation_contract": {
                "qualified_pure_oracle_controls": qualification["mutations_rejected"],
                "future_HLT14_runtime_requirements": {
                    "raw_bundle_byte_and_hash_recomputation": {
                        "required": True, "completed": False,
                    },
                    "external_Git_authority_blob_and_live_import_validation": {
                        "required": True, "completed": False,
                    },
                    "semantic_historical_journal_payload_replay": {
                        "required": True, "completed": False,
                    },
                    "namespace_reuse_and_foreign_store_rejection": {
                        "required": True, "completed": False,
                    },
                },
            },
            "decision": {
                "PROTO17_frozen": True,
                "successor_runtime_owner": "FGC-1-HLT14-MON14",
                "concrete_GenesisSpec_and_physical_bundle_remain_future_HLT14_adapter_inputs": True,
                "runtime_implemented": False,
                "pretrajectory_authorized": False,
                "calibration_authorized": False,
                "candidate_execution_authorized": False,
            },
        },
        "gate_status": config["claims"],
        "nonclaims": {key: False for key in config["claims"] if key not in TRUE},
    }


def verify(path: Path = OUTPUT) -> None:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_reject_duplicates)
    if source != _canonical(value) or value != record(require_absent_namespace=False):
        raise ValueError("PROTO17 canonical result differs")


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
    parser.add_argument("--verify-freeze-boundary", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.verify:
        verify(args.output)
        print(f"verified {args.output}")
    elif args.verify_freeze_boundary:
        record(require_absent_namespace=True)
        print("verified PROTO17 freeze namespace boundary")
    else:
        _atomic_write(args.output, _canonical(record(require_absent_namespace=True)).encode("utf-8"))
        print(f"wrote {args.output}")
