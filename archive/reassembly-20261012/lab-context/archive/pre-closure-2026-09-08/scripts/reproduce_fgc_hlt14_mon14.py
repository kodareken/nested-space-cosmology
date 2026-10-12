#!/usr/bin/env python3
"""Reproduce HLT14's completed synthetic-only PROTO17 qualification.

The reproducer reads immutable compact lineage and creates only temporary test
Git repositories, NPZ bundles, and lifecycle stores.  It never opens a real
campaign bundle or materializes a production ``runs/`` namespace.
"""

from __future__ import annotations

import argparse
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

from recursive_horizons.fgc.evolution.proto17_hlt14_qualification import (  # noqa: E402
    SEALED_HASHES,
    synthetic_qualification,
)


ARTIFACT_ID = "FGC-1-HLT14-MON14"
PROJECT_VERSION = "0.11.0"
BASE = "4da668fad42d4006d3e09f511eafa0ef0e69d7ed"
CONFIG = ROOT / "configs/fgc/fgc-1-hlt14-mon14.toml"
OUTPUT = ROOT / "results/fgc-1-hlt14-mon14.json"
DOC = ROOT / "docs/fgc-hlt14-mon14.md"

MEMBERS = [
    "RK4-2049", "RK4-4097", "RK4-8193",
    "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
]
LINEAGE = {
    "pref25_config_sha256": (
        "configs/fgc/fgc-1-pro17-pref25.toml",
        "60e78ac3778bd0411bc52dc8756396a18b8f3f1be54b18938fa5cdcb83c59279",
    ),
    "pref25_result_sha256": (
        "results/fgc-1-pro17-pref25.json",
        "7bba3af7b18798496800cdc024168c5593d6cb39ac643c5933d4df101191b199",
    ),
    "pref25_reproducer_sha256": (
        "scripts/reproduce_fgc_pro17_pref25.py",
        "dcd0c0eee47724133b91e66b03d844c24fb23b489f05ebbf7c73438510431279",
    ),
    "pref25_theorem_sha256": (
        "src/recursive_horizons/fgc/evolution/proto17_construction_theorem.py",
        "59b1eec6bcbde0bded065b281f6f2dc94615b0fd86e94e15e96a9ff756eb18ca",
    ),
    "pref25_document_sha256": (
        "docs/fgc-pro17-pref25.md",
        "ca9ae473fffb268de12a3542211868a88cfcdb5150a15a85a65ef3cb4fcd3e04",
    ),
    "pref25_test_sha256": (
        "tests/test_fgc_pro17_pref25_reproduction.py",
        "c1b2891459850dfae285982781be27e9d0490ff73f1275767c2c14fd2de3fcca",
    ),
    "proto17_protocol_sha256": (
        "configs/fgc/fgc-2-sf1-protocol-v17.toml",
        "520466fb2f798159f7d0183938327eec0fa1ccf8b49c41101924f1ff01fcd484",
    ),
    "proto17_freeze_config_sha256": (
        "configs/fgc/fgc-1-pro17-frz1.toml",
        "6ff3cfad5eebb7824434e3b7f97da3ce7a8610312523e3461c8999532678d2ec",
    ),
    "proto17_freeze_result_sha256": (
        "results/fgc-1-pro17-frz1.json",
        "8e4b3f77015374cf825c141dff34d343eb3d5f4ccebfed5082bc283c1d414fea",
    ),
    "proto17_freeze_reproducer_sha256": (
        "scripts/reproduce_fgc_pro17_frz1.py",
        "1481d8da60d2684a692fe9526236a0b82e6d63299e922545c9b9ee75b9c1a3cf",
    ),
    "proto17_protocol_core_sha256": (
        "src/recursive_horizons/fgc/evolution/protocol_v17.py",
        "2933d60d90607e0361d7deaec5395d8ad70a15ec41ab98ba2dd6aedbe9c5085c",
    ),
    "proto17_construction_sha256": (
        "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py",
        "5fbd7f990d9f821344bbaae80840070386b3ae483fc278a85ae2a3ee5de8fa1d",
    ),
}
FUTURE = {
    "raw_bundle_byte_and_hash_recomputation": {"required": True, "completed": False},
    "external_Git_authority_blob_and_live_import_validation": {"required": True, "completed": False},
    "semantic_historical_journal_payload_replay": {"required": True, "completed": False},
    "namespace_reuse_and_foreign_store_rejection": {"required": True, "completed": False},
}
SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO17",
    "runtime_owner": ARTIFACT_ID,
    "role": "completed_synthetic_only_PROTO17_authority_bundle_and_durable_runtime_qualification",
    "sealed_compact_inputs_read_from_immutable_commit": True,
    "synthetic_real_NPZ_bundles_opened": True,
    "synthetic_Git_authority_preflight_exercised": True,
    "synthetic_temporary_store_materialized": True,
    "real_production_state_or_raw_bundle_opened": False,
    "real_historical_journal_payload_semantics_replayed": False,
    "external_production_launch_authority_instance_validated": False,
    "production_namespace_materialized": False,
    "pretrajectory_operation": False,
    "trajectory_read": False,
    "mechanism_question_answered": False,
}
ADAPTER_CONTRACT = {
    "canonical_member_order": MEMBERS,
    "generation_zero_event_index": 23,
    "generation_zero_accepted_boundary_time": "23/16",
    "generation_zero_active_target_time": "3/2",
    "common_event_completed_index": 24,
    "common_event_next_target_time": "25/16",
    "expected_generation_zero_checkpoint_sha256": SEALED_HASHES["generation_zero_checkpoint_sha256"],
    "expected_common_event_receipt_sha256": SEALED_HASHES["common_event_receipt_sha256"],
    "expected_successor_checkpoint_sha256": SEALED_HASHES["successor_checkpoint_sha256"],
    "sealed_abstract_fixture_must_match_pure_constructor_byte_for_byte": True,
    "real_NPZ_fixture_hashes_are_input_dependent": True,
    "all_six_raw_bundle_bytes_and_array_hashes_must_be_recomputed": True,
    "authority_commit_blob_and_three_source_roles_must_be_checked": True,
    "import_origin_check_is_preflight_not_TOCTOU_closure": True,
    "synthetic_store_must_be_temporary_and_isolated": True,
    "real_runs_namespace_must_not_be_created": True,
}
QUALIFICATION_CONTRACT = {
    "qualification_case_count": 72,
    "nominal_case_count": 4,
    "bundle_admission_case_count": 6,
    "mutation_case_count": 18,
    "in_process_exception_fault_case_count": 44,
    "qualification_case_digest": "445765193eddba247b7452635ca382388950a55a3a2bfb128e0b125f530b6baf",
    "complete_visible_synthetic_genesis_and_store_plan_verified": True,
    "content_addressed_state_store_and_checkpoint_reload_verified": True,
    "structural_six_record_accepted_fine_prefix_verified": True,
    "predecessor_only_common_event_and_complete_successor_checkpoint_verified": True,
    "exact_lone_receipt_recovery_only_verified": True,
    "fully_rehashed_semantic_plan_and_high_water_forgeries_rejected": True,
    "raw_hash_array_shape_metadata_canonicality_and_symlink_mutations_rejected": True,
    "authority_commit_schema_source_drift_and_import_origin_preflight_exercised": True,
    "accepted_prefix_common_event_and_generation_zero_exception_windows_exercised": True,
    "preexisting_root_foreign_store_and_production_path_guards_exercised": True,
    "sealed_PREF25_hash_vector_reproduced": True,
    "all_qualification_cases_passed": True,
    "faults_are_in_process_exceptions_not_power_loss_proof": True,
}
AUTHORIZATION = {
    "HLT14_implementation_and_synthetic_qualification_authorized": True,
    "HLT14_runtime_implemented": True,
    "HLT14_synthetic_qualification_passed": True,
    "production_prelaunch_successor_design_authorized": True,
    "pretrajectory_operation_authorized": False,
    "production_namespace_creation_authorized": False,
    "fresh_GR0_calibration_authorized": False,
    "candidate_execution_authorized": False,
    "physical_claim_authorized": False,
}
TRUE_CLAIMS = {
    "PROTO17_Pref25_authorization_consumed",
    "HLT14_synthetic_runtime_contract_frozen",
    "HLT14_implementation_and_synthetic_qualification_authorized",
    "HLT14_runtime_implemented",
    "HLT14_synthetic_qualification_passed",
    "HLT14_synthetic_Git_authority_and_import_origin_preflight_exercised",
    "HLT14_synthetic_all_six_raw_bundle_admission_exercised",
    "HLT14_synthetic_semantic_replay_and_exception_recovery_qualified",
    "HLT14_sealed_PREF25_hash_vector_reproduced",
    "production_prelaunch_successor_design_authorized",
}
FALSE_CLAIMS = {
    "HLT14_real_production_state_loaded",
    "HLT14_real_raw_bundle_bytes_recomputed",
    "HLT14_real_historical_journal_payload_semantics_replayed",
    "HLT14_external_production_launch_authority_validated",
    "HLT14_production_namespace_reuse_or_foreign_store_rejection_validated",
    "HLT14_production_namespace_materialized",
    "HLT14_production_filesystem_durability_or_restart_qualified",
    "HLT14_power_loss_durability_qualified",
    "HLT14_pretrajectory_operation_authorized",
    "HLT14_trajectory_read",
    "HLT14_fresh_GR0_dynamic_calibration_authorized",
    "HLT14_calibration_completed",
    "GR0_case_eligible",
    "classical_spherical_diagnostic_authorized",
    "SGBL_execution_authorized",
    "FGCQR_holdout_execution_authorized",
    "DEF1_execution_authorized",
    "retained_EFT_evolution_authorized",
    "physical_transition_claim_authorized",
    "FGCQR_mechanism_rejected",
    "general_gradient_route_rejected",
    "singularity_resolution_derived",
    "child_domain_or_topology_derived",
    "dark_sector_mechanism_derived",
    "varying_locally_measured_c_derived",
}
CLAIMS = {**{key: True for key in TRUE_CLAIMS}, **{key: False for key in FALSE_CLAIMS}}
ROOT_KEYS = {
    "schema_version", "artifact_id", "project_version", "protocol_config",
    "protocol_freeze_config", "protocol_freeze_result", "pref25_config",
    "pref25_result", "input_adapter_module", "runtime_module",
    "qualification_module", "scope", "immutable_lineage",
    "synthetic_adapter_contract", "synthetic_qualification",
    "future_HLT14_production_instance_requirements", "authorization", "claims",
}
IMPLEMENTATION_PATHS = (
    "scripts/reproduce_fgc_hlt14_mon14.py",
    "src/recursive_horizons/fgc/evolution/proto17_hlt14_inputs.py",
    "src/recursive_horizons/fgc/evolution/proto17_hlt14_runtime.py",
    "src/recursive_horizons/fgc/evolution/proto17_hlt14_qualification.py",
    "tests/test_fgc_proto17_hlt14_inputs.py",
    "tests/test_fgc_proto17_hlt14_runtime.py",
    "tests/test_fgc_hlt14_mon14_reproduction.py",
)


def _relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False,
    ) + "\n"


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def validate_config_mapping(config: Mapping[str, Any]) -> dict[str, Any]:
    """Enforce HLT14's entire authorization and nonclaim boundary exactly."""
    if not isinstance(config, Mapping) or set(config) != ROOT_KEYS:
        raise ValueError("HLT14 config root differs")
    value = dict(config)
    if (
        value["schema_version"], value["artifact_id"], value["project_version"]
    ) != (1, ARTIFACT_ID, PROJECT_VERSION):
        raise ValueError("HLT14 config identity differs")
    expected_paths = {
        "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v17.toml",
        "protocol_freeze_config": "configs/fgc/fgc-1-pro17-frz1.toml",
        "protocol_freeze_result": "results/fgc-1-pro17-frz1.json",
        "pref25_config": "configs/fgc/fgc-1-pro17-pref25.toml",
        "pref25_result": "results/fgc-1-pro17-pref25.json",
        "input_adapter_module": "src/recursive_horizons/fgc/evolution/proto17_hlt14_inputs.py",
        "runtime_module": "src/recursive_horizons/fgc/evolution/proto17_hlt14_runtime.py",
        "qualification_module": "src/recursive_horizons/fgc/evolution/proto17_hlt14_qualification.py",
    }
    if any(value[key] != expected for key, expected in expected_paths.items()):
        raise ValueError("HLT14 input paths differ")
    if value["scope"] != SCOPE:
        raise ValueError("HLT14 scope boundary differs")
    expected_lineage = {
        "checkpoint_commit": BASE,
        "checkpoint_must_be_ancestor_of_HEAD": True,
        "all_bound_blobs_must_match_checkpoint_commit_and_worktree": True,
        **{key: digest for key, (_path, digest) in LINEAGE.items()},
    }
    if value["immutable_lineage"] != expected_lineage:
        raise ValueError("HLT14 immutable lineage declaration differs")
    if value["synthetic_adapter_contract"] != ADAPTER_CONTRACT:
        raise ValueError("HLT14 adapter contract differs")
    if value["synthetic_qualification"] != QUALIFICATION_CONTRACT:
        raise ValueError("HLT14 qualification contract differs")
    if value["future_HLT14_production_instance_requirements"] != FUTURE:
        raise ValueError("HLT14 future production duties differ")
    if value["authorization"] != AUTHORIZATION:
        raise ValueError("HLT14 authorization boundary differs")
    if value["claims"] != CLAIMS:
        raise ValueError("HLT14 claim boundary differs")
    return value


def load_config(path: Path = CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    return validate_config_mapping(value)


def _git_show_sha(commit: str, relative: str) -> str:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{relative}"], cwd=ROOT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if completed.returncode:
        raise ValueError(f"missing immutable HLT14 input: {relative}")
    return sha256(completed.stdout).hexdigest()


def _sealed_lineage(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", BASE, "HEAD"],
        cwd=ROOT, check=False,
    ).returncode:
        raise ValueError("HLT14 checkpoint is not an ancestor of HEAD")
    records: list[dict[str, Any]] = []
    for key, (relative, expected) in LINEAGE.items():
        if config["immutable_lineage"][key] != expected:
            raise ValueError(f"HLT14 lineage declaration differs: {relative}")
        live = ROOT / relative
        if (
            not live.is_file()
            or _sha(live) != expected
            or _git_show_sha(BASE, relative) != expected
        ):
            raise ValueError(f"HLT14 sealed/live lineage differs: {relative}")
        records.append({
            "path": relative,
            "sha256": expected,
            "matches_checkpoint_and_worktree": True,
        })
    return records


def _qualification_summary(config: Mapping[str, Any]) -> dict[str, Any]:
    summary = synthetic_qualification()
    required = {
        "case_names", "case_count", "case_digest", "all_cases_passed",
        "sealed_PREF25_hash_vector", "sealed_PREF25_hash_vector_reproduced",
        "real_NPZ_authority_path_has_input_dependent_hashes",
        "temporary_directories_only",
        "in_process_exception_faults_not_power_loss_proof",
        "real_production_instance_requirements_completed",
    }
    if not isinstance(summary, dict) or set(summary) != required:
        raise ValueError("HLT14 qualification summary schema differs")
    names = summary["case_names"]
    expected = config["synthetic_qualification"]
    if (
        not isinstance(names, list)
        or names != sorted(names)
        or len(names) != len(set(names))
        or summary["case_count"] != expected["qualification_case_count"]
        or summary["case_count"] != len(names)
        or summary["case_digest"] != expected["qualification_case_digest"]
        or sum(name.startswith("nominal:") for name in names) != expected["nominal_case_count"]
        or sum(name.startswith("bundle:") for name in names) != expected["bundle_admission_case_count"]
        or sum(name.startswith("mutation:") for name in names) != expected["mutation_case_count"]
        or sum(name.startswith("fault:") for name in names) != expected["in_process_exception_fault_case_count"]
        or summary["sealed_PREF25_hash_vector"] != SEALED_HASHES
        or summary["sealed_PREF25_hash_vector_reproduced"] is not True
        or summary["real_NPZ_authority_path_has_input_dependent_hashes"] is not True
        or summary["temporary_directories_only"] is not True
        or summary["in_process_exception_faults_not_power_loss_proof"] is not True
        or summary["real_production_instance_requirements_completed"] is not False
        or summary["all_cases_passed"] is not True
    ):
        raise ValueError("HLT14 executed qualification differs")
    return summary


def record(config_path: Path = CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _sealed_lineage(config)
    summary = _qualification_summary(config)
    claims = dict(config["claims"])
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "completed_synthetic_only_PROTO17_authority_bundle_and_durable_runtime_qualification_without_production_prelaunch_or_trajectory",
        "generated_by": _relative(Path(__file__)),
        "derivation_document": _relative(DOC),
        "derivation_document_sha256": _sha(DOC),
        "source_config_sha256": {
            _relative(config_path): _sha(config_path),
            config["protocol_config"]: _sha(ROOT / config["protocol_config"]),
            config["protocol_freeze_config"]: _sha(ROOT / config["protocol_freeze_config"]),
            config["pref25_config"]: _sha(ROOT / config["pref25_config"]),
        },
        "implementation_sha256": {
            relative: _sha(ROOT / relative) for relative in IMPLEMENTATION_PATHS
        },
        "scope_bindings": dict(config["scope"]),
        "artifact_payload": {
            "immutable_lineage": {
                "evidence_checkpoint_commit": BASE,
                "checkpoint_is_ancestor_of_HEAD": True,
                "tracked_compact_records": lineage,
            },
            "synthetic_adapter_contract": dict(config["synthetic_adapter_contract"]),
            "synthetic_qualification_contract": dict(config["synthetic_qualification"]),
            "synthetic_qualification_summary": summary,
            "evidence_partition": {
                "sealed_abstract_fixture_hash_vector_reproduced": True,
                "real_NPZ_fixture_hashes_are_input_dependent": True,
                "synthetic_Git_authority_and_import_origin_preflight_exercised": True,
                "external_production_launch_authority_validated": False,
                "in_process_exception_recovery_exercised": True,
                "power_loss_durability_qualified": False,
            },
            "future_HLT14_boundary": FUTURE,
            "decision": {
                "HLT14_runtime_implemented": True,
                "HLT14_synthetic_qualification_passed": True,
                "production_prelaunch_successor_design_authorized": True,
                "pretrajectory_operation_authorized": False,
                "production_namespace_creation_authorized": False,
                "fresh_GR0_calibration_authorized": False,
                "candidate_execution_authorized": False,
                "physical_claim_authorized": False,
            },
        },
        "gate_status": claims,
        "nonclaims": {key: False for key in FALSE_CLAIMS},
    }


def validate_result_mapping(
    observed: Mapping[str, Any], config_path: Path = CONFIG,
) -> dict[str, Any]:
    if not isinstance(observed, Mapping) or dict(observed) != record(config_path):
        raise ValueError("HLT14 result differs from fresh synthetic reproduction")
    return dict(observed)


def load_canonical_result(path: Path = OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    observed = json.loads(source, object_pairs_hook=_reject_duplicates)
    if (
        not isinstance(observed, dict)
        or observed.get("artifact_id") != ARTIFACT_ID
        or source != _canonical(observed)
    ):
        raise ValueError("HLT14 canonical result differs")
    return observed


def verify(path: Path = OUTPUT, config_path: Path = CONFIG) -> None:
    validate_result_mapping(load_canonical_result(path), config_path)


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent,
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        parent = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(parent)
        finally:
            os.close(parent)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    if arguments.verify:
        verify(arguments.output, arguments.config)
        print(f"verified {arguments.output}")
    else:
        _atomic_write(
            arguments.output,
            _canonical(record(arguments.config)).encode("utf-8"),
        )
        print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
