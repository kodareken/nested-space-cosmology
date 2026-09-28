#!/usr/bin/env python3
"""Reproduce the raw-dependent, read-only FGC-1-PRO18-PREF26 binder."""

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
sys.path[:0] = [str(ROOT / "src")]
from recursive_horizons.fgc.evolution.proto18_pref26_binder import (  # noqa: E402
    Proto18Pref26BinderError,
    bind_pref26,
    validate_pref26_contract,
)


ARTIFACT = "FGC-1-PRO18-PREF26"
SEALED_COMMIT = "8a41a70606dbdc0092e45eed215fd108b35fac9b"
CONFIG = ROOT / "configs/fgc/fgc-1-pro18-pref26.toml"
PRO18_CONFIG = ROOT / "configs/fgc/fgc-1-pro18-frz1.toml"
OUTPUT = ROOT / "results/fgc-1-pro18-pref26.json"
DOC = ROOT / "docs/fgc-pro18-pref26.md"
SCRIPT = ROOT / "scripts/reproduce_fgc_pro18_pref26.py"
IMPLEMENTATION_PATHS = (
    "src/recursive_horizons/fgc/evolution/proto18_production_inputs.py",
    "src/recursive_horizons/fgc/evolution/proto18_historical_replay.py",
    "src/recursive_horizons/fgc/evolution/proto18_pref26_binder.py",
)
SEALED_BLOBS = {
    "configs/fgc/fgc-1-pro18-frz1.toml": "f0997e28bd92c013dac34e8c4c0bc39de4ac7db882a6dc40708724bf58df4a60",
    "results/fgc-1-pro18-frz1.json": "d9aaad258ee9a6faddef445a998fda0faeb925351b412f2ae714b1ff41238698",
    "scripts/reproduce_fgc_pro18_frz1.py": "2159834037ced5c57fbe7ced91ccfb83eec5b98d62232bef15a2d56014029d60",
    "src/recursive_horizons/fgc/evolution/proto18_prelaunch_contract.py": "e7a86c5d52ce6e020ae5316038518b960d608e8467622dd8cb129a85cc5bfa21",
    "docs/fgc-pro18-frz1.md": "4fc2b5850918c3d1ab7555646304374210e7e2a6118bbfabd578d544bd3cc2fb",
    "configs/fgc/fgc-2-sf1-protocol-v17.toml": "520466fb2f798159f7d0183938327eec0fa1ccf8b49c41101924f1ff01fcd484",
    "src/recursive_horizons/fgc/evolution/protocol_v17.py": "2933d60d90607e0361d7deaec5395d8ad70a15ec41ab98ba2dd6aedbe9c5085c",
    "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py": "5fbd7f990d9f821344bbaae80840070386b3ae483fc278a85ae2a3ee5de8fa1d",
    "results/fgc-1-pro15-frz1.json": "99d24425cb5d40da2e248605265b83bcdc82a2126cfa4b5bce4239c6ef30635c",
}
TRUE_CLAIMS = {
    "PRO18_predecessor_bound",
    "PREF26_contract_frozen",
    "PREF26_read_only_source_execution_authorized",
    "PREF26_completed",
    "actual_raw_bundle_byte_and_hash_recomputed",
    "actual_legacy_history_semantically_replayed",
    "AUTH1_input_evidence_derived",
    "AUTH1_implementation_or_design_authorized",
}
CONFIG_TRUE_CLAIMS = {
    "PRO18_predecessor_bound",
    "PREF26_contract_frozen",
    "PREF26_read_only_source_execution_authorized",
}
FALSE_CLAIMS = {
    "PROTO17_genesis_spec_constructed",
    "AUTH1_authorized", "AUTH1_committed", "HLT15_GEN1_authorized",
    "PREF27_authorized", "future_output_roots_reobserved",
    "future_namespace_created", "pretrajectory_operation_authorized",
    "trajectory_read", "fresh_GR0_calibration_authorized",
    "fresh_GR0_dynamic_calibration_completed", "GR0_case_eligible",
    "candidate_execution_authorized", "SGBL_execution_authorized",
    "FGCQR_holdout_execution_authorized", "DEF1_execution_authorized",
    "retained_EFT_evolution_authorized", "physical_transition_claim_authorized",
    "global_continuation_authorized", "FGCQR_mechanism_rejected",
    "general_gradient_route_rejected", "singularity_resolution_derived",
    "child_domain_or_topology_derived", "dark_sector_mechanism_derived",
    "varying_locally_measured_c_derived",
}


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
        + "\n"
    ).encode("utf-8")


def _duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
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
        raise Proto18Pref26BinderError(
            "PREF26_INVALID_OR_NONCONVERGED_BINDER", f"{path.name} is not a mapping",
        )
    return value


def _git(*arguments: str, check: bool = True) -> bytes:
    try:
        return subprocess.run(
            ["git", *arguments], cwd=ROOT, check=check,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        raise Proto18Pref26BinderError(
            "PREF26_PREDECESSOR_LINEAGE_DRIFT", "sealed Git operation failed",
        ) from error


def _verify_sealed_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    predecessor = config.get("predecessor")
    construction = config.get("construction_inputs")
    expected_predecessor = {
        "sealed_commit": SEALED_COMMIT,
        "sealed_commit_must_be_ancestor_of_HEAD": True,
        "predecessor_artifact": "FGC-1-PRO18-FRZ1",
        "target_protocol": "FGC-2-SF1-PROTO17",
        "predecessor_config": "configs/fgc/fgc-1-pro18-frz1.toml",
        "predecessor_config_sha256": SEALED_BLOBS["configs/fgc/fgc-1-pro18-frz1.toml"],
        "predecessor_result": "results/fgc-1-pro18-frz1.json",
        "predecessor_result_sha256": SEALED_BLOBS["results/fgc-1-pro18-frz1.json"],
        "predecessor_reproducer": "scripts/reproduce_fgc_pro18_frz1.py",
        "predecessor_reproducer_sha256": SEALED_BLOBS["scripts/reproduce_fgc_pro18_frz1.py"],
        "predecessor_contract": "src/recursive_horizons/fgc/evolution/proto18_prelaunch_contract.py",
        "predecessor_contract_sha256": SEALED_BLOBS["src/recursive_horizons/fgc/evolution/proto18_prelaunch_contract.py"],
        "predecessor_document": "docs/fgc-pro18-frz1.md",
        "predecessor_document_sha256": SEALED_BLOBS["docs/fgc-pro18-frz1.md"],
        "authority_identity": "FGC-1-PRO18-AUTH1",
    }
    expected_construction = {
        "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v17.toml",
        "protocol_config_sha256": SEALED_BLOBS["configs/fgc/fgc-2-sf1-protocol-v17.toml"],
        "protocol_schema_module": "src/recursive_horizons/fgc/evolution/protocol_v17.py",
        "protocol_schema_module_sha256": SEALED_BLOBS["src/recursive_horizons/fgc/evolution/protocol_v17.py"],
        "protocol_construction_module": "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py",
        "protocol_construction_module_sha256": SEALED_BLOBS["src/recursive_horizons/fgc/evolution/proto17_pure_construction.py"],
        "restart_lineage_result": "results/fgc-1-pro15-frz1.json",
        "restart_lineage_result_sha256": SEALED_BLOBS["results/fgc-1-pro15-frz1.json"],
        "initial_TDG6_ledger_policy": "exact_zero_initialization_at_23_over_16",
        "final_PROTO17_GenesisSpec_construction_deferred_to_AUTH1": True,
    }
    if predecessor != expected_predecessor or construction != expected_construction:
        raise Proto18Pref26BinderError(
            "PREF26_PREDECESSOR_LINEAGE_DRIFT", "sealed lineage contract differs",
        )
    _git("cat-file", "-e", f"{SEALED_COMMIT}^{{commit}}")
    status = subprocess.run(
        ["git", "merge-base", "--is-ancestor", SEALED_COMMIT, "HEAD"], cwd=ROOT,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    ).returncode
    if status:
        raise Proto18Pref26BinderError(
            "PREF26_PREDECESSOR_LINEAGE_DRIFT", "sealed commit is not an ancestor",
        )
    receipts: dict[str, Any] = {}
    for relative, expected in SEALED_BLOBS.items():
        historical = _git("show", f"{SEALED_COMMIT}:{relative}")
        live = ROOT / relative
        if (sha256(historical).hexdigest() != expected or not live.is_file()
                or live.is_symlink() or _sha(live) != expected):
            raise Proto18Pref26BinderError(
                "PREF26_PREDECESSOR_LINEAGE_DRIFT", f"sealed blob differs: {relative}",
            )
        receipts[relative] = expected
    return {
        "sealed_commit": SEALED_COMMIT,
        "sealed_commit_is_ancestor_of_HEAD": True,
        "tracked_blob_sha256": receipts,
    }


def _validate_config(config: Mapping[str, Any]) -> None:
    expected_top = {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "scope", "predecessor", "construction_inputs", "source_contract",
        "resource_bounds", "semantic_replay", "future_chain", "typed_stops",
        "duty_transition_on_success", "authorization_on_success", "claims",
    }
    if set(config) != expected_top:
        raise Proto18Pref26BinderError(
            "PREF26_INVALID_OR_NONCONVERGED_BINDER", "configuration fields differ",
        )
    validate_pref26_contract(config)
    future = {
        "PREF26": "FGC-1-PRO18-PREF26", "AUTH1": "FGC-1-PRO18-AUTH1",
        "HLT15_GEN1": "FGC-1-HLT15-GEN1", "PREF27": "FGC-1-PRO18-PREF27",
        "required_order": ["PREF26", "AUTH1", "HLT15_GEN1", "PREF27"],
        "PREF26_must_not_authorize_AUTH1_or_later": True,
        "AUTH1_must_be_committed_before_HLT15_GEN1": True,
        "HLT15_GEN1_must_recheck_authority_and_import_in_one_isolated_invocation": True,
        "PREF27_owns_actual_namespace_reuse_and_foreign_store_rejection": True,
    }
    transition = {
        "raw_bundle_byte_and_hash_recomputation_required": True,
        "raw_bundle_byte_and_hash_recomputation_completed": True,
        "semantic_historical_journal_payload_replay_required": True,
        "semantic_historical_journal_payload_replay_completed": True,
        "external_Git_authority_blob_and_live_import_validation_required": True,
        "external_Git_authority_blob_and_live_import_validation_completed": False,
        "namespace_reuse_and_foreign_store_rejection_required": True,
        "namespace_reuse_and_foreign_store_rejection_completed": False,
    }
    authorization = {
        "PREF26_completed": True, "AUTH1_implementation_or_design_authorized": True,
        "AUTH1_authorized": False, "HLT15_GEN1_authorized": False,
        "PREF27_authorized": False, "fresh_GR0_calibration_authorized": False,
    }
    claims = config.get("claims")
    # The two actual-evidence fields and completion fields are all false before
    # the operation; the result promotes only the explicitly frozen success transition.
    expected_claims = {
        key: key in CONFIG_TRUE_CLAIMS
        for key in sorted(TRUE_CLAIMS | FALSE_CLAIMS)
    }
    if (config.get("future_chain") != future
            or config.get("duty_transition_on_success") != transition
            or config.get("authorization_on_success") != authorization
            or not isinstance(claims, Mapping) or dict(claims) != expected_claims
            or set(claims) != TRUE_CLAIMS | FALSE_CLAIMS):
        raise Proto18Pref26BinderError(
            "PREF26_AUTHORITY_OR_LAUNCH_PROMOTION_FORBIDDEN",
            "future chain, transition, or claim boundary differs",
        )


def reproduce() -> dict[str, Any]:
    config = _load_toml(CONFIG)
    _validate_config(config)
    lineage = _verify_sealed_lineage(config)
    pro18 = _load_toml(PRO18_CONFIG)
    evidence = bind_pref26(ROOT, config, pro18)
    true_claims = {key: True for key in sorted(TRUE_CLAIMS)}
    false_claims = {key: False for key in sorted(FALSE_CLAIMS)}
    decision = {**true_claims, **false_claims}
    implementation = {
        str(SCRIPT.relative_to(ROOT)): _sha(SCRIPT),
        **{relative: _sha(ROOT / relative) for relative in IMPLEMENTATION_PATHS},
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT,
        "project_version": "0.11.0",
        "classification": "completed_read_only_actual_two_container_restart_and_legacy_history_binding_without_authority_namespace_or_trajectory",
        "source_config_sha256": {str(CONFIG.relative_to(ROOT)): _sha(CONFIG)},
        "implementation_sha256": implementation,
        "derivation_document": str(DOC.relative_to(ROOT)),
        "derivation_document_sha256": _sha(DOC),
        "generated_by": str(SCRIPT.relative_to(ROOT)),
        "scope_bindings": {
            "target_protocol": "FGC-2-SF1-PROTO17",
            "role": "completed_read_only_actual_source_container_and_legacy_history_binder",
            "raw_source_checkpoints_opened_read_only": True,
            "historical_event_logs_opened_read_only": True,
            "future_output_roots_observed": False,
            "future_namespace_created": False,
            "pretrajectory_operation": False,
            "trajectory_read": False,
            "mechanism_question_answered": False,
        },
        "artifact_payload": {
            "decision": decision,
            "immutable_lineage": lineage,
            "source_archive_manifest": evidence.source_archive_manifest,
            "historical_replay": evidence.historical_replay,
            "AUTH1_input_evidence": evidence.auth1_input_evidence,
            "HLT14_production_duties": evidence.duty_status,
            "authority_boundary": {
                "AUTH1_identity": "FGC-1-PRO18-AUTH1",
                "AUTH1_is_separate_future_committed_artifact": True,
                "PREF26_result_is_not_launch_authority": True,
                "PROTO17_genesis_spec_constructed": False,
                "future_runtime_adapter_and_runner_pins_deferred_to_AUTH1": True,
                "external_Git_trust_is_not_cryptographic_authentication": True,
            },
        },
        "gate_status": decision,
        "nonclaims": false_claims,
    }


def _load_result(path: Path = OUTPUT) -> dict[str, Any]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Proto18Pref26BinderError(
            "PREF26_INVALID_OR_NONCONVERGED_BINDER", "stored result is malformed",
        ) from error
    if not isinstance(value, dict) or _canonical(value) != raw:
        raise Proto18Pref26BinderError(
            "PREF26_INVALID_OR_NONCONVERGED_BINDER", "stored result is noncanonical",
        )
    return value


def _atomic_write(payload: bytes) -> None:
    """Atomically replace only PREF26's canonical compact-result path."""
    path = OUTPUT
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent,
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    if not arguments.write and not arguments.verify:
        arguments.verify = True
    try:
        result = reproduce()
        encoded = _canonical(result)
        if arguments.write:
            _atomic_write(encoded)
            print(f"wrote {OUTPUT}")
        if arguments.verify:
            if _load_result() != result:
                raise Proto18Pref26BinderError(
                    "PREF26_INVALID_OR_NONCONVERGED_BINDER",
                    "stored result differs from complete read-only reproduction",
                )
            print("FGC-1-PRO18-PREF26 verification passed")
        return 0
    except (Proto18Pref26BinderError, OSError, ValueError) as error:
        print(f"PREF26 failed closed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
