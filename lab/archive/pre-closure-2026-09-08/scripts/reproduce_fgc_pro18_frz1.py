#!/usr/bin/env python3
"""Reproduce the compact-only PRO18 production-prelaunch freeze."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
from recursive_horizons.fgc.evolution.proto18_prelaunch_contract import (  # noqa: E402
    PRO18_FALSE_CLAIMS, PRO18_TRUE_CLAIMS, Proto18PrelaunchStop, canonical_bytes,
    duplicate_safe_json, inspect_future_output_roots_only, require_exact_claims,
    verify_immutable_compact_files,
)

ARTIFACT = "FGC-1-PRO18-FRZ1"
CONFIG = ROOT / "configs/fgc/fgc-1-pro18-frz1.toml"
OUTPUT = ROOT / "results/fgc-1-pro18-frz1.json"
DOC = ROOT / "docs/fgc-pro18-frz1.md"
MODULE = ROOT / "src/recursive_horizons/fgc/evolution/proto18_prelaunch_contract.py"
EXPECTED_ORDER = ["RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"]
EXPECTED_IDENTITIES = [
    ("RK4-2049", "PROTO12", 0, "4bb04ae97c09f66881084b5816d14ce300c54152ae11c182f53b64d5e887c072", "3bfc9a137f20353e9ead2f742f7780b049a11ab576f6e4de98fe4e5d0ad5c0d3", 342, 1710),
    ("RK4-4097", "PROTO12", 0, "cc6d3b6fbb3a5b6abd93f69f7f59ba67c38e9b7d332fac991a88847970b57aac", "ffbde443e99535e74900b1e04c1e60cc27da105dd5d605286e77a6c77620c0bb", 673, 3365),
    ("RK4-8193", "PROTO12", 0, "b95ea0be0779f712fb0c354f9c10f28b9ccc6c737fb403d5710dbd07e0e506f9", "076bebfc71482793017d2d306001b3416f5ad7a4088f35351d4d26322a09d12c", 976, 4880),
    ("SSPRK3-4097", "PROTO12", 0, "dc056db806e5d983852135138dd48dbb9b995a8585f56900258bece7e16e7016", "3cc3dc3a262f83cc19f07c49d477f5e43309f0c9e0f16bdeaf6d0c7dc02bd48f", 631, 2524),
    ("SSPRK3-8193", "PROTO12", 0, "b262bf4bd180717d42a7fa4d502f006bac0950669bc9e86c474ad6e4a21ce7ad", "d4bd39d36dfe749cc6ec04e5f112fc5544e3defdce67315ee871d2d8cfc75b5f", 987, 3948),
    ("SSPRK3-16385", "RSP2", 1, "62784c8cc63eedc210366f146f28d535343b966e1f988a4298894e78ec560884", "2946239e5d23b4e7a3cd74419a0b9c214eb89195f72bc84d4b89f39b6ef01963", 1771, 7084),
]
EXPECTED_COMPACT_PATHS = {
    "configs/fgc/fgc-1-hlt14-mon14.toml",
    "results/fgc-1-hlt14-mon14.json",
    "configs/fgc/fgc-2-sf1-protocol-v17.toml",
    "configs/fgc/fgc-1-pro17-frz1.toml",
    "results/fgc-1-pro17-frz1.json",
    "configs/fgc/fgc-1-pro17-pref25.toml",
    "results/fgc-1-pro17-pref25.json",
    "configs/fgc/fgc-1-pro15-frz1.toml",
    "src/recursive_horizons/fgc/evolution/proto15_runtime.py",
    "configs/fgc/fgc-1-cal9-pref13.toml",
    "results/fgc-1-cal9-pref13.json",
    "configs/fgc/fgc-1-rsp2-frz1.toml",
    "results/fgc-1-rsp2-frz1.json",
    "configs/fgc/fgc-1-rsp2-pref14.toml",
    "results/fgc-1-rsp2-pref14.json",
    "configs/fgc/fgc-1-pro13-frz1.toml",
    "results/fgc-1-pro13-frz1.json",
    "results/fgc-1-pro15-frz1.json",
}
EXPECTED_SOURCE_MAP = {
    "container_format": "NPZ",
    "source_container_count": 2,
    "source_event_log_name": "events.jsonl",
    "raw_archive_or_event_payload_must_not_be_opened_by_PRO18": True,
    "shared_containers": [
        "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
        "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
    ],
    "legacy_event_logs": [
        "runs/fgc-2-sf1/proto12/calibration/events.jsonl",
        "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/events.jsonl",
    ],
    "source_container_grammar": (
        "exactly_two_shared_NPZ_containers_plus_matching_legacy_events_jsonl; "
        "six selectors may share a container but never a selector identity"
    ),
}
EXPECTED_FUTURE_ROOTS = {
    "calibration": "runs/fgc-2-sf1/proto17/calibration",
    "holdout": "runs/fgc-2-sf1/proto17/holdout",
    "temporal_absence_observation_at_freeze": True,
    "only_output_root_existence_may_be_inspected": True,
}
EXPECTED_FUTURE_CHAIN = {
    "PREF26": "FGC-1-PRO18-PREF26",
    "AUTH1": "FGC-1-PRO18-AUTH1",
    "HLT15_GEN1": "FGC-1-HLT15-GEN1",
    "PREF27": "FGC-1-PRO18-PREF27",
    "required_order": ["PREF26", "AUTH1", "HLT15_GEN1", "PREF27"],
}
EXPECTED_TYPED_STOPS = {
    "raw_archive_access": "PRO18_RAW_ARCHIVE_ACCESS_FORBIDDEN",
    "historical_event_access": "PRO18_HISTORICAL_EVENT_ACCESS_FORBIDDEN",
    "output_root_present": "PRO18_FUTURE_NAMESPACE_ALREADY_PRESENT",
    "immutable_blob_drift": "PRO18_IMMUTABLE_COMPACT_BLOB_DRIFT",
    "selector_drift": "PRO18_SELECTOR_IDENTITY_DRIFT",
    "authority_or_image_drift": "PRO18_AUTHORITY_OR_DETACHED_IMAGE_DRIFT",
    "unsafe_output_root": "PRO18_UNSAFE_OUTPUT_ROOT",
}
EXPECTED_AUTHORIZATION = {
    "PREF26_implementation_or_design_authorized": True,
    "AUTH1_authorized": False,
    "HLT15_GEN1_authorized": False,
    "PREF27_authorized": False,
}
EXPECTED_DUTIES = {
    name: {"required": True, "completed": False}
    for name in (
        "raw_bundle_byte_and_hash_recomputation",
        "external_Git_authority_blob_and_live_import_validation",
        "semantic_historical_journal_payload_replay",
        "namespace_reuse_and_foreign_store_rejection",
    )
}

def sha(path: Path) -> str: return sha256(path.read_bytes()).hexdigest()

def load_config(path: Path = CONFIG) -> dict:
    raw = path.read_bytes(); value = tomllib.loads(raw.decode("utf-8"))
    if set(value) != {"schema_version", "artifact_id", "project_version", "target_protocol", "scope", "immutable_lineage", "source_map", "selectors", "future_output_roots", "future_chain", "typed_stops", "claims", "authorization", "HLT14_production_duties"}:
        raise Proto18PrelaunchStop("PRO18 config schema differs")
    if value["schema_version"] != 1 or value["project_version"] != "0.11.0" or value["artifact_id"] != ARTIFACT or value["target_protocol"] != "FGC-2-SF1-PROTO17":
        raise Proto18PrelaunchStop("PRO18 identity differs")
    if value["scope"] != {"role": "static_production_prelaunch_overlay_without_raw_history_or_namespace_access", "raw_archive_opened": False, "historical_event_log_opened": False, "future_namespace_created": False, "trajectory_read": False, "mechanism_question_answered": False}:
        raise Proto18PrelaunchStop("PRO18 compact scope differs")
    lineage = value["immutable_lineage"]
    expected_lineage_keys = {
        "checkpoint_commit", "checkpoint_must_be_ancestor_of_HEAD", "authority_identity",
        "detached_launch_image_required", "isolated_detached_launch_image_required",
        "same_invocation_authority_and_import_recheck_required",
        "external_Git_trust_is_not_cryptographically_authenticated",
        "TOCTOU_fully_eliminated", "authority_recheck_required_immediately_before_atomic_genesis",
        "compact_files",
    }
    if set(lineage) != expected_lineage_keys:
        raise Proto18PrelaunchStop("PRO18 immutable-lineage schema differs")
    expected_lineage_values = {
        "checkpoint_commit": "4b27a48facf75b1c3186efed6cdfa752ce160e93",
        "checkpoint_must_be_ancestor_of_HEAD": True,
        "authority_identity": "FGC-1-PRO18-AUTH1",
        "detached_launch_image_required": True,
        "isolated_detached_launch_image_required": True,
        "same_invocation_authority_and_import_recheck_required": True,
        "external_Git_trust_is_not_cryptographically_authenticated": True,
        "TOCTOU_fully_eliminated": False,
        "authority_recheck_required_immediately_before_atomic_genesis": True,
    }
    if any(lineage[key] != expected for key, expected in expected_lineage_values.items()):
        raise Proto18PrelaunchStop("PRO18 immutable authority contract differs")
    if set(lineage["compact_files"]) != EXPECTED_COMPACT_PATHS:
        raise Proto18PrelaunchStop("PRO18 compact-lineage file set differs")
    if any(
        not isinstance(digest, str)
        or len(digest) != 64
        or digest.lower() != digest
        or any(character not in "0123456789abcdef" for character in digest)
        for digest in lineage["compact_files"].values()
    ):
        raise Proto18PrelaunchStop("PRO18 compact-lineage digest is malformed")
    if value["source_map"] != EXPECTED_SOURCE_MAP:
        raise Proto18PrelaunchStop("PRO18 source map differs")
    if value["future_output_roots"] != EXPECTED_FUTURE_ROOTS:
        raise Proto18PrelaunchStop("PRO18 future output-root contract differs")
    if value["future_chain"] != EXPECTED_FUTURE_CHAIN:
        raise Proto18PrelaunchStop("PRO18 future chain ordering differs")
    if value["typed_stops"] != EXPECTED_TYPED_STOPS:
        raise Proto18PrelaunchStop("PRO18 typed-stop contract differs")
    if set(value["selectors"]) != {
        "canonical_order", "restart_coordinate_time", "accepted_boundary_time",
        "active_target_time", "member",
    }:
        raise Proto18PrelaunchStop("PRO18 selector-root schema differs")
    members = value["selectors"].get("member")
    if not isinstance(members, list) or len(members) != 6:
        raise Proto18PrelaunchStop("exactly six selectors required")
    order = value["selectors"]["canonical_order"]
    if [m.get("key") for m in members] != order or len(set(order)) != 6:
        raise Proto18PrelaunchStop("selector order or uniqueness differs")
    if order != EXPECTED_ORDER or value["selectors"].get("restart_coordinate_time") != "23/16" or value["selectors"].get("accepted_boundary_time") != "23/16" or value["selectors"].get("active_target_time") != "3/2":
        raise Proto18PrelaunchStop("selector time identity differs")
    for item in members:
        if set(item) != {"key", "source", "archive_index", "method", "point_count", "array_prefix", "state_sha256", "physical_state_sha256", "restart_sha256", "input_sha256", "coordinate_time", "source_retry_count", "CFL_retry_count", "accepted_stage_count", "event_sample_count", "tracer_count", "state_shape", "step_index", "transaction_serial"}:
            raise Proto18PrelaunchStop("selector schema differs")
        if item["archive_index"] not in (0, 1) or item["step_index"] <= 0 or item["transaction_serial"] <= 0:
            raise Proto18PrelaunchStop("selector identity differs")
    found = [(m["key"], m["source"], m["archive_index"], m["restart_sha256"], m["input_sha256"], m["step_index"], m["transaction_serial"]) for m in members]
    if found != EXPECTED_IDENTITIES:
        raise Proto18PrelaunchStop("PRO18_SELECTOR_IDENTITY_DRIFT")
    if value["authorization"] != EXPECTED_AUTHORIZATION:
        raise Proto18PrelaunchStop("PRO18 authorization map differs")
    if value["HLT14_production_duties"] != EXPECTED_DUTIES:
        raise Proto18PrelaunchStop("HLT14 production duties differ")
    require_exact_claims(value["claims"])
    _cross_check_selectors_from_sealed_compact(value)
    return value

def _cross_check_selectors_from_sealed_compact(config: dict) -> None:
    """Cross-check selector detail from sealed compact PRO13/PROTO15 records only."""
    def parse(path: str) -> list[dict]:
        value = json.loads((ROOT / path).read_text())
        return value["artifact_payload"]["restart_members"]
    pro13 = parse("results/fgc-1-pro13-frz1.json")
    pro15 = parse("results/fgc-1-pro15-frz1.json")
    fields = ("key", "method", "point_count", "state_sha256", "restart_payload_sha256", "input_hash", "coordinate_time", "source_retry_count", "CFL_retry_count", "accepted_stage_count", "event_sample_count", "tracer_count", "state_shape", "step_index", "transaction_serial")
    for selector, thirteen, fifteen in zip(config["selectors"]["member"], pro13, pro15, strict=True):
        if any(thirteen[name] != fifteen[name] for name in fields):
            raise Proto18PrelaunchStop("sealed PRO13/PROTO15 selector records disagree")
        expected = {"source": thirteen["source_checkpoint"], "method": thirteen["method"], "point_count": thirteen["point_count"], "state_sha256": thirteen["state_sha256"], "physical_state_sha256": thirteen["state_sha256"], "restart_sha256": thirteen["restart_payload_sha256"], "input_sha256": thirteen["input_hash"], "coordinate_time": "23/16", "source_retry_count": thirteen["source_retry_count"], "CFL_retry_count": thirteen["CFL_retry_count"], "accepted_stage_count": thirteen["accepted_stage_count"], "event_sample_count": thirteen["event_sample_count"], "tracer_count": thirteen["tracer_count"], "state_shape": thirteen["state_shape"], "step_index": thirteen["step_index"], "transaction_serial": thirteen["transaction_serial"], "array_prefix": thirteen["key"].replace("-", "_")}
        if any(selector[name] != expected[name] for name in expected):
            raise Proto18PrelaunchStop("PRO18_SELECTOR_IDENTITY_DRIFT")

def build(*, observe_prelaunch_roots: bool) -> dict:
    config = load_config()
    lineage = config["immutable_lineage"]
    hashes = verify_immutable_compact_files(ROOT, lineage["checkpoint_commit"], lineage["compact_files"])
    root_state = None
    if observe_prelaunch_roots:
        root_state = inspect_future_output_roots_only(ROOT, config["future_output_roots"])
        if any(root_state.values()): raise Proto18PrelaunchStop("PRO18_FUTURE_NAMESPACE_ALREADY_PRESENT")
    payload = {
        "scope": config["scope"], "immutable_lineage": {"checkpoint_commit": lineage["checkpoint_commit"], "authority_identity": lineage["authority_identity"], "compact_files": hashes, "isolated_detached_launch_image_required": lineage["isolated_detached_launch_image_required"], "same_invocation_authority_and_import_recheck_required": lineage["same_invocation_authority_and_import_recheck_required"], "external_Git_trust_is_not_cryptographically_authenticated": lineage["external_Git_trust_is_not_cryptographically_authenticated"], "TOCTOU_fully_eliminated": lineage["TOCTOU_fully_eliminated"]},
        "source_map": config["source_map"], "selectors": config["selectors"],
        "future_output_roots": {
            "recorded_absent_at_freeze": {
                name: True for name in (root_state or {"calibration": False, "holdout": False})
            },
            "absence_was_observed_when_result_was_written": True,
            "ordinary_verification_does_not_reobserve_current_roots": True,
            "temporal_observation_only": True,
        },
        "future_chain": config["future_chain"], "typed_stops": config["typed_stops"], "claims": config["claims"], "authorization": config["authorization"], "HLT14_production_duties": config["HLT14_production_duties"],
    }
    return {
        "schema_version": 1, "artifact_id": ARTIFACT, "project_version": config["project_version"],
        "classification": "static_compact_only_production_prelaunch_overlay_without_raw_history_or_namespace_access",
        "generated_by": "scripts/reproduce_fgc_pro18_frz1.py", "source_config_sha256": {str(CONFIG.relative_to(ROOT)): sha(CONFIG)},
        "implementation_sha256": {str(MODULE.relative_to(ROOT)): sha(MODULE), str(Path(__file__).relative_to(ROOT)): sha(Path(__file__))},
        "derivation_document": str(DOC.relative_to(ROOT)), "derivation_document_sha256": sha(DOC),
        "scope_bindings": config["scope"], "gate_status": config["claims"],
        "nonclaims": {name: False for name in sorted(PRO18_FALSE_CLAIMS)},
        "artifact_payload": payload,
    }

def verify_result(result: dict) -> None:
    expected = build(observe_prelaunch_roots=False)
    if result != expected:
        raise Proto18PrelaunchStop("canonical PRO18 result differs from current compact construction")

def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true"); parser.add_argument("--verify", action="store_true"); parser.add_argument("--verify-prelaunch", action="store_true")
    args = parser.parse_args()
    if args.verify_prelaunch:
        build(observe_prelaunch_roots=True)
    result = build(observe_prelaunch_roots=args.write)
    if args.write: OUTPUT.write_bytes(canonical_bytes(result) + b"\n")
    if args.verify:
        stored = duplicate_safe_json(OUTPUT.read_bytes().rstrip(b"\n"), "PRO18 result")
        verify_result(stored)
    if not args.write and not args.verify and not args.verify_prelaunch: print(json.dumps(result, indent=2, sort_keys=True))
    return 0
if __name__ == "__main__": raise SystemExit(main())
