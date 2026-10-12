"""Static contract validator for the prospective first PROTO18 GR-0 edge."""
from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import stat
from typing import Any, Mapping


class Proto19FreezeError(ValueError):
    pass


MEMBERS = ("RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385")
STATE_HASHES = (
    "b19d22516a4c8a517e20e651a6f7a15fd68e11e7bfc74fe4d25eba13e4d8d40d",
    "61a19b6758c13105268fcca090fc29139e8b238cb4ceb8227c7d25bfde455e15",
    "671f49ea5820c797305cc3fba1a7f4ea1b4a9d29d50e2caf9aff7722d41a7fb1",
    "1ab7f9fbd4e75a0409c700674d81bb536ce30f5ff25f354b4cd40e0bfbad4b6c",
    "bca065fd543ff8c4629584d64f8a1ec49d2ab5a653e73bfad08b1fb5ad66d25d",
    "783844bc83ce2b55c9f16b343d3f081b381159a70cd16578f04692b0273a2447",
)


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def digest(value: object) -> str:
    return sha256(canonical(value)).hexdigest()


def _require_hash(root: Path, entry: Mapping[str, Any]) -> None:
    if set(entry) != {"path", "sha256"} or not isinstance(entry["path"], str):
        raise Proto19FreezeError("bound artifact schema differs")
    relative = Path(entry["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise Proto19FreezeError("bound artifact path is unsafe")
    current = root
    for part in relative.parts[:-1]:
        current /= part
        try: info = current.lstat()
        except OSError as error: raise Proto19FreezeError("bound artifact ancestor is absent") from error
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            raise Proto19FreezeError("bound artifact ancestor is unsafe")
    path = root / relative
    try: before = path.lstat()
    except OSError as error: raise Proto19FreezeError("bound artifact is absent") from error
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise Proto19FreezeError("bound artifact is unsafe")
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        try:
            payload = os.read(descriptor, before.st_size + 1)
            after = os.fstat(descriptor)
        finally: os.close(descriptor)
    except OSError as error: raise Proto19FreezeError("bound artifact cannot be snapshotted") from error
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) or len(payload) != before.st_size or sha256(payload).hexdigest() != entry["sha256"]:
        raise Proto19FreezeError(f"bound artifact differs: {entry['path']}")


def build_freeze(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    """Validate only compact frozen inputs; no store/raw-array/run access."""
    if (config.get("schema_version"), config.get("artifact_id"), config.get("protocol_artifact_id"), config.get("frozen")) != (1, "FGC-1-PRO19-FRZ1", "FGC-2-SF1-PROTO18", True):
        raise Proto19FreezeError("PROTO18 identity differs")
    scope = config.get("scope")
    if not isinstance(scope, Mapping) or any(scope.get(key) is not False for key in (
        "arrays_inspected", "namespace_created", "state_advanced", "candidate_execution", "holdout_opened", "trajectory_read"
    )):
        raise Proto19FreezeError("premise-only scope differs")
    expected_paths = {
        "results/fgc-1-pro18-pref27.json", "configs/fgc/fgc-2-sf1-protocol-v18.toml", "configs/fgc/fgc-2-sf1-protocol-v17.toml",
        "configs/fgc/fgc-1-pro15-frz1.toml", "results/fgc-1-hlt13-mon13.json", "results/fgc-1-tdg6-imp2.json", "results/fgc-1-tdg7-imp3.json",
        "src/recursive_horizons/fgc/evolution/proto15_runtime.py", "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
        "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_runtime.py", "src/recursive_horizons/fgc/evolution/gr0_calibration.py",
        "src/recursive_horizons/fgc/evolution/calibration_runtime.py", "src/recursive_horizons/fgc/evolution/health_monitor.py", "src/recursive_horizons/fgc/evolution/numerical_engine.py",
    }
    bindings = config.get("bindings")
    if not isinstance(bindings, list) or {entry.get("path") for entry in bindings if isinstance(entry, Mapping)} != expected_paths or len(bindings) != len(expected_paths):
        raise Proto19FreezeError("exact binding inventory differs")
    for entry in bindings:
        _require_hash(root, entry)
    if config.get("sealed_foundation_commit") != "142733cc2a26407e7150962879fb6a760001ea18":
        raise Proto19FreezeError("sealed PREF27 foundation commit differs")
    facts = config.get("sealed_genesis")
    compact = config.get("compact_input_evidence")
    if not isinstance(facts, Mapping) or not isinstance(compact, Mapping):
        raise Proto19FreezeError("sealed compact facts are absent")
    if facts != {
        "campaign_id": "FGC-2-SF1-PROTO17-GR0-A3-CAL-290a65bcd6a2ba820683", "namespace": "runs/fgc-2-sf1/proto17/calibration",
        "checkpoint_sha256": "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d", "receipt_sha256": "78ae89a9f74193658003cf51f8567c4088db80fd6fb86867c2d581a0c16d5baf",
        "store_tree_sha256": "95a1bb96570216b8007e32a0cc05f4c31f420052f1d90c1b7134fa27ce958585", "generation": 0,
        "journal_tip_sha256": "0" * 64, "terminal_lock": False, "cursor_mode": "FRESH_READY",
        "zero_accepted_macro_steps": True, "zero_source_retries": True, "zero_CFL_retries": True,
        "zero_temporal_retries": True, "zero_debits": True,
    }:
        raise Proto19FreezeError("sealed genesis facts differ")
    if compact != {
        "auth1_result": {"path": "results/fgc-1-pro18-auth1.json", "sha256": "dd967568fbae0bd19e063af94ea8e2d53480c98afec5c78bd86601c20faf4efa"},
        "pref26_evidence_sha256": "64d903b4e0e95eefef95d60df4cafeadb709ab6200831fdc0b025c5313842968",
        "proto12_container_sha256": "c784d4706911029883d14e1763c2d36c5dbe0e619a4e50416cc6d666d5b7ea17",
        "rsp2_container_sha256": "000546dcf3726c7882e80b7e70b64a5d2e060e2e0f7a567e174f714432ae7f34",
        "python_version": "3.14.3", "numpy_version": "2.5.1", "system": "Darwin", "machine": "arm64",
    }:
        raise Proto19FreezeError("compact input evidence differs")
    auth = compact.get("auth1_result")
    if not isinstance(auth, Mapping):
        raise Proto19FreezeError("AUTH1 compact binding is absent")
    _require_hash(root, auth)
    pref27 = json.loads((root / "results/fgc-1-pro18-pref27.json").read_text("utf-8"))
    payload = pref27.get("artifact_payload", {})
    if (pref27.get("artifact_id") != "FGC-1-PRO18-PREF27"
            or payload.get("checkpoint_sha256") != facts.get("checkpoint_sha256")
            or payload.get("receipt_sha256") != facts.get("receipt_sha256")
            or payload.get("store_tree_sha256") != facts.get("store_tree_sha256")
            or payload.get("claims", {}).get("state_advanced") is not False):
        raise Proto19FreezeError("PREF27 compact payload differs")
    auth_payload = json.loads((root / auth["path"]).read_text("utf-8")).get("artifact_payload", {})
    genesis = auth_payload.get("genesis_spec", {})
    evidence = auth_payload.get("auth1_input_evidence", {})
    if (genesis.get("campaign_id") != facts.get("campaign_id") or genesis.get("namespace") != facts.get("namespace")
            or genesis.get("common_event_index") != 23 or genesis.get("restart_coordinate_time", {}).get("rational") != "23/16"
            or genesis.get("active_event_target_time", {}).get("rational") != "3/2"
            or auth_payload.get("auth1_input_evidence_sha256") != compact.get("pref26_evidence_sha256")
            or evidence.get("physical_input_manifest", {}).get("bundle_members", [{}, {}])[0].get("raw_file_sha256") != compact.get("proto12_container_sha256")
            or evidence.get("physical_input_manifest", {}).get("bundle_members", [{}, {}, {}, {}, {}, {}])[-1].get("raw_file_sha256") != compact.get("rsp2_container_sha256")):
        raise Proto19FreezeError("AUTH1/PREF26 compact evidence differs")
    environment = auth_payload.get("numerical_environment", {}).get("contract", {})
    if tuple(environment.get(key) for key in ("python_version", "numpy_version", "system", "machine")) != tuple(compact.get(key) for key in ("python_version", "numpy_version", "system", "machine")):
        raise Proto19FreezeError("pinned runtime environment differs")
    checkpoint = auth_payload.get("genesis_spec", {})
    if checkpoint.get("genesis_checkpoint_sha256") != facts.get("checkpoint_sha256"):
        raise Proto19FreezeError("sealed checkpoint identity differs")
    for descriptor in genesis.get("member_descriptors", []):
        state = descriptor.get("state_object", {})
        ledger = descriptor.get("initial_TDG6_ledger", {})
        if (state.get("source_retry_count"), ledger.get("accepted_macro_step_count"), ledger.get("cumulative_temporal_retry_count"), ledger.get("current_macro_step_temporal_retry_count")) != (0, 0, 0, 0) or any(value != "0x0.0p+0" for value in ledger.get("accumulated_debit_vector_hex", [])):
            raise Proto19FreezeError("nonzero genesis retry or debit state")
    edge = config.get("first_edge")
    if not isinstance(edge, Mapping) or dict(edge) != {
        "branch": "GR-0", "amplitude": "3", "start_time": "23/16", "target_time": "3/2",
        "common_event_index": 23, "event_count": 1, "members": list(MEMBERS),
        "state_hashes": list(STATE_HASHES), "candidate_forbidden": True,
        "raw_import_permitted": False, "run_permitted": False,
    }:
        raise Proto19FreezeError("first progression edge differs")
    claims = config.get("claims")
    required = {
        "PROTO18_first_edge_frozen": True, "PROTO18_runtime_implemented": False,
        "PROTO18_preflight_authorized": False, "PROTO18_state_advance_authorized": False,
        "fresh_GR0_dynamic_calibration_authorized": False, "GR0_case_eligible": False,
        "SGBL_execution_authorized": False, "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False, "ROB1_authorized": False,
        "candidate_execution_authorized": False, "physical_transition_claim_authorized": False,
    }
    if claims != required:
        raise Proto19FreezeError("claim boundary differs")
    return {
        "artifact_id": "FGC-1-PRO19-FRZ1",
        "artifact_payload": {"first_edge": dict(edge), "bindings": list(config["bindings"]), "claims": dict(claims)},
    }
