"""Pure compact qualification record for the HLT16 durable-runtime layer.

This builder inventories tracked implementation surfaces only.  It has no
array decoder, store path, runner entrypoint, or numerical dependency: MON16
is deliberately unable to inspect or advance an experiment.
"""
from __future__ import annotations

from hashlib import sha256
import os
from pathlib import Path
import stat
from typing import Any, Mapping


ARTIFACT_ID = "FGC-1-HLT16-MON16"
PENDING_SHA256 = "PENDING_COORDINATOR_REFRESH"
PENDING_CLASSIFICATION = "implemented_qualified_runtime_inventory_pending_integration_refresh"
SEALED_CLASSIFICATION = "sealed_durable_runtime_qualified_nonexecuting"
_ROLES = frozenset({"source", "test", "status", "runner", "manifest"})
_CLAIMS = {
    "durable_runtime_implemented": True,
    "durable_runtime_qualified": True,
    "one_event_preflight_authorized": False,
    "state_advance_authorized": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_transition_claim_authorized": False,
}


class HLT16Mon16ArtifactError(ValueError):
    """The compact MON16 contract or its source inventory is malformed."""


def _sha(value: object, label: str, *, allow_pending: bool) -> str:
    if allow_pending and value == PENDING_SHA256:
        return PENDING_SHA256
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise HLT16Mon16ArtifactError(f"{label} differs")
    try:
        int(value, 16)
    except ValueError as error:
        raise HLT16Mon16ArtifactError(f"{label} differs") from error
    return value


def _relative_path(value: object) -> str:
    if not isinstance(value, str) or not value or value.startswith("/"):
        raise HLT16Mon16ArtifactError("inventory path differs")
    path = Path(value)
    if ".." in path.parts or path.parts[0] not in {"src", "tests", "scripts", "configs", "results"}:
        raise HLT16Mon16ArtifactError("inventory path differs")
    # MON16 cannot reach mutable payloads, arrays, or run namespaces even if a
    # caller tried to smuggle one in through the inventory.
    if "runs" in path.parts or path.suffix in {".npy", ".npz"}:
        raise HLT16Mon16ArtifactError("inventory attempts to bind runtime payloads")
    return path.as_posix()


def _nofollow_regular_bytes(repository: Path, relative: str) -> bytes:
    """Read one unique regular source file without trusting path traversal.

    MON16 binds source text, so a source-tree symlink or a raced replacement is
    an authority failure rather than an alternate way to obtain the same bytes.
    This implementation is deliberately local rather than shared with launch
    authority: the two gates must not inherit one traversal defect.
    """
    parts = Path(relative).parts
    directory_fd = -1
    leaf_fd = -1
    try:
        directory_fd = os.open(
            repository,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        root_stat = os.fstat(directory_fd)
        if not stat.S_ISDIR(root_stat.st_mode) or root_stat.st_nlink < 1:
            raise HLT16Mon16ArtifactError("repository root is unsafe")
        for component in parts[:-1]:
            child_fd = os.open(
                component,
                os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=directory_fd,
            )
            child_stat = os.fstat(child_fd)
            if not stat.S_ISDIR(child_stat.st_mode):
                os.close(child_fd)
                raise HLT16Mon16ArtifactError("required inventory parent is unsafe")
            os.close(directory_fd)
            directory_fd = child_fd
        leaf_fd = os.open(
            parts[-1],
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
        before = os.fstat(leaf_fd)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise HLT16Mon16ArtifactError("required inventory path is absent or unsafe")
        chunks: list[bytes] = []
        while True:
            block = os.read(leaf_fd, 1024 * 1024)
            if not block:
                break
            chunks.append(block)
        after = os.fstat(leaf_fd)
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
        ):
            raise HLT16Mon16ArtifactError("required inventory path changed during read")
        return b"".join(chunks)
    except FileNotFoundError as error:
        raise HLT16Mon16ArtifactError("required inventory path is absent or unsafe") from error
    except OSError as error:
        raise HLT16Mon16ArtifactError("required inventory path is absent or unsafe") from error
    finally:
        if leaf_fd != -1:
            os.close(leaf_fd)
        if directory_fd != -1:
            os.close(directory_fd)


def validate_mon16_config(config: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(config, Mapping) or {
        "schema_version", "artifact_id", "target_protocol", "predecessors",
        "auth1_environment", "claims", "inventory",
    } != set(config):
        raise HLT16Mon16ArtifactError("MON16 config fields differ")
    if (config["schema_version"], config["artifact_id"], config["target_protocol"]) != (
        1, ARTIFACT_ID, "FGC-2-SF1-PROTO18",
    ):
        raise HLT16Mon16ArtifactError("MON16 identity differs")
    claims = config["claims"]
    if not isinstance(claims, Mapping) or dict(claims) != _CLAIMS:
        raise HLT16Mon16ArtifactError("MON16 claim boundary differs")
    predecessors = config["predecessors"]
    if not isinstance(predecessors, list) or len(predecessors) != 2:
        raise HLT16Mon16ArtifactError("MON16 predecessor inventory differs")
    observed_predecessors: list[dict[str, str]] = []
    required_paths = {"configs/fgc/fgc-1-pro19-frz1.toml", "results/fgc-1-pro18-pref27.json"}
    for item in predecessors:
        if not isinstance(item, Mapping) or set(item) != {"artifact_id", "path", "sha256"}:
            raise HLT16Mon16ArtifactError("MON16 predecessor differs")
        path = _relative_path(item["path"])
        if not isinstance(item["artifact_id"], str) or not item["artifact_id"]:
            raise HLT16Mon16ArtifactError("MON16 predecessor differs")
        observed_predecessors.append({
            "artifact_id": item["artifact_id"], "path": path,
            "sha256": _sha(item["sha256"], "predecessor hash", allow_pending=False),
        })
    if {item["path"] for item in observed_predecessors} != required_paths:
        raise HLT16Mon16ArtifactError("MON16 predecessor paths differ")
    environment = config["auth1_environment"]
    if not isinstance(environment, Mapping) or set(environment) != {
        "python_version", "numpy_version", "system", "machine",
    } or any(not isinstance(value, str) or not value for value in environment.values()):
        raise HLT16Mon16ArtifactError("AUTH1 environment contract differs")
    inventory = config["inventory"]
    if not isinstance(inventory, list) or not inventory:
        raise HLT16Mon16ArtifactError("MON16 inventory differs")
    normalized: list[dict[str, str]] = []
    seen: set[str] = set()
    for item in inventory:
        if not isinstance(item, Mapping) or set(item) != {"role", "path", "sha256"}:
            raise HLT16Mon16ArtifactError("MON16 inventory entry differs")
        role, path = item["role"], _relative_path(item["path"])
        if role not in _ROLES or path in seen:
            raise HLT16Mon16ArtifactError("MON16 inventory role/path differs")
        seen.add(path)
        normalized.append({"role": role, "path": path,
                           "sha256": _sha(item["sha256"], "inventory hash", allow_pending=True)})
    if not {item["role"] for item in normalized}.issuperset(_ROLES):
        raise HLT16Mon16ArtifactError("MON16 inventory omits a required role")
    return {
        "predecessors": sorted(observed_predecessors, key=lambda item: item["path"]),
        "auth1_environment": dict(environment),
        "inventory": sorted(normalized, key=lambda item: item["path"]),
    }


def _observed_binding(repository: Path, item: Mapping[str, str], *, allow_pending: bool) -> dict[str, str]:
    observed = sha256(_nofollow_regular_bytes(repository, item["path"])).hexdigest()
    configured = item["sha256"]
    if configured != PENDING_SHA256 and observed != configured:
        raise HLT16Mon16ArtifactError(f"configured hash differs: {item['path']}")
    if configured == PENDING_SHA256 and not allow_pending:
        raise HLT16Mon16ArtifactError(f"predecessor binding is not sealed: {item['path']}")
    return {
        **dict(item), "observed_sha256": observed,
        "binding_status": "pending_coordinator_refresh" if configured == PENDING_SHA256 else "sealed",
    }


def build_mon16(repository: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    """Build a compact, read-only MON16 record from source-text bindings."""
    checked = validate_mon16_config(config)
    root = Path(repository)
    predecessors = [_observed_binding(root, item, allow_pending=False) for item in checked["predecessors"]]
    inventory = [_observed_binding(root, item, allow_pending=True) for item in checked["inventory"]]
    pending = [item["path"] for item in inventory if item["binding_status"] != "sealed"]
    return {
        "artifact_id": ARTIFACT_ID,
        "classification": SEALED_CLASSIFICATION if not pending else PENDING_CLASSIFICATION,
        "artifact_payload": {
            "predecessors": predecessors,
            "auth1_environment": checked["auth1_environment"],
            "inventory": inventory,
            "coordinator_refresh_required_paths": pending,
            "claims": dict(_CLAIMS),
            "explicit_nonexecution": {
                "arrays_opened": False,
                "live_run_store_opened": False,
                "runner_invoked": False,
                "one_event_preflight_completed": False,
            },
        },
    }


__all__ = [
    "ARTIFACT_ID", "HLT16Mon16ArtifactError", "PENDING_CLASSIFICATION",
    "PENDING_SHA256", "SEALED_CLASSIFICATION", "build_mon16",
    "validate_mon16_config",
]
