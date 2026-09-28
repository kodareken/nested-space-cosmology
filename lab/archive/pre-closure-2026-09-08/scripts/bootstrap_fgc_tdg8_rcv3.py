#!/usr/bin/env python3
"""Install the exact FRZ1-authorized RCV3 generation-nine projection.

Only the standard library is imported until the tracked execution image has
been authenticated. This command performs no evolution: after compact FRZ1
validation it invokes only the fixed production prefix installer.
"""

from __future__ import annotations

from hashlib import sha256
import importlib
import json
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
import sys
import tomllib
from typing import Any, Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_ID = "FGC-1-TDG8-RCV3-FRZ1"
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-rcv3-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg8-rcv3-frz1.json"
FREEZE_MODULE_PATH = "src/recursive_horizons/fgc/evolution/tdg8_rcv3_freeze.py"
RUNTIME_PATH = "src/recursive_horizons/fgc/evolution/tdg8_rcv3_fork_runtime.py"
RUNTIME_SHA256 = "9c7e40dd76d126fa4ce708025483d0c15a65e9f36a80c5919335bd53c73dca36"
SCRIPT_PATH = "scripts/bootstrap_fgc_tdg8_rcv3.py"
PRODUCTION_INSTALL_FUNCTION = "tdg8_rcv3_fork_runtime.install_recovery_fork"
PRODUCTION_SPEC = "tdg8_rcv3_fork_runtime.PRODUCTION_SPEC"


class TDG8RCV3BootstrapError(RuntimeError):
    """The compact authority, execution image, Git state, or destination differs."""


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _safe_parts(relative: str) -> tuple[str, ...]:
    path = PurePosixPath(relative)
    if (
        not relative
        or path.is_absolute()
        or str(path) != relative
        or any(part in ("", ".", "..") for part in path.parts)
    ):
        raise TDG8RCV3BootstrapError("bootstrap path is unsafe")
    return path.parts


def _read_nofollow(repository: Path, relative: str) -> bytes:
    parts = _safe_parts(relative)
    descriptors: list[int] = []
    try:
        root = os.open(
            Path(repository),
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        descriptors.append(root)
        parent = root
        for part in parts[:-1]:
            child = os.open(
                part,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=parent,
            )
            if not stat.S_ISDIR(os.fstat(child).st_mode):
                os.close(child)
                raise TDG8RCV3BootstrapError("bootstrap parent is unsafe")
            descriptors.append(child)
            parent = child
        leaf = os.open(
            parts[-1],
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent,
        )
        descriptors.append(leaf)
        metadata = os.fstat(leaf)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise TDG8RCV3BootstrapError("bootstrap leaf is unsafe")
        chunks: list[bytes] = []
        while True:
            block = os.read(leaf, 1024 * 1024)
            if not block:
                break
            chunks.append(block)
        raw = b"".join(chunks)
        if len(raw) != metadata.st_size:
            raise TDG8RCV3BootstrapError("bootstrap leaf changed during read")
        return raw
    except TDG8RCV3BootstrapError:
        raise
    except OSError as exc:
        raise TDG8RCV3BootstrapError("bootstrap image cannot be read safely") from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _git(repository: Path, *arguments: str) -> bytes:
    process = subprocess.run(
        ["git", *arguments],
        cwd=repository,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if process.returncode != 0:
        raise TDG8RCV3BootstrapError("cannot authenticate committed bootstrap image")
    return process.stdout


def _json(raw: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in items:
            if key in value:
                raise ValueError(key)
            value[key] = item
        return value
    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG8RCV3BootstrapError(f"{label} is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8RCV3BootstrapError(f"{label} is not object-valued")
    return value


def _canonical_result(value: object) -> bytes:
    try:
        return (
            json.dumps(
                value,
                sort_keys=True,
                indent=2,
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n"
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG8RCV3BootstrapError("FRZ1 result is noncanonical") from exc


def _prevalidate_authority_bytes(raw_by_path: Mapping[str, bytes]) -> dict[str, Any]:
    try:
        config = tomllib.loads(raw_by_path[CONFIG_PATH].decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG8RCV3BootstrapError("FRZ1 config is malformed") from exc
    if (
        not isinstance(config, dict)
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("classification")
        != "premise_only_generation9_prefix_projection_freeze"
    ):
        raise TDG8RCV3BootstrapError("FRZ1 config identity differs")
    expected_bootstrap = {
        "runtime_path": RUNTIME_PATH,
        "runtime_sha256": RUNTIME_SHA256,
        "script_path": SCRIPT_PATH,
        "script_sha256": _sha(raw_by_path[SCRIPT_PATH]),
        "production_install_function": PRODUCTION_INSTALL_FUNCTION,
        "production_spec": PRODUCTION_SPEC,
        "candidate_or_evolution_entrypoints": False,
    }
    if config.get("bootstrap") != expected_bootstrap:
        raise TDG8RCV3BootstrapError("bootstrap authority contract differs")
    result = _json(raw_by_path[RESULT_PATH], "FRZ1 result")
    if (
        _canonical_result(result) != raw_by_path[RESULT_PATH]
        or result.get("artifact_id") != ARTIFACT_ID
        or result.get("gate_status") != "pass"
        or result.get("source_config_sha256") != _sha(raw_by_path[CONFIG_PATH])
        or result.get("artifact_payload", {}).get("bootstrap") != expected_bootstrap
    ):
        raise TDG8RCV3BootstrapError("FRZ1 result/config/bootstrap binding differs")
    return config


def preauthenticate_execution_image(repository: Path) -> dict[str, object]:
    """Authenticate tracked bytes and Git state before importing project code."""
    root = Path(repository)
    paths = (
        CONFIG_PATH,
        RESULT_PATH,
        FREEZE_MODULE_PATH,
        RUNTIME_PATH,
        SCRIPT_PATH,
    )
    raw_by_path = {path: _read_nofollow(root, path) for path in paths}
    for path in paths:
        _git(root, "ls-files", "--error-unmatch", "--", path)
    if _git(root, "status", "--porcelain=v1", "--untracked-files=no"):
        raise TDG8RCV3BootstrapError("tracked Git state is dirty")
    for path, raw in raw_by_path.items():
        if _git(root, "show", f"HEAD:{path}") != raw:
            raise TDG8RCV3BootstrapError("live bootstrap bytes differ from HEAD")
    if _sha(raw_by_path[RUNTIME_PATH]) != RUNTIME_SHA256:
        raise TDG8RCV3BootstrapError("bound RCV3 runtime hash differs")
    config = _prevalidate_authority_bytes(raw_by_path)
    head = _git(root, "rev-parse", "HEAD").decode("ascii").strip()
    if len(head) != 40 or head.lower() != head:
        raise TDG8RCV3BootstrapError("RCV3 authority commit identity differs")
    try:
        int(head, 16)
    except ValueError as exc:
        raise TDG8RCV3BootstrapError("RCV3 authority commit identity differs") from exc
    return {
        "head": head,
        "config": config,
        "config_raw": raw_by_path[CONFIG_PATH],
        "result_raw": raw_by_path[RESULT_PATH],
        "config_sha256": _sha(raw_by_path[CONFIG_PATH]),
        "result_sha256": _sha(raw_by_path[RESULT_PATH]),
        "freeze_module_sha256": _sha(raw_by_path[FREEZE_MODULE_PATH]),
        "runtime_sha256": RUNTIME_SHA256,
        "script_sha256": _sha(raw_by_path[SCRIPT_PATH]),
    }


def bootstrap(repository: Path = ROOT) -> dict[str, object]:
    """Validate FRZ1 and invoke only the fixed production prefix installer."""
    root = Path(repository)
    image = preauthenticate_execution_image(root)

    sys.path.insert(0, str(root / "src"))
    freeze = importlib.import_module(
        "recursive_horizons.fgc.evolution.tdg8_rcv3_freeze"
    )
    runtime = importlib.import_module(
        "recursive_horizons.fgc.evolution.tdg8_rcv3_fork_runtime"
    )
    result = freeze.validate_compact(image["config_raw"], image["result_raw"])
    config = freeze.parse_config(image["config_raw"])
    if config["bootstrap"] != freeze.expected_bootstrap_contract():
        raise TDG8RCV3BootstrapError("full FRZ1 bootstrap contract differs")

    wrapper_absent, store_absent = freeze.destination_absence_status(root)
    if wrapper_absent != store_absent:
        raise TDG8RCV3BootstrapError("RCV3 destination is partial")
    destination_mode = "absent" if wrapper_absent else "existing_requires_exact_validation"

    final_image = preauthenticate_execution_image(root)
    if any(
        final_image[key] != image[key]
        for key in (
            "head", "config_sha256", "result_sha256", "freeze_module_sha256",
            "runtime_sha256", "script_sha256",
        )
    ):
        raise TDG8RCV3BootstrapError("bootstrap image changed before dispatch")

    external_authority = runtime.RecoveryForkAuthority(
        authority_commit_sha=image["head"],
        frz1_result_sha256=image["result_sha256"],
        runtime_sha256=image["runtime_sha256"],
        bootstrap_script_sha256=image["script_sha256"],
    )
    try:
        installed = runtime.install_recovery_fork(
            root,
            runtime.PRODUCTION_SPEC,
            authority=external_authority,
        )
    except runtime.TDG8RCV3ForkRuntimeError as exc:
        raise TDG8RCV3BootstrapError("RCV3 production install failed closed") from exc

    return {
        "artifact_id": ARTIFACT_ID,
        "projection_id": installed.projection_id,
        "freeze_result_sha256": image["result_sha256"],
        "freeze_config_sha256": image["config_sha256"],
        "execution_head": image["head"],
        "runtime_sha256": image["runtime_sha256"],
        "script_sha256": image["script_sha256"],
        "destination_mode_before_call": destination_mode,
        "installed_now": installed.installed_now,
        "receipt_relative": installed.receipt_relative,
        "receipt_sha256": installed.receipt_sha256,
        "anchor_checkpoint_sha256": installed.anchor_checkpoint_sha256,
        "bounded_execution_authorized": False,
        "candidate_branch_opened": False,
        "physical_result_earned": False,
        "compact_freeze_gate": result["gate_status"],
    }


def main() -> int:
    if not sys.flags.isolated or not sys.dont_write_bytecode:
        raise TDG8RCV3BootstrapError(
            "RCV3 bootstrap requires isolated no-bytecode launch: python -I -B"
        )
    print(json.dumps(bootstrap(), sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
