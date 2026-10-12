"""Exact-delta authority for the one prospective HLT17 SRCQ1 run.

The module is implementation-image code (commit A2). A later freeze commit B
adds only the config and updates its owner document. Validation authenticates
that exact relationship, environment, source closure, caps, resources and
absent output namespace. It never calls the physical source.
"""

from __future__ import annotations

from hashlib import sha256
import os
from pathlib import Path, PurePosixPath
import platform
import re
import subprocess
import sys
import tomllib
from typing import Mapping

import numpy._core._multiarray_umath as numpy_extension

from recursive_horizons.evidence_io import read_regular_file

from .hlt17_admission_runtime import hlt17_implementation_identity
from .hlt17_srcq1 import (
    HLT17SRCQ1AuthoritySeamError,
    MEMBER_KEYS,
    observe_environment,
    recommended_resource_ceilings,
)


ARTIFACT_ID = "FGC-1-HLT17-SRCQ1-FRZ1"
CONFIG_PATH = "configs/fgc/fgc-1-hlt17-srcq1-frz1.toml"
OWNER_PATH = "docs/fgc-hlt17-srcq1-frz1.md"
HARNESS_PATH = "src/recursive_horizons/fgc/evolution/hlt17_srcq1.py"
RUNNER_PATH = "scripts/qualify_fgc_hlt17_srcq1.py"
OUTPUT_NAMESPACE = "runs/fgc-2-sf1/hlt17-srcq1-source-origin-qualification"
ORIGIN_CAPTURE_SHA256 = "2824af06c674336d479ec53faa4269f02fb598055860ec22f2ad84602a83fe72"
ALLOWED_AUTHORITY_DELTA = (CONFIG_PATH, OWNER_PATH)
_COMMIT = re.compile(r"\A(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_SOURCE_DOMAIN = b"FGC-1-HLT17-SRCQ1-SOURCE-CLOSURE-v1\n"
_ENV = {"LC_ALL": "C", "LANG": "C", "PATH": os.environ.get("PATH", "")}


class HLT17SRCQ1AuthorityError(HLT17SRCQ1AuthoritySeamError):
    """Freeze B or the live implementation/environment differs."""


def _root(value: Path) -> Path:
    root = Path(value)
    if not root.is_absolute() or root.resolve(strict=True) != root:
        raise HLT17SRCQ1AuthorityError("repository root is not canonical")
    return root


def _commit(value: object, label: str) -> str:
    if type(value) is not str or _COMMIT.fullmatch(value) is None:
        raise HLT17SRCQ1AuthorityError(f"{label} is not an exact commit")
    return value


def _git(root: Path, *arguments: str) -> bytes:
    completed = subprocess.run(
        ["git", "--no-replace-objects", *arguments],
        cwd=root,
        env=_ENV,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise HLT17SRCQ1AuthorityError("Git authority query failed")
    return completed.stdout


def _git_text(root: Path, *arguments: str) -> str:
    try:
        return _git(root, *arguments).decode("ascii").strip()
    except UnicodeDecodeError as error:
        raise HLT17SRCQ1AuthorityError("Git authority output is not ASCII") from error


def _sha_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def environment_identity() -> dict[str, str]:
    observed = observe_environment()
    extension = Path(numpy_extension.__file__).resolve(strict=True)
    executable = Path(sys.executable).resolve(strict=True)
    return {
        **observed,
        "executable_realpath": str(executable),
        "executable_sha256": _sha_file(executable),
        "numpy_extension_realpath": str(extension),
        "numpy_extension_sha256": _sha_file(extension),
        "kernel_release": platform.release(),
    }


def _tree_paths(root: Path, commit: str) -> tuple[str, ...]:
    raw = _git(root, "ls-tree", "-r", "--name-only", "-z", commit, "--", "src/recursive_horizons")
    try:
        values = tuple(item.decode("utf-8") for item in raw.split(b"\0") if item)
    except UnicodeDecodeError as error:
        raise HLT17SRCQ1AuthorityError("source closure path is not UTF-8") from error
    python = tuple(path for path in values if path.endswith((".py", ".pyi")))
    required = (RUNNER_PATH,)
    return tuple(sorted(set(python + required)))


def source_closure_identity(root: Path, implementation_commit: str) -> dict[str, object]:
    repository = _root(root)
    commit = _commit(implementation_commit, "implementation commit")
    records: list[tuple[str, str]] = []
    digest = sha256(_SOURCE_DOMAIN)
    for path in _tree_paths(repository, commit):
        pure = PurePosixPath(path)
        if pure.is_absolute() or ".." in pure.parts:
            raise HLT17SRCQ1AuthorityError("source closure path is unsafe")
        raw = _git(repository, "show", f"{commit}:{path}")
        value = sha256(raw).hexdigest()
        records.append((path, value))
        digest.update(f"{path}\0{value}\n".encode("ascii"))
    environment = environment_identity()
    for name in ("executable_sha256", "numpy_extension_sha256"):
        digest.update(f"environment:{name}\0{environment[name]}\n".encode("ascii"))
    return {
        "sha256": digest.hexdigest(),
        "file_count": len(records),
        "files": tuple(records),
        "executable_sha256": environment["executable_sha256"],
        "numpy_extension_sha256": environment["numpy_extension_sha256"],
    }


def _parse_config(raw: bytes) -> dict[str, object]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise HLT17SRCQ1AuthorityError("SRCQ1 freeze config is malformed") from error
    expected = {
        "artifact_id", "schema_version", "implementation_commit",
        "implementation_sha256", "source_closure_sha256",
        "origin_capture_sha256", "output_namespace", "resource",
        "environment", "implementation_identity", "member_cap",
    }
    if type(value) is not dict or set(value) != expected:
        raise HLT17SRCQ1AuthorityError("SRCQ1 freeze config fields differ")
    if value["artifact_id"] != ARTIFACT_ID or value["schema_version"] != 1:
        raise HLT17SRCQ1AuthorityError("SRCQ1 freeze config identity differs")
    return value


def _cap_mapping(config: Mapping[str, object]) -> dict[str, str]:
    rows = config.get("member_cap")
    if type(rows) is not list or len(rows) != len(MEMBER_KEYS):
        raise HLT17SRCQ1AuthorityError("SRCQ1 member-cap inventory differs")
    result: dict[str, str] = {}
    for expected, row in zip(MEMBER_KEYS, rows, strict=True):
        if type(row) is not dict or set(row) != {
            "member_key", "grid_spacing_hex", "previous_speed_upper_hex",
            "cfl_maximum_hex", "requested_cap_hex",
        }:
            raise HLT17SRCQ1AuthorityError("SRCQ1 member-cap fields differ")
        if row["member_key"] != expected:
            raise HLT17SRCQ1AuthorityError("SRCQ1 member-cap order differs")
        try:
            spacing = float.fromhex(row["grid_spacing_hex"])
            speed = float.fromhex(row["previous_speed_upper_hex"])
            cfl = float.fromhex(row["cfl_maximum_hex"])
            requested = float.fromhex(row["requested_cap_hex"])
        except (TypeError, ValueError) as error:
            raise HLT17SRCQ1AuthorityError("SRCQ1 member cap is not binary64 hex") from error
        for name in ("grid_spacing_hex", "previous_speed_upper_hex", "cfl_maximum_hex", "requested_cap_hex"):
            if type(row[name]) is not str or row[name].lower() != row[name]:
                raise HLT17SRCQ1AuthorityError("SRCQ1 member cap is not canonical lowercase hex")
        if min(spacing, speed, cfl, requested) <= 0.0 or (cfl * spacing / speed).hex() != requested.hex():
            raise HLT17SRCQ1AuthorityError("SRCQ1 member-cap formula differs")
        if requested.hex() != row["requested_cap_hex"]:
            raise HLT17SRCQ1AuthorityError("SRCQ1 requested-cap encoding differs")
        result[expected] = row["requested_cap_hex"]
    return result


def validate_authority_delta(
    *, repository_root: Path, authority_commit: str, implementation_sha256: str
) -> dict[str, object]:
    root = _root(Path(repository_root))
    authority = _commit(authority_commit, "authority commit")
    if _git_text(root, "rev-parse", "HEAD") != authority:
        raise HLT17SRCQ1AuthorityError("live HEAD is not the SRCQ1 authority commit")
    status = _git(root, "status", "--porcelain=v1", "-z", "--untracked-files=all")
    if status:
        raise HLT17SRCQ1AuthorityError("SRCQ1 authority requires a clean worktree")
    parents = _git_text(root, "rev-list", "--parents", "-n", "1", authority).split()
    if len(parents) != 2:
        raise HLT17SRCQ1AuthorityError("SRCQ1 authority must have one implementation parent")
    implementation = parents[1]
    config_raw = read_regular_file(root, CONFIG_PATH)
    committed_config = _git(root, "show", f"{authority}:{CONFIG_PATH}")
    if config_raw != committed_config:
        raise HLT17SRCQ1AuthorityError("live SRCQ1 config differs from authority commit")
    config = _parse_config(config_raw)
    if config["implementation_commit"] != implementation:
        raise HLT17SRCQ1AuthorityError("SRCQ1 config implementation parent differs")
    supplied_sha = str(config["implementation_sha256"])
    if supplied_sha != implementation_sha256 or supplied_sha != sha256(
        _git(root, "show", f"{implementation}:{HARNESS_PATH}")
    ).hexdigest():
        raise HLT17SRCQ1AuthorityError("SRCQ1 harness implementation hash differs")
    delta = tuple(
        item.decode("utf-8")
        for item in _git(root, "diff", "--name-only", "-z", implementation, authority).split(b"\0")
        if item
    )
    if tuple(sorted(delta)) != tuple(sorted(ALLOWED_AUTHORITY_DELTA)):
        raise HLT17SRCQ1AuthorityError("SRCQ1 authority delta paths differ")
    source = source_closure_identity(root, implementation)
    if config["source_closure_sha256"] != source["sha256"]:
        raise HLT17SRCQ1AuthorityError("SRCQ1 source closure differs")
    if config["origin_capture_sha256"] != ORIGIN_CAPTURE_SHA256:
        raise HLT17SRCQ1AuthorityError("SRCQ1 origin capture differs")
    if config["output_namespace"] != OUTPUT_NAMESPACE:
        raise HLT17SRCQ1AuthorityError("SRCQ1 output namespace differs")
    if (root / OUTPUT_NAMESPACE).exists():
        raise HLT17SRCQ1AuthorityError("SRCQ1 output namespace already exists")
    observed_environment = environment_identity()
    if config["environment"] != observed_environment:
        raise HLT17SRCQ1AuthorityError("SRCQ1 environment differs")
    if config["implementation_identity"] != hlt17_implementation_identity():
        raise HLT17SRCQ1AuthorityError("SRCQ1 C1R1 implementation identity differs")
    resources = recommended_resource_ceilings().as_mapping()
    if config["resource"] != resources:
        raise HLT17SRCQ1AuthorityError("SRCQ1 resource ceilings differ")
    caps = _cap_mapping(config)
    return {
        "artifact_id": ARTIFACT_ID,
        "authority_commit": authority,
        "implementation_commit": implementation,
        "implementation_sha256": supplied_sha,
        "exact_delta_paths": list(ALLOWED_AUTHORITY_DELTA),
        "source_closure_sha256": source["sha256"],
        "source_closure_file_count": source["file_count"],
        "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
        "output_namespace": OUTPUT_NAMESPACE,
        "environment": observed_environment,
        "resource_ceilings": resources,
        "implementation_identity": hlt17_implementation_identity(),
        "prospective_requested_cap_hex_by_member": caps,
        "boolean_flag_accepted": False,
        "campaign_execution_authorized": False,
    }


__all__ = [
    "ALLOWED_AUTHORITY_DELTA", "ARTIFACT_ID", "CONFIG_PATH",
    "HLT17SRCQ1AuthorityError", "ORIGIN_CAPTURE_SHA256", "OUTPUT_NAMESPACE",
    "environment_identity", "source_closure_identity", "validate_authority_delta",
]
