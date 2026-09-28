"""Exact-delta authority for the one prospective HLT17 SRCQ1-REC1 run.

This is implementation A2. A later freeze commit B adds only the freeze
config, updates the owner document, and updates current-foundation Make
routing. Validation authenticates that exact relationship, environment,
source closure, caps, resources, Attempt2 identities, and absent output
namespaces. It never calls the physical source and does not import
historical SRCQ1 authority decision code.
"""

from __future__ import annotations

from hashlib import sha256
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import re
import stat
import subprocess
import sys
import tomllib
from types import MappingProxyType
from typing import Mapping, NoReturn

from recursive_horizons.evidence_io import (
    DEFAULT_MAX_BYTES,
    CommitParents,
    EvidenceIOError,
    InspectDelta,
    InspectTree,
    InspectWorktree,
    ReadBlob,
    ResolveCommit,
    canonical_json_bytes,
    git_read,
    read_regular_file,
)

from .hlt17_admission_runtime import hlt17_implementation_identity
from .hlt17_srcq1_rec1 import (
    ATTEMPT2_COMPACT_RELATIVE,
    ATTEMPT2_DERIVATION_DOCUMENT,
    ATTEMPT2_METADATA_LINKED_COMPACT_COMMIT,
    ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
    ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256,
    ATTEMPT2_RAW_LEAF_BYTES,
    ATTEMPT2_RAW_LEAF_RELATIVE,
    ATTEMPT2_RAW_LEAF_SHA256,
    ATTEMPT2_SEALED_COMPACT_COMMIT,
    ATTEMPT2_SEALED_COMPACT_SHA256,
    FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER,
    HLT17SRCQ1REC1AuthoritySeamError,
    IMPLEMENTATION_RELATIVE,
    ORIGIN_CAPTURE_SHA256,
    observe_environment,
    recommended_resource_ceilings,
)
from .pro20_origin import MEMBER_KEYS


ARTIFACT_ID = "FGC-1-HLT17-SRCQ1-REC1-FRZ1"
CONFIG_PATH = "configs/fgc/fgc-1-hlt17-srcq1-rec1-frz1.toml"
OWNER_PATH = "docs/fgc-hlt17-srcq1-rec1-frz1.md"
MAKE_PATH = "mk/current-foundation.mk"
HARNESS_PATH = IMPLEMENTATION_RELATIVE
RUNNER_PATH = "scripts/qualify_fgc_hlt17_srcq1_rec1.py"
AUTHORITY_PATH = "src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1_auth1.py"
OUTPUT_NAMESPACE = "runs/fgc-2-sf1/hlt17-srcq1-rec1-source-origin-qualification"
ATTEMPT2_RAW_DIRECTORY = (
    "runs/fgc-2-sf1/hlt17-srcq1-source-origin-qualification"
)
PRO20_EVENT_NAMESPACE = "runs/fgc-2-sf1/pro20-event1/calibration"
ALLOWED_AUTHORITY_DELTA = (CONFIG_PATH, OWNER_PATH, MAKE_PATH)
CFL_MAXIMUM_HEX = "0x1.0000000000000p-3"
_COMMIT = re.compile(r"\A(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_SHA = re.compile(r"\A[0-9a-f]{64}\Z")
_SOURCE_DOMAIN = b"FGC-1-HLT17-SRCQ1-REC1-SOURCE-CLOSURE-v1\n"
_SOURCE_TREE = "src/recursive_horizons"
_STAGING_PREFIX = ".hlt17-srcq1-rec1-source-origin-qualification.evidence-io-stage-"
_PS_MAX_BYTES = DEFAULT_MAX_BYTES
_REGULAR_GIT_MODES = frozenset({"100644", "100755"})
_FORBIDDEN_PROCESS = (
    "scripts/run_fgc_",
    "scripts/qualify_fgc_",
    "run_fgc_test_partition",
    "-m unittest",
    "unittest.main",
    "-m pytest",
    "pytest",
    "codex review",
    "review-agent",
)
ENVIRONMENT_KEYS = (
    "python_implementation",
    "python_version",
    "numpy_version",
    "blas_name",
    "blas_version",
    "system",
    "machine",
    "byteorder",
    "executable",
    "executable_realpath",
    "executable_sha256",
    "numpy_extension_realpath",
    "numpy_extension_sha256",
    "kernel_release",
)
MEMBER_CAP_INPUTS = (
    (
        "RK4-2049",
        "0x1.0000000000000p-4",
        "0x1.33345e70b5188p+0",
        "0x1.aaa90b0fb5c26p-8",
    ),
    (
        "RK4-4097",
        "0x1.0000000000000p-5",
        "0x1.3333333aecf3ep+0",
        "0x1.aaaaaa9fefc9cp-9",
    ),
    (
        "RK4-8193",
        "0x1.0000000000000p-6",
        "0x1.3333333333653p+0",
        "0x1.aaaaaaaaaa654p-10",
    ),
    (
        "SSPRK3-4097",
        "0x1.0000000000000p-5",
        "0x1.3333333fe51c4p+0",
        "0x1.aaaaaa9908e70p-9",
    ),
    (
        "SSPRK3-8193",
        "0x1.0000000000000p-6",
        "0x1.33333333349e5p+0",
        "0x1.aaaaaaaaa8b26p-10",
    ),
    (
        "SSPRK3-16385",
        "0x1.0000000000000p-7",
        "0x1.3333333333346p+0",
        "0x1.aaaaaaaaaaa91p-11",
    ),
)
SOURCE_CONFIGURATION_SHA256_BY_MEMBER = MappingProxyType(
    {
        "RK4-2049": (
            "76e0cae19347d29f3958decf115df9a20953b0d831cda904b41aa14f017e7565"
        ),
        "RK4-4097": (
            "8abbdc5978e989262855281ba37e783642adac354d85fca57aa931107dbbd860"
        ),
        "RK4-8193": (
            "55b8247845791d5aa5c5131a521f668fa830628583059fb5ea55e2fa1aa82eb5"
        ),
        "SSPRK3-4097": (
            "d39b7bbb33dd3c8ffa8985d76d9b1deb5249a725f73ed158ca46e01b105c1dd8"
        ),
        "SSPRK3-8193": (
            "25671dbf22da99f1662e4fc209d2680dc0adc8f2e58b0f3d983999d705ba5f9c"
        ),
        "SSPRK3-16385": (
            "3f97c57ec9514fa7bef695a8f967110bd1dfa8141d1be8acdbc07190d3ccae61"
        ),
    }
)
_CONFIG_TOP_LEVEL = (
    "artifact_id",
    "schema_version",
    "implementation_commit",
    "implementation_sha256",
    "runner_sha256",
    "authority_sha256",
    "source_closure_sha256",
    "origin_capture_sha256",
    "output_namespace",
    "accepted_state_advanced",
    "campaign_store_written",
    "endpoint_adopted_or_serialized",
    "pro20_namespace_absent",
    "retry_authorized",
    "campaign_execution_authorized",
    "resource",
    "environment",
    "implementation_identity",
    "attempt2",
    "source_configuration",
    "member_cap",
)
_ATTEMPT2_FIELDS = (
    "sealed_compact_commit",
    "sealed_compact_sha256",
    "metadata_linked_compact_commit",
    "metadata_linked_compact_sha256",
    "derivation_document",
    "raw_leaf_relative",
    "raw_leaf_sha256",
    "raw_directory_payload_sha256",
    "raw_leaf_bytes",
)


class HLT17SRCQ1REC1AuthorityError(HLT17SRCQ1REC1AuthoritySeamError):
    """Freeze B or the live implementation/environment differs."""


def _fail(message: str) -> NoReturn:
    raise HLT17SRCQ1REC1AuthorityError(message)


def _wrap(error: Exception, message: str) -> NoReturn:
    raise HLT17SRCQ1REC1AuthorityError(message) from error


def _root(value: Path) -> Path:
    if not isinstance(value, Path):
        value = Path(value)
    supplied = Path(os.fspath(value))
    if not supplied.is_absolute():
        _fail("repository root is not canonical")
    try:
        resolved = supplied.resolve(strict=True)
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError(
            "repository root cannot be resolved"
        ) from error
    if resolved != supplied:
        _fail("repository root traverses a symlink")
    return supplied


def _commit(value: object, label: str) -> str:
    if type(value) is not str or _COMMIT.fullmatch(value) is None:
        _fail(f"{label} is not an exact commit")
    return value


def _sha(value: object, label: str) -> str:
    if type(value) is not str or _SHA.fullmatch(value) is None:
        _fail(f"{label} is not a lowercase SHA-256 digest")
    return value


def _git_query(root: Path, operation: object) -> object:
    try:
        return git_read(root, operation)
    except EvidenceIOError as error:
        _wrap(error, "committed image cannot be inspected")


def _blob(root: Path, commit: str, path: str) -> bytes:
    raw = _git_query(root, ReadBlob(commit=commit, path=path))
    if type(raw) is not bytes:
        _fail(f"{path} blob is not bytes")
    return raw


def _sha256_hex(raw: bytes) -> str:
    if type(raw) is not bytes:
        _fail("payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _stat_identity(metadata: os.stat_result) -> tuple[int, ...]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_mode,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
        metadata.st_nlink,
    )


def _hash_software_file(path: str) -> str:
    if type(path) is not str or not path or "\0" in path or len(path) > 1024:
        _fail("software path is unsafe")
    try:
        resolved = os.path.realpath(path)
        before = os.lstat(resolved)
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError(
            "software path cannot be read"
        ) from error
    if (
        type(resolved) is not str
        or not resolved
        or "\0" in resolved
        or stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_size < 0
    ):
        _fail("software path is not a regular file")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(resolved, flags)
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError(
            "software path cannot be opened"
        ) from error
    try:
        opened = os.fstat(descriptor)
        if _stat_identity(opened) != _stat_identity(before):
            _fail("software path changed before read")
        digest = sha256()
        remaining = int(before.st_size)
        while remaining > 0:
            chunk = os.read(descriptor, min(remaining, 1 << 20))
            if not chunk:
                break
            digest.update(chunk)
            remaining -= len(chunk)
        after = os.fstat(descriptor)
        path_after = os.lstat(resolved)
        if (
            remaining != 0
            or _stat_identity(after) != _stat_identity(before)
            or _stat_identity(path_after) != _stat_identity(before)
            or not stat.S_ISREG(after.st_mode)
        ):
            _fail("software path changed during read")
        return digest.hexdigest()
    finally:
        os.close(descriptor)


def _runtime_regular_path(path: str) -> str:
    if type(path) is not str or not path or "\0" in path:
        _fail("runtime path is unsafe")
    try:
        resolved = os.path.realpath(path)
        metadata = os.lstat(resolved)
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError(
            "runtime path cannot be resolved"
        ) from error
    if (
        type(resolved) is not str
        or not resolved
        or "\0" in resolved
        or len(resolved) > 1024
        or stat.S_ISLNK(metadata.st_mode)
        or not stat.S_ISREG(metadata.st_mode)
    ):
        _fail("runtime path is not a regular file")
    return resolved


def environment_identity() -> dict[str, str]:
    """Expanded live arithmetic image, including executable and NumPy hashes."""

    observed = observe_environment()
    executable = _runtime_regular_path(sys.executable)
    import numpy._core._multiarray_umath as numpy_extension

    extension_path = getattr(numpy_extension, "__file__", None)
    extension = _runtime_regular_path(str(extension_path))
    record = {
        **observed,
        "executable_realpath": executable,
        "executable_sha256": _hash_software_file(executable),
        "numpy_extension_realpath": extension,
        "numpy_extension_sha256": _hash_software_file(extension),
        "kernel_release": platform.release(),
    }
    if tuple(record) != ENVIRONMENT_KEYS:
        _fail("expanded environment keys differ")
    for key in ENVIRONMENT_KEYS:
        value = record[key]
        if type(value) is not str or not value:
            _fail(f"environment {key} is empty")
    return record


def _tree_python_paths(root: Path, commit: str) -> tuple[str, ...]:
    tree = _git_query(root, InspectTree(commit=commit, path=_SOURCE_TREE))
    paths: list[str] = []
    seen: set[str] = set()
    for entry in tree:
        path = getattr(entry, "path", None)
        kind = getattr(entry, "kind", None)
        mode = getattr(entry, "mode", None)
        if type(path) is not str or type(kind) is not str or type(mode) is not str:
            _fail("source closure tree entry differs")
        if not path.endswith((".py", ".pyi")):
            continue
        if kind != "blob" or mode not in _REGULAR_GIT_MODES:
            _fail("source closure contains a nonregular Python path")
        pure = PurePosixPath(path)
        if pure.is_absolute() or ".." in pure.parts:
            _fail("source closure path is unsafe")
        if path in seen:
            _fail("source closure repeats a path")
        seen.add(path)
        paths.append(path)
    required = (RUNNER_PATH,)
    for path in required:
        if path not in seen:
            paths.append(path)
            seen.add(path)
    return tuple(sorted(paths))


def source_closure_identity(
    root: Path,
    implementation_commit: str,
    *,
    environment: Mapping[str, str] | None = None,
) -> dict[str, object]:
    repository = _root(root)
    commit = _commit(implementation_commit, "implementation commit")
    records: list[tuple[str, str]] = []
    digest = sha256(_SOURCE_DOMAIN)
    for path in _tree_python_paths(repository, commit):
        raw = _blob(repository, commit, path)
        value = _sha256_hex(raw)
        records.append((path, value))
        digest.update(f"{path}\0{value}\n".encode("ascii"))
    observed_environment = (
        environment_identity() if environment is None else dict(environment)
    )
    if tuple(observed_environment) != ENVIRONMENT_KEYS:
        _fail("source-closure environment keys differ")
    for name in ("executable_sha256", "numpy_extension_sha256"):
        digest.update(
            f"environment:{name}\0{_sha(observed_environment[name], name)}\n".encode(
                "ascii"
            )
        )
    return {
        "sha256": digest.hexdigest(),
        "file_count": len(records),
        "files": tuple(records),
        "executable_sha256": observed_environment["executable_sha256"],
        "numpy_extension_sha256": observed_environment["numpy_extension_sha256"],
    }


def _toml_literal(value: object) -> str:
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int and not isinstance(value, bool):
        return str(value)
    if type(value) is float:
        if not math.isfinite(value):
            _fail("config float is not finite")
        return json.dumps(value, ensure_ascii=True, allow_nan=False)
    if type(value) is str:
        return json.dumps(value, ensure_ascii=True)
    _fail("unsupported TOML scalar")


def render_config(config: Mapping[str, object]) -> bytes:
    """Deterministic TOML for the closed REC1 freeze schema; no file writes."""

    if type(config) is not dict:
        _fail("emitted config is not a mapping")
    lines: list[str] = []
    for key in _CONFIG_TOP_LEVEL:
        if key not in config:
            _fail(f"emitted config omits {key}")
        item = config[key]
        if key in {
            "resource",
            "environment",
            "implementation_identity",
            "attempt2",
            "source_configuration",
            "member_cap",
        }:
            continue
        lines.append(f"{key} = {_toml_literal(item)}")
    resource = config["resource"]
    if type(resource) is not dict:
        _fail("emitted resource table differs")
    lines.extend(("", "[resource]"))
    for key, item in resource.items():
        if key in {"max_wall_seconds_per_member", "production_enclosure_limits"}:
            continue
        lines.append(f"{key} = {_toml_literal(item)}")
    walls = resource["max_wall_seconds_per_member"]
    if type(walls) is not dict:
        _fail("emitted wall table differs")
    lines.extend(("", "[resource.max_wall_seconds_per_member]"))
    for key in MEMBER_KEYS:
        lines.append(f"{key} = {_toml_literal(walls[key])}")
    limits = resource["production_enclosure_limits"]
    if type(limits) is not list:
        _fail("emitted enclosure inventory differs")
    for row in limits:
        if type(row) is not dict:
            _fail("emitted enclosure row differs")
        lines.extend(("", "[[resource.production_enclosure_limits]]"))
        for key, item in row.items():
            lines.append(f"{key} = {_toml_literal(item)}")
    for table in ("environment", "implementation_identity", "attempt2"):
        payload = config[table]
        if type(payload) is not dict:
            _fail(f"emitted {table} table differs")
        lines.extend(("", f"[{table}]"))
        for key, item in payload.items():
            lines.append(f"{key} = {_toml_literal(item)}")
    rows = config["source_configuration"]
    if type(rows) is not list:
        _fail("emitted source-configuration inventory differs")
    for row in rows:
        if type(row) is not dict:
            _fail("emitted source-configuration row differs")
        lines.extend(("", "[[source_configuration]]"))
        for key, item in row.items():
            lines.append(f"{key} = {_toml_literal(item)}")
    caps = config["member_cap"]
    if type(caps) is not list:
        _fail("emitted member-cap inventory differs")
    for row in caps:
        if type(row) is not dict:
            _fail("emitted member-cap row differs")
        lines.extend(("", "[[member_cap]]"))
        for key, item in row.items():
            lines.append(f"{key} = {_toml_literal(item)}")
    return ("\n".join(lines) + "\n").encode("ascii")


def _cap_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key, spacing, speed, requested in MEMBER_CAP_INPUTS:
        if key not in MEMBER_KEYS:
            _fail("member-cap inventory differs")
        if FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[key] != requested:
            _fail("frozen requested-cap hex differs")
        rows.append(
            {
                "member_key": key,
                "grid_spacing_hex": spacing,
                "previous_speed_upper_hex": speed,
                "cfl_maximum_hex": CFL_MAXIMUM_HEX,
                "requested_cap_hex": requested,
            }
        )
    if tuple(row["member_key"] for row in rows) != MEMBER_KEYS:
        _fail("member-cap order differs")
    return rows


def _source_configuration_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key in MEMBER_KEYS:
        rows.append(
            {
                "member_key": key,
                "sha256": SOURCE_CONFIGURATION_SHA256_BY_MEMBER[key],
            }
        )
    return rows


def _attempt2_config() -> dict[str, object]:
    return {
        "sealed_compact_commit": ATTEMPT2_SEALED_COMPACT_COMMIT,
        "sealed_compact_sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
        "metadata_linked_compact_commit": ATTEMPT2_METADATA_LINKED_COMPACT_COMMIT,
        "metadata_linked_compact_sha256": ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
        "derivation_document": ATTEMPT2_DERIVATION_DOCUMENT,
        "raw_leaf_relative": ATTEMPT2_RAW_LEAF_RELATIVE,
        "raw_leaf_sha256": ATTEMPT2_RAW_LEAF_SHA256,
        "raw_directory_payload_sha256": ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256,
        "raw_leaf_bytes": ATTEMPT2_RAW_LEAF_BYTES,
    }


def _false_flags() -> dict[str, object]:
    return {
        "accepted_state_advanced": False,
        "campaign_store_written": False,
        "endpoint_adopted_or_serialized": False,
        "pro20_namespace_absent": True,
        "retry_authorized": False,
        "campaign_execution_authorized": False,
    }


def _live_implementation_hashes(root: Path) -> dict[str, str]:
    hashes = {
        "implementation_sha256": _sha256_hex(read_regular_file(root, HARNESS_PATH)),
        "runner_sha256": _sha256_hex(read_regular_file(root, RUNNER_PATH)),
        "authority_sha256": _sha256_hex(read_regular_file(root, AUTHORITY_PATH)),
    }
    for label, value in hashes.items():
        _sha(value, label)
    return hashes


def _committed_implementation_hashes(root: Path, commit: str) -> dict[str, str]:
    hashes = {
        "implementation_sha256": _sha256_hex(_blob(root, commit, HARNESS_PATH)),
        "runner_sha256": _sha256_hex(_blob(root, commit, RUNNER_PATH)),
        "authority_sha256": _sha256_hex(_blob(root, commit, AUTHORITY_PATH)),
    }
    for label, value in hashes.items():
        _sha(value, label)
    return hashes


def expected_config(
    *,
    implementation_commit: str,
    implementation_hashes: Mapping[str, str],
    source_closure_sha256: str,
    environment: Mapping[str, str],
) -> dict[str, object]:
    resources = recommended_resource_ceilings().as_mapping()
    if (
        resources["max_rss_bytes"] != 4 * 1024 * 1024 * 1024
        or resources["max_total_wall_seconds"] != 3600.0
        or resources["scientific_threshold_fitted"] is not False
        or resources["source_retry_cap"] != 32
        or resources["cfl_retry_cap"] != 32
    ):
        _fail("recommended resource ceilings differ")
    walls = resources["max_wall_seconds_per_member"]
    if type(walls) is not dict or any(walls[key] != 900.0 for key in MEMBER_KEYS):
        _fail("recommended per-member wall ceilings differ")
    config: dict[str, object] = {
        "artifact_id": ARTIFACT_ID,
        "schema_version": 1,
        "implementation_commit": implementation_commit,
        "implementation_sha256": implementation_hashes["implementation_sha256"],
        "runner_sha256": implementation_hashes["runner_sha256"],
        "authority_sha256": implementation_hashes["authority_sha256"],
        "source_closure_sha256": source_closure_sha256,
        "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
        "output_namespace": OUTPUT_NAMESPACE,
        **_false_flags(),
        "resource": resources,
        "environment": dict(environment),
        "implementation_identity": hlt17_implementation_identity(),
        "attempt2": _attempt2_config(),
        "source_configuration": _source_configuration_rows(),
        "member_cap": _cap_rows(),
    }
    if tuple(config) != _CONFIG_TOP_LEVEL:
        _fail("emitted config field order differs")
    return config


def emit_config_bytes(
    root: Path, *, implementation_commit: str | None = None
) -> bytes:
    """Deterministic freeze-config TOML from live files and environment.

    This does not authorize a run, write a file, or construct a source.
    """

    repository = _root(root)
    if implementation_commit is not None:
        commit = _commit(implementation_commit, "implementation commit")
    else:
        head = _git_query(repository, ResolveCommit("HEAD"))
        commit = _commit(head, "HEAD")
    try:
        environment = environment_identity()
        hashes = _live_implementation_hashes(repository)
        source = source_closure_identity(
            repository, commit, environment=environment
        )
        config = expected_config(
            implementation_commit=commit,
            implementation_hashes=hashes,
            source_closure_sha256=str(source["sha256"]),
            environment=environment,
        )
        return render_config(config)
    except EvidenceIOError as error:
        _wrap(error, "freeze config cannot be emitted")


def check_tracked_config(root: Path) -> dict[str, object]:
    """Compare a later tracked freeze config to a fresh emission. No authorize."""

    repository = _root(root)
    try:
        raw = read_regular_file(repository, CONFIG_PATH)
    except EvidenceIOError as error:
        _wrap(error, "tracked REC1 freeze config cannot be read")
    parsed = _parse_config(raw)
    expected = emit_config_bytes(
        repository,
        implementation_commit=str(parsed["implementation_commit"]),
    )
    if raw != expected:
        _fail("tracked REC1 freeze config differs from emitter")
    return {
        "config_sha256": _sha256_hex(raw),
        "matched": True,
        "authorized": False,
        "executed": False,
        "physical_source_qualification_executed": False,
        "campaign_execution_authorized": False,
    }


def _parse_config(raw: bytes) -> dict[str, object]:
    if type(raw) is not bytes:
        _fail("REC1 freeze config is not bytes")
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise HLT17SRCQ1REC1AuthorityError(
            "REC1 freeze config is malformed"
        ) from error
    if type(value) is not dict or tuple(value) != _CONFIG_TOP_LEVEL:
        _fail("REC1 freeze config fields differ")
    if value["artifact_id"] != ARTIFACT_ID or value["schema_version"] != 1:
        _fail("REC1 freeze config identity differs")
    if value["output_namespace"] != OUTPUT_NAMESPACE:
        _fail("REC1 output namespace differs")
    if value["origin_capture_sha256"] != ORIGIN_CAPTURE_SHA256:
        _fail("REC1 origin capture differs")
    flags = _false_flags()
    for name, expected in flags.items():
        if value[name] is not expected:
            _fail(f"REC1 {name} flag differs")
    return value


def _cap_mapping(config: Mapping[str, object]) -> dict[str, str]:
    rows = config.get("member_cap")
    expected_rows = _cap_rows()
    if type(rows) is not list or len(rows) != len(MEMBER_KEYS):
        _fail("REC1 member-cap inventory differs")
    result: dict[str, str] = {}
    for expected, row in zip(MEMBER_KEYS, rows, strict=True):
        if type(row) is not dict or set(row) != {
            "member_key",
            "grid_spacing_hex",
            "previous_speed_upper_hex",
            "cfl_maximum_hex",
            "requested_cap_hex",
        }:
            _fail("REC1 member-cap fields differ")
        if row["member_key"] != expected:
            _fail("REC1 member-cap order differs")
        wanted = expected_rows[len(result)]
        for name in (
            "member_key",
            "grid_spacing_hex",
            "previous_speed_upper_hex",
            "cfl_maximum_hex",
        ):
            if row[name] != wanted[name]:
                _fail("REC1 member-cap values differ")
        try:
            spacing = float.fromhex(row["grid_spacing_hex"])
            speed = float.fromhex(row["previous_speed_upper_hex"])
            cfl = float.fromhex(row["cfl_maximum_hex"])
            requested = float.fromhex(row["requested_cap_hex"])
        except (TypeError, ValueError) as error:
            raise HLT17SRCQ1REC1AuthorityError(
                "REC1 member cap is not binary64 hex"
            ) from error
        for name in (
            "grid_spacing_hex",
            "previous_speed_upper_hex",
            "cfl_maximum_hex",
            "requested_cap_hex",
        ):
            if type(row[name]) is not str or row[name].lower() != row[name]:
                _fail("REC1 member cap is not canonical lowercase hex")
        if min(spacing, speed, cfl, requested) <= 0.0:
            _fail("REC1 member-cap formula differs")
        if (cfl * spacing / speed).hex() != requested.hex():
            _fail("REC1 member-cap formula differs")
        if requested.hex() != row["requested_cap_hex"]:
            _fail("REC1 requested-cap encoding differs")
        if row["requested_cap_hex"] != FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[
            expected
        ]:
            _fail("REC1 requested cap is not the frozen SRCQ1 hex")
        result[expected] = row["requested_cap_hex"]
    return result


def _source_configuration_mapping(config: Mapping[str, object]) -> dict[str, str]:
    rows = config.get("source_configuration")
    if type(rows) is not list or len(rows) != len(MEMBER_KEYS):
        _fail("REC1 source-configuration inventory differs")
    result: dict[str, str] = {}
    for expected, row in zip(MEMBER_KEYS, rows, strict=True):
        if type(row) is not dict or set(row) != {"member_key", "sha256"}:
            _fail("REC1 source-configuration fields differ")
        if row["member_key"] != expected:
            _fail("REC1 source-configuration order differs")
        digest = _sha(row["sha256"], f"{expected} source configuration")
        if digest != SOURCE_CONFIGURATION_SHA256_BY_MEMBER[expected]:
            _fail("REC1 source-configuration mapping differs")
        result[expected] = digest
    return result


def _json_object(raw: bytes, label: str) -> dict[str, object]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise HLT17SRCQ1REC1AuthorityError(f"{label} is not JSON") from error
    if type(value) is not dict:
        _fail(f"{label} is not an object")
    return value


def _require_committed_image(root: Path, authority_commit: str) -> str:
    head = _git_query(root, ResolveCommit("HEAD"))
    if head != authority_commit:
        _fail("live HEAD is not the REC1 authority commit")
    parents = _git_query(root, CommitParents(authority_commit))
    if type(parents) is not tuple or len(parents) != 1:
        _fail("REC1 authority must have one implementation parent")
    implementation = parents[0]
    if type(implementation) is not str or _COMMIT.fullmatch(implementation) is None:
        _fail("REC1 implementation parent is not an exact commit")
    delta = _git_query(
        root, InspectDelta(commit=authority_commit, against=implementation)
    )
    paths = tuple(sorted(item.path for item in delta))
    if paths != tuple(sorted(ALLOWED_AUTHORITY_DELTA)):
        _fail("REC1 authority delta paths differ")
    worktree = _git_query(root, InspectWorktree())
    if not getattr(worktree, "clean", False):
        _fail("REC1 authority requires a clean worktree")
    for path in ALLOWED_AUTHORITY_DELTA:
        try:
            live = read_regular_file(root, path)
            committed = _blob(root, authority_commit, path)
        except EvidenceIOError as error:
            _wrap(error, f"freeze path cannot be read: {path}")
        if live != committed:
            _fail(f"live freeze path differs from authority commit: {path}")
    return implementation


def _lexists(path: Path) -> bool:
    try:
        path.lstat()
    except FileNotFoundError:
        return False
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError("namespace path cannot be read") from error
    return True


def _require_relative(relative: str, label: str) -> Path:
    value = Path(relative)
    if value.is_absolute() or any(part in {"", ".", ".."} for part in value.parts):
        _fail(f"{label} is unsafe")
    return value


def require_output_absent(root: Path) -> None:
    """Require the REC1 output namespace, staging names, and symlinks absent."""

    repository = _root(root)
    relative = _require_relative(OUTPUT_NAMESPACE, "REC1 output namespace")
    current = repository
    for part in relative.parent.parts:
        current = current / part
        if not _lexists(current):
            return
        try:
            metadata = current.lstat()
        except OSError as error:
            raise HLT17SRCQ1REC1AuthorityError(
                "output parent cannot be read"
            ) from error
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _fail("output parent is unsafe")
    target = repository / relative
    if _lexists(target):
        _fail("REC1 output namespace already exists")
    if not _lexists(current):
        return
    try:
        names = os.listdir(current)
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError("output parent cannot be listed") from error
    if any(
        name.startswith(_STAGING_PREFIX) or name.startswith(f".{relative.name}.")
        for name in names
    ):
        _fail("staging namespace already exists")


def require_pro20_namespace_absent(root: Path) -> None:
    repository = _root(root)
    relative = _require_relative(PRO20_EVENT_NAMESPACE, "PRO20 event namespace")
    target = repository / relative
    if _lexists(target):
        _fail("PRO20 event namespace already exists")
    current = repository
    for part in relative.parts:
        current = current / part
        if not _lexists(current):
            return
        try:
            metadata = current.lstat()
        except OSError as error:
            raise HLT17SRCQ1REC1AuthorityError(
                "PRO20 namespace cannot be read"
            ) from error
        if stat.S_ISLNK(metadata.st_mode):
            _fail("PRO20 event namespace is unsafe")


def _directory_payload_digest(files: Mapping[str, bytes]) -> str:
    if type(files) is not dict or len(files) != 1:
        _fail("Attempt2 raw directory must contain exactly one leaf")
    inventory = {path: _sha256_hex(payload) for path, payload in files.items()}
    return _sha256_hex(canonical_json_bytes(inventory))


def authenticate_attempt2_compacts(root: Path) -> dict[str, object]:
    """Independently bind sealed Git and metadata-linked live compact blobs."""

    repository = _root(root)
    sealed = _blob(
        repository, ATTEMPT2_SEALED_COMPACT_COMMIT, ATTEMPT2_COMPACT_RELATIVE
    )
    linked = _blob(
        repository,
        ATTEMPT2_METADATA_LINKED_COMPACT_COMMIT,
        ATTEMPT2_COMPACT_RELATIVE,
    )
    try:
        live = read_regular_file(repository, ATTEMPT2_COMPACT_RELATIVE)
    except EvidenceIOError as error:
        _wrap(error, "live Attempt2 compact cannot be read")
    if _sha256_hex(sealed) != ATTEMPT2_SEALED_COMPACT_SHA256:
        _fail("sealed Attempt2 compact SHA-256 differs")
    if _sha256_hex(linked) != ATTEMPT2_METADATA_LINKED_COMPACT_SHA256:
        _fail("metadata-linked Attempt2 compact SHA-256 differs")
    if live != linked or _sha256_hex(live) != ATTEMPT2_METADATA_LINKED_COMPACT_SHA256:
        _fail("live Attempt2 compact is not the metadata-linked blob")
    sealed_obj = _json_object(sealed, "sealed Attempt2 compact")
    linked_obj = _json_object(linked, "metadata-linked Attempt2 compact")
    if "derivation_document" in sealed_obj:
        _fail("sealed Attempt2 compact must not carry derivation_document")
    expected = dict(sealed_obj)
    expected["derivation_document"] = ATTEMPT2_DERIVATION_DOCUMENT
    if linked_obj != expected:
        _fail("metadata-linked compact is not sealed payload plus derivation_document")
    if linked_obj.get("derivation_document") != ATTEMPT2_DERIVATION_DOCUMENT:
        _fail("metadata-linked derivation_document differs")
    return {
        "sealed_historical_blob": {
            "commit": ATTEMPT2_SEALED_COMPACT_COMMIT,
            "sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
            "derivation_document": None,
        },
        "metadata_linked_blob": {
            "commit": ATTEMPT2_METADATA_LINKED_COMPACT_COMMIT,
            "sha256": ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
            "derivation_document": ATTEMPT2_DERIVATION_DOCUMENT,
        },
        "compact": linked_obj,
    }


def authenticate_attempt2_raw(root: Path) -> dict[str, object]:
    """Recompute the Attempt2 raw directory digest with exactly one no-follow leaf."""

    repository = _root(root)
    relative = _require_relative(ATTEMPT2_RAW_DIRECTORY, "Attempt2 raw directory")
    directory = repository / relative
    try:
        metadata = directory.lstat()
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError(
            "Attempt2 raw directory is absent"
        ) from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("Attempt2 raw directory is unsafe")
    flags = (
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        descriptor = os.open(directory, flags)
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError(
            "Attempt2 raw directory cannot be opened"
        ) from error
    try:
        opened_directory = os.fstat(descriptor)
        if _stat_identity(opened_directory) != _stat_identity(metadata):
            _fail("Attempt2 raw directory changed before read")
        names = tuple(os.listdir(descriptor))
        if names != ("qualification.json",):
            _fail("Attempt2 raw directory leaf count or name differs")
        leaf_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
        try:
            leaf = os.open("qualification.json", leaf_flags, dir_fd=descriptor)
        except OSError as error:
            raise HLT17SRCQ1REC1AuthorityError(
                "Attempt2 raw leaf cannot be opened"
            ) from error
        try:
            leaf_before = os.fstat(leaf)
            if (
                not stat.S_ISREG(leaf_before.st_mode)
                or leaf_before.st_size != ATTEMPT2_RAW_LEAF_BYTES
            ):
                _fail("Attempt2 raw leaf is not the expected regular file")
            chunks: list[bytes] = []
            remaining = ATTEMPT2_RAW_LEAF_BYTES
            while remaining > 0:
                chunk = os.read(leaf, min(remaining, 1 << 20))
                if not chunk:
                    break
                chunks.append(chunk)
                remaining -= len(chunk)
            raw = b"".join(chunks)
            leaf_after = os.fstat(leaf)
            if (
                remaining != 0
                or _stat_identity(leaf_after) != _stat_identity(leaf_before)
            ):
                _fail("Attempt2 raw leaf changed during read")
        finally:
            os.close(leaf)
        if tuple(os.listdir(descriptor)) != ("qualification.json",):
            _fail("Attempt2 raw directory changed during read")
        directory_after = os.fstat(descriptor)
        path_after = directory.lstat()
        if (
            _stat_identity(directory_after) != _stat_identity(opened_directory)
            or _stat_identity(path_after) != _stat_identity(metadata)
        ):
            _fail("Attempt2 raw directory changed during read")
    finally:
        os.close(descriptor)
    if len(raw) != ATTEMPT2_RAW_LEAF_BYTES:
        _fail("Attempt2 raw leaf byte length differs")
    if _sha256_hex(raw) != ATTEMPT2_RAW_LEAF_SHA256:
        _fail("Attempt2 raw leaf SHA-256 differs")
    digest = _directory_payload_digest({"qualification.json": raw})
    if digest != ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256:
        _fail("Attempt2 raw directory payload SHA-256 differs")
    payload = _json_object(raw, "Attempt2 raw terminal")
    if payload.get("accepted_state_advanced") is not False:
        _fail("Attempt2 raw accepted_state_advanced differs")
    if payload.get("endpoint_adopted") is not False:
        _fail("Attempt2 raw endpoint_adopted differs")
    if payload.get("store_published") is not False:
        _fail("Attempt2 raw store_published differs")
    if payload.get("campaign_execution_authorized") is not False:
        _fail("Attempt2 raw campaign_execution_authorized differs")
    observed = payload.get("source_configuration")
    if type(observed) is not dict:
        _fail("Attempt2 raw source-configuration mapping differs")
    expected = dict(SOURCE_CONFIGURATION_SHA256_BY_MEMBER)
    if observed != expected:
        _fail("Attempt2 raw source-configuration mapping differs")
    return {
        "raw_leaf_relative": ATTEMPT2_RAW_LEAF_RELATIVE,
        "raw_leaf_sha256": ATTEMPT2_RAW_LEAF_SHA256,
        "raw_directory_payload_sha256": ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256,
        "raw_leaf_bytes": ATTEMPT2_RAW_LEAF_BYTES,
        "source_configuration": expected,
        "leaf_count": 1,
    }


def _process_rows() -> tuple[tuple[int, int, str], ...]:
    from recursive_horizons.evidence_io import _communicate_bounded

    try:
        proc = subprocess.Popen(
            ["/bin/ps", "-axo", "pid=,ppid=,command="],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            start_new_session=True,
        )
    except OSError as error:
        raise HLT17SRCQ1REC1AuthorityError(
            "process table cannot be observed"
        ) from error
    assert proc.stdout is not None
    try:
        raw = _communicate_bounded(proc, max_bytes=_PS_MAX_BYTES, timeout=5)
        code = proc.returncode
    except Exception as error:
        proc.kill()
        proc.wait()
        raise HLT17SRCQ1REC1AuthorityError(
            "process table cannot be observed"
        ) from error
    if len(raw) > _PS_MAX_BYTES or code != 0:
        _fail("process table is unbounded or unavailable")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise HLT17SRCQ1REC1AuthorityError("process table is not text") from error
    rows: list[tuple[int, int, str]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split(None, 2)
        if len(parts) < 2:
            _fail("process table is malformed")
        try:
            pid = int(parts[0])
            ppid = int(parts[1])
        except ValueError:
            _fail("process table is malformed")
        command = parts[2] if len(parts) == 3 else ""
        rows.append((pid, ppid, command))
    return tuple(rows)


def _protected_pids(rows: tuple[tuple[int, int, str], ...]) -> set[int]:
    parents = {pid: ppid for pid, ppid, _command in rows}
    protected = {os.getpid(), os.getppid()}
    current = os.getpid()
    seen: set[int] = set()
    while current and current not in seen:
        seen.add(current)
        protected.add(current)
        current = parents.get(current, 0)
    return protected


def _require_process_preflight() -> None:
    rows = _process_rows()
    protected = _protected_pids(rows)
    for pid, _ppid, command in rows:
        if pid in protected:
            continue
        if any(marker in command for marker in _FORBIDDEN_PROCESS):
            _fail(
                "another scientific runner, test partition, or review is active"
            )


def _require_attempt2_config(config: Mapping[str, object]) -> None:
    payload = config.get("attempt2")
    expected = _attempt2_config()
    if type(payload) is not dict or tuple(payload) != _ATTEMPT2_FIELDS:
        _fail("REC1 Attempt2 config fields differ")
    for key, value in expected.items():
        if payload.get(key) != value:
            _fail(f"REC1 Attempt2 {key} differs")


def validate_authority_delta(
    *,
    repository_root: Path,
    authority_commit: str,
    implementation_sha256: str,
) -> dict[str, object]:
    """Freeze-B exact-delta receipt. Never constructs source or runtime."""

    root = _root(Path(repository_root))
    authority = _commit(authority_commit, "authority commit")
    supplied_sha = _sha(implementation_sha256, "implementation_sha256")
    implementation = _require_committed_image(root, authority)
    try:
        config_raw = read_regular_file(root, CONFIG_PATH)
    except EvidenceIOError as error:
        _wrap(error, "live REC1 freeze config cannot be read")
    config = _parse_config(config_raw)
    if config["implementation_commit"] != implementation:
        _fail("REC1 config implementation parent differs")
    committed_hashes = _committed_implementation_hashes(root, implementation)
    live_hashes = _live_implementation_hashes(root)
    if committed_hashes != live_hashes:
        _fail("live implementation hashes differ from the implementation parent")
    if (
        supplied_sha != committed_hashes["implementation_sha256"]
        or config["implementation_sha256"] != committed_hashes["implementation_sha256"]
        or config["runner_sha256"] != committed_hashes["runner_sha256"]
        or config["authority_sha256"] != committed_hashes["authority_sha256"]
    ):
        _fail("REC1 core/runner/auth implementation hash differs")
    _require_attempt2_config(config)
    compact = authenticate_attempt2_compacts(root)
    raw = authenticate_attempt2_raw(root)
    configs = _source_configuration_mapping(config)
    if configs != raw["source_configuration"]:
        _fail("REC1 source-configuration mapping differs")
    observed_environment = environment_identity()
    if config["environment"] != observed_environment:
        _fail("REC1 environment differs")
    source = source_closure_identity(
        root, implementation, environment=observed_environment
    )
    if config["source_closure_sha256"] != source["sha256"]:
        _fail("REC1 source closure differs")
    identity = hlt17_implementation_identity()
    if config["implementation_identity"] != identity:
        _fail("REC1 C1R1 implementation identity differs")
    resources = recommended_resource_ceilings().as_mapping()
    if config["resource"] != resources:
        _fail("REC1 resource ceilings differ")
    caps = _cap_mapping(config)
    require_output_absent(root)
    require_pro20_namespace_absent(root)
    _require_process_preflight()
    return {
        "artifact_id": ARTIFACT_ID,
        "authority_commit": authority,
        "implementation_commit": implementation,
        "implementation_sha256": committed_hashes["implementation_sha256"],
        "runner_sha256": committed_hashes["runner_sha256"],
        "authority_sha256": committed_hashes["authority_sha256"],
        "exact_delta_paths": list(ALLOWED_AUTHORITY_DELTA),
        "source_closure_sha256": source["sha256"],
        "source_closure_file_count": source["file_count"],
        "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
        "output_namespace": OUTPUT_NAMESPACE,
        "environment": observed_environment,
        "resource_ceilings": resources,
        "implementation_identity": identity,
        "prospective_requested_cap_hex_by_member": caps,
        "source_configuration": configs,
        "attempt2": {
            "sealed_historical_blob": compact["sealed_historical_blob"],
            "metadata_linked_blob": compact["metadata_linked_blob"],
            "raw_leaf_relative": raw["raw_leaf_relative"],
            "raw_leaf_sha256": raw["raw_leaf_sha256"],
            "raw_directory_payload_sha256": raw["raw_directory_payload_sha256"],
            "raw_leaf_bytes": raw["raw_leaf_bytes"],
        },
        "accepted_state_advanced": False,
        "campaign_store_written": False,
        "endpoint_adopted_or_serialized": False,
        "pro20_namespace_absent": True,
        "retry_authorized": False,
        "boolean_flag_accepted": False,
        "campaign_execution_authorized": False,
        "store_publication_authorized": False,
        "endpoint_adoption_authorized": False,
        "physical_source_qualification_executed": False,
        "rk_members_remeasured": False,
    }


__all__ = [
    "ALLOWED_AUTHORITY_DELTA",
    "ARTIFACT_ID",
    "AUTHORITY_PATH",
    "CONFIG_PATH",
    "HARNESS_PATH",
    "HLT17SRCQ1REC1AuthorityError",
    "MAKE_PATH",
    "ORIGIN_CAPTURE_SHA256",
    "OUTPUT_NAMESPACE",
    "OWNER_PATH",
    "PRO20_EVENT_NAMESPACE",
    "RUNNER_PATH",
    "SOURCE_CONFIGURATION_SHA256_BY_MEMBER",
    "authenticate_attempt2_compacts",
    "authenticate_attempt2_raw",
    "check_tracked_config",
    "emit_config_bytes",
    "environment_identity",
    "expected_config",
    "render_config",
    "require_output_absent",
    "require_pro20_namespace_absent",
    "source_closure_identity",
    "validate_authority_delta",
]
