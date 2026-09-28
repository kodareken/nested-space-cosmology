"""Exact-delta authority for the prospective RSRC1 common-event run.

This is the implementation-parent authority. A later freeze commit may add
only the declared freeze delta. Validation authenticates that relationship,
the clean Git image, config/source/runner/authority hashes, source closure,
REC1-PREF1 compact license, historical HLT15/source tree, static inputs,
environment, absent output/staging namespace, process isolation, Planck and
runs baselines, and resource/disk premises. It never constructs the physical
source, never seeds the production store, and does not import historical
PRO19/HLT16 runner or recovery writers.
"""

from __future__ import annotations

from hashlib import sha256
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import re
import resource
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
from .hlt17_imp1_cursor import (
    HLT17_CFL_RETRY_CAP,
    HLT17_PRODUCTION_MEMBERS,
    HLT17_SOURCE_RETRY_CAP,
)
from .pro20_origin import (
    MEMBER_KEYS,
    PREF27_COMPACT_RELATIVE,
    SOURCE_TREE_RELATIVE,
    Pro20OriginError,
    capture_pro20_historical_origin,
    inventory_source_tree,
)
from .proto19_gr0_static_factory import STATIC_INPUT_PATHS
from .pro20_rsrc1_isolation import (
    CHILD_MAX_PEAK_RSS_BYTES,
    HOST_RESERVE_BYTES,
    MIN_AVAILABLE_MEMORY_BYTES,
    PARENT_MAX_CURRENT_RSS_BYTES,
)
from .pro20_rsrc1_runtime import (
    IMPLEMENTATION_PATHS,
    available_memory_bytes,
    current_parent_rss_bytes,
    implementation_identity,
)


ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1-FRZ1"
CONFIG_PATH = "configs/fgc/fgc-1-pro20-ev1-rsrc1-frz1.toml"
OWNER_PATH = "docs/fgc-pro20-ev1-rsrc1-frz1.md"
MAKE_PATH = "mk/current-foundation.mk"
CATALOG_BUILDER_PATH = "scripts/build_artifact_catalog.py"
CATALOG_OUTPUT_PATH = "configs/fgc/artifact-catalog.json"
CATALOG_INDEX_PATH = "results/README.md"
CATALOG_CHECKER_PATH = "scripts/repo_checks/catalog.py"
CATALOG_TEST_PATH = "tests/test_phase_minus1_artifact_catalog.py"
MAKE_TEST_PATH = "tests/test_phase_minus1_make_routing.py"
README_PATH = "README.md"
CHANGELOG_PATH = "CHANGELOG.md"
ACTIVE_MAP_PATH = "docs/active-code-map.md"
CLAIM_LEDGER_PATH = "docs/claim-ledger.md"
RUNTIME_MATRIX_PATH = "docs/fgc-runtime-matrix.md"
ROADMAP_PATH = "docs/research-roadmap.md"
HARNESS_PATH = "src/recursive_horizons/fgc/evolution/pro20_rsrc1_runtime.py"
RUNNER_PATH = "scripts/run_fgc_pro20_rsrc1.py"
AUTHORITY_PATH = "src/recursive_horizons/fgc/evolution/pro20_rsrc1_authority.py"
OUTPUT_NAMESPACE = "runs/fgc-2-sf1/pro20-rsrc1-event1/event"
PRODUCTION_CONTAINER = "runs/fgc-2-sf1/pro20-rsrc1-event1"
REC1_PREF1_COMPACT_PATH = "results/fgc-1-hlt17-srcq1-rec1-pref1.json"
REC1_PREF1_CONFIG_PATH = "configs/fgc/fgc-1-hlt17-srcq1-rec1-pref1.toml"
REC1_OUTPUT_NAMESPACE = (
    "runs/fgc-2-sf1/hlt17-srcq1-rec1-source-origin-qualification"
)
ALLOWED_FREEZE_DELTA = (
    CONFIG_PATH,
    OWNER_PATH,
)
CFL_MAXIMUM_HEX = "0x1.0000000000000p-3"
ORIGIN_CAPTURE_SHA256 = (
    "2824af06c674336d479ec53faa4269f02fb598055860ec22f2ad84602a83fe72"
)
REC1_PREF1_COMPACT_SHA256 = (
    "ea69bc1fd7ee31c41982f7eabfc3bf3e93c7e237765e71e5b112dadf69b34e75"
)
PREF27_COMPACT_RAW_SHA256 = (
    "102b385d058066d668e9e8541c974c736c583f6d5c4afaf9b99e95a16aabeb83"
)
PREF27_TREE_SHA256 = (
    "95a1bb96570216b8007e32a0cc05f4c31f420052f1d90c1b7134fa27ce958585"
)
HLT15_CHECKPOINT_RAW_SHA256 = (
    "9bcbc68ca2f168d9f7ae87c37346d9ebf90aeaf2568f26546441634734511b3b"
)
HLT15_RECEIPT_RAW_SHA256 = (
    "7e99fb5fb70b6d2d590232d21222e9aa510b823dfdb3342aa86489b21d4d2f4a"
)
GENERATION1_BRIDGE_CONTENT_ID = (
    "23f4a1ac64af5e1463862bf67ca240980921611a4a77d31050b0f92c81d44895"
)
GENERATION1_BRIDGE_RAW_SHA256 = (
    "b760ee80dca4c1334c45dd9dcab6669a7e538c7881b497222a50da0d720c8af1"
)
PROTO17_LEAF_COUNT = 56
PROTO17_BYTE_COUNT = 6054860
PROTO17_CLEANUP_DIGEST = (
    "8e0b1f3ba70a4de7061d322ed041ff0328ae4e6ca7ed78da79421883f47ad777"
)
HISTORICAL_RUNS_ROOT = "runs"
HISTORICAL_RUNS_FILE_COUNT = 1381
HISTORICAL_RUNS_BYTE_COUNT = 305779519
HISTORICAL_RUNS_DIGEST = (
    "5878f23a69d55776e3a42e8aef86b4a764a77b49e03d46c3614a9d66d9969b12"
)
PLANCK_FILES = MappingProxyType(
    {
        "collected-data/smica_2048.fits": (
            2013312960,
            "60952c645eb33d151905ddf5837477e15ca02a3261feca4ceca3c3feece4f9ac",
        ),
        "collected-data/mask_common.fits": (
            201335040,
            "23f49a3479073408b83a7f1bb2cca490d9be7f934b2e0e891524b6a1d06a5fc4",
        ),
    }
)
STATIC_INPUT_SHA256 = MappingProxyType(
    {
        "configs/fgc/fgc-1-cal9-run1.toml": (
            "8e6fafff638ffd49549efab8b6f969a01964c915116db9ce3698953201172176"
        ),
        "configs/fgc/fgc-1-rsp2-run1.toml": (
            "8213341f6cc5dd1d044a1ee9ae68b9ce7c7445ef9b933bcf8bcc944d30e4abd1"
        ),
        "results/fgc-1-hlt10-mon10.json": (
            "630d0d1843130180cb49dd237fe32e83e0336e7e44e7c684e4baf406b1883fe5"
        ),
        "results/fgc-1-rsp2-frz1.json": (
            "a6c147d7b369387e70318ee42eb56c97332702a227fd93d7ad20a1cfe2c1aa22"
        ),
    }
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
MAX_TOTAL_WALL_SECONDS = 86400.0
MAX_RSS_BYTES = CHILD_MAX_PEAK_RSS_BYTES
MAX_NAMESPACE_BYTES = 17179869184
MIN_FREE_DISK_BYTES = 34359738368
MAX_PUBLISHED_ATTEMPTS = 4096
PLANNING_ACCEPTED_STEP_ESTIMATE = 224
_COMMIT = re.compile(r"\A(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_SHA = re.compile(r"\A[0-9a-f]{64}\Z")
_SOURCE_DOMAIN = b"FGC-1-PRO20-EV1-RSRC1-SOURCE-CLOSURE-v1\n"
_SOURCE_TREE = "src/recursive_horizons"
_STAGING_PREFIX = ".pro20-rsrc1-event1.evidence-io-stage-"
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
_CONFIG_TOP_LEVEL = (
    "artifact_id",
    "schema_version",
    "implementation_commit",
    "implementation_sha256",
    "runner_sha256",
    "authority_sha256",
    "source_closure_sha256",
    "origin_capture_sha256",
    "rec1_pref1_compact_sha256",
    "output_namespace",
    "campaign_execution_authorized",
    "calibration_eligible",
    "resume_authorized",
    "takeover_authorized",
    "rsrc1_namespace_absent",
    "resource",
    "environment",
    "implementation_identity",
    "source_configuration",
    "member_cap",
    "historical",
    "planning",
)
_FALSE_FLAGS = (
    "campaign_execution_authorized",
    "calibration_eligible",
    "resume_authorized",
    "takeover_authorized",
)


class PRO20RSRC1AuthorityError(ValueError):
    """Freeze successor or the live implementation/environment differs."""

    physical_source_constructed = False


def _fail(message: str) -> NoReturn:
    raise PRO20RSRC1AuthorityError(message)


def _wrap(error: Exception, message: str) -> NoReturn:
    raise PRO20RSRC1AuthorityError(message) from error


def _root(value: Path) -> Path:
    if not isinstance(value, Path):
        value = Path(value)
    supplied = Path(os.fspath(value))
    if not supplied.is_absolute():
        _fail("repository root is not canonical")
    try:
        resolved = supplied.resolve(strict=True)
    except OSError as error:
        raise PRO20RSRC1AuthorityError("repository root cannot be resolved") from error
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
        raise PRO20RSRC1AuthorityError("software path cannot be read") from error
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
        raise PRO20RSRC1AuthorityError("software path cannot be opened") from error
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
        raise PRO20RSRC1AuthorityError("runtime path cannot be resolved") from error
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


def observe_environment() -> dict[str, str]:
    import numpy as np

    configuration = np.show_config(mode="dicts")
    if type(configuration) is not dict:
        _fail("NumPy build configuration is unavailable")
    dependencies = configuration.get("Build Dependencies", {})
    blas = dependencies.get("blas", {}) if type(dependencies) is dict else {}
    if type(blas) is not dict:
        _fail("NumPy BLAS configuration is unavailable")
    observed = {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": str(np.__version__),
        "blas_name": str(blas.get("name") or ""),
        "blas_version": str(blas.get("version") or "unknown"),
        "system": platform.system(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
        "executable": sys.executable,
    }
    for key, value in observed.items():
        if type(value) is not str or not value:
            _fail(f"environment {key} is empty")
    return observed


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


def environment_sha256(environment: Mapping[str, str]) -> str:
    if tuple(environment) != ENVIRONMENT_KEYS:
        _fail("environment keys differ")
    return _sha256_hex(canonical_json_bytes(dict(environment)))


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
    for path in (RUNNER_PATH, AUTHORITY_PATH, HARNESS_PATH):
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
    """Deterministic TOML for the closed PRO20 freeze schema; no file writes."""

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
            "source_configuration",
            "member_cap",
            "historical",
            "planning",
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
    for table in ("environment", "implementation_identity"):
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
    historical = config["historical"]
    if type(historical) is not dict:
        _fail("emitted historical table differs")
    lines.extend(("", "[historical]"))
    nested = {
        "static_inputs",
        "planck",
        "historical_runs",
    }
    for key, item in historical.items():
        if key in nested:
            continue
        lines.append(f"{key} = {_toml_literal(item)}")
    static_inputs = historical["static_inputs"]
    if type(static_inputs) is not list:
        _fail("emitted static-input inventory differs")
    for row in static_inputs:
        if type(row) is not dict:
            _fail("emitted static-input row differs")
        lines.extend(("", "[[historical.static_inputs]]"))
        for key, item in row.items():
            lines.append(f"{key} = {_toml_literal(item)}")
    planck = historical["planck"]
    if type(planck) is not list:
        _fail("emitted Planck inventory differs")
    for row in planck:
        if type(row) is not dict:
            _fail("emitted Planck row differs")
        lines.extend(("", "[[historical.planck]]"))
        for key, item in row.items():
            lines.append(f"{key} = {_toml_literal(item)}")
    runs = historical["historical_runs"]
    if type(runs) is not dict:
        _fail("emitted historical-runs table differs")
    lines.extend(("", "[historical.historical_runs]"))
    for key, item in runs.items():
        lines.append(f"{key} = {_toml_literal(item)}")
    planning = config["planning"]
    if type(planning) is not dict:
        _fail("emitted planning table differs")
    lines.extend(("", "[planning]"))
    for key, item in planning.items():
        lines.append(f"{key} = {_toml_literal(item)}")
    return ("\n".join(lines) + "\n").encode("ascii")


def recommended_resource_ceilings() -> dict[str, object]:
    return {
        "max_rss_bytes": MAX_RSS_BYTES,
        "child_max_peak_rss_bytes": CHILD_MAX_PEAK_RSS_BYTES,
        "parent_max_current_rss_bytes": PARENT_MAX_CURRENT_RSS_BYTES,
        "host_reserve_bytes": HOST_RESERVE_BYTES,
        "min_available_memory_bytes": MIN_AVAILABLE_MEMORY_BYTES,
        "max_total_wall_seconds": MAX_TOTAL_WALL_SECONDS,
        "max_wall_seconds_per_member": {
            key: MAX_TOTAL_WALL_SECONDS for key in MEMBER_KEYS
        },
        "max_namespace_bytes": MAX_NAMESPACE_BYTES,
        "min_free_disk_bytes": MIN_FREE_DISK_BYTES,
        "max_published_attempts": MAX_PUBLISHED_ATTEMPTS,
        "source_retry_cap": HLT17_SOURCE_RETRY_CAP,
        "cfl_retry_cap": HLT17_CFL_RETRY_CAP,
        "scientific_threshold_fitted": False,
        "coordinator_freeze_required": True,
        "production_enclosure_limits": [
            {
                "member_key": item.member_key,
                "point_count": item.point_count,
                "owned_rows": item.owned_rows,
                "fallback_d01": item.fallback_d01,
                "fallback_d12": item.fallback_d12,
            }
            for item in HLT17_PRODUCTION_MEMBERS
        ],
    }


def _cap_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key, spacing, speed, requested in MEMBER_CAP_INPUTS:
        if key not in MEMBER_KEYS:
            _fail("member-cap inventory differs")
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
    return [
        {"member_key": key, "sha256": SOURCE_CONFIGURATION_SHA256_BY_MEMBER[key]}
        for key in MEMBER_KEYS
    ]


def _static_input_rows() -> list[dict[str, str]]:
    return [
        {"path": path, "sha256": STATIC_INPUT_SHA256[path]}
        for path in STATIC_INPUT_PATHS
    ]


def _planck_rows() -> list[dict[str, object]]:
    return [
        {"path": path, "byte_count": size, "sha256": digest}
        for path, (size, digest) in PLANCK_FILES.items()
    ]


def _historical_config() -> dict[str, object]:
    return {
        "pref27_compact_path": PREF27_COMPACT_RELATIVE,
        "pref27_compact_raw_sha256": PREF27_COMPACT_RAW_SHA256,
        "pref27_tree_sha256": PREF27_TREE_SHA256,
        "hlt15_checkpoint_raw_sha256": HLT15_CHECKPOINT_RAW_SHA256,
        "hlt15_receipt_raw_sha256": HLT15_RECEIPT_RAW_SHA256,
        "generation1_bridge_content_id": GENERATION1_BRIDGE_CONTENT_ID,
        "generation1_bridge_raw_sha256": GENERATION1_BRIDGE_RAW_SHA256,
        "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
        "proto17_relative": SOURCE_TREE_RELATIVE,
        "proto17_leaf_count": PROTO17_LEAF_COUNT,
        "proto17_byte_count": PROTO17_BYTE_COUNT,
        "proto17_cleanup_digest": PROTO17_CLEANUP_DIGEST,
        "rec1_pref1_compact_path": REC1_PREF1_COMPACT_PATH,
        "rec1_pref1_compact_sha256": REC1_PREF1_COMPACT_SHA256,
        "static_inputs": _static_input_rows(),
        "planck": _planck_rows(),
        "historical_runs": {
            "file_count": HISTORICAL_RUNS_FILE_COUNT,
            "byte_count": HISTORICAL_RUNS_BYTE_COUNT,
            "digest": HISTORICAL_RUNS_DIGEST,
            "excluded_namespace": REC1_OUTPUT_NAMESPACE,
        },
    }


def expected_config(
    *,
    implementation_commit: str,
    implementation_hashes: Mapping[str, str],
    source_closure_sha256: str,
    environment: Mapping[str, str],
) -> dict[str, object]:
    resources = recommended_resource_ceilings()
    if (
        resources["max_rss_bytes"] != MAX_RSS_BYTES
        or resources["max_total_wall_seconds"] != MAX_TOTAL_WALL_SECONDS
        or resources["scientific_threshold_fitted"] is not False
        or resources["source_retry_cap"] != 32
        or resources["cfl_retry_cap"] != 32
        or resources["max_published_attempts"] != MAX_PUBLISHED_ATTEMPTS
    ):
        _fail("recommended resource ceilings differ")
    config: dict[str, object] = {
        "artifact_id": ARTIFACT_ID,
        "schema_version": 1,
        "implementation_commit": implementation_commit,
        "implementation_sha256": implementation_hashes["implementation_sha256"],
        "runner_sha256": implementation_hashes["runner_sha256"],
        "authority_sha256": implementation_hashes["authority_sha256"],
        "source_closure_sha256": source_closure_sha256,
        "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
        "rec1_pref1_compact_sha256": REC1_PREF1_COMPACT_SHA256,
        "output_namespace": OUTPUT_NAMESPACE,
        "campaign_execution_authorized": False,
        "calibration_eligible": False,
        "resume_authorized": False,
        "takeover_authorized": False,
        "rsrc1_namespace_absent": True,
        "resource": resources,
        "environment": dict(environment),
        "implementation_identity": hlt17_implementation_identity(),
        "source_configuration": _source_configuration_rows(),
        "member_cap": _cap_rows(),
        "historical": _historical_config(),
        "planning": {
            "accepted_step_estimate": PLANNING_ACCEPTED_STEP_ESTIMATE,
            "accepted_step_estimate_is_not_a_gate": True,
        },
    }
    if tuple(config) != _CONFIG_TOP_LEVEL:
        _fail("emitted config field order differs")
    return config


def _live_implementation_hashes(root: Path) -> dict[str, str]:
    aggregate = implementation_identity(root)
    hashes = {
        "implementation_sha256": str(aggregate["implementation_sha256"]),
        "runner_sha256": _sha256_hex(read_regular_file(root, RUNNER_PATH)),
        "authority_sha256": _sha256_hex(read_regular_file(root, AUTHORITY_PATH)),
    }
    for label, value in hashes.items():
        _sha(value, label)
    return hashes


def _committed_implementation_hashes(root: Path, commit: str) -> dict[str, str]:
    path_hashes = {
        relative: _sha256_hex(_blob(root, commit, relative))
        for relative in IMPLEMENTATION_PATHS
    }
    aggregate = _sha256_hex(
        "".join(
            f"{path}\0{path_hashes[path]}\n" for path in IMPLEMENTATION_PATHS
        ).encode("ascii")
    )
    hashes = {
        "implementation_sha256": aggregate,
        "runner_sha256": _sha256_hex(_blob(root, commit, RUNNER_PATH)),
        "authority_sha256": _sha256_hex(_blob(root, commit, AUTHORITY_PATH)),
    }
    for label, value in hashes.items():
        _sha(value, label)
    return hashes


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
        _wrap(error, "tracked PRO20 freeze config cannot be read")
    parsed = _parse_config(raw)
    expected = emit_config_bytes(
        repository,
        implementation_commit=str(parsed["implementation_commit"]),
    )
    if raw != expected:
        _fail("tracked PRO20 freeze config differs from emitter")
    return {
        "config_sha256": _sha256_hex(raw),
        "matched": True,
        "authorized": False,
        "executed": False,
        "physical_source_constructed": False,
        "campaign_execution_authorized": False,
        "calibration_eligible": False,
    }


def _parse_config(raw: bytes) -> dict[str, object]:
    if type(raw) is not bytes:
        _fail("PRO20 freeze config is not bytes")
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise PRO20RSRC1AuthorityError("PRO20 freeze config is malformed") from error
    if type(value) is not dict or tuple(value) != _CONFIG_TOP_LEVEL:
        _fail("PRO20 freeze config fields differ")
    if value["artifact_id"] != ARTIFACT_ID or value["schema_version"] != 1:
        _fail("PRO20 freeze config identity differs")
    if value["output_namespace"] != OUTPUT_NAMESPACE:
        _fail("PRO20 output namespace differs")
    if value["origin_capture_sha256"] != ORIGIN_CAPTURE_SHA256:
        _fail("PRO20 origin capture differs")
    if value["rec1_pref1_compact_sha256"] != REC1_PREF1_COMPACT_SHA256:
        _fail("REC1-PREF1 compact license hash differs")
    if value["rsrc1_namespace_absent"] is not True:
        _fail("PRO20 namespace-absent flag differs")
    for name in _FALSE_FLAGS:
        if value[name] is not False:
            _fail(f"PRO20 {name} flag differs")
    return value


def _cap_mapping(config: Mapping[str, object]) -> dict[str, str]:
    rows = config.get("member_cap")
    expected_rows = _cap_rows()
    if type(rows) is not list or len(rows) != len(MEMBER_KEYS):
        _fail("PRO20 member-cap inventory differs")
    result: dict[str, str] = {}
    for expected, row in zip(MEMBER_KEYS, rows, strict=True):
        if type(row) is not dict or set(row) != {
            "member_key",
            "grid_spacing_hex",
            "previous_speed_upper_hex",
            "cfl_maximum_hex",
            "requested_cap_hex",
        }:
            _fail("PRO20 member-cap fields differ")
        if row["member_key"] != expected:
            _fail("PRO20 member-cap order differs")
        wanted = expected_rows[len(result)]
        for name in (
            "member_key",
            "grid_spacing_hex",
            "previous_speed_upper_hex",
            "cfl_maximum_hex",
        ):
            if row[name] != wanted[name]:
                _fail("PRO20 member-cap values differ")
        try:
            spacing = float.fromhex(row["grid_spacing_hex"])
            speed = float.fromhex(row["previous_speed_upper_hex"])
            cfl = float.fromhex(row["cfl_maximum_hex"])
            requested = float.fromhex(row["requested_cap_hex"])
        except (TypeError, ValueError) as error:
            raise PRO20RSRC1AuthorityError("PRO20 member cap is not binary64 hex") from error
        if min(spacing, speed, cfl, requested) <= 0.0:
            _fail("PRO20 member-cap formula differs")
        if (cfl * spacing / speed).hex() != requested.hex():
            _fail("PRO20 member-cap formula differs")
        if requested.hex() != row["requested_cap_hex"]:
            _fail("PRO20 requested-cap encoding differs")
        if row["requested_cap_hex"] != wanted["requested_cap_hex"]:
            _fail("PRO20 requested cap is not the frozen SRCQ1 hex")
        result[expected] = row["requested_cap_hex"]
    return result


def _source_configuration_mapping(config: Mapping[str, object]) -> dict[str, str]:
    rows = config.get("source_configuration")
    if type(rows) is not list or len(rows) != len(MEMBER_KEYS):
        _fail("PRO20 source-configuration inventory differs")
    result: dict[str, str] = {}
    for expected, row in zip(MEMBER_KEYS, rows, strict=True):
        if type(row) is not dict or set(row) != {"member_key", "sha256"}:
            _fail("PRO20 source-configuration fields differ")
        if row["member_key"] != expected:
            _fail("PRO20 source-configuration order differs")
        digest = _sha(row["sha256"], f"{expected} source configuration")
        if digest != SOURCE_CONFIGURATION_SHA256_BY_MEMBER[expected]:
            _fail("PRO20 source-configuration mapping differs")
        result[expected] = digest
    return result


def _json_object(raw: bytes, label: str) -> dict[str, object]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PRO20RSRC1AuthorityError(f"{label} is not JSON") from error
    if type(value) is not dict:
        _fail(f"{label} is not an object")
    return value


def _require_committed_image(root: Path, authority_commit: str) -> str:
    head = _git_query(root, ResolveCommit("HEAD"))
    if head != authority_commit:
        _fail("live HEAD is not the PRO20 authority commit")
    parents = _git_query(root, CommitParents(authority_commit))
    if type(parents) is not tuple or len(parents) != 1:
        _fail("PRO20 authority must have one implementation parent")
    implementation = parents[0]
    if type(implementation) is not str or _COMMIT.fullmatch(implementation) is None:
        _fail("PRO20 implementation parent is not an exact commit")
    delta = _git_query(
        root, InspectDelta(commit=authority_commit, against=implementation)
    )
    paths = tuple(sorted(item.path for item in delta))
    if paths != tuple(sorted(ALLOWED_FREEZE_DELTA)):
        _fail("PRO20 authority delta paths differ")
    worktree = _git_query(root, InspectWorktree())
    if not getattr(worktree, "clean", False):
        _fail("PRO20 authority requires a clean worktree")
    for path in ALLOWED_FREEZE_DELTA:
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
        raise PRO20RSRC1AuthorityError("namespace path cannot be read") from error
    return True


def _require_relative(relative: str, label: str) -> Path:
    value = Path(relative)
    if value.is_absolute() or any(part in {"", ".", ".."} for part in value.parts):
        _fail(f"{label} is unsafe")
    return value


def require_output_absent(root: Path) -> None:
    """Require the PRO20 container, calibration store, and staging names absent."""

    repository = _root(root)
    relative = _require_relative(PRODUCTION_CONTAINER, "PRO20 production container")
    current = repository
    for part in relative.parent.parts:
        current = current / part
        if not _lexists(current):
            return
        try:
            metadata = current.lstat()
        except OSError as error:
            raise PRO20RSRC1AuthorityError("output parent cannot be read") from error
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _fail("output parent is unsafe")
    container = repository / relative
    store = repository / _require_relative(OUTPUT_NAMESPACE, "PRO20 output namespace")
    if _lexists(container) or _lexists(store):
        _fail("PRO20 output namespace already exists")
    parent = container.parent
    if not _lexists(parent):
        return
    try:
        names = os.listdir(parent)
    except OSError as error:
        raise PRO20RSRC1AuthorityError("output parent cannot be listed") from error
    if any(
        name.startswith(_STAGING_PREFIX)
        for name in names
    ):
        _fail("staging namespace already exists")


def authenticate_rec1_pref1_compact(root: Path) -> dict[str, object]:
    """Authenticate the independent REC1-PREF1 compact license. No source call."""

    repository = _root(root)
    try:
        raw = read_regular_file(repository, REC1_PREF1_COMPACT_PATH)
    except EvidenceIOError as error:
        _wrap(error, "REC1-PREF1 compact cannot be read")
    digest = _sha256_hex(raw)
    if digest != REC1_PREF1_COMPACT_SHA256:
        _fail("REC1-PREF1 compact SHA-256 differs")
    payload = _json_object(raw, "REC1-PREF1 compact")
    inner = payload.get("artifact_payload", payload)
    if type(inner) is not dict:
        _fail("REC1-PREF1 compact payload differs")
    if payload.get("artifact_id") != "FGC-1-HLT17-SRCQ1-REC1-PREF1":
        _fail("REC1-PREF1 compact identity differs")
    if (
        inner.get("classification")
        != "independently_bound_all_six_source_c1r1_qualified_no_state_advance"
    ):
        _fail("REC1-PREF1 compact classification differs")
    if inner.get("pro20_event_execution_authorized") is not False:
        _fail("REC1-PREF1 licenses PRO20 execution")
    claims = inner.get("claims")
    if type(claims) is not dict:
        _fail("REC1-PREF1 claims are absent")
    if claims.get("pro20_event_execution_authorized") is not False:
        _fail("REC1-PREF1 claim licenses PRO20 execution")
    if claims.get("pro20_first_event_authority_licensed_separately") is not True:
        _fail("REC1-PREF1 does not license a separate PRO20 authority")
    if claims.get("independently_bound_all_six_source_c1r1_qualified") is not True:
        _fail("REC1-PREF1 all-six source/C1R1 qualification differs")
    origin = None
    authority = inner.get("authority")
    if type(authority) is dict:
        origin = authority.get("origin_capture_sha256")
    if origin is not None and origin != ORIGIN_CAPTURE_SHA256:
        _fail("REC1-PREF1 origin capture differs")
    return {
        "path": REC1_PREF1_COMPACT_PATH,
        "sha256": digest,
        "artifact_id": "FGC-1-HLT17-SRCQ1-REC1-PREF1",
        "pro20_event_execution_authorized": False,
        "physical_source_constructed": False,
    }


def authenticate_historical_source_tree(root: Path) -> dict[str, object]:
    """Hash the preserved proto17 tree and PREF27/HLT15 leaves. No RHS call."""

    repository = _root(root)
    try:
        compact = read_regular_file(repository, PREF27_COMPACT_RELATIVE)
    except EvidenceIOError as error:
        _wrap(error, "PREF27 compact cannot be read")
    if _sha256_hex(compact) != PREF27_COMPACT_RAW_SHA256:
        _fail("PREF27 compact raw SHA-256 differs")
    try:
        inventory = inventory_source_tree(repository, SOURCE_TREE_RELATIVE)
    except Pro20OriginError as error:
        _wrap(error, "historical HLT15/source tree cannot be inventoried")
    if (
        inventory.leaf_count != PROTO17_LEAF_COUNT
        or inventory.byte_count != PROTO17_BYTE_COUNT
        or inventory.cleanup_digest != PROTO17_CLEANUP_DIGEST
    ):
        _fail("historical HLT15/source tree differs")
    try:
        capture = capture_pro20_historical_origin(repository)
    except Pro20OriginError as error:
        _wrap(error, "historical HLT15 generation-one capture differs")
    if capture.capture_sha256 != ORIGIN_CAPTURE_SHA256:
        _fail("historical origin capture SHA-256 differs")
    if capture.generation1_bridge_content_id != GENERATION1_BRIDGE_CONTENT_ID:
        _fail("historical generation-one bridge content ID differs")
    if (
        capture.generation1_checkpoint.content_id != GENERATION1_BRIDGE_CONTENT_ID
        or capture.generation1_checkpoint.raw_sha256
        != GENERATION1_BRIDGE_RAW_SHA256
    ):
        _fail("historical generation-one checkpoint differs")
    if (
        capture.pref27_compact.raw_sha256 != PREF27_COMPACT_RAW_SHA256
        or capture.pref27_tree_sha256 != PREF27_TREE_SHA256
    ):
        _fail("historical PREF27 capture differs")
    originals = {item.content_id: item for item in capture.original_leaves}
    checkpoint = originals.get(
        "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d"
    )
    receipt = originals.get(
        "78ae89a9f74193658003cf51f8567c4088db80fd6fb86867c2d581a0c16d5baf"
    )
    if checkpoint is None or checkpoint.raw_sha256 != HLT15_CHECKPOINT_RAW_SHA256:
        _fail("HLT15 checkpoint raw SHA-256 differs")
    if receipt is None or receipt.raw_sha256 != HLT15_RECEIPT_RAW_SHA256:
        _fail("HLT15 receipt raw SHA-256 differs")
    return {
        "proto17_leaf_count": inventory.leaf_count,
        "proto17_byte_count": inventory.byte_count,
        "proto17_cleanup_digest": inventory.cleanup_digest,
        "pref27_compact_raw_sha256": PREF27_COMPACT_RAW_SHA256,
        "hlt15_checkpoint_raw_sha256": HLT15_CHECKPOINT_RAW_SHA256,
        "hlt15_receipt_raw_sha256": HLT15_RECEIPT_RAW_SHA256,
        "generation1_bridge_content_id": GENERATION1_BRIDGE_CONTENT_ID,
        "generation1_bridge_raw_sha256": GENERATION1_BRIDGE_RAW_SHA256,
        "origin_capture_sha256": capture.capture_sha256,
        "physical_source_constructed": False,
    }


def authenticate_static_inputs(root: Path) -> dict[str, str]:
    repository = _root(root)
    observed: dict[str, str] = {}
    for relative in STATIC_INPUT_PATHS:
        try:
            raw = read_regular_file(repository, relative)
        except EvidenceIOError as error:
            _wrap(error, f"static input cannot be read: {relative}")
        digest = _sha256_hex(raw)
        if digest != STATIC_INPUT_SHA256[relative]:
            _fail(f"static input identity differs: {relative}")
        observed[relative] = digest
    if tuple(observed) != STATIC_INPUT_PATHS:
        _fail("static input order differs")
    return observed


def authenticate_planck(root: Path) -> dict[str, object]:
    repository = _root(root)
    records: dict[str, object] = {}
    for relative, (size, digest) in PLANCK_FILES.items():
        path = repository / relative
        try:
            metadata = path.lstat()
        except FileNotFoundError:
            records[relative] = {
                "present": False,
                "byte_count": size,
                "sha256": digest,
            }
            continue
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
            _fail(f"Planck path is unsafe: {relative}")
        if metadata.st_size != size:
            _fail(f"Planck byte count differs: {relative}")
        if _hash_software_file(str(path)) != digest:
            _fail(f"Planck SHA-256 differs: {relative}")
        records[relative] = {
            "present": True,
            "byte_count": size,
            "sha256": digest,
        }
    return records


def historical_runs_digest(
    records: Mapping[str, tuple[int, str]] | list[tuple[str, int, str]]
) -> str:
    if isinstance(records, Mapping):
        rows = [
            (path, int(item[0]), _sha(item[1], "runs leaf"))
            for path, item in records.items()
        ]
    else:
        rows = list(records)
    lines = "".join(
        f"{digest}  {path}\n" for path, _size, digest in sorted(rows, key=lambda item: item[0])
    )
    return _sha256_hex(lines.encode("ascii"))


def inventory_historical_runs(root: Path) -> dict[str, object]:
    repository = _root(root)
    relative = _require_relative(HISTORICAL_RUNS_ROOT, "runs root")
    base = repository / relative
    try:
        info = base.lstat()
    except OSError as error:
        raise PRO20RSRC1AuthorityError("historical runs tree is absent") from error
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        _fail("historical runs tree is unsafe")
    pending = [str(relative)]
    files: list[str] = []
    excluded_prefix = REC1_OUTPUT_NAMESPACE + "/"
    pro20_prefix = PRODUCTION_CONTAINER + "/"
    while pending:
        current_relative = pending.pop()
        current = repository / current_relative
        try:
            with os.scandir(current) as iterator:
                children = sorted(iterator, key=lambda item: item.name)
        except OSError as error:
            raise PRO20RSRC1AuthorityError("historical runs tree cannot be scanned") from error
        for child in children:
            child_relative = f"{current_relative}/{child.name}"
            try:
                child_info = child.stat(follow_symlinks=False)
            except OSError as error:
                raise PRO20RSRC1AuthorityError(
                    f"historical runs path is unsafe: {child_relative}"
                ) from error
            if stat.S_ISLNK(child_info.st_mode):
                _fail(f"historical runs path is unsafe: {child_relative}")
            if stat.S_ISDIR(child_info.st_mode):
                pending.append(child_relative)
                continue
            if not stat.S_ISREG(child_info.st_mode) or child_info.st_nlink != 1:
                _fail(f"historical runs path is unsafe: {child_relative}")
            if child_relative.startswith(excluded_prefix) or child_relative.startswith(
                pro20_prefix
            ):
                continue
            files.append(child_relative)
    files = sorted(files)
    rows: list[tuple[str, int, str]] = []
    total = 0
    for path in files:
        absolute = repository / path
        try:
            metadata = absolute.lstat()
        except OSError as error:
            raise PRO20RSRC1AuthorityError(f"historical runs path vanished: {path}") from error
        digest = _hash_software_file(str(absolute))
        rows.append((path, metadata.st_size, digest))
        total += metadata.st_size
    digest = historical_runs_digest(rows)
    if (
        len(rows) != HISTORICAL_RUNS_FILE_COUNT
        or total != HISTORICAL_RUNS_BYTE_COUNT
        or digest != HISTORICAL_RUNS_DIGEST
    ):
        _fail("historical runs baseline differs")
    return {
        "file_count": len(rows),
        "byte_count": total,
        "digest": digest,
        "excluded_namespace": REC1_OUTPUT_NAMESPACE,
    }


def current_rss_bytes() -> int:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    rss = int(usage.ru_maxrss)
    if rss < 0:
        _fail("RSS observation is negative")
    if sys.platform == "darwin":
        return rss
    return rss * 1024


def free_disk_bytes(path: Path) -> int:
    try:
        stats = os.statvfs(os.fspath(path))
    except OSError as error:
        raise PRO20RSRC1AuthorityError("free disk cannot be observed") from error
    available = int(stats.f_bavail) * int(stats.f_frsize)
    if available < 0:
        _fail("free disk observation is negative")
    return available


def _require_disk_and_rss(root: Path) -> dict[str, int]:
    rss = current_parent_rss_bytes()
    if rss > PARENT_MAX_CURRENT_RSS_BYTES:
        _fail("parent current RSS exceeds the operational ceiling")
    available = available_memory_bytes()
    if available < MIN_AVAILABLE_MEMORY_BYTES:
        _fail("available memory is below the operational reserve")
    disk = free_disk_bytes(root)
    if disk < MIN_FREE_DISK_BYTES:
        _fail("free disk is below the operational ceiling")
    return {
        "parent_current_rss_bytes": rss,
        "available_memory_bytes": available,
        "free_disk_bytes": disk,
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
        raise PRO20RSRC1AuthorityError("process table cannot be observed") from error
    assert proc.stdout is not None
    try:
        raw = _communicate_bounded(proc, max_bytes=_PS_MAX_BYTES, timeout=5)
        code = proc.returncode
    except Exception as error:
        proc.kill()
        proc.wait()
        raise PRO20RSRC1AuthorityError("process table cannot be observed") from error
    if len(raw) > _PS_MAX_BYTES or code != 0:
        _fail("process table is unbounded or unavailable")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PRO20RSRC1AuthorityError("process table is not text") from error
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
            _fail("another scientific runner, test partition, or review is active")


def _receipt_identities(
    *,
    implementation_hashes: Mapping[str, str],
    config_sha256: str,
    source_closure_sha256: str,
    environment: Mapping[str, str],
) -> dict[str, str]:
    identities = {
        "authority_sha256": implementation_hashes["authority_sha256"],
        "implementation_sha256": implementation_hashes["implementation_sha256"],
        "config_sha256": config_sha256,
        "source_sha256": source_closure_sha256,
        "origin_sha256": ORIGIN_CAPTURE_SHA256,
        "environment_sha256": environment_sha256(environment),
    }
    if len(set(identities.values())) != 6:
        _fail("authority/implementation/config/source/origin/environment identities collided")
    return identities


def validate_authority_delta(
    *,
    repository_root: Path,
    authority_commit: str,
    implementation_sha256: str,
) -> dict[str, object]:
    """Freeze exact-delta receipt. Never constructs source or runtime."""

    root = _root(Path(repository_root))
    authority = _commit(authority_commit, "authority commit")
    supplied_sha = _sha(implementation_sha256, "implementation_sha256")
    implementation = _require_committed_image(root, authority)
    try:
        config_raw = read_regular_file(root, CONFIG_PATH)
    except EvidenceIOError as error:
        _wrap(error, "live PRO20 freeze config cannot be read")
    config = _parse_config(config_raw)
    if config["implementation_commit"] != implementation:
        _fail("PRO20 config implementation parent differs")
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
        _fail("PRO20 core/runner/auth implementation hash differs")
    rec1 = authenticate_rec1_pref1_compact(root)
    historical = authenticate_historical_source_tree(root)
    static_inputs = authenticate_static_inputs(root)
    planck = authenticate_planck(root)
    runs = inventory_historical_runs(root)
    configs = _source_configuration_mapping(config)
    observed_environment = environment_identity()
    if config["environment"] != observed_environment:
        _fail("PRO20 environment differs")
    source = source_closure_identity(
        root, implementation, environment=observed_environment
    )
    if config["source_closure_sha256"] != source["sha256"]:
        _fail("PRO20 source closure differs")
    identity = hlt17_implementation_identity()
    if config["implementation_identity"] != identity:
        _fail("PRO20 C1R1 implementation identity differs")
    resources = recommended_resource_ceilings()
    if config["resource"] != resources:
        _fail("PRO20 resource ceilings differ")
    caps = _cap_mapping(config)
    if config["historical"] != _historical_config():
        _fail("PRO20 historical pins differ")
    planning = config.get("planning")
    if type(planning) is not dict:
        _fail("PRO20 planning table differs")
    if planning.get("accepted_step_estimate") != PLANNING_ACCEPTED_STEP_ESTIMATE:
        _fail("PRO20 planning estimate differs")
    if planning.get("accepted_step_estimate_is_not_a_gate") is not True:
        _fail("planning estimate must not be a pass/fail gate")
    require_output_absent(root)
    _require_process_preflight()
    disk = _require_disk_and_rss(root)
    config_digest = _sha256_hex(config_raw)
    identities = _receipt_identities(
        implementation_hashes=committed_hashes,
        config_sha256=config_digest,
        source_closure_sha256=str(source["sha256"]),
        environment=observed_environment,
    )
    return {
        "artifact_id": ARTIFACT_ID,
        "authority_commit": authority,
        "implementation_commit": implementation,
        "implementation_sha256": committed_hashes["implementation_sha256"],
        "runner_sha256": committed_hashes["runner_sha256"],
        "authority_sha256": committed_hashes["authority_sha256"],
        "config_sha256": config_digest,
        "source_closure_sha256": source["sha256"],
        "source_sha256": source["sha256"],
        "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
        "origin_sha256": ORIGIN_CAPTURE_SHA256,
        "environment_sha256": identities["environment_sha256"],
        "receipt_identities": identities,
        "exact_delta_paths": list(ALLOWED_FREEZE_DELTA),
        "source_closure_file_count": source["file_count"],
        "output_namespace": OUTPUT_NAMESPACE,
        "environment": observed_environment,
        "resource_ceilings": resources,
        "implementation_identity": identity,
        "prospective_requested_cap_hex_by_member": caps,
        "source_configuration": configs,
        "static_inputs": static_inputs,
        "rec1_pref1": rec1,
        "historical_source_tree": historical,
        "planck": planck,
        "historical_runs": runs,
        "disk": disk,
        "campaign_id": "FGC-2-SF1-PRO20-RSRC1-EVENT1",
        "campaign_execution_authorized": False,
        "calibration_eligible": False,
        "resume_authorized": False,
        "takeover_authorized": False,
        "rsrc1_namespace_absent": True,
        "boolean_flag_accepted": False,
        "physical_source_constructed": False,
        "physical_source_qualification_executed": False,
        "planning_accepted_step_estimate": PLANNING_ACCEPTED_STEP_ESTIMATE,
        "planning_accepted_step_estimate_is_not_a_gate": True,
    }


__all__ = [
    "ALLOWED_FREEZE_DELTA",
    "ARTIFACT_ID",
    "AUTHORITY_PATH",
    "CONFIG_PATH",
    "HARNESS_PATH",
    "ORIGIN_CAPTURE_SHA256",
    "OUTPUT_NAMESPACE",
    "OWNER_PATH",
    "PLANCK_FILES",
    "PRO20RSRC1AuthorityError",
    "REC1_PREF1_COMPACT_SHA256",
    "RUNNER_PATH",
    "SOURCE_CONFIGURATION_SHA256_BY_MEMBER",
    "authenticate_historical_source_tree",
    "authenticate_planck",
    "authenticate_rec1_pref1_compact",
    "authenticate_static_inputs",
    "check_tracked_config",
    "emit_config_bytes",
    "environment_identity",
    "expected_config",
    "inventory_historical_runs",
    "render_config",
    "require_output_absent",
    "source_closure_identity",
    "validate_authority_delta",
]
