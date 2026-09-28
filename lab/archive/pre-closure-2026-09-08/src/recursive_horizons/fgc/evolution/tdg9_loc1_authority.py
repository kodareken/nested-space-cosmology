"""Prospective, premise-only authority for the bounded TDG9 LOC1 diagnostic."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
import platform
from pathlib import Path
import re
import stat
import subprocess
import sys
import tomllib
from typing import NoReturn

import numpy as np


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-LOC1-FRZ1"
CLASSIFICATION = "premise_only_exact_failed_cubic_localization_authority"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
PROJECT_VERSION = "0.11.0"
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-loc1-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg9-loc1-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-loc1-frz1.md"
OUTPUT_NAMESPACE = "runs/fgc-2-sf1/tdg9-loc1/rcv3-retries-3-5"
STAGING_PREFIX = ".rcv3-retries-3-5.stage-"

BASE_COMMIT = "e6052fe62c4e8eb54b64f92fe501a923f3d0e087"
PREF1_COMMIT = BASE_COMMIT
PREF1_CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ar1-pref1.toml"
PREF1_CONFIG_SHA256 = "c9a87856988f807c60bd44ac5eca1428402e7b2aa2eff51984b9899ce878dd79"
PREF1_RESULT_PATH = "results/fgc-1-tdg9-ar1-pref1.json"
PREF1_RESULT_SHA256 = "e061bcbcebf0919ef237cb1e757b0fd9c40c0b37348b3bfb3409354c97b84bdf"
PREF1_BINDER_PATH = "src/recursive_horizons/fgc/evolution/tdg9_ar1_pref1_binder.py"
PREF1_BINDER_SHA256 = "59c9d1f52263fa58775e228f62a868d427e661b5d2c29e3ab5628f2abd61b4a6"
PREF1_OWNER_PATH = "docs/fgc-tdg9-ar1-pref1.md"
PREF1_OWNER_SHA256 = "735dfbec447a25e8046ad27d7578f13cd9056b1900266be7f7e0b75c1300b8a1"

METHOD = "RK4"
MEMBER_KEY = "RK4-2049"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
GRID_SPACING_HEX = "0x1.0000000000000p-4"
OUTER_RADIUS_HEX = "0x1.0000000000000p+7"
RETRIES = (3, 4, 5)
FAILED_OCCURRENCES = (
    (3, "u:alpha"), (3, "u:R"),
    (4, "u:alpha"), (4, "u:lambda"), (4, "u:R"),
    (5, "u:alpha"), (5, "u:v"), (5, "u:lambda"), (5, "u:R"), (5, "q:R"),
)
POLYNOMIALS_PER_OCCURRENCE = 6 * OWNED_ROW_COUNT
TOTAL_POLYNOMIALS = len(FAILED_OCCURRENCES) * POLYNOMIALS_PER_OCCURRENCE
COMPONENT_CUBIC_COUNT = TOTAL_POLYNOMIALS
EXECUTED_COMPONENT_CUBIC_COUNT = 3 * TOTAL_POLYNOMIALS
COMPLETE_C_MAXIMUM_CANDIDATES = 4 * TOTAL_POLYNOMIALS
VALUE_V_MAXIMUM_CANDIDATES = 4 * TOTAL_POLYNOMIALS
SLOPE_S_MAXIMUM_CANDIDATES = 4 * TOTAL_POLYNOMIALS
TOTAL_VSC_MAXIMUM_CANDIDATES = (
    COMPLETE_C_MAXIMUM_CANDIDATES
    + VALUE_V_MAXIMUM_CANDIDATES
    + SLOPE_S_MAXIMUM_CANDIDATES
)
MAXIMUM_CANDIDATES = COMPLETE_C_MAXIMUM_CANDIDATES
PRIMARY_REFINEMENT_DEPTH = 160
INDEPENDENT_REFINEMENT_BITS = 192
ROW_REGIONS = (
    (0, 0, "centre_parity_derivative_row"),
    (1, 2, "inner_stencil_transition"),
    (3, 2040, "smooth_interior"),
    (2041, 2041, "outer_dissipation_transition"),
    (2042, 2043, "outer_projector_influenced"),
)
ENVIRONMENT = {
    "python_implementation": "CPython", "python_version": "3.14.3",
    "numpy_version": "2.5.1", "system": "Darwin", "sys_platform": "darwin",
    "machine": "arm64", "byteorder": "little",
}

DELTA_PATHS = tuple(sorted((
    "Makefile", "README.md", CONFIG_PATH, "docs/claim-ledger.md",
    "docs/fgc-runtime-matrix.md", OWNER_DOCUMENT, "docs/research-roadmap.md",
    "results/README.md", RESULT_PATH, "scripts/check_repo.py",
    "scripts/reproduce_fgc_tdg9_loc1_frz1.py", "scripts/run_fgc_tdg9_loc1.py",
    "src/recursive_horizons/fgc/evolution/tdg9_loc1_authority.py",
    "src/recursive_horizons/fgc/evolution/tdg9_local_extrema.py",
    "src/recursive_horizons/fgc/evolution/tdg9_local_extrema_independent.py",
    "tests/test_check_repo_tdg9_loc1_frz1.py",
    "tests/test_fgc_tdg9_loc1_authority.py",
    "tests/test_fgc_tdg9_loc1_runner.py",
    "tests/test_fgc_tdg9_local_extrema.py",
)))
ALLOWED_UNTRACKED = (".qdrant-initialized",)
SCOPE = {
    "live_store_read_only": True, "compact_verifier_store_blind": True,
    "exact_retries_3_4_5_only": True, "exact_failed_occurrences_only": True,
    "historical_threshold_changed": False, "ULP_used_as_tolerance": False,
    "fourth_width_authorized": False, "resource_escalation_authorized": False,
    "campaign_store_mutation_authorized": False, "PDE_state_commit_authorized": False,
    "continuation_authorized": False, "candidate_branches_authorized": False,
    "stage2_source_decomposition_authorized": False,
    "SSPRK3_comparator_authorized": False,
}
CLAIMS = {
    "PREF1_terminal_bound": True, "localization_execution_authorized": True,
    "localization_result_earned": False, "common_event_completed": False,
    "GR0_calibration_completed": False, "candidate_execution_authorized": False,
    "mechanism_result_earned": False, "physical_result_earned": False,
}

_SHA = re.compile(r"[0-9a-f]{64}\Z")
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class LOC1AuthorityError(RuntimeError):
    pass


def _fail(message: str) -> NoReturn:
    raise LOC1AuthorityError(message)


def canonical_pretty(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("duplicate compact JSON key")
        answer[key] = value
    return answer


def _read(root: Path, relative: str, maximum: int = 16 * 1024 * 1024) -> bytes:
    candidate = Path(relative)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        _fail("unsafe authority path")
    path = root / candidate
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_size > maximum:
        _fail("unsafe authority leaf")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        active = os.fstat(descriptor)
        identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        if identity != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns):
            _fail("authority leaf raced")
        raw = b""
        while len(raw) <= maximum:
            block = os.read(descriptor, min(1 << 20, maximum + 1 - len(raw)))
            if not block:
                break
            raw += block
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) or len(raw) != before.st_size:
        _fail("authority leaf changed")
    return raw


def _git(root: Path, *arguments: str) -> bytes:
    environment = dict(os.environ, LC_ALL="C", LANG="C", GIT_NO_REPLACE_OBJECTS="1", GIT_CONFIG_NOSYSTEM="1")
    try:
        return subprocess.run(("git", "--no-replace-objects", "--no-optional-locks", "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false", *arguments), cwd=root, env=environment, check=True, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    except subprocess.CalledProcessError as exc:
        raise LOC1AuthorityError("Git image differs") from exc


def _paths(root: Path, *arguments: str) -> tuple[str, ...]:
    raw = _git(root, *arguments)
    if not raw:
        return ()
    chunks = raw.split(b"\0")
    if chunks[-1] != b"":
        _fail("Git path stream differs")
    return tuple(item.decode() for item in chunks[:-1])


def _head(root: Path) -> str:
    answer = _git(root, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if not _COMMIT.fullmatch(answer):
        _fail("HEAD identity differs")
    return answer


def _parents(root: Path, commit: str) -> tuple[str, ...]:
    fields = _git(root, "rev-list", "--parents", "-n", "1", commit).decode().split()
    return tuple(fields[1:])


def _working(root: Path) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    return (_paths(root, "diff", "--name-only", "-z", "--"), _paths(root, "diff", "--cached", "--name-only", "-z", "--"), _paths(root, "ls-files", "--others", "--exclude-standard", "-z", "--"))


def _environment() -> dict[str, str]:
    observed = {"python_implementation": platform.python_implementation(), "python_version": platform.python_version(), "numpy_version": np.__version__, "system": platform.system(), "sys_platform": sys.platform, "machine": platform.machine(), "byteorder": sys.byteorder}
    if observed != ENVIRONMENT:
        _fail("execution environment differs")
    return observed


def _expected_config() -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION, "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION, "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL, "owner_document": OWNER_DOCUMENT,
        "output_namespace": OUTPUT_NAMESPACE,
        "predecessor": {"commit": PREF1_COMMIT, "config_path": PREF1_CONFIG_PATH, "config_sha256": PREF1_CONFIG_SHA256, "result_path": PREF1_RESULT_PATH, "result_sha256": PREF1_RESULT_SHA256, "binder_path": PREF1_BINDER_PATH, "binder_sha256": PREF1_BINDER_SHA256, "owner_path": PREF1_OWNER_PATH, "owner_sha256": PREF1_OWNER_SHA256},
        "selection": {"method": METHOD, "member_key": MEMBER_KEY, "point_count": POINT_COUNT, "owned_row_count": OWNED_ROW_COUNT, "grid_spacing_hex": GRID_SPACING_HEX, "outer_radius_hex": OUTER_RADIUS_HEX, "retries": list(RETRIES), "failed_occurrences": [{"retry": retry, "channel": channel} for retry, channel in FAILED_OCCURRENCES]},
        "work_budget": {"polynomials_per_occurrence": POLYNOMIALS_PER_OCCURRENCE, "base_cubic_count": TOTAL_POLYNOMIALS, "component_polynomials_per_base_cubic": 3, "value_V_component_cubic_count": COMPONENT_CUBIC_COUNT, "slope_S_component_cubic_count": COMPONENT_CUBIC_COUNT, "complete_C_component_cubic_count": COMPONENT_CUBIC_COUNT, "executed_component_cubic_count": EXECUTED_COMPONENT_CUBIC_COUNT, "complete_C_maximum_candidates": COMPLETE_C_MAXIMUM_CANDIDATES, "value_V_maximum_candidates": VALUE_V_MAXIMUM_CANDIDATES, "slope_S_maximum_candidates": SLOPE_S_MAXIMUM_CANDIDATES, "total_VSC_maximum_candidates": TOTAL_VSC_MAXIMUM_CANDIDATES, "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH, "independent_refinement_bits": INDEPENDENT_REFINEMENT_BITS},
        "row_regions": [{"owned_start": start, "owned_end": end, "name": name} for start, end, name in ROW_REGIONS],
        "environment": dict(ENVIRONMENT), "committed_image": {"base_commit": BASE_COMMIT, "delta_paths": list(DELTA_PATHS), "allowed_untracked_paths": list(ALLOWED_UNTRACKED), "single_direct_successor_required": True},
        "scope": dict(SCOPE), "claims": dict(CLAIMS),
    }


def parse_config(raw: bytes) -> dict[str, object]:
    try:
        value = tomllib.loads(raw.decode())
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise LOC1AuthorityError("invalid LOC1 TOML") from exc
    if value != _expected_config():
        _fail("LOC1 config differs from frozen contract")
    return value


def _require_predecessor(root: Path) -> None:
    for path, digest in ((PREF1_CONFIG_PATH, PREF1_CONFIG_SHA256), (PREF1_RESULT_PATH, PREF1_RESULT_SHA256), (PREF1_BINDER_PATH, PREF1_BINDER_SHA256), (PREF1_OWNER_PATH, PREF1_OWNER_SHA256)):
        raw = _git(root, "show", f"{PREF1_COMMIT}:{path}")
        if sha256(raw).hexdigest() != digest:
            _fail("sealed PREF1 blob differs")


def _require_absent(root: Path) -> None:
    path = root / OUTPUT_NAMESPACE
    current = root
    for part in path.relative_to(root).parent.parts:
        current /= part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            return
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _fail("LOC1 output parent is unsafe")
    try:
        path.lstat()
    except FileNotFoundError:
        pass
    else:
        _fail("LOC1 output namespace already exists")
    if any(item.name.startswith(STAGING_PREFIX) for item in os.scandir(path.parent)):
        _fail("LOC1 staging namespace already exists")


def build_prelaunch(config_raw: bytes, root: Path) -> dict[str, object]:
    repository = root.resolve()
    parse_config(config_raw)
    _environment()
    if _head(repository) != BASE_COMMIT:
        _fail("prelaunch HEAD is not sealed PREF1")
    unstaged, staged, untracked = _working(repository)
    observed = tuple(sorted((*unstaged, *(item for item in untracked if item not in ALLOWED_UNTRACKED))))
    if staged or observed != DELTA_PATHS or any(item not in ALLOWED_UNTRACKED and item not in DELTA_PATHS for item in untracked):
        _fail("prospective LOC1 delta differs")
    for relative in DELTA_PATHS:
        _read(repository, relative)
    _require_predecessor(repository)
    _require_absent(repository)
    return {"schema_version": SCHEMA_VERSION, "artifact_id": ARTIFACT_ID, "classification": CLASSIFICATION, "gate_status": "pass", "source_config_sha256": sha256(config_raw).hexdigest(), "artifact_payload": {**_expected_config(), "output_namespace_absent_at_prelaunch": True, "PREF1_commit_bound": True, "diagnostic_executed": False}}


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    parse_config(config_raw)
    try:
        result = json.loads(result_raw, object_pairs_hook=_unique)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise LOC1AuthorityError("invalid LOC1 compact result") from exc
    expected = {"schema_version": SCHEMA_VERSION, "artifact_id": ARTIFACT_ID, "classification": CLASSIFICATION, "gate_status": "pass", "source_config_sha256": sha256(config_raw).hexdigest(), "artifact_payload": {**_expected_config(), "output_namespace_absent_at_prelaunch": True, "PREF1_commit_bound": True, "diagnostic_executed": False}}
    if result != expected or canonical_pretty(result) != result_raw:
        _fail("LOC1 compact result differs")
    return result


@dataclass(frozen=True, slots=True)
class LOC1Authority:
    authority_commit: str
    output_namespace: str = OUTPUT_NAMESPACE
    total_polynomials: int = TOTAL_POLYNOMIALS
    maximum_candidates: int = COMPLETE_C_MAXIMUM_CANDIDATES
    total_VSC_maximum_candidates: int = TOTAL_VSC_MAXIMUM_CANDIDATES
    execution_authorized: bool = True
    store_mutation_authorized: bool = False
    candidate_branches_authorized: bool = False


def authorize(root: Path, config_raw: bytes, result_raw: bytes, authority_commit: str) -> LOC1Authority:
    repository = root.resolve()
    if not _COMMIT.fullmatch(authority_commit):
        _fail("authority commit malformed")
    validate_compact(config_raw, result_raw)
    _environment()
    if _head(repository) != authority_commit or _parents(repository, authority_commit) != (BASE_COMMIT,):
        _fail("authority is not exact direct successor")
    unstaged, staged, untracked = _working(repository)
    if unstaged or staged or any(item not in ALLOWED_UNTRACKED for item in untracked):
        _fail("authority working image is not clean")
    delta = tuple(sorted(_paths(repository, "diff", "--name-only", "-z", BASE_COMMIT, authority_commit, "--")))
    if delta != DELTA_PATHS:
        _fail("committed LOC1 delta differs")
    if _git(repository, "show", f"{authority_commit}:{CONFIG_PATH}") != config_raw or _git(repository, "show", f"{authority_commit}:{RESULT_PATH}") != result_raw:
        _fail("committed LOC1 authority bytes differ")
    _require_predecessor(repository)
    return LOC1Authority(authority_commit)


__all__ = ["ARTIFACT_ID", "CONFIG_PATH", "RESULT_PATH", "OUTPUT_NAMESPACE", "STAGING_PREFIX", "BASE_COMMIT", "FAILED_OCCURRENCES", "TOTAL_POLYNOMIALS", "COMPONENT_CUBIC_COUNT", "EXECUTED_COMPONENT_CUBIC_COUNT", "MAXIMUM_CANDIDATES", "COMPLETE_C_MAXIMUM_CANDIDATES", "VALUE_V_MAXIMUM_CANDIDATES", "SLOPE_S_MAXIMUM_CANDIDATES", "TOTAL_VSC_MAXIMUM_CANDIDATES", "ROW_REGIONS", "LOC1Authority", "LOC1AuthorityError", "authorize", "build_prelaunch", "canonical_pretty", "parse_config", "validate_compact"]
