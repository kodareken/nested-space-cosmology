"""Prospective authority for the bounded TDG9 LOC2 diagnostic."""

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
ARTIFACT_ID = "FGC-1-TDG9-LOC2-FRZ1"
CLASSIFICATION = "premise_only_scale_invariant_independent_localization_authority"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-loc2-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg9-loc2-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-loc2-frz1.md"
OUTPUT_NAMESPACE = "runs/fgc-2-sf1/tdg9-loc2/rcv3-retries-3-5"
STAGING_PREFIX = ".rcv3-retries-3-5.stage-"
BASE_COMMIT = "14e26c43edbc88003d68cbcf24681bb090d98534"

ERR1_RESULT_PATH = "results/fgc-1-tdg9-loc1-err1.json"
ERR1_RESULT_SHA256 = "07f85b3692b22235380965a9cd6656f1dbdb4948386956d4769a58ea45ef0890"
ERR1_REPRODUCER_PATH = "scripts/reproduce_fgc_tdg9_loc1_err1.py"
ERR1_REPRODUCER_SHA256 = "7ad83e4ce8797f3ea88109259c252ecb1e602a1e52dbe03d6dfbdf63b8b18cfd"
SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SNAPSHOT_SHA256 = "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
EVALUATOR_PATH = "src/recursive_horizons/fgc/evolution/tdg9_local_extrema_independent_v2.py"
EVALUATOR_SHA256 = "750c38bb6bad4b37d946f66708055f3dc0835c807244e9c25c6e87e4e87e0511"
EVALUATOR_ID = "tdg9_loc2_endpoint_deflated_adaptive_discriminant_v2"
RUNNER_PATH = "scripts/run_fgc_tdg9_loc2.py"
RUNNER_SHA256 = "e3ad84f818fe33afce707bcbec1e21a99e0892b3d924536549e0ead740351122"

METHOD = "RK4"
MEMBER_KEY = "RK4-2049"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
GRID_SPACING_HEX = "0x1.0000000000000p-4"
OUTER_RADIUS_HEX = "0x1.0000000000000p+7"
FAILED_OCCURRENCES = (
    (3, "u:alpha"), (3, "u:R"),
    (4, "u:alpha"), (4, "u:lambda"), (4, "u:R"),
    (5, "u:alpha"), (5, "u:v"), (5, "u:lambda"), (5, "u:R"), (5, "q:R"),
)
BASE_CUBIC_COUNT = 122_640
COMPONENT_CUBIC_COUNT = 367_920
FAMILY_CANDIDATE_CEILING = 490_560
AGGREGATE_CANDIDATE_CEILING = 1_471_680
PRIMARY_REFINEMENT_DEPTH = 160
INITIAL_REFINEMENT_BITS = 192
REFINEMENT_SCHEDULE = "doubling_plus_exact_proof_cap"
GLOBAL_PROOF_BIT_CEILING = 9216

COEFFICIENT_DEPTH_AUDIT = {
    "counts": {"base_cubics": 122_640, "component_cubics": 367_920, "endpoint_deflated": 278_706, "positive_rational_square": 278_706, "positive_nonsquare": 89_207, "negative_discriminant": 7, "zero_discriminant": 0},
    "maxima": {"coefficient_numerator_bits": 95, "coefficient_denominator_bits": 427, "derivative_numerator_bits": 96, "derivative_denominator_bits": 427, "discriminant_numerator_bits": 190, "discriminant_denominator_bits": 851, "primitive_height_bits": 96, "proof_bits": 200},
    "maximum_proof_owner": [3, "u:alpha", "D12", "complete_C", 458],
    "initial_refinement_bits": INITIAL_REFINEMENT_BITS,
    "refinement_schedule": REFINEMENT_SCHEDULE,
    "global_proof_bit_ceiling": GLOBAL_PROOF_BIT_CEILING,
    "all_proof_caps_within_global_ceiling": True,
    "localization_executed": False,
}
ENVIRONMENT = {"python_implementation": "CPython", "python_version": "3.14.3", "numpy_version": "2.5.1", "system": "Darwin", "sys_platform": "darwin", "machine": "arm64", "byteorder": "little"}

DELTA_PATHS = tuple(sorted((
    "Makefile", "README.md", CONFIG_PATH, "docs/claim-ledger.md",
    "docs/fgc-runtime-matrix.md", "docs/fgc-tdg9-loc1-err1.md",
    OWNER_DOCUMENT, "docs/research-roadmap.md", "results/README.md",
    ERR1_RESULT_PATH, RESULT_PATH, "scripts/check_repo.py",
    ERR1_REPRODUCER_PATH, "scripts/reproduce_fgc_tdg9_loc2_frz1.py",
    RUNNER_PATH, EVALUATOR_PATH,
    "tests/test_check_repo_tdg9_loc2_frz1.py",
    "tests/test_fgc_tdg9_loc1_err1.py",
    "tests/test_fgc_tdg9_loc2_authority.py",
    "tests/test_fgc_tdg9_loc2_runner.py",
    "tests/test_fgc_tdg9_local_extrema_v2.py",
    "src/recursive_horizons/fgc/evolution/tdg9_loc2_authority.py",
)))
ALLOWED_UNTRACKED = (".qdrant-initialized",)
SCOPE = {"live_store_read_only": True, "compact_verifier_store_blind": True, "historical_threshold_changed": False, "ULP_used_as_tolerance": False, "fourth_width_authorized": False, "stage2_source_decomposition_authorized": False, "SSPRK3_comparator_authorized": False, "campaign_store_mutation_authorized": False, "PDE_state_commit_authorized": False, "continuation_authorized": False, "candidate_branches_authorized": False, "mechanism_result_authorized": False, "physical_result_authorized": False}
CLAIMS = {"ERR1_diagnosis_bound": True, "LOC2_execution_authorized": True, "LOC2_result_earned": False, "common_event_completed": False, "GR0_calibration_completed": False, "candidate_execution_authorized": False, "mechanism_result_earned": False, "physical_result_earned": False}

_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class LOC2AuthorityError(RuntimeError):
    pass


def _fail(message: str) -> NoReturn:
    raise LOC2AuthorityError(message)


def canonical_pretty(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("duplicate compact JSON key")
        answer[key] = value
    return answer


def _git(root: Path, *arguments: str) -> bytes:
    environment = dict(os.environ, LC_ALL="C", LANG="C", GIT_NO_REPLACE_OBJECTS="1", GIT_CONFIG_NOSYSTEM="1")
    try:
        return subprocess.run(("git", "--no-replace-objects", "--no-optional-locks", "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false", *arguments), cwd=root, env=environment, check=True, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    except subprocess.CalledProcessError as exc:
        raise LOC2AuthorityError("Git image differs") from exc


def _paths(root: Path, *arguments: str) -> tuple[str, ...]:
    raw = _git(root, *arguments)
    if not raw:
        return ()
    chunks = raw.split(b"\0")
    if chunks[-1] != b"":
        _fail("Git path stream differs")
    return tuple(item.decode() for item in chunks[:-1])


def _read(root: Path, relative: str) -> bytes:
    candidate = Path(relative)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        _fail("unsafe authority path")
    path = root / candidate
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        _fail("unsafe authority leaf")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        active = os.fstat(descriptor)
        raw = b""
        while block := os.read(descriptor, 1 << 20):
            raw += block
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    if identity != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns) or identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) or len(raw) != before.st_size:
        _fail("authority leaf changed")
    return raw


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
        "predecessor": {"commit": BASE_COMMIT, "ERR1_result_path": ERR1_RESULT_PATH, "ERR1_result_sha256": ERR1_RESULT_SHA256, "ERR1_reproducer_path": ERR1_REPRODUCER_PATH, "ERR1_reproducer_sha256": ERR1_REPRODUCER_SHA256, "sealed_store_leaf_count": SEALED_STORE_LEAF_COUNT, "sealed_store_snapshot_sha256": SEALED_STORE_SNAPSHOT_SHA256},
        "implementation": {"evaluator_path": EVALUATOR_PATH, "evaluator_sha256": EVALUATOR_SHA256, "evaluator_id": EVALUATOR_ID, "runner_path": RUNNER_PATH, "runner_sha256": RUNNER_SHA256, "initial_refinement_bits": INITIAL_REFINEMENT_BITS, "refinement_schedule": REFINEMENT_SCHEDULE, "proof_bit_formula": "2*h+8", "global_proof_bit_ceiling": GLOBAL_PROOF_BIT_CEILING, "typed_root_stops": ["root_proof_cap_exceeded", "root_location_inconclusive", "root_separation_inconclusive", "root_count_invariant_failure"], "boundary_clipping_authorized": False},
        "selection": {"method": METHOD, "member_key": MEMBER_KEY, "point_count": POINT_COUNT, "owned_row_count": OWNED_ROW_COUNT, "grid_spacing_hex": GRID_SPACING_HEX, "outer_radius_hex": OUTER_RADIUS_HEX, "retries": [3, 4, 5], "failed_occurrences": [{"retry": retry, "channel": channel} for retry, channel in FAILED_OCCURRENCES]},
        "work_budget": {"base_cubic_count": BASE_CUBIC_COUNT, "component_cubic_count": COMPONENT_CUBIC_COUNT, "value_V_cubic_count": BASE_CUBIC_COUNT, "slope_S_cubic_count": BASE_CUBIC_COUNT, "complete_C_cubic_count": BASE_CUBIC_COUNT, "family_candidate_ceiling": FAMILY_CANDIDATE_CEILING, "aggregate_candidate_ceiling": AGGREGATE_CANDIDATE_CEILING, "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH},
        "coefficient_depth_audit": COEFFICIENT_DEPTH_AUDIT,
        "environment": ENVIRONMENT, "scope": SCOPE, "claims": CLAIMS,
    }


def parse_config(raw: bytes) -> dict[str, object]:
    try:
        value = tomllib.loads(raw.decode())
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise LOC2AuthorityError("invalid LOC2 TOML") from exc
    if value != _expected_config():
        _fail("LOC2 config differs from frozen contract")
    return value


def _require_absent(root: Path) -> None:
    target = root / OUTPUT_NAMESPACE
    current = root
    for part in target.relative_to(root).parent.parts:
        current /= part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            return
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _fail("LOC2 output parent is unsafe")
    try:
        target.lstat()
    except FileNotFoundError:
        pass
    else:
        _fail("LOC2 output namespace already exists")
    if any(item.name.startswith(STAGING_PREFIX) for item in os.scandir(target.parent)):
        _fail("LOC2 staging namespace already exists")


def _verify_committed_bindings(root: Path, authority_commit: str) -> None:
    for path, digest in (
        (ERR1_RESULT_PATH, ERR1_RESULT_SHA256),
        (ERR1_REPRODUCER_PATH, ERR1_REPRODUCER_SHA256),
        (EVALUATOR_PATH, EVALUATOR_SHA256),
        (RUNNER_PATH, RUNNER_SHA256),
    ):
        raw = _git(root, "show", f"{authority_commit}:{path}")
        if sha256(raw).hexdigest() != digest:
            _fail(f"committed LOC2 binding differs: {path}")


def build_prelaunch(config_raw: bytes, root: Path, coefficient_audit: dict[str, object]) -> dict[str, object]:
    repository = root.resolve()
    parse_config(config_raw)
    _environment()
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != BASE_COMMIT:
        _fail("prelaunch HEAD is not sealed LOC1 authority")
    unstaged, staged, untracked = _working(repository)
    observed = tuple(sorted((*unstaged, *(item for item in untracked if item not in ALLOWED_UNTRACKED))))
    if staged or observed != DELTA_PATHS:
        _fail("prospective LOC2 delta differs")
    for relative in DELTA_PATHS:
        _read(repository, relative)
    for path, digest in ((ERR1_RESULT_PATH, ERR1_RESULT_SHA256), (ERR1_REPRODUCER_PATH, ERR1_REPRODUCER_SHA256), (EVALUATOR_PATH, EVALUATOR_SHA256), (RUNNER_PATH, RUNNER_SHA256)):
        if sha256(_read(repository, path)).hexdigest() != digest:
            _fail("LOC2 implementation hash differs")
    if coefficient_audit != COEFFICIENT_DEPTH_AUDIT:
        _fail("coefficient-depth audit differs")
    _require_absent(repository)
    return {"schema_version": SCHEMA_VERSION, "artifact_id": ARTIFACT_ID, "classification": CLASSIFICATION, "gate_status": "pass", "source_config_sha256": sha256(config_raw).hexdigest(), "artifact_payload": {**_expected_config(), "coefficient_depth_audit_reproduced": True, "output_namespace_absent_at_prelaunch": True, "diagnostic_executed": False}}


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    parse_config(config_raw)
    try:
        result = json.loads(result_raw, object_pairs_hook=_unique)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise LOC2AuthorityError("invalid LOC2 compact result") from exc
    expected = {"schema_version": SCHEMA_VERSION, "artifact_id": ARTIFACT_ID, "classification": CLASSIFICATION, "gate_status": "pass", "source_config_sha256": sha256(config_raw).hexdigest(), "artifact_payload": {**_expected_config(), "coefficient_depth_audit_reproduced": True, "output_namespace_absent_at_prelaunch": True, "diagnostic_executed": False}}
    if result != expected or result_raw != canonical_pretty(result):
        _fail("LOC2 compact result differs")
    return result


@dataclass(frozen=True, slots=True)
class LOC2Authority:
    authority_commit: str
    base_cubic_count: int = BASE_CUBIC_COUNT
    execution_authorized: bool = True
    store_mutation_authorized: bool = False
    candidate_branches_authorized: bool = False


def authorize(root: Path, config_raw: bytes, result_raw: bytes, authority_commit: str) -> LOC2Authority:
    repository = root.resolve()
    if not _COMMIT.fullmatch(authority_commit):
        _fail("authority commit malformed")
    validate_compact(config_raw, result_raw)
    _environment()
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    parents = _git(repository, "rev-list", "--parents", "-n", "1", authority_commit).decode().split()[1:]
    if head != authority_commit or parents != [BASE_COMMIT]:
        _fail("authority is not exact direct successor")
    unstaged, staged, untracked = _working(repository)
    if unstaged or staged or any(item not in ALLOWED_UNTRACKED for item in untracked):
        _fail("authority working image is not clean")
    delta = tuple(sorted(_paths(repository, "diff", "--name-only", "-z", BASE_COMMIT, authority_commit, "--")))
    if delta != DELTA_PATHS:
        _fail("committed LOC2 delta differs")
    _verify_committed_bindings(repository, authority_commit)
    _require_absent(repository)
    return LOC2Authority(authority_commit)


__all__ = ["ARTIFACT_ID", "CONFIG_PATH", "RESULT_PATH", "OUTPUT_NAMESPACE", "STAGING_PREFIX", "BASE_COMMIT", "ERR1_RESULT_SHA256", "SEALED_STORE_LEAF_COUNT", "SEALED_STORE_SNAPSHOT_SHA256", "FAILED_OCCURRENCES", "BASE_CUBIC_COUNT", "COMPONENT_CUBIC_COUNT", "FAMILY_CANDIDATE_CEILING", "AGGREGATE_CANDIDATE_CEILING", "OWNED_ROW_COUNT", "PRIMARY_REFINEMENT_DEPTH", "LOC2Authority", "LOC2AuthorityError", "authorize", "build_prelaunch", "canonical_pretty", "parse_config", "validate_compact"]
