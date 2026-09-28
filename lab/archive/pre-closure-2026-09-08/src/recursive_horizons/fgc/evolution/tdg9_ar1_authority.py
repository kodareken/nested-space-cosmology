"""Prospective authority for one read-only TDG9 exact-arithmetic audit.

AR1 may deterministically replay exactly three already-persisted GR-0 TDG6
rejections (retry depths 3, 4, and 5) and apply the two implementation-only
exact Hermite evaluators.  It cannot commit a PDE state, alter the RCV3 store,
open a candidate branch, select a runtime remedy, or authorize continuation.

The one-time prelaunch builder authenticates the immutable PREF2 terminal and
observes that the raw AR1 output namespace is absent.  Compact verification is
deliberately store-blind after that observation is committed.
"""
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
from typing import Any, Mapping, NoReturn

import numpy as np

from . import tdg8_rcv3_pref2_binder as pref2
from . import tdg9_exact_temporal_arithmetic as exact_primary
from . import tdg9_exact_temporal_arithmetic_independent as exact_independent


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-AR1-AUTH1"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_read_only_three_width_exact_arithmetic_authority"
GATE_STATUS = "pass"
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ar1-auth1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ar1-auth1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ar1-auth1.md"
OUTPUT_NAMESPACE = "runs/fgc-2-sf1/tdg9-ar1/rcv3-retries-3-5"
STAGING_PREFIX = ".rcv3-retries-3-5.stage-"

PYTHON_IMPLEMENTATION = "CPython"
PYTHON_VERSION = "3.14.3"
NUMPY_VERSION = "2.5.1"
SYSTEM = "Darwin"
SYS_PLATFORM = "darwin"
MACHINE = "arm64"
BYTEORDER = "little"

CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"
CAMPAIGN_PLAN_SHA256 = "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"
CAMPAIGN_AUTHORIZATION_COMMIT = "6df67967e6e0bc63eca4dc6951a60d3787368f26"
BRANCH = "GR-0"
AMPLITUDE = "3"
EVENT = 23
TARGET_RATIONAL = "3/2"
TARGET_BINARY64_HEX = "0x1.8000000000000p+0"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044

PREF2_COMMIT = "8c3a3bbe3476284ed46b50bc363eec7cf40405b3"
PREF2_CONFIG_PATH = pref2.CONFIG_PATH
PREF2_CONFIG_SHA256 = "9957ef6901dc2e332acd07a5812d975e08651d3ea7d80a1b19e9469c89c79f72"
PREF2_RESULT_PATH = pref2.RESULT_PATH
PREF2_RESULT_SHA256 = "7e676cb40fd83deeff0912603641fd93ec139f274287b1d092aa8eccd2c1231f"
PREF2_BINDER_PATH = "src/recursive_horizons/fgc/evolution/tdg8_rcv3_pref2_binder.py"
PREF2_BINDER_SHA256 = "4fedf666dfe5dfc34e17a2a55495b0109fde36d6a90747f4dc5a48921262ea4a"
PREF2_STORE_MANIFEST_SHA256 = "46b57bc38f3bbfbec70c4b6b11a089cb799dc73135af0d9dac0e3dc289b952d4"
PREF2_STORE_LEAF_COUNT = 115
PREF2_STORE_BYTE_COUNT = 25_734_805
PREF2_STORE_PATH = pref2.STORE_PATH
PREF2_PHYSICAL_STATE_SHA256 = pref2.PHYSICAL_STATE_SHA256

EVALUATOR_COMMIT = "3e37e3baa5a7cf862c4c0c004c7f43dcfa1313b2"
PRIMARY_EVALUATOR_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg9_exact_temporal_arithmetic.py"
)
PRIMARY_EVALUATOR_SHA256 = "2c9a0c1a9c6560ed787a605ef508ede7471953ceb222527889c5622ea6d8e82d"
INDEPENDENT_EVALUATOR_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg9_exact_temporal_arithmetic_independent.py"
)
INDEPENDENT_EVALUATOR_SHA256 = "851c1c69daa7366c90798b7fe7ba11e2e5150c06450ae0fc209a926cc7e62850"
EVALUATOR_TEST_PATH = "tests/test_fgc_tdg9_exact_temporal_arithmetic.py"
EVALUATOR_TEST_SHA256 = "7d44f1bc4c5929520bd0ef13ca23294aca39e7c9f34cd72dc620c8c78f6445ec"

METHOD = "RK4"
MEMBER_KEY = "RK4-2049"
CHANNEL_ORDER = tuple(
    f"{block}:{field}"
    for block in ("u", "p", "q")
    for field in ("alpha", "v", "lambda", "R", "phi", "chi")
)
ORDER_THRESHOLD = "3/2"
ORDER_SQUARED_MULTIPLIER = 8
MAX_DEPTH = 64
MAX_NODES_PER_CHANNEL_PER_EVALUATOR = 4_194_304
MAX_EVALUATOR_CALLS = 108

ENVIRONMENT = {
    "python_implementation": PYTHON_IMPLEMENTATION,
    "python_version": PYTHON_VERSION,
    "numpy_version": NUMPY_VERSION,
    "system": SYSTEM,
    "sys_platform": SYS_PLATFORM,
    "machine": MACHINE,
    "byteorder": BYTEORDER,
}
SELECTION = {
    "protocol": TARGET_PROTOCOL,
    "campaign_id": CAMPAIGN_ID,
    "campaign_plan_sha256": CAMPAIGN_PLAN_SHA256,
    "campaign_authorization_commit": CAMPAIGN_AUTHORIZATION_COMMIT,
    "branch": BRANCH,
    "amplitude": AMPLITUDE,
    "method": METHOD,
    "member_key": MEMBER_KEY,
    "point_count": POINT_COUNT,
    "owned_row_count": OWNED_ROW_COUNT,
    "event": EVENT,
    "target_rational": TARGET_RATIONAL,
    "target_binary64_hex": TARGET_BINARY64_HEX,
}

AUTHORITY_DELTA_PATHS = (
    "Makefile",
    "README.md",
    CONFIG_PATH,
    "docs/claim-ledger.md",
    "docs/fgc-runtime-matrix.md",
    OWNER_DOCUMENT,
    "docs/research-roadmap.md",
    "paper/fgc-local-defocusing/README.md",
    "results/README.md",
    RESULT_PATH,
    "scripts/check_repo.py",
    "scripts/reproduce_fgc_tdg9_ar1_auth1.py",
    "scripts/run_fgc_tdg9_ar1.py",
    "src/recursive_horizons/fgc/evolution/tdg9_ar1_authority.py",
    "tests/test_check_repo_tdg9_ar1_auth1.py",
    "tests/test_fgc_tdg9_ar1_authority.py",
    "tests/test_fgc_tdg9_ar1_runner.py",
)
ALLOWED_UNTRACKED_PATHS = (".qdrant-initialized",)
COMMITTED_IMAGE = {
    "base_commit": EVALUATOR_COMMIT,
    "authority_delta_paths": list(AUTHORITY_DELTA_PATHS),
    "allowed_untracked_paths": list(ALLOWED_UNTRACKED_PATHS),
    "prelaunch_requires_base_commit_at_HEAD": True,
    "execution_requires_authority_commit_at_clean_tracked_HEAD": True,
    "authority_commit_must_be_single_direct_successor_of_base": True,
    "predecessor_commits_must_be_ancestors": True,
    "replace_refs_forbidden": True,
    "exact_delta_locks_unlisted_transitive_modules": True,
}

REPLAYS = (
    {
        "retry": 3,
        "prior_retry_count": 2,
        "predecessor_generation": 9,
        "checkpoint_sha256": "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56",
        "checkpoint_raw_sha256": "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846",
        "journal_sequence": 11,
        "journal_sha256": "341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a",
        "journal_raw_sha256": "e4f9ca42f2d014d4bdbd875a0f24a9bbb9c04159fb0d9e0e7de8d0d72db61fdf",
        "attempted_width_hex": "0x1.aaa9612df8000p-11",
    },
    {
        "retry": 4,
        "prior_retry_count": 3,
        "predecessor_generation": 10,
        "checkpoint_sha256": "6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187",
        "checkpoint_raw_sha256": "9abb59999809d22354d6b5f810bc592cd2673e5e4d354d0ec58e9ad53522dc4b",
        "journal_sequence": 13,
        "journal_sha256": "1bd21553f3b6e0af613997da094a46d4dcd966a0e4346692bae62b05046862ca",
        "journal_raw_sha256": "26542403a17b05664d5a447b88860186d5f554b1f4d555624493ed40b941694c",
        "attempted_width_hex": "0x1.aaa9612df8000p-12",
    },
    {
        "retry": 5,
        "prior_retry_count": 4,
        "predecessor_generation": 11,
        "checkpoint_sha256": "11ec8a80d300de72075661a97b616aecead8ec78208bb541b93e62ca4c8124c8",
        "checkpoint_raw_sha256": "786bdcb835e043c65ddd33fad3a2af0e3d33e768dc26a9f9f00d2757b0da3ea5",
        "journal_sequence": 15,
        "journal_sha256": "475014d918952adef032d71b81073773364a37fa795491e739aae609f77e8537",
        "journal_raw_sha256": "1997342415f74f330ce6b45ca6fd550ad0e5e674b85148fea5047ff82ac3f152",
        "attempted_width_hex": "0x1.aaa9612df0000p-13",
    },
)

SCOPE = {
    "live_store_read_only": True,
    "compact_verifier_store_blind": True,
    "authenticated_deterministic_replay_only": True,
    "exact_three_retries_only": True,
    "historical_proposal_arrays_persisted": False,
    "PDE_state_commit_authorized": False,
    "campaign_store_mutation_authorized": False,
    "automatic_resource_escalation_authorized": False,
    "fourth_retry_width_authorized": False,
    "continuation_authorized": False,
    "candidate_branches_authorized": False,
}
CLAIMS = {
    "PREF2_terminal_bound": True,
    "exact_evaluator_implementation_bound": True,
    "read_only_AR1_execution_authorized": True,
    "raw_AR1_result_earned": False,
    "AR1_binder_completed": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "mechanism_result_earned": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}
NONCLAIMS = (
    "AR1 certifies only the finite-dimensional Hermite diagnostic reconstructed from authenticated deterministic binary64 replays; historical proposal arrays were not persisted.",
    "AR1 does not certify exact PDE history, continuum convergence, an asymptotic regime, GR-0 calibration, candidate dynamics, a mechanism, or physics.",
    "The historical TDG6 rejection records and threshold remain unchanged; AR1 neither repairs nor replaces the runtime gate.",
    "No AR1 outcome authorizes a campaign-state commit, continuation, a fourth width, automatic resource escalation, SGB-L, FGC-QR, retained-EFT evolution, or a transition claim.",
)

IMPLEMENTATION_INVENTORY = (
    ("authority", "src/recursive_horizons/fgc/evolution/tdg9_ar1_authority.py"),
    ("runner", "scripts/run_fgc_tdg9_ar1.py"),
    ("reproducer", "scripts/reproduce_fgc_tdg9_ar1_auth1.py"),
    ("owner_document", OWNER_DOCUMENT),
    ("exact_evaluator", PRIMARY_EVALUATOR_PATH),
    ("exact_evaluator", INDEPENDENT_EVALUATOR_PATH),
    ("predecessor_binder", PREF2_BINDER_PATH),
    ("store", "src/recursive_horizons/fgc/evolution/hlt16_campaign_store.py"),
    ("store", "src/recursive_horizons/fgc/evolution/hlt16_campaign_schema.py"),
    ("store", "src/recursive_horizons/fgc/evolution/hlt16_member_codec.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_lifecycle.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/proto15_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/proto19_gr0_static_factory.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_design.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg7_binary64_subdivision_lattice.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg8_persisted_retry_replay.py"),
    ("publication", "src/recursive_horizons/fgc/evolution/proto17_hlt15_runtime.py"),
    ("numerical", "src/recursive_horizons/fgc/evolution/numerical_engine.py"),
    ("numerical", "src/recursive_horizons/fgc/evolution/proto5_runtime.py"),
    ("numerical", "src/recursive_horizons/fgc/evolution/proto6_runtime.py"),
    ("numerical", "src/recursive_horizons/fgc/evolution/proto7_runtime.py"),
    ("source", "src/recursive_horizons/fgc/evolution/proto11_runtime.py"),
    ("source", "src/recursive_horizons/fgc/evolution/proto12_runtime.py"),
    ("source", "src/recursive_horizons/fgc/evolution/src4_vectorized_reference_source.py"),
    ("source", "src/recursive_horizons/fgc/evolution/vectorized_source.py"),
)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_MAX_CONFIG_BYTES = 2 * 1024 * 1024
_MAX_RESULT_BYTES = 8 * 1024 * 1024


class TDG9AR1AuthorityError(RuntimeError):
    """The prospective exact-arithmetic authority differs from its freeze."""


@dataclass(frozen=True, slots=True)
class TDG9AR1Authority:
    authority_commit: str
    pref2_commit: str
    evaluator_commit: str
    store_manifest_sha256: str
    output_namespace: str
    retries: tuple[int, int, int]
    retry_width_hex: tuple[str, str, str]
    max_depth: int
    max_nodes_per_channel_per_evaluator: int
    max_evaluator_calls: int
    channel_order: tuple[str, ...]
    environment: tuple[tuple[str, str], ...]
    selection: tuple[tuple[str, object], ...]
    authority_delta_paths: tuple[str, ...]
    config_sha256: str
    result_sha256: str
    implementation_inventory: tuple[tuple[str, str, str], ...]
    diagnostic_execution_authorized: bool = True
    campaign_state_mutation_authorized: bool = False
    continuation_authorized: bool = False
    candidate_branches_authorized: bool = False


def canonical(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def canonical_pretty(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _fail(message: str) -> NoReturn:
    raise TDG9AR1AuthorityError(message)


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _observe_environment() -> dict[str, str]:
    observed = {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "system": platform.system(),
        "sys_platform": sys.platform,
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
    }
    if observed != ENVIRONMENT:
        _fail("AR1 execution environment differs from its freeze")
    return observed


def _git_environment() -> dict[str, str]:
    path = os.environ.get("PATH")
    if not path:
        _fail("AR1 Git PATH is absent")
    return {
        "PATH": path,
        "LC_ALL": "C",
        "LANG": "C",
        "GIT_NO_REPLACE_OBJECTS": "1",
        "GIT_CONFIG_NOSYSTEM": "1",
    }


def _read_regular(root: Path, relative: str, maximum: int) -> bytes:
    candidate = Path(relative)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        _fail("authority path is unsafe")
    path = root / candidate
    try:
        before = path.lstat()
    except OSError as exc:
        raise TDG9AR1AuthorityError("authority leaf is absent") from exc
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_size > maximum:
        _fail("authority leaf is unsafe")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        active = os.fstat(descriptor)
        identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        if identity != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns):
            _fail("authority leaf raced")
        chunks: list[bytes] = []
        size = 0
        while size <= maximum:
            block = os.read(descriptor, min(1 << 20, maximum + 1 - size))
            if not block:
                break
            chunks.append(block)
            size += len(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        _fail("authority leaf changed during read")
    raw = b"".join(chunks)
    if len(raw) != before.st_size or len(raw) > maximum:
        _fail("authority leaf size differs")
    return raw


def _expected_predecessor() -> dict[str, object]:
    return {
        "pref2_commit": PREF2_COMMIT,
        "pref2_config_path": PREF2_CONFIG_PATH,
        "pref2_config_sha256": PREF2_CONFIG_SHA256,
        "pref2_result_path": PREF2_RESULT_PATH,
        "pref2_result_sha256": PREF2_RESULT_SHA256,
        "pref2_binder_path": PREF2_BINDER_PATH,
        "pref2_binder_sha256": PREF2_BINDER_SHA256,
        "store_path": PREF2_STORE_PATH,
        "store_manifest_sha256": PREF2_STORE_MANIFEST_SHA256,
        "store_leaf_count": PREF2_STORE_LEAF_COUNT,
        "store_byte_count": PREF2_STORE_BYTE_COUNT,
        "accepted_physical_state_sha256": PREF2_PHYSICAL_STATE_SHA256,
    }


def _expected_exact_arithmetic() -> dict[str, object]:
    return {
        "evaluator_commit": EVALUATOR_COMMIT,
        "primary_path": PRIMARY_EVALUATOR_PATH,
        "primary_sha256": PRIMARY_EVALUATOR_SHA256,
        "primary_id": exact_primary.TDG9_EXACT_ARITHMETIC_EVALUATOR_ID,
        "independent_path": INDEPENDENT_EVALUATOR_PATH,
        "independent_sha256": INDEPENDENT_EVALUATOR_SHA256,
        "independent_id": exact_independent.TDG9_INDEPENDENT_EVALUATOR_ID,
        "test_path": EVALUATOR_TEST_PATH,
        "test_sha256": EVALUATOR_TEST_SHA256,
        "order_threshold": ORDER_THRESHOLD,
        "squared_multiplier": ORDER_SQUARED_MULTIPLIER,
        "max_depth": MAX_DEPTH,
        "max_nodes_per_channel_per_evaluator": MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
        "max_evaluator_calls": MAX_EVALUATOR_CALLS,
        "automatic_resource_escalation": False,
    }


def _parse_config(raw: bytes) -> dict[str, Any]:
    if len(raw) > _MAX_CONFIG_BYTES:
        _fail("AR1 config is too large")
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG9AR1AuthorityError("AR1 config is invalid TOML") from exc
    required = {
        "schema_version", "artifact_id", "classification", "project_version",
        "target_protocol", "owner_document", "output_namespace", "predecessor",
        "exact_arithmetic", "environment", "selection", "committed_image",
        "replays", "channel_order", "scope", "claims", "implementation_inventory",
    }
    if set(value) != required:
        _fail("AR1 config top-level keys differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["classification"] != CLASSIFICATION
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["owner_document"] != OWNER_DOCUMENT
        or value["output_namespace"] != OUTPUT_NAMESPACE
        or value["predecessor"] != _expected_predecessor()
        or value["exact_arithmetic"] != _expected_exact_arithmetic()
        or value["environment"] != ENVIRONMENT
        or value["selection"] != SELECTION
        or value["committed_image"] != COMMITTED_IMAGE
        or value["replays"] != list(REPLAYS)
        or value["channel_order"] != list(CHANNEL_ORDER)
        or value["scope"] != SCOPE
        or value["claims"] != CLAIMS
    ):
        _fail("AR1 config differs from the frozen contract")
    inventory = value["implementation_inventory"]
    if not isinstance(inventory, list) or len(inventory) != len(IMPLEMENTATION_INVENTORY):
        _fail("AR1 implementation inventory length differs")
    expected_pairs = list(IMPLEMENTATION_INVENTORY)
    observed_pairs: list[tuple[str, str]] = []
    for item in inventory:
        if not isinstance(item, dict) or set(item) != {"role", "path", "sha256"}:
            _fail("AR1 implementation inventory entry differs")
        if not isinstance(item["sha256"], str) or not _SHA256.fullmatch(item["sha256"]):
            _fail("AR1 implementation inventory digest differs")
        observed_pairs.append((item["role"], item["path"]))
    if observed_pairs != expected_pairs:
        _fail("AR1 implementation inventory order differs")
    return value


def _inventory(root: Path, config: Mapping[str, Any]) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    for item in config["implementation_inventory"]:
        raw = _read_regular(root, item["path"], 16 * 1024 * 1024)
        digest = _sha(raw)
        if digest != item["sha256"]:
            _fail(f"AR1 implementation inventory drifted: {item['path']}")
        result.append({"role": item["role"], "path": item["path"], "sha256": digest})
    return result


def _git(root: Path, *args: str) -> bytes:
    try:
        return subprocess.run(
            (
                "git", "--no-replace-objects", "--no-optional-locks",
                "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false",
                *args,
            ),
            cwd=root,
            env=_git_environment(),
            stdin=subprocess.DEVNULL,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except subprocess.CalledProcessError as exc:
        raise TDG9AR1AuthorityError("AR1 Git image differs") from exc


def _git_paths(root: Path, *args: str) -> tuple[str, ...]:
    raw = _git(root, *args)
    if not raw:
        return ()
    entries = raw.split(b"\0")
    if entries[-1] != b"":
        _fail("AR1 Git path stream is not NUL terminated")
    result: list[str] = []
    for item in entries[:-1]:
        try:
            path = item.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise TDG9AR1AuthorityError("AR1 Git path is not UTF-8") from exc
        candidate = Path(path)
        if (
            not path
            or candidate.is_absolute()
            or any(part in {"", ".", ".."} for part in candidate.parts)
            or candidate.as_posix() != path
        ):
            _fail("AR1 Git path is unsafe")
        result.append(path)
    if len(set(result)) != len(result):
        _fail("AR1 Git path stream contains duplicates")
    return tuple(result)


def _require_no_replace_refs(root: Path) -> None:
    if _git(root, "for-each-ref", "--format=%(refname)", "refs/replace").strip():
        _fail("AR1 Git replace refs are forbidden")


def _head(root: Path) -> str:
    try:
        value = _git(root, "rev-parse", "--verify", "HEAD^{commit}").decode("ascii").strip()
    except UnicodeDecodeError as exc:
        raise TDG9AR1AuthorityError("AR1 HEAD is not ASCII") from exc
    if not _COMMIT.fullmatch(value):
        _fail("AR1 HEAD identity differs")
    return value


def _require_ancestor(root: Path, ancestor: str, descendant: str) -> None:
    try:
        subprocess.run(
            (
                "git", "--no-replace-objects", "--no-optional-locks",
                "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false",
                "merge-base", "--is-ancestor", ancestor, descendant,
            ),
            cwd=root,
            env=_git_environment(),
            stdin=subprocess.DEVNULL,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        _fail("AR1 predecessor commit is not an ancestor")


def _commit_parents(root: Path, commit: str) -> tuple[str, ...]:
    try:
        fields = _git(root, "rev-list", "--parents", "-n", "1", commit).decode(
            "ascii"
        ).split()
    except UnicodeDecodeError as exc:
        raise TDG9AR1AuthorityError("AR1 commit ancestry is not ASCII") from exc
    if not fields or fields[0] != commit or any(
        not _COMMIT.fullmatch(value) for value in fields
    ):
        _fail("AR1 commit ancestry differs")
    return tuple(fields[1:])


def _working_image(root: Path) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    unstaged = _git_paths(root, "diff", "--name-only", "-z", "--")
    staged = _git_paths(root, "diff", "--cached", "--name-only", "-z", "--")
    untracked = _git_paths(
        root, "ls-files", "--others", "--exclude-standard", "-z", "--"
    )
    return unstaged, staged, untracked


def _require_prospective_image(root: Path) -> None:
    _require_no_replace_refs(root)
    if _head(root) != EVALUATOR_COMMIT:
        _fail("AR1 prelaunch HEAD is not the exact evaluator commit")
    unstaged, staged, untracked = _working_image(root)
    if staged:
        _fail("AR1 prelaunch index is not clean")
    task_untracked = tuple(
        path for path in untracked if path not in ALLOWED_UNTRACKED_PATHS
    )
    foreign_untracked = tuple(
        path
        for path in untracked
        if path in ALLOWED_UNTRACKED_PATHS and path != ".qdrant-initialized"
    )
    observed = tuple(sorted((*unstaged, *task_untracked)))
    if (
        foreign_untracked
        or len(set((*unstaged, *task_untracked))) != len(observed)
        or observed != AUTHORITY_DELTA_PATHS
    ):
        _fail("AR1 prospective authority delta differs")
    for path in AUTHORITY_DELTA_PATHS:
        _read_regular(root, path, 16 * 1024 * 1024)


def _require_runtime_image(root: Path, authority_commit: str) -> None:
    _require_no_replace_refs(root)
    if _head(root) != authority_commit:
        _fail("AR1 authority commit is not exact HEAD")
    if _commit_parents(root, authority_commit) != (EVALUATOR_COMMIT,):
        _fail("AR1 authority commit is not the single direct successor")
    unstaged, staged, untracked = _working_image(root)
    if unstaged or staged or any(
        path not in ALLOWED_UNTRACKED_PATHS for path in untracked
    ):
        _fail("AR1 runtime tracked image is not clean")
    _require_ancestor(root, PREF2_COMMIT, authority_commit)
    _require_ancestor(root, EVALUATOR_COMMIT, authority_commit)
    delta = tuple(
        sorted(
            _git_paths(
                root,
                "diff", "--name-only", "-z", EVALUATOR_COMMIT,
                authority_commit, "--",
            )
        )
    )
    if delta != AUTHORITY_DELTA_PATHS:
        _fail("AR1 committed authority delta differs")


def _require_commit_blob(root: Path, commit: str, path: str, digest: str) -> None:
    if not _COMMIT.fullmatch(commit) or not _SHA256.fullmatch(digest):
        _fail("AR1 commit/blob identity is malformed")
    raw = _git(root, "show", f"{commit}:{path}")
    if _sha(raw) != digest:
        _fail(f"AR1 committed predecessor blob differs: {path}")


def _validate_predecessors(root: Path) -> None:
    for path, digest in (
        (PREF2_CONFIG_PATH, PREF2_CONFIG_SHA256),
        (PREF2_RESULT_PATH, PREF2_RESULT_SHA256),
        (PREF2_BINDER_PATH, PREF2_BINDER_SHA256),
    ):
        _require_commit_blob(root, PREF2_COMMIT, path, digest)
    for path, digest in (
        (PRIMARY_EVALUATOR_PATH, PRIMARY_EVALUATOR_SHA256),
        (INDEPENDENT_EVALUATOR_PATH, INDEPENDENT_EVALUATOR_SHA256),
        (EVALUATOR_TEST_PATH, EVALUATOR_TEST_SHA256),
    ):
        _require_commit_blob(root, EVALUATOR_COMMIT, path, digest)


def _require_namespace_absent(root: Path) -> None:
    candidate = root / OUTPUT_NAMESPACE
    parent = candidate.parent
    current = root
    for part in Path(OUTPUT_NAMESPACE).parent.parts:
        current = current / part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            return
        except OSError as exc:
            raise TDG9AR1AuthorityError("AR1 output namespace cannot be inspected") from exc
        if stat.S_ISLNK(metadata.st_mode):
            _fail("AR1 output namespace is redirected")
        if not stat.S_ISDIR(metadata.st_mode):
            _fail("AR1 output namespace parent is not a directory")
    try:
        with os.scandir(parent) as entries:
            names = tuple(sorted(item.name for item in entries))
    except OSError as exc:
        raise TDG9AR1AuthorityError("AR1 output parent cannot be inspected") from exc
    if any(name.startswith(STAGING_PREFIX) for name in names):
        _fail("AR1 private output staging already exists")
    try:
        metadata = candidate.lstat()
    except FileNotFoundError:
        return
    except OSError as exc:
        raise TDG9AR1AuthorityError("AR1 output namespace cannot be inspected") from exc
    if stat.S_ISLNK(metadata.st_mode):
        _fail("AR1 output namespace is redirected")
    _fail("AR1 output namespace already exists")


def _authenticate_pref2_live(root: Path) -> None:
    config_raw = _read_regular(root, PREF2_CONFIG_PATH, _MAX_CONFIG_BYTES)
    result_raw = _read_regular(root, PREF2_RESULT_PATH, _MAX_RESULT_BYTES)
    if _sha(config_raw) != PREF2_CONFIG_SHA256 or _sha(result_raw) != PREF2_RESULT_SHA256:
        _fail("PREF2 tracked bundle differs")
    try:
        regenerated = pref2.build_pref2_result(config_raw, root)
    except Exception as exc:
        raise TDG9AR1AuthorityError("PREF2 live terminal authentication failed") from exc
    if pref2.canonical_result(regenerated) != result_raw:
        _fail("PREF2 live terminal differs from its compact result")


def _payload(config: Mapping[str, Any], inventory: list[dict[str, str]]) -> dict[str, object]:
    return {
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": TARGET_PROTOCOL,
        "predecessor": dict(config["predecessor"]),
        "exact_arithmetic": dict(config["exact_arithmetic"]),
        "environment": dict(config["environment"]),
        "selection": dict(config["selection"]),
        "committed_image": dict(config["committed_image"]),
        "replays": [dict(item) for item in config["replays"]],
        "channel_order": list(CHANNEL_ORDER),
        "output_namespace": OUTPUT_NAMESPACE,
        "output_namespace_absent_at_prelaunch": True,
        "PREF2_live_terminal_authenticated_at_prelaunch": True,
        "implementation_inventory": inventory,
        "scope": dict(SCOPE),
        "claims": dict(CLAIMS),
        "nonclaims": list(NONCLAIMS),
    }


def build_prelaunch(config_raw: bytes, root: Path) -> dict[str, object]:
    repository = Path(root).resolve()
    config = _parse_config(config_raw)
    _observe_environment()
    _require_prospective_image(repository)
    _validate_predecessors(repository)
    inventory = _inventory(repository, config)
    _authenticate_pref2_live(repository)
    _require_namespace_absent(repository)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "gate_status": GATE_STATUS,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(config, inventory),
    }


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    config = _parse_config(config_raw)
    try:
        result = json.loads(result_raw, object_pairs_hook=_unique)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TDG9AR1AuthorityError("AR1 compact result is invalid JSON") from exc
    expected = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "gate_status": GATE_STATUS,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(config, [dict(item) for item in config["implementation_inventory"]]),
    }
    if result != expected or canonical_pretty(result) != result_raw:
        _fail("AR1 compact result differs")
    return result


def _unique(pairs: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError("duplicate JSON key")
        answer[key] = value
    return answer


def authorize_diagnostic(
    root: Path,
    config_raw: bytes,
    result_raw: bytes,
    authority_commit: str,
) -> TDG9AR1Authority:
    repository = Path(root).resolve()
    if not isinstance(authority_commit, str) or not _COMMIT.fullmatch(authority_commit):
        _fail("AR1 authority commit is malformed")
    _observe_environment()
    _require_runtime_image(repository, authority_commit)
    config = _parse_config(config_raw)
    result = validate_compact(config_raw, result_raw)
    _validate_predecessors(repository)
    inventory = _inventory(repository, config)
    if _git(repository, "show", f"{authority_commit}:{CONFIG_PATH}") != config_raw:
        _fail("AR1 committed config differs")
    if _git(repository, "show", f"{authority_commit}:{RESULT_PATH}") != result_raw:
        _fail("AR1 committed result differs")
    for item in inventory:
        _require_commit_blob(repository, authority_commit, item["path"], item["sha256"])
    payload = result["artifact_payload"]
    if not isinstance(payload, Mapping) or payload.get("claims") != CLAIMS:
        _fail("AR1 typed claim surface differs")
    return TDG9AR1Authority(
        authority_commit=authority_commit,
        pref2_commit=PREF2_COMMIT,
        evaluator_commit=EVALUATOR_COMMIT,
        store_manifest_sha256=PREF2_STORE_MANIFEST_SHA256,
        output_namespace=OUTPUT_NAMESPACE,
        retries=tuple(item["retry"] for item in REPLAYS),
        retry_width_hex=tuple(item["attempted_width_hex"] for item in REPLAYS),
        max_depth=MAX_DEPTH,
        max_nodes_per_channel_per_evaluator=MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
        max_evaluator_calls=MAX_EVALUATOR_CALLS,
        channel_order=CHANNEL_ORDER,
        environment=tuple(ENVIRONMENT.items()),
        selection=tuple(SELECTION.items()),
        authority_delta_paths=AUTHORITY_DELTA_PATHS,
        config_sha256=_sha(config_raw),
        result_sha256=_sha(result_raw),
        implementation_inventory=tuple(
            (item["role"], item["path"], item["sha256"]) for item in inventory
        ),
    )


__all__ = [
    "ARTIFACT_ID", "CLASSIFICATION", "CONFIG_PATH", "RESULT_PATH",
    "OUTPUT_NAMESPACE", "STAGING_PREFIX", "REPLAYS", "CHANNEL_ORDER", "MAX_DEPTH",
    "MAX_NODES_PER_CHANNEL_PER_EVALUATOR", "MAX_EVALUATOR_CALLS",
    "TDG9AR1Authority", "TDG9AR1AuthorityError", "authorize_diagnostic",
    "build_prelaunch", "canonical", "canonical_pretty", "validate_compact",
]
