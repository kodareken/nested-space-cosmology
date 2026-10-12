"""Outcome-neutral, read-only binder for the PRO19 generation-nine terminal.

This module authenticates only durable evidence.  It does not import or call
the event runner, recovery coordinator, progression attempt, proposal code, or
candidate runtimes.  The mutable campaign and the three SID3 authority inputs
are read through no-follow descriptors; Git is queried only for immutable
commit/tree/blob facts.  No writer or recovery API is exposed here.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import tomllib
from typing import Any, Iterable, Mapping
import zipfile

import numpy as np


ARTIFACT_ID = "FGC-1-PRO19-PREF28"
SCHEMA_VERSION = 1
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "outcome_neutral_generation9_invalid_terminal_binder"
CONFIG_PATH = "configs/fgc/fgc-1-pro19-pref28.toml"
RESULT_PATH = "results/fgc-1-pro19-pref28.json"
STORE_PATH = "runs/fgc-2-sf1/proto17/calibration"

AUTHORITY_COMMIT = "b8fea44384be62918670ccbab5bd1a2d92057224"
IMPLEMENTATION_COMMIT = "4d7133c0737524712220f77ee5ae4d9c3550c75a"
CAMPAIGN_ID = "FGC-2-SF1-PROTO17-GR0-A3-CAL-290a65bcd6a2ba820683"
ORIGINAL_AUTHORIZATION_COMMIT = "c11ba422ce49ddd6f6175de1a9d2da675b97f2de"
PLAN_SHA256 = "f3b365978aeb2ae6063742807a3161bf5d6aa0c8392cd283fdec2a79d3e46060"
GEN8_CHECKPOINT_SHA256 = "6d2b6e37e7c7cfd4ce64504353d886099003f959ea2aaee1a14fc778ce14fa46"
GEN8_CHECKPOINT_RAW_SHA256 = "dca398efa21d699ed09060c2afd1a13b7b3672ca62573c24bd846a0ff6a3ef7d"
SEQ8_JOURNAL_SHA256 = "949b29e46081bf5a42c34461b3894419db86aeca6c733b7a2636c91c6cb2ee4a"
GEN9_CHECKPOINT_SHA256 = "770728faa3b0a7233e755ceadb2cc10a7109d231989e3824a40004da7ded5ea7"
GEN9_CHECKPOINT_RAW_SHA256 = "3fd5eb4e79ffaaa0e62bbd28597fc551fb53a97dc35170a4e4bfe3854ed43421"
SEQ9_JOURNAL_SHA256 = "9e8b8a58aea5a1e615e9f316e5c8d1872b6d4e813b09765f8da6b339607688bb"
SEQ9_JOURNAL_RAW_SHA256 = "d21f0eaf765c307bd0377503c510cae7e593a5912539e5318c17560d59f11e18"
TERMINAL_LOCK_SHA256 = "1726b50aa3c7921819bcd077939ddb4f82c2dc994247a6bcc68600cd61920d57"
TERMINAL_LOCK_RAW_SHA256 = "2c7434d40cecaa67317fad59a2fedf323044509622285a56618e602b2a35c19c"
MEMBER_MAP_SHA256 = "4dac7e098558313bc812832ec60abd8ec9b055fed6c70cc0e34aaaeb4c5541f3"

MEMBER_KEYS = (
    "RK4-2049",
    "RK4-4097",
    "RK4-8193",
    "SSPRK3-4097",
    "SSPRK3-8193",
    "SSPRK3-16385",
)
PHYSICAL_MEMBER = "RK4-2049"
PHYSICAL_DESCRIPTOR = "6b0dd985ab95ac52950348df5079a6b8fc4f5e7a78cf95031237cb0155c49c56"
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"

EXPECTED_PREDECESSOR = {
    "artifact_id": "FGC-1-PRO19-SID3-AUTH1",
    "authority_commit": AUTHORITY_COMMIT,
    "implementation_commit": IMPLEMENTATION_COMMIT,
    "config_path": "configs/fgc/fgc-1-pro19-sid3-auth1.toml",
    "config_sha256": "54c6f972d2b80d15e5b5a2307c1c69135e0c0be23b63ca1d412d8b20ef1898e5",
    "execution_closure_path": "configs/fgc/fgc-1-pro19-sid3-execution-closure.json",
    "execution_closure_sha256": "494842c38faa4337544259b0f204dc285a1ff3b7be50061187b959cb5ca3b136",
    "result_path": "results/fgc-1-pro19-sid3-auth1.json",
    "result_sha256": "0d0371f74d4bcb503482047b94cdb4a2b269eeaf275ed80bb0a6317420478500",
    "hlt16_config_path": "configs/fgc/fgc-1-hlt16-mon16.toml",
    "hlt16_config_sha256": "8ee09fbdd7a750849858d7fce252a01195f0f34862ebc892f2535dd5f8c0e526",
    "hlt16_result_path": "results/fgc-1-hlt16-mon16.json",
    "hlt16_result_sha256": "f6b8a64afc569af31ca4df943582db520b8eb714e38091fb6087c8ed50bacc52",
    "launch_manifest_path": "configs/fgc/fgc-1-pro19-launch-authority.toml",
    "launch_manifest_sha256": "0ba671505e5d5a791d19f0d9c48880d95a79014874a73e03e9a27a4bd98f9786",
}
EXPECTED_CAMPAIGN = {
    "original_authorization_commit": ORIGINAL_AUTHORIZATION_COMMIT,
    "original_plan_sha256": PLAN_SHA256,
    "campaign_id": CAMPAIGN_ID,
    "branch": "GR-0",
    "amplitude": "3",
    "event": 23,
    "target_rational": "3/2",
}
EXPECTED_GENERATION8 = {
    "generation": 8,
    "checkpoint_sha256": GEN8_CHECKPOINT_SHA256,
    "checkpoint_raw_sha256": GEN8_CHECKPOINT_RAW_SHA256,
    "journal_sequence": 8,
    "journal_sha256": SEQ8_JOURNAL_SHA256,
    "disposition": "nonterminal",
}
EXPECTED_GENERATION9 = {
    "generation": 9,
    "checkpoint_sha256": GEN9_CHECKPOINT_SHA256,
    "checkpoint_raw_sha256": GEN9_CHECKPOINT_RAW_SHA256,
    "parent_checkpoint_sha256": GEN8_CHECKPOINT_SHA256,
    "journal_sequence": 9,
    "journal_sha256": SEQ9_JOURNAL_SHA256,
    "journal_raw_sha256": SEQ9_JOURNAL_RAW_SHA256,
    "journal_parent_sha256": SEQ8_JOURNAL_SHA256,
    "disposition": "invalid_terminal",
    "event": 23,
    "target_rational": "3/2",
}
EXPECTED_TERMINAL = {
    "member_key": PHYSICAL_MEMBER,
    "classification": "invalid",
    "owner": "invalid",
    "reason": "invalid",
    "exception_type": "HLT16AttemptError",
    "evidence_kind": "invalid",
    "rejection_record_sha256": "",
    "lock_path": "locks/terminal.lock",
    "lock_sha256": TERMINAL_LOCK_SHA256,
    "lock_raw_sha256": TERMINAL_LOCK_RAW_SHA256,
}
ADDITIVE_SUFFIX_PATHS = (
    f"checkpoints/{9:020d}-{GEN9_CHECKPOINT_SHA256}.json",
    f"journal/{9:020d}-{SEQ9_JOURNAL_SHA256}.journal",
    "locks/.active-write.lock.hlt16-quarantine-33ee84d9183e76f764f67680f3455416c496f431d59d33f568fda926f96ce043",
    "locks/terminal.lock",
)
EXPECTED_STORE = {
    "path": STORE_PATH,
    "regular_leaf_count": 56,
    "lexical_inventory_sha256": "c574d7f92403f8d1dfcaae880edda34fa8de9f3809a508c13db0e58f89f9a216",
    "generation8_baseline_leaf_count": 52,
    "generation8_baseline_inventory_sha256": "8196b088a19f70ae1016a3d9dddb508add1ce1ce5cd227d105c954bd0c151c27",
    "additive_suffix_paths": list(ADDITIVE_SUFFIX_PATHS),
    "active_writer_absent": True,
    "terminal_lock_present": True,
    "staging_absent": True,
    "orphan_absent": True,
    "fork_absent": True,
    "symlink_and_special_file_absent": True,
}
EXPECTED_UNCHANGED_STATE = {
    "member_map_sha256": MEMBER_MAP_SHA256,
    "member_descriptor_count": 6,
    "new_state_or_payload_leaf_count": 0,
    "affected_member_key": PHYSICAL_MEMBER,
    "accepted_time_hex": ACCEPTED_TIME_HEX,
    "descriptor_sha256": PHYSICAL_DESCRIPTOR,
    "u_bytes_sha256": "856fc327d3163d038f9fdbfc660968509dd7360faa0bb9d158721eb3049c2fe4",
    "p_bytes_sha256": "32939115b86938cf8a813bbf9a50a8c90bbaceaee161fa65695deb435d1243b7",
    "q_bytes_sha256": "48ec07bcace3935173c06e7b41915cacb17f14e1c9d3500522d4eb59a205c4a3",
}
EXPECTED_SCOPE = {
    "live_store_read_only": True,
    "PDE_proposal_execution": False,
    "accepted_state_mutation": False,
    "event_commit": False,
    "trajectory_resume": False,
    "candidate_branch_opened": False,
    "root_cause_diagnosis": False,
}
EXPECTED_CLAIMS = {
    "generation9_terminal_boundary_bound": True,
    "invalid_runtime_terminal_bound": True,
    "accepted_member_map_unchanged": True,
    "terminal_lock_valid": True,
    "root_cause_localized": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "candidate_runtime_configuration_state_or_outcome_opened": False,
    "physical_result_earned": False,
}
EXPECTED_NONCLAIMS = [
    "PREF28 binds one fail-closed invalid terminal record; it does not identify the unpersisted HLT16AttemptError message or diagnose its root cause.",
    "PREF28 proves that no accepted member state changed and that event 23 remains incomplete; it is not a GR-0 calibration or physical result.",
    "PREF28 authorizes no resume, rollback, repair, new campaign, SGB-L, FGC-QR, DEF1, retained-EFT, transition, or physical claim.",
]

_EXPECTED_DIRECTORIES = ("checkpoints", "journal", "locks", "payloads", "receipts", "states")
_CHECKPOINT_NAME = re.compile(r"(?P<number>[0-9]{20})-(?P<digest>[0-9a-f]{64})\.json\Z")
_JOURNAL_NAME = re.compile(r"(?P<number>[0-9]{20})-(?P<digest>[0-9a-f]{64})\.journal\Z")
_STATE_NAME = re.compile(r"(?P<digest>[0-9a-f]{64})\.json\Z")
_PAYLOAD_NAME = re.compile(r"(?P<digest>[0-9a-f]{64})\.npz\Z")
_QUARANTINE_NAME = re.compile(r"\.active-write\.lock\.hlt16-quarantine-(?P<digest>[0-9a-f]{64})\Z")
_MAX_LEAF_BYTES = 96 * 1024 * 1024
_MAX_TREE_BYTES = 384 * 1024 * 1024


class PREF28BinderError(ValueError):
    """The PREF28 config, authority, store, or compact record differs."""

    def __init__(self, stop_id: str, detail: str) -> None:
        super().__init__(f"{stop_id}: {detail}")
        self.stop_id = stop_id
        self.detail = detail


def _stop(stop_id: str, detail: str) -> None:
    raise PREF28BinderError(stop_id, detail)


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as error:
        _stop("PREF28_CANONICAL_DRIFT", "value is not canonical JSON")
        raise AssertionError from error


def canonical_result(value: object) -> bytes:
    """Encode one deterministic tracked PREF28 result."""
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
    except (TypeError, ValueError) as error:
        _stop("PREF28_COMPACT_DRIFT", "result is not canonical JSON")
        raise AssertionError from error


def _digest(value: object) -> str:
    return sha256(_canonical(value)).hexdigest()


def _duplicates(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError(key)
        answer[key] = value
    return answer


def _json(raw: bytes, label: str, *, pretty: bool = False) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        _stop("PREF28_CANONICAL_DRIFT", f"{label} is malformed")
        raise AssertionError from error
    expected = canonical_result(value) if pretty else _canonical(value)
    if not isinstance(value, dict) or expected != raw:
        _stop("PREF28_CANONICAL_DRIFT", f"{label} is noncanonical")
    return value


def _object_address(value: Mapping[str, Any], field: str, address: str, label: str) -> None:
    if value.get(field) != address:
        _stop("PREF28_HASH_DRIFT", f"{label} embedded address differs")
    body = dict(value)
    body.pop(field, None)
    if _digest(body) != address:
        _stop("PREF28_HASH_DRIFT", f"{label} content address differs")


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        _stop("PREF28_CONFIG_DRIFT", "PREF28 config is malformed")
        raise AssertionError from error
    if not isinstance(value, dict) or set(value) != {
        "schema_version",
        "artifact_id",
        "project_version",
        "target_protocol",
        "classification",
        "nonclaims",
        "predecessor",
        "campaign",
        "generation8",
        "generation9",
        "terminal",
        "store",
        "unchanged_state",
        "scope",
        "claims",
    }:
        _stop("PREF28_CONFIG_DRIFT", "PREF28 config fields differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != EXPECTED_NONCLAIMS
        or value["predecessor"] != EXPECTED_PREDECESSOR
        or value["campaign"] != EXPECTED_CAMPAIGN
        or value["generation8"] != EXPECTED_GENERATION8
        or value["generation9"] != EXPECTED_GENERATION9
        or value["terminal"] != EXPECTED_TERMINAL
        or value["store"] != EXPECTED_STORE
        or value["unchanged_state"] != EXPECTED_UNCHANGED_STATE
        or value["scope"] != EXPECTED_SCOPE
        or value["claims"] != EXPECTED_CLAIMS
    ):
        _stop("PREF28_CONFIG_DRIFT", "PREF28 exact contract differs")
    return value


@dataclass(frozen=True, slots=True)
class _Leaf:
    relative: str
    byte_size: int
    raw_sha256: str
    raw: bytes


def _directory_flags() -> int:
    return os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)


def _open_directory_chain(root: Path, parts: tuple[str, ...], label: str) -> int:
    try:
        root_info = root.lstat()
    except OSError as error:
        _stop("PREF28_PATH_UNSAFE", f"{label} root is absent")
        raise AssertionError from error
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        _stop("PREF28_PATH_UNSAFE", f"{label} root is unsafe")
    try:
        current = os.open(root, _directory_flags())
    except OSError as error:
        _stop("PREF28_PATH_UNSAFE", f"{label} root cannot be opened")
        raise AssertionError from error
    try:
        for component in parts:
            if component in {"", ".", ".."}:
                _stop("PREF28_PATH_UNSAFE", f"{label} component is unsafe")
            before = os.stat(component, dir_fd=current, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                _stop("PREF28_PATH_UNSAFE", f"{label} ancestor is unsafe: {component}")
            child = os.open(component, _directory_flags(), dir_fd=current)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                os.close(child)
                _stop("PREF28_PATH_UNSAFE", f"{label} ancestor raced: {component}")
            os.close(current)
            current = child
        return current
    except OSError as error:
        os.close(current)
        _stop("PREF28_PATH_UNSAFE", f"{label} ancestor cannot be opened safely")
        raise AssertionError from error
    except Exception:
        os.close(current)
        raise


def _read_leaf_fd(parent: int, name: str, before: os.stat_result, label: str) -> bytes:
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > _MAX_LEAF_BYTES
    ):
        _stop("PREF28_PATH_UNSAFE", f"{label} is not a single-linked bounded regular file")
    descriptor = -1
    try:
        descriptor = os.open(name, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0), dir_fd=parent)
        active = os.fstat(descriptor)
        identity = (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
            before.st_nlink,
        )
        if identity != (
            active.st_dev,
            active.st_ino,
            active.st_size,
            active.st_mtime_ns,
            active.st_ctime_ns,
            active.st_nlink,
        ):
            _stop("PREF28_PATH_UNSAFE", f"{label} raced before read")
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _stop("PREF28_PATH_UNSAFE", f"{label} ended early")
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _stop("PREF28_PATH_UNSAFE", f"{label} grew during read")
        raw = b"".join(chunks)
        after = os.fstat(descriptor)
        if identity != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
            after.st_nlink,
        ):
            _stop("PREF28_PATH_UNSAFE", f"{label} changed during read")
        return raw
    except OSError as error:
        _stop("PREF28_PATH_UNSAFE", f"{label} cannot be read safely")
        raise AssertionError from error
    finally:
        if descriptor != -1:
            os.close(descriptor)


def _read_repository_leaf(repository: Path, relative: str, label: str) -> bytes:
    path = Path(relative)
    if path.is_absolute() or not path.parts or any(part in {"", ".", ".."} for part in path.parts):
        _stop("PREF28_PATH_UNSAFE", f"{label} relative path is unsafe")
    parent = _open_directory_chain(Path(repository), tuple(path.parts[:-1]), label)
    try:
        before = os.stat(path.parts[-1], dir_fd=parent, follow_symlinks=False)
        return _read_leaf_fd(parent, path.parts[-1], before, label)
    except FileNotFoundError as error:
        _stop("PREF28_PATH_UNSAFE", f"{label} is absent")
        raise AssertionError from error
    finally:
        os.close(parent)


def _git(repository: Path, *args: str) -> bytes:
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"})
    try:
        return subprocess.run(
            [
                "git",
                "--no-replace-objects",
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                *args,
            ],
            cwd=repository,
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        _stop("PREF28_AUTHORITY_GIT_INVALID", "Git authority lookup failed")
        raise AssertionError from error


def _bind_authority(repository: Path, predecessor: Mapping[str, Any]) -> dict[str, Any]:
    """Authenticate authority C, its sole parent, and its exact SID3 blobs."""
    repository = Path(repository)
    root_fd = _open_directory_chain(repository, (), "repository")
    os.close(root_fd)
    try:
        top = Path(_git(repository, "rev-parse", "--show-toplevel").decode("utf-8").strip())
        if not os.path.samefile(repository, top):
            _stop("PREF28_AUTHORITY_GIT_INVALID", "repository is not the Git worktree root")
    except (UnicodeDecodeError, OSError) as error:
        _stop("PREF28_AUTHORITY_GIT_INVALID", "Git worktree root differs")
        raise AssertionError from error

    resolved = _git(repository, "rev-parse", "--verify", f"{AUTHORITY_COMMIT}^{{commit}}")
    if resolved != f"{AUTHORITY_COMMIT}\n".encode("ascii"):
        _stop("PREF28_AUTHORITY_GIT_INVALID", "authority commit differs")
    parents = _git(repository, "rev-list", "--parents", "-n", "1", AUTHORITY_COMMIT)
    if parents != f"{AUTHORITY_COMMIT} {IMPLEMENTATION_COMMIT}\n".encode("ascii"):
        _stop("PREF28_AUTHORITY_GIT_INVALID", "authority parent relation differs")
    _git(repository, "merge-base", "--is-ancestor", AUTHORITY_COMMIT, "HEAD")

    raw_by_path: dict[str, bytes] = {}
    for path_key, hash_key in (
        ("config_path", "config_sha256"),
        ("execution_closure_path", "execution_closure_sha256"),
        ("result_path", "result_sha256"),
        ("hlt16_config_path", "hlt16_config_sha256"),
        ("hlt16_result_path", "hlt16_result_sha256"),
        ("launch_manifest_path", "launch_manifest_sha256"),
    ):
        relative = str(predecessor[path_key])
        expected_hash = str(predecessor[hash_key])
        tree = _git(repository, "ls-tree", "-z", AUTHORITY_COMMIT, "--", relative)
        entries = [entry for entry in tree.split(b"\0") if entry]
        if len(entries) != 1 or not entries[0].startswith(b"100644 blob "):
            _stop("PREF28_AUTHORITY_GIT_INVALID", f"SID3 blob is not tracked: {relative}")
        try:
            _, recorded_path = entries[0].split(b"\t", 1)
        except ValueError:
            _stop("PREF28_AUTHORITY_GIT_INVALID", f"SID3 tree entry differs: {relative}")
        if recorded_path.decode("utf-8") != relative:
            _stop("PREF28_AUTHORITY_GIT_INVALID", f"SID3 tree path differs: {relative}")
        committed = _git(repository, "show", f"{AUTHORITY_COMMIT}:{relative}")
        if sha256(committed).hexdigest() != expected_hash:
            _stop("PREF28_AUTHORITY_DRIFT", f"SID3 committed bytes differ: {relative}")
        if _read_repository_leaf(repository, relative, f"SID3 {relative}") != committed:
            _stop("PREF28_AUTHORITY_DRIFT", f"SID3 worktree bytes differ: {relative}")
        raw_by_path[relative] = committed

    try:
        sid3_config = tomllib.loads(raw_by_path[str(predecessor["config_path"])].decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        _stop("PREF28_AUTHORITY_DRIFT", "SID3 config is malformed")
        raise AssertionError from error
    closure = _json(
        raw_by_path[str(predecessor["execution_closure_path"])],
        "SID3 execution closure",
        pretty=True,
    )
    result = _json(raw_by_path[str(predecessor["result_path"])], "SID3 result", pretty=True)
    execution = sid3_config.get("execution_closure")
    image = result.get("artifact_payload", {}).get("execution_image") if isinstance(result.get("artifact_payload"), Mapping) else None
    if (
        sid3_config.get("artifact_id") != predecessor["artifact_id"]
        or not isinstance(execution, Mapping)
        or execution.get("path") != predecessor["execution_closure_path"]
        or execution.get("raw_sha256") != predecessor["execution_closure_sha256"]
        or execution.get("implementation_commit") != IMPLEMENTATION_COMMIT
        or closure.get("schema") != "FGC-1-PRO19-execution-source-closure-v1"
        or closure.get("implementation_commit") != IMPLEMENTATION_COMMIT
        or result.get("artifact_id") != predecessor["artifact_id"]
        or result.get("gate_status") != "pass"
        or result.get("source_config_sha256") != predecessor["config_sha256"]
        or result.get("execution_closure_raw_sha256") != predecessor["execution_closure_sha256"]
        or not isinstance(image, Mapping)
        or image.get("implementation_commit") != IMPLEMENTATION_COMMIT
    ):
        _stop("PREF28_AUTHORITY_DRIFT", "SID3 authority cross-binding differs")
    return dict(predecessor)


def _scan_tree_fd(root_fd: int) -> tuple[dict[str, _Leaf], tuple[str, ...]]:
    leaves: dict[str, _Leaf] = {}
    directories: list[str] = []
    total = 0

    def descend(directory_fd: int, prefix: str) -> None:
        nonlocal total
        before_directory = os.fstat(directory_fd)
        try:
            with os.scandir(directory_fd) as scan:
                entries = sorted(scan, key=lambda entry: entry.name)
        except OSError as error:
            _stop("PREF28_PATH_UNSAFE", "store cannot be scanned safely")
            raise AssertionError from error
        names = tuple(entry.name for entry in entries)
        for entry in entries:
            name = entry.name
            relative = f"{prefix}/{name}" if prefix else name
            try:
                name.encode("utf-8")
                before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            except (UnicodeEncodeError, OSError) as error:
                _stop("PREF28_PATH_UNSAFE", f"store node cannot be statted: {relative}")
                raise AssertionError from error
            if stat.S_ISLNK(before.st_mode):
                _stop("PREF28_PATH_UNSAFE", f"symlink appears in store: {relative}")
            if stat.S_ISDIR(before.st_mode):
                try:
                    child = os.open(name, _directory_flags(), dir_fd=directory_fd)
                except OSError as error:
                    _stop("PREF28_PATH_UNSAFE", f"directory cannot be opened: {relative}")
                    raise AssertionError from error
                try:
                    after = os.fstat(child)
                    if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                        _stop("PREF28_PATH_UNSAFE", f"directory raced: {relative}")
                    directories.append(relative)
                    descend(child, relative)
                finally:
                    os.close(child)
                continue
            if not stat.S_ISREG(before.st_mode):
                _stop("PREF28_PATH_UNSAFE", f"special file appears in store: {relative}")
            raw = _read_leaf_fd(directory_fd, name, before, relative)
            total += len(raw)
            if total > _MAX_TREE_BYTES:
                _stop("PREF28_PATH_UNSAFE", "store exceeds the read-only byte budget")
            leaves[relative] = _Leaf(relative, len(raw), sha256(raw).hexdigest(), raw)

        try:
            with os.scandir(directory_fd) as scan:
                after_names = tuple(sorted(entry.name for entry in scan))
        except OSError as error:
            _stop("PREF28_PATH_UNSAFE", "store cannot be rescanned safely")
            raise AssertionError from error
        after_directory = os.fstat(directory_fd)
        if names != after_names or (
            before_directory.st_dev,
            before_directory.st_ino,
            before_directory.st_mtime_ns,
            before_directory.st_ctime_ns,
        ) != (
            after_directory.st_dev,
            after_directory.st_ino,
            after_directory.st_mtime_ns,
            after_directory.st_ctime_ns,
        ):
            _stop("PREF28_PATH_UNSAFE", f"store directory raced: {prefix or '.'}")

    descend(root_fd, "")
    return leaves, tuple(sorted(directories))


def _stage_siblings(parent: int, final: str) -> tuple[str, ...]:
    try:
        with os.scandir(parent) as scan:
            return tuple(
                sorted(
                    entry.name
                    for entry in scan
                    if entry.name.startswith(f".{final}.hlt15-stage-")
                    or entry.name.startswith(f".{final}.hlt16-stage-")
                )
            )
    except OSError as error:
        _stop("PREF28_PATH_UNSAFE", "store parent cannot be scanned")
        raise AssertionError from error


def _scan_store(repository: Path) -> tuple[dict[str, _Leaf], tuple[str, ...]]:
    parts = Path(STORE_PATH).parts
    parent = _open_directory_chain(Path(repository), parts[:-1], "repository/store")
    store_fd = -1
    try:
        final = parts[-1]
        if _stage_siblings(parent, final):
            _stop("PREF28_TREE_DRIFT", "campaign staging sibling is present")
        before = os.stat(final, dir_fd=parent, follow_symlinks=False)
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
            _stop("PREF28_PATH_UNSAFE", "campaign store root is unsafe")
        store_fd = os.open(final, _directory_flags(), dir_fd=parent)
        after = os.fstat(store_fd)
        if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
            _stop("PREF28_PATH_UNSAFE", "campaign store root raced")
        first = _scan_tree_fd(store_fd)
        second = _scan_tree_fd(store_fd)
        if first != second or _stage_siblings(parent, final):
            _stop("PREF28_PATH_UNSAFE", "campaign store changed during authentication")
        return first
    except FileNotFoundError as error:
        _stop("PREF28_PATH_UNSAFE", "campaign store is absent")
        raise AssertionError from error
    except OSError as error:
        _stop("PREF28_PATH_UNSAFE", "campaign store cannot be opened safely")
        raise AssertionError from error
    finally:
        if store_fd != -1:
            os.close(store_fd)
        os.close(parent)


def _inventory_sha256(leaves: Mapping[str, _Leaf]) -> str:
    digestor = sha256()
    for relative in sorted(leaves):
        leaf = leaves[relative]
        digestor.update(relative.encode("utf-8"))
        digestor.update(b"\0")
        digestor.update(str(leaf.byte_size).encode("ascii"))
        digestor.update(b"\0")
        digestor.update(leaf.raw_sha256.encode("ascii"))
        digestor.update(b"\n")
    return digestor.hexdigest()


def _validate_tree_grammar(leaves: Mapping[str, _Leaf], directories: tuple[str, ...]) -> None:
    if directories != _EXPECTED_DIRECTORIES:
        _stop("PREF28_TREE_DRIFT", "store directory grammar differs")
    for relative, leaf in leaves.items():
        parts = Path(relative).parts
        if len(parts) != 2 or parts[0] not in _EXPECTED_DIRECTORIES:
            _stop("PREF28_TREE_DRIFT", f"nested or foreign leaf appears: {relative}")
        directory, name = parts
        if ".hlt16-stage-" in name or ".hlt15-stage-" in name:
            _stop("PREF28_TREE_DRIFT", "publication staging residue is present")
        valid = False
        if directory == "checkpoints":
            valid = _CHECKPOINT_NAME.fullmatch(name) is not None
        elif directory == "journal":
            valid = _JOURNAL_NAME.fullmatch(name) is not None
        elif directory == "states":
            valid = _STATE_NAME.fullmatch(name) is not None
        elif directory == "payloads":
            valid = _PAYLOAD_NAME.fullmatch(name) is not None
        elif directory == "receipts":
            valid = name == "generation-zero.json"
        elif directory == "locks":
            if name == "active-write.lock":
                _stop("PREF28_TREE_DRIFT", "active writer lease is present")
            valid = (
                name in {"bootstrap.guard", "writer.guard", "terminal.lock"}
                or _QUARANTINE_NAME.fullmatch(name) is not None
            )
            if name in {"bootstrap.guard", "writer.guard"} and leaf.raw != b"":
                _stop("PREF28_TREE_DRIFT", f"advisory guard bytes differ: {name}")
        if not valid:
            _stop("PREF28_TREE_DRIFT", f"foreign store leaf appears: {relative}")


def _indexed_objects(
    leaves: Mapping[str, _Leaf],
    directory: str,
    pattern: re.Pattern[str],
    address_field: str,
) -> dict[int, dict[str, Any]]:
    answer: dict[int, dict[str, Any]] = {}
    for relative in sorted(path for path in leaves if path.startswith(f"{directory}/")):
        match = pattern.fullmatch(Path(relative).name)
        if match is None:
            _stop("PREF28_TREE_DRIFT", f"{directory} filename differs")
        index = int(match.group("number"))
        address = match.group("digest")
        if index in answer:
            _stop("PREF28_TREE_DRIFT", f"{directory} fork appears at {index}")
        value = _json(leaves[relative].raw, f"{directory} {index}")
        _object_address(value, address_field, address, f"{directory} {index}")
        answer[index] = value
    return answer


def _validate_retired_leases(leaves: Mapping[str, _Leaf]) -> None:
    for relative, leaf in leaves.items():
        if not relative.startswith("locks/.active-write.lock.hlt16-quarantine-"):
            continue
        match = _QUARANTINE_NAME.fullmatch(Path(relative).name)
        assert match is not None
        if leaf.raw_sha256 != match.group("digest"):
            _stop("PREF28_TREE_DRIFT", "retired writer content address differs")
        lease = _json(leaf.raw, "retired writer lease")
        if set(lease) != {
            "authorization_commit",
            "checkpoint_sha256",
            "host",
            "owner_token",
            "pid",
            "plan_sha256",
            "schema",
        } or (
            lease["schema"] != "FGC-1-HLT16-active-writer-v1"
            or lease["authorization_commit"] != ORIGINAL_AUTHORIZATION_COMMIT
            or lease["plan_sha256"] != PLAN_SHA256
            or not isinstance(lease["host"], str)
            or not lease["host"]
            or isinstance(lease["pid"], bool)
            or not isinstance(lease["pid"], int)
            or lease["pid"] <= 0
            or not isinstance(lease["owner_token"], str)
            or len(lease["owner_token"]) != 64
        ):
            _stop("PREF28_TREE_DRIFT", "retired writer lease differs")


def _validate_physical_payload(leaves: Mapping[str, _Leaf]) -> dict[str, str]:
    descriptor_leaf = leaves.get(f"states/{PHYSICAL_DESCRIPTOR}.json")
    if descriptor_leaf is None:
        _stop("PREF28_STATE_DRIFT", "affected member descriptor is absent")
    descriptor = _json(descriptor_leaf.raw, "affected member descriptor")
    _object_address(descriptor, "descriptor_sha256", PHYSICAL_DESCRIPTOR, "affected member descriptor")
    manifest = descriptor.get("arrays")
    semantic = descriptor.get("semantic_sha256")
    archive_hash = descriptor.get("raw_archive_sha256")
    if not isinstance(manifest, list) or not isinstance(semantic, str) or _digest(manifest) != semantic:
        _stop("PREF28_STATE_DRIFT", "affected member manifest differs")
    payload_leaf = leaves.get(f"payloads/{semantic}.npz")
    if payload_leaf is None or payload_leaf.raw_sha256 != archive_hash:
        _stop("PREF28_STATE_DRIFT", "affected member payload differs")
    by_name = {
        item.get("name"): item
        for item in manifest
        if isinstance(item, Mapping) and isinstance(item.get("name"), str)
    }
    observed: dict[str, str] = {}
    try:
        with zipfile.ZipFile(BytesIO(payload_leaf.raw), "r") as archive:
            if len(archive.infolist()) != len({item.filename for item in archive.infolist()}):
                _stop("PREF28_STATE_DRIFT", "affected member ZIP has duplicate entries")
            for name in ("u", "p", "q"):
                info = archive.getinfo(f"{name}.npy")
                if (
                    info.is_dir()
                    or info.flag_bits & 0x1
                    or info.file_size <= 0
                    or info.file_size > 32 * 1024 * 1024
                    or Path(info.filename).name != info.filename
                ):
                    _stop("PREF28_STATE_DRIFT", f"affected member {name} archive entry is unsafe")
                encoded = archive.read(info)
                value = np.load(BytesIO(encoded), allow_pickle=False)
                if (
                    not isinstance(value, np.ndarray)
                    or value.dtype != np.dtype("<f8")
                    or not value.flags.c_contiguous
                    or not value.shape
                    or not np.isfinite(value).all()
                ):
                    _stop("PREF28_STATE_DRIFT", f"affected member {name} array differs")
                digest = sha256(value.tobytes(order="C")).hexdigest()
                record = by_name.get(name)
                if not isinstance(record, Mapping) or record.get("bytes_sha256") != digest:
                    _stop("PREF28_STATE_DRIFT", f"affected member {name} manifest hash differs")
                observed[f"{name}_bytes_sha256"] = digest
    except (KeyError, OSError, ValueError, EOFError, zipfile.BadZipFile) as error:
        _stop("PREF28_STATE_DRIFT", "affected member payload is malformed")
        raise AssertionError from error
    return observed


def _bind_terminal_edge(
    leaves: Mapping[str, _Leaf],
    checkpoints: Mapping[int, Mapping[str, Any]],
    journals: Mapping[int, Mapping[str, Any]],
) -> None:
    if set(checkpoints) != set(range(10)) or set(journals) != set(range(10)):
        _stop("PREF28_LINEAGE_DRIFT", "checkpoint/journal sequence has a gap or fork")
    before, after, record = checkpoints[8], checkpoints[9], journals[9]
    if (
        before.get("checkpoint_sha256") != GEN8_CHECKPOINT_SHA256
        or before.get("generation") != 8
        or before.get("journal_sequence") != 8
        or before.get("journal_tip_sha256") != SEQ8_JOURNAL_SHA256
        or before.get("disposition") != "nonterminal"
        or before.get("terminal") is not None
        or after.get("checkpoint_sha256") != GEN9_CHECKPOINT_SHA256
        or after.get("generation") != 9
        or after.get("parent_sha256") != GEN8_CHECKPOINT_SHA256
        or after.get("journal_sequence") != 9
        or after.get("journal_tip_sha256") != SEQ9_JOURNAL_SHA256
        or after.get("disposition") != "invalid_terminal"
    ):
        _stop("PREF28_LINEAGE_DRIFT", "generation-eight/nine checkpoint relation differs")
    for value in (before, after):
        if (
            value.get("authorization_commit") != ORIGINAL_AUTHORIZATION_COMMIT
            or value.get("plan_sha256") != PLAN_SHA256
            or value.get("protocol") != TARGET_PROTOCOL
            or value.get("campaign_id") != CAMPAIGN_ID
            or value.get("event") != 23
            or value.get("target")
            != {"binary64_hex": "0x1.8000000000000p+0", "rational": "3/2"}
        ):
            _stop("PREF28_LINEAGE_DRIFT", "campaign identity differs")
    before_members, after_members = before.get("members"), after.get("members")
    if (
        not isinstance(before_members, Mapping)
        or set(before_members) != set(MEMBER_KEYS)
        or after_members != before_members
        or _digest(before_members) != MEMBER_MAP_SHA256
    ):
        _stop("PREF28_STATE_DRIFT", "accepted member map advanced or differs")
    affected = before_members.get(PHYSICAL_MEMBER)
    if not isinstance(affected, Mapping):
        _stop("PREF28_STATE_DRIFT", "affected member differs")
    cursor = affected.get("cursor")
    ledger = affected.get("ledger")
    if (
        affected.get("descriptor_sha256") != PHYSICAL_DESCRIPTOR
        or not isinstance(cursor, Mapping)
        or cursor.get("accepted_state_sha256") != PHYSICAL_DESCRIPTOR
        or cursor.get("accepted_boundary_time", {}).get("binary64_hex") != ACCEPTED_TIME_HEX
        or not isinstance(ledger, Mapping)
        or ledger.get("last_accepted_time_hex") != ACCEPTED_TIME_HEX
    ):
        _stop("PREF28_STATE_DRIFT", "affected accepted identity differs")

    if (
        set(record) != {
            "schema",
            "sequence",
            "campaign_id",
            "generation",
            "event",
            "previous_record_sha256",
            "kind",
            "payload",
            "record_sha256",
        }
        or record.get("schema") != "FGC-1-HLT16-campaign-journal-v1"
        or record.get("sequence") != 9
        or record.get("campaign_id") != CAMPAIGN_ID
        or record.get("generation") != 9
        or record.get("event") != 23
        or record.get("previous_record_sha256") != SEQ8_JOURNAL_SHA256
        or record.get("kind") != "terminal_lock"
        or record.get("record_sha256") != SEQ9_JOURNAL_SHA256
    ):
        _stop("PREF28_LINEAGE_DRIFT", "sequence-nine terminal record differs")
    payload = record.get("payload")
    expected_payload = {
        "classification": "invalid",
        "cursor_sha256": cursor.get("cursor_chain_sha256"),
        "descriptor_sha256": PHYSICAL_DESCRIPTOR,
        "evidence": {"exception_type": "HLT16AttemptError", "kind": "invalid"},
        "member_key": PHYSICAL_MEMBER,
        "owner": "invalid",
        "reason": "invalid",
        "rejection_record_sha256": None,
    }
    if payload != expected_payload:
        _stop("PREF28_TERMINAL_DRIFT", "persisted terminal payload differs")
    expected_terminal = {
        "classification": "invalid",
        "event": 23,
        "member_key": PHYSICAL_MEMBER,
        "owner": "invalid",
        "reason": "invalid",
        "target": {"binary64_hex": "0x1.8000000000000p+0", "rational": "3/2"},
    }
    if after.get("terminal") != expected_terminal:
        _stop("PREF28_TERMINAL_DRIFT", "generation-nine terminal body differs")

    lock_leaf = leaves.get("locks/terminal.lock")
    if lock_leaf is None or lock_leaf.raw_sha256 != TERMINAL_LOCK_RAW_SHA256:
        _stop("PREF28_TERMINAL_DRIFT", "terminal lock raw hash differs")
    lock = _json(lock_leaf.raw, "terminal lock")
    if set(lock) != {"schema", "checkpoint_sha256", "journal_tip_sha256", "lock_sha256"}:
        _stop("PREF28_TERMINAL_DRIFT", "terminal lock fields differ")
    _object_address(lock, "lock_sha256", TERMINAL_LOCK_SHA256, "terminal lock")
    if (
        lock["schema"] != "FGC-1-HLT16-terminal-lock-v1"
        or lock["checkpoint_sha256"] != GEN9_CHECKPOINT_SHA256
        or lock["journal_tip_sha256"] != SEQ9_JOURNAL_SHA256
    ):
        _stop("PREF28_TERMINAL_DRIFT", "terminal lock relation differs")


def expected_boundary(config_raw: bytes) -> dict[str, Any]:
    """Construct the exact expected terminal anchor without live I/O."""
    config = _config(config_raw)
    return {
        "generation8": dict(config["generation8"]),
        "generation9": dict(config["generation9"]),
        "terminal": dict(config["terminal"]),
        "store": dict(config["store"]),
        "unchanged_state": dict(config["unchanged_state"]),
    }


def bind_terminal_store(config_raw: bytes, repository: Path) -> dict[str, Any]:
    """Authenticate authority C and the exact generation-eight/nine terminal edge."""
    config = _config(config_raw)
    _bind_authority(Path(repository), config["predecessor"])
    leaves, directories = _scan_store(Path(repository))
    _validate_tree_grammar(leaves, directories)
    _validate_retired_leases(leaves)

    inventory = _inventory_sha256(leaves)
    if len(leaves) != EXPECTED_STORE["regular_leaf_count"] or inventory != EXPECTED_STORE["lexical_inventory_sha256"]:
        _stop("PREF28_TREE_DRIFT", f"terminal inventory differs: leaves={len(leaves)} sha256={inventory}")
    if any(path not in leaves for path in ADDITIVE_SUFFIX_PATHS):
        _stop("PREF28_TREE_DRIFT", "additive terminal suffix is incomplete")
    baseline = {path: leaf for path, leaf in leaves.items() if path not in ADDITIVE_SUFFIX_PATHS}
    if (
        len(baseline) != EXPECTED_STORE["generation8_baseline_leaf_count"]
        or _inventory_sha256(baseline) != EXPECTED_STORE["generation8_baseline_inventory_sha256"]
        or any(path.startswith(("states/", "payloads/")) for path in ADDITIVE_SUFFIX_PATHS)
    ):
        _stop("PREF28_TREE_DRIFT", "generation-eight additive suffix proof differs")

    checkpoint8_path = f"checkpoints/{8:020d}-{GEN8_CHECKPOINT_SHA256}.json"
    checkpoint9_path = ADDITIVE_SUFFIX_PATHS[0]
    journal9_path = ADDITIVE_SUFFIX_PATHS[1]
    if (
        leaves[checkpoint8_path].raw_sha256 != GEN8_CHECKPOINT_RAW_SHA256
        or leaves[checkpoint9_path].raw_sha256 != GEN9_CHECKPOINT_RAW_SHA256
        or leaves[journal9_path].raw_sha256 != SEQ9_JOURNAL_RAW_SHA256
    ):
        _stop("PREF28_HASH_DRIFT", "generation-eight/nine raw evidence differs")

    checkpoints = _indexed_objects(leaves, "checkpoints", _CHECKPOINT_NAME, "checkpoint_sha256")
    journals = _indexed_objects(leaves, "journal", _JOURNAL_NAME, "record_sha256")
    _bind_terminal_edge(leaves, checkpoints, journals)
    physical_hashes = _validate_physical_payload(leaves)
    if physical_hashes != {
        key: EXPECTED_UNCHANGED_STATE[key]
        for key in ("u_bytes_sha256", "p_bytes_sha256", "q_bytes_sha256")
    }:
        _stop("PREF28_STATE_DRIFT", "preserved physical array hashes differ")
    return expected_boundary(config_raw)


def _result_from_config(config_raw: bytes, config: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "artifact_payload": {
            "predecessor": dict(config["predecessor"]),
            "campaign": dict(config["campaign"]),
            "live_anchor": expected_boundary(config_raw),
            "scope": dict(config["scope"]),
            "claims": dict(config["claims"]),
            "nonclaims": list(config["nonclaims"]),
        },
    }


def build_pref28_result(config_raw: bytes, repository: Path) -> dict[str, Any]:
    """Build the compact result only after authenticating the read-only boundary."""
    config = _config(config_raw)
    observed = bind_terminal_store(config_raw, Path(repository))
    if observed != expected_boundary(config_raw):
        _stop("PREF28_LINEAGE_DRIFT", "live store/config cross-binding differs")
    return _result_from_config(config_raw, config)


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    """Validate tracked bytes without Git access or reopening the ignored store."""
    config = _config(config_raw)
    try:
        value = json.loads(result_raw.decode("ascii"), object_pairs_hook=_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        _stop("PREF28_COMPACT_DRIFT", "PREF28 result is malformed")
        raise AssertionError from error
    if not isinstance(value, dict) or canonical_result(value) != result_raw:
        _stop("PREF28_COMPACT_DRIFT", "PREF28 result is noncanonical")
    if value != _result_from_config(config_raw, config):
        _stop("PREF28_COMPACT_DRIFT", "PREF28 compact result differs")
    return value


__all__ = [
    "ARTIFACT_ID",
    "CONFIG_PATH",
    "RESULT_PATH",
    "STORE_PATH",
    "PREF28BinderError",
    "bind_terminal_store",
    "build_pref28_result",
    "canonical_result",
    "expected_boundary",
    "validate_compact_result",
]
