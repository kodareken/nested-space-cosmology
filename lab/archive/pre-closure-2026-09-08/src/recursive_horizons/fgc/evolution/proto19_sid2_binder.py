"""Independent, read-only binder for the completed PRO19 SID2 recovery edge.

The binder deliberately does not import the event runner, SID1's live-anchor
reader, the HLT16 recovery coordinator, or any proposal/candidate runtime.
It traverses the persisted store with ``os.scandir`` and no-follow file
descriptors, parses the durable bytes independently, and returns only a
deterministic metadata certificate.  There is no writer or recovery API here.
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
import struct
import tomllib
from typing import Any, Iterable, Mapping
import zipfile

import numpy as np


ARTIFACT_ID = "FGC-1-PRO19-SID2-PREF1"
SCHEMA_VERSION = 1
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "post_recovery_generation8_metadata_only_binder"
CONFIG_PATH = "configs/fgc/fgc-1-pro19-sid2-pref1.toml"
RESULT_PATH = "results/fgc-1-pro19-sid2-pref1.json"
STORE_PATH = "runs/fgc-2-sf1/proto17/calibration"

MEMBER_KEYS = (
    "RK4-2049",
    "RK4-4097",
    "RK4-8193",
    "SSPRK3-4097",
    "SSPRK3-8193",
    "SSPRK3-16385",
)
ARRAY_NAMES = (
    "u",
    "p",
    "q",
    "grid_coordinates",
    "tracer_labels",
    "tracer_positions",
    "tracer_proper_times",
    "event_proper_times",
    "event_fields",
)
ZERO_SHA256 = "0" * 64
CAMPAIGN_ID = "FGC-2-SF1-PROTO17-GR0-A3-CAL-290a65bcd6a2ba820683"
AUTHORIZATION_COMMIT = "c11ba422ce49ddd6f6175de1a9d2da675b97f2de"
PLAN_SHA256 = "f3b365978aeb2ae6063742807a3161bf5d6aa0c8392cd283fdec2a79d3e46060"
GEN7_PARENT_SHA256 = "6277b3be4ffe1bd79b858e2b4770c244b88b9971cbff67ffb1a40d40f14e2e63"
GEN7_CHECKPOINT_SHA256 = "c0ebfe09ef4058e117e46d34424edc8139419bfb4680cec979a1f5ab64d37be9"
GEN7_CHECKPOINT_RAW_SHA256 = "209f673ad78d96553e26965c7ecf00953b7e9b57ece0b14d00e92e20a0623b21"
SEQ6_CFL_SHA256 = "78e65c08380a4fc979269e21db314232fd3fdc0e428c7e259f2031e66bd30bf7"
SEQ7_REJECTION_SHA256 = "e23d81705f1038f8240249658d0160cd5ca8ea8e33969cb4b32de3f03d426473"
SEQ7_REJECTION_RAW_SHA256 = "d1d929378bb0d9ca5146a9f83acf2d4f6be986503720596a93a5ae2b5f81292f"
GEN8_CHECKPOINT_SHA256 = "6d2b6e37e7c7cfd4ce64504353d886099003f959ea2aaee1a14fc778ce14fa46"
GEN8_CHECKPOINT_RAW_SHA256 = "dca398efa21d699ed09060c2afd1a13b7b3672ca62573c24bd846a0ff6a3ef7d"
SEQ8_TRANSITION_SHA256 = "949b29e46081bf5a42c34461b3894419db86aeca6c733b7a2636c91c6cb2ee4a"
SEQ8_TRANSITION_RAW_SHA256 = "8facb8cb67e6b9a891f1a8cb63056cc0b00f142592085ef21d3c9f1ced7c6d3c"
GEN7_CURSOR_SHA256 = "dee74d3243c16d1f987ec7c20b7bf6da40f407507716b65a86249605f6be3dd0"
GEN8_CURSOR_SHA256 = "6d594999f6aaebb5d56f0fd6ff12499e2b22a47e3daf74ac8b13b15d273fa177"
GEN7_LEDGER_SHA256 = "2763399f9824fcb476dc2e617d35cf1510ef08f11e65757ed5ef7e7840fd4984"
GEN8_LEDGER_SHA256 = "bc8a19299a3eb7f9ff57939a9c9c3c59943387a6644fc2c3be903258527eb0bc"
PHYSICAL_MEMBER = "RK4-2049"
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
REJECTED_MACRO_WIDTH_HEX = "0x1.aaa9612df9000p-9"
PENDING_CAP_HEX = "0x1.aaa9612df9000p-10"
SUCCESSOR_MACRO_WIDTH_HEX = "0x1.aaa9612df8000p-10"

EXPECTED_MEMBER_DESCRIPTORS = {
    "RK4-2049": "6b0dd985ab95ac52950348df5079a6b8fc4f5e7a78cf95031237cb0155c49c56",
    "RK4-4097": "540ed681955f244c82b9f91538f9d797e5237b04fcdd7011826efcade53f9388",
    "RK4-8193": "53b8413090d5af595d2bbcbc58a7ead07a9ab8c37aa90847dc669bd3e4bfa058",
    "SSPRK3-4097": "d4376e1025c1e5579d2816c87fbd32b46c3768a72c6e7753dd1c83d3d8ec14d5",
    "SSPRK3-8193": "070427ff80ee72eb16396f6d651e020ad6423b5068da8388e5e6fed994608f32",
    "SSPRK3-16385": "e7e2d42df4aa77301fe7b42fc8c64158bd4b3327ad402d55f3bd5f9ae7969564",
}
EXPECTED_PHYSICAL_IDENTITY = {
    "member_key": PHYSICAL_MEMBER,
    "accepted_time_hex": ACCEPTED_TIME_HEX,
    "descriptor_sha256": EXPECTED_MEMBER_DESCRIPTORS[PHYSICAL_MEMBER],
    "semantic_payload_sha256": "0e48e0e2396ea7d9ef79d1641e6944a2dfdaa3716bb65e0b882b97156ff9372f",
    "raw_archive_sha256": "f47a4f533eb8e0ed9f8ea337d520fa6f4e0bb6c1b27a4ef8023ee78bc830f5a4",
    "evolution_state_sha256": "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a",
    "u_bytes_sha256": "856fc327d3163d038f9fdbfc660968509dd7360faa0bb9d158721eb3049c2fe4",
    "p_bytes_sha256": "32939115b86938cf8a813bbf9a50a8c90bbaceaee161fa65695deb435d1243b7",
    "q_bytes_sha256": "48ec07bcace3935173c06e7b41915cacb17f14e1c9d3500522d4eb59a205c4a3",
}
EXPECTED_RETRY_STATE = {
    "member_key": PHYSICAL_MEMBER,
    "cursor_sha256": GEN8_CURSOR_SHA256,
    "mode": "RETRY_PENDING",
    "owner": "temporal",
    "pending_cap_hex": PENDING_CAP_HEX,
    "successor_macro_width_hex": SUCCESSOR_MACRO_WIDTH_HEX,
    "rejected_macro_width_hex": REJECTED_MACRO_WIDTH_HEX,
    "cfl_current": 1,
    "cfl_total": 1,
    "temporal_retry_total": 1,
    "TDG6_ledger_sha256": GEN8_LEDGER_SHA256,
}
EXPECTED_PREDECESSOR = {
    "artifact_id": "FGC-1-PRO19-SID1-FRZ1",
    "authority_commit": "c83ef958283be5e36ca31dd1f13df4a1b74f8696",
    "config_path": "configs/fgc/fgc-1-pro19-sid1-frz1.toml",
    "config_sha256": "5019067dd56445b33f270a2b6531572ee2e8e3b2318dc5f82e70d0f554b3c394",
    "result_path": "results/fgc-1-pro19-sid1-frz1.json",
    "result_sha256": "e8f937f8513c0ccb5d7c02bb1d83c4320c7ef0ded26a5461a212abf766882465",
    "historical_launch_manifest_path": "configs/fgc/fgc-1-pro19-launch-authority.toml",
    "historical_launch_manifest_sha256": "c3084ce9024319cde337603a00f052a6c7aefb04ad4dbc8cf475f5b7997fbd7d",
    "original_authorization_commit": AUTHORIZATION_COMMIT,
    "original_plan_sha256": PLAN_SHA256,
    "campaign_id": CAMPAIGN_ID,
}
EXPECTED_STORE = {
    "path": STORE_PATH,
    "regular_leaf_count": 52,
    "lexical_inventory_sha256": "8196b088a19f70ae1016a3d9dddb508add1ce1ce5cd227d105c954bd0c151c27",
    "nonfollowing_traversal_required": True,
    "active_writer_absent": True,
    "terminal_lock_absent": True,
    "staging_absent": True,
    "orphan_absent": True,
    "fork_absent": True,
    "symlink_and_special_file_absent": True,
}
EXPECTED_GENERATION7 = {
    "generation": 7,
    "checkpoint_sha256": GEN7_CHECKPOINT_SHA256,
    "checkpoint_raw_sha256": GEN7_CHECKPOINT_RAW_SHA256,
    "journal_sequence": 7,
    "journal_sha256": SEQ7_REJECTION_SHA256,
    "journal_raw_sha256": SEQ7_REJECTION_RAW_SHA256,
}
EXPECTED_GENERATION8 = {
    "generation": 8,
    "checkpoint_sha256": GEN8_CHECKPOINT_SHA256,
    "checkpoint_raw_sha256": GEN8_CHECKPOINT_RAW_SHA256,
    "parent_checkpoint_sha256": GEN7_CHECKPOINT_SHA256,
    "journal_sequence": 8,
    "journal_sha256": SEQ8_TRANSITION_SHA256,
    "journal_raw_sha256": SEQ8_TRANSITION_RAW_SHA256,
    "journal_parent_sha256": SEQ7_REJECTION_SHA256,
    "event": 23,
    "target_rational": "3/2",
    "disposition": "nonterminal",
    "suffix_classification": "clean",
}
EXPECTED_SCOPE = {
    "raw_store_mutation": False,
    "writer_lease_acquisition": False,
    "PDE_proposal_execution": False,
    "GEN0_reimport": False,
    "trajectory_resume": False,
    "candidate_branch_opened": False,
}
EXPECTED_NONCLAIMS = [
    "SID2 binds one completed metadata-only recovery; it does not advance a physical state or replay the rejected proposal.",
    "SID2 does not authorize trajectory resume by itself; a later committed-image authority and read-only resume preflight remain mandatory.",
    "SID2 does not complete a common event or GR-0 calibration and does not authorize SGB-L, FGC-QR, DEF1, retained-EFT, transition, or physical claims.",
]

# Updated by the integration contract before sealing.  In particular, the
# binder freezes the recovery boundary; it does not freeze a resume contract.
EXPECTED_CLAIMS = {
    "post_recovery_boundary_bound": True,
    "metadata_only_recovery_bound": True,
    "physical_state_advanced": False,
    "rejected_proposal_replayed": False,
    "recovery_boundary_frozen": True,
    "trajectory_resume_authorized": False,
    "safe_to_resume_trajectory": False,
    "committed_resume_authority_present": False,
    "committed_resume_image_authenticated": False,
    "committed_resume_image_authentication_required": True,
    "resume_authorized_by_artifact_alone": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_result_earned": False,
}

_EXPECTED_DIRECTORIES = (
    "checkpoints",
    "journal",
    "locks",
    "payloads",
    "receipts",
    "states",
)
_CHECKPOINT_NAME = re.compile(r"(?P<number>[0-9]{20})-(?P<digest>[0-9a-f]{64})\.json\Z")
_JOURNAL_NAME = re.compile(r"(?P<number>[0-9]{20})-(?P<digest>[0-9a-f]{64})\.journal\Z")
_STATE_NAME = re.compile(r"(?P<digest>[0-9a-f]{64})\.json\Z")
_PAYLOAD_NAME = re.compile(r"(?P<digest>[0-9a-f]{64})\.npz\Z")
_QUARANTINE_NAME = re.compile(r"\.active-write\.lock\.hlt16-quarantine-(?P<digest>[0-9a-f]{64})\Z")
_CHECKPOINT_FIELDS = frozenset({
    "authorization_commit", "campaign_id", "checkpoint_sha256", "disposition",
    "event", "generation", "journal_sequence", "journal_tip_sha256", "members",
    "parent_sha256", "plan_sha256", "protocol", "target", "terminal",
})
_JOURNAL_FIELDS = frozenset({
    "schema", "sequence", "campaign_id", "generation", "event",
    "previous_record_sha256", "kind", "payload", "record_sha256",
})
_CURSOR_FIELDS = frozenset({
    "protocol_artifact_id", "campaign_id", "member_key", "method", "point_count",
    "committed_common_event_index", "accepted_boundary_time", "accepted_state_sha256",
    "previous_step_index", "previous_transaction_serial", "TDG6_ledger_sha256",
    "cursor_generation", "event_target_time", "attempt_serial", "attempt_id", "mode",
    "retry_successor_payload_or_none", "journal_tip_sha256", "cursor_chain_parent_sha256",
    "cursor_chain_sha256",
})
_RETRY_FIELDS = frozenset({
    "predecessor_plan", "successor_plan", "TDG6_rejection_evidence",
    "complete_rejection_prefix", "durable_rejection_record_sha256",
    "predecessor_journal_tip_sha256", "predecessor_attempt_id",
    "successor_attempt_id", "member_key", "method", "point_count",
    "accepted_state_sha256", "accepted_boundary_time", "previous_step_index",
    "previous_transaction_serial", "retry_count", "half_cap", "minimum_macro_step",
})
_MAX_LEAF_BYTES = 96 * 1024 * 1024
_MAX_TREE_BYTES = 384 * 1024 * 1024
_MAX_ARRAY_MEMBER_BYTES = 32 * 1024 * 1024
_MAX_TOTAL_UNCOMPRESSED_BYTES = 96 * 1024 * 1024


class SID2BinderError(ValueError):
    """The SID2 config, compact record, or post-recovery store differs."""

    def __init__(self, stop_id: str, detail: str) -> None:
        super().__init__(f"{stop_id}: {detail}")
        self.stop_id = stop_id
        self.detail = detail


def _stop(stop_id: str, detail: str) -> None:
    raise SID2BinderError(stop_id, detail)


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
        _stop("SID2_CANONICAL_DRIFT", "value is not canonical JSON")
        raise AssertionError from error


def canonical_result(value: object) -> bytes:
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
        _stop("SID2_COMPACT_DRIFT", "result is not canonical JSON")
        raise AssertionError from error


def _digest(value: object) -> str:
    return sha256(_canonical(value)).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        _stop("SID2_HASH_DRIFT", f"{label} is not SHA-256")
    try:
        int(value, 16)
    except ValueError:
        _stop("SID2_HASH_DRIFT", f"{label} is not SHA-256")
    return value


def _json(raw: bytes, label: str) -> dict[str, Any]:
    def pairs(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
        answer: dict[str, Any] = {}
        for key, value in items:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer

    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=pairs)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        _stop("SID2_CANONICAL_DRIFT", f"{label} is malformed")
        raise AssertionError from error
    if not isinstance(value, dict) or _canonical(value) != raw:
        _stop("SID2_CANONICAL_DRIFT", f"{label} is noncanonical")
    return value


def _object_address(value: Mapping[str, Any], field: str, address: str, label: str) -> None:
    if value.get(field) != address:
        _stop("SID2_HASH_DRIFT", f"{label} embedded address differs")
    body = dict(value)
    body.pop(field, None)
    if _digest(body) != address:
        _stop("SID2_HASH_DRIFT", f"{label} content address differs")


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
        _stop("SID2_PATH_UNSAFE", f"{label} root is absent")
        raise AssertionError from error
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        _stop("SID2_PATH_UNSAFE", f"{label} root is unsafe")
    try:
        current = os.open(root, _directory_flags())
    except OSError as error:
        _stop("SID2_PATH_UNSAFE", f"{label} root cannot be opened")
        raise AssertionError from error
    try:
        for component in parts:
            before = os.stat(component, dir_fd=current, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                _stop("SID2_PATH_UNSAFE", f"{label} ancestor is unsafe: {component}")
            child = os.open(component, _directory_flags(), dir_fd=current)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                os.close(child)
                _stop("SID2_PATH_UNSAFE", f"{label} ancestor raced: {component}")
            os.close(current)
            current = child
        return current
    except OSError as error:
        os.close(current)
        _stop("SID2_PATH_UNSAFE", f"{label} ancestor cannot be opened safely")
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
        _stop("SID2_PATH_UNSAFE", f"{label} is not a single-linked bounded regular file")
    descriptor = -1
    try:
        descriptor = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent,
        )
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
            _stop("SID2_PATH_UNSAFE", f"{label} raced before read")
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _stop("SID2_PATH_UNSAFE", f"{label} ended early")
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _stop("SID2_PATH_UNSAFE", f"{label} grew during read")
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
            _stop("SID2_PATH_UNSAFE", f"{label} changed during read")
        return raw
    except OSError as error:
        _stop("SID2_PATH_UNSAFE", f"{label} cannot be read safely")
        raise AssertionError from error
    finally:
        if descriptor != -1:
            os.close(descriptor)


def _scan_tree_fd(root_fd: int) -> tuple[dict[str, _Leaf], tuple[str, ...]]:
    leaves: dict[str, _Leaf] = {}
    directories: list[str] = []
    total = 0

    def descend(directory_fd: int, prefix: str) -> None:
        nonlocal total
        try:
            with os.scandir(directory_fd) as scan:
                entries = sorted(scan, key=lambda entry: entry.name)
        except OSError as error:
            _stop("SID2_PATH_UNSAFE", "store cannot be scanned safely")
            raise AssertionError from error
        for entry in entries:
            name = entry.name
            relative = f"{prefix}/{name}" if prefix else name
            try:
                name.encode("utf-8")
                before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            except (UnicodeEncodeError, OSError) as error:
                _stop("SID2_PATH_UNSAFE", f"store node cannot be statted: {relative}")
                raise AssertionError from error
            if stat.S_ISLNK(before.st_mode):
                _stop("SID2_PATH_UNSAFE", f"symlink appears in store: {relative}")
            if stat.S_ISDIR(before.st_mode):
                try:
                    child = os.open(name, _directory_flags(), dir_fd=directory_fd)
                except OSError as error:
                    _stop("SID2_PATH_UNSAFE", f"directory cannot be opened: {relative}")
                    raise AssertionError from error
                try:
                    after = os.fstat(child)
                    if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                        _stop("SID2_PATH_UNSAFE", f"directory raced: {relative}")
                    directories.append(relative)
                    descend(child, relative)
                finally:
                    os.close(child)
                continue
            if not stat.S_ISREG(before.st_mode):
                _stop("SID2_PATH_UNSAFE", f"special file appears in store: {relative}")
            raw = _read_leaf_fd(directory_fd, name, before, relative)
            total += len(raw)
            if total > _MAX_TREE_BYTES:
                _stop("SID2_PATH_UNSAFE", "store exceeds the read-only byte budget")
            leaves[relative] = _Leaf(relative, len(raw), sha256(raw).hexdigest(), raw)

    descend(root_fd, "")
    return leaves, tuple(sorted(directories))


def _scan_store(repository: Path) -> tuple[dict[str, _Leaf], tuple[str, ...]]:
    parts = Path(STORE_PATH).parts
    parent = _open_directory_chain(Path(repository), parts[:-1], "repository/store")
    store_fd = -1
    try:
        try:
            with os.scandir(parent) as scan:
                siblings = sorted(scan, key=lambda entry: entry.name)
        except OSError as error:
            _stop("SID2_PATH_UNSAFE", "store parent cannot be scanned")
            raise AssertionError from error
        final = parts[-1]
        for entry in siblings:
            if entry.name.startswith(f".{final}.hlt15-stage-") or entry.name.startswith(
                f".{final}.hlt16-stage-"
            ):
                _stop("SID2_TREE_DRIFT", "campaign staging sibling is present")
        before = os.stat(final, dir_fd=parent, follow_symlinks=False)
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
            _stop("SID2_PATH_UNSAFE", "campaign store root is unsafe")
        store_fd = os.open(final, _directory_flags(), dir_fd=parent)
        after = os.fstat(store_fd)
        if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
            _stop("SID2_PATH_UNSAFE", "campaign store root raced")
        return _scan_tree_fd(store_fd)
    except FileNotFoundError as error:
        _stop("SID2_PATH_UNSAFE", "campaign store is absent")
        raise AssertionError from error
    except OSError as error:
        _stop("SID2_PATH_UNSAFE", "campaign store cannot be opened safely")
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
        _stop("SID2_TREE_DRIFT", "store directory grammar differs")
    for relative, leaf in leaves.items():
        parts = Path(relative).parts
        if len(parts) != 2 or parts[0] not in _EXPECTED_DIRECTORIES:
            _stop("SID2_TREE_DRIFT", f"nested or foreign leaf appears: {relative}")
        directory, name = parts
        if ".hlt16-stage-" in name or ".hlt15-stage-" in name:
            _stop("SID2_TREE_DRIFT", "publication staging residue is present")
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
            if name in {"active-write.lock", "terminal.lock"}:
                _stop("SID2_TREE_DRIFT", f"forbidden live lock is present: {name}")
            valid = name in {"bootstrap.guard", "writer.guard"} or _QUARANTINE_NAME.fullmatch(name) is not None
            if name in {"bootstrap.guard", "writer.guard"} and leaf.raw != b"":
                _stop("SID2_TREE_DRIFT", f"advisory guard bytes differ: {name}")
        if not valid:
            _stop("SID2_TREE_DRIFT", f"foreign store leaf appears: {relative}")


def _indexed_objects(
    leaves: Mapping[str, _Leaf],
    directory: str,
    pattern: re.Pattern[str],
    address_field: str,
) -> dict[int, dict[str, Any]]:
    answer: dict[int, dict[str, Any]] = {}
    for relative in sorted(path for path in leaves if path.startswith(f"{directory}/")):
        name = Path(relative).name
        match = pattern.fullmatch(name)
        if match is None:
            _stop("SID2_TREE_DRIFT", f"{directory} filename differs")
        index = int(match.group("number"))
        address = match.group("digest")
        if index in answer:
            _stop("SID2_TREE_DRIFT", f"{directory} fork appears at {index}")
        value = _json(leaves[relative].raw, f"{directory} {index}")
        _object_address(value, address_field, address, f"{directory} {index}")
        answer[index] = value
    return answer


def _validate_cursor(value: object, member_key: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != _CURSOR_FIELDS:
        _stop("SID2_RETRY_DRIFT", f"{member_key} cursor fields differ")
    cursor = dict(value)
    observed = _require_sha(cursor.get("cursor_chain_sha256"), f"{member_key} cursor")
    bare = dict(cursor)
    bare.pop("cursor_chain_sha256")
    if _digest(bare) != observed:
        _stop("SID2_RETRY_DRIFT", f"{member_key} cursor hash differs")
    try:
        method, points = member_key.split("-", 1)
        point_count = int(points)
    except ValueError:
        _stop("SID2_RETRY_DRIFT", f"{member_key} is malformed")
    if (
        cursor["member_key"] != member_key
        or cursor["method"] != method
        or cursor["point_count"] != point_count
        or cursor["campaign_id"] != CAMPAIGN_ID
        or cursor["protocol_artifact_id"] != TARGET_PROTOCOL
        or cursor["attempt_id"]
        != f"{CAMPAIGN_ID}:{member_key}:{cursor['committed_common_event_index']}:{cursor['attempt_serial']}"
        or cursor["mode"] not in {"FRESH_READY", "RETRY_PENDING"}
    ):
        _stop("SID2_RETRY_DRIFT", f"{member_key} cursor identity differs")
    for field in (
        "accepted_state_sha256",
        "TDG6_ledger_sha256",
        "journal_tip_sha256",
        "cursor_chain_parent_sha256",
    ):
        _require_sha(cursor[field], f"{member_key} {field}")
    retry = cursor["retry_successor_payload_or_none"]
    if cursor["mode"] == "FRESH_READY" and retry is not None:
        _stop("SID2_RETRY_DRIFT", f"{member_key} fresh cursor has retry payload")
    if cursor["mode"] == "RETRY_PENDING" and (
        not isinstance(retry, Mapping) or set(retry) != _RETRY_FIELDS
    ):
        _stop("SID2_RETRY_DRIFT", f"{member_key} pending cursor payload differs")
    return cursor


def _validate_runtime_checkpoints(
    checkpoints: Mapping[int, Mapping[str, Any]],
    journals: Mapping[int, Mapping[str, Any]],
) -> None:
    if set(checkpoints) != set(range(9)) or set(journals) != set(range(9)):
        _stop("SID2_TREE_DRIFT", "checkpoint/journal sequence has a gap or fork")
    genesis = checkpoints[0]
    if (
        genesis.get("campaign_generation") != 0
        or genesis.get("checkpoint_sha256")
        != "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d"
        or genesis.get("parent_checkpoint_sha256") != ZERO_SHA256
        or genesis.get("terminal_lock") is not False
    ):
        _stop("SID2_LINEAGE_DRIFT", "generation-zero checkpoint differs")

    previous = genesis["checkpoint_sha256"]
    for generation in range(1, 9):
        checkpoint = checkpoints[generation]
        if set(checkpoint) != _CHECKPOINT_FIELDS:
            _stop("SID2_LINEAGE_DRIFT", f"generation {generation} checkpoint fields differ")
        if (
            checkpoint["generation"] != generation
            or checkpoint["parent_sha256"] != previous
            or checkpoint["authorization_commit"] != AUTHORIZATION_COMMIT
            or checkpoint["plan_sha256"] != PLAN_SHA256
            or checkpoint["protocol"] != TARGET_PROTOCOL
            or checkpoint["campaign_id"] != CAMPAIGN_ID
            or checkpoint["event"] != 23
            or checkpoint["target"]
            != {"binary64_hex": "0x1.8000000000000p+0", "rational": "3/2"}
            or checkpoint["disposition"] != "nonterminal"
            or checkpoint["terminal"] is not None
            or set(checkpoint["members"]) != set(MEMBER_KEYS)
        ):
            _stop("SID2_LINEAGE_DRIFT", f"generation {generation} checkpoint relation differs")
        sequence = checkpoint["journal_sequence"]
        if sequence not in journals or checkpoint["journal_tip_sha256"] != journals[sequence]["record_sha256"]:
            _stop("SID2_LINEAGE_DRIFT", f"generation {generation} checkpoint tip differs")
        for key in MEMBER_KEYS:
            member = checkpoint["members"][key]
            if not isinstance(member, Mapping):
                _stop("SID2_LINEAGE_DRIFT", f"generation {generation} member differs")
            cursor = _validate_cursor(member.get("cursor"), key)
            ledger = member.get("ledger")
            if not isinstance(ledger, Mapping) or _digest(ledger) != cursor["TDG6_ledger_sha256"]:
                _stop("SID2_RETRY_DRIFT", f"generation {generation} {key} ledger differs")
            if (
                member.get("descriptor_sha256") != cursor["accepted_state_sha256"]
                or ledger.get("last_accepted_time_hex")
                != cursor.get("accepted_boundary_time", {}).get("binary64_hex")
            ):
                _stop("SID2_LINEAGE_DRIFT", f"generation {generation} {key} accepted state differs")
        previous = checkpoint["checkpoint_sha256"]

    parent = ZERO_SHA256
    for sequence in range(9):
        record = journals[sequence]
        if set(record) != _JOURNAL_FIELDS:
            _stop("SID2_LINEAGE_DRIFT", f"journal {sequence} fields differ")
        if (
            record["schema"] != "FGC-1-HLT16-campaign-journal-v1"
            or record["sequence"] != sequence
            or record["campaign_id"] != CAMPAIGN_ID
            or record["event"] != 23
            or record["previous_record_sha256"] != parent
        ):
            _stop("SID2_LINEAGE_DRIFT", f"journal {sequence} chain differs")
        parent = record["record_sha256"]


def _descriptor_objects(
    leaves: Mapping[str, _Leaf],
    checkpoints: Mapping[int, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    actual_states = {
        match.group("digest")
        for relative in leaves
        if relative.startswith("states/")
        for match in [_STATE_NAME.fullmatch(Path(relative).name)]
        if match is not None
    }
    genesis_states = checkpoints[0].get("states")
    if not isinstance(genesis_states, Mapping) or set(genesis_states) != set(MEMBER_KEYS):
        _stop("SID2_DESCRIPTOR_DRIFT", "generation-zero state set differs")
    legacy = {_require_sha(value, "generation-zero state") for value in genesis_states.values()}
    runtime = {
        _require_sha(checkpoints[generation]["members"][key]["descriptor_sha256"], "runtime descriptor")
        for generation in range(1, 9)
        for key in MEMBER_KEYS
    }
    if actual_states != legacy | runtime or legacy & runtime:
        _stop("SID2_TREE_DRIFT", "state descriptor orphan or omission appears")
    for address in legacy:
        leaf = leaves[f"states/{address}.json"]
        _json(leaf.raw, f"legacy state {address}")
        if leaf.raw_sha256 != address:
            _stop("SID2_DESCRIPTOR_DRIFT", "legacy state address differs")

    descriptors: dict[str, dict[str, Any]] = {}
    for address in runtime:
        leaf = leaves[f"states/{address}.json"]
        value = _json(leaf.raw, f"runtime descriptor {address}")
        if set(value) != {
            "schema", "semantic_sha256", "raw_archive_sha256", "arrays", "metadata",
            "descriptor_sha256",
        } or value["schema"] != "FGC-1-HLT16-member-state-v1":
            _stop("SID2_DESCRIPTOR_DRIFT", "runtime descriptor schema differs")
        _object_address(value, "descriptor_sha256", address, "runtime descriptor")
        _require_sha(value["semantic_sha256"], "semantic payload")
        _require_sha(value["raw_archive_sha256"], "raw archive")
        descriptors[address] = value
    return descriptors


def _decode_npz(raw: bytes, descriptor: Mapping[str, Any], label: str) -> dict[str, np.ndarray]:
    try:
        archive = zipfile.ZipFile(BytesIO(raw), "r")
        infos = archive.infolist()
    except (OSError, zipfile.BadZipFile) as error:
        _stop("SID2_PAYLOAD_DRIFT", f"{label} ZIP is malformed")
        raise AssertionError from error
    arrays: dict[str, np.ndarray] = {}
    try:
        expected = tuple(f"{name}.npy" for name in ARRAY_NAMES)
        if tuple(info.filename for info in infos) != expected or len({info.filename for info in infos}) != len(infos):
            _stop("SID2_PAYLOAD_DRIFT", f"{label} ZIP inventory/order differs")
        total = 0
        for info in infos:
            mode = (info.external_attr >> 16) & 0o170000
            if (
                info.is_dir()
                or info.flag_bits & 0x1
                or info.file_size <= 0
                or info.file_size > _MAX_ARRAY_MEMBER_BYTES + (1 << 20)
                or info.compress_size < 0
                or mode not in {0, stat.S_IFREG}
                or Path(info.filename).name != info.filename
            ):
                _stop("SID2_PAYLOAD_DRIFT", f"{label} ZIP member is unsafe")
            total += info.file_size
        if total > _MAX_TOTAL_UNCOMPRESSED_BYTES + len(infos) * (1 << 20):
            _stop("SID2_PAYLOAD_DRIFT", f"{label} ZIP expansion budget differs")
        for name, info in zip(ARRAY_NAMES, infos, strict=True):
            encoded = archive.read(info)
            if len(encoded) != info.file_size:
                _stop("SID2_PAYLOAD_DRIFT", f"{label} NPY size differs")
            try:
                value = np.load(BytesIO(encoded), allow_pickle=False)
            except (OSError, ValueError, EOFError) as error:
                _stop("SID2_PAYLOAD_DRIFT", f"{label} NPY member is malformed")
                raise AssertionError from error
            if (
                not isinstance(value, np.ndarray)
                or value.dtype != np.dtype("<f8")
                or not value.flags.c_contiguous
                or not value.shape
                or not np.isfinite(value).all()
            ):
                _stop("SID2_PAYLOAD_DRIFT", f"{label} NPY dtype/layout/value differs")
            arrays[name] = value.copy(order="C")
    except (RuntimeError, zipfile.BadZipFile) as error:
        _stop("SID2_PAYLOAD_DRIFT", f"{label} ZIP member cannot be decoded")
        raise AssertionError from error
    finally:
        archive.close()

    manifest = []
    for name in ARRAY_NAMES:
        value = arrays[name]
        manifest.append({
            "name": name,
            "dtype": "<f8",
            "shape": list(value.shape),
            "order": "C",
            "byte_count": int(value.nbytes),
            "bytes_sha256": sha256(value.tobytes(order="C")).hexdigest(),
        })
    if manifest != descriptor.get("arrays") or _digest(manifest) != descriptor.get("semantic_sha256"):
        _stop("SID2_PAYLOAD_DRIFT", f"{label} array manifest differs")
    return arrays


def _load_payloads(
    leaves: Mapping[str, _Leaf],
    descriptors: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, np.ndarray]]:
    semantics = {str(value["semantic_sha256"]) for value in descriptors.values()}
    actual = {
        match.group("digest")
        for relative in leaves
        if relative.startswith("payloads/")
        for match in [_PAYLOAD_NAME.fullmatch(Path(relative).name)]
        if match is not None
    }
    if actual != semantics:
        _stop("SID2_TREE_DRIFT", "payload orphan or omission appears")
    by_semantic: dict[str, dict[str, np.ndarray]] = {}
    expected_raw: dict[str, str] = {}
    for address, descriptor in descriptors.items():
        semantic = str(descriptor["semantic_sha256"])
        raw_hash = str(descriptor["raw_archive_sha256"])
        if semantic in expected_raw and expected_raw[semantic] != raw_hash:
            _stop("SID2_PAYLOAD_DRIFT", "shared semantic payload has conflicting raw hashes")
        expected_raw[semantic] = raw_hash
        if semantic not in by_semantic:
            leaf = leaves[f"payloads/{semantic}.npz"]
            if leaf.raw_sha256 != raw_hash:
                _stop("SID2_PAYLOAD_DRIFT", f"payload raw hash differs: {semantic}")
            by_semantic[semantic] = _decode_npz(leaf.raw, descriptor, f"payload {semantic}")
        else:
            arrays = by_semantic[semantic]
            manifest = [
                {
                    "name": name,
                    "dtype": "<f8",
                    "shape": list(arrays[name].shape),
                    "order": "C",
                    "byte_count": int(arrays[name].nbytes),
                    "bytes_sha256": sha256(arrays[name].tobytes(order="C")).hexdigest(),
                }
                for name in ARRAY_NAMES
            ]
            if manifest != descriptor.get("arrays"):
                _stop("SID2_PAYLOAD_DRIFT", f"descriptor {address} shares incompatible payload")
    return by_semantic


def _evolution_state_sha256(*arrays: np.ndarray) -> str:
    digestor = sha256()
    for value in arrays:
        normalized = np.asarray(value, dtype="<f8", order="C")
        digestor.update(struct.pack("<I", normalized.ndim))
        digestor.update(struct.pack(f"<{normalized.ndim}Q", *normalized.shape))
        digestor.update(normalized.tobytes(order="C"))
    return digestor.hexdigest()


def _bind_recovery_edge(
    checkpoints: Mapping[int, Mapping[str, Any]],
    journals: Mapping[int, Mapping[str, Any]],
    descriptors: Mapping[str, Mapping[str, Any]],
    payloads: Mapping[str, Mapping[str, np.ndarray]],
) -> tuple[dict[str, str], dict[str, Any]]:
    before, after = checkpoints[7], checkpoints[8]
    rejection, transition = journals[7], journals[8]
    if (
        before["checkpoint_sha256"] != GEN7_CHECKPOINT_SHA256
        or before["parent_sha256"] != GEN7_PARENT_SHA256
        or before["journal_sequence"] != 6
        or before["journal_tip_sha256"] != SEQ6_CFL_SHA256
        or after["checkpoint_sha256"] != GEN8_CHECKPOINT_SHA256
        or after["parent_sha256"] != GEN7_CHECKPOINT_SHA256
        or after["journal_sequence"] != 8
        or after["journal_tip_sha256"] != SEQ8_TRANSITION_SHA256
    ):
        _stop("SID2_LINEAGE_DRIFT", "generation-seven/eight checkpoint identity differs")
    if (
        rejection["generation"] != 8
        or rejection["kind"] != "tdg6_rejection"
        or rejection["previous_record_sha256"] != SEQ6_CFL_SHA256
        or rejection["record_sha256"] != SEQ7_REJECTION_SHA256
        or transition["generation"] != 8
        or transition["kind"] != "cursor_transition"
        or transition["previous_record_sha256"] != SEQ7_REJECTION_SHA256
        or transition["record_sha256"] != SEQ8_TRANSITION_SHA256
    ):
        _stop("SID2_LINEAGE_DRIFT", "generation-eight journal edge differs")

    before_descriptors = {key: before["members"][key]["descriptor_sha256"] for key in MEMBER_KEYS}
    after_descriptors = {key: after["members"][key]["descriptor_sha256"] for key in MEMBER_KEYS}
    if before_descriptors != EXPECTED_MEMBER_DESCRIPTORS or after_descriptors != before_descriptors:
        _stop("SID2_DESCRIPTOR_DRIFT", "six-member descriptor set advanced or differs")
    for key in MEMBER_KEYS:
        if key != PHYSICAL_MEMBER and before["members"][key] != after["members"][key]:
            _stop("SID2_LINEAGE_DRIFT", f"unchanged member moved during recovery: {key}")

    old_member, new_member = before["members"][PHYSICAL_MEMBER], after["members"][PHYSICAL_MEMBER]
    old_cursor = _validate_cursor(old_member["cursor"], PHYSICAL_MEMBER)
    new_cursor = _validate_cursor(new_member["cursor"], PHYSICAL_MEMBER)
    old_ledger, new_ledger = old_member["ledger"], new_member["ledger"]
    rejection_payload, transition_payload = rejection["payload"], transition["payload"]
    if not isinstance(rejection_payload, Mapping) or not isinstance(transition_payload, Mapping):
        _stop("SID2_RETRY_DRIFT", "recovery record payload differs")
    if (
        old_cursor["cursor_chain_sha256"] != GEN7_CURSOR_SHA256
        or new_cursor["cursor_chain_sha256"] != GEN8_CURSOR_SHA256
        or new_cursor["cursor_chain_parent_sha256"] != GEN7_CURSOR_SHA256
        or new_cursor["cursor_generation"] != old_cursor["cursor_generation"] + 1
        or new_cursor["attempt_serial"] != old_cursor["attempt_serial"] + 1
        or old_cursor["accepted_state_sha256"] != new_cursor["accepted_state_sha256"]
        or old_cursor["accepted_boundary_time"] != new_cursor["accepted_boundary_time"]
        or new_cursor["journal_tip_sha256"] != SEQ7_REJECTION_SHA256
        or new_cursor["mode"] != "RETRY_PENDING"
        or new_member["descriptor_sha256"] != old_member["descriptor_sha256"]
        or new_member["source_current"] != old_member["source_current"]
        or new_member["source_total"] != old_member["source_total"]
        or new_member["cfl_current"] != old_member["cfl_current"]
        or new_member["cfl_total"] != old_member["cfl_total"]
        or new_member["pending_owner"] != "temporal"
        or new_member["pending_cap_hex"] != PENDING_CAP_HEX
    ):
        _stop("SID2_RETRY_DRIFT", "recovery cursor/state transition differs")
    if _digest(old_ledger) != GEN7_LEDGER_SHA256 or _digest(new_ledger) != GEN8_LEDGER_SHA256:
        _stop("SID2_RETRY_DRIFT", "recovery ledger hash differs")
    evidence = rejection_payload.get("evidence")
    retry = new_cursor["retry_successor_payload_or_none"]
    if not isinstance(evidence, Mapping) or not isinstance(retry, Mapping):
        _stop("SID2_RETRY_DRIFT", "TDG6 retry evidence differs")
    evidence_text = _canonical(evidence).decode("ascii")
    if (
        rejection_payload.get("member_key") != PHYSICAL_MEMBER
        or rejection_payload.get("predecessor_descriptor_sha256")
        != EXPECTED_MEMBER_DESCRIPTORS[PHYSICAL_MEMBER]
        or rejection_payload.get("predecessor_cursor_sha256") != GEN7_CURSOR_SHA256
        or transition_payload.get("member_key") != PHYSICAL_MEMBER
        or transition_payload.get("predecessor_cursor_sha256") != GEN7_CURSOR_SHA256
        or transition_payload.get("rejection_record_sha256") != SEQ7_REJECTION_SHA256
        or transition_payload.get("successor_cursor_sha256") != GEN8_CURSOR_SHA256
        or transition_payload.get("successor_cursor") != new_cursor
        or retry["TDG6_rejection_evidence"] != evidence
        or retry["complete_rejection_prefix"] != [evidence_text]
        or retry["durable_rejection_record_sha256"] != SEQ7_REJECTION_SHA256
        or retry["predecessor_journal_tip_sha256"] != SEQ6_CFL_SHA256
        or retry["accepted_state_sha256"] != EXPECTED_MEMBER_DESCRIPTORS[PHYSICAL_MEMBER]
        or retry["accepted_boundary_time"]["binary64_hex"] != ACCEPTED_TIME_HEX
        or retry["retry_count"] != 1
        or retry["half_cap"] != PENDING_CAP_HEX
        or retry["predecessor_plan"]["macro_width_hex"] != REJECTED_MACRO_WIDTH_HEX
        or retry["successor_plan"]["macro_width_hex"] != SUCCESSOR_MACRO_WIDTH_HEX
    ):
        _stop("SID2_RETRY_DRIFT", "durable TDG6 retry relation differs")
    if (
        evidence.get("initial_state_sha256") != EXPECTED_PHYSICAL_IDENTITY["evolution_state_sha256"]
        or not isinstance(evidence.get("initial_time"), float)
        or float(evidence["initial_time"]).hex() != ACCEPTED_TIME_HEX
        or not isinstance(evidence.get("attempted_macro_step_size"), float)
        or float(evidence["attempted_macro_step_size"]).hex() != REJECTED_MACRO_WIDTH_HEX
        or evidence.get("retry_count_for_current_macro_step") != 1
        or evidence.get("cumulative_temporal_retry_count") != 1
        or evidence.get("failed_channels") != ["u:R"]
        or evidence.get("accepted_state_bitwise_preserved") is not True
        or evidence.get("accumulated_debit_unchanged") is not True
    ):
        _stop("SID2_RETRY_DRIFT", "TDG6 rejection no-commit evidence differs")
    if (
        new_ledger.get("accepted_macro_step_count") != old_ledger.get("accepted_macro_step_count")
        or new_ledger.get("accumulated_debit_vector_hex") != old_ledger.get("accumulated_debit_vector_hex")
        or new_ledger.get("last_accepted_time_hex") != old_ledger.get("last_accepted_time_hex")
        or old_ledger.get("current_macro_step_temporal_retry_count") != 0
        or old_ledger.get("cumulative_temporal_retry_count") != 0
        or new_ledger.get("current_macro_step_temporal_retry_count") != 1
        or new_ledger.get("cumulative_temporal_retry_count") != 1
        or new_ledger.get("last_accepted_macro_step_temporal_retry_count") != 0
        or new_ledger.get("serialized_temporal_rejections") != [evidence_text]
    ):
        _stop("SID2_RETRY_DRIFT", "TDG6 retry ledger update differs")

    descriptor = descriptors[EXPECTED_MEMBER_DESCRIPTORS[PHYSICAL_MEMBER]]
    semantic = str(descriptor["semantic_sha256"])
    arrays = payloads[semantic]
    physical = {
        "member_key": PHYSICAL_MEMBER,
        "accepted_time_hex": new_cursor["accepted_boundary_time"]["binary64_hex"],
        "descriptor_sha256": new_member["descriptor_sha256"],
        "semantic_payload_sha256": semantic,
        "raw_archive_sha256": descriptor["raw_archive_sha256"],
        "evolution_state_sha256": _evolution_state_sha256(arrays["u"], arrays["p"], arrays["q"]),
        "u_bytes_sha256": sha256(arrays["u"].tobytes(order="C")).hexdigest(),
        "p_bytes_sha256": sha256(arrays["p"].tobytes(order="C")).hexdigest(),
        "q_bytes_sha256": sha256(arrays["q"].tobytes(order="C")).hexdigest(),
    }
    if physical != EXPECTED_PHYSICAL_IDENTITY:
        _stop("SID2_DESCRIPTOR_DRIFT", "preserved physical u,p,q identity differs")
    retry_state = {
        "member_key": PHYSICAL_MEMBER,
        "cursor_sha256": new_cursor["cursor_chain_sha256"],
        "mode": new_cursor["mode"],
        "owner": new_member["pending_owner"],
        "pending_cap_hex": new_member["pending_cap_hex"],
        "successor_macro_width_hex": retry["successor_plan"]["macro_width_hex"],
        "rejected_macro_width_hex": retry["predecessor_plan"]["macro_width_hex"],
        "cfl_current": new_member["cfl_current"],
        "cfl_total": new_member["cfl_total"],
        "temporal_retry_total": new_ledger["cumulative_temporal_retry_count"],
        "TDG6_ledger_sha256": new_cursor["TDG6_ledger_sha256"],
    }
    if retry_state != EXPECTED_RETRY_STATE:
        _stop("SID2_RETRY_DRIFT", "post-recovery retry state differs")
    return physical, retry_state


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        _stop("SID2_CONFIG_DRIFT", "SID2 config is malformed")
        raise AssertionError from error
    if not isinstance(value, dict) or set(value) != {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "classification", "nonclaims", "predecessor", "store", "generation7",
        "generation8", "member_descriptors", "physical_identity", "retry_state",
        "scope", "claims",
    }:
        _stop("SID2_CONFIG_DRIFT", "SID2 config fields differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != EXPECTED_NONCLAIMS
        or value["predecessor"] != EXPECTED_PREDECESSOR
        or value["store"] != EXPECTED_STORE
        or value["generation7"] != EXPECTED_GENERATION7
        or value["generation8"] != EXPECTED_GENERATION8
        or value["physical_identity"] != EXPECTED_PHYSICAL_IDENTITY
        or value["retry_state"] != EXPECTED_RETRY_STATE
        or value["scope"] != EXPECTED_SCOPE
        or value["claims"] != EXPECTED_CLAIMS
    ):
        _stop("SID2_CONFIG_DRIFT", "SID2 exact contract differs")
    configured_members = value["member_descriptors"]
    expected_config_members = {key.replace("-", "_"): digest for key, digest in EXPECTED_MEMBER_DESCRIPTORS.items()}
    if configured_members != expected_config_members:
        _stop("SID2_CONFIG_DRIFT", "SID2 six-member descriptor contract differs")
    return value


def expected_anchor(config_raw: bytes) -> dict[str, Any]:
    """Construct the exact expected live anchor without opening the live store."""
    config = _config(config_raw)
    return {
        "store": dict(config["store"]),
        "generation7": dict(config["generation7"]),
        "generation8": dict(config["generation8"]),
        "member_descriptors": {
            key: config["member_descriptors"][key.replace("-", "_")] for key in MEMBER_KEYS
        },
        "physical_identity": dict(config["physical_identity"]),
        "retry_state": dict(config["retry_state"]),
    }


def bind_live_store(repository: Path) -> dict[str, Any]:
    """Independently authenticate the exact clean generation-eight boundary."""
    leaves, directories = _scan_store(Path(repository))
    _validate_tree_grammar(leaves, directories)
    checkpoints = _indexed_objects(leaves, "checkpoints", _CHECKPOINT_NAME, "checkpoint_sha256")
    journals = _indexed_objects(leaves, "journal", _JOURNAL_NAME, "record_sha256")
    _validate_runtime_checkpoints(checkpoints, journals)
    descriptors = _descriptor_objects(leaves, checkpoints)
    payloads = _load_payloads(leaves, descriptors)
    physical, retry_state = _bind_recovery_edge(checkpoints, journals, descriptors, payloads)

    for relative, leaf in leaves.items():
        if relative.startswith("locks/.active-write.lock.hlt16-quarantine-"):
            match = _QUARANTINE_NAME.fullmatch(Path(relative).name)
            assert match is not None
            if leaf.raw_sha256 != match.group("digest"):
                _stop("SID2_TREE_DRIFT", "retired writer content address differs")
            retired = _json(leaf.raw, "retired writer lease")
            if retired.get("schema") != "FGC-1-HLT16-active-writer-v1":
                _stop("SID2_TREE_DRIFT", "retired writer schema differs")

    inventory = _inventory_sha256(leaves)
    if len(leaves) != EXPECTED_STORE["regular_leaf_count"] or inventory != EXPECTED_STORE["lexical_inventory_sha256"]:
        _stop(
            "SID2_TREE_DRIFT",
            f"lexical inventory differs: leaves={len(leaves)} sha256={inventory}",
        )
    if (
        leaves[f"checkpoints/{7:020d}-{GEN7_CHECKPOINT_SHA256}.json"].raw_sha256
        != GEN7_CHECKPOINT_RAW_SHA256
        or leaves[f"checkpoints/{8:020d}-{GEN8_CHECKPOINT_SHA256}.json"].raw_sha256
        != GEN8_CHECKPOINT_RAW_SHA256
        or leaves[f"journal/{7:020d}-{SEQ7_REJECTION_SHA256}.journal"].raw_sha256
        != SEQ7_REJECTION_RAW_SHA256
        or leaves[f"journal/{8:020d}-{SEQ8_TRANSITION_SHA256}.journal"].raw_sha256
        != SEQ8_TRANSITION_RAW_SHA256
    ):
        _stop("SID2_HASH_DRIFT", "generation-seven/eight raw leaf hash differs")
    return {
        "store": dict(EXPECTED_STORE),
        "generation7": dict(EXPECTED_GENERATION7),
        "generation8": dict(EXPECTED_GENERATION8),
        "member_descriptors": dict(EXPECTED_MEMBER_DESCRIPTORS),
        "physical_identity": physical,
        "retry_state": retry_state,
    }


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
            "live_anchor": expected_anchor(config_raw),
            "scope": dict(config["scope"]),
            "claims": dict(config["claims"]),
            "nonclaims": list(config["nonclaims"]),
        },
    }


def build_sid2_result(config_raw: bytes, repository: Path) -> dict[str, Any]:
    """Build the deterministic result only after authenticating the live store."""
    config = _config(config_raw)
    observed = bind_live_store(Path(repository))
    if observed != expected_anchor(config_raw):
        _stop("SID2_LINEAGE_DRIFT", "live store/config cross-binding differs")
    return _result_from_config(config_raw, config)


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    """Validate tracked compact bytes without reopening the ignored live store."""
    config = _config(config_raw)
    try:
        value = json.loads(
            result_raw.decode("ascii"),
            object_pairs_hook=lambda items: _reject_result_duplicates(items),
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        _stop("SID2_COMPACT_DRIFT", "SID2 result is malformed")
        raise AssertionError from error
    if not isinstance(value, dict) or canonical_result(value) != result_raw:
        _stop("SID2_COMPACT_DRIFT", "SID2 result is noncanonical")
    if value != _result_from_config(config_raw, config):
        _stop("SID2_COMPACT_DRIFT", "SID2 compact result differs")
    return value


def _reject_result_duplicates(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError(key)
        answer[key] = value
    return answer


__all__ = [
    "ARTIFACT_ID",
    "CONFIG_PATH",
    "RESULT_PATH",
    "STORE_PATH",
    "SID2BinderError",
    "bind_live_store",
    "build_sid2_result",
    "canonical_result",
    "expected_anchor",
    "validate_compact_result",
]
