"""Externally identified RCV3 projection of the closed generation-nine lineage.

The original PROTO19 store is an immutable invalid-terminal history.  RCV3
does not rewrite, truncate, or resume that store.  It authenticates the exact
generation-ten terminal, selects the already closed generation-nine ancestor,
and installs a byte-identical fifty-leaf lineage projection in a separately
identified wrapper.  The HLT16 campaign id, plan hash, and authorization
commit inside the projected store remain unchanged.

This module contains no PDE proposal or candidate-branch entry point.  The
installed projection is merely a restartable finite-state premise for a later,
separately authorized runner.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import stat
import tempfile
from types import MappingProxyType
from typing import Callable, Iterable, Mapping

from . import proto17_hlt15_runtime as hlt15
from .hlt16_campaign_schema import CampaignCheckpoint
from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError


PROJECTION_ID = "FGC-1-TDG8-RCV3-PROJECTION-1"
SOURCE_RELATIVE = Path("runs/fgc-2-sf1/proto19/calibration")
WRAPPER_RELATIVE = Path("runs/fgc-2-sf1/tdg8-rcv3")
PROJECTED_STORE_RELATIVE = WRAPPER_RELATIVE / "calibration"
RECEIPT_RELATIVE = WRAPPER_RELATIVE / "bootstrap-receipt.json"

INTERNAL_CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"
INTERNAL_AUTHORIZATION_COMMIT = "6df67967e6e0bc63eca4dc6951a60d3787368f26"
INTERNAL_PLAN_SHA256 = (
    "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"
)
TERMINAL_GENERATION = 10
TERMINAL_CHECKPOINT_SHA256 = (
    "291eb466d5cb566b1f03822bf109c0c9e0f6ff930e9be496d26309b53b5005ea"
)
TERMINAL_JOURNAL_SEQUENCE = 11
TERMINAL_JOURNAL_SHA256 = (
    "394c803df15f0b73827425dec4c914a1dd998aad8de23d8a8ecf3d637e1cf340"
)
ANCHOR_GENERATION = 9
ANCHOR_CHECKPOINT_SHA256 = (
    "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56"
)
ANCHOR_JOURNAL_SEQUENCE = 10
ANCHOR_JOURNAL_SHA256 = (
    "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff"
)
EVENT = 23
TARGET = MappingProxyType(
    {"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"}
)

SOURCE_LEAF_COUNT = 59
PROJECTED_LEAF_COUNT = 50
SOURCE_TREE_SHA256 = (
    "b076b5233e5e2ee9898444f6cade495e2e0ed8d24697348f640485104d2cf145"
)
PROJECTED_TREE_SHA256 = (
    "d84e85c3dd29ac1554111ed00b37289d2b89273e7a7c2c47e3d21f5b9fda9e0b"
)
STORE_DIRECTORIES = (
    "checkpoints", "journal", "locks", "payloads", "receipts", "states",
)

EXCLUDED_SOURCE_PATHS = (
    "checkpoints/00000000000000000010-"
    f"{TERMINAL_CHECKPOINT_SHA256}.json",
    "journal/00000000000000000011-"
    f"{TERMINAL_JOURNAL_SHA256}.journal",
    "locks/.active-write.lock.hlt16-quarantine-"
    "1d1ebcfbcb57ecf91adb191815acd57c1137609cd50f35fa98e7d18d2c91086c",
    "locks/.active-write.lock.hlt16-quarantine-"
    "2f76eeac85746761d38617804da02f0d520676a86a243a5e29b2a9f0fdfa71e9",
    "locks/.active-write.lock.hlt16-quarantine-"
    "6c409d46d6c049cced6a26d9ad733911070dad9d53e966674a7fdce4ba372946",
    "locks/.active-write.lock.hlt16-quarantine-"
    "f474863d36b1c0312ff677560f1385c2913add7353bde694f08e69161bcef5ec",
    "locks/bootstrap.guard",
    "locks/terminal.lock",
    "locks/writer.guard",
)

_MAX_LEAF_BYTES = 64 * 1024 * 1024
_MAX_TREE_BYTES = 128 * 1024 * 1024
FaultHook = Callable[[str], None]


class TDG8RCV3ForkRuntimeError(ValueError):
    """The terminal source or its external recovery projection is unsafe."""


@dataclass(frozen=True, slots=True)
class RecoveryForkSpec:
    projection_id: str
    source_relative: str
    wrapper_relative: str
    projected_store_relative: str
    receipt_relative: str

    def validate(self) -> None:
        expected = (
            PROJECTION_ID,
            SOURCE_RELATIVE.as_posix(),
            WRAPPER_RELATIVE.as_posix(),
            PROJECTED_STORE_RELATIVE.as_posix(),
            RECEIPT_RELATIVE.as_posix(),
        )
        observed = (
            self.projection_id,
            self.source_relative,
            self.wrapper_relative,
            self.projected_store_relative,
            self.receipt_relative,
        )
        if observed != expected:
            raise TDG8RCV3ForkRuntimeError("RCV3 production projection spec differs")
        for value in observed[1:]:
            _safe_relative(value, "RCV3 projection path")


PRODUCTION_SPEC = RecoveryForkSpec(
    projection_id=PROJECTION_ID,
    source_relative=SOURCE_RELATIVE.as_posix(),
    wrapper_relative=WRAPPER_RELATIVE.as_posix(),
    projected_store_relative=PROJECTED_STORE_RELATIVE.as_posix(),
    receipt_relative=RECEIPT_RELATIVE.as_posix(),
)


@dataclass(frozen=True, slots=True)
class RecoveryForkAuthority:
    """External committed authority supplied by the frozen bootstrap wrapper.

    The runtime intentionally knows no eventual commit or artifact digest.
    The prospective bootstrap script authenticates those bytes and passes the
    resulting tuple here; this type enforces its canonical shape and binds it
    into the installed receipt.
    """

    authority_commit_sha: str
    frz1_result_sha256: str
    runtime_sha256: str
    bootstrap_script_sha256: str

    def validate(self) -> None:
        _hex(self.authority_commit_sha, 40, "RCV3 authority commit")
        _hex(self.frz1_result_sha256, 64, "RCV3 FRZ1 result")
        _hex(self.runtime_sha256, 64, "RCV3 runtime")
        _hex(self.bootstrap_script_sha256, 64, "RCV3 bootstrap script")

    def as_mapping(self) -> dict[str, str]:
        self.validate()
        return {
            "authority_commit_sha": self.authority_commit_sha,
            "frz1_result_sha256": self.frz1_result_sha256,
            "runtime_sha256": self.runtime_sha256,
            "bootstrap_script_sha256": self.bootstrap_script_sha256,
        }


@dataclass(frozen=True, slots=True)
class TreeLeaf:
    path: str
    raw_sha256: str
    size: int

    def as_mapping(self) -> dict[str, object]:
        return {
            "path": self.path,
            "byte_count": self.size,
            "sha256": self.raw_sha256,
        }


@dataclass(frozen=True, slots=True)
class _TreeSnapshot:
    leaves: tuple[TreeLeaf, ...]
    directories: tuple[str, ...]
    raw_by_path: Mapping[str, bytes]
    tree_sha256: str


@dataclass(frozen=True, slots=True)
class AuthenticatedRecoveryForkSource:
    spec: RecoveryForkSpec
    terminal_checkpoint: CampaignCheckpoint
    anchor_checkpoint: CampaignCheckpoint
    source_leaves: tuple[TreeLeaf, ...]
    projected_leaves: tuple[TreeLeaf, ...]
    raw_by_path: Mapping[str, bytes]
    source_tree_sha256: str
    projected_tree_sha256: str


@dataclass(frozen=True, slots=True)
class RecoveryForkInstallation:
    projection_id: str
    wrapper_relative: str
    projected_store_relative: str
    receipt_relative: str
    receipt_sha256: str
    anchor_checkpoint_sha256: str
    installed_now: bool


def canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG8RCV3ForkRuntimeError("RCV3 object is not canonical") from exc


def _safe_relative(value: object, label: str) -> tuple[str, ...]:
    if not isinstance(value, str):
        raise TDG8RCV3ForkRuntimeError(f"{label} is not a path")
    path = Path(value)
    if (
        path.is_absolute()
        or "\\" in value
        or not path.parts
        or path.as_posix() != value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise TDG8RCV3ForkRuntimeError(f"{label} is unsafe")
    return path.parts


def _hex(value: object, length: int, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != length
        or value.lower() != value
    ):
        raise TDG8RCV3ForkRuntimeError(f"{label} is not hexadecimal")
    try:
        int(value, 16)
    except ValueError as exc:
        raise TDG8RCV3ForkRuntimeError(f"{label} is not hexadecimal") from exc
    return value


def _repository(root: Path) -> Path:
    repository = Path(os.path.abspath(os.fspath(root)))
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = -1
    try:
        before = repository.lstat()
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
            raise TDG8RCV3ForkRuntimeError("repository root is unsafe")
        descriptor = os.open(repository, flags)
        after = os.fstat(descriptor)
        if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
            raise TDG8RCV3ForkRuntimeError("repository root raced")
    except OSError as exc:
        raise TDG8RCV3ForkRuntimeError("repository root cannot be opened safely") from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)
    return repository


def _open_relative_directory(root: Path, relative: str) -> int:
    parts = _safe_relative(relative, "tree root")
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    parent = -1
    try:
        root_before = root.lstat()
        if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
            raise TDG8RCV3ForkRuntimeError("tree repository root is unsafe")
        parent = os.open(root, flags)
        root_after = os.fstat(parent)
        if (root_before.st_dev, root_before.st_ino) != (
            root_after.st_dev,
            root_after.st_ino,
        ):
            raise TDG8RCV3ForkRuntimeError("tree repository root raced")
        for component in parts:
            before = os.stat(component, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise TDG8RCV3ForkRuntimeError("tree parent is unsafe")
            child = os.open(component, flags, dir_fd=parent)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                os.close(child)
                raise TDG8RCV3ForkRuntimeError("tree parent raced")
            os.close(parent)
            parent = child
        return parent
    except OSError as exc:
        if parent != -1:
            os.close(parent)
        raise TDG8RCV3ForkRuntimeError("tree cannot be opened safely") from exc
    except BaseException:
        if parent != -1:
            os.close(parent)
        raise


def _read_leaf(parent: int, name: str, relative: str) -> bytes:
    leaf = -1
    try:
        before = os.stat(name, dir_fd=parent, follow_symlinks=False)
        if (
            stat.S_ISLNK(before.st_mode)
            or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_size > _MAX_LEAF_BYTES
        ):
            raise TDG8RCV3ForkRuntimeError(f"unsafe RCV3 leaf: {relative}")
        leaf = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent,
        )
        active = os.fstat(leaf)
        identity = (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
        )
        if identity != (
            active.st_dev,
            active.st_ino,
            active.st_size,
            active.st_mtime_ns,
            active.st_ctime_ns,
        ):
            raise TDG8RCV3ForkRuntimeError(f"RCV3 leaf raced: {relative}")
        chunks: list[bytes] = []
        remaining = _MAX_LEAF_BYTES + 1
        while remaining > 0:
            block = os.read(leaf, min(1 << 20, remaining))
            if not block:
                break
            chunks.append(block)
            remaining -= len(block)
        raw = b"".join(chunks)
        after = os.fstat(leaf)
        path_after = os.stat(name, dir_fd=parent, follow_symlinks=False)
        after_identity = (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        )
        if (
            len(raw) != before.st_size
            or after_identity != identity
            or (path_after.st_dev, path_after.st_ino) != (before.st_dev, before.st_ino)
        ):
            raise TDG8RCV3ForkRuntimeError(f"RCV3 leaf changed during read: {relative}")
        return raw
    except OSError as exc:
        raise TDG8RCV3ForkRuntimeError(f"RCV3 leaf cannot be read: {relative}") from exc
    finally:
        if leaf != -1:
            os.close(leaf)


def _snapshot_tree(root: Path, relative: str) -> _TreeSnapshot:
    root_fd = _open_relative_directory(root, relative)
    leaves: list[TreeLeaf] = []
    directories: list[str] = []
    raw_by_path: dict[str, bytes] = {}
    total = 0
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)

    def descend(directory_fd: int, prefix: str) -> None:
        nonlocal total
        try:
            with os.scandir(directory_fd) as scan:
                entries = sorted(scan, key=lambda item: item.name)
        except OSError as exc:
            raise TDG8RCV3ForkRuntimeError("RCV3 tree cannot be scanned") from exc
        for entry in entries:
            name = entry.name
            if name in {"", ".", ".."} or "/" in name or "\\" in name:
                raise TDG8RCV3ForkRuntimeError("RCV3 tree name is unsafe")
            path = f"{prefix}/{name}" if prefix else name
            try:
                before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            except OSError as exc:
                raise TDG8RCV3ForkRuntimeError(f"RCV3 node raced: {path}") from exc
            if stat.S_ISLNK(before.st_mode):
                raise TDG8RCV3ForkRuntimeError(f"unsafe RCV3 node: {path}")
            if stat.S_ISREG(before.st_mode):
                raw = _read_leaf(directory_fd, name, path)
                total += len(raw)
                if total > _MAX_TREE_BYTES:
                    raise TDG8RCV3ForkRuntimeError("RCV3 tree exceeds its size bound")
                leaf = TreeLeaf(path, sha256(raw).hexdigest(), len(raw))
                leaves.append(leaf)
                raw_by_path[path] = raw
                continue
            if not stat.S_ISDIR(before.st_mode):
                raise TDG8RCV3ForkRuntimeError(f"unsafe RCV3 node: {path}")
            child = -1
            try:
                child = os.open(name, flags, dir_fd=directory_fd)
                after = os.fstat(child)
                if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                    raise TDG8RCV3ForkRuntimeError(f"RCV3 directory raced: {path}")
                directories.append(path)
                descend(child, path)
            except OSError as exc:
                raise TDG8RCV3ForkRuntimeError(
                    f"RCV3 directory cannot be opened: {path}"
                ) from exc
            finally:
                if child != -1:
                    os.close(child)

    try:
        descend(root_fd, "")
    finally:
        os.close(root_fd)
    ordered = tuple(sorted(leaves, key=lambda leaf: leaf.path))
    if tuple(leaves) != ordered:
        leaves = list(ordered)
    digest = sha256(canonical([leaf.as_mapping() for leaf in ordered])).hexdigest()
    return _TreeSnapshot(
        leaves=ordered,
        directories=tuple(sorted(directories)),
        raw_by_path=MappingProxyType(dict(raw_by_path)),
        tree_sha256=digest,
    )


def _source_identity(checkpoint: CampaignCheckpoint, *, terminal: bool) -> None:
    expected_generation = TERMINAL_GENERATION if terminal else ANCHOR_GENERATION
    expected_sha = TERMINAL_CHECKPOINT_SHA256 if terminal else ANCHOR_CHECKPOINT_SHA256
    expected_sequence = TERMINAL_JOURNAL_SEQUENCE if terminal else ANCHOR_JOURNAL_SEQUENCE
    expected_journal = TERMINAL_JOURNAL_SHA256 if terminal else ANCHOR_JOURNAL_SHA256
    if (
        checkpoint.generation != expected_generation
        or checkpoint.sha256 != expected_sha
        or checkpoint.journal_sequence != expected_sequence
        or checkpoint.journal_tip_sha256 != expected_journal
        or checkpoint.protocol != "FGC-2-SF1-PROTO18"
        or checkpoint.campaign_id != INTERNAL_CAMPAIGN_ID
        or checkpoint.authorization_commit != INTERNAL_AUTHORIZATION_COMMIT
        or checkpoint.plan_sha256 != INTERNAL_PLAN_SHA256
        or checkpoint.event != EVENT
        or checkpoint.target != dict(TARGET)
    ):
        raise TDG8RCV3ForkRuntimeError("RCV3 HLT16 checkpoint identity differs")
    if terminal:
        if (
            checkpoint.disposition != "invalid_terminal"
            or checkpoint.terminal
            != {
                "classification": "invalid",
                "event": EVENT,
                "member_key": "RK4-2049",
                "owner": "invalid",
                "reason": "invalid",
                "target": dict(TARGET),
            }
        ):
            raise TDG8RCV3ForkRuntimeError("RCV3 terminal checkpoint differs")
    elif (
        checkpoint.disposition != "nonterminal"
        or checkpoint.terminal is not None
    ):
        raise TDG8RCV3ForkRuntimeError("RCV3 anchor checkpoint differs")


def _authenticate_source_once(
    repository: Path, spec: RecoveryForkSpec
) -> AuthenticatedRecoveryForkSource:
    tree = _snapshot_tree(repository, spec.source_relative)
    if (
        tree.directories != STORE_DIRECTORIES
        or len(tree.leaves) != SOURCE_LEAF_COUNT
        or tree.tree_sha256 != SOURCE_TREE_SHA256
    ):
        raise TDG8RCV3ForkRuntimeError("RCV3 source tree identity differs")
    source_paths = {leaf.path for leaf in tree.leaves}
    if not set(EXCLUDED_SOURCE_PATHS).issubset(source_paths):
        raise TDG8RCV3ForkRuntimeError("RCV3 source exclusions are absent")
    projected = tuple(
        leaf for leaf in tree.leaves if leaf.path not in EXCLUDED_SOURCE_PATHS
    )
    if (
        len(projected) != PROJECTED_LEAF_COUNT
        or sha256(canonical([leaf.as_mapping() for leaf in projected])).hexdigest()
        != PROJECTED_TREE_SHA256
    ):
        raise TDG8RCV3ForkRuntimeError("RCV3 projected lineage differs")

    store = HLT16CampaignStore(repository / spec.source_relative)
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
        anchor = store.authenticated_checkpoint_at_generation(ANCHOR_GENERATION)
    except HLT16CampaignStoreError as exc:
        raise TDG8RCV3ForkRuntimeError("RCV3 source HLT16 validation failed") from exc
    terminal = snapshot.checkpoint
    _source_identity(terminal, terminal=True)
    _source_identity(anchor, terminal=False)
    if (
        snapshot.suffix_records
        or snapshot.suffix_kinds
        or snapshot.suffix_classification != "clean"
        or snapshot.orphan_descriptor_sha256 is not None
        or snapshot.orphan_payload_semantic_sha256 is not None
        or snapshot.staging_paths
        or snapshot.terminal_lock_present is not True
        or status.checkpoint_generation != TERMINAL_GENERATION
        or status.checkpoint_sha256 != TERMINAL_CHECKPOINT_SHA256
        or status.last_complete_journal_sequence != TERMINAL_JOURNAL_SEQUENCE
        or status.journal_tip_sha256 != TERMINAL_JOURNAL_SHA256
        or status.state != "terminal"
        or status.safe_to_restart is not False
        or status.terminal is not True
        or status.active_write is not False
        or status.writer_state is not None
    ):
        raise TDG8RCV3ForkRuntimeError("RCV3 source is not the exact closed terminal")
    return AuthenticatedRecoveryForkSource(
        spec=spec,
        terminal_checkpoint=terminal,
        anchor_checkpoint=anchor,
        source_leaves=tree.leaves,
        projected_leaves=projected,
        raw_by_path=tree.raw_by_path,
        source_tree_sha256=tree.tree_sha256,
        projected_tree_sha256=PROJECTED_TREE_SHA256,
    )


def authenticate_recovery_fork_source(
    root: Path, spec: RecoveryForkSpec = PRODUCTION_SPEC
) -> AuthenticatedRecoveryForkSource:
    """Authenticate the exact source twice and return its closed gen-nine cut."""
    repository = _repository(root)
    spec.validate()
    first = _authenticate_source_once(repository, spec)
    second = _authenticate_source_once(repository, spec)
    if (
        first.source_leaves != second.source_leaves
        or first.raw_by_path != second.raw_by_path
        or first.terminal_checkpoint != second.terminal_checkpoint
        or first.anchor_checkpoint != second.anchor_checkpoint
    ):
        raise TDG8RCV3ForkRuntimeError("RCV3 source changed during authentication")
    return first


def _receipt(
    source: AuthenticatedRecoveryForkSource,
    authority: RecoveryForkAuthority,
) -> dict[str, object]:
    body: dict[str, object] = {
        "schema": "FGC-1-TDG8-RCV3-external-recovery-fork-v1",
        "projection_id": source.spec.projection_id,
        "kind": "externally_identified_recovery_fork_projection",
        "source_terminal_store": source.spec.source_relative,
        "wrapper_root": source.spec.wrapper_relative,
        "projected_store": source.spec.projected_store_relative,
        "external_install_authority": authority.as_mapping(),
        "source_terminal": {
            "generation": TERMINAL_GENERATION,
            "checkpoint_sha256": TERMINAL_CHECKPOINT_SHA256,
            "journal_sequence": TERMINAL_JOURNAL_SEQUENCE,
            "journal_sha256": TERMINAL_JOURNAL_SHA256,
            "disposition": "invalid_terminal",
            "terminal_lock_authenticated": True,
            "writer_absent": True,
            "suffix_absent": True,
        },
        "selected_lineage_anchor": {
            "generation": ANCHOR_GENERATION,
            "checkpoint_sha256": ANCHOR_CHECKPOINT_SHA256,
            "journal_sequence": ANCHOR_JOURNAL_SEQUENCE,
            "journal_sha256": ANCHOR_JOURNAL_SHA256,
            "disposition": "nonterminal",
            "event": EVENT,
            "target": dict(TARGET),
        },
        "internal_hlt16_identity": {
            "campaign_id": INTERNAL_CAMPAIGN_ID,
            "authorization_commit": INTERNAL_AUTHORIZATION_COMMIT,
            "plan_sha256": INTERNAL_PLAN_SHA256,
            "byte_identical_to_anchor": True,
            "new_internal_campaign_created": False,
        },
        "tree_projection": {
            "source_leaf_count": SOURCE_LEAF_COUNT,
            "source_tree_sha256": source.source_tree_sha256,
            "projected_leaf_count": PROJECTED_LEAF_COUNT,
            "projected_tree_sha256": source.projected_tree_sha256,
            "excluded_paths": list(EXCLUDED_SOURCE_PATHS),
        },
        "claims": {
            "source_store_mutated": False,
            "accepted_physical_state_advanced": False,
            "common_event_completed": False,
            "GR0_calibration_completed": False,
            "candidate_execution_authorized": False,
            "SGBL_or_FGCQR_opened": False,
            "mechanism_or_physical_result_earned": False,
            "retained_EFT_evolution_authorized": False,
            "physical_transition_claim_authorized": False,
        },
    }
    return {**body, "receipt_sha256": sha256(canonical(body)).hexdigest()}


def _write_noreplace(root: Path, relative: str, payload: bytes, hook: FaultHook | None) -> None:
    parts = _safe_relative(relative, "RCV3 staged leaf")
    parent_relative = "/".join(parts[:-1])
    parent = _open_relative_directory(root, parent_relative) if parent_relative else os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    descriptor = -1
    try:
        if hook is not None:
            hook(f"before_leaf:{relative}")
        descriptor = os.open(
            parts[-1],
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o600,
            dir_fd=parent,
        )
        view = memoryview(payload)
        offset = 0
        while offset < len(view):
            written = os.write(descriptor, view[offset:])
            if written <= 0:
                raise TDG8RCV3ForkRuntimeError("RCV3 staged write made no progress")
            offset += written
        os.fsync(descriptor)
        observed = os.fstat(descriptor)
        if observed.st_size != len(payload) or not stat.S_ISREG(observed.st_mode):
            raise TDG8RCV3ForkRuntimeError("RCV3 staged leaf differs after fsync")
        os.close(descriptor)
        descriptor = -1
        os.fsync(parent)
        if hook is not None:
            hook(f"after_leaf:{relative}")
    except OSError as exc:
        raise TDG8RCV3ForkRuntimeError(f"RCV3 staged leaf cannot be published: {relative}") from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)
        os.close(parent)


def _mkdir_noreplace(root: Path, relative: str) -> None:
    parts = _safe_relative(relative, "RCV3 staged directory")
    parent_relative = "/".join(parts[:-1])
    parent = _open_relative_directory(root, parent_relative) if parent_relative else os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        os.mkdir(parts[-1], mode=0o700, dir_fd=parent)
        os.fsync(parent)
    except OSError as exc:
        raise TDG8RCV3ForkRuntimeError(
            f"RCV3 staged directory cannot be created: {relative}"
        ) from exc
    finally:
        os.close(parent)


def _wrapper_expectations(
    source: AuthenticatedRecoveryForkSource,
    authority: RecoveryForkAuthority,
) -> tuple[dict[str, bytes], tuple[str, ...], dict[str, object]]:
    receipt = _receipt(source, authority)
    leaves = {
        f"calibration/{leaf.path}": source.raw_by_path[leaf.path]
        for leaf in source.projected_leaves
    }
    leaves["bootstrap-receipt.json"] = canonical(receipt)
    directories = ("calibration",) + tuple(
        f"calibration/{name}" for name in STORE_DIRECTORIES
    )
    return leaves, directories, receipt


def _validate_wrapper(
    repository: Path,
    source: AuthenticatedRecoveryForkSource,
    authority: RecoveryForkAuthority,
    *,
    wrapper_relative: str | None = None,
    projected_store_relative: str | None = None,
) -> dict[str, object]:
    expected_leaves, expected_directories, receipt = _wrapper_expectations(
        source, authority
    )
    wrapper = source.spec.wrapper_relative if wrapper_relative is None else wrapper_relative
    projected_store = (
        source.spec.projected_store_relative
        if projected_store_relative is None
        else projected_store_relative
    )
    observed = _snapshot_tree(repository, wrapper)
    if observed.directories != expected_directories:
        raise TDG8RCV3ForkRuntimeError("RCV3 installed directory inventory differs")
    if set(observed.raw_by_path) != set(expected_leaves):
        raise TDG8RCV3ForkRuntimeError("RCV3 installed leaf inventory differs")
    for path, raw in expected_leaves.items():
        if observed.raw_by_path[path] != raw:
            raise TDG8RCV3ForkRuntimeError(f"RCV3 installed bytes differ: {path}")
    store = HLT16CampaignStore(repository / projected_store)
    try:
        anchor = store.authenticated_checkpoint_at_generation(ANCHOR_GENERATION)
        status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise TDG8RCV3ForkRuntimeError("RCV3 projected HLT16 store differs") from exc
    _source_identity(anchor, terminal=False)
    if (
        status.checkpoint_generation != ANCHOR_GENERATION
        or status.checkpoint_sha256 != ANCHOR_CHECKPOINT_SHA256
        or status.last_complete_journal_sequence != ANCHOR_JOURNAL_SEQUENCE
        or status.journal_tip_sha256 != ANCHOR_JOURNAL_SHA256
        or status.state != "clean_checkpoint"
        or status.safe_to_restart is not True
        or status.terminal is not False
        or status.active_write is not False
        or status.writer_state is not None
    ):
        raise TDG8RCV3ForkRuntimeError("RCV3 projected HLT16 restart boundary differs")
    return receipt


def _same_source(
    before: AuthenticatedRecoveryForkSource,
    after: AuthenticatedRecoveryForkSource,
) -> None:
    if (
        before.source_leaves != after.source_leaves
        or before.raw_by_path != after.raw_by_path
        or before.terminal_checkpoint != after.terminal_checkpoint
        or before.anchor_checkpoint != after.anchor_checkpoint
    ):
        raise TDG8RCV3ForkRuntimeError("RCV3 source bytes changed during projection")


def _installation(
    source: AuthenticatedRecoveryForkSource,
    authority: RecoveryForkAuthority,
    *,
    installed_now: bool,
) -> RecoveryForkInstallation:
    receipt = _receipt(source, authority)
    return RecoveryForkInstallation(
        projection_id=source.spec.projection_id,
        wrapper_relative=source.spec.wrapper_relative,
        projected_store_relative=source.spec.projected_store_relative,
        receipt_relative=source.spec.receipt_relative,
        receipt_sha256=str(receipt["receipt_sha256"]),
        anchor_checkpoint_sha256=ANCHOR_CHECKPOINT_SHA256,
        installed_now=installed_now,
    )


def _cleanup_private_stage(stage: Path) -> None:
    try:
        info = stage.lstat()
    except FileNotFoundError:
        return
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise TDG8RCV3ForkRuntimeError("RCV3 private stage became unsafe")
    shutil.rmtree(stage)


def install_recovery_fork(
    root: Path,
    spec: RecoveryForkSpec = PRODUCTION_SPEC,
    *,
    authority: RecoveryForkAuthority,
    fault_hook: FaultHook | None = None,
) -> RecoveryForkInstallation:
    """Install or exactly re-open the external RCV3 recovery projection.

    The wrapper is staged on the destination filesystem and adopted with
    Darwin ``RENAME_EXCL``.  An existing destination is accepted only when its
    complete bytes, receipt, and HLT16 generation-nine state are exact.
    """
    repository = _repository(root)
    spec.validate()
    if not isinstance(authority, RecoveryForkAuthority):
        raise TDG8RCV3ForkRuntimeError("RCV3 install authority type differs")
    authority.validate()
    source = authenticate_recovery_fork_source(repository, spec)
    if fault_hook is not None:
        fault_hook("after_source_authentication")

    try:
        target = hlt15._safe_parent_chain(repository, spec.wrapper_relative)
    except hlt15.Proto17HLT15Error as exc:
        raise TDG8RCV3ForkRuntimeError("RCV3 wrapper parent is unsafe") from exc
    if target != repository / spec.wrapper_relative:
        raise TDG8RCV3ForkRuntimeError("RCV3 wrapper resolution differs")

    try:
        target.lstat()
    except FileNotFoundError:
        pass
    except OSError as exc:
        raise TDG8RCV3ForkRuntimeError("RCV3 destination cannot be inspected") from exc
    else:
        _validate_wrapper(repository, source, authority)
        after = authenticate_recovery_fork_source(repository, spec)
        _same_source(source, after)
        return _installation(source, authority, installed_now=False)

    stage = Path(
        tempfile.mkdtemp(
            prefix=f".tdg8-rcv3-{sha256(canonical(source.projected_tree_sha256)).hexdigest()[:12]}-",
            dir=target.parent,
        )
    )
    adopted = False
    try:
        _mkdir_noreplace(stage, "calibration")
        for name in STORE_DIRECTORIES:
            _mkdir_noreplace(stage, f"calibration/{name}")
        if fault_hook is not None:
            fault_hook("after_stage_directories")

        for leaf in source.projected_leaves:
            _write_noreplace(
                stage,
                f"calibration/{leaf.path}",
                source.raw_by_path[leaf.path],
                fault_hook,
            )
        if fault_hook is not None:
            fault_hook("after_projected_leaves")
        _write_noreplace(
            stage,
            "bootstrap-receipt.json",
            canonical(_receipt(source, authority)),
            fault_hook,
        )
        if fault_hook is not None:
            fault_hook("after_receipt")
        stage_fd = os.open(
            stage,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            os.fsync(stage_fd)
        finally:
            os.close(stage_fd)
        if fault_hook is not None:
            fault_hook("after_stage_fsync")

        # Validate the stage through the unchanged HLT16 reader, independently
        # of the copy loop and external receipt comparison.
        staged_repository = stage.parent
        _validate_wrapper(
            staged_repository,
            source,
            authority,
            wrapper_relative=stage.name,
            projected_store_relative=f"{stage.name}/calibration",
        )
        if fault_hook is not None:
            fault_hook("after_staged_validation")
        current = authenticate_recovery_fork_source(repository, spec)
        _same_source(source, current)
        if fault_hook is not None:
            fault_hook("after_source_revalidation")
            fault_hook("before_exclusive_rename")
        try:
            hlt15._rename_exclusive(stage, target)
        except hlt15.Proto17HLT15Error as exc:
            # A concurrent exact installer wins harmlessly.  Any other object
            # remains untouched and is rejected as foreign.
            try:
                _validate_wrapper(repository, source, authority)
            except TDG8RCV3ForkRuntimeError:
                raise TDG8RCV3ForkRuntimeError(
                    "RCV3 destination appeared or exclusive install failed"
                ) from exc
            current = authenticate_recovery_fork_source(repository, spec)
            _same_source(source, current)
            return _installation(source, authority, installed_now=False)
        adopted = True
        if fault_hook is not None:
            fault_hook("after_exclusive_rename")
        parent = os.open(
            target.parent,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            os.fsync(parent)
        finally:
            os.close(parent)
        if fault_hook is not None:
            fault_hook("after_parent_fsync")
        _validate_wrapper(repository, source, authority)
        after = authenticate_recovery_fork_source(repository, spec)
        _same_source(source, after)
        if fault_hook is not None:
            fault_hook("after_final_validation")
        return _installation(source, authority, installed_now=True)
    finally:
        if not adopted:
            _cleanup_private_stage(stage)


__all__ = [
    "ANCHOR_CHECKPOINT_SHA256",
    "ANCHOR_GENERATION",
    "AuthenticatedRecoveryForkSource",
    "EXCLUDED_SOURCE_PATHS",
    "PROJECTED_LEAF_COUNT",
    "PROJECTED_STORE_RELATIVE",
    "PRODUCTION_SPEC",
    "PROJECTION_ID",
    "RECEIPT_RELATIVE",
    "RecoveryForkInstallation",
    "RecoveryForkAuthority",
    "RecoveryForkSpec",
    "SOURCE_LEAF_COUNT",
    "SOURCE_RELATIVE",
    "TDG8RCV3ForkRuntimeError",
    "TERMINAL_CHECKPOINT_SHA256",
    "TERMINAL_GENERATION",
    "WRAPPER_RELATIVE",
    "authenticate_recovery_fork_source",
    "install_recovery_fork",
]
