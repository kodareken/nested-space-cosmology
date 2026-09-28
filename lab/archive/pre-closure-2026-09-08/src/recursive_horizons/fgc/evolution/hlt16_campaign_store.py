"""Durable HLT16 campaign adoption, journal, checkpoint, and recovery.

The sealed PROTO17 generation-zero tree is a read-only predecessor. HLT16
adds only new sibling objects, so its eight legacy leaves can never be renamed,
overwritten, or mistaken for restartable array payloads.
"""
from __future__ import annotations

from dataclasses import dataclass
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import secrets
import stat
from typing import Any, Callable, Iterable, Mapping

from .hlt16_campaign_schema import (
    CampaignCheckpoint, HLT16CampaignSchemaError, MemberCampaignState,
    ROOT_SHA256, journal_envelope, validate_checkpoint_successor,
    validate_generation_one_bridge, validate_journal,
)
from .hlt16_member_codec import HLT16MemberSnapshot, validate_codec_metadata
from .hlt16_state_store import (
    AuthenticatedGen0Import, HLT16StateStore, HLT16StateStoreError,
    SEALED_GEN0_CHECKPOINT_SHA256, SEALED_GEN0_RECEIPT_SHA256, canonical,
    read_nofollow, STATE_SCHEMA, atomic_publish_noreplace, decode_payload,
    digest, encode_payload,
    _read_leaf_fd,
)
from .proto17_pure_construction import MEMBER_KEYS
from .proto17_pure_construction import validate_checkpoint as validate_proto17_checkpoint
from . import hlt16_lifecycle as lifecycle
from . import proto15_runtime as p15
from .numerical_engine import array_content_sha256

FaultHook = Callable[[str], None]
_RUNTIME_DIRS = ("payloads", "states", "journal", "checkpoints", "locks")


class HLT16CampaignStoreError(ValueError):
    """A campaign store is malformed, unsafe, incomplete, or forked."""


@dataclass(frozen=True, slots=True)
class WriterCapability:
    """Unforgeable-in-process handle for one active HLT16 writer lease.

    It is issued only after the exact active-lock bytes have been installed.
    Every post-genesis mutator revalidates it against the on-disk lease while
    holding ``writer.guard``; holding a stale Python object alone is never
    authority to append campaign evidence.
    """

    authorization_commit: str
    plan_sha256: str
    owner_token: str
    host: str
    pid: int


@dataclass(frozen=True, slots=True)
class _ValidatedLinkedStage:
    """An after-link stage proven to be the exact final object's twin.

    This is deliberately private and short lived.  It is not a relaxation of
    ``read_nofollow``: only a read-only store inspection may use it, and only
    for the exact final name whose task-owned stage has already been proved to
    be the same two-link inode.
    """

    stage_relative: str
    target_relative: str
    expected_sha256: str
    link_count: int
    identity: tuple[int, int, int, int, int]


def _stage_target(name: str) -> tuple[str, str] | None:
    marker = ".hlt16-stage-"
    if not name.startswith(".") or marker not in name:
        return None
    final, digest_value = name[1:].rsplit(marker, 1)
    if not final or len(digest_value) != 64 or digest_value.lower() != digest_value:
        return None
    try:
        int(digest_value, 16)
    except ValueError:
        return None
    return final, digest_value


@dataclass(frozen=True, slots=True)
class RecoveryStatus:
    authorization_commit: str
    plan_sha256: str
    checkpoint_generation: int
    checkpoint_sha256: str
    member_times: Mapping[str, str]
    member_modes: Mapping[str, str]
    active_target: Mapping[str, str]
    last_complete_journal_sequence: int
    journal_tip_sha256: str
    state: str
    safe_to_restart: bool
    terminal: bool
    active_write: bool
    writer_state: str | None


@dataclass(frozen=True, slots=True)
class AuthenticatedCampaignSnapshot:
    """Read-only facts from one completely validated campaign-store view.

    This deliberately says nothing about which suffix should be replayed,
    committed, or discarded.  A recovery coordinator may consume these facts,
    but it must make and durably bind its own transition decision.
    """

    checkpoint: CampaignCheckpoint
    suffix_records: tuple[Mapping[str, Any], ...]
    suffix_kinds: tuple[str, ...]
    suffix_classification: str
    orphan_descriptor_sha256: str | None
    orphan_payload_semantic_sha256: str | None
    staging_paths: tuple[str, ...]
    terminal_lock_present: bool


@dataclass(frozen=True, slots=True)
class _ExpectedBootstrapBridge:
    """One entirely precomputed GEN0-to-generation-one publication sequence."""

    publications: tuple[tuple[str, bytes], ...]
    members: Mapping[str, MemberCampaignState]
    record: Mapping[str, Any]
    checkpoint: CampaignCheckpoint


def _decode(raw: bytes, label: str) -> dict[str, Any]:
    def pairs(items: Iterable[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in items:
            if key in result:
                raise ValueError(key)
            result[key] = value
        return result
    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=pairs)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise HLT16CampaignStoreError(f"{label} is malformed") from exc
    if not isinstance(value, dict) or canonical(value) != raw:
        raise HLT16CampaignStoreError(f"{label} is noncanonical")
    return value


def _walk_tree(root: Path) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """FD-relative lexical traversal that never follows a replaced directory."""
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        root_fd = os.open(root, flags)
    except OSError as exc:
        raise HLT16CampaignStoreError("store root absent or unsafe") from exc
    leaves: list[str] = []
    directories: list[str] = []
    hardlinks: dict[tuple[int, int], list[str]] = {}

    def descend(directory_fd: int, prefix: str) -> None:
        try:
            with os.scandir(directory_fd) as scan:
                entries = sorted(scan, key=lambda item: item.name)
        except OSError as exc:
            raise HLT16CampaignStoreError("store directory cannot be scanned safely") from exc
        for item in entries:
            relative = f"{prefix}/{item.name}" if prefix else item.name
            try:
                before = os.stat(item.name, dir_fd=directory_fd, follow_symlinks=False)
            except OSError as exc:
                raise HLT16CampaignStoreError(f"node raced: {relative}") from exc
            if stat.S_ISLNK(before.st_mode):
                raise HLT16CampaignStoreError(f"unsafe node: {relative}")
            if stat.S_ISREG(before.st_mode):
                if before.st_nlink not in {1, 2}:
                    raise HLT16CampaignStoreError(f"hard-linked leaf: {relative}")
                leaves.append(relative)
                if before.st_nlink == 2:
                    hardlinks.setdefault((before.st_dev, before.st_ino), []).append(relative)
                continue
            if not stat.S_ISDIR(before.st_mode):
                raise HLT16CampaignStoreError(f"unsafe node: {relative}")
            try:
                child = os.open(item.name, flags, dir_fd=directory_fd)
            except OSError as exc:
                raise HLT16CampaignStoreError(f"directory raced: {relative}") from exc
            try:
                after = os.fstat(child)
                if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                    raise HLT16CampaignStoreError(f"directory raced: {relative}")
                directories.append(relative)
                descend(child, relative)
            finally:
                os.close(child)

    try:
        descend(root_fd, "")
    finally:
        os.close(root_fd)
    for paths in hardlinks.values():
        if len(paths) != 2:
            raise HLT16CampaignStoreError("hard-link inventory differs")
        staged = [path for path in paths if _stage_target(Path(path).name) is not None]
        if len(staged) != 1:
            raise HLT16CampaignStoreError("hard link is not a publication stage")
        target_name, _digest_value = _stage_target(Path(staged[0]).name)  # type: ignore[misc]
        target = f"{Path(staged[0]).parent.as_posix()}/{target_name}"
        if target not in paths:
            raise HLT16CampaignStoreError("publication stage hard link target differs")
    return tuple(sorted(leaves)), tuple(sorted(directories))


def _walk(root: Path) -> tuple[str, ...]:
    return _walk_tree(root)[0]


def _safe_read(root: Path, relative: str, label: str) -> bytes:
    try:
        return read_nofollow(root, relative, label)
    except HLT16StateStoreError as exc:
        raise HLT16CampaignStoreError(f"{label} cannot be read safely") from exc


def _read_validated_linked_final(
    root: Path,
    relative: str,
    label: str,
    stage: _ValidatedLinkedStage,
) -> bytes:
    """Read one already-validated final inode without following a stage.

    ``read_nofollow`` remains the default one-link reader.  This exception is
    available only while a read-only snapshot is observing the short window
    after ``link(stage, final)`` and before stage cleanup.  It repeats the
    inode-pair proof immediately before and after the final read so a renamed,
    replaced, or unrelated hard link cannot inherit that authorization.
    """
    if relative != stage.target_relative:
        raise HLT16CampaignStoreError("linked-stage reader target differs")
    target_parts = Path(relative).parts
    stage_parts = Path(stage.stage_relative).parts
    if (
        len(target_parts) != 2
        or len(stage_parts) != 2
        or target_parts[0] != stage_parts[0]
    ):
        raise HLT16CampaignStoreError("linked-stage reader path differs")
    directory, final_name = target_parts
    _stage_directory, stage_name = stage_parts
    parsed = _stage_target(stage_name)
    if parsed != (final_name, stage.expected_sha256):
        raise HLT16CampaignStoreError("linked-stage reader name differs")
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    root_fd = child = descriptor = -1
    try:
        root_fd = os.open(root, flags)
        child = os.open(directory, flags, dir_fd=root_fd)

        def pair() -> tuple[os.stat_result, os.stat_result]:
            final_info = os.stat(final_name, dir_fd=child, follow_symlinks=False)
            stage_info = os.stat(stage_name, dir_fd=child, follow_symlinks=False)
            for info in (final_info, stage_info):
                if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode) or info.st_nlink != 2:
                    raise HLT16CampaignStoreError("linked-stage reader leaf is unsafe")
            observed = (
                final_info.st_dev,
                final_info.st_ino,
                final_info.st_size,
                final_info.st_mtime_ns,
                final_info.st_ctime_ns,
            )
            stage_identity = (
                stage_info.st_dev,
                stage_info.st_ino,
                stage_info.st_size,
                stage_info.st_mtime_ns,
                stage_info.st_ctime_ns,
            )
            if observed != stage.identity or stage_identity != stage.identity:
                raise HLT16CampaignStoreError("linked-stage reader inode differs")
            return final_info, stage_info

        before, _stage_before = pair()
        descriptor = os.open(final_name, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0), dir_fd=child)
        active = os.fstat(descriptor)
        active_identity = (
            active.st_dev, active.st_ino, active.st_size,
            active.st_mtime_ns, active.st_ctime_ns,
        )
        if active_identity != stage.identity:
            raise HLT16CampaignStoreError("linked-stage reader raced")
        chunks: list[bytes] = []
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            chunks.append(block)
        raw = b"".join(chunks)
        after = os.fstat(descriptor)
        after_identity = (
            after.st_dev, after.st_ino, after.st_size,
            after.st_mtime_ns, after.st_ctime_ns,
        )
        if len(raw) != before.st_size or after_identity != stage.identity:
            raise HLT16CampaignStoreError("linked-stage reader changed during read")
        pair()
        if hashlib.sha256(raw).hexdigest() != stage.expected_sha256:
            raise HLT16CampaignStoreError("linked-stage reader hash differs")
        return raw
    except OSError as exc:
        raise HLT16CampaignStoreError("linked-stage reader cannot inspect safely") from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)
        if child != -1:
            os.close(child)
        if root_fd != -1:
            os.close(root_fd)


def _safe_read_snapshot(
    root: Path,
    relative: str,
    label: str,
    stages_by_target: Mapping[str, _ValidatedLinkedStage],
) -> bytes:
    """Use the narrowly validated two-link reader only for its exact target."""
    stage = stages_by_target.get(relative)
    if stage is None:
        return _safe_read(root, relative, label)
    return _read_validated_linked_final(root, relative, label, stage)


def _valid_runtime_leaf(directory: str, name: str) -> bool:
    if directory == "payloads":
        return len(name) == 68 and name.endswith(".npz")
    if directory == "journal":
        return len(name) == 93 and name.endswith(".journal")
    if directory == "states":
        return len(name) == 69 and name.endswith(".json")
    if directory == "checkpoints":
        return len(name) == 90 and name.endswith(".json")
    if directory == "locks":
        return name in {
            "active-write.lock", "terminal.lock", "writer.guard",
            # This is an advisory flock inode, not a scientific record.  It
            # serializes the one irreversible GEN0 -> HLT16 adoption before a
            # generation-one checkpoint exists to which an active writer can
            # be bound.
            "bootstrap.guard",
        } or (
            name.startswith(".active-write.lock.hlt16-quarantine-")
            and len(name) == len(".active-write.lock.hlt16-quarantine-") + 64
            and all(char in "0123456789abcdef" for char in name.rsplit("-", 1)[1])
        )
    return directory == "receipts" and name == "generation-zero.json"


def _retired_writer_digest(name: str) -> str | None:
    prefix = ".active-write.lock.hlt16-quarantine-"
    if not name.startswith(prefix) or len(name) != len(prefix) + 64:
        return None
    digest_value = name[len(prefix):]
    try:
        int(digest_value, 16)
    except ValueError:
        return None
    return digest_value if digest_value == digest_value.lower() else None


def _validated_stage(root: Path, relative: str) -> _ValidatedLinkedStage:
    """Validate one exact task-stage and return its intended final path/hash."""
    parts = Path(relative).parts
    if len(parts) != 2:
        raise HLT16CampaignStoreError("nested publication stage differs")
    directory, stage_name = parts
    parsed = _stage_target(stage_name)
    if parsed is None:
        raise HLT16CampaignStoreError("publication stage name differs")
    final, expected_hash = parsed
    if not _valid_runtime_leaf(directory, final) or directory in {"receipts"}:
        raise HLT16CampaignStoreError("publication stage target differs")
    directory_flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    parent = os.open(root, directory_flags)
    child = -1
    descriptor = -1
    try:
        child = os.open(directory, directory_flags, dir_fd=parent)
        before = os.stat(stage_name, dir_fd=child, follow_symlinks=False)
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_nlink not in {1, 2}:
            raise HLT16CampaignStoreError("publication stage is unsafe")
        descriptor = os.open(
            stage_name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=child,
        )
        active = os.fstat(descriptor)
        if (before.st_dev, before.st_ino, before.st_size) != (active.st_dev, active.st_ino, active.st_size):
            raise HLT16CampaignStoreError("publication stage raced")
        digestor = hashlib.sha256()
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            digestor.update(block)
        observed_hash = digestor.hexdigest()
        try:
            final_info = os.stat(final, dir_fd=child, follow_symlinks=False)
        except FileNotFoundError:
            final_info = None
        if before.st_nlink == 2:
            if (
                final_info is None
                or not stat.S_ISREG(final_info.st_mode)
                or (before.st_dev, before.st_ino) != (final_info.st_dev, final_info.st_ino)
                or observed_hash != expected_hash
            ):
                raise HLT16CampaignStoreError("linked publication stage differs")
        elif final_info is not None:
            raise HLT16CampaignStoreError("unlinked publication stage forks a final object")
        # A one-link stage may be a complete pre-link object or an interrupted
        # partial write.  Both are a recognized, nonrestartable suffix; the
        # exact publication call owns repair or rollback.
        identity = (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
        )
        return _ValidatedLinkedStage(
            stage_relative=relative,
            target_relative=f"{directory}/{final}",
            expected_sha256=expected_hash,
            link_count=before.st_nlink,
            identity=identity,
        )
    except OSError as exc:
        raise HLT16CampaignStoreError("publication stage cannot be inspected safely") from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)
        if child != -1:
            os.close(child)
        os.close(parent)


def _validate_stage(root: Path, relative: str) -> tuple[str, str]:
    """Compatibility view of one validated publication stage."""
    stage = _validated_stage(root, relative)
    return stage.target_relative, stage.expected_sha256


def _descriptor_only(
    root: Path,
    address: str,
    *,
    snapshot_stages: Mapping[str, _ValidatedLinkedStage] | None = None,
) -> dict[str, Any]:
    relative = f"states/{address}.json"
    if snapshot_stages is None:
        raw = _safe_read(root, relative, "state descriptor")
    else:
        raw = _safe_read_snapshot(root, relative, "state descriptor", snapshot_stages)
    value = _decode(raw, "state descriptor")
    observed = value.get("descriptor_sha256")
    bare = dict(value)
    bare.pop("descriptor_sha256", None)
    if (
        value.get("schema") != STATE_SCHEMA
        or observed != address
        or hashlib.sha256(canonical(bare)).hexdigest() != address
    ):
        raise HLT16CampaignStoreError("state descriptor address differs")
    try:
        metadata = validate_codec_metadata(value["metadata"])
    except Exception as exc:
        raise HLT16CampaignStoreError("state descriptor metadata differs") from exc
    return {**value, "metadata": metadata}


class HLT16CampaignStore:
    """In-place HLT16 extension bound to an immutable PROTO17 GEN0 root."""

    def __init__(self, legacy_root: Path, *, fault_hook: FaultHook | None = None) -> None:
        self.legacy_root = Path(legacy_root)
        self.root = self.legacy_root
        self.fault_hook = fault_hook
        self._states = HLT16StateStore(self.root, fault_hook=self._state_fault)
        self._writer_capability: WriterCapability | None = None

    def _cut(self, phase: str) -> None:
        if self.fault_hook is not None:
            self.fault_hook(phase)

    def _state_fault(self, phase: str) -> None:
        self._cut(f"state:{phase}")

    def _repair_linked_stages(self) -> None:
        """Remove only a proven two-link stage whose final inode is present."""
        leaves, _directories = _walk_tree(self.root)
        candidates = [item for item in leaves if _stage_target(Path(item).name) is not None]
        for relative in candidates:
            _target, _observed = _validate_stage(self.root, relative)
            directory, name = Path(relative).parts
            flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
            root_fd = os.open(self.root, flags)
            child = -1
            try:
                child = os.open(directory, flags, dir_fd=root_fd)
                info = os.stat(name, dir_fd=child, follow_symlinks=False)
                if info.st_nlink == 2:
                    try:
                        os.unlink(name, dir_fd=child)
                    except FileNotFoundError:
                        pass
                    os.fsync(child)
            finally:
                if child != -1:
                    os.close(child)
                os.close(root_fd)

    def _runtime_inventory(self, *, repair_stages: bool = True) -> tuple[str, ...]:
        if repair_stages:
            self._repair_linked_stages()
        observed = _walk(self.root)
        return observed

    def _legacy_inventory(
        self,
        authorities: Mapping[str, AuthenticatedGen0Import] | None = None,
        *,
        repair_stages: bool = True,
        snapshot_stages: Mapping[str, _ValidatedLinkedStage] | None = None,
    ) -> dict[str, str]:
        """Validate the exact eight-leaf source tree HLT16 is permitted to bridge."""
        if repair_stages:
            self._repair_linked_stages()
        checkpoint_rel = f"checkpoints/00000000000000000000-{SEALED_GEN0_CHECKPOINT_SHA256}.json"
        receipt_rel = "receipts/generation-zero.json"
        read = (
            _safe_read
            if snapshot_stages is None
            else lambda root, relative, label: _safe_read_snapshot(
                root, relative, label, snapshot_stages
            )
        )
        receipt = _decode(read(self.legacy_root, receipt_rel, "GEN0 receipt"), "GEN0 receipt")
        checkpoint = _decode(read(self.legacy_root, checkpoint_rel, "GEN0 checkpoint"), "GEN0 checkpoint")
        try:
            checked = validate_proto17_checkpoint(checkpoint)
        except Exception as exc:
            raise HLT16CampaignStoreError("sealed GEN0 checkpoint schema differs") from exc
        if (checked.get("checkpoint_sha256") != SEALED_GEN0_CHECKPOINT_SHA256
                or receipt.get("receipt_sha256") != SEALED_GEN0_RECEIPT_SHA256
                or receipt.get("checkpoint_sha256") != SEALED_GEN0_CHECKPOINT_SHA256
                or receipt.get("state_object_set_sha256") != checked.get("accepted_state_set_sha256")):
            raise HLT16CampaignStoreError("sealed GEN0 receipt/checkpoint binding differs")
        states = checked.get("states")
        cursors = checked.get("cursors")
        if not isinstance(states, Mapping) or not isinstance(cursors, Mapping) or set(states) != set(MEMBER_KEYS):
            raise HLT16CampaignStoreError("sealed GEN0 member inventory differs")
        source = {key: str(states[key]) for key in MEMBER_KEYS}
        if len(set(source.values())) != len(source):
            raise HLT16CampaignStoreError("GEN0 source state identities are not distinct")
        expected = {*(f"states/{value}.json" for value in source.values()), checkpoint_rel, receipt_rel}
        leaves, directories = _walk_tree(self.legacy_root)
        allowed_directories = {"checkpoints", "receipts", "states", "payloads", "journal", "locks"}
        if (not {"checkpoints", "receipts", "states"}.issubset(directories)
                or any("/" in item or item not in allowed_directories for item in directories)):
            raise HLT16CampaignStoreError("legacy GEN0 directory inventory differs")
        observed = set(leaves)
        allowed_new = {item for item in observed if item.startswith(("payloads/", "journal/", "locks/"))
                       or (item.startswith("checkpoints/") and not Path(item).name.startswith("00000000000000000000-"))
                       or (item.startswith("states/") and Path(item).name not in {f"{value}.json" for value in source.values()})}
        if observed - allowed_new != expected:
            raise HLT16CampaignStoreError("legacy GEN0 eight-leaf inventory differs")
        for key, digest in source.items():
            payload = _decode(
                read(self.legacy_root, f"states/{digest}.json", f"GEN0 state {key}"),
                f"GEN0 state {key}",
            )
            if hashlib.sha256(canonical(payload)).hexdigest() != digest:
                raise HLT16CampaignStoreError("GEN0 source state identity differs")
        if authorities is not None:
            if set(authorities) != set(MEMBER_KEYS):
                raise HLT16CampaignStoreError("GEN0 authority member inventory differs")
            for key in MEMBER_KEYS:
                authority = authorities[key]
                if not isinstance(authority, AuthenticatedGen0Import):
                    raise HLT16CampaignStoreError("GEN0 authority type differs")
                authority.validate()
                if authority.source_state_object_sha256 != source[key]:
                    raise HLT16CampaignStoreError("GEN0 authority/state relation differs")
        return source

    def _make_runtime_dirs_exclusive(self) -> None:
        flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
        try:
            root_fd = os.open(self.root, flags)
        except OSError as exc:
            raise HLT16CampaignStoreError("legacy root unsafe") from exc
        try:
            for name in ("payloads", "journal", "locks", "checkpoints"):
                try:
                    os.mkdir(name, mode=0o755, dir_fd=root_fd)
                    os.fsync(root_fd)
                except FileExistsError:
                    before = os.stat(name, dir_fd=root_fd, follow_symlinks=False)
                    if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                        raise HLT16CampaignStoreError("HLT16 runtime directory unsafe")
                    child = os.open(name, flags, dir_fd=root_fd)
                    try:
                        after = os.fstat(child)
                        if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                            raise HLT16CampaignStoreError("HLT16 runtime directory raced")
                    finally:
                        os.close(child)
        finally:
            os.close(root_fd)

    def _load_state_snapshot(
        self,
        address: str,
        *,
        snapshot_stages: Mapping[str, _ValidatedLinkedStage] | None,
    ) -> Mapping[str, Any]:
        """Validate one state for a read-only snapshot without changing loads.

        The ordinary state-store API intentionally insists on single-link
        leaves.  A snapshot may instead consume the already proved stage map
        to inspect a final descriptor or payload in the tiny after-link
        window.  This is validation-only; it never becomes a restart API.
        """
        if snapshot_stages is None:
            return self._states.load_state(address).arrays
        descriptor = _descriptor_only(
            self.root, address, snapshot_stages=snapshot_stages
        )
        semantic = descriptor["semantic_sha256"]
        raw = _safe_read_snapshot(
            self.root,
            f"payloads/{semantic}.npz",
            "state payload",
            snapshot_stages,
        )
        if hashlib.sha256(raw).hexdigest() != descriptor["raw_archive_sha256"]:
            raise HLT16CampaignStoreError("state payload raw hash differs")
        try:
            arrays, observed, _manifest = decode_payload(
                raw,
                expected_semantic_sha256=semantic,
                expected_manifest=descriptor["arrays"],
            )
        except HLT16StateStoreError as exc:
            raise HLT16CampaignStoreError("state payload differs") from exc
        if observed != semantic:
            raise HLT16CampaignStoreError("state payload semantic identity differs")
        return arrays

    def _publish(self, relative: str, payload: bytes, *, phase: str) -> None:
        try:
            atomic_publish_noreplace(
                self.root,
                relative,
                payload,
                fault_hook=lambda cut: self._cut(f"{phase}:{cut}"),
            )
        except HLT16StateStoreError as exc:
            raise HLT16CampaignStoreError("runtime publication failed safely") from exc

    def _journal(
        self,
        *,
        repair_stages: bool = True,
        snapshot_stages: Mapping[str, _ValidatedLinkedStage] | None = None,
    ) -> tuple[dict[str, Any], ...]:
        records: list[dict[str, Any]] = []
        for relative in sorted(
            item for item in self._runtime_inventory(repair_stages=repair_stages)
            if item.startswith("journal/") and _stage_target(Path(item).name) is None
        ):
            name = Path(relative).name
            if not name.endswith(".journal") or len(name) != 93:
                raise HLT16CampaignStoreError("journal filename differs")
            sequence, observed = name[:-8].split("-", 1)
            if sequence != f"{len(records):020d}":
                raise HLT16CampaignStoreError("journal sequence/fork differs")
            try:
                raw = (
                    _safe_read(self.root, relative, "journal")
                    if snapshot_stages is None
                    else _safe_read_snapshot(self.root, relative, "journal", snapshot_stages)
                )
                record = validate_journal(_decode(raw, "journal"))
            except HLT16CampaignSchemaError as exc:
                raise HLT16CampaignStoreError("journal schema differs") from exc
            parent = ROOT_SHA256 if not records else records[-1]["record_sha256"]
            if observed != record["record_sha256"] or record["previous_record_sha256"] != parent:
                raise HLT16CampaignStoreError("journal hash-chain differs")
            records.append(record)
        return tuple(records)

    def _checkpoints(
        self,
        *,
        repair_stages: bool = True,
        snapshot_stages: Mapping[str, _ValidatedLinkedStage] | None = None,
    ) -> tuple[CampaignCheckpoint, ...]:
        result: list[CampaignCheckpoint] = []
        for relative in sorted(
            item for item in self._runtime_inventory(repair_stages=repair_stages)
            if item.startswith("checkpoints/")
            and _stage_target(Path(item).name) is None
            and not Path(item).name.startswith("00000000000000000000-")
        ):
            name = Path(relative).name
            if not name.endswith(".json") or len(name) != 90:
                raise HLT16CampaignStoreError("checkpoint filename differs")
            generation, observed = name[:-5].split("-", 1)
            if generation != f"{len(result) + 1:020d}":
                raise HLT16CampaignStoreError("checkpoint generation/fork differs")
            try:
                raw = (
                    _safe_read(self.root, relative, "checkpoint")
                    if snapshot_stages is None
                    else _safe_read_snapshot(self.root, relative, "checkpoint", snapshot_stages)
                )
                checkpoint = CampaignCheckpoint.from_mapping(_decode(raw, "checkpoint"))
            except HLT16CampaignSchemaError as exc:
                raise HLT16CampaignStoreError("checkpoint schema differs") from exc
            if observed != checkpoint.sha256 or checkpoint.generation != len(result) + 1:
                raise HLT16CampaignStoreError("checkpoint identity differs")
            result.append(checkpoint)
        return tuple(result)

    def _validate_full_store(
        self,
        authorities: Mapping[str, AuthenticatedGen0Import] | None = None,
        *,
        repair_stages: bool = True,
        include_snapshot_facts: bool = False,
    ) -> (
        tuple[tuple[dict[str, Any], ...], tuple[CampaignCheckpoint, ...], str | None]
        | tuple[tuple[dict[str, Any], ...], tuple[CampaignCheckpoint, ...], str | None, AuthenticatedCampaignSnapshot]
    ):
        """Validate the complete store plus one narrowly recognized crash suffix."""
        # A read-only snapshot must be able to inspect the exact short window
        # after ``link(stage, final)``.  Prove every such pair before any
        # target is read, then pass the proof only to snapshot reads.  Normal
        # recovery/publication still uses the one-link reader after repairing
        # these task-owned stages.
        pre_inventory, _pre_directories = _walk_tree(self.root)
        pre_stage_paths = {
            item for item in pre_inventory if _stage_target(Path(item).name) is not None
        }
        snapshot_stages: dict[str, _ValidatedLinkedStage] = {}
        if not repair_stages:
            for relative in sorted(pre_stage_paths):
                stage = _validated_stage(self.root, relative)
                if stage.target_relative in snapshot_stages:
                    raise HLT16CampaignStoreError("publication stages share a final target")
                snapshot_stages[stage.target_relative] = stage
        source = self._legacy_inventory(
            authorities,
            repair_stages=repair_stages,
            snapshot_stages=snapshot_stages if not repair_stages else None,
        )
        records = self._journal(
            repair_stages=repair_stages,
            snapshot_stages=snapshot_stages if not repair_stages else None,
        )
        checkpoints = self._checkpoints(
            repair_stages=repair_stages,
            snapshot_stages=snapshot_stages if not repair_stages else None,
        )
        inventory, directories = _walk_tree(self.root)
        expected_directories = (
            "checkpoints", "journal", "locks", "payloads", "receipts", "states",
        )
        if directories != expected_directories:
            raise HLT16CampaignStoreError("store directory inventory differs")
        if not checkpoints:
            raise HLT16CampaignStoreError("HLT16 bridge checkpoint absent")

        previous: CampaignCheckpoint | None = None
        start = 0
        references: set[str] = set()
        baselines: dict[str, tuple[int, int]] = {}
        descriptor_cache: dict[str, dict[str, Any]] = {}
        evolution_state_hash_cache: dict[str, str] = {}

        def descriptor_for(address: str) -> dict[str, Any]:
            if address not in descriptor_cache:
                descriptor_cache[address] = _descriptor_only(
                    self.root,
                    address,
                    snapshot_stages=snapshot_stages if not repair_stages else None,
                )
            return descriptor_cache[address]

        def evolution_state_hash_for(address: str) -> str:
            if address not in evolution_state_hash_cache:
                arrays = self._load_state_snapshot(
                    address,
                    snapshot_stages=snapshot_stages if not repair_stages else None,
                )
                evolution_state_hash_cache[address] = array_content_sha256(
                    arrays["u"], arrays["p"], arrays["q"]
                )
            return evolution_state_hash_cache[address]

        for checkpoint in checkpoints:
            if checkpoint.journal_sequence < start or checkpoint.journal_sequence >= len(records):
                raise HLT16CampaignStoreError("checkpoint journal span differs")
            transition = records[start:checkpoint.journal_sequence + 1]
            try:
                if previous is None:
                    validate_generation_one_bridge(
                        checkpoint,
                        records[0],
                        external_checkpoint_sha256=SEALED_GEN0_CHECKPOINT_SHA256,
                    )
                else:
                    validate_checkpoint_successor(previous, checkpoint, transition)
            except HLT16CampaignSchemaError as exc:
                raise HLT16CampaignStoreError("checkpoint successor relation differs") from exc
            accepted = transition[0]["payload"] if previous is not None and len(transition) == 1 and transition[0]["kind"] == "accepted_fine" else None
            for key in MEMBER_KEYS:
                address = checkpoint.members[key].descriptor_sha256
                descriptor = descriptor_for(address)
                metadata = descriptor["metadata"]
                runtime = metadata["runtime_identity"]
                cursor = checkpoint.members[key].cursor
                ledger = checkpoint.members[key].ledger
                extension = metadata["tdg6"]["extension"]
                debit = metadata["tdg6"]["accumulated_debit_hex"]
                encoded_ledger = {
                    "last_accepted_time_hex": float(extension["last_accepted_time"]).hex(),
                    "accepted_macro_step_count": extension["accepted_macro_step_count"],
                    "cumulative_temporal_retry_count": extension["cumulative_temporal_retry_count"],
                    "current_macro_step_temporal_retry_count": extension["current_macro_step_temporal_retry_count"],
                    "last_accepted_macro_step_temporal_retry_count": extension["last_accepted_macro_step_temporal_retry_count"],
                    "accumulated_debit_vector_hex": list(debit),
                    "serialized_temporal_rejections": extension["serialized_temporal_rejections"],
                }
                if (
                    runtime["member_key"] != key
                    or runtime["campaign_id"] != checkpoint.campaign_id
                    or runtime["accepted_time"] != cursor["accepted_boundary_time"]
                    or metadata["counters"]["step_index"] != cursor["previous_step_index"]
                    or metadata["counters"]["transaction_serial"] != cursor["previous_transaction_serial"]
                    or runtime["generation"] > checkpoint.generation
                ):
                    raise HLT16CampaignStoreError("checkpoint descriptor member differs")
                if cursor["mode"] == "RETRY_PENDING":
                    try:
                        lifecycle.validate_retry_payload(
                            p15.Proto15Cursor(cursor),
                            p15._ledger_from_mapping(ledger),
                            evolution_state_sha256=evolution_state_hash_for(address),
                        )
                    except Exception as exc:
                        raise HLT16CampaignStoreError(
                            "checkpoint retry state identity differs"
                        ) from exc
                if checkpoint.generation == 1:
                    if (
                        runtime["generation"] != 0
                        or metadata["provenance_cursor"]["source_protocol"] != "FGC-2-SF1-PROTO17"
                        # A GEN0 descriptor deliberately retains the immutable
                        # PROTO17 predecessor cursor.  The bridge cursor is a
                        # distinct PROTO18 object in the genesis journal; it
                        # must point back to that descriptor provenance rather
                        # than be substituted into the historical metadata.
                        or metadata["provenance_cursor"]["payload"].get("cursor_chain_sha256")
                        != records[0]["payload"]["bridge_cursors"][key].get("cursor_chain_parent_sha256")
                        or encoded_ledger != ledger
                    ):
                        raise HLT16CampaignStoreError("generation-one descriptor provenance differs")
                    baselines[key] = (
                        metadata["counters"]["source_retry_count"],
                        metadata["counters"]["CFL_retry_count"],
                    )
                if accepted is not None and accepted["member_key"] == key:
                    before = previous.members[key]
                    if (
                        runtime["generation"] != checkpoint.generation
                        or metadata["provenance_cursor"]["source_protocol"] != "FGC-2-SF1-PROTO18"
                        or metadata["provenance_cursor"]["payload"] != before.cursor
                        or metadata["accepted_plan"] != accepted["executed_plan"]
                        or encoded_ledger != ledger
                        or metadata["counters"]["source_retry_count"] != baselines[key][0] + checkpoint.members[key].source_total
                        or metadata["counters"]["CFL_retry_count"] != baselines[key][1] + checkpoint.members[key].cfl_total
                    ):
                        raise HLT16CampaignStoreError("accepted descriptor lineage differs")
                references.add(address)
            previous = checkpoint
            start = checkpoint.journal_sequence + 1

        latest = checkpoints[-1]
        suffix = records[latest.journal_sequence + 1:]
        suffix_kinds = tuple(record["kind"] for record in suffix)
        permitted_suffixes = {
            (),
            ("accepted_fine",),
            ("source_rejection",),
            ("cfl_rejection",),
            ("tdg6_rejection",),
            ("tdg6_rejection", "cursor_transition"),
            ("terminal_lock",),
            ("source_rejection", "terminal_lock"),
            ("cfl_rejection", "terminal_lock"),
            ("tdg6_rejection", "terminal_lock"),
            ("common_event_commit",),
        }
        if suffix_kinds not in permitted_suffixes:
            raise HLT16CampaignStoreError("unrecognized journal crash suffix")
        for record in suffix:
            if (
                record["campaign_id"] != latest.campaign_id
                or record["generation"] != latest.generation + 1
                or record["event"] != latest.event
            ):
                raise HLT16CampaignStoreError("journal crash suffix authority differs")
        if suffix_kinds and suffix_kinds[0] == "tdg6_rejection":
            rejection = suffix[0]["payload"]
            key = rejection["member_key"]
            member = latest.members[key]
            physical = evolution_state_hash_for(member.descriptor_sha256)
            try:
                typed = lifecycle.typed_retry_from_evidence(
                    p15.Proto15Cursor(member.cursor),
                    p15._ledger_from_mapping(member.ledger),
                    rejection["evidence"],
                    evolution_state_sha256=physical,
                )
                if suffix_kinds == ("tdg6_rejection", "cursor_transition"):
                    lifecycle.validate_retry_payload(
                        p15.Proto15Cursor(
                            suffix[1]["payload"]["successor_cursor"]
                        ),
                        typed.updated_ledger,
                        evolution_state_sha256=physical,
                    )
            except Exception as exc:
                raise HLT16CampaignStoreError(
                    "TDG6 suffix evidence reconstruction differs"
                ) from exc

        inventory_set = set(inventory)
        stage_paths = {item for item in inventory if _stage_target(Path(item).name) is not None}
        if not repair_stages and stage_paths != pre_stage_paths:
            raise HLT16CampaignStoreError("publication stage inventory raced")
        for item in stage_paths:
            stage = _validated_stage(self.root, item)
            if not repair_stages and snapshot_stages.get(stage.target_relative) != stage:
                raise HLT16CampaignStoreError("publication stage validation raced")

        legacy_paths = {
            "receipts/generation-zero.json",
            f"checkpoints/00000000000000000000-{SEALED_GEN0_CHECKPOINT_SHA256}.json",
            *(f"states/{value}.json" for value in source.values()),
        }
        journal_paths = {
            f"journal/{record['sequence']:020d}-{record['record_sha256']}.journal"
            for record in records
        }
        checkpoint_paths = {
            f"checkpoints/{checkpoint.generation:020d}-{checkpoint.sha256}.json"
            for checkpoint in checkpoints
        }
        actual_descriptor_paths = {
            item for item in inventory
            if item.startswith("states/")
            and item not in legacy_paths
            and item not in stage_paths
        }
        actual_descriptors = {Path(item).stem for item in actual_descriptor_paths}
        extra_descriptors = actual_descriptors - references
        if len(extra_descriptors) > 1:
            raise HLT16CampaignStoreError("multiple uncheckpointed descriptors differ")

        if suffix_kinds == ("accepted_fine",):
            accepted_hash = suffix[0]["payload"]["descriptor_sha256"]
            if extra_descriptors != {accepted_hash}:
                raise HLT16CampaignStoreError("accepted journal/descriptor suffix differs")
        elif suffix_kinds and extra_descriptors:
            raise HLT16CampaignStoreError("nonaccepted suffix owns an evolved descriptor")

        if extra_descriptors:
            address = next(iter(extra_descriptors))
            descriptor = descriptor_for(address)
            metadata = descriptor["metadata"]
            key = metadata["runtime_identity"]["member_key"]
            before = latest.members[key]
            if (
                metadata["runtime_identity"]["generation"] != latest.generation + 1
                or metadata["runtime_identity"]["campaign_id"] != latest.campaign_id
                or metadata["provenance_cursor"]["source_protocol"] != "FGC-2-SF1-PROTO18"
                or metadata["provenance_cursor"]["payload"] != before.cursor
            ):
                raise HLT16CampaignStoreError("uncheckpointed descriptor lineage differs")
            if suffix_kinds == ("accepted_fine",) and metadata["accepted_plan"] != suffix[0]["payload"]["executed_plan"]:
                raise HLT16CampaignStoreError("uncheckpointed descriptor plan differs")
            self._load_state_snapshot(
                address,
                snapshot_stages=snapshot_stages if not repair_stages else None,
            )
            references.add(address)

        expected_payloads = {
            str(descriptor_for(address)["semantic_sha256"]) for address in references
        }
        actual_payload_paths = {
            item for item in inventory
            if item.startswith("payloads/") and item not in stage_paths
        }
        actual_payloads = {Path(item).stem for item in actual_payload_paths}
        orphan_payloads = actual_payloads - expected_payloads
        if len(orphan_payloads) > 1:
            raise HLT16CampaignStoreError("multiple orphan payloads differ")
        for semantic in orphan_payloads:
            relative = f"payloads/{semantic}.npz"
            raw = (
                _safe_read(self.root, relative, "orphan payload")
                if repair_stages
                else _safe_read_snapshot(self.root, relative, "orphan payload", snapshot_stages)
            )
            _arrays, observed, _manifest = decode_payload(
                raw,
                expected_semantic_sha256=semantic,
            )
            if observed != semantic:
                raise HLT16CampaignStoreError("orphan payload identity differs")

        # A restart is allowed only after all six currently selected payloads
        # have been independently decoded and hash-checked.  Historical payload
        # arrays remain content-addressed but are not re-inflated at every step.
        for key in MEMBER_KEYS:
            self._load_state_snapshot(
                latest.members[key].descriptor_sha256,
                snapshot_stages=snapshot_stages if not repair_stages else None,
            )

        lock_paths = {
            item for item in inventory
            if item.startswith("locks/") and item not in stage_paths
        }
        if any(not _valid_runtime_leaf("locks", Path(item).name) for item in lock_paths):
            raise HLT16CampaignStoreError("lock inventory differs")
        trusted_checkpoint_sha256s = {item.sha256 for item in checkpoints}
        for relative in sorted(lock_paths):
            name = Path(relative).name
            digest_value = _retired_writer_digest(name)
            if digest_value is None:
                continue
            raw = (
                _safe_read(self.root, relative, "retired writer lease")
                if repair_stages
                else _safe_read_snapshot(self.root, relative, "retired writer lease", snapshot_stages)
            )
            if hashlib.sha256(raw).hexdigest() != digest_value:
                raise HLT16CampaignStoreError("retired writer lease filename hash differs")
            value = _decode(raw, "retired writer lease")
            if (
                set(value) != {"schema", "authorization_commit", "plan_sha256", "checkpoint_sha256", "host", "pid", "owner_token"}
                or value["schema"] != "FGC-1-HLT16-active-writer-v1"
                or value["authorization_commit"] != latest.authorization_commit
                or value["plan_sha256"] != latest.plan_sha256
                or value["checkpoint_sha256"] not in trusted_checkpoint_sha256s
                or not isinstance(value["host"], str) or not value["host"]
                or isinstance(value["pid"], bool) or not isinstance(value["pid"], int) or value["pid"] <= 0
            ):
                raise HLT16CampaignStoreError("retired writer lease authority differs")
            self._owner_token(value["owner_token"])
        expected = (
            legacy_paths
            | journal_paths
            | checkpoint_paths
            | {f"states/{address}.json" for address in references}
            | {f"payloads/{semantic}.npz" for semantic in expected_payloads | orphan_payloads}
            | lock_paths
            | stage_paths
        )
        if inventory_set != expected:
            raise HLT16CampaignStoreError("store leaf inventory differs")

        suffix_state: str | None = None
        if stage_paths:
            suffix_state = "recognized_staging_suffix"
        elif suffix:
            suffix_state = "recognized_uncheckpointed_journal_suffix"
        elif extra_descriptors or orphan_payloads:
            suffix_state = "recognized_uncheckpointed_state_suffix"
        if include_snapshot_facts:
            if stage_paths:
                classification = "staging_suffix"
            elif suffix_kinds:
                classification = "+".join(suffix_kinds)
            elif extra_descriptors:
                classification = "descriptor_orphan"
            elif orphan_payloads:
                classification = "payload_orphan"
            else:
                classification = "clean"
            snapshot = AuthenticatedCampaignSnapshot(
                checkpoint=latest,
                suffix_records=tuple(dict(record) for record in suffix),
                suffix_kinds=suffix_kinds,
                suffix_classification=classification,
                orphan_descriptor_sha256=next(iter(extra_descriptors), None),
                orphan_payload_semantic_sha256=next(iter(orphan_payloads), None),
                staging_paths=tuple(sorted(stage_paths)),
                terminal_lock_present=self._terminal_lock(
                    latest,
                    snapshot_stages=snapshot_stages,
                ),
            )
            return records, checkpoints, suffix_state, snapshot
        return records, checkpoints, suffix_state

    def authenticated_snapshot(
        self,
        *,
        authorities: Mapping[str, AuthenticatedGen0Import] | None = None,
    ) -> AuthenticatedCampaignSnapshot:
        """Return validated current facts without repairing or deciding a suffix.

        In particular, a completed hard-link publication stage remains visible
        as staging here.  Only an explicit writer/recovery operation is allowed
        to remove that task-owned stage after deciding the next action.
        """
        result = self._validate_full_store(
            authorities,
            repair_stages=False,
            include_snapshot_facts=True,
        )
        return result[3]  # type: ignore[index,return-value]

    def authenticated_checkpoint_at_generation(
        self,
        generation: int,
        *,
        authorities: Mapping[str, AuthenticatedGen0Import] | None = None,
    ) -> CampaignCheckpoint:
        """Return one closed historical checkpoint after validating the full tree.

        A successor campaign occasionally needs an immutable accepted boundary
        that predates the latest terminal checkpoint.  Selecting that boundary
        must not weaken validation to a single JSON leaf, repair a publication
        suffix, or make the selected ancestor writable.  This accessor first
        authenticates the complete checkpoint/journal/state lineage with
        repairs disabled and then returns the unique requested generation.
        """
        if isinstance(generation, bool) or not isinstance(generation, int) or generation < 1:
            raise HLT16CampaignStoreError("historical checkpoint generation differs")
        _records, checkpoints, suffix = self._validate_full_store(
            authorities,
            repair_stages=False,
        )
        if suffix is not None:
            raise HLT16CampaignStoreError(
                "historical checkpoint selection requires a closed store"
            )
        selected = tuple(item for item in checkpoints if item.generation == generation)
        if len(selected) != 1:
            raise HLT16CampaignStoreError("historical checkpoint is unavailable")
        return selected[0]

    def inspect_recovery(
        self,
        *,
        authorities: Mapping[str, AuthenticatedGen0Import] | None = None,
    ) -> RecoveryStatus:
        """Inspect restart eligibility without repairing any crash suffix.

        This is deliberately distinct from :meth:`recover`.  It validates the
        complete store through :meth:`authenticated_snapshot`, which proves a
        completed publication stage is safe to read *without unlinking it*.
        Consequently it is suitable for a preflight certificate: it neither
        acquires a writer lease nor repairs, publishes, removes, or renames a
        store leaf.  A caller wanting to resume must still invoke the explicit
        recovery/writer path after its authority has been checked.
        """
        snapshot = self.authenticated_snapshot(authorities=authorities)
        latest = snapshot.checkpoint
        writer = self._active_writer()
        active = writer is not None
        trusted_checkpoint_sha256s = {latest.sha256}
        # The snapshot has already verified the entire checkpoint chain.  A
        # lease may bind an ancestor, but this read-only view need only reject
        # a lease whose declared authority cannot possibly bind this lineage.
        if writer is not None and (
            writer["authorization_commit"] != latest.authorization_commit
            or writer["plan_sha256"] != latest.plan_sha256
        ):
            raise HLT16CampaignStoreError("active writer lock authority differs")
        if (
            writer is not None
            and writer["checkpoint_sha256"] not in trusted_checkpoint_sha256s
        ):
            # ``authenticated_snapshot`` exposes only the latest checkpoint;
            # use the full non-mutating validation to admit a verified ancestor
            # lease without repairing publication stages.
            records, checkpoints, _suffix = self._validate_full_store(
                authorities,
                repair_stages=False,
            )
            del records
            trusted_checkpoint_sha256s = {item.sha256 for item in checkpoints}
            if writer["checkpoint_sha256"] not in trusted_checkpoint_sha256s:
                raise HLT16CampaignStoreError("active writer lock authority differs")
        checkpoint_is_terminal = latest.disposition in {
            "scientific_terminal", "invalid_terminal"
        }
        terminal = snapshot.terminal_lock_present
        if terminal and not checkpoint_is_terminal:
            raise HLT16CampaignStoreError("terminal lock/checkpoint differs")
        writer_state = self._writer_liveness(writer) if writer is not None else None
        suffix_state = (
            "recognized_staging_suffix"
            if snapshot.staging_paths
            else "recognized_uncheckpointed_journal_suffix"
            if snapshot.suffix_records
            else "recognized_uncheckpointed_state_suffix"
            if snapshot.orphan_descriptor_sha256 is not None
            or snapshot.orphan_payload_semantic_sha256 is not None
            else None
        )
        state = suffix_state or (
            "terminal"
            if terminal
            else "terminal_checkpoint_unlocked"
            if checkpoint_is_terminal
            else writer_state or "clean_checkpoint"
        )
        return RecoveryStatus(
            authorization_commit=latest.authorization_commit,
            plan_sha256=latest.plan_sha256,
            checkpoint_generation=latest.generation,
            checkpoint_sha256=latest.sha256,
            member_times={
                key: latest.members[key].cursor["accepted_boundary_time"]["binary64_hex"]
                for key in MEMBER_KEYS
            },
            member_modes={key: latest.members[key].cursor["mode"] for key in MEMBER_KEYS},
            active_target=dict(latest.target),
            last_complete_journal_sequence=latest.journal_sequence,
            journal_tip_sha256=latest.journal_tip_sha256,
            state=state,
            safe_to_restart=state == "clean_checkpoint",
            terminal=terminal,
            active_write=active,
            writer_state=writer_state,
        )

    def runtime_checkpoint_present(self) -> bool:
        """Report whether the store contains any post-GEN0 checkpoint leaf.

        This is deliberately an inventory discriminator, not authentication:
        callers must immediately validate a positive result with
        :meth:`authenticated_snapshot`.  Its purpose is to make reopening the
        historical raw import impossible once HLT16 has begun publishing its
        own checkpoint lineage, including if that lineage is malformed.
        """
        leaves, _directories = _walk_tree(self.root)
        prefix = "checkpoints/"
        sealed = f"{SEALED_GEN0_CHECKPOINT_SHA256}.json"
        return any(
            item.startswith(prefix)
            and Path(item).name != f"00000000000000000000-{sealed}"
            for item in leaves
        )

    def _expected_bootstrap_bridge(
        self,
        *,
        plan_sha256: str,
        authorization_commit: str,
        campaign_id: str,
        snapshots: Mapping[str, HLT16MemberSnapshot],
        member_factory: Callable[[str, str], MemberCampaignState],
        target: Mapping[str, str],
        event: int,
    ) -> _ExpectedBootstrapBridge:
        """Construct every bootstrap byte before the first bridge write."""
        publications: list[tuple[str, bytes]] = []
        publication_bytes: dict[str, bytes] = {}

        def publish_once(path: str, raw: bytes) -> None:
            prior = publication_bytes.get(path)
            if prior is None:
                publication_bytes[path] = raw
                publications.append((path, raw))
            elif prior != raw:
                raise HLT16CampaignStoreError("bootstrap duplicate publication bytes differ")
        members: dict[str, MemberCampaignState] = {}
        for key in MEMBER_KEYS:
            snapshot = snapshots[key]
            metadata = validate_codec_metadata(snapshot.metadata)
            if metadata["runtime_identity"]["member_key"] != key:
                raise HLT16CampaignStoreError("GEN0 snapshot member differs")
            payload, semantic, manifest = encode_payload(snapshot.arrays)
            raw_hash = hashlib.sha256(payload).hexdigest()
            body = {
                "schema": STATE_SCHEMA,
                "semantic_sha256": semantic,
                "raw_archive_sha256": raw_hash,
                "arrays": list(manifest),
                "metadata": metadata,
            }
            address = digest(body)
            descriptor = {**body, "descriptor_sha256": address}
            publish_once(f"payloads/{semantic}.npz", payload)
            publish_once(f"states/{address}.json", canonical(descriptor))
            member = member_factory(key, address)
            member.validate()
            members[key] = member
        payload = {
            "external_protocol": "FGC-2-SF1-PROTO17",
            "external_checkpoint_sha256": SEALED_GEN0_CHECKPOINT_SHA256,
            "plan_sha256": plan_sha256,
            "member_descriptors": {key: members[key].descriptor_sha256 for key in MEMBER_KEYS},
            "bridge_cursors": {key: dict(members[key].cursor) for key in MEMBER_KEYS},
        }
        record = journal_envelope(
            sequence=0, campaign_id=campaign_id, generation=1, event=event,
            previous_record_sha256=ROOT_SHA256, kind="genesis_bridge", payload=payload,
        )
        publish_once(
            f"journal/{0:020d}-{record['record_sha256']}.journal", canonical(record),
        )
        checkpoint = CampaignCheckpoint(
            plan_sha256=plan_sha256,
            authorization_commit=authorization_commit,
            protocol="FGC-2-SF1-PROTO18",
            campaign_id=campaign_id,
            event=event,
            target=dict(target),
            members=members,
            journal_tip_sha256=record["record_sha256"],
            parent_sha256=SEALED_GEN0_CHECKPOINT_SHA256,
            generation=1,
            journal_sequence=0,
        )
        publish_once(
            f"checkpoints/{checkpoint.generation:020d}-{checkpoint.sha256}.json",
            canonical(checkpoint.mapping()),
        )
        return _ExpectedBootstrapBridge(
            publications=tuple(publications), members=members,
            record=record, checkpoint=checkpoint,
        )

    def _validate_bootstrap_prefix(
        self,
        source: Mapping[str, str],
        expected: _ExpectedBootstrapBridge,
    ) -> None:
        """Prove existing bridge leaves are one exact deterministic prefix.

        This is invoked while ``bootstrap.guard`` is held and before any
        payload, descriptor, journal, or checkpoint publication.  A one-link
        stage is deliberately *not* adopted: without its final hard link there
        is no store-internal fact proving that the current process created the
        inode, even when its deterministic name and bytes happen to match.
        Such a suffix is retained as ambiguous evidence and the live store
        fails closed.  A two-link stage must already be byte-identical to its
        expected final object and may therefore be completed deterministically.
        """
        leaves, directories = _walk_tree(self.root)
        expected_gen0 = {
            "receipts/generation-zero.json",
            f"checkpoints/00000000000000000000-{SEALED_GEN0_CHECKPOINT_SHA256}.json",
            *(f"states/{address}.json" for address in source.values()),
        }
        allowed_directories = {"checkpoints", "receipts", "states", "payloads", "journal", "locks"}
        if (not {"checkpoints", "receipts", "states"}.issubset(directories)
                or not set(directories).issubset(allowed_directories)):
            raise HLT16CampaignStoreError("bootstrap prefix directory inventory differs")
        expected_bytes = dict(expected.publications)
        order = tuple(path for path, _raw in expected.publications)
        observed = set(leaves)
        stages: dict[str, _ValidatedLinkedStage] = {}
        for relative in observed - expected_gen0:
            path = Path(relative)
            if len(path.parts) != 2:
                raise HLT16CampaignStoreError("bootstrap prefix leaf differs")
            directory, name = path.parts
            if directory == "locks":
                if name != "bootstrap.guard":
                    raise HLT16CampaignStoreError("bootstrap prefix owns a non-bootstrap lease")
                continue
            parsed = _stage_target(name)
            if parsed is not None:
                final, digest_value = parsed
                target = f"{directory}/{final}"
                if target not in expected_bytes or hashlib.sha256(expected_bytes[target]).hexdigest() != digest_value:
                    raise HLT16CampaignStoreError("bootstrap stage target differs")
                if target in stages:
                    raise HLT16CampaignStoreError("bootstrap has duplicate publication stages")
                try:
                    stages[target] = _validated_stage(self.root, relative)
                except HLT16CampaignStoreError:
                    raise
                continue
            if relative not in expected_bytes:
                raise HLT16CampaignStoreError("bootstrap prefix owns a foreign bridge leaf")

        finals = {path for path in order if path in observed}
        prefix_length = 0
        for path in order:
            if path in finals:
                if prefix_length != order.index(path):
                    raise HLT16CampaignStoreError("bootstrap final publications are not a prefix")
                prefix_length += 1
            elif any(later in finals for later in order[order.index(path) + 1:]):
                raise HLT16CampaignStoreError("bootstrap final publications are not a prefix")
            else:
                break
        if len(stages) > 1:
            raise HLT16CampaignStoreError("bootstrap has multiple publication stages")
        if stages:
            target, stage = next(iter(stages.items()))
            index = order.index(target)
            if stage.link_count == 2:
                if index != prefix_length - 1:
                    raise HLT16CampaignStoreError("linked bootstrap stage is not the last publication")
            elif stage.link_count == 1:
                # A deterministic filename is not an ownership credential.
                # Before link(stage, final), an interrupted task stage and a
                # same-named foreign inode are observationally identical in a
                # same-user writable directory.  Never delete, overwrite, or
                # adopt that inode.  Recovery must start from a fresh verified
                # copy of the immutable GEN0 authority while preserving this
                # tree as failure evidence.
                raise HLT16CampaignStoreError(
                    "unlinked bootstrap stage ownership is unprovable"
                )
            else:  # construction is exhaustive, retained as a fail-closed guard
                raise HLT16CampaignStoreError("bootstrap stage link count differs")
        for path in finals:
            if path in stages:
                # _validated_stage already tied this two-link inode to the
                # expected raw hash.  A one-link stage cannot coexist with a
                # final object.
                continue
            raw = _safe_read(self.root, path, "bootstrap final")
            if raw != expected_bytes[path]:
                raise HLT16CampaignStoreError("bootstrap bridge leaf bytes differ")

    def adopt_generation_zero(self, *, plan_sha256: str, authorization_commit: str, campaign_id: str,
                              snapshots: Mapping[str, HLT16MemberSnapshot], authorities: Mapping[str, AuthenticatedGen0Import],
                              member_factory: Callable[[str, str], MemberCampaignState], target: Mapping[str, str], event: int = 23) -> CampaignCheckpoint:
        """Persist six authenticated payloads, then bind one generation-one bridge."""
        # Everything before this point is read-only.  In particular, do not
        # create payload, descriptor, journal, or checkpoint objects until a
        # store-owned exclusive bootstrap lease is held.  GEN0 has no HLT16
        # checkpoint yet, so the ordinary checkpoint-bound writer lease is
        # intentionally unavailable at this boundary.
        # Do not let an ordinary legacy preflight repair a staged publication:
        # a bootstrap caller owns no exclusive lease yet.  An interrupted
        # stage therefore fails closed here and is handled only after the
        # bootstrap guard has been acquired.
        source = self._legacy_inventory(authorities, repair_stages=False)
        if set(snapshots) != set(MEMBER_KEYS) or not callable(member_factory):
            raise HLT16CampaignStoreError("GEN0 bridge inputs differ")
        with self._bootstrap_guard():
            # A checkpoint-bound lease can only be legitimate after genesis;
            # seeing one while deciding genesis is a foreign or stale store
            # mutation, never authority to proceed around it.
            if self._active_writer() is not None:
                raise HLT16CampaignStoreError("GEN0 bridge is already writer-owned")
            expected = self._expected_bootstrap_bridge(
                plan_sha256=plan_sha256,
                authorization_commit=authorization_commit,
                campaign_id=campaign_id,
                snapshots=snapshots,
                member_factory=member_factory,
                target=target,
                event=event,
            )
            leaves, _directories = _walk_tree(self.root)
            has_complete_checkpoint = any(
                item.startswith("checkpoints/")
                and not Path(item).name.startswith("00000000000000000000-")
                and _stage_target(Path(item).name) is None
                for item in leaves
            )
            has_stage = any(_stage_target(Path(item).name) is not None for item in leaves)
            if not (has_complete_checkpoint and not has_stage):
                self._validate_bootstrap_prefix(source, expected)
            # Only after expected-byte membership has been proved may a
            # completed bridge inspect/repair its own deterministic stages.
            # A foreign two-link stage must never disappear before that proof.
            existing = self._existing_bridge_or_none(authorities)
            if existing is not None:
                if (
                    existing.plan_sha256 != plan_sha256
                    or existing.authorization_commit != authorization_commit
                    or existing.campaign_id != campaign_id
                    or existing.event != event
                    or existing.target != dict(target)
                ):
                    raise HLT16CampaignStoreError("generation-one bridge authority differs")
                return existing
            self._make_runtime_dirs_exclusive()
            for key in MEMBER_KEYS:
                snapshot = snapshots[key]
                persisted = self._states.persist_snapshot(snapshot, gen0_authority=authorities[key])
                expected_address = expected.members[key].descriptor_sha256
                expected_raw = dict(expected.publications)[f"states/{expected_address}.json"]
                if persisted.descriptor_sha256 != expected_address or canonical(persisted.descriptor) != expected_raw:
                    raise HLT16CampaignStoreError("GEN0 persistence differs from precomputed bridge")
            record = expected.record
            checkpoint = expected.checkpoint
            self._publish(
                f"journal/{0:020d}-{record['record_sha256']}.journal",
                dict(expected.publications)[f"journal/{0:020d}-{record['record_sha256']}.journal"],
                phase="journal",
            )
            self._publish_checkpoint(checkpoint, previous=None, records=(record,))
            _records, observed, suffix = self._validate_full_store(authorities)
            if observed[-1].sha256 != checkpoint.sha256 or suffix is not None:
                raise HLT16CampaignStoreError("generation-one bridge did not close cleanly")
            return checkpoint

    # The name records the recovery guarantee explicitly.  Repeating this
    # exact call after a publication cut reuses byte-identical state objects,
    # the bridge record, and the checkpoint; any different bytes are a fork.
    initialize_or_recover_bridge = adopt_generation_zero

    def _require_writer_capability(
        self,
        capability: WriterCapability,
        checkpoint: CampaignCheckpoint,
    ) -> None:
        """Prove a mutator belongs to the live lease and current predecessor."""
        if not isinstance(capability, WriterCapability) or capability is not self._writer_capability:
            raise HLT16CampaignStoreError("writer capability is absent or foreign")
        if capability.host != os.uname().nodename or capability.pid != os.getpid():
            raise HLT16CampaignStoreError("writer capability process differs")
        writer = self._active_writer()
        if writer is None or any(
            writer[field] != getattr(capability, field)
            for field in ("authorization_commit", "plan_sha256", "owner_token", "host", "pid")
        ):
            raise HLT16CampaignStoreError("writer capability lease differs")
        if (checkpoint.authorization_commit != capability.authorization_commit
                or checkpoint.plan_sha256 != capability.plan_sha256):
            raise HLT16CampaignStoreError("writer capability authority differs")

    def active_writer_capability(self, checkpoint: CampaignCheckpoint) -> WriterCapability:
        """Return the live capability for the durable predecessor checkpoint.

        A recognized publication suffix is precisely when recovery needs this
        capability.  The durable checkpoint must still be the latest closed
        checkpoint, but an authenticated payload/journal/stage suffix does not
        invalidate the lease that owns its completion.
        """
        with self._writer_guard():
            _records, checkpoints, _suffix = self._validate_full_store()
            if checkpoints[-1].sha256 != checkpoint.sha256:
                raise HLT16CampaignStoreError("writer capability checkpoint is not current")
            capability = self._writer_capability
            if capability is None:
                raise HLT16CampaignStoreError("writer capability is absent")
            self._require_writer_capability(capability, checkpoint)
            return capability

    def persist_snapshot(
        self,
        snapshot: HLT16MemberSnapshot,
        *,
        capability: WriterCapability,
        checkpoint: CampaignCheckpoint,
    ):
        """Persist an evolved payload only under the current active lease."""
        with self._writer_guard():
            _records, checkpoints, suffix = self._validate_full_store()
            if suffix is not None or checkpoints[-1].sha256 != checkpoint.sha256:
                raise HLT16CampaignStoreError("writer payload predecessor is not current")
            self._require_writer_capability(capability, checkpoint)
            try:
                return self._states._persist_snapshot_internal(snapshot)
            except HLT16StateStoreError as exc:
                raise HLT16CampaignStoreError("writer payload publication failed safely") from exc

    def load_state(self, descriptor_sha256: str):
        """Read one already committed payload; this API has no write surface."""
        try:
            return self._states.load_state(descriptor_sha256)
        except HLT16StateStoreError as exc:
            raise HLT16CampaignStoreError("campaign state cannot be restored") from exc

    def publish_journal(self, record: Mapping[str, Any], *, capability: WriterCapability) -> str:
        """Durably append one already schema-valid, next-sequence record."""
        with self._writer_guard():
            existing, checkpoints, _suffix = self._validate_full_store()
            latest = checkpoints[-1]
            self._require_writer_capability(capability, latest)
            try:
                checked = validate_journal(record)
            except HLT16CampaignSchemaError as exc:
                raise HLT16CampaignStoreError("journal record schema differs") from exc
        # A coordinator may crash after fsyncing its rejection record but
        # before returning to its caller.  Retrying that *exact* append is a
        # recovery operation, not a second event.  Any same-sequence variant
        # remains a fork and is rejected before it can become a convenient
        # alternative retry history.
            if checked["sequence"] < len(existing):
                prior = existing[checked["sequence"]]
                if canonical(prior) != canonical(checked):
                    raise HLT16CampaignStoreError("journal replay forks existing record")
                return checked["record_sha256"]
            expected_parent = ROOT_SHA256 if not existing else existing[-1]["record_sha256"]
            if checked["sequence"] != len(existing) or checked["previous_record_sha256"] != expected_parent:
                raise HLT16CampaignStoreError("journal append relation differs")
            if (checked["campaign_id"] != latest.campaign_id
                    or checked["generation"] != latest.generation + 1
                    or checked["event"] != latest.event):
                raise HLT16CampaignStoreError("journal append authority differs")
            self._publish(f"journal/{checked['sequence']:020d}-{checked['record_sha256']}.journal", canonical(checked), phase="journal")
            return checked["record_sha256"]

    def publish_checkpoint(self, checkpoint: CampaignCheckpoint, *, previous: CampaignCheckpoint | None,
                           records: tuple[Mapping[str, Any], ...], capability: WriterCapability) -> None:
        """Narrow public primitive for a coordinator's schema-valid successor."""
        with self._writer_guard():
            _all_records, checkpoints, _suffix = self._validate_full_store()
            if previous is None or checkpoints[-1].sha256 != previous.sha256:
                raise HLT16CampaignStoreError("checkpoint predecessor is not current")
            self._require_writer_capability(capability, previous)
            self._publish_checkpoint(checkpoint, previous=previous, records=records)
            _records, observed, suffix = self._validate_full_store()
            if observed[-1].sha256 != checkpoint.sha256 or suffix is not None:
                raise HLT16CampaignStoreError("published checkpoint did not close its suffix")

    def _publish_checkpoint(self, checkpoint: CampaignCheckpoint, *, previous: CampaignCheckpoint | None, records: tuple[Mapping[str, Any], ...]) -> None:
        try:
            if previous is None:
                validate_generation_one_bridge(checkpoint, records[0], external_checkpoint_sha256=SEALED_GEN0_CHECKPOINT_SHA256)
            else:
                validate_checkpoint_successor(previous, checkpoint, records)
        except HLT16CampaignSchemaError as exc:
            raise HLT16CampaignStoreError("checkpoint successor differs") from exc
        for key in MEMBER_KEYS:
            self._states.load_state(checkpoint.members[key].descriptor_sha256)
        self._publish(f"checkpoints/{checkpoint.generation:020d}-{checkpoint.sha256}.json", canonical(checkpoint.mapping()), phase="checkpoint")

    @contextmanager
    def _bootstrap_guard(self):
        """Serialize GEN0 adoption before a checkpoint-bound lease exists.

        The guard is deliberately distinct from ``active-write.lock``:
        authorizing that lock requires a generation-one checkpoint, while this
        guard is what makes the first such checkpoint unique.  It remains an
        empty advisory-lock inode after release so later validators can see a
        stable, non-following inventory without treating a process-local flock
        as scientific evidence.
        """
        with self._lease_guard("bootstrap.guard"):
            yield

    @contextmanager
    def _writer_guard(self):
        """Serialize acquisition, owner-bound release, and stale takeover."""
        with self._lease_guard("writer.guard"):
            yield

    @contextmanager
    def _lease_guard(self, name: str):
        """Take one store-owned, non-following advisory lock inode."""
        if name not in {"bootstrap.guard", "writer.guard"}:
            raise HLT16CampaignStoreError("guard name differs")
        directory_flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
        root_fd = os.open(self.root, directory_flags)
        locks_fd = -1
        descriptor = -1
        try:
            # Sealed GEN0 predates HLT16 and may not contain the empty locks
            # directory.  Creating that directory is the bootstrap primitive
            # itself, not a bridge object: the following non-following flock
            # is acquired before any payload/descriptor/journal/checkpoint
            # mutation.  Concurrent mkdir is harmless; every contender then
            # serializes on this one inode.
            try:
                os.mkdir("locks", mode=0o755, dir_fd=root_fd)
                os.fsync(root_fd)
            except FileExistsError:
                pass
            locks_fd = os.open("locks", directory_flags, dir_fd=root_fd)
            descriptor = os.open(
                name,
                os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0),
                0o600,
                dir_fd=locks_fd,
            )
            info = os.fstat(descriptor)
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise HLT16CampaignStoreError("store guard is unsafe")
            fcntl.flock(descriptor, fcntl.LOCK_EX)
            yield
        finally:
            try:
                if descriptor != -1:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)
            finally:
                if descriptor != -1:
                    os.close(descriptor)
                if locks_fd != -1:
                    os.close(locks_fd)
                os.close(root_fd)

    def _existing_bridge_or_none(
        self,
        authorities: Mapping[str, AuthenticatedGen0Import],
    ) -> CampaignCheckpoint | None:
        """Return a closed generation-one bridge, or prove no bridge exists.

        This runs only while the bootstrap guard is held.  A journal/state
        suffix with no checkpoint is an interrupted copy of the *same* bridge:
        exact content-addressed publication may deterministically resume it.
        A checkpoint, however, is already a durable genesis and must validate
        completely before it can be reused.
        """
        checkpoints = self._checkpoints(repair_stages=True)
        if not checkpoints:
            return None
        records, observed, suffix = self._validate_full_store(authorities)
        if suffix is not None or len(observed) != 1:
            raise HLT16CampaignStoreError("existing generation-one bridge is not closed")
        return observed[0]

    @staticmethod
    def _owner_token(value: object) -> str:
        if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
            raise HLT16CampaignStoreError("writer owner token differs")
        try:
            int(value, 16)
        except ValueError as exc:
            raise HLT16CampaignStoreError("writer owner token differs") from exc
        return value

    def acquire_active_writer(
        self,
        checkpoint: CampaignCheckpoint,
        *,
        host: str,
        pid: int | None = None,
        owner_token: str | None = None,
    ) -> str:
        """Atomically lease only the exact current nonterminal checkpoint."""
        checkpoint.validate()
        local = os.uname().nodename
        if host != local:
            raise HLT16CampaignStoreError("writer host is not this machine")
        owner = os.getpid() if pid is None else pid
        if isinstance(owner, bool) or not isinstance(owner, int) or owner <= 0:
            raise HLT16CampaignStoreError("writer pid differs")
        token = secrets.token_hex(32) if owner_token is None else self._owner_token(owner_token)
        with self._writer_guard():
            _records, checkpoints, suffix = self._validate_full_store()
            latest = checkpoints[-1]
            if (latest.sha256 != checkpoint.sha256 or suffix is not None
                    or latest.disposition in {"scientific_terminal", "invalid_terminal"}
                    or self._terminal_lock(latest) or self._active_writer() is not None):
                raise HLT16CampaignStoreError("writer checkpoint is not exclusively runnable")
            payload = {
                "schema": "FGC-1-HLT16-active-writer-v1",
                "authorization_commit": checkpoint.authorization_commit,
                "plan_sha256": checkpoint.plan_sha256,
                "checkpoint_sha256": checkpoint.sha256,
                "host": host,
                "pid": owner,
                "owner_token": token,
            }
            self._publish("locks/active-write.lock", canonical(payload), phase="active_lock")
            self._writer_capability = WriterCapability(
                authorization_commit=checkpoint.authorization_commit,
                plan_sha256=checkpoint.plan_sha256,
                owner_token=token, host=host, pid=owner,
            )
        return token

    def recover_interrupted_writer_acquisition(
        self,
        checkpoint: CampaignCheckpoint,
        *,
        host: str,
    ) -> str:
        """Recover only the deterministic stage for ``active-write.lock``.

        A process can die while publishing its lease, before the final hard
        link exists.  A complete one-link stage is replayed byte-for-byte so
        ordinary liveness/takeover rules can decide its owner.  A provably
        partial task-owned stage is removed because no complete lease exists.
        No other publication stage is touched by this narrow operation.
        """
        checkpoint.validate()
        if host != os.uname().nodename:
            raise HLT16CampaignStoreError("writer-recovery host is not this machine")
        with self._writer_guard():
            _records, checkpoints, suffix = self._validate_full_store()
            latest = checkpoints[-1]
            if latest.sha256 != checkpoint.sha256:
                raise HLT16CampaignStoreError("writer-recovery checkpoint is not current")
            if self._active_writer() is not None:
                return "active_writer_present"
            snapshot = self.authenticated_snapshot()
            if not snapshot.staging_paths:
                return "no_writer_stage"
            if len(snapshot.staging_paths) != 1:
                raise HLT16CampaignStoreError("campaign has multiple publication stages")
            relative = snapshot.staging_paths[0]
            parsed = _stage_target(Path(relative).name)
            if (parsed is None or parsed[0] != "active-write.lock"
                    or Path(relative).parent.as_posix() != "locks"):
                raise HLT16CampaignStoreError("campaign has a non-writer publication stage")
            _validated_stage(self.root, relative)
            _final, expected = parsed
            staged_raw = _safe_read(self.root, relative, "writer publication stage")
            observed = hashlib.sha256(staged_raw).hexdigest()
            path = self.root / relative
            if observed != expected:
                # Recheck the exact one-link inode immediately before removing
                # only this deterministic, incomplete task stage.
                before = path.lstat()
                if (stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode)
                        or before.st_nlink != 1):
                    raise HLT16CampaignStoreError("partial writer stage is unsafe")
                directory = os.open(
                    path.parent,
                    os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
                )
                try:
                    current = os.stat(path.name, dir_fd=directory, follow_symlinks=False)
                    if (before.st_dev, before.st_ino, before.st_size,
                            before.st_mtime_ns, before.st_ctime_ns) != (
                            current.st_dev, current.st_ino, current.st_size,
                            current.st_mtime_ns, current.st_ctime_ns):
                        raise HLT16CampaignStoreError("partial writer stage raced")
                    os.unlink(path.name, dir_fd=directory)
                    os.fsync(directory)
                finally:
                    os.close(directory)
                return "discarded_partial_writer_stage"
            raw = staged_raw
            value = _decode(raw, "complete writer stage")
            trusted = {item.sha256 for item in checkpoints}
            if (set(value) != {
                    "schema", "authorization_commit", "plan_sha256",
                    "checkpoint_sha256", "host", "pid", "owner_token",
                }
                    or value["schema"] != "FGC-1-HLT16-active-writer-v1"
                    or value["authorization_commit"] != latest.authorization_commit
                    or value["plan_sha256"] != latest.plan_sha256
                    or value["checkpoint_sha256"] not in trusted
                    or value["host"] != host
                    or isinstance(value["pid"], bool)
                    or not isinstance(value["pid"], int)
                    or value["pid"] <= 0):
                raise HLT16CampaignStoreError("complete writer stage authority differs")
            self._owner_token(value["owner_token"])
            self._publish("locks/active-write.lock", raw, phase="active_lock_recovery")
            return "replayed_complete_writer_stage"

    def _retire_verified_active_writer(self, writer: Mapping[str, Any]) -> None:
        """Atomically retire, never unlink, a writer lease under ``writer.guard``.

        POSIX has no compare-and-unlink operation.  Moving the active pathname
        to a deterministic content-addressed retired leaf means a replacement
        at the final pathname is preserved as evidence: it is moved, reread,
        and rejected rather than deleted.  The proven retired lease is also
        deliberately retained, so there is no second path-based deletion race.
        """
        flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
        root_fd = locks_fd = -1
        try:
            root_fd = os.open(self.root, flags)
            locks_fd = os.open("locks", flags, dir_fd=root_fd)
            before = os.stat("active-write.lock", dir_fd=locks_fd, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
                raise HLT16CampaignStoreError("active writer lock changed before retirement")
            raw = _read_leaf_fd(locks_fd, "active-write.lock", "active writer lock")
            observed = _decode(raw, "active writer lock")
            if dict(observed) != dict(writer):
                raise HLT16CampaignStoreError("active writer lock changed before retirement")
            current = os.stat("active-write.lock", dir_fd=locks_fd, follow_symlinks=False)
            if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
                current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns, current.st_ctime_ns,
            ):
                raise HLT16CampaignStoreError("active writer lock raced before retirement")
            retired = ".active-write.lock.hlt16-quarantine-" + hashlib.sha256(raw).hexdigest()
            try:
                os.stat(retired, dir_fd=locks_fd, follow_symlinks=False)
            except FileNotFoundError:
                pass
            else:
                raise HLT16CampaignStoreError("writer retirement identity already exists")
            self._cut("before_writer_retirement_move")
            os.rename("active-write.lock", retired, src_dir_fd=locks_fd, dst_dir_fd=locks_fd)
            os.fsync(locks_fd)
            retired_raw = _read_leaf_fd(locks_fd, retired, "retired writer lease")
            if retired_raw != raw or _decode(retired_raw, "retired writer lease") != dict(writer):
                raise HLT16CampaignStoreError("retired writer lease differs; retained fail-closed")
            try:
                os.stat("active-write.lock", dir_fd=locks_fd, follow_symlinks=False)
            except FileNotFoundError:
                return
            raise HLT16CampaignStoreError("replacement active writer lease retained fail-closed")
        except OSError as exc:
            raise HLT16CampaignStoreError("active writer lock retirement failed safely") from exc
        finally:
            if locks_fd != -1:
                os.close(locks_fd)
            if root_fd != -1:
                os.close(root_fd)

    def release_active_writer(self, *, owner_token: str, host: str, pid: int | None = None) -> None:
        token = self._owner_token(owner_token)
        owner = os.getpid() if pid is None else pid
        with self._writer_guard():
            writer = self._active_writer()
            if writer is None:
                raise HLT16CampaignStoreError("active writer lock is absent")
            if writer["owner_token"] != token or writer["host"] != host or writer["pid"] != owner:
                raise HLT16CampaignStoreError("active writer ownership differs")
            self._retire_verified_active_writer(writer)
            self._writer_capability = None

    def _active_writer(self) -> Mapping[str, Any] | None:
        path = self.root / "locks" / "active-write.lock"
        try:
            info = path.lstat()
        except FileNotFoundError:
            return None
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
            raise HLT16CampaignStoreError("active writer lock unsafe")
        value = _decode(_safe_read(self.root, "locks/active-write.lock", "active writer lock"), "active writer lock")
        if (set(value) != {"schema", "authorization_commit", "plan_sha256", "checkpoint_sha256", "host", "pid", "owner_token"}
                or value["schema"] != "FGC-1-HLT16-active-writer-v1"
                or not isinstance(value["host"], str) or not value["host"]
                or isinstance(value["pid"], bool) or not isinstance(value["pid"], int) or value["pid"] <= 0):
            raise HLT16CampaignStoreError("active writer lock schema differs")
        self._owner_token(value["owner_token"])
        return value

    def _writer_liveness(self, writer: Mapping[str, Any]) -> str:
        """Only a missing PID on this exact host is a verified stale lease."""
        local_host = os.uname().nodename
        if writer["host"] != local_host:
            return "foreign_host_writer"
        try:
            os.kill(int(writer["pid"]), 0)
        except ProcessLookupError:
            return "stale_verified_lock"
        except PermissionError:
            return "live_local_writer"
        return "live_local_writer"

    def take_over_stale_writer(
        self,
        checkpoint: CampaignCheckpoint,
        *,
        host: str,
        pid: int | None = None,
        owner_token: str | None = None,
    ) -> str:
        """Replace only a same-host dead-PID lease anchored in this exact chain.

        A writer lease is acquired at one authenticated checkpoint and remains
        constant while that owner appends successor checkpoints.  Its bound
        checkpoint must therefore be an ancestor present in the validated
        append-only chain, not necessarily the latest checkpoint.
        """
        local = os.uname().nodename
        if host != local:
            raise HLT16CampaignStoreError("takeover host is not this machine")
        owner = os.getpid() if pid is None else pid
        token = secrets.token_hex(32) if owner_token is None else self._owner_token(owner_token)
        with self._writer_guard():
            _records, checkpoints, suffix = self._validate_full_store()
            writer = self._active_writer()
            # A dead process may have left an independently recognized
            # publication suffix.  Recovery must acquire exclusive ownership
            # before replaying that suffix, so its very presence cannot forbid
            # a same-host stale-lease takeover.  The full store validator above
            # has already rejected every ambiguous, unsafe, or forked suffix.
            if (checkpoints[-1].sha256 != checkpoint.sha256
                    or writer is None or self._writer_liveness(writer) != "stale_verified_lock"):
                raise HLT16CampaignStoreError("writer lease is not verified stale/current")
            trusted_checkpoint_sha256s = {item.sha256 for item in checkpoints}
            if (writer["authorization_commit"] != checkpoint.authorization_commit
                    or writer["plan_sha256"] != checkpoint.plan_sha256
                    or writer["checkpoint_sha256"] not in trusted_checkpoint_sha256s):
                raise HLT16CampaignStoreError("stale writer authority differs")
            self._retire_verified_active_writer(writer)
            payload = {
                "schema": "FGC-1-HLT16-active-writer-v1",
                "authorization_commit": checkpoint.authorization_commit,
                "plan_sha256": checkpoint.plan_sha256,
                "checkpoint_sha256": checkpoint.sha256,
                "host": host,
                "pid": owner,
                "owner_token": token,
            }
            self._publish("locks/active-write.lock", canonical(payload), phase="active_lock")
            self._writer_capability = WriterCapability(
                authorization_commit=checkpoint.authorization_commit,
                plan_sha256=checkpoint.plan_sha256,
                owner_token=token, host=host, pid=owner,
            )
        return token

    def publish_terminal_lock(self, checkpoint: CampaignCheckpoint, *, capability: WriterCapability) -> None:
        """Bind a terminal lease to the already durable terminal checkpoint."""
        if checkpoint.disposition not in {"scientific_terminal", "invalid_terminal"}:
            raise HLT16CampaignStoreError("terminal checkpoint disposition differs")
        with self._writer_guard():
            records, checkpoints, suffix = self._validate_full_store()
            if not checkpoints or checkpoints[-1].sha256 != checkpoint.sha256 or suffix is not None:
                raise HLT16CampaignStoreError("terminal checkpoint is not current")
            self._require_writer_capability(capability, checkpoint)
            body = {"schema": "FGC-1-HLT16-terminal-lock-v1", "checkpoint_sha256": checkpoint.sha256,
                    "journal_tip_sha256": checkpoint.journal_tip_sha256}
            payload = {**body, "lock_sha256": hashlib.sha256(canonical(body)).hexdigest()}
            self._publish("locks/terminal.lock", canonical(payload), phase="terminal_lock")

    def _terminal_lock(
        self,
        latest: CampaignCheckpoint,
        *,
        snapshot_stages: Mapping[str, _ValidatedLinkedStage] | None = None,
    ) -> bool:
        path = self.root / "locks" / "terminal.lock"
        try:
            info = path.lstat()
        except FileNotFoundError:
            return False
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
            raise HLT16CampaignStoreError("terminal lock unsafe")
        raw = (
            _safe_read(self.root, "locks/terminal.lock", "terminal lock")
            if snapshot_stages is None
            else _safe_read_snapshot(
                self.root,
                "locks/terminal.lock",
                "terminal lock",
                snapshot_stages,
            )
        )
        value = _decode(raw, "terminal lock")
        body = dict(value); observed = body.pop("lock_sha256", None)
        if (set(body) != {"schema", "checkpoint_sha256", "journal_tip_sha256"}
                or body["schema"] != "FGC-1-HLT16-terminal-lock-v1"
                or observed != hashlib.sha256(canonical(body)).hexdigest()
                or body["checkpoint_sha256"] != latest.sha256
                or body["journal_tip_sha256"] != latest.journal_tip_sha256):
            raise HLT16CampaignStoreError("terminal lock binding differs")
        return True

    def recover(self, *, authorities: Mapping[str, AuthenticatedGen0Import] | None = None) -> RecoveryStatus:
        records, checkpoints, suffix = self._validate_full_store(authorities)
        latest = checkpoints[-1]
        terminal = self._terminal_lock(latest)
        writer = self._active_writer()
        active = writer is not None
        trusted_checkpoint_sha256s = {item.sha256 for item in checkpoints}
        if writer is not None and (writer["authorization_commit"] != latest.authorization_commit
                                   or writer["plan_sha256"] != latest.plan_sha256
                                   or writer["checkpoint_sha256"] not in trusted_checkpoint_sha256s):
            raise HLT16CampaignStoreError("active writer lock authority differs")
        checkpoint_is_terminal = latest.disposition in {
            "scientific_terminal", "invalid_terminal"
        }
        if terminal and not checkpoint_is_terminal:
            raise HLT16CampaignStoreError("terminal lock/checkpoint differs")
        writer_state = self._writer_liveness(writer) if writer is not None else None
        state = suffix or (
            "terminal"
            if terminal
            else "terminal_checkpoint_unlocked"
            if checkpoint_is_terminal
            else writer_state or "clean_checkpoint"
        )
        return RecoveryStatus(authorization_commit=latest.authorization_commit, plan_sha256=latest.plan_sha256,
            checkpoint_generation=latest.generation, checkpoint_sha256=latest.sha256,
            member_times={key: latest.members[key].cursor["accepted_boundary_time"]["binary64_hex"] for key in MEMBER_KEYS},
            member_modes={key: latest.members[key].cursor["mode"] for key in MEMBER_KEYS}, active_target=dict(latest.target),
            last_complete_journal_sequence=latest.journal_sequence, journal_tip_sha256=latest.journal_tip_sha256,
            state=state, safe_to_restart=state == "clean_checkpoint", terminal=terminal,
            active_write=active, writer_state=writer_state)


__all__ = [
    "AuthenticatedCampaignSnapshot",
    "HLT16CampaignStore",
    "HLT16CampaignStoreError",
    "RecoveryStatus",
]
