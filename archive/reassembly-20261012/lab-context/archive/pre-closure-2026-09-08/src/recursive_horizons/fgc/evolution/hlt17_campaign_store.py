"""Append-only HLT17/PROTO19 first-event publisher.

The store consumes validated member bundles of three immutable byte-string
components and publishes each generation as one exclusive no-replace
directory.  It does not authenticate a physical origin, source, or
environment, authorize production writes, or run a campaign.  A free OS lock
is not authority to take over an unclosed session.  Reads authenticate the
unique chain and every referenced leaf without mutation.  Complete visible
bytes after publication uncertainty may be read as evidence; they do not
authorize resuming an unclosed or uncertain writer.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import fcntl
import os
from pathlib import Path
import stat
import threading
from types import MappingProxyType
from typing import Mapping

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    DEFAULT_MAX_BYTES,
    PostpublicationUncertainty,
    PrepublicationError,
    UnsafePathError,
    canonical_json_bytes,
    load_canonical_json,
    publish_exclusive_directory,
    publish_exclusive_file,
    read_regular_file,
)

from .hlt17_member_checkpoint import HLT17GenerationBundle
from .protocol_v19 import (
    EVENT_ORIGIN_HEX,
    EVENT_TARGET_HEX,
    FUTURE_ENVIRONMENT_GATE_SEAM,
    FUTURE_INDEPENDENT_BINDER_SEAM,
    FUTURE_LIVE_RUNNER_SEAM,
    FUTURE_PHYSICAL_ORIGIN_SEAM,
    FUTURE_RESOURCE_GATE_SEAM,
    FUTURE_SOURCE_MANIFEST_SEAM_TEXT,
    HLT17ProtocolError,
    HLT17ValidatedGeneration,
    PROTOCOL_ARTIFACT_ID,
    authenticate_generation_chain,
    build_root_manifest,
    build_seed_generation,
    build_successor_generation,
    encode_session_close,
    encode_session_open,
    encode_terminal_lock,
    generation_directory_name,
    parse_generation_directory_name,
    parse_session_name,
    session_close_name,
    session_open_name,
    validate_root_manifest,
    validate_session_close,
    validate_session_open,
    validate_terminal_lock,
)


ROOT_MANIFEST_NAME = "root_manifest.json"
WRITER_EXCLUSION_NAME = "writer_exclusion.lock"
TERMINAL_LOCK_NAME = "terminal.lock"
WRITER_EXCLUSION_SCHEMA = "FGC-1-HLT17-writer-exclusion-v1"
PRODUCTION_NAMESPACE_MARKERS = ("/runs/", "\\runs\\", "/runs\\", "\\runs/")
PRODUCTION_EVENT_PATH = "runs/fgc-2-sf1/pro20-event1/calibration"


class HLT17CampaignStoreError(ValueError):
    """A PROTO19 campaign store is unsafe, incomplete, forked, or unauthorized."""

    published = False


class HLT17WriterConflict(HLT17CampaignStoreError):
    """Process or thread exclusion refused a second writer."""


class HLT17UnclosedSession(HLT17CampaignStoreError):
    """A free OS lock does not authorize takeover of an unclosed session."""


class HLT17TerminalStore(HLT17CampaignStoreError):
    """A typed terminal preserves history and forbids further appends."""


class HLT17PostpublicationUncertainty(HLT17CampaignStoreError):
    """A generation directory became visible, but a later durability step failed.

    Complete visible bytes may be readable as evidence.  They are not
    permission to resume an unclosed or uncertain writer.
    """

    published = True


def _fail(message: str) -> None:
    raise HLT17CampaignStoreError(message)


def _digest(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _canonical_root(root: Path) -> Path:
    if not isinstance(root, Path):
        raise TypeError("store root must be a pathlib.Path")
    try:
        supplied = Path(os.fspath(root))
    except TypeError as error:
        raise HLT17CampaignStoreError("store root is unsafe") from error
    if not supplied.is_absolute() or any(part in {".", ".."} for part in supplied.parts):
        _fail("store root must be an absolute canonical path")
    try:
        resolved = supplied.resolve(strict=True)
    except OSError as error:
        raise HLT17CampaignStoreError("store root cannot be resolved safely") from error
    if resolved != supplied:
        _fail("store root path traverses a symlinked ancestor")
    text = resolved.as_posix()
    if PRODUCTION_EVENT_PATH in text or any(marker in text for marker in PRODUCTION_NAMESPACE_MARKERS):
        _fail("production namespace writes are not authorized")
    return resolved


def _exclusion_payload() -> bytes:
    return canonical_json_bytes(
        {
            "schema": WRITER_EXCLUSION_SCHEMA,
            "schema_version": 1,
            "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
            "role": "process_thread_exclusion",
            "production_write_authorized": False,
            "campaign_execution_authorized": False,
            "automatic_stale_writer_recovery": False,
        }
    )


def _flag(name: str) -> int:
    value = int(getattr(os, name, 0) or 0)
    if value == 0:
        raise HLT17CampaignStoreError(f"{name} is required")
    return value


def _scandir_fd(directory_fd: int):
    return os.scandir(directory_fd)


def _classify_root_entry(name: str, info: os.stat_result) -> str:
    if stat.S_ISLNK(info.st_mode):
        _fail(f"unsafe symlink: {name}")
    if stat.S_ISREG(info.st_mode):
        if info.st_nlink != 1:
            _fail(f"hard-linked leaf: {name}")
        if name in {ROOT_MANIFEST_NAME, WRITER_EXCLUSION_NAME, TERMINAL_LOCK_NAME}:
            return "file"
        try:
            parse_session_name(name)
        except HLT17ProtocolError as error:
            raise HLT17CampaignStoreError(f"unexpected file: {name}") from error
        return "session"
    if not stat.S_ISDIR(info.st_mode):
        _fail(f"unsafe node: {name}")
    try:
        parse_generation_directory_name(name)
    except HLT17ProtocolError as error:
        raise HLT17CampaignStoreError(f"unexpected directory: {name}") from error
    return "generation"


def _open_root_fd(root: Path) -> int:
    flags = os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW")
    try:
        before = root.lstat()
    except OSError as error:
        raise HLT17CampaignStoreError("store root absent or unsafe") from error
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
        _fail("store root is unsafe")
    descriptor = os.open(root, flags)
    try:
        after = os.fstat(descriptor)
        if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
            _fail("store root raced")
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _list_root(root: Path) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    descriptor = _open_root_fd(root)
    files: list[str] = []
    sessions: list[str] = []
    generations: list[str] = []
    try:
        with _scandir_fd(descriptor) as scan:
            entries = sorted(scan, key=lambda item: item.name)
        for item in entries:
            try:
                info = os.stat(item.name, dir_fd=descriptor, follow_symlinks=False)
            except OSError as error:
                raise HLT17CampaignStoreError(f"node raced: {item.name}") from error
            kind = _classify_root_entry(item.name, info)
            if kind == "file":
                files.append(item.name)
            elif kind == "session":
                sessions.append(item.name)
            else:
                generations.append(item.name)
    finally:
        os.close(descriptor)
    return tuple(files), tuple(sessions), tuple(generations)


def _read_json(root: Path, relative: str) -> tuple[bytes, dict[str, object]]:
    try:
        raw = read_regular_file(root, relative)
        value = load_canonical_json(raw)
    except (UnsafePathError, CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17CampaignStoreError(f"{relative} cannot be read safely") from error
    if type(value) is not dict:
        _fail(f"{relative} is not a JSON object")
    return raw, value


def _read_blob(root: Path, relative: str) -> bytes:
    try:
        return read_regular_file(root, relative, max_bytes=DEFAULT_MAX_BYTES)
    except (UnsafePathError, TypeError, ValueError) as error:
        raise HLT17CampaignStoreError(f"{relative} cannot be read safely") from error


def _generation_files(root: Path, directory: str) -> dict[str, bytes]:
    descriptor = _open_root_fd(root)
    files: dict[str, bytes] = {}
    child = None
    try:
        try:
            before = os.stat(directory, dir_fd=descriptor, follow_symlinks=False)
        except OSError as error:
            raise HLT17CampaignStoreError(f"generation absent: {directory}") from error
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
            _fail(f"generation is not a directory: {directory}")
        child = os.open(
            directory,
            os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW"),
            dir_fd=descriptor,
        )
        after = os.fstat(child)
        if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
            _fail(f"generation raced: {directory}")
        with _scandir_fd(child) as scan:
            entries = sorted(scan, key=lambda item: item.name)
        blob_dir = False
        for item in entries:
            info = os.stat(item.name, dir_fd=child, follow_symlinks=False)
            if stat.S_ISLNK(info.st_mode):
                _fail(f"unsafe symlink: {directory}/{item.name}")
            if item.name in {"generation_manifest.json", "checkpoint.json", "journal.json"}:
                if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                    _fail(f"unsafe generation leaf: {directory}/{item.name}")
                files[item.name] = _read_blob(root, f"{directory}/{item.name}")
                continue
            if item.name == "blobs":
                if not stat.S_ISDIR(info.st_mode):
                    _fail(f"{directory}/blobs is not a directory")
                blob_dir = True
                continue
            _fail(f"unexpected generation entry: {directory}/{item.name}")
        if blob_dir:
            blob_fd = os.open(
                "blobs",
                os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW"),
                dir_fd=child,
            )
            try:
                with _scandir_fd(blob_fd) as scan:
                    blobs = sorted(scan, key=lambda item: item.name)
                for item in blobs:
                    info = os.stat(item.name, dir_fd=blob_fd, follow_symlinks=False)
                    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                        _fail(f"unsafe blob: {directory}/blobs/{item.name}")
                    relative = f"{directory}/blobs/{item.name}"
                    files[f"blobs/{item.name}"] = _read_blob(root, relative)
            finally:
                os.close(blob_fd)
    finally:
        if child is not None:
            os.close(child)
        os.close(descriptor)
    required = {"generation_manifest.json", "checkpoint.json", "journal.json"}
    if not required.issubset(files):
        _fail(f"partial generation: {directory}")
    return files


@dataclass(frozen=True, slots=True)
class HLT17WriterCapability:
    """In-process handle issued only after session publication and exclusion."""

    campaign_id: str
    session_index: int
    owner_token: str
    host: str
    pid: int
    thread_id: int
    session_open_sha256: str


@dataclass(frozen=True, slots=True)
class HLT17AuthenticatedView:
    """Read-only reconstruction of a unique complete published chain."""

    root: Path
    root_manifest_raw: bytes
    generations: tuple[HLT17ValidatedGeneration, ...]
    unclosed_session: bool

    def __post_init__(self) -> None:
        if type(self.root_manifest_raw) is not bytes:
            raise TypeError("root manifest raw must be bytes")
        object.__setattr__(self, "generations", tuple(self.generations))
        try:
            validate_root_manifest(load_canonical_json(self.root_manifest_raw))
        except (CanonicalJSONError, HLT17ProtocolError, TypeError, ValueError) as error:
            raise HLT17CampaignStoreError("root manifest differs") from error
        if not self.generations:
            raise HLT17CampaignStoreError("campaign store has no published generations")
        last = self.generations[-1]
        last.__post_init__()

    @property
    def root_manifest(self) -> Mapping[str, object]:
        mapping = validate_root_manifest(load_canonical_json(self.root_manifest_raw))
        return MappingProxyType(mapping)

    @property
    def bundles(self) -> Mapping[str, HLT17GenerationBundle]:
        last = self.generations[-1]
        last.__post_init__()
        return last.bundles

    @property
    def terminal(self) -> bool:
        last = self.generations[-1]
        last.__post_init__()
        return last.terminal

    @property
    def first_event_complete(self) -> bool:
        last = self.generations[-1]
        last.__post_init__()
        return last.first_event_complete

    @property
    def checkpoint_sha256(self) -> str:
        last = self.generations[-1]
        last.__post_init__()
        return last.checkpoint_sha256

    @property
    def journal_sha256(self) -> str:
        last = self.generations[-1]
        last.__post_init__()
        return last.journal_sha256


def _dir_identity(root: Path) -> tuple[int, int]:
    try:
        info = root.lstat()
    except OSError as error:
        raise HLT17CampaignStoreError("store root absent or unsafe") from error
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        _fail("store root is unsafe")
    return (info.st_dev, info.st_ino)


def _lock_path_identity(root: Path) -> tuple[int, int]:
    root_fd = _open_root_fd(root)
    try:
        try:
            info = os.stat(WRITER_EXCLUSION_NAME, dir_fd=root_fd, follow_symlinks=False)
        except OSError as error:
            raise HLT17CampaignStoreError("writer exclusion lock cannot be read") from error
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            _fail("writer exclusion lock is unsafe")
        return (info.st_dev, info.st_ino)
    finally:
        os.close(root_fd)


class HLT17CampaignStore:
    """Synthetic first-event publisher. Not a production runner or origin proof."""

    def __init__(self, root: Path, *, _manifest: dict[str, object]) -> None:
        self.root = _canonical_root(root)
        self._manifest = validate_root_manifest(_manifest)
        self._thread_lock = threading.Lock()
        self._lock_fd: int | None = None
        self._issued: HLT17WriterCapability | None = None
        self._lock_held_here = False
        self._publication_uncertain = False
        self._root_identity = _dir_identity(self.root)
        self._root_manifest_raw = _read_blob(self.root, ROOT_MANIFEST_NAME)
        if self._root_manifest_raw != canonical_json_bytes(self._manifest):
            _fail("root manifest bytes differ from the validated mapping")
        self._lock_raw = _read_blob(self.root, WRITER_EXCLUSION_NAME)
        if self._lock_raw != _exclusion_payload():
            _fail("writer exclusion lock bytes differ")
        self._lock_identity = _lock_path_identity(self.root)

    @property
    def campaign_id(self) -> str:
        return str(self._manifest["campaign_id"])

    @property
    def member_keys(self) -> tuple[str, ...]:
        return tuple(self._manifest["member_keys"])

    @property
    def remaining_seams(self) -> dict[str, str]:
        return {
            "physical_origin": FUTURE_PHYSICAL_ORIGIN_SEAM,
            "source_manifest": FUTURE_SOURCE_MANIFEST_SEAM_TEXT,
            "environment_gate": FUTURE_ENVIRONMENT_GATE_SEAM,
            "resource_gate": FUTURE_RESOURCE_GATE_SEAM,
            "independent_binder": FUTURE_INDEPENDENT_BINDER_SEAM,
            "live_runner": FUTURE_LIVE_RUNNER_SEAM,
        }

    def _reauthenticate_identity(self) -> tuple[bytes, dict[str, object]]:
        """Re-read immutable root/lock bytes and inode identity before use."""

        root = _canonical_root(self.root)
        if root != self.root:
            _fail("store root path raced")
        if _dir_identity(root) != self._root_identity:
            _fail("store root inode was replaced")
        raw, mapping = _read_json(root, ROOT_MANIFEST_NAME)
        try:
            manifest = validate_root_manifest(mapping)
        except HLT17ProtocolError as error:
            raise HLT17CampaignStoreError("root manifest differs") from error
        if raw != self._root_manifest_raw:
            _fail("root manifest bytes drifted")
        if canonical_json_bytes(manifest) != raw:
            _fail("root manifest is not a redundant canonical mapping")
        lock = _read_blob(root, WRITER_EXCLUSION_NAME)
        if lock != self._lock_raw or lock != _exclusion_payload():
            _fail("writer exclusion lock bytes drifted")
        lock_identity = _lock_path_identity(root)
        if lock_identity != self._lock_identity:
            _fail("writer exclusion lock inode was replaced")
        if self._lock_fd is not None:
            try:
                held = os.fstat(self._lock_fd)
            except OSError as error:
                raise HLT17CampaignStoreError("held exclusion lock cannot be read") from error
            if (held.st_dev, held.st_ino) != lock_identity:
                _fail("held exclusion inode is not the current lock file")
            if not stat.S_ISREG(held.st_mode) or held.st_nlink != 1:
                _fail("held exclusion lock is unsafe")
        self._manifest = manifest
        return raw, manifest

    def _refuse_uncertain_mutation(self) -> None:
        if self._publication_uncertain:
            raise HLT17PostpublicationUncertainty(
                "uncertain publication forbids further mutation"
            )

    def _poison_writer(self, message: str) -> HLT17PostpublicationUncertainty:
        self._publication_uncertain = True
        return HLT17PostpublicationUncertainty(message)

    @classmethod
    def create(
        cls,
        root: Path,
        *,
        campaign_id: str,
        member_keys: object,
    ) -> "HLT17CampaignStore":
        root = _canonical_root(root)
        files, sessions, generations = _list_root(root)
        if files or sessions or generations:
            _fail("synthetic store root must be empty")
        manifest = build_root_manifest(campaign_id=campaign_id, member_keys=member_keys)
        try:
            publish_exclusive_file(root, ROOT_MANIFEST_NAME, canonical_json_bytes(manifest))
        except PrepublicationError as error:
            raise HLT17CampaignStoreError("root manifest was not published") from error
        except PostpublicationUncertainty as error:
            raise HLT17PostpublicationUncertainty(
                "root manifest became visible but was not confirmed"
            ) from error
        try:
            publish_exclusive_file(root, WRITER_EXCLUSION_NAME, _exclusion_payload())
        except (PrepublicationError, PostpublicationUncertainty) as error:
            # The manifest is already visible. Failure of the second identity
            # leaf is therefore uncertainty after publication for the store
            # creation operation, even when the lock leaf itself is absent.
            raise HLT17PostpublicationUncertainty(
                "root manifest was published but writer exclusion was not confirmed"
            ) from error
        return cls(root, _manifest=manifest)

    @classmethod
    def open(cls, root: Path) -> "HLT17CampaignStore":
        root = _canonical_root(root)
        raw, mapping = _read_json(root, ROOT_MANIFEST_NAME)
        del raw
        try:
            manifest = validate_root_manifest(mapping)
        except HLT17ProtocolError as error:
            raise HLT17CampaignStoreError("root manifest differs") from error
        lock = _read_blob(root, WRITER_EXCLUSION_NAME)
        if lock != _exclusion_payload():
            _fail("writer exclusion lock bytes differ")
        return cls(root, _manifest=manifest)

    def _session_state(self) -> tuple[int, bool, dict[str, object] | None]:
        _files, sessions, _ignored = _list_root(self.root)
        opened: dict[int, tuple[bytes, dict[str, object]]] = {}
        closed: dict[int, tuple[bytes, dict[str, object]]] = {}
        for name in sessions:
            index, kind = parse_session_name(name)
            raw, mapping = _read_json(self.root, name)
            if kind == "open":
                item = validate_session_open(mapping)
                if item["campaign_id"] != self.campaign_id:
                    _fail("session campaign_id differs")
                if item["exclusion_lock_sha256"] != _digest(_exclusion_payload()):
                    _fail("session does not bind the exclusion lock")
                opened[index] = (raw, item)
            else:
                item = validate_session_close(mapping)
                if item["campaign_id"] != self.campaign_id:
                    _fail("session close campaign_id differs")
                closed[index] = (raw, item)
        if opened and set(opened) != set(range(max(opened) + 1)):
            _fail("writer session sequence has a hole")
        for index, (raw, item) in closed.items():
            if index not in opened:
                _fail("session close without a matching open")
            expected = encode_session_close(
                campaign_id=self.campaign_id,
                session_index=index,
                session_open_sha256=_digest(opened[index][0]),
                owner_token=str(opened[index][1]["owner_token"]),
            )
            if item != expected or raw != canonical_json_bytes(expected):
                _fail("session close does not bind the open session")
        unclosed = [index for index in sorted(opened) if index not in closed]
        if len(unclosed) > 1:
            _fail("multiple unclosed writer sessions")
        next_index = 0 if not opened else max(opened) + 1
        if unclosed:
            return unclosed[0], True, opened[unclosed[0]][1]
        return next_index, False, None

    def _load_generations(self) -> tuple[HLT17ValidatedGeneration, ...]:
        self._reauthenticate_identity()
        files, sessions, directories = _list_root(self.root)
        del sessions
        required = {ROOT_MANIFEST_NAME, WRITER_EXCLUSION_NAME}
        extra = [name for name in files if name not in required | {TERMINAL_LOCK_NAME}]
        if extra:
            _fail(f"unexpected file: {extra[0]}")
        ordered = sorted(directories, key=parse_generation_directory_name)
        expected = [generation_directory_name(index) for index in range(len(ordered))]
        if list(ordered) != expected:
            _fail("generation directories are not a unique complete sequence")
        payloads = [_generation_files(self.root, name) for name in ordered]
        if not payloads:
            return ()
        try:
            return authenticate_generation_chain(root=self._manifest, generations=payloads)
        except HLT17ProtocolError as error:
            raise HLT17CampaignStoreError(str(error)) from error

    def _terminal_lock(self, generation: HLT17ValidatedGeneration | None) -> dict[str, object] | None:
        files, _sessions, _directories = _list_root(self.root)
        if TERMINAL_LOCK_NAME not in files:
            return None
        _raw, mapping = _read_json(self.root, TERMINAL_LOCK_NAME)
        try:
            item = validate_terminal_lock(mapping)
        except HLT17ProtocolError as error:
            raise HLT17CampaignStoreError("terminal lock differs") from error
        if item["campaign_id"] != self.campaign_id:
            _fail("terminal lock campaign_id differs")
        if generation is None:
            _fail("terminal lock without a published generation")
        expected = encode_terminal_lock(
            campaign_id=self.campaign_id,
            checkpoint_sha256=generation.checkpoint_sha256,
            journal_record_sha256=generation.journal_sha256,
            store_generation=generation.store_generation,
            disposition=generation.disposition,
            changed_member_key=generation.changed_member_key,
        )
        if item != expected:
            _fail("terminal lock does not bind the published checkpoint")
        return item

    def authenticate(self) -> HLT17AuthenticatedView:
        """Read-only unique-chain reconstruction. Does not mutate the store."""

        manifest_raw, _manifest = self._reauthenticate_identity()
        generations = self._load_generations()
        last = generations[-1] if generations else None
        lock = self._terminal_lock(last)
        if last is None:
            if lock is not None:
                _fail("terminal lock without a published generation")
            _fail("campaign store has no published generations")
        last.__post_init__()
        if lock is not None and not last.terminal:
            _fail("terminal lock present on a nonterminal checkpoint")
        _index, unclosed, _session = self._session_state()
        del _index
        return HLT17AuthenticatedView(
            root=self.root,
            root_manifest_raw=manifest_raw,
            generations=generations,
            unclosed_session=unclosed,
        )

    @classmethod
    def authenticate_read_only(cls, root: Path) -> HLT17AuthenticatedView:
        return cls.open(root).authenticate()

    def _require_capability(
        self,
        capability: HLT17WriterCapability,
        *,
        allow_uncertain: bool = False,
    ) -> None:
        if capability is not self._issued:
            _fail("writer capability was not issued by this store session")
        if capability.pid != os.getpid() or capability.thread_id != threading.get_ident():
            raise HLT17WriterConflict("writer capability does not own this process/thread")
        if self._lock_fd is None or not self._lock_held_here:
            raise HLT17WriterConflict("writer exclusion is not held")
        self._reauthenticate_identity()
        if not allow_uncertain:
            self._refuse_uncertain_mutation()
        _index, unclosed, session = self._session_state()
        if not unclosed or session is None:
            raise HLT17UnclosedSession("writer session is not open")
        if int(session["session_index"]) != capability.session_index:
            _fail("writer capability session differs")
        if str(session["owner_token"]) != capability.owner_token:
            _fail("writer capability owner differs")
        open_raw, _mapping = _read_json(self.root, session_open_name(capability.session_index))
        if _digest(open_raw) != capability.session_open_sha256:
            _fail("writer session bytes differ")
        del _index

    def _lock_exclusion(self) -> None:
        self._reauthenticate_identity()
        flags = os.O_RDONLY | _flag("O_NOFOLLOW")
        root_fd = _open_root_fd(self.root)
        try:
            try:
                lock_fd = os.open(WRITER_EXCLUSION_NAME, flags, dir_fd=root_fd)
            except OSError as error:
                raise HLT17CampaignStoreError("writer exclusion lock cannot be opened") from error
        finally:
            os.close(root_fd)
        try:
            info = os.fstat(lock_fd)
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise HLT17CampaignStoreError("writer exclusion lock is unsafe")
            current = _lock_path_identity(self.root)
            if (info.st_dev, info.st_ino) != current:
                raise HLT17CampaignStoreError(
                    "held exclusion inode is not the current lock file"
                )
            fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            os.close(lock_fd)
            raise HLT17WriterConflict("writer exclusion is held") from error
        except Exception:
            os.close(lock_fd)
            raise
        self._lock_fd = lock_fd
        self._lock_held_here = True
        self._lock_identity = (info.st_dev, info.st_ino)

    def _unlock_exclusion(self) -> None:
        lock_fd = self._lock_fd
        self._lock_fd = None
        self._lock_held_here = False
        if lock_fd is None:
            return
        try:
            fcntl.flock(lock_fd, fcntl.LOCK_UN)
        finally:
            os.close(lock_fd)

    def acquire_writer(self, *, owner_token: str, host: str) -> HLT17WriterCapability:
        """Publish a durable session after exclusion. Does not take over an unclosed session."""

        self._reauthenticate_identity()
        if not self._thread_lock.acquire(blocking=False):
            raise HLT17WriterConflict("in-process writer exclusion is held")
        try:
            if self._issued is not None:
                raise HLT17WriterConflict("this process already holds the writer")
            self._lock_exclusion()
            try:
                next_index, unclosed, _session = self._session_state()
                if unclosed:
                    raise HLT17UnclosedSession(
                        "a free OS lock is not authority to take over an unclosed session"
                    )
                generations = self._load_generations()
                last = generations[-1] if generations else None
                lock = self._terminal_lock(last)
                if last is not None:
                    last.__post_init__()
                    if last.terminal or last.first_event_complete:
                        raise HLT17TerminalStore("typed terminal forbids further appends")
                if lock is not None:
                    raise HLT17TerminalStore("typed terminal forbids further appends")
                payload = encode_session_open(
                    campaign_id=self.campaign_id,
                    session_index=next_index,
                    owner_token=owner_token,
                    host=host,
                    pid=os.getpid(),
                    thread_id=threading.get_ident(),
                    exclusion_lock_sha256=_digest(_exclusion_payload()),
                )
                raw = canonical_json_bytes(payload)
                publish_exclusive_file(self.root, session_open_name(next_index), raw)
                capability = HLT17WriterCapability(
                    campaign_id=self.campaign_id,
                    session_index=next_index,
                    owner_token=str(payload["owner_token"]),
                    host=str(payload["host"]),
                    pid=int(payload["pid"]),
                    thread_id=int(payload["thread_id"]),
                    session_open_sha256=_digest(raw),
                )
                self._issued = capability
                return capability
            except Exception:
                self._unlock_exclusion()
                raise
        except Exception:
            if self._issued is None:
                self._thread_lock.release()
            raise

    def close_writer(self, capability: HLT17WriterCapability) -> None:
        self._require_capability(capability)
        self._refuse_uncertain_mutation()
        payload = encode_session_close(
            campaign_id=self.campaign_id,
            session_index=capability.session_index,
            session_open_sha256=capability.session_open_sha256,
            owner_token=capability.owner_token,
        )
        try:
            publish_exclusive_file(
                self.root,
                session_close_name(capability.session_index),
                canonical_json_bytes(payload),
            )
        except PrepublicationError as error:
            raise HLT17CampaignStoreError("writer session could not be closed") from error
        except PostpublicationUncertainty as error:
            raise HLT17PostpublicationUncertainty(
                "writer session close became visible but was not confirmed"
            ) from error
        self._issued = None
        self._unlock_exclusion()
        self._thread_lock.release()

    def abandon_writer(self, capability: HLT17WriterCapability) -> None:
        """Release process/thread exclusion without closing durable session ownership.

        This is resource cleanup only.  It does not reconcile, repair, retry, or
        publish session_close.  A later acquire still refuses an unclosed session.
        """

        if capability is not self._issued:
            _fail("writer capability was not issued by this store session")
        if capability.pid != os.getpid() or capability.thread_id != threading.get_ident():
            raise HLT17WriterConflict("writer capability does not own this process/thread")
        if self._lock_fd is None or not self._lock_held_here:
            raise HLT17WriterConflict("writer exclusion is not held")
        self._issued = None
        self._unlock_exclusion()
        self._thread_lock.release()

    def _publish_generation(
        self,
        generation: HLT17ValidatedGeneration,
        *,
        fault_hook=None,
    ) -> HLT17ValidatedGeneration:
        relative = generation.directory_name
        try:
            publish_exclusive_directory(
                self.root,
                relative,
                dict(generation.files),
                _fault_hook=fault_hook,
            )
        except PrepublicationError as error:
            raise HLT17CampaignStoreError(
                "generation was not published"
            ) from error
        except PostpublicationUncertainty as error:
            raise self._poison_writer(
                "generation became visible but a later durability step failed"
            ) from error
        if generation.terminal:
            lock = encode_terminal_lock(
                campaign_id=self.campaign_id,
                checkpoint_sha256=generation.checkpoint_sha256,
                journal_record_sha256=generation.journal_sha256,
                store_generation=generation.store_generation,
                disposition=generation.disposition,
                changed_member_key=generation.changed_member_key,
            )
            try:
                publish_exclusive_file(
                    self.root,
                    TERMINAL_LOCK_NAME,
                    canonical_json_bytes(lock),
                )
            except PrepublicationError as error:
                raise self._poison_writer(
                    "terminal generation was published but terminal.lock was not"
                ) from error
            except PostpublicationUncertainty as error:
                raise self._poison_writer(
                    "terminal lock became visible but was not confirmed"
                ) from error
        return generation

    def publish_seed(
        self,
        capability: HLT17WriterCapability,
        bundles: Mapping[str, HLT17GenerationBundle],
        *,
        _fault_hook=None,
    ) -> HLT17ValidatedGeneration:
        self._require_capability(capability)
        generations = self._load_generations()
        if generations:
            _fail("seed generation already exists")
        try:
            generation = build_seed_generation(root=self._manifest, bundles=bundles)
        except HLT17ProtocolError as error:
            raise HLT17CampaignStoreError(str(error)) from error
        return self._publish_generation(generation, fault_hook=_fault_hook)

    def publish_attempt(
        self,
        capability: HLT17WriterCapability,
        *,
        kind: str,
        member_key: str,
        successor_bundle: HLT17GenerationBundle,
        terminal_evidence: Mapping[str, object] | None = None,
        _fault_hook=None,
    ) -> HLT17ValidatedGeneration:
        self._require_capability(capability)
        generations = self._load_generations()
        if not generations:
            _fail("attempt publication requires a seeded predecessor")
        predecessor = generations[-1]
        predecessor.__post_init__()
        lock = self._terminal_lock(predecessor)
        if predecessor.terminal or predecessor.first_event_complete or lock is not None:
            raise HLT17TerminalStore("typed terminal forbids further appends")
        try:
            generation = build_successor_generation(
                root=self._manifest,
                predecessor=predecessor,
                predecessor_bundles=predecessor.bundles,
                kind=kind,
                member_key=member_key,
                successor_bundle=successor_bundle,
                terminal_evidence=terminal_evidence,
            )
        except HLT17ProtocolError as error:
            raise HLT17CampaignStoreError(str(error)) from error
        return self._publish_generation(generation, fault_hook=_fault_hook)


__all__ = [
    "EVENT_ORIGIN_HEX",
    "EVENT_TARGET_HEX",
    "HLT17AuthenticatedView",
    "HLT17CampaignStore",
    "HLT17CampaignStoreError",
    "HLT17PostpublicationUncertainty",
    "HLT17TerminalStore",
    "HLT17UnclosedSession",
    "HLT17WriterCapability",
    "HLT17WriterConflict",
    "PRODUCTION_EVENT_PATH",
    "ROOT_MANIFEST_NAME",
    "TERMINAL_LOCK_NAME",
    "WRITER_EXCLUSION_NAME",
]
