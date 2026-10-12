"""Atomic prospective RSRC1 fresh first-event publisher.

This store is a new namespace/schema sibling, not a switch in or append path
to the sealed PRO20 store. It consumes validated live C1R1 member bundles and an externally
supplied immutable authority receipt, then publishes an entire seed store
atomically to the production namespace shape.  Generation-level authority and
claim booleans remain false.  A free OS lock is not authority to take over an
unclosed session.  Reads authenticate the unique chain without mutation.
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
from .pro20_rsrc1_protocol import (
    EXACT_DELTA_AUTHORITY_SEAM,
    FREEZE_ARTIFACT_ID,
    INDEPENDENT_BINDER_SEAM,
    PRODUCTION_NAMESPACE,
    PROTOCOL_ARTIFACT_ID,
    PRO20RSRC1AuthorityReceipt,
    PRO20RSRC1ProtocolError,
    PRO20RSRC1ValidatedGeneration,
    RESOURCE_GATE_SEAM,
    RUNNER_SEAM,
    authenticate_generation_chain,
    authenticate_generation_chain_latest,
    authority_receipt_from_root,
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
WRITER_EXCLUSION_SCHEMA = "FGC-1-PRO20-EV1-RSRC1-writer-exclusion-v1"
ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1"
PRODUCTION_CONTAINER = "runs/fgc-2-sf1/pro20-rsrc1-event1"
PRODUCTION_STORE_LEAF = "event"

if PRODUCTION_NAMESPACE != f"{PRODUCTION_CONTAINER}/{PRODUCTION_STORE_LEAF}":
    raise RuntimeError("PRO20 RSRC1 production namespace decomposition differs")


class PRO20RSRC1CampaignStoreError(ValueError):
    """A PRO20-EV1 campaign store is unsafe, incomplete, forked, or unauthorized."""

    published = False


class PRO20RSRC1WriterConflict(PRO20RSRC1CampaignStoreError):
    """Process or thread exclusion refused a second writer."""


class PRO20RSRC1UnclosedSession(PRO20RSRC1CampaignStoreError):
    """A free OS lock does not authorize takeover of an unclosed session."""


class PRO20RSRC1TerminalStore(PRO20RSRC1CampaignStoreError):
    """A typed or successful-event terminal preserves history and forbids appends."""


class PRO20RSRC1PostpublicationUncertainty(PRO20RSRC1CampaignStoreError):
    """A generation or seed store became visible, but a later durability step failed.

    Complete visible bytes may be readable as evidence.  They are not
    permission to resume an unclosed or uncertain writer.
    """

    published = True


def _fail(message: str) -> None:
    raise PRO20RSRC1CampaignStoreError(message)


def _digest(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _canonical_repository(root: Path) -> Path:
    if not isinstance(root, Path):
        raise TypeError("repository root must be a pathlib.Path")
    try:
        supplied = Path(os.fspath(root))
    except TypeError as error:
        raise PRO20RSRC1CampaignStoreError("repository root is unsafe") from error
    if not supplied.is_absolute() or any(part in {".", ".."} for part in supplied.parts):
        _fail("repository root must be an absolute canonical path")
    try:
        resolved = supplied.resolve(strict=True)
    except OSError as error:
        raise PRO20RSRC1CampaignStoreError("repository root cannot be resolved safely") from error
    if resolved != supplied:
        _fail("repository root path traverses a symlinked ancestor")
    return resolved


def _namespace_root(repository: Path) -> Path:
    repo = _canonical_repository(repository)
    store = repo.joinpath(*PRODUCTION_NAMESPACE.split("/"))
    if store.as_posix() != f"{repo.as_posix()}/{PRODUCTION_NAMESPACE}":
        _fail("production namespace path differs")
    return store


def _exclusion_payload(receipt: PRO20RSRC1AuthorityReceipt) -> bytes:
    receipt.__post_init__()
    return canonical_json_bytes(
        {
            "schema": WRITER_EXCLUSION_SCHEMA,
            "schema_version": 1,
            "protocol_artifact_id": PROTOCOL_ARTIFACT_ID,
            "freeze_artifact_id": FREEZE_ARTIFACT_ID,
            "role": "process_thread_exclusion",
            "namespace": PRODUCTION_NAMESPACE,
            "authority_receipt_sha256": receipt.sha256,
            "production_write_authorized": False,
            "campaign_execution_authorized": False,
            "automatic_stale_writer_recovery": False,
        }
    )


def _flag(name: str) -> int:
    value = int(getattr(os, name, 0) or 0)
    if value == 0:
        raise PRO20RSRC1CampaignStoreError(f"{name} is required")
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
        except PRO20RSRC1ProtocolError as error:
            raise PRO20RSRC1CampaignStoreError(f"unexpected file: {name}") from error
        return "session"
    if not stat.S_ISDIR(info.st_mode):
        _fail(f"unsafe node: {name}")
    try:
        parse_generation_directory_name(name)
    except PRO20RSRC1ProtocolError as error:
        raise PRO20RSRC1CampaignStoreError(f"unexpected directory: {name}") from error
    return "generation"


def _open_root_fd(root: Path) -> int:
    flags = os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW")
    try:
        before = root.lstat()
    except OSError as error:
        raise PRO20RSRC1CampaignStoreError("store root absent or unsafe") from error
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
                raise PRO20RSRC1CampaignStoreError(f"node raced: {item.name}") from error
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
        raise PRO20RSRC1CampaignStoreError(f"{relative} cannot be read safely") from error
    if type(value) is not dict:
        _fail(f"{relative} is not a JSON object")
    return raw, value


def _read_blob(root: Path, relative: str) -> bytes:
    try:
        return read_regular_file(root, relative, max_bytes=DEFAULT_MAX_BYTES)
    except (UnsafePathError, TypeError, ValueError) as error:
        raise PRO20RSRC1CampaignStoreError(f"{relative} cannot be read safely") from error


def _generation_files(root: Path, directory: str) -> dict[str, bytes]:
    descriptor = _open_root_fd(root)
    files: dict[str, bytes] = {}
    child = None
    try:
        try:
            before = os.stat(directory, dir_fd=descriptor, follow_symlinks=False)
        except OSError as error:
            raise PRO20RSRC1CampaignStoreError(f"generation absent: {directory}") from error
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
class PRO20RSRC1WriterCapability:
    """In-process handle issued only after session publication and exclusion."""

    campaign_id: str
    session_index: int
    owner_token: str
    host: str
    pid: int
    thread_id: int
    session_open_sha256: str
    authority_receipt_sha256: str
    namespace: str


@dataclass(frozen=True, slots=True)
class PRO20RSRC1AuthenticatedView:
    """Read-only reconstruction of a unique complete published production chain."""

    repository_root: Path
    store_root: Path
    root_manifest_raw: bytes
    generations: tuple[PRO20RSRC1ValidatedGeneration, ...]
    unclosed_session: bool

    def __post_init__(self) -> None:
        if type(self.root_manifest_raw) is not bytes:
            raise TypeError("root manifest raw must be bytes")
        object.__setattr__(self, "generations", tuple(self.generations))
        try:
            validate_root_manifest(load_canonical_json(self.root_manifest_raw))
        except (CanonicalJSONError, PRO20RSRC1ProtocolError, TypeError, ValueError) as error:
            raise PRO20RSRC1CampaignStoreError("root manifest differs") from error
        if not self.generations:
            raise PRO20RSRC1CampaignStoreError("campaign store has no published generations")
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
        raise PRO20RSRC1CampaignStoreError("store root absent or unsafe") from error
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        _fail("store root is unsafe")
    return (info.st_dev, info.st_ino)


def _lock_path_identity(root: Path) -> tuple[int, int]:
    root_fd = _open_root_fd(root)
    try:
        try:
            info = os.stat(WRITER_EXCLUSION_NAME, dir_fd=root_fd, follow_symlinks=False)
        except OSError as error:
            raise PRO20RSRC1CampaignStoreError("writer exclusion lock cannot be read") from error
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            _fail("writer exclusion lock is unsafe")
        return (info.st_dev, info.st_ino)
    finally:
        os.close(root_fd)


class PRO20RSRC1CampaignStore:
    """Production first-event publisher. Not a runner, freeze, or origin proof."""

    def __init__(
        self,
        repository_root: Path,
        *,
        _manifest: dict[str, object],
        _receipt: PRO20RSRC1AuthorityReceipt,
    ) -> None:
        self.repository_root = _canonical_repository(repository_root)
        self.root = _namespace_root(self.repository_root)
        self._manifest = validate_root_manifest(_manifest)
        self._receipt = _receipt
        self._receipt.__post_init__()
        if self._receipt.sha256 != self._manifest["authority_receipt_sha256"]:
            _fail("store receipt digest disagrees with the root manifest")
        self._thread_lock = threading.Lock()
        self._lock_fd: int | None = None
        self._issued: PRO20RSRC1WriterCapability | None = None
        self._lock_held_here = False
        self._publication_uncertain = False
        self._root_identity = _dir_identity(self.root)
        self._root_manifest_raw = _read_blob(self.root, ROOT_MANIFEST_NAME)
        if self._root_manifest_raw != canonical_json_bytes(self._manifest):
            _fail("root manifest bytes differ from the validated mapping")
        self._lock_raw = _read_blob(self.root, WRITER_EXCLUSION_NAME)
        if self._lock_raw != _exclusion_payload(self._receipt):
            _fail("writer exclusion lock bytes differ")
        self._lock_identity = _lock_path_identity(self.root)

    @property
    def campaign_id(self) -> str:
        return str(self._manifest["campaign_id"])

    @property
    def member_keys(self) -> tuple[str, ...]:
        return tuple(self._manifest["member_keys"])

    @property
    def authority_receipt_sha256(self) -> str:
        return self._receipt.sha256

    @property
    def remaining_seams(self) -> dict[str, str]:
        return {
            "runner": RUNNER_SEAM,
            "exact_delta_authority": EXACT_DELTA_AUTHORITY_SEAM,
            "resource_gate": RESOURCE_GATE_SEAM,
            "independent_binder": INDEPENDENT_BINDER_SEAM,
        }

    def _reauthenticate_identity(self) -> tuple[bytes, dict[str, object]]:
        repository = _canonical_repository(self.repository_root)
        if repository != self.repository_root:
            _fail("repository root path raced")
        root = _namespace_root(repository)
        if root != self.root:
            _fail("store root path raced")
        if _dir_identity(root) != self._root_identity:
            _fail("store root inode was replaced")
        raw, mapping = _read_json(root, ROOT_MANIFEST_NAME)
        try:
            manifest = validate_root_manifest(mapping)
        except PRO20RSRC1ProtocolError as error:
            raise PRO20RSRC1CampaignStoreError("root manifest differs") from error
        if raw != self._root_manifest_raw:
            _fail("root manifest bytes drifted")
        if canonical_json_bytes(manifest) != raw:
            _fail("root manifest is not a redundant canonical mapping")
        lock = _read_blob(root, WRITER_EXCLUSION_NAME)
        if lock != self._lock_raw or lock != _exclusion_payload(self._receipt):
            _fail("writer exclusion lock bytes drifted")
        lock_identity = _lock_path_identity(root)
        if lock_identity != self._lock_identity:
            _fail("writer exclusion lock inode was replaced")
        if self._lock_fd is not None:
            try:
                held = os.fstat(self._lock_fd)
            except OSError as error:
                raise PRO20RSRC1CampaignStoreError("held exclusion lock cannot be read") from error
            if (held.st_dev, held.st_ino) != lock_identity:
                _fail("held exclusion inode is not the current lock file")
            if not stat.S_ISREG(held.st_mode) or held.st_nlink != 1:
                _fail("held exclusion lock is unsafe")
        self._manifest = manifest
        return raw, manifest

    def _refuse_uncertain_mutation(self) -> None:
        if self._publication_uncertain:
            raise PRO20RSRC1PostpublicationUncertainty(
                "uncertain publication forbids further mutation"
            )

    def _poison_writer(self, message: str) -> PRO20RSRC1PostpublicationUncertainty:
        self._publication_uncertain = True
        return PRO20RSRC1PostpublicationUncertainty(message)

    @classmethod
    def publish_seed_store(
        cls,
        repository_root: Path,
        *,
        receipt: object,
        bundles: Mapping[str, HLT17GenerationBundle],
        _fault_hook=None,
    ) -> "PRO20RSRC1CampaignStore":
        """Publish the entire seed store atomically to the production namespace."""

        repository = _canonical_repository(repository_root)
        _namespace_root(repository)
        container = repository.joinpath(*PRODUCTION_CONTAINER.split("/"))
        parent = container.parent
        try:
            parent_info = parent.lstat()
        except OSError as error:
            raise PRO20RSRC1CampaignStoreError(
                "production container parent is absent"
            ) from error
        if stat.S_ISLNK(parent_info.st_mode) or not stat.S_ISDIR(parent_info.st_mode):
            _fail("production container parent is unsafe")
        if not isinstance(receipt, PRO20RSRC1AuthorityReceipt):
            receipt = PRO20RSRC1AuthorityReceipt.from_mapping(receipt)
        receipt.__post_init__()
        try:
            manifest = build_root_manifest(receipt)
            generation = build_seed_generation(root=manifest, bundles=bundles)
        except PRO20RSRC1ProtocolError as error:
            raise PRO20RSRC1CampaignStoreError(str(error)) from error
        files: dict[str, bytes] = {
            f"{PRODUCTION_STORE_LEAF}/{ROOT_MANIFEST_NAME}": canonical_json_bytes(manifest),
            f"{PRODUCTION_STORE_LEAF}/{WRITER_EXCLUSION_NAME}": _exclusion_payload(receipt),
        }
        prefix = f"{PRODUCTION_STORE_LEAF}/{generation.directory_name}"
        for relative, payload in generation.files.items():
            files[f"{prefix}/{relative}"] = payload
        try:
            publish_exclusive_directory(
                repository,
                PRODUCTION_CONTAINER,
                files,
                _fault_hook=_fault_hook,
            )
        except PrepublicationError as error:
            raise PRO20RSRC1CampaignStoreError("seed store was not published") from error
        except PostpublicationUncertainty as error:
            raise PRO20RSRC1PostpublicationUncertainty(
                "seed store became visible but a later durability step failed"
            ) from error
        return cls(repository, _manifest=manifest, _receipt=receipt)

    @classmethod
    def open(cls, repository_root: Path) -> "PRO20RSRC1CampaignStore":
        repository = _canonical_repository(repository_root)
        root = _namespace_root(repository)
        raw, mapping = _read_json(root, ROOT_MANIFEST_NAME)
        del raw
        try:
            manifest = validate_root_manifest(mapping)
            receipt = authority_receipt_from_root(manifest)
        except PRO20RSRC1ProtocolError as error:
            raise PRO20RSRC1CampaignStoreError("root manifest differs") from error
        if receipt.sha256 != manifest["authority_receipt_sha256"]:
            _fail("root authority receipt digest differs")
        lock = _read_blob(root, WRITER_EXCLUSION_NAME)
        if lock != _exclusion_payload(receipt):
            _fail("writer exclusion lock bytes differ")
        return cls(repository, _manifest=manifest, _receipt=receipt)

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
                if item["authority_receipt_sha256"] != self._receipt.sha256:
                    _fail("session does not bind the authority receipt")
                if item["exclusion_lock_sha256"] != _digest(_exclusion_payload(self._receipt)):
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
                authority_receipt_sha256=self._receipt.sha256,
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

    def _ordered_generation_directories(self) -> tuple[str, ...]:
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
        return tuple(ordered)

    def _load_generations(self) -> tuple[PRO20RSRC1ValidatedGeneration, ...]:
        ordered = self._ordered_generation_directories()
        payloads = [_generation_files(self.root, name) for name in ordered]
        if not payloads:
            return ()
        try:
            return authenticate_generation_chain(root=self._manifest, generations=payloads)
        except PRO20RSRC1ProtocolError as error:
            raise PRO20RSRC1CampaignStoreError(str(error)) from error

    def _load_latest_generation(self) -> PRO20RSRC1ValidatedGeneration:
        """Authenticate the complete chain with bounded live history state."""

        ordered = self._ordered_generation_directories()
        try:
            return authenticate_generation_chain_latest(
                root=self._manifest,
                generations=(_generation_files(self.root, name) for name in ordered),
            )
        except PRO20RSRC1ProtocolError as error:
            raise PRO20RSRC1CampaignStoreError(str(error)) from error

    def _terminal_lock(
        self, generation: PRO20RSRC1ValidatedGeneration | None
    ) -> dict[str, object] | None:
        files, _sessions, _directories = _list_root(self.root)
        if TERMINAL_LOCK_NAME not in files:
            return None
        _raw, mapping = _read_json(self.root, TERMINAL_LOCK_NAME)
        try:
            item = validate_terminal_lock(mapping)
        except PRO20RSRC1ProtocolError as error:
            raise PRO20RSRC1CampaignStoreError("terminal lock differs") from error
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
            first_event_complete=generation.first_event_complete,
        )
        if item != expected:
            _fail("terminal lock does not bind the published checkpoint")
        return item

    def authenticate(self) -> PRO20RSRC1AuthenticatedView:
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
        return PRO20RSRC1AuthenticatedView(
            repository_root=self.repository_root,
            store_root=self.root,
            root_manifest_raw=manifest_raw,
            generations=generations,
            unclosed_session=unclosed,
        )

    @classmethod
    def authenticate_read_only(cls, repository_root: Path) -> PRO20RSRC1AuthenticatedView:
        return cls.open(repository_root).authenticate()

    def authenticate_latest(self) -> PRO20RSRC1ValidatedGeneration:
        """Read-only complete-chain authentication retaining only the latest item."""

        last = self._load_latest_generation()
        lock = self._terminal_lock(last)
        last.__post_init__()
        if lock is not None and not last.terminal:
            _fail("terminal lock present on a nonterminal checkpoint")
        return last

    def _require_capability(
        self,
        capability: PRO20RSRC1WriterCapability,
        *,
        allow_uncertain: bool = False,
    ) -> None:
        if capability is not self._issued:
            _fail("writer capability was not issued by this store session")
        if capability.pid != os.getpid() or capability.thread_id != threading.get_ident():
            raise PRO20RSRC1WriterConflict("writer capability does not own this process/thread")
        if capability.authority_receipt_sha256 != self._receipt.sha256:
            _fail("writer capability receipt differs")
        if self._lock_fd is None or not self._lock_held_here:
            raise PRO20RSRC1WriterConflict("writer exclusion is not held")
        self._reauthenticate_identity()
        if not allow_uncertain:
            self._refuse_uncertain_mutation()
        _index, unclosed, session = self._session_state()
        if not unclosed or session is None:
            raise PRO20RSRC1UnclosedSession("writer session is not open")
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
                raise PRO20RSRC1CampaignStoreError(
                    "writer exclusion lock cannot be opened"
                ) from error
        finally:
            os.close(root_fd)
        try:
            info = os.fstat(lock_fd)
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise PRO20RSRC1CampaignStoreError("writer exclusion lock is unsafe")
            current = _lock_path_identity(self.root)
            if (info.st_dev, info.st_ino) != current:
                raise PRO20RSRC1CampaignStoreError(
                    "held exclusion inode is not the current lock file"
                )
            fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            os.close(lock_fd)
            raise PRO20RSRC1WriterConflict("writer exclusion is held") from error
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

    def acquire_writer(self, *, owner_token: str, host: str) -> PRO20RSRC1WriterCapability:
        """Publish a durable session after exclusion. Does not take over an unclosed session."""

        self._reauthenticate_identity()
        if not self._thread_lock.acquire(blocking=False):
            raise PRO20RSRC1WriterConflict("in-process writer exclusion is held")
        try:
            if self._issued is not None:
                raise PRO20RSRC1WriterConflict("this process already holds the writer")
            self._lock_exclusion()
            try:
                next_index, unclosed, _session = self._session_state()
                if unclosed:
                    raise PRO20RSRC1UnclosedSession(
                        "a free OS lock is not authority to take over an unclosed session"
                    )
                last = self._load_latest_generation()
                lock = self._terminal_lock(last)
                last.__post_init__()
                if last.terminal or last.first_event_complete:
                    raise PRO20RSRC1TerminalStore("typed terminal forbids further appends")
                if lock is not None:
                    raise PRO20RSRC1TerminalStore("typed terminal forbids further appends")
                payload = encode_session_open(
                    campaign_id=self.campaign_id,
                    session_index=next_index,
                    owner_token=owner_token,
                    host=host,
                    pid=os.getpid(),
                    thread_id=threading.get_ident(),
                    exclusion_lock_sha256=_digest(_exclusion_payload(self._receipt)),
                    authority_receipt_sha256=self._receipt.sha256,
                )
                raw = canonical_json_bytes(payload)
                publish_exclusive_file(self.root, session_open_name(next_index), raw)
                capability = PRO20RSRC1WriterCapability(
                    campaign_id=self.campaign_id,
                    session_index=next_index,
                    owner_token=str(payload["owner_token"]),
                    host=str(payload["host"]),
                    pid=int(payload["pid"]),
                    thread_id=int(payload["thread_id"]),
                    session_open_sha256=_digest(raw),
                    authority_receipt_sha256=self._receipt.sha256,
                    namespace=PRODUCTION_NAMESPACE,
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

    def close_writer(self, capability: PRO20RSRC1WriterCapability) -> None:
        self._require_capability(capability)
        self._refuse_uncertain_mutation()
        payload = encode_session_close(
            campaign_id=self.campaign_id,
            session_index=capability.session_index,
            session_open_sha256=capability.session_open_sha256,
            owner_token=capability.owner_token,
            authority_receipt_sha256=capability.authority_receipt_sha256,
        )
        try:
            publish_exclusive_file(
                self.root,
                session_close_name(capability.session_index),
                canonical_json_bytes(payload),
            )
        except PrepublicationError as error:
            raise PRO20RSRC1CampaignStoreError("writer session could not be closed") from error
        except PostpublicationUncertainty as error:
            raise PRO20RSRC1PostpublicationUncertainty(
                "writer session close became visible but was not confirmed"
            ) from error
        self._issued = None
        self._unlock_exclusion()
        self._thread_lock.release()

    def abandon_writer(self, capability: PRO20RSRC1WriterCapability) -> None:
        """Release process/thread exclusion without closing durable session ownership.

        This is resource cleanup only.  It does not reconcile, repair, retry, or
        publish session_close.  A later acquire still refuses an unclosed session.
        """

        if capability is not self._issued:
            _fail("writer capability was not issued by this store session")
        if capability.pid != os.getpid() or capability.thread_id != threading.get_ident():
            raise PRO20RSRC1WriterConflict("writer capability does not own this process/thread")
        if self._lock_fd is None or not self._lock_held_here:
            raise PRO20RSRC1WriterConflict("writer exclusion is not held")
        self._issued = None
        self._unlock_exclusion()
        self._thread_lock.release()

    def _publish_generation(
        self,
        generation: PRO20RSRC1ValidatedGeneration,
        *,
        fault_hook=None,
    ) -> PRO20RSRC1ValidatedGeneration:
        relative = generation.directory_name
        try:
            publish_exclusive_directory(
                self.root,
                relative,
                dict(generation.files),
                _fault_hook=fault_hook,
            )
        except PrepublicationError as error:
            raise PRO20RSRC1CampaignStoreError("generation was not published") from error
        except PostpublicationUncertainty as error:
            raise self._poison_writer(
                "generation became visible but a later durability step failed"
            ) from error
        if generation.terminal or generation.first_event_complete:
            lock = encode_terminal_lock(
                campaign_id=self.campaign_id,
                checkpoint_sha256=generation.checkpoint_sha256,
                journal_record_sha256=generation.journal_sha256,
                store_generation=generation.store_generation,
                disposition=generation.disposition,
                changed_member_key=generation.changed_member_key,
                first_event_complete=generation.first_event_complete,
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

    def publish_attempt(
        self,
        capability: PRO20RSRC1WriterCapability,
        *,
        kind: str,
        member_key: str,
        successor_bundle: HLT17GenerationBundle,
        terminal_evidence: Mapping[str, object] | None = None,
        _fault_hook=None,
    ) -> PRO20RSRC1ValidatedGeneration:
        self._require_capability(capability)
        predecessor = self._load_latest_generation()
        predecessor.__post_init__()
        lock = self._terminal_lock(predecessor)
        if predecessor.terminal or predecessor.first_event_complete or lock is not None:
            raise PRO20RSRC1TerminalStore("typed terminal forbids further appends")
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
        except PRO20RSRC1ProtocolError as error:
            raise PRO20RSRC1CampaignStoreError(str(error)) from error
        return self._publish_generation(generation, fault_hook=_fault_hook)


__all__ = [
    "PRODUCTION_NAMESPACE",
    "PRO20RSRC1AuthenticatedView",
    "PRO20RSRC1CampaignStore",
    "PRO20RSRC1CampaignStoreError",
    "PRO20RSRC1PostpublicationUncertainty",
    "PRO20RSRC1TerminalStore",
    "PRO20RSRC1UnclosedSession",
    "PRO20RSRC1WriterCapability",
    "PRO20RSRC1WriterConflict",
    "PRODUCTION_CONTAINER",
    "ROOT_MANIFEST_NAME",
    "TERMINAL_LOCK_NAME",
    "WRITER_EXCLUSION_NAME",
]
