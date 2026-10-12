"""Content-addressed finite-state persistence for HLT16.

There is one production descriptor schema and no journal/checkpoint authority
in this module.  Campaign ordering belongs to ``hlt16_campaign_store``.
Generation-zero import requires an authority bound to the sealed PREF27
receipt/checkpoint/member identity; later states must never use that route.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
import stat
from typing import Any, Callable, Iterable, Mapping
import zipfile

import numpy as np

from .hlt16_member_codec import (
    ARRAY_NAMES,
    GEN0_PROTOCOL,
    HLT16MemberSnapshot,
    validate_codec_metadata,
)


STATE_SCHEMA = "FGC-1-HLT16-member-state-v1"
GEN0_AUTHORITY_SCHEMA = "FGC-1-HLT16-authenticated-gen0-import-v1"
SEALED_GEN0_RECEIPT_SHA256 = (
    "78ae89a9f74193658003cf51f8567c4088db80fd6fb86867c2d581a0c16d5baf"
)
SEALED_GEN0_CHECKPOINT_SHA256 = (
    "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d"
)
MAX_RAW_ARCHIVE_BYTES = 64 * 1024 * 1024
MAX_ARRAY_MEMBER_BYTES = 32 * 1024 * 1024
MAX_TOTAL_UNCOMPRESSED_BYTES = 96 * 1024 * 1024
FaultHook = Callable[[str], None]


class HLT16StateStoreError(ValueError):
    """A content-addressed state object is unsafe, malformed, or inconsistent."""


def canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as error:
        raise HLT16StateStoreError("object is not canonical JSON") from error


def digest(value: object) -> str:
    return sha256(canonical(value)).hexdigest()


def _sha(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise HLT16StateStoreError(f"{label} is not SHA-256")
    try:
        int(value, 16)
    except ValueError as error:
        raise HLT16StateStoreError(f"{label} is not SHA-256") from error
    return value


def _safe_relative(value: str) -> tuple[str, ...]:
    path = Path(value)
    if (
        path.is_absolute()
        or "\\" in value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise HLT16StateStoreError("unsafe relative path")
    return path.parts


def _decode_json(raw: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(items: Iterable[tuple[str, object]]) -> dict[str, object]:
        answer: dict[str, object] = {}
        for key, value in items:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer

    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise HLT16StateStoreError(f"{label} is malformed") from error
    if not isinstance(value, dict) or canonical(value) != raw:
        raise HLT16StateStoreError(f"{label} is noncanonical")
    return value


def read_nofollow(root: Path, relative: str, label: str) -> bytes:
    """Read a single-linked regular leaf through an all-nofollow parent chain."""
    parts = _safe_relative(relative)
    directory_flags = (
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    )
    leaf_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    opened: list[int] = []
    leaf: int | None = None
    try:
        root_info = root.lstat()
        if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
            raise HLT16StateStoreError(f"{label} root is unsafe")
        parent = os.open(root, directory_flags)
        opened.append(parent)
        for component in parts[:-1]:
            before = os.stat(component, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise HLT16StateStoreError(f"{label} parent is unsafe")
            child = os.open(component, directory_flags, dir_fd=parent)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                raise HLT16StateStoreError(f"{label} parent raced")
            opened.append(child)
            parent = child
        before = os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
        if (
            stat.S_ISLNK(before.st_mode)
            or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
        ):
            raise HLT16StateStoreError(f"{label} leaf is unsafe")
        leaf = os.open(parts[-1], leaf_flags, dir_fd=parent)
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
            raise HLT16StateStoreError(f"{label} leaf raced")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(leaf, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        raw = b"".join(chunks)
        after = os.fstat(leaf)
        if len(raw) != before.st_size or identity != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise HLT16StateStoreError(f"{label} leaf changed during read")
        return raw
    except OSError as error:
        raise HLT16StateStoreError(f"{label} cannot be read safely") from error
    finally:
        if leaf is not None:
            os.close(leaf)
        for descriptor in reversed(opened):
            os.close(descriptor)


def _read_leaf_fd(parent: int, name: str, label: str, *, links: tuple[int, ...] = (1,)) -> bytes:
    """Read one FD-relative regular leaf while proving its identity is stable."""
    before = os.stat(name, dir_fd=parent, follow_symlinks=False)
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink not in links
    ):
        raise HLT16StateStoreError(f"{label} leaf is unsafe")
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=parent,
    )
    try:
        active = os.fstat(descriptor)
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
            raise HLT16StateStoreError(f"{label} leaf raced")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        raw = b"".join(chunks)
        after = os.fstat(descriptor)
        if len(raw) != before.st_size or identity != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise HLT16StateStoreError(f"{label} leaf changed during read")
        return raw
    finally:
        os.close(descriptor)


def _fsync_leaf_fd(parent: int, name: str, label: str) -> None:
    """Make an already validated staging inode durable before linking it."""
    before = os.stat(name, dir_fd=parent, follow_symlinks=False)
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
        raise HLT16StateStoreError(f"{label} leaf is unsafe")
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=parent,
    )
    try:
        active = os.fstat(descriptor)
        if (before.st_dev, before.st_ino) != (active.st_dev, active.st_ino):
            raise HLT16StateStoreError(f"{label} leaf raced")
        os.fsync(descriptor)
        after = os.fstat(descriptor)
        if (active.st_dev, active.st_ino, active.st_size) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
        ):
            raise HLT16StateStoreError(f"{label} leaf changed during fsync")
    finally:
        os.close(descriptor)


def _open_parent_nofollow(root: Path, parts: tuple[str, ...], label: str) -> int:
    """Open an existing parent through an FD-relative all-no-follow chain."""
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    parent = os.open(root, flags)
    try:
        for component in parts:
            before = os.stat(component, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise HLT16StateStoreError(f"{label} parent is unsafe")
            child = os.open(component, flags, dir_fd=parent)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                os.close(child)
                raise HLT16StateStoreError(f"{label} parent raced")
            os.close(parent)
            parent = child
        return parent
    except Exception:
        os.close(parent)
        raise


def atomic_publish_noreplace(
    root: Path,
    relative: str,
    payload: bytes,
    *,
    fault_hook: FaultHook | None = None,
) -> None:
    """Publish exact bytes without replacement and recover our deterministic stage.

    A process kill may leave the narrowly named stage before or after the
    no-replace hard link.  Repeating the exact publication repairs only that
    stage, verifies any linked final inode, and resumes deterministically.
    """
    parts = _safe_relative(relative)
    if not isinstance(payload, bytes):
        raise HLT16StateStoreError("publication payload is not bytes")
    parent = _open_parent_nofollow(Path(root), parts[:-1], "publication")
    final = parts[-1]
    stage = f".{final}.hlt16-stage-{sha256(payload).hexdigest()}"

    def cut(name: str) -> None:
        if fault_hook is not None:
            fault_hook(name)

    try:
        try:
            staged = os.stat(stage, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            staged = None
        if staged is not None:
            if stat.S_ISLNK(staged.st_mode) or not stat.S_ISREG(staged.st_mode) or staged.st_nlink not in {1, 2}:
                raise HLT16StateStoreError("publication stage is unsafe")
            staged_raw = _read_leaf_fd(parent, stage, "publication stage", links=(1, 2))
            try:
                final_info = os.stat(final, dir_fd=parent, follow_symlinks=False)
            except FileNotFoundError:
                final_info = None
            if staged.st_nlink == 2:
                if (
                    final_info is None
                    or not stat.S_ISREG(final_info.st_mode)
                    or (staged.st_dev, staged.st_ino) != (final_info.st_dev, final_info.st_ino)
                    or staged_raw != payload
                ):
                    raise HLT16StateStoreError("linked publication stage differs")
                os.unlink(stage, dir_fd=parent)
                os.fsync(parent)
                staged = None
            elif staged_raw != payload:
                # Only this exact task-owned deterministic staging name is
                # recoverably discarded after an interrupted partial write.
                os.unlink(stage, dir_fd=parent)
                os.fsync(parent)
                staged = None

        try:
            final_info = os.stat(final, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            final_info = None
        if final_info is not None:
            if (
                stat.S_ISLNK(final_info.st_mode)
                or not stat.S_ISREG(final_info.st_mode)
                or final_info.st_nlink != 1
                or _read_leaf_fd(parent, final, "existing publication") != payload
            ):
                raise HLT16StateStoreError("publication collision")
            if staged is not None:
                os.unlink(stage, dir_fd=parent)
                os.fsync(parent)
            return

        if staged is None:
            descriptor = os.open(
                stage,
                os.O_WRONLY
                | os.O_CREAT
                | os.O_EXCL
                | getattr(os, "O_NOFOLLOW", 0),
                0o600,
                dir_fd=parent,
            )
            try:
                cut("before_write")
                view = memoryview(payload)
                while view:
                    count = os.write(descriptor, view)
                    if count <= 0:
                        raise HLT16StateStoreError("short publication write")
                    view = view[count:]
                cut("after_write")
                os.fsync(descriptor)
                cut("after_file_fsync")
            finally:
                os.close(descriptor)
        else:
            # A previous process may have died after completing the write but
            # before fsyncing it.  Exact bytes in page cache are not durability
            # evidence, so an adopted one-link stage is fsynced again here.
            _fsync_leaf_fd(parent, stage, "publication stage")
        os.link(
            stage,
            final,
            src_dir_fd=parent,
            dst_dir_fd=parent,
            follow_symlinks=False,
        )
        cut("after_link")
        os.fsync(parent)
        cut("after_parent_fsync")
        try:
            os.unlink(stage, dir_fd=parent)
        except FileNotFoundError:
            # A concurrent read-only recovery may have removed the proven
            # two-link publication stage after the final inode became durable.
            pass
        os.fsync(parent)
        cut("after_stage_cleanup")
    except FileExistsError as error:
        raise HLT16StateStoreError("publication collision") from error
    except OSError as error:
        raise HLT16StateStoreError("publication failed safely") from error
    finally:
        os.close(parent)


def _manifest(arrays: Mapping[str, np.ndarray]) -> tuple[dict[str, Any], ...]:
    if not isinstance(arrays, Mapping) or tuple(arrays) != ARRAY_NAMES:
        raise HLT16StateStoreError("array inventory/order differs")
    answer: list[dict[str, Any]] = []
    for name in ARRAY_NAMES:
        value = arrays[name]
        if (
            not isinstance(value, np.ndarray)
            or value.dtype != np.dtype("<f8")
            or not value.flags.c_contiguous
            or not value.shape
            or not np.isfinite(value).all()
        ):
            raise HLT16StateStoreError(f"array {name} is invalid")
        size = int(value.nbytes)
        if size > MAX_ARRAY_MEMBER_BYTES:
            raise HLT16StateStoreError(f"array {name} exceeds the frozen byte budget")
        answer.append(
            {
                "name": name,
                "dtype": "<f8",
                "shape": list(value.shape),
                "order": "C",
                "byte_count": size,
                "bytes_sha256": sha256(value.tobytes(order="C")).hexdigest(),
            }
        )
    if sum(item["byte_count"] for item in answer) > MAX_TOTAL_UNCOMPRESSED_BYTES:
        raise HLT16StateStoreError("state exceeds the frozen total byte budget")
    return tuple(answer)


def semantic_state_sha256(arrays: Mapping[str, np.ndarray]) -> str:
    return digest(list(_manifest(arrays)))


def encode_payload(
    arrays: Mapping[str, np.ndarray],
) -> tuple[bytes, str, tuple[dict[str, Any], ...]]:
    manifest = _manifest(arrays)
    stream = BytesIO()
    with zipfile.ZipFile(stream, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in ARRAY_NAMES:
            member = BytesIO()
            np.save(member, arrays[name], allow_pickle=False)
            info = zipfile.ZipInfo(f"{name}.npy", (1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100600 << 16
            archive.writestr(info, member.getvalue(), compresslevel=9)
    raw = stream.getvalue()
    if len(raw) > MAX_RAW_ARCHIVE_BYTES:
        raise HLT16StateStoreError("payload exceeds the frozen archive byte budget")
    return raw, digest(list(manifest)), manifest


def decode_payload(
    raw: bytes,
    *,
    expected_semantic_sha256: str | None = None,
    expected_manifest: Iterable[Mapping[str, Any]] | None = None,
) -> tuple[dict[str, np.ndarray], str, tuple[dict[str, Any], ...]]:
    """Decode only after bounding every ZIP member before decompression."""
    if not isinstance(raw, bytes) or len(raw) > MAX_RAW_ARCHIVE_BYTES:
        raise HLT16StateStoreError("payload archive byte budget differs")
    try:
        archive = zipfile.ZipFile(BytesIO(raw), "r")
        infos = archive.infolist()
    except (OSError, zipfile.BadZipFile) as error:
        raise HLT16StateStoreError("payload ZIP is invalid") from error
    expected_names = tuple(f"{name}.npy" for name in ARRAY_NAMES)
    try:
        observed_names = tuple(info.filename for info in infos)
        if observed_names != expected_names or len(set(observed_names)) != len(observed_names):
            raise HLT16StateStoreError("payload ZIP inventory/order differs")
        total = 0
        for info in infos:
            mode = (info.external_attr >> 16) & 0o170000
            if (
                info.is_dir()
                or info.flag_bits & 0x1
                or info.file_size <= 0
                or info.file_size > MAX_ARRAY_MEMBER_BYTES + (1 << 20)
                or info.compress_size < 0
                or mode not in {0, stat.S_IFREG}
                or Path(info.filename).name != info.filename
            ):
                raise HLT16StateStoreError("payload ZIP member is unsafe")
            total += info.file_size
        if total > MAX_TOTAL_UNCOMPRESSED_BYTES + len(infos) * (1 << 20):
            raise HLT16StateStoreError("payload ZIP expansion budget differs")
        data: dict[str, np.ndarray] = {}
        for name, info in zip(ARRAY_NAMES, infos, strict=True):
            encoded = archive.read(info)
            if len(encoded) != info.file_size:
                raise HLT16StateStoreError("payload ZIP member size differs")
            try:
                value = np.load(BytesIO(encoded), allow_pickle=False)
            except (OSError, ValueError) as error:
                raise HLT16StateStoreError("payload NPY member is invalid") from error
            if (
                not isinstance(value, np.ndarray)
                or value.dtype != np.dtype("<f8")
                or not value.flags.c_contiguous
                or not value.shape
            ):
                raise HLT16StateStoreError("payload dtype/layout differs")
            data[name] = value.copy(order="C")
    finally:
        archive.close()
    manifest = _manifest(data)
    if expected_manifest is not None and list(manifest) != list(expected_manifest):
        raise HLT16StateStoreError("payload manifest differs")
    semantic = digest(list(manifest))
    if expected_semantic_sha256 is not None and semantic != _sha(
        expected_semantic_sha256, "semantic state"
    ):
        raise HLT16StateStoreError("payload semantic identity differs")
    return data, semantic, manifest


@dataclass(frozen=True, slots=True)
class AuthenticatedGen0Import:
    schema: str
    campaign_id: str
    member_key: str
    source_receipt_sha256: str
    source_checkpoint_sha256: str
    source_state_object_sha256: str
    authority_sha256: str

    def body(self) -> dict[str, str]:
        return {
            "schema": self.schema,
            "campaign_id": self.campaign_id,
            "member_key": self.member_key,
            "source_receipt_sha256": self.source_receipt_sha256,
            "source_checkpoint_sha256": self.source_checkpoint_sha256,
            "source_state_object_sha256": self.source_state_object_sha256,
        }

    def validate(self) -> None:
        if (
            self.schema != GEN0_AUTHORITY_SCHEMA
            or not self.campaign_id
            or not self.member_key
            or self.source_receipt_sha256 != SEALED_GEN0_RECEIPT_SHA256
            or self.source_checkpoint_sha256 != SEALED_GEN0_CHECKPOINT_SHA256
            or digest(self.body()) != self.authority_sha256
        ):
            raise HLT16StateStoreError("GEN0 import authority differs")
        _sha(self.source_state_object_sha256, "GEN0 source state")


def authorize_gen0_import(snapshot: HLT16MemberSnapshot) -> AuthenticatedGen0Import:
    """Create a narrow authority from an already validated sealed-GEN0 snapshot."""
    try:
        metadata = validate_codec_metadata(snapshot.metadata)
    except Exception as error:
        raise HLT16StateStoreError("GEN0 snapshot metadata differs") from error
    runtime = metadata["runtime_identity"]
    provenance = metadata["provenance_cursor"]
    if runtime["generation"] != 0 or provenance["source_protocol"] != GEN0_PROTOCOL:
        raise HLT16StateStoreError("GEN0 authority cannot cover an evolved state")
    body = {
        "schema": GEN0_AUTHORITY_SCHEMA,
        "campaign_id": runtime["campaign_id"],
        "member_key": runtime["member_key"],
        "source_receipt_sha256": SEALED_GEN0_RECEIPT_SHA256,
        "source_checkpoint_sha256": SEALED_GEN0_CHECKPOINT_SHA256,
        "source_state_object_sha256": provenance["payload"]["accepted_state_sha256"],
    }
    authority = AuthenticatedGen0Import(**body, authority_sha256=digest(body))
    authority.validate()
    return authority


@dataclass(frozen=True, slots=True)
class PersistedState:
    semantic_sha256: str
    descriptor_sha256: str
    raw_archive_sha256: str
    descriptor: Mapping[str, Any]
    arrays: Mapping[str, np.ndarray]


class HLT16StateStore:
    def __init__(self, root: Path, *, fault_hook: FaultHook | None = None) -> None:
        self.root = Path(root)
        self.fault_hook = fault_hook

    def _cut(self, phase: str) -> None:
        if self.fault_hook is not None:
            self.fault_hook(phase)

    def _directory(self, relative: str) -> Path:
        current = self.root
        try:
            root_info = current.lstat()
        except FileNotFoundError:
            current.mkdir(mode=0o755)
            root_info = current.lstat()
        if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
            raise HLT16StateStoreError("store root is unsafe")
        for component in _safe_relative(relative):
            child = current / component
            try:
                info = child.lstat()
            except FileNotFoundError:
                parent_fd = os.open(
                    current,
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                )
                try:
                    os.mkdir(component, mode=0o755, dir_fd=parent_fd)
                finally:
                    os.close(parent_fd)
                info = child.lstat()
            if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
                raise HLT16StateStoreError("store directory is unsafe")
            current = child
        return current

    def _publish(self, relative: str, payload: bytes) -> None:
        parts = _safe_relative(relative)
        self._directory("/".join(parts[:-1]))
        atomic_publish_noreplace(
            self.root,
            relative,
            payload,
            fault_hook=self._cut,
        )

    def persist_snapshot(
        self,
        snapshot: HLT16MemberSnapshot,
        *,
        gen0_authority: AuthenticatedGen0Import | None = None,
    ) -> PersistedState:
        """Import only authenticated generation-zero state into a campaign.

        Evolved persistence is intentionally not a public low-level API: it
        must pass through ``HLT16CampaignStore.persist_snapshot`` where the
        active writer capability and predecessor checkpoint are checked.
        """
        if not isinstance(snapshot, HLT16MemberSnapshot):
            raise HLT16StateStoreError("snapshot type differs")
        try:
            generation = validate_codec_metadata(snapshot.metadata)["runtime_identity"]["generation"]
        except Exception as error:
            raise HLT16StateStoreError("codec metadata differs") from error
        if generation != 0:
            raise HLT16StateStoreError("evolved persistence requires campaign writer authority")
        return self._persist_snapshot_internal(snapshot, gen0_authority=gen0_authority)

    def _persist_snapshot_internal(
        self,
        snapshot: HLT16MemberSnapshot,
        *,
        gen0_authority: AuthenticatedGen0Import | None = None,
    ) -> PersistedState:
        """Canonical publisher owned by the campaign store after GEN0."""
        if not isinstance(snapshot, HLT16MemberSnapshot):
            raise HLT16StateStoreError("snapshot type differs")
        try:
            metadata = validate_codec_metadata(snapshot.metadata)
        except Exception as error:
            raise HLT16StateStoreError("codec metadata differs") from error
        generation = metadata["runtime_identity"]["generation"]
        if generation == 0:
            if not isinstance(gen0_authority, AuthenticatedGen0Import):
                raise HLT16StateStoreError("GEN0 import authority is absent")
            gen0_authority.validate()
            provenance = metadata["provenance_cursor"]["payload"]
            if (
                gen0_authority.campaign_id != metadata["runtime_identity"]["campaign_id"]
                or gen0_authority.member_key != metadata["runtime_identity"]["member_key"]
                or gen0_authority.source_state_object_sha256
                != provenance["accepted_state_sha256"]
            ):
                raise HLT16StateStoreError("GEN0 authority/snapshot relation differs")
        elif gen0_authority is not None:
            raise HLT16StateStoreError("raw GEN0 authority is forbidden after generation zero")

        payload, semantic, manifest = encode_payload(snapshot.arrays)
        raw_hash = sha256(payload).hexdigest()
        body = {
            "schema": STATE_SCHEMA,
            "semantic_sha256": semantic,
            "raw_archive_sha256": raw_hash,
            "arrays": list(manifest),
            "metadata": metadata,
        }
        address = digest(body)
        descriptor = {**body, "descriptor_sha256": address}
        self._publish(f"payloads/{semantic}.npz", payload)
        self._publish(f"states/{address}.json", canonical(descriptor))
        return PersistedState(
            semantic,
            address,
            raw_hash,
            descriptor,
            {name: value.copy(order="C") for name, value in snapshot.arrays.items()},
        )

    def load_state(self, address: str) -> PersistedState:
        address = _sha(address, "descriptor")
        descriptor = _decode_json(
            read_nofollow(self.root, f"states/{address}.json", "descriptor"),
            "descriptor",
        )
        if descriptor.get("schema") != STATE_SCHEMA:
            raise HLT16StateStoreError("descriptor schema differs")
        observed = descriptor.get("descriptor_sha256")
        bare = dict(descriptor)
        bare.pop("descriptor_sha256", None)
        if observed != address or digest(bare) != address:
            raise HLT16StateStoreError("descriptor address differs")
        try:
            metadata = validate_codec_metadata(descriptor["metadata"])
        except Exception as error:
            raise HLT16StateStoreError("descriptor codec metadata differs") from error
        raw = read_nofollow(
            self.root,
            f"payloads/{descriptor['semantic_sha256']}.npz",
            "payload",
        )
        if sha256(raw).hexdigest() != descriptor.get("raw_archive_sha256"):
            raise HLT16StateStoreError("payload raw hash differs")
        arrays, semantic, manifest = decode_payload(
            raw,
            expected_semantic_sha256=descriptor["semantic_sha256"],
            expected_manifest=descriptor["arrays"],
        )
        return PersistedState(
            semantic,
            address,
            sha256(raw).hexdigest(),
            {**descriptor, "metadata": metadata},
            arrays,
        )


__all__ = [
    "AuthenticatedGen0Import",
    "GEN0_AUTHORITY_SCHEMA",
    "HLT16StateStore",
    "HLT16StateStoreError",
    "MAX_ARRAY_MEMBER_BYTES",
    "MAX_RAW_ARCHIVE_BYTES",
    "MAX_TOTAL_UNCOMPRESSED_BYTES",
    "PersistedState",
    "SEALED_GEN0_CHECKPOINT_SHA256",
    "SEALED_GEN0_RECEIPT_SHA256",
    "STATE_SCHEMA",
    "authorize_gen0_import",
    "atomic_publish_noreplace",
    "canonical",
    "decode_payload",
    "digest",
    "encode_payload",
    "read_nofollow",
    "semantic_state_sha256",
]
