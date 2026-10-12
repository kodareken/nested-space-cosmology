"""Read-only admission of the two real PRO18 restart-source containers.

Unlike HLT14's synthetic six-file adapter, the production evidence consists of
two legacy shared NPZ containers.  This module never creates or inspects a
PROTO17 output namespace.  It takes one immutable snapshot of each source
file, hashes that exact byte string, and parses only that byte string.  This
prevents a path being hashed and a later, different path image being parsed.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from io import BytesIO
import json
import math
import os
from pathlib import Path
import stat
import struct
from types import MappingProxyType
from typing import Any, Mapping, Sequence
import zipfile

import numpy as np

from .numerical_engine import array_content_sha256


ARRAY_SUFFIXES = (
    "u", "p", "q", "tracer_positions", "tracer_proper_times",
    "event_proper_times", "event_fields",
)
_ARRAY_SHAPES = {
    "u": lambda points: (points, 6),
    "p": lambda points: (points, 6),
    "q": lambda points: (points, 6),
    "tracer_positions": lambda _points: (48,),
    "tracer_proper_times": lambda _points: (48,),
    "event_proper_times": lambda _points: (24, 48),
    "event_fields": lambda _points: (24, 48, 6),
}
_MAX_SOURCE_BYTES = 128 << 20
_MAX_METADATA_BYTES = 2 << 20
_MAX_NPY_HEADER_BYTES = 65536
_HEX64 = frozenset("0123456789abcdef")
_PROTO12_MEMBER_FIELDS = frozenset({
    "CFL_retry_count", "causal_state", "input_hash", "runtime_monitor_state",
    "source_retry_count", "step_index", "time", "transaction_serial",
})
_RSP2_MEMBER_FIELDS = _PROTO12_MEMBER_FIELDS | {"state_sha256"}


class Proto18ProductionInputError(ValueError):
    """A real source-container premise failed before any launch operation."""


class Proto18ProductionPathError(Proto18ProductionInputError):
    """A source path is unsafe, symlinked, mutable during snapshot, or absent."""


class Proto18ProductionArchiveError(Proto18ProductionInputError):
    """A source archive differs from the frozen shared-container grammar."""


class Proto18ProductionResourceError(Proto18ProductionInputError):
    """A source or archive member exceeds a declared admission budget."""


class Proto18ProductionSourceIdentityError(Proto18ProductionInputError):
    """A safely captured source byte image differs from its frozen identity."""


class Proto18ProductionSelectorError(Proto18ProductionInputError):
    """One logical restart member differs from its frozen compact identity."""


@dataclass(frozen=True, slots=True)
class SourceSnapshot:
    """One byte image and its hash; all parsing must consume ``payload`` only."""

    relative_path: str
    payload: bytes
    sha256: str
    device: int
    inode: int
    size: int
    mtime_ns: int
    ctime_ns: int


@dataclass(frozen=True, slots=True)
class MemberSelector:
    key: str
    source: str
    archive_index: int
    method: str
    point_count: int
    array_prefix: str
    state_sha256: str
    restart_sha256: str
    input_sha256: str
    coordinate_time: str
    source_retry_count: int
    cfl_retry_count: int
    accepted_stage_count: int
    event_sample_count: int
    tracer_count: int
    state_shape: tuple[int, int]
    step_index: int
    transaction_serial: int


@dataclass(frozen=True, slots=True)
class SourceArchiveSpec:
    source: str
    archive_index: int
    archive_path: str
    event_log_path: str
    selectors: tuple[MemberSelector, ...]


@dataclass(frozen=True, slots=True)
class VerifiedProductionArchive:
    spec: SourceArchiveSpec
    archive_snapshot: SourceSnapshot
    event_snapshot: SourceSnapshot
    archive_member_names: tuple[str, ...]
    metadata: Mapping[str, Any]
    metadata_sha256: str
    embedded_event_log_payload: bytes
    embedded_event_log_sha256: str
    total_uncompressed_bytes: int
    largest_uncompressed_member_bytes: int


@dataclass(frozen=True, slots=True)
class VerifiedProductionMember:
    selector: MemberSelector
    source_archive_sha256: str
    source_event_log_sha256: str
    array_sha256: Mapping[str, str]
    physical_state_sha256: str
    restart_payload_sha256: str
    metadata: Mapping[str, Any]
    arrays: Mapping[str, np.ndarray]


@dataclass(frozen=True, slots=True)
class VerifiedProductionCorpus:
    archives: Mapping[str, VerifiedProductionArchive]
    members: Mapping[str, VerifiedProductionMember]


def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and set(value) <= _HEX64


def _finite_int(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise Proto18ProductionSelectorError(f"{label} is not a nonnegative integer")
    return value


def _safe_relative(relative: object, *, label: str) -> str:
    if not isinstance(relative, str) or not relative:
        raise Proto18ProductionPathError(f"{label} is not nonempty text")
    path = Path(relative)
    if path.is_absolute() or "\\" in relative or any(part in {"", ".", ".."} for part in path.parts):
        raise Proto18ProductionPathError(f"{label} is not a safe relative path")
    if not relative.startswith("runs/fgc-2-sf1/"):
        raise Proto18ProductionPathError(f"{label} escapes the frozen raw-source corpus")
    return relative


def _regular_snapshot(root: Path, relative: object, *, label: str) -> SourceSnapshot:
    """Read one stable regular-file image without following a symlink.

    The returned payload is the sole object passed to ZIP/JSON/NPY parsing.
    We compare lstat before and after the descriptor read so a path replacement
    is a typed provenance stop rather than an unobserved mixed-image read.
    """
    rel = _safe_relative(relative, label=label)
    root_path = Path(root)
    if not root_path.is_absolute():
        raise Proto18ProductionPathError("repository root must be absolute")
    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_DIRECTORY", 0)
    )
    file_flags = (
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    )
    directory_descriptors: list[int] = []
    descriptor: int | None = None
    try:
        root_before = root_path.lstat()
        if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
            raise Proto18ProductionPathError("repository root is unsafe")
        root_descriptor = os.open(root_path, directory_flags)
        directory_descriptors.append(root_descriptor)
        root_open = os.fstat(root_descriptor)
        if ((root_open.st_dev, root_open.st_ino) != (root_before.st_dev, root_before.st_ino)
                or not stat.S_ISDIR(root_open.st_mode)):
            raise Proto18ProductionPathError("repository root changed before snapshot")

        parent = root_descriptor
        chain: list[tuple[int, str, int, int]] = []
        parts = Path(rel).parts
        for part in parts[:-1]:
            before_directory = os.stat(part, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before_directory.st_mode) or not stat.S_ISDIR(before_directory.st_mode):
                raise Proto18ProductionPathError(
                    f"{label} has a symlinked or non-directory ancestor"
                )
            child = os.open(part, directory_flags, dir_fd=parent)
            directory_descriptors.append(child)
            opened_directory = os.fstat(child)
            if ((opened_directory.st_dev, opened_directory.st_ino)
                    != (before_directory.st_dev, before_directory.st_ino)
                    or not stat.S_ISDIR(opened_directory.st_mode)):
                raise Proto18ProductionPathError(f"{label} ancestor changed during traversal")
            chain.append((parent, part, before_directory.st_dev, before_directory.st_ino))
            parent = child

        leaf = parts[-1]
        before = os.stat(leaf, dir_fd=parent, follow_symlinks=False)
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_size < 0:
            raise Proto18ProductionPathError(f"{label} is not a safe regular file")
        if before.st_size > _MAX_SOURCE_BYTES:
            raise Proto18ProductionResourceError(f"{label} exceeds the source-byte budget")
        descriptor = os.open(leaf, file_flags, dir_fd=parent)
        opened = os.fstat(descriptor)
        if (not stat.S_ISREG(opened.st_mode)
                or (opened.st_dev, opened.st_ino) != (before.st_dev, before.st_ino)):
            raise Proto18ProductionPathError(f"{label} changed before snapshot")
        chunks: list[bytes] = []
        remaining = _MAX_SOURCE_BYTES + 1
        while remaining:
            chunk = os.read(descriptor, min(1 << 20, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        payload = b"".join(chunks)
        after_open = os.fstat(descriptor)
        after_path = os.stat(leaf, dir_fd=parent, follow_symlinks=False)
        for ancestor_parent, part, device, inode in chain:
            after_directory = os.stat(part, dir_fd=ancestor_parent, follow_symlinks=False)
            if ((after_directory.st_dev, after_directory.st_ino) != (device, inode)
                    or stat.S_ISLNK(after_directory.st_mode)):
                raise Proto18ProductionPathError(f"{label} ancestor changed during snapshot")
        root_after = root_path.lstat()
        if ((root_after.st_dev, root_after.st_ino) != (root_before.st_dev, root_before.st_ino)
                or stat.S_ISLNK(root_after.st_mode)):
            raise Proto18ProductionPathError("repository root changed during snapshot")
    except Proto18ProductionInputError:
        raise
    except OSError as error:
        raise Proto18ProductionPathError(f"{label} cannot be traversed or read safely") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
        for directory_descriptor in reversed(directory_descriptors):
            os.close(directory_descriptor)
    if len(payload) > _MAX_SOURCE_BYTES:
        raise Proto18ProductionResourceError(f"{label} exceeds the source-byte budget")
    if (len(payload) != before.st_size
            or (after_open.st_dev, after_open.st_ino, after_open.st_size,
                after_open.st_mtime_ns, after_open.st_ctime_ns)
            != (before.st_dev, before.st_ino, before.st_size,
                before.st_mtime_ns, before.st_ctime_ns)
            or (after_path.st_dev, after_path.st_ino, after_path.st_size,
                after_path.st_mtime_ns, after_path.st_ctime_ns)
            != (before.st_dev, before.st_ino, before.st_size,
                before.st_mtime_ns, before.st_ctime_ns)):
        raise Proto18ProductionPathError(f"{label} changed during snapshot")
    return SourceSnapshot(
        rel, payload, sha256(payload).hexdigest(), before.st_dev, before.st_ino,
        len(payload), before.st_mtime_ns, before.st_ctime_ns,
    )


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate key {key}")
        result[key] = value
    return result


def _legacy_json(raw: bytes, *, label: str) -> dict[str, Any]:
    def reject_constant(value: str) -> None:
        raise ValueError(f"nonfinite JSON constant {value}")
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=_no_duplicates, parse_constant=reject_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Proto18ProductionArchiveError(f"{label} is malformed legacy JSON") from error
    if not isinstance(value, dict):
        raise Proto18ProductionArchiveError(f"{label} is not an object")
    def require_finite(item: object) -> None:
        if isinstance(item, float) and not math.isfinite(item):
            raise Proto18ProductionArchiveError(f"{label} contains a nonfinite number")
        if isinstance(item, Mapping):
            for child in item.values():
                require_finite(child)
        elif isinstance(item, list):
            for child in item:
                require_finite(child)
    require_finite(value)
    return value


def _freeze(value: Any) -> Any:
    """Recursively freeze admitted JSON-like evidence after validation."""
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(child) for child in value)
    return value


def _selector(value: Mapping[str, Any]) -> MemberSelector:
    required = {
        "key", "source", "archive_index", "method", "point_count", "array_prefix",
        "state_sha256", "physical_state_sha256", "restart_sha256", "input_sha256",
        "coordinate_time", "source_retry_count", "CFL_retry_count", "accepted_stage_count",
        "event_sample_count", "tracer_count", "state_shape", "step_index", "transaction_serial",
    }
    if set(value) != required:
        raise Proto18ProductionSelectorError("selector schema differs")
    if value["state_sha256"] != value["physical_state_sha256"]:
        raise Proto18ProductionSelectorError("selector physical and state hashes differ")
    strings = ("key", "source", "method", "array_prefix", "coordinate_time")
    if any(not isinstance(value[name], str) or not value[name] for name in strings):
        raise Proto18ProductionSelectorError("selector textual identity differs")
    if value["source"] not in {"PROTO12", "RSP2"} or value["method"] not in {"RK4", "SSPRK3"}:
        raise Proto18ProductionSelectorError("selector source or method differs")
    for name in ("state_sha256", "physical_state_sha256", "restart_sha256", "input_sha256"):
        if not _is_sha256(value[name]):
            raise Proto18ProductionSelectorError(f"selector {name} is not SHA-256")
    points = _finite_int(value["point_count"], "selector point count")
    shape = value["state_shape"]
    if (not isinstance(shape, list) or len(shape) != 2 or shape != [points, 6]):
        raise Proto18ProductionSelectorError("selector state shape differs")
    for name in ("archive_index", "source_retry_count", "CFL_retry_count", "accepted_stage_count", "event_sample_count", "tracer_count", "step_index", "transaction_serial"):
        _finite_int(value[name], f"selector {name}")
    if value["coordinate_time"] != "23/16" or value["event_sample_count"] != 24 or value["tracer_count"] != 48:
        raise Proto18ProductionSelectorError("selector restart/event/tracer contract differs")
    return MemberSelector(
        key=value["key"], source=value["source"], archive_index=value["archive_index"],
        method=value["method"], point_count=points, array_prefix=value["array_prefix"],
        state_sha256=value["state_sha256"], restart_sha256=value["restart_sha256"],
        input_sha256=value["input_sha256"], coordinate_time=value["coordinate_time"],
        source_retry_count=value["source_retry_count"], cfl_retry_count=value["CFL_retry_count"],
        accepted_stage_count=value["accepted_stage_count"], event_sample_count=value["event_sample_count"],
        tracer_count=value["tracer_count"], state_shape=(points, 6),
        step_index=value["step_index"], transaction_serial=value["transaction_serial"],
    )


def load_source_specs(config: Mapping[str, Any]) -> tuple[SourceArchiveSpec, ...]:
    """Load the frozen two-container/six-selector map without touching raw bytes."""
    if not isinstance(config, Mapping):
        raise Proto18ProductionSelectorError("PRO18 config is not a mapping")
    source_map, selectors = config.get("source_map"), config.get("selectors")
    if not isinstance(source_map, Mapping) or not isinstance(selectors, Mapping):
        raise Proto18ProductionSelectorError("PRO18 source map is absent")
    paths, logs = source_map.get("shared_containers"), source_map.get("legacy_event_logs")
    members = selectors.get("member")
    if (source_map.get("container_format") != "NPZ" or source_map.get("source_container_count") != 2
            or not isinstance(paths, list) or not isinstance(logs, list) or len(paths) != 2 or len(logs) != 2
            or not isinstance(members, list) or len(members) != 6):
        raise Proto18ProductionSelectorError("PRO18 two-container source grammar differs")
    expected_paths = [
        "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
        "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
    ]
    expected_logs = [
        "runs/fgc-2-sf1/proto12/calibration/events.jsonl",
        "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/events.jsonl",
    ]
    if paths != expected_paths or logs != expected_logs:
        raise Proto18ProductionSelectorError("PRO18 exact source paths differ")
    parsed = tuple(_selector(item) for item in members if isinstance(item, Mapping))
    if len(parsed) != 6 or len({item.key for item in parsed}) != 6 or len({item.array_prefix for item in parsed}) != 6:
        raise Proto18ProductionSelectorError("PRO18 selector identities are not one-to-one")
    if tuple(item.key for item in parsed) != ("RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"):
        raise Proto18ProductionSelectorError("PRO18 canonical selector order differs")
    for item in parsed:
        method, points = item.key.split("-", 1)
        if (item.method != method or item.point_count != int(points)
                or item.array_prefix != item.key.replace("-", "_")
                or item.archive_index != (1 if item.source == "RSP2" else 0)):
            raise Proto18ProductionSelectorError("PRO18 selector method/grid/prefix ownership differs")
    specs: list[SourceArchiveSpec] = []
    for index, source in enumerate(("PROTO12", "RSP2")):
        chosen = tuple(item for item in parsed if item.archive_index == index)
        if not chosen or any(item.source != source for item in chosen):
            raise Proto18ProductionSelectorError("PRO18 source/archive ownership differs")
        specs.append(SourceArchiveSpec(source, index, _safe_relative(paths[index], label="source archive"), _safe_relative(logs[index], label="source event log"), chosen))
    if {item.key for item in specs[0].selectors} != {"RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193"} or {item.key for item in specs[1].selectors} != {"SSPRK3-16385"}:
        raise Proto18ProductionSelectorError("PROTO12/RSP2 source corpus differs")
    return tuple(specs)


def _npy_member_bytes(archive: zipfile.ZipFile, info: zipfile.ZipInfo, *, expected_size: int, label: str) -> bytes:
    if info.file_size > expected_size + _MAX_NPY_HEADER_BYTES or info.file_size <= 0:
        raise Proto18ProductionResourceError(f"{label} NPY member size exceeds its declaration")
    if info.compress_size <= 0 or info.file_size > info.compress_size * 256 + _MAX_NPY_HEADER_BYTES:
        raise Proto18ProductionResourceError(f"{label} NPY compression ratio is unsafe")
    try:
        payload = archive.read(info)
    except (KeyError, OSError, RuntimeError, zipfile.BadZipFile) as error:
        raise Proto18ProductionArchiveError(f"{label} NPY member is unreadable") from error
    if len(payload) != info.file_size:
        raise Proto18ProductionArchiveError(f"{label} NPY member was truncated")
    return payload


def _load_array(member: bytes, *, expected_shape: tuple[int, ...], label: str) -> np.ndarray:
    """Validate NPY header before NumPy is allowed to allocate the payload."""
    stream = BytesIO(member)
    try:
        major, minor = np.lib.format.read_magic(stream)
        if (major, minor) != (1, 0):
            raise ValueError("unsupported NPY version")
        shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
    except (ValueError, EOFError) as error:
        raise Proto18ProductionArchiveError(f"{label} has malformed NPY header") from error
    if dtype != np.dtype("<f8") or fortran or tuple(shape) != expected_shape:
        raise Proto18ProductionArchiveError(f"{label} NPY dtype/layout/shape differs")
    offset = stream.tell()
    if offset > _MAX_NPY_HEADER_BYTES or len(member) != offset + math.prod(expected_shape) * 8:
        raise Proto18ProductionArchiveError(f"{label} NPY byte layout differs")
    try:
        array = np.frombuffer(
            member,
            dtype=np.dtype("<f8"),
            count=math.prod(expected_shape),
            offset=offset,
        ).reshape(expected_shape)
    except ValueError as error:
        raise Proto18ProductionArchiveError(f"{label} NPY payload is unreadable") from error
    if array.dtype != np.dtype("<f8") or not array.flags.c_contiguous or tuple(array.shape) != expected_shape or not np.isfinite(array).all():
        raise Proto18ProductionArchiveError(f"{label} decoded array differs")
    return array


def _load_metadata(member: bytes, *, label: str) -> tuple[dict[str, Any], str]:
    stream = BytesIO(member)
    try:
        major, minor = np.lib.format.read_magic(stream)
        if (major, minor) != (1, 0):
            raise ValueError("unsupported NPY version")
        shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
    except (ValueError, EOFError) as error:
        raise Proto18ProductionArchiveError(f"{label} metadata header is malformed") from error
    offset = stream.tell()
    if len(shape) != 1:
        raise Proto18ProductionArchiveError(f"{label} metadata NPY shape differs")
    if shape[0] > _MAX_METADATA_BYTES:
        raise Proto18ProductionResourceError(f"{label} metadata exceeds its declaration")
    if dtype != np.dtype("uint8") or fortran or offset > _MAX_NPY_HEADER_BYTES or len(member) != offset + shape[0]:
        raise Proto18ProductionArchiveError(f"{label} metadata NPY layout differs")
    raw = member[offset:]
    return _legacy_json(raw, label=label), sha256(raw).hexdigest()


def _expected_member_names(spec: SourceArchiveSpec) -> set[str]:
    prefixes = {item.array_prefix for item in spec.selectors}
    if spec.source == "PROTO12":
        prefixes.add("SSPRK3_2049")
    return {"metadata_utf8.npy", "event_log_utf8.npy"} | {
        f"{prefix}_{suffix}.npy" for prefix in prefixes for suffix in ARRAY_SUFFIXES
    }


def _archive_members(snapshot: SourceSnapshot, spec: SourceArchiveSpec) -> tuple[zipfile.ZipFile, dict[str, zipfile.ZipInfo]]:
    try:
        archive = zipfile.ZipFile(BytesIO(snapshot.payload))
        infos = archive.infolist()
    except (EOFError, OSError, ValueError, zipfile.BadZipFile) as error:
        raise Proto18ProductionArchiveError(f"{spec.source} source archive is not ZIP/NPZ") from error
    names = [item.filename for item in infos]
    expected = _expected_member_names(spec)
    if (len(names) != len(set(names)) or set(names) != expected
            or any(item.filename.startswith(".") or "/" in item.filename or "\\" in item.filename
                   or item.flag_bits != 0 or item.is_dir()
                   or item.compress_type not in {zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED}
                   or stat.S_IFMT(item.external_attr >> 16) == stat.S_IFLNK for item in infos)):
        archive.close()
        raise Proto18ProductionArchiveError(f"{spec.source} shared-container ZIP inventory differs")
    return archive, {item.filename: item for item in infos}


def _validate_archive_metadata(spec: SourceArchiveSpec, metadata: Mapping[str, Any]) -> None:
    common = {"schema_version", "campaign_id", "amplitude", "amplitude_index", "completed_amplitude_records", "completed_common_event_index", "consecutive_qualified_events", "event_log", "members", "terminal"}
    expected = common if spec.source == "PROTO12" else common | {"manifest_sha256", "study_id", "terminal_result"}
    if set(metadata) != expected or metadata.get("schema_version") != 1 or metadata.get("completed_common_event_index") != 23 or metadata.get("terminal") is not True:
        raise Proto18ProductionArchiveError(f"{spec.source} archive metadata grammar differs")
    members = metadata.get("members")
    expected_keys = ({item.key for item in spec.selectors} | ({"SSPRK3-2049"} if spec.source == "PROTO12" else set()))
    if not isinstance(members, Mapping) or set(members) != expected_keys:
        raise Proto18ProductionArchiveError(f"{spec.source} metadata member inventory differs")
    expected_fields = _PROTO12_MEMBER_FIELDS if spec.source == "PROTO12" else _RSP2_MEMBER_FIELDS
    for key, value in members.items():
        if not isinstance(value, Mapping) or set(value) != expected_fields:
            raise Proto18ProductionArchiveError(f"{spec.source} metadata grammar differs for {key}")


def _enforce_frozen_source_identity(
    spec: SourceArchiveSpec,
    archive_snapshot: SourceSnapshot,
    event_snapshot: SourceSnapshot,
    expected: Mapping[str, Any] | None,
) -> None:
    """Reject byte drift before ZIP/NPY decompression or semantic parsing."""
    if expected is None:
        return
    if not isinstance(expected, Mapping):
        raise Proto18ProductionSourceIdentityError(
            f"{spec.source} frozen source identity is not a mapping"
        )
    observed = {
        "source_id": spec.source,
        "checkpoint": spec.archive_path,
        "checkpoint_sha256": archive_snapshot.sha256,
        "checkpoint_byte_count": archive_snapshot.size,
        "event_log": spec.event_log_path,
        "event_log_sha256": event_snapshot.sha256,
        "event_log_byte_count": event_snapshot.size,
    }
    if any(expected.get(key) != value for key, value in observed.items()):
        raise Proto18ProductionSourceIdentityError(
            f"{spec.source} safely captured source byte identity differs"
        )


def read_source_archive(
    root: Path,
    spec: SourceArchiveSpec,
    *,
    expected_source_identity: Mapping[str, Any] | None = None,
) -> VerifiedProductionArchive:
    """Snapshot, hash, and structurally admit one actual shared source container."""
    archive_snapshot = _regular_snapshot(root, spec.archive_path, label=f"{spec.source} archive")
    event_snapshot = _regular_snapshot(root, spec.event_log_path, label=f"{spec.source} event log")
    _enforce_frozen_source_identity(
        spec, archive_snapshot, event_snapshot, expected_source_identity,
    )
    archive, infos = _archive_members(archive_snapshot, spec)
    try:
        metadata_member = _npy_member_bytes(archive, infos["metadata_utf8.npy"], expected_size=_MAX_METADATA_BYTES, label=f"{spec.source} metadata")
        mutable_metadata, metadata_hash = _load_metadata(
            metadata_member, label=f"{spec.source} metadata",
        )
        _validate_archive_metadata(spec, mutable_metadata)
        metadata = _freeze(mutable_metadata)
        event_member = _npy_member_bytes(archive, infos["event_log_utf8.npy"], expected_size=_MAX_SOURCE_BYTES, label=f"{spec.source} embedded event log")
        stream = BytesIO(event_member)
        try:
            major, minor = np.lib.format.read_magic(stream)
            if (major, minor) != (1, 0): raise ValueError("unsupported NPY version")
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
        except (ValueError, EOFError) as error:
            raise Proto18ProductionArchiveError(f"{spec.source} embedded event-log header is malformed") from error
        offset = stream.tell()
        if len(shape) != 1:
            raise Proto18ProductionArchiveError(
                f"{spec.source} embedded event-log shape differs"
            )
        if shape[0] > _MAX_SOURCE_BYTES:
            raise Proto18ProductionResourceError(
                f"{spec.source} embedded event log exceeds its declaration"
            )
        if dtype != np.dtype("uint8") or fortran or offset > _MAX_NPY_HEADER_BYTES or len(event_member) != offset + shape[0]:
            raise Proto18ProductionArchiveError(f"{spec.source} embedded event-log layout differs")
        embedded = event_member[offset:]
        if embedded != event_snapshot.payload:
            raise Proto18ProductionArchiveError(f"{spec.source} archive and external event-log bytes differ")
        sizes = [item.file_size for item in infos.values()]
        return VerifiedProductionArchive(
            spec, archive_snapshot, event_snapshot, tuple(sorted(infos)), metadata,
            metadata_hash, embedded, sha256(embedded).hexdigest(), sum(sizes), max(sizes),
        )
    finally:
        archive.close()


def verify_selected_member(archive: VerifiedProductionArchive, selector: MemberSelector) -> VerifiedProductionMember:
    if selector.source != archive.spec.source or selector not in archive.spec.selectors:
        raise Proto18ProductionSelectorError("selector does not belong to its archive")
    zip_archive, infos = _archive_members(archive.archive_snapshot, archive.spec)
    try:
        arrays: dict[str, np.ndarray] = {}
        array_hashes: dict[str, str] = {}
        for suffix in ARRAY_SUFFIXES:
            expected_shape = _ARRAY_SHAPES[suffix](selector.point_count)
            member_name = f"{selector.array_prefix}_{suffix}.npy"
            payload = _npy_member_bytes(zip_archive, infos[member_name], expected_size=math.prod(expected_shape) * 8, label=f"{selector.key} {suffix}")
            array = _load_array(payload, expected_shape=expected_shape, label=f"{selector.key} {suffix}")
            arrays[suffix] = array
            array_hashes[suffix] = sha256(array.tobytes(order="C")).hexdigest()
        metadata = archive.metadata["members"].get(selector.key)
        if not isinstance(metadata, Mapping):
            raise Proto18ProductionSelectorError("selected member metadata is absent")
        values = {
            "input_hash": selector.input_sha256, "time": 23.0 / 16.0,
            "source_retry_count": selector.source_retry_count, "CFL_retry_count": selector.cfl_retry_count,
            "step_index": selector.step_index, "transaction_serial": selector.transaction_serial,
        }
        if any(metadata.get(key) != value for key, value in values.items()):
            raise Proto18ProductionSelectorError(f"selected metadata differs for {selector.key}")
        monitor, causal = metadata.get("runtime_monitor_state"), metadata.get("causal_state")
        if (not isinstance(monitor, Mapping) or monitor.get("accepted_stage_count") != selector.accepted_stage_count
                or monitor.get("last_accepted_time") != 23.0 / 16.0
                or monitor.get("last_transaction_serial") != selector.transaction_serial - 1
                or monitor.get("first_failed_premise") is not None
                or not isinstance(causal, Mapping) or causal.get("accepted_time") != 23.0 / 16.0):
            raise Proto18ProductionSelectorError(f"selected monitor/causal metadata differs for {selector.key}")
        physical = array_content_sha256(arrays["u"], arrays["p"], arrays["q"])
        restart = array_content_sha256(*(arrays[name] for name in ARRAY_SUFFIXES))
        if physical != selector.state_sha256 or restart != selector.restart_sha256:
            raise Proto18ProductionSelectorError(f"selected array hash differs for {selector.key}")
        if selector.source == "RSP2" and metadata.get("state_sha256") != physical:
            raise Proto18ProductionSelectorError("RSP2 metadata physical-state hash differs")
        return VerifiedProductionMember(
            selector,
            archive.archive_snapshot.sha256,
            archive.event_snapshot.sha256,
            MappingProxyType(dict(array_hashes)),
            physical,
            restart,
            metadata,
            MappingProxyType(dict(arrays)),
        )
    finally:
        zip_archive.close()


def verify_production_corpus(
    root: Path,
    config: Mapping[str, Any],
    *,
    source_identities: Mapping[str, Mapping[str, Any]] | None = None,
) -> VerifiedProductionCorpus:
    """Admit exactly the frozen PROTO12/RSP2 six-member source corpus."""
    specs = load_source_specs(config)
    if source_identities is not None and set(source_identities) != {"PROTO12", "RSP2"}:
        raise Proto18ProductionSourceIdentityError(
            "frozen source identity inventory differs"
        )
    archives = {
        spec.source: read_source_archive(
            root,
            spec,
            expected_source_identity=(
                source_identities[spec.source] if source_identities is not None else None
            ),
        )
        for spec in specs
    }
    members: dict[str, VerifiedProductionMember] = {}
    for spec in specs:
        for selector in spec.selectors:
            if selector.key in members:
                raise Proto18ProductionSelectorError("selected member appears twice")
            members[selector.key] = verify_selected_member(archives[spec.source], selector)
    if tuple(members) != ("RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"):
        raise Proto18ProductionSelectorError("complete selected corpus order differs")
    return VerifiedProductionCorpus(
        MappingProxyType(dict(archives)), MappingProxyType(dict(members)),
    )


__all__ = [
    "ARRAY_SUFFIXES", "MemberSelector", "Proto18ProductionArchiveError",
    "Proto18ProductionInputError", "Proto18ProductionPathError",
    "Proto18ProductionResourceError",
    "Proto18ProductionSelectorError", "Proto18ProductionSourceIdentityError",
    "SourceArchiveSpec", "SourceSnapshot",
    "VerifiedProductionArchive", "VerifiedProductionCorpus", "VerifiedProductionMember",
    "load_source_specs", "read_source_archive", "verify_production_corpus",
    "verify_selected_member",
]
