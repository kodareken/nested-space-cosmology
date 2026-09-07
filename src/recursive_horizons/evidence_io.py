"""Future-only neutral evidence I/O helpers.

This module is infrastructure.  It does not classify scientific results,
authorize a campaign, bind a compact artifact, or share decision logic with
any historical runner or binder.  Callers that need those meanings must own
them separately.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import ctypes
import errno
import json
import math
import os
from pathlib import Path
import re
import secrets
import selectors
import signal
import stat
import subprocess
import sys
import time
from typing import Any, Callable, Union


DEFAULT_MAX_BYTES = 16 * 1024 * 1024
GIT_TIMEOUT_SECONDS = 30
MAX_DIRECTORY_ENTRIES = 1024
MAX_GIT_ENTRIES = 16384
MAX_JSON_DEPTH = 64
MAX_PATH_PARTS = 64
MAX_PATH_PART_BYTES = 180
_COMMIT_ID = re.compile(r"\A(?:HEAD|[0-9a-f]{40}|[0-9a-f]{64})\Z")
_OBJECT_ID = re.compile(r"\A(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_GIT_PATH_FORBIDDEN = frozenset("*?[]:")
_GIT_EXECUTABLE_CANDIDATES = (
    "/usr/bin/git",
    "/opt/homebrew/bin/git",
    "/usr/local/bin/git",
)
_GIT_MINIMAL_PATH = "/usr/bin:/bin"
_RENAME_EXCL_DARWIN = 0x0004
_RENAME_NOREPLACE_LINUX = 0x0001
_CHUNK = 65536
FaultHook = Callable[[str], None]


class EvidenceIOError(ValueError):
    """Base failure for the future-only evidence I/O utility."""

    published = False


class CanonicalJSONError(EvidenceIOError):
    """Bytes or values are not canonical, finite, duplicate-free JSON."""


class UnsafePathError(EvidenceIOError):
    """A path escaped its root, followed a link, or failed an identity check."""


class GitQueryError(EvidenceIOError):
    """A closed read-only Git operation was rejected or failed."""


class PrepublicationError(EvidenceIOError):
    """Exclusive publication failed before a destination became visible.

    Staging is removed only when this call created the staging inode and that
    inode is still present under the staging name.  A foreign, colliding, or
    substituted staging object is retained.
    """

    published = False

    def __init__(self, message: str = "", *, staging_retained: bool = False) -> None:
        super().__init__(message)
        self.staging_retained = staging_retained


class UnsupportedPublication(PrepublicationError):
    """The platform cannot offer exclusive no-replace publication."""


class PostpublicationUncertainty(EvidenceIOError):
    """The destination was published, but a later durability or identity step failed.

    The published result is left in place.  This utility does not delete or
    retry it.  Uncertain leftover staging is also retained.
    """

    published = True

    def __init__(self, message: str, *, relative_path: str) -> None:
        super().__init__(message)
        self.relative_path = relative_path


@dataclass(frozen=True, slots=True)
class PublicationReceipt:
    """Record that an exclusive publication completed, fsynced, and revalidated."""

    relative_path: str
    kind: str
    payload_sha256: str


@dataclass(frozen=True, slots=True)
class ResolveCommit:
    """Resolve a revision to a commit object name."""

    revision: str


@dataclass(frozen=True, slots=True)
class ReadBlob:
    """Read one regular blob at a repository-relative path."""

    commit: str
    path: str


@dataclass(frozen=True, slots=True)
class InspectTree:
    """List recursive tree entries at an optional repository-relative path."""

    commit: str
    path: str | None = None


@dataclass(frozen=True, slots=True)
class InspectDelta:
    """List the name-status delta of a commit, optionally against another."""

    commit: str
    against: str | None = None


@dataclass(frozen=True, slots=True)
class CommitParents:
    """Read the ordered parent object names of one commit, including merges."""

    commit: str


@dataclass(frozen=True, slots=True)
class InspectWorktree:
    """Inspect changes and index flags without refreshing or writing the index."""


@dataclass(frozen=True, slots=True)
class WorktreeChange:
    index_status: str
    worktree_status: str
    path: str


@dataclass(frozen=True, slots=True)
class WorktreeInspection:
    changes: tuple[WorktreeChange, ...]
    hidden_index_paths: tuple[str, ...]

    @property
    def clean(self) -> bool:
        return not self.changes and not self.hidden_index_paths


@dataclass(frozen=True, slots=True)
class TreeEntry:
    mode: str
    kind: str
    object_id: str
    path: str


@dataclass(frozen=True, slots=True)
class DeltaEntry:
    status: str
    path: str


GitOperation = Union[
    ResolveCommit, ReadBlob, InspectTree, InspectDelta, CommitParents, InspectWorktree
]


__all__ = (
    "CanonicalJSONError",
    "CommitParents",
    "DEFAULT_MAX_BYTES",
    "DeltaEntry",
    "EvidenceIOError",
    "GitQueryError",
    "InspectDelta",
    "InspectTree",
    "InspectWorktree",
    "MAX_GIT_ENTRIES",
    "MAX_JSON_DEPTH",
    "PostpublicationUncertainty",
    "PrepublicationError",
    "PublicationReceipt",
    "ReadBlob",
    "ResolveCommit",
    "TreeEntry",
    "UnsafePathError",
    "UnsupportedPublication",
    "WorktreeChange",
    "WorktreeInspection",
    "canonical_json_bytes",
    "git_read",
    "load_canonical_json",
    "publish_exclusive_directory",
    "publish_exclusive_file",
    "read_regular_file",
)


def _digest(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _devino(metadata: os.stat_result) -> tuple[int, int]:
    return (metadata.st_dev, metadata.st_ino)


def _leaf_identity(metadata: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
        metadata.st_nlink,
    )


def _require_jsonable(
    value: object,
    *,
    depth: int = 0,
    seen: set[int] | None = None,
) -> None:
    if depth > MAX_JSON_DEPTH:
        raise CanonicalJSONError("JSON value exceeds depth bound")
    if seen is None:
        seen = set()
    if value is None:
        return
    value_type = type(value)
    if value_type is bool or value_type is str:
        return
    if value_type is int:
        return
    if value_type is float:
        if not math.isfinite(value):
            raise CanonicalJSONError("nonfinite JSON number")
        return
    if value_type is list or value_type is dict:
        marker = id(value)
        if marker in seen:
            raise CanonicalJSONError("cyclic JSON value")
        seen.add(marker)
        try:
            if value_type is list:
                for item in value:
                    _require_jsonable(item, depth=depth + 1, seen=seen)
                return
            for key, item in value.items():
                if type(key) is not str:
                    raise CanonicalJSONError("JSON object keys must be strings")
                _require_jsonable(item, depth=depth + 1, seen=seen)
            return
        finally:
            seen.discard(marker)
    raise CanonicalJSONError(f"unsupported JSON value: {value_type.__name__}")


def canonical_json_bytes(value: object) -> bytes:
    """Encode a JSON value as sorted compact ASCII bytes.

    The encoder rejects non-finite numbers, subclasses, cycles, over-deep
    nesting, and values outside the JSON type set.  This is a byte contract,
    not a scientific reduction.
    """

    _require_jsonable(value)
    try:
        raw = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError, UnicodeEncodeError, RecursionError) as error:
        raise CanonicalJSONError("value is not canonical JSON") from error
    if len(raw) > DEFAULT_MAX_BYTES:
        raise CanonicalJSONError("JSON value exceeds the byte bound")
    return raw


def load_canonical_json(raw: bytes) -> Any:
    """Parse canonical JSON bytes, rejecting duplicates and non-finite tokens.

    The input must already match :func:`canonical_json_bytes`.  Pretty-printed
    JSON, duplicate keys, ``NaN``/``Infinity``, and unsupported values fail.
    """

    if type(raw) not in (bytes, bytearray):
        raise CanonicalJSONError("JSON input is not bytes")
    if len(raw) > DEFAULT_MAX_BYTES:
        raise CanonicalJSONError("JSON input exceeds the byte bound")

    def reject_duplicates(
        items: list[tuple[str, object]],
    ) -> dict[str, object]:
        mapping: dict[str, object] = {}
        for key, item in items:
            if key in mapping:
                raise CanonicalJSONError(f"duplicate JSON key {key!r}")
            mapping[key] = item
        return mapping

    def reject_constant(name: str) -> None:
        raise CanonicalJSONError(f"nonfinite JSON token {name}")

    try:
        text = bytes(raw).decode("ascii")
        value = json.loads(
            text,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except CanonicalJSONError:
        raise
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        ValueError,
        TypeError,
        RecursionError,
    ) as error:
        raise CanonicalJSONError("JSON input is malformed") from error
    _require_jsonable(value)
    if canonical_json_bytes(value) != bytes(raw):
        raise CanonicalJSONError("JSON input is not canonical")
    return value


def _flag(name: str) -> int:
    value = int(getattr(os, name, 0) or 0)
    if value == 0:
        raise UnsafePathError(f"{name} is required")
    return value


def _relative_parts(relative: str, *, git: bool = False) -> tuple[str, ...]:
    error_type: type[EvidenceIOError] = GitQueryError if git else UnsafePathError
    if type(relative) is not str or not relative:
        raise error_type("path is empty")
    if "\x00" in relative or "\\" in relative:
        raise error_type("path contains a forbidden character")
    if relative.startswith("/") or relative.endswith("/") or "//" in relative:
        raise error_type("path is not a safe relative POSIX path")
    candidate = Path(relative)
    if candidate.is_absolute() or candidate.anchor:
        raise error_type("path is absolute")
    if candidate.as_posix() != relative:
        raise error_type("path is not a safe relative POSIX path")
    parts = candidate.parts
    if not parts or len(parts) > MAX_PATH_PARTS:
        raise error_type("path has an unsafe number of components")
    for part in parts:
        encoded = part.encode("utf-8")
        if (
            part in {"", ".", ".."}
            or "/" in part
            or "\\" in part
            or len(encoded) > MAX_PATH_PART_BYTES
        ):
            raise error_type("path component is unsafe")
        if git and (any(character in part for character in _GIT_PATH_FORBIDDEN)):
            raise GitQueryError("Git path contains pathspec magic")
    return parts


def _canonical_root(root: Path) -> Path:
    """Require an existing absolute root with no symlinked ancestors."""

    if not isinstance(root, Path):
        raise UnsafePathError("root is not a path")
    original = os.fspath(root)
    if not isinstance(original, str) or not original or "\x00" in original:
        raise UnsafePathError("root is unsafe")
    supplied = Path(original)
    if not supplied.is_absolute() or any(
        part in {".", ".."} for part in supplied.parts
    ):
        raise UnsafePathError("root must be an absolute canonical path")
    try:
        resolved = supplied.resolve(strict=True)
    except OSError as error:
        raise UnsafePathError("root cannot be resolved safely") from error
    if resolved != supplied:
        raise UnsafePathError("root path traverses a symlinked ancestor")
    return supplied


def _close_fds(descriptors: list[int]) -> None:
    for descriptor in reversed(descriptors):
        try:
            os.close(descriptor)
        except OSError:
            pass


def _open_root(root: Path) -> int:
    root = _canonical_root(root)
    directory = _flag("O_DIRECTORY")
    nofollow = _flag("O_NOFOLLOW")
    try:
        before = root.lstat()
    except OSError as error:
        raise UnsafePathError("root cannot be read safely") from error
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
        raise UnsafePathError("root is unsafe")
    descriptor = os.open(root, os.O_RDONLY | directory | nofollow)
    try:
        after = os.fstat(descriptor)
        if _devino(before) != _devino(after):
            raise UnsafePathError("root raced")
        if stat.S_ISLNK(after.st_mode) or not stat.S_ISDIR(after.st_mode):
            raise UnsafePathError("root is unsafe")
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _open_child_directory(parent: int, name: str) -> int:
    directory = _flag("O_DIRECTORY")
    nofollow = _flag("O_NOFOLLOW")
    try:
        before = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except OSError as error:
        raise UnsafePathError("parent cannot be read safely") from error
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
        raise UnsafePathError("parent is unsafe")
    child = os.open(name, os.O_RDONLY | directory | nofollow, dir_fd=parent)
    try:
        after = os.fstat(child)
        if _devino(before) != _devino(after):
            raise UnsafePathError("parent raced")
        if stat.S_ISLNK(after.st_mode) or not stat.S_ISDIR(after.st_mode):
            raise UnsafePathError("parent is unsafe")
        return child
    except Exception:
        os.close(child)
        raise


def _open_chain(
    root: Path, parts: tuple[str, ...]
) -> tuple[list[int], list[tuple[int, int]]]:
    """Open root plus each relative parent, keeping the whole chain live."""

    descriptors: list[int] = []
    identities: list[tuple[int, int]] = []
    try:
        root_fd = _open_root(root)
        descriptors.append(root_fd)
        identities.append(_devino(os.fstat(root_fd)))
        for component in parts:
            child = _open_child_directory(descriptors[-1], component)
            descriptors.append(child)
            identities.append(_devino(os.fstat(child)))
        return descriptors, identities
    except Exception:
        _close_fds(descriptors)
        raise


def _revalidate_chain(
    root: Path,
    parts: tuple[str, ...],
    expected: list[tuple[int, int]],
    live: list[int],
) -> None:
    if len(live) != len(expected):
        raise UnsafePathError("parent chain raced")
    for descriptor, identity in zip(live, expected, strict=True):
        active = os.fstat(descriptor)
        if _devino(active) != identity or not stat.S_ISDIR(active.st_mode):
            raise UnsafePathError("parent chain raced")
    fresh, identities = _open_chain(root, parts)
    try:
        if identities != expected:
            raise UnsafePathError("parent chain raced")
    finally:
        _close_fds(fresh)


def _read_regular_from_parent(
    parent_fd: int,
    name: str,
    *,
    max_bytes: int,
    expected_identity: tuple[int, int] | None = None,
) -> tuple[bytes, tuple[int, int]]:
    before = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > max_bytes
    ):
        raise UnsafePathError("leaf is unsafe")
    identity = _devino(before)
    if expected_identity is not None and identity != expected_identity:
        raise UnsafePathError("leaf raced")
    flags = os.O_RDONLY | _flag("O_NOFOLLOW") | getattr(os, "O_NONBLOCK", 0)
    leaf = os.open(name, flags, dir_fd=parent_fd)
    try:
        active = os.fstat(leaf)
        if (
            _leaf_identity(active) != _leaf_identity(before)
            or not stat.S_ISREG(active.st_mode)
            or active.st_nlink != 1
            or active.st_size != before.st_size
        ):
            raise UnsafePathError("leaf raced")
        chunks: list[bytes] = []
        remaining = int(before.st_size)
        while remaining > 0:
            chunk = os.read(leaf, min(remaining, 1 << 20))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        raw = b"".join(chunks)
        after = os.fstat(leaf)
        named_after = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if (
            len(raw) != before.st_size
            or _leaf_identity(after) != _leaf_identity(before)
            or _leaf_identity(named_after) != _leaf_identity(before)
            or not stat.S_ISREG(after.st_mode)
            or after.st_nlink != 1
            or after.st_size != before.st_size
        ):
            raise UnsafePathError("leaf changed during read")
        return raw, identity
    finally:
        os.close(leaf)


def read_regular_file(
    root: Path,
    relative: str,
    *,
    max_bytes: int = DEFAULT_MAX_BYTES,
) -> bytes:
    """Read a single-linked regular file through a no-follow parent chain.

    The read is bounded, refuses traversal, leaf/ancestor symlinks, hard
    links, and special files, keeps the parent-chain identities live, and
    re-walks the caller-visible path before returning bytes.
    """

    if type(max_bytes) is not int or max_bytes < 0 or max_bytes > DEFAULT_MAX_BYTES:
        raise UnsafePathError("byte bound is unsafe")
    root = _canonical_root(Path(root))
    parts = _relative_parts(relative)
    descriptors, identities = _open_chain(root, parts[:-1])
    try:
        raw, leaf_identity = _read_regular_from_parent(
            descriptors[-1],
            parts[-1],
            max_bytes=max_bytes,
        )
        _revalidate_chain(root, parts[:-1], identities, descriptors)
        fresh, fresh_ids = _open_chain(root, parts[:-1])
        try:
            if fresh_ids != identities:
                raise UnsafePathError("parent chain raced")
            observed, observed_identity = _read_regular_from_parent(
                fresh[-1],
                parts[-1],
                max_bytes=max_bytes,
                expected_identity=leaf_identity,
            )
        finally:
            _close_fds(fresh)
        if observed != raw or observed_identity != leaf_identity:
            raise UnsafePathError("leaf raced")
        _revalidate_chain(root, parts[:-1], identities, descriptors)
        return raw
    except UnsafePathError:
        raise
    except OSError as error:
        raise UnsafePathError("leaf cannot be read safely") from error
    finally:
        _close_fds(descriptors)


def _git_environment() -> dict[str, str]:
    return {
        "LC_ALL": "C",
        "LANG": "C",
        "PATH": _GIT_MINIMAL_PATH,
        "GIT_ATTR_NOSYSTEM": "1",
        "GIT_CONFIG_COUNT": "0",
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_CONFIG_SYSTEM": os.devnull,
        "GIT_NO_REPLACE_OBJECTS": "1",
        "GIT_OPTIONAL_LOCKS": "0",
        "GIT_TERMINAL_PROMPT": "0",
    }


def _git_executable() -> str:
    for candidate in _GIT_EXECUTABLE_CANDIDATES:
        try:
            info = Path(candidate).lstat()
        except OSError:
            continue
        if stat.S_ISREG(info.st_mode) and os.access(candidate, os.X_OK):
            return candidate
        if stat.S_ISLNK(info.st_mode) and os.access(candidate, os.X_OK):
            return candidate
    raise GitQueryError("trusted Git executable is unavailable")


def _validate_revision(revision: str) -> str:
    if not isinstance(revision, str) or not _COMMIT_ID.fullmatch(revision):
        raise GitQueryError("Git revision is not a closed commit identifier")
    return revision


def _communicate_bounded(
    proc: subprocess.Popen[bytes],
    *,
    max_bytes: int,
    timeout: float,
) -> bytes:
    if proc.stdout is None or proc.stderr is None:
        raise GitQueryError("Git capture pipes are absent")
    stdout_chunks: list[bytes] = []
    sizes = {"out": 0, "err": 0}
    deadline = time.monotonic() + timeout
    selector = selectors.DefaultSelector()

    def terminate() -> None:
        try:
            if os.getpgid(proc.pid) == proc.pid:
                os.killpg(proc.pid, signal.SIGKILL)
            else:
                proc.kill()
        except ProcessLookupError:
            pass

    try:
        for stream, kind in ((proc.stdout, "out"), (proc.stderr, "err")):
            os.set_blocking(stream.fileno(), False)
            selector.register(stream.fileno(), selectors.EVENT_READ, kind)
        while selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise GitQueryError("Git read timed out")
            for key, _mask in selector.select(min(0.1, remaining)):
                data = os.read(key.fd, _CHUNK)
                if not data:
                    selector.unregister(key.fd)
                    continue
                kind = key.data
                sizes[kind] += len(data)
                if sizes[kind] > max_bytes:
                    raise GitQueryError("Git read exceeded the byte bound")
                if kind == "out":
                    stdout_chunks.append(data)
        proc.wait(timeout=max(0.01, deadline - time.monotonic()))
        if proc.returncode != 0:
            raise GitQueryError("Git read failed")
        return b"".join(stdout_chunks)
    except Exception as error:
        terminate()
        proc.wait(timeout=1)
        if isinstance(error, GitQueryError):
            raise
        raise GitQueryError("Git read failed") from error
    finally:
        selector.close()
        proc.stdout.close()
        proc.stderr.close()


def _run_git(repository: Path, arguments: tuple[str, ...]) -> bytes:
    if not arguments or any(
        not isinstance(item, str) or item == "" for item in arguments
    ):
        raise GitQueryError("Git argument is not an allowed read-only operand")
    root = _canonical_root(Path(repository))
    root_fd = _open_root(root)
    root_identity = _devino(os.fstat(root_fd))
    os.close(root_fd)
    git_path = _git_executable()
    command = (
        git_path,
        "-C",
        os.fspath(root),
        "--no-replace-objects",
        "--no-optional-locks",
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "-c",
        "core.pager=",
        "-c",
        "core.useReplaceRefs=false",
        "-c",
        "diff.external=",
        "-c",
        "diff.textconv=",
        "-c",
        "alias.rev-parse=",
        "-c",
        "alias.cat-file=",
        "-c",
        "alias.ls-tree=",
        "-c",
        "alias.diff-tree=",
        "-c",
        "alias.status=",
        "-c",
        "alias.ls-files=",
        *arguments,
    )
    try:
        proc = subprocess.Popen(
            command,
            cwd=root,
            env=_git_environment(),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            start_new_session=True,
        )
    except OSError as error:
        raise GitQueryError("Git read failed") from error
    payload = _communicate_bounded(
        proc, max_bytes=DEFAULT_MAX_BYTES, timeout=GIT_TIMEOUT_SECONDS
    )
    checked_root = _open_root(root)
    try:
        if _devino(os.fstat(checked_root)) != root_identity:
            raise GitQueryError("repository root changed during Git query")
    finally:
        os.close(checked_root)
    return payload


def _require_worktree_root(repository: Path) -> None:
    repository = _canonical_root(Path(repository))
    raw = _run_git(repository, ("rev-parse", "--show-toplevel"))
    try:
        top = Path(raw.decode("utf-8").strip())
        if _devino(repository.lstat()) != _devino(top.lstat()):
            raise GitQueryError("repository is not the Git worktree root")
    except (OSError, UnicodeDecodeError) as error:
        raise GitQueryError("repository is not the Git worktree root") from error


def _parse_ls_tree(raw: bytes) -> tuple[TreeEntry, ...]:
    if raw == b"":
        return ()
    if not raw.endswith(b"\0"):
        raise GitQueryError("tree listing is malformed")
    entries: list[TreeEntry] = []
    for item in raw.split(b"\0"):
        if item == b"":
            continue
        try:
            meta, path_raw = item.split(b"\t", 1)
            mode_b, kind_b, object_b = meta.split(b" ")
            path = path_raw.decode("utf-8")
            mode = mode_b.decode("ascii")
            kind = kind_b.decode("ascii")
            object_id = object_b.decode("ascii")
        except (ValueError, UnicodeDecodeError) as error:
            raise GitQueryError("tree listing is malformed") from error
        if (mode, kind) not in {
            ("100644", "blob"),
            ("100755", "blob"),
            ("120000", "blob"),
            ("040000", "tree"),
            ("160000", "commit"),
        } or not _OBJECT_ID.fullmatch(object_id):
            raise GitQueryError("tree entry identity is unsupported")
        _relative_parts(path, git=True)
        entries.append(TreeEntry(mode=mode, kind=kind, object_id=object_id, path=path))
        if len(entries) > MAX_GIT_ENTRIES:
            raise GitQueryError("tree listing exceeded the entry bound")
    return tuple(entries)


def _parse_name_status(raw: bytes) -> tuple[DeltaEntry, ...]:
    if raw == b"":
        return ()
    if not raw.endswith(b"\0"):
        raise GitQueryError("delta listing is malformed")
    parts = raw.split(b"\0")
    entries: list[DeltaEntry] = []
    index = 0
    while index < len(parts):
        status_raw = parts[index]
        if status_raw == b"":
            index += 1
            continue
        if index + 1 >= len(parts):
            raise GitQueryError("delta listing is malformed")
        try:
            status = status_raw.decode("ascii")
            path = parts[index + 1].decode("utf-8")
        except UnicodeDecodeError as error:
            raise GitQueryError("delta listing is malformed") from error
        if status not in {"A", "D", "M", "T", "U", "X", "B"}:
            raise GitQueryError("delta status is unsafe")
        _relative_parts(path, git=True)
        entries.append(DeltaEntry(status=status, path=path))
        if len(entries) > MAX_GIT_ENTRIES:
            raise GitQueryError("delta listing exceeded the entry bound")
        index += 2
    return tuple(entries)


def _parse_worktree(status: bytes, flags: bytes) -> WorktreeInspection:
    changes: list[WorktreeChange] = []
    hidden: list[str] = []
    for raw, label in ((status, "status"), (flags, "index flags")):
        if raw and not raw.endswith(b"\0"):
            raise GitQueryError(f"worktree {label} is not NUL-terminated")
    status_fields = status.split(b"\0")[:-1] if status else []
    flag_fields = flags.split(b"\0")[:-1] if flags else []
    if max(len(status_fields), len(flag_fields)) > MAX_GIT_ENTRIES:
        raise GitQueryError("worktree listing exceeded the entry bound")
    for field in status_fields:
        if len(field) < 4 or field[2:3] != b" ":
            raise GitQueryError("worktree status is malformed")
        try:
            pair, path = field[:2].decode("ascii"), field[3:].decode("utf-8")
        except UnicodeDecodeError as error:
            raise GitQueryError("worktree status is malformed") from error
        # Rename detection is disabled in the query: no second pathname may
        # be silently consumed as a status entry.
        if pair != "??" and (
            any(char not in " MADTU" for char in pair) or pair == "  "
        ):
            raise GitQueryError("worktree status code is unsupported")
        _relative_parts(path, git=True)
        changes.append(WorktreeChange(pair[0], pair[1], path))
    for field in flag_fields:
        if (
            len(field) < 3
            or field[1:2] != b" "
            or field[:1] not in b"HSMRCK?hs mrck".replace(b" ", b"")
        ):
            raise GitQueryError("worktree index flags are malformed")
        try:
            path = field[2:].decode("utf-8")
        except UnicodeDecodeError as error:
            raise GitQueryError("worktree index path is malformed") from error
        _relative_parts(path, git=True)
        if field[:1] != b"H":
            hidden.append(path)
    if len({item.path for item in changes}) != len(changes):
        raise GitQueryError("worktree status repeats a path")
    return WorktreeInspection(tuple(changes), tuple(sorted(set(hidden))))


def _parse_commit_parents(raw: bytes) -> tuple[str, ...]:
    header, separator, _message = raw.partition(b"\n\n")
    lines = header.split(b"\n")
    if not separator or not lines or not lines[0].startswith(b"tree "):
        raise GitQueryError("commit header is malformed")
    try:
        tree = lines[0][5:].decode("ascii")
    except UnicodeDecodeError as error:
        raise GitQueryError("commit tree is malformed") from error
    if not _OBJECT_ID.fullmatch(tree):
        raise GitQueryError("commit tree is malformed")
    parents: list[str] = []
    parents_closed = False
    for line in lines[1:]:
        if not line.startswith(b"parent "):
            parents_closed = True
            continue
        if parents_closed:
            raise GitQueryError("commit parent headers are not contiguous")
        try:
            parent = line[7:].decode("ascii")
        except UnicodeDecodeError as error:
            raise GitQueryError("commit parent is malformed") from error
        if not _OBJECT_ID.fullmatch(parent) or len(parents) >= MAX_GIT_ENTRIES:
            raise GitQueryError("commit parent is invalid or exceeds the bound")
        parents.append(parent)
    return tuple(parents)


def git_read(repository: Path, operation: GitOperation) -> object:
    """Run one closed, sanitized, read-only Git query.

    The process uses a trusted absolute Git executable and a minimal
    environment that does not inherit ``PATH``, ``LD_*``, or ``DYLD_*``.
    stdout and stderr are captured with a hard byte cap.  There is no generic
    mutation or exec command.
    """

    if not isinstance(
        operation,
        (
            ResolveCommit,
            ReadBlob,
            InspectTree,
            InspectDelta,
            CommitParents,
            InspectWorktree,
        ),
    ):
        raise GitQueryError("Git operation is not in the closed read-only set")
    repository = _canonical_root(Path(repository))
    _require_worktree_root(repository)
    if isinstance(operation, InspectWorktree):
        return _parse_worktree(
            _run_git(
                repository,
                (
                    "status",
                    "--porcelain=v1",
                    "-z",
                    "--untracked-files=all",
                    "--ignore-submodules=none",
                    "--no-renames",
                ),
            ),
            _run_git(repository, ("ls-files", "-v", "-z")),
        )
    if isinstance(operation, ResolveCommit):
        revision = _validate_revision(operation.revision)
        raw = _run_git(
            repository,
            ("rev-parse", "--verify", "--end-of-options", f"{revision}^{{commit}}"),
        )
        try:
            resolved = raw.decode("ascii").strip()
        except UnicodeDecodeError as error:
            raise GitQueryError("commit resolution is malformed") from error
        if not _OBJECT_ID.fullmatch(resolved):
            raise GitQueryError("commit resolution is not an object identity")
        return resolved
    commit = _validate_revision(operation.commit)
    if isinstance(operation, CommitParents):
        return _parse_commit_parents(
            _run_git(repository, ("cat-file", "commit", commit))
        )
    if isinstance(operation, ReadBlob):
        parts = _relative_parts(operation.path, git=True)
        path = "/".join(parts)
        listing = _parse_ls_tree(
            _run_git(repository, ("ls-tree", "-z", commit, "--", path))
        )
        if (
            len(listing) != 1
            or listing[0].path != path
            or listing[0].kind != "blob"
            or listing[0].mode not in {"100644", "100755"}
        ):
            raise GitQueryError("blob path is not a regular Git blob")
        raw = _run_git(repository, ("cat-file", "blob", listing[0].object_id))
        if len(raw) > DEFAULT_MAX_BYTES:
            raise GitQueryError("blob exceeded the byte bound")
        return raw
    if isinstance(operation, InspectTree):
        arguments: tuple[str, ...]
        if operation.path is None:
            arguments = ("ls-tree", "-z", "-r", commit)
        else:
            path = "/".join(_relative_parts(operation.path, git=True))
            arguments = ("ls-tree", "-z", "-r", commit, "--", path)
        return _parse_ls_tree(_run_git(repository, arguments))
    against = (
        None if operation.against is None else _validate_revision(operation.against)
    )
    if against is None:
        arguments = (
            "diff-tree",
            "--root",
            "--no-commit-id",
            "--name-status",
            "--no-ext-diff",
            "--no-renames",
            "-r",
            "-z",
            commit,
        )
    else:
        arguments = (
            "diff-tree",
            "--no-commit-id",
            "--name-status",
            "--no-ext-diff",
            "--no-renames",
            "-r",
            "-z",
            against,
            commit,
        )
    return _parse_name_status(_run_git(repository, arguments))


def _load_exclusive_rename() -> tuple[Any, int] | tuple[None, None]:
    if sys.platform == "darwin":
        try:
            libc = ctypes.CDLL("/usr/lib/libSystem.B.dylib", use_errno=True)
            func = libc.renameatx_np
        except (AttributeError, OSError):
            return None, None
        func.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        func.restype = ctypes.c_int
        return func, _RENAME_EXCL_DARWIN
    if sys.platform.startswith("linux"):
        try:
            libc = ctypes.CDLL(None, use_errno=True)
            func = libc.renameat2
        except (AttributeError, OSError):
            return None, None
        func.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        func.restype = ctypes.c_int
        return func, _RENAME_NOREPLACE_LINUX
    return None, None


_RENAME_FUNC, _RENAME_FLAG = _load_exclusive_rename()


def _require_exclusive_rename_support() -> None:
    if _RENAME_FUNC is None or _RENAME_FLAG is None:
        raise UnsupportedPublication(
            "exclusive no-replace rename is unsupported on this platform"
        )


def _exclusive_rename(
    src_dir_fd: int, src_name: str, dst_dir_fd: int, dst_name: str
) -> None:
    _require_exclusive_rename_support()
    ctypes.set_errno(0)
    result = _RENAME_FUNC(
        src_dir_fd,
        os.fsencode(src_name),
        dst_dir_fd,
        os.fsencode(dst_name),
        _RENAME_FLAG,
    )
    if result != 0:
        error = ctypes.get_errno() or errno.EIO
        raise OSError(error, os.strerror(error), dst_name)


def _stage_name(final: str) -> str:
    name = f".{final}.evidence-io-stage-{secrets.token_hex(16)}"
    if name in {".", ".."} or "/" in name or len(name.encode("utf-8")) > 255:
        raise PrepublicationError("staging name is unsafe")
    return name


def _name_exists(parent_fd: int, name: str) -> bool:
    try:
        os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return False
    return True


def _remove_owned_tree(
    parent_fd: int,
    name: str,
    owned: tuple[int, int],
    *,
    expected_file: bytes | None,
    expected_directory: dict[str, bytes] | None,
) -> bool:
    """Remove only an owned, still-exact staging tree; retain unknown entries."""

    try:
        info = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return True
    if _devino(info) != owned:
        return False
    if stat.S_ISREG(info.st_mode) and expected_file is not None:
        try:
            raw, identity = _read_regular_from_parent(
                parent_fd, name, max_bytes=DEFAULT_MAX_BYTES, expected_identity=owned
            )
        except (OSError, EvidenceIOError):
            return False
        if raw != expected_file or identity != owned:
            return False
        os.unlink(name, dir_fd=parent_fd)
        return True
    if not stat.S_ISDIR(info.st_mode) or expected_directory is None:
        return False
    child = os.open(
        name,
        os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW"),
        dir_fd=parent_fd,
    )
    try:
        if _devino(os.fstat(child)) != owned:
            return False
        expected_names = {path.split("/", 1)[0] for path in expected_directory}
        if set(os.listdir(child)) != expected_names:
            return False
        for entry in sorted(expected_names):
            entry_info = os.stat(entry, dir_fd=child, follow_symlinks=False)
            prefix = entry + "/"
            subtree = {
                path[len(prefix) :]: payload
                for path, payload in expected_directory.items()
                if path.startswith(prefix)
            }
            if subtree:
                removed = _remove_owned_tree(
                    child,
                    entry,
                    _devino(entry_info),
                    expected_file=None,
                    expected_directory=subtree,
                )
            else:
                removed = _remove_owned_tree(
                    child,
                    entry,
                    _devino(entry_info),
                    expected_file=expected_directory.get(entry),
                    expected_directory=None,
                )
            if not removed:
                return False
    finally:
        os.close(child)
    try:
        still = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return True
    if _devino(still) != owned:
        return False
    os.rmdir(name, dir_fd=parent_fd)
    return True


def _snapshot_file_payload(payload: object) -> bytes:
    if type(payload) not in (bytes, bytearray):
        raise PrepublicationError("publication payload is not bytes")
    if len(payload) > DEFAULT_MAX_BYTES:
        raise PrepublicationError("publication payload exceeds the byte bound")
    return bytes(payload)


def _snapshot_directory_payload(files: object) -> dict[str, bytes]:
    if type(files) is not dict:
        raise PrepublicationError("directory payload is unsafe")
    snapshot: dict[str, bytes] = {}
    total = 0
    count = 0
    for relative, payload in files.items():
        count += 1
        if count > MAX_DIRECTORY_ENTRIES:
            raise PrepublicationError("directory payload exceeds the entry bound")
        if type(relative) is not str:
            raise PrepublicationError("directory payload path is unsafe")
        if type(payload) not in (bytes, bytearray):
            raise PrepublicationError("publication payload is not bytes")
        size = len(payload)
        if size > DEFAULT_MAX_BYTES:
            raise PrepublicationError("publication payload exceeds the byte bound")
        path = "/".join(_relative_parts(relative))
        if path in snapshot:
            raise PrepublicationError("directory payload paths overlap")
        copied = bytes(payload)
        total += len(copied)
        if len(copied) > DEFAULT_MAX_BYTES or total > DEFAULT_MAX_BYTES:
            raise PrepublicationError("directory payload exceeds the byte bound")
        snapshot[path] = copied
    ordered = sorted(snapshot)
    for index, path in enumerate(ordered):
        prefix = path + "/"
        for other in ordered[index + 1 :]:
            if other.startswith(prefix):
                raise PrepublicationError("directory payload paths overlap")
    return snapshot


def _write_owned_file(parent_fd: int, name: str, payload: bytes) -> tuple[int, int]:
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | _flag("O_NOFOLLOW")
        | getattr(os, "O_NONBLOCK", 0)
    )
    descriptor = os.open(name, flags, 0o600, dir_fd=parent_fd)
    try:
        view = memoryview(payload)
        while view:
            count = os.write(descriptor, view)
            if count <= 0:
                raise PrepublicationError("short publication write")
            view = view[count:]
        os.fsync(descriptor)
        after = os.fstat(descriptor)
        if (
            not stat.S_ISREG(after.st_mode)
            or after.st_nlink != 1
            or after.st_size != len(payload)
        ):
            raise PrepublicationError("staging file is unsafe")
        return _devino(after)
    finally:
        os.close(descriptor)


def _mkdir_exclusive(parent_fd: int, name: str) -> tuple[int, int]:
    os.mkdir(name, 0o700, dir_fd=parent_fd)
    info = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise PrepublicationError("staging directory is unsafe")
    descriptor = _open_child_directory(parent_fd, name)
    try:
        active = os.fstat(descriptor)
        named = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if _devino(active) != _devino(info) or _devino(named) != _devino(info):
            raise PrepublicationError("staging directory changed during creation")
        return _devino(active)
    finally:
        os.close(descriptor)


def _ensure_owned_directory(parent_fd: int, name: str) -> int:
    try:
        os.mkdir(name, 0o700, dir_fd=parent_fd)
    except FileExistsError:
        info = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            raise PrepublicationError("staging parent is unsafe")
    return _open_child_directory(parent_fd, name)


def _write_directory_payload(stage_fd: int, files: dict[str, bytes]) -> None:
    for relative, payload in files.items():
        parts = _relative_parts(relative)
        parent_fd = stage_fd
        opened: list[int] = []
        try:
            for component in parts[:-1]:
                child = _ensure_owned_directory(parent_fd, component)
                opened.append(child)
                parent_fd = child
            _write_owned_file(parent_fd, parts[-1], payload)
        finally:
            _close_fds(opened)


def _directory_names(files: dict[str, bytes]) -> frozenset[str]:
    names: set[str] = set()
    for path in files:
        parts = path.split("/")
        names.update("/".join(parts[:index]) for index in range(1, len(parts)))
    return frozenset(names)


def _collect_inventory(
    dir_fd: int,
    expected_directories: frozenset[str],
    prefix: str = "",
) -> dict[str, bytes]:
    inventory: dict[str, bytes] = {}
    for name in os.listdir(dir_fd):
        if name in {".", ".."}:
            continue
        info = os.stat(name, dir_fd=dir_fd, follow_symlinks=False)
        relative = name if prefix == "" else f"{prefix}/{name}"
        _relative_parts(relative)
        if stat.S_ISLNK(info.st_mode):
            raise PrepublicationError("staging tree contains a symlink")
        if stat.S_ISDIR(info.st_mode):
            if relative not in expected_directories:
                raise PrepublicationError(
                    "staging tree contains an unexpected directory"
                )
            child = os.open(
                name,
                os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW"),
                dir_fd=dir_fd,
            )
            try:
                nested = _collect_inventory(child, expected_directories, relative)
            finally:
                os.close(child)
            inventory.update(nested)
        elif stat.S_ISREG(info.st_mode):
            raw, _ignored = _read_regular_from_parent(
                dir_fd, name, max_bytes=DEFAULT_MAX_BYTES
            )
            inventory[relative] = raw
        else:
            raise PrepublicationError("staging tree contains a special file")
        if len(inventory) > MAX_DIRECTORY_ENTRIES:
            raise PrepublicationError("staging tree exceeds the payload bound")
        if sum(len(item) for item in inventory.values()) > DEFAULT_MAX_BYTES:
            raise PrepublicationError("staging tree exceeds the payload bound")
    return inventory


def _fsync_directory_tree(
    dir_fd: int,
    expected_directories: frozenset[str],
    prefix: str = "",
) -> None:
    for name in os.listdir(dir_fd):
        info = os.stat(name, dir_fd=dir_fd, follow_symlinks=False)
        relative = name if not prefix else f"{prefix}/{name}"
        _relative_parts(relative)
        if stat.S_ISDIR(info.st_mode):
            if relative not in expected_directories:
                raise PrepublicationError(
                    "staging tree contains an unexpected directory"
                )
            child = _open_child_directory(dir_fd, name)
            try:
                _fsync_directory_tree(child, expected_directories, relative)
            finally:
                os.close(child)
        elif stat.S_ISREG(info.st_mode) and info.st_nlink == 1:
            _fsync_descriptor(dir_fd, name, directory=False)
        else:
            raise PrepublicationError("staging tree contains an unsafe leaf")
    os.fsync(dir_fd)


def _fsync_descriptor(parent_fd: int, name: str, *, directory: bool) -> None:
    flags = os.O_RDONLY | _flag("O_NOFOLLOW")
    if directory:
        flags |= _flag("O_DIRECTORY")
    descriptor = os.open(name, flags, dir_fd=parent_fd)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _inventory_digest(files: dict[str, bytes]) -> str:
    return _digest(
        canonical_json_bytes(
            {path: _digest(payload) for path, payload in files.items()}
        )
    )


def _verify_stage(
    parent_fd: int,
    name: str,
    *,
    kind: str,
    owned: tuple[int, int],
    expected_file: bytes | None,
    expected_directory: dict[str, bytes] | None,
) -> None:
    info = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    if _devino(info) != owned:
        raise PrepublicationError(
            "staging identity is not owned by this call",
            staging_retained=True,
        )
    if kind == "file":
        if not stat.S_ISREG(info.st_mode):
            raise PrepublicationError("staging file is unsafe", staging_retained=True)
        raw, identity = _read_regular_from_parent(
            parent_fd, name, max_bytes=DEFAULT_MAX_BYTES, expected_identity=owned
        )
        if identity != owned or raw != expected_file:
            raise PrepublicationError("staging bytes do not match the snapshot")
        return
    if not stat.S_ISDIR(info.st_mode):
        raise PrepublicationError("staging directory is unsafe", staging_retained=True)
    stage_fd = os.open(
        name,
        os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW"),
        dir_fd=parent_fd,
    )
    try:
        if _devino(os.fstat(stage_fd)) != owned:
            raise PrepublicationError(
                "staging identity is not owned by this call",
                staging_retained=True,
            )
        try:
            observed = _collect_inventory(
                stage_fd, _directory_names(expected_directory or {})
            )
        except UnsafePathError as error:
            raise PrepublicationError(
                "staging inventory cannot be read safely"
            ) from error
    finally:
        os.close(stage_fd)
    if observed != expected_directory:
        raise PrepublicationError("staging inventory does not match the snapshot")


def _verify_destination(
    parent_fd: int,
    name: str,
    *,
    relative_path: str,
    kind: str,
    owned: tuple[int, int],
    expected_file: bytes | None,
    expected_directory: dict[str, bytes] | None,
) -> str:
    info = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    if _devino(info) != owned:
        raise PostpublicationUncertainty(
            "published identity does not match the owned staging inode",
            relative_path=relative_path,
        )
    if kind == "file":
        raw, identity = _read_regular_from_parent(
            parent_fd, name, max_bytes=DEFAULT_MAX_BYTES, expected_identity=owned
        )
        if identity != owned or expected_file is None or raw != expected_file:
            raise PostpublicationUncertainty(
                "published bytes do not match the snapshot",
                relative_path=relative_path,
            )
        return _digest(raw)
    if not stat.S_ISDIR(info.st_mode):
        raise PostpublicationUncertainty(
            "published destination is not a directory",
            relative_path=relative_path,
        )
    dest_fd = os.open(
        name,
        os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW"),
        dir_fd=parent_fd,
    )
    try:
        if _devino(os.fstat(dest_fd)) != owned:
            raise PostpublicationUncertainty(
                "published identity does not match the owned staging inode",
                relative_path=relative_path,
            )
        observed = _collect_inventory(
            dest_fd, _directory_names(expected_directory or {})
        )
    finally:
        os.close(dest_fd)
    if expected_directory is None or observed != expected_directory:
        raise PostpublicationUncertainty(
            "published inventory does not match the snapshot",
            relative_path=relative_path,
        )
    return _inventory_digest(observed)


def _publish(
    root: Path,
    relative: str,
    *,
    kind: str,
    expected_file: bytes | None,
    expected_directory: dict[str, bytes] | None,
    populate: Callable[[int, str], tuple[int, int]],
    fault_hook: FaultHook | None,
) -> PublicationReceipt:
    _require_exclusive_rename_support()
    root = _canonical_root(Path(root))
    parts = _relative_parts(relative)
    descriptors, identities = _open_chain(root, parts[:-1])
    parent_fd = descriptors[-1]
    final = parts[-1]
    stage = _stage_name(final)
    published = False
    created = False
    owned: tuple[int, int] | None = None

    def fault(cut: str) -> None:
        if fault_hook is not None:
            fault_hook(cut)

    try:
        if _name_exists(parent_fd, final):
            raise PrepublicationError("destination exists")
        if _name_exists(parent_fd, stage):
            raise PrepublicationError(
                "staging name is not exclusively available",
                staging_retained=True,
            )
        fault("before_write")
        owned = populate(parent_fd, stage)
        created = True
        fault("after_write")
        if kind == "file":
            _fsync_descriptor(parent_fd, stage, directory=False)
        else:
            stage_fd = os.open(
                stage,
                os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW"),
                dir_fd=parent_fd,
            )
            try:
                _fsync_directory_tree(
                    stage_fd, _directory_names(expected_directory or {})
                )
            finally:
                os.close(stage_fd)
        fault("after_file_fsync")
        fault("before_publish")
        _revalidate_chain(root, parts[:-1], identities, descriptors)
        _verify_stage(
            parent_fd,
            stage,
            kind=kind,
            owned=owned,
            expected_file=expected_file,
            expected_directory=expected_directory,
        )
        try:
            _exclusive_rename(parent_fd, stage, parent_fd, final)
        except FileExistsError as error:
            raise PrepublicationError("destination exists") from error
        published = True
        fault("after_publish")
        os.fsync(parent_fd)
        fault("after_parent_fsync")
        _revalidate_chain(root, parts[:-1], identities, descriptors)
        payload_sha256 = _verify_destination(
            parent_fd,
            final,
            relative_path=relative,
            kind=kind,
            owned=owned,
            expected_file=expected_file,
            expected_directory=expected_directory,
        )
        _revalidate_chain(root, parts[:-1], identities, descriptors)
        return PublicationReceipt(
            relative_path=relative,
            kind=kind,
            payload_sha256=payload_sha256,
        )
    except PostpublicationUncertainty:
        raise
    except Exception as error:
        if published:
            raise PostpublicationUncertainty(
                "destination was published but a later durability step failed",
                relative_path=relative,
            ) from error
        retained = True
        if created and owned is not None:
            try:
                retained = not _remove_owned_tree(
                    parent_fd,
                    stage,
                    owned,
                    expected_file=expected_file,
                    expected_directory=expected_directory,
                )
            except (OSError, EvidenceIOError):
                retained = True
        elif not created:
            retained = _name_exists(parent_fd, stage)
        else:
            retained = True
        if isinstance(error, EvidenceIOError):
            if isinstance(error, PrepublicationError):
                error.staging_retained = bool(error.staging_retained or retained)
            raise
        raise PrepublicationError(
            "publication failed before the destination was visible",
            staging_retained=retained,
        ) from error
    finally:
        _close_fds(descriptors)


def publish_exclusive_file(
    root: Path,
    relative: str,
    payload: bytes,
    *,
    _fault_hook: FaultHook | None = None,
) -> PublicationReceipt:
    """Publish exact bytes to a new regular file without replacing a target.

    The payload is snapshotted before staging.  Staging is created exclusively
    by this call; a colliding name is left untouched.  After pre-publish hooks
    the owned inode and exact bytes are revalidated, then independently
    revalidated at the destination before a success receipt is returned.
    """

    snapshot = _snapshot_file_payload(payload)

    def populate(parent_fd: int, stage: str) -> tuple[int, int]:
        return _write_owned_file(parent_fd, stage, snapshot)

    return _publish(
        root,
        relative,
        kind="file",
        expected_file=snapshot,
        expected_directory=None,
        populate=populate,
        fault_hook=_fault_hook,
    )


def publish_exclusive_directory(
    root: Path,
    relative: str,
    files: dict[str, bytes | bytearray],
    *,
    _fault_hook: FaultHook | None = None,
) -> PublicationReceipt:
    """Publish a new directory of regular files without replacing a target.

    Nested paths are snapshotted into plain ``bytes`` with entry and total
    bounds before any staging.  Extra, mutated, or replaced leaves after
    pre-publish hooks fail closed.  Platforms without exclusive no-replace
    rename fail before the destination becomes visible.
    """

    snapshot = _snapshot_directory_payload(files)

    def populate(parent_fd: int, stage: str) -> tuple[int, int]:
        owned = _mkdir_exclusive(parent_fd, stage)
        stage_fd = os.open(
            stage,
            os.O_RDONLY | _flag("O_DIRECTORY") | _flag("O_NOFOLLOW"),
            dir_fd=parent_fd,
        )
        try:
            if _devino(os.fstat(stage_fd)) != owned:
                raise PrepublicationError(
                    "staging identity is not owned by this call",
                    staging_retained=True,
                )
            _write_directory_payload(stage_fd, snapshot)
        finally:
            os.close(stage_fd)
        return owned

    return _publish(
        root,
        relative,
        kind="directory",
        expected_file=None,
        expected_directory=snapshot,
        populate=populate,
        fault_hook=_fault_hook,
    )
