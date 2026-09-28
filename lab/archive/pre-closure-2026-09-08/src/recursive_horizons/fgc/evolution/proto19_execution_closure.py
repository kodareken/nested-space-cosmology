"""Standard-library execution-source closure for a later PRO19 bootstrap.

This module is deliberately independent of the Recursive Horizons package and
of every numerical dependency.  It can be loaded from bytes obtained with
``git show`` and run under ``python -I -B`` before a repository package is
imported.  Its only responsibilities are to authenticate an implementation
image and to constrain subsequent repository-local imports.  It does not
authorize execution, inspect a campaign, acquire a writer, or mutate state.

The closure uses a two-commit model.  ``implementation_commit`` owns every
executable/source/input blob.  A later externally supplied authority commit may
change exactly ``authority_delta_paths`` (normally prospective config, result,
documentation, and manifest files), but it may not change Python/importable
source.  Every bound path must retain the implementation blob at the authority
commit and in the live worktree.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from hashlib import sha256
import importlib
import importlib.abc
import importlib.machinery
import json
import keyword
import os
from pathlib import Path, PurePosixPath
import platform
import stat
import subprocess
import sys
import sysconfig
from types import ModuleType
from typing import Iterable, Iterator, Mapping, Sequence


SCHEMA = "FGC-1-PRO19-execution-source-closure-v1"
FILE_KINDS = frozenset({"python", "source", "input"})
IMPORT_ROOTS = ("src", "")
_SOURCE_SUFFIXES = (".py", ".pyi", ".pyw", ".pyc", ".pyo") + tuple(
    importlib.machinery.EXTENSION_SUFFIXES
)
_HEX = frozenset("0123456789abcdef")
_IDENTITY_FIELDS = frozenset(
    {
        "implementation",
        "implementation_version",
        "cache_tag",
        "hexversion",
        "executable",
        "executable_realpath",
        "executable_sha256",
        "prefix",
        "base_prefix",
        "exec_prefix",
        "base_exec_prefix",
        "platform",
        "platform_release",
        "machine",
        "byteorder",
        "filesystem_encoding",
        "filesystem_errors",
        "sysconfig_platform",
        "soabi",
        "flags",
        "base_sys_path",
    }
)


class ExecutionClosureError(RuntimeError):
    """The prospective closure or observed execution image is not exact."""


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise ExecutionClosureError("closure record is not canonicalizable") from exc


def _sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise ExecutionClosureError(f"{label} is not a sequence")
    return value


def _hex(value: object, length: int, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != length
        or value.lower() != value
        or any(character not in _HEX for character in value)
    ):
        raise ExecutionClosureError(f"{label} is not lowercase hexadecimal")
    return value


def _relative(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value or "\x00" in value:
        raise ExecutionClosureError(f"{label} is not a portable repository path")
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or path.as_posix() != value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ExecutionClosureError(f"{label} escapes or aliases the repository")
    return value


def _prefix_path(value: object, label: str) -> str:
    # Prefixes may intentionally end in a partial file stem, but still must be
    # rooted in a canonical repository directory.
    path = _relative(value, label)
    if path.endswith("/"):
        raise ExecutionClosureError(f"{label} is noncanonical")
    return path


def _module(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ExecutionClosureError(f"{label} is absent")
    parts = value.split(".")
    if any(not part.isidentifier() or keyword.iskeyword(part) for part in parts):
        raise ExecutionClosureError(f"{label} is not a canonical module name")
    return value


def _matches_module(name: str, prefix: str) -> bool:
    return name == prefix or name.startswith(prefix + ".")


def _matches_path(path: str, prefix: str) -> bool:
    return path == prefix or path.startswith(prefix)


def _is_importable_source(path: str) -> bool:
    lowered = path.lower()
    return any(lowered.endswith(suffix.lower()) for suffix in _SOURCE_SUFFIXES)


@dataclass(frozen=True, slots=True)
class InterpreterIdentity:
    """Stable, non-secret identity observed before repository imports."""

    implementation: str
    implementation_version: str
    cache_tag: str
    hexversion: int
    executable: str
    executable_realpath: str
    executable_sha256: str
    prefix: str
    base_prefix: str
    exec_prefix: str
    base_exec_prefix: str
    platform: str
    platform_release: str
    machine: str
    byteorder: str
    filesystem_encoding: str
    filesystem_errors: str
    sysconfig_platform: str
    soabi: str
    flags: tuple[tuple[str, int | bool], ...]
    base_sys_path: tuple[str, ...]

    def to_mapping(self) -> dict[str, object]:
        return {
            "implementation": self.implementation,
            "implementation_version": self.implementation_version,
            "cache_tag": self.cache_tag,
            "hexversion": self.hexversion,
            "executable": self.executable,
            "executable_realpath": self.executable_realpath,
            "executable_sha256": self.executable_sha256,
            "prefix": self.prefix,
            "base_prefix": self.base_prefix,
            "exec_prefix": self.exec_prefix,
            "base_exec_prefix": self.base_exec_prefix,
            "platform": self.platform,
            "platform_release": self.platform_release,
            "machine": self.machine,
            "byteorder": self.byteorder,
            "filesystem_encoding": self.filesystem_encoding,
            "filesystem_errors": self.filesystem_errors,
            "sysconfig_platform": self.sysconfig_platform,
            "soabi": self.soabi,
            "flags": [[name, value] for name, value in self.flags],
            "base_sys_path": list(self.base_sys_path),
        }

    @classmethod
    def from_mapping(cls, value: object) -> "InterpreterIdentity":
        if not isinstance(value, Mapping) or set(value) != _IDENTITY_FIELDS:
            raise ExecutionClosureError("interpreter identity fields differ")
        flags: list[tuple[str, int | bool]] = []
        for item in _sequence(value["flags"], "interpreter flags"):
            pair = _sequence(item, "interpreter flag")
            if (
                len(pair) != 2
                or not isinstance(pair[0], str)
                or not isinstance(pair[1], (int, bool))
            ):
                raise ExecutionClosureError("interpreter flag record differs")
            flags.append((pair[0], pair[1]))
        paths = _sequence(value["base_sys_path"], "base sys.path")
        if any(not isinstance(item, str) or not item for item in paths):
            raise ExecutionClosureError("base sys.path record differs")
        text_fields = _IDENTITY_FIELDS - {"hexversion", "flags", "base_sys_path"}
        if any(not isinstance(value[name], str) for name in text_fields):
            raise ExecutionClosureError("interpreter text identity differs")
        if not isinstance(value["hexversion"], int):
            raise ExecutionClosureError("interpreter hexversion differs")
        identity = cls(
            **{name: value[name] for name in text_fields},  # type: ignore[arg-type]
            hexversion=value["hexversion"],
            flags=tuple(flags),
            base_sys_path=tuple(paths),  # type: ignore[arg-type]
        )
        if tuple(sorted(identity.flags)) != identity.flags:
            raise ExecutionClosureError("interpreter flags are not canonical")
        _hex(identity.executable_sha256, 64, "interpreter executable SHA-256")
        return identity


@dataclass(frozen=True, slots=True)
class FilePin:
    path: str
    kind: str
    git_blob_oid: str
    sha256: str

    def to_mapping(self) -> dict[str, str]:
        return {
            "path": self.path,
            "kind": self.kind,
            "git_blob_oid": self.git_blob_oid,
            "sha256": self.sha256,
        }


@dataclass(frozen=True, slots=True)
class ImportPin:
    module: str
    path: str
    sha256: str

    def to_mapping(self) -> dict[str, str]:
        return {"module": self.module, "path": self.path, "sha256": self.sha256}


@dataclass(frozen=True, slots=True)
class NamespacePin:
    module: str
    path: str

    def to_mapping(self) -> dict[str, str]:
        return {"module": self.module, "path": self.path}


@dataclass(frozen=True, slots=True)
class ExecutionClosureRecord:
    """Prospective content record; it is not an execution authorization."""

    implementation_commit: str
    git_object_format: str
    files: tuple[FilePin, ...]
    interpreter: InterpreterIdentity
    imports: tuple[ImportPin, ...]
    namespaces: tuple[NamespacePin, ...]
    forbidden_module_prefixes: tuple[str, ...]
    forbidden_script_prefixes: tuple[str, ...]
    authority_delta_paths: tuple[str, ...]
    schema: str = SCHEMA

    def to_mapping(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "implementation_commit": self.implementation_commit,
            "git_object_format": self.git_object_format,
            "files": [item.to_mapping() for item in self.files],
            "interpreter": self.interpreter.to_mapping(),
            "imports": [item.to_mapping() for item in self.imports],
            "namespaces": [item.to_mapping() for item in self.namespaces],
            "forbidden_module_prefixes": list(self.forbidden_module_prefixes),
            "forbidden_script_prefixes": list(self.forbidden_script_prefixes),
            "authority_delta_paths": list(self.authority_delta_paths),
        }

    @property
    def canonical_sha256(self) -> str:
        return sha256(_canonical(self.to_mapping())).hexdigest()


@dataclass(frozen=True, slots=True)
class CapturedFile:
    pin: FilePin
    content: bytes
    device: int
    inode: int


@dataclass(frozen=True, slots=True)
class VerifiedExecutionClosure:
    """Immutable byte capture returned by the read-only pre-import verifier."""

    root: Path
    authority_commit: str
    record: ExecutionClosureRecord
    files: tuple[CapturedFile, ...]
    import_path: tuple[str, str]
    require_detached: bool

    def bytes_for_path(self, path: str) -> bytes:
        for item in self.files:
            if item.pin.path == path:
                return item.content
        raise ExecutionClosureError(f"path is outside the verified closure: {path}")


def parse_execution_closure_record(value: object) -> ExecutionClosureRecord:
    """Strictly parse the future config/result mapping without reading a repo."""
    fields = {
        "schema",
        "implementation_commit",
        "git_object_format",
        "files",
        "interpreter",
        "imports",
        "namespaces",
        "forbidden_module_prefixes",
        "forbidden_script_prefixes",
        "authority_delta_paths",
    }
    if (
        not isinstance(value, Mapping)
        or set(value) != fields
        or value.get("schema") != SCHEMA
    ):
        raise ExecutionClosureError("execution closure record fields or schema differ")
    object_format = value["git_object_format"]
    if object_format not in {"sha1", "sha256"}:
        raise ExecutionClosureError("Git object format differs")
    oid_length = 40 if object_format == "sha1" else 64
    files: list[FilePin] = []
    for raw in _sequence(value["files"], "closure files"):
        if not isinstance(raw, Mapping) or set(raw) != {
            "path",
            "kind",
            "git_blob_oid",
            "sha256",
        }:
            raise ExecutionClosureError("closure file pin fields differ")
        path = _relative(raw["path"], "closure file")
        kind = raw["kind"]
        if kind not in FILE_KINDS:
            raise ExecutionClosureError("closure file kind differs")
        files.append(
            FilePin(
                path=path,
                kind=str(kind),
                git_blob_oid=_hex(raw["git_blob_oid"], oid_length, f"Git blob {path}"),
                sha256=_hex(raw["sha256"], 64, f"live SHA-256 {path}"),
            )
        )
    imports: list[ImportPin] = []
    for raw in _sequence(value["imports"], "closure imports"):
        if not isinstance(raw, Mapping) or set(raw) != {"module", "path", "sha256"}:
            raise ExecutionClosureError("import pin fields differ")
        imports.append(
            ImportPin(
                module=_module(raw["module"], "import module"),
                path=_relative(raw["path"], "import origin"),
                sha256=_hex(raw["sha256"], 64, "import SHA-256"),
            )
        )
    namespaces: list[NamespacePin] = []
    for raw in _sequence(value["namespaces"], "closure namespaces"):
        if not isinstance(raw, Mapping) or set(raw) != {"module", "path"}:
            raise ExecutionClosureError("namespace pin fields differ")
        namespaces.append(
            NamespacePin(
                module=_module(raw["module"], "namespace module"),
                path=_relative(raw["path"], "namespace path"),
            )
        )

    def prefixes(name: str, *, modules: bool) -> tuple[str, ...]:
        output: list[str] = []
        for item in _sequence(value[name], name):
            output.append(_module(item, name) if modules else _prefix_path(item, name))
        return tuple(output)

    record = ExecutionClosureRecord(
        implementation_commit=str(value["implementation_commit"]),
        git_object_format=str(object_format),
        files=tuple(files),
        interpreter=InterpreterIdentity.from_mapping(value["interpreter"]),
        imports=tuple(imports),
        namespaces=tuple(namespaces),
        forbidden_module_prefixes=prefixes("forbidden_module_prefixes", modules=True),
        forbidden_script_prefixes=prefixes("forbidden_script_prefixes", modules=False),
        authority_delta_paths=tuple(
            _relative(item, "authority delta path")
            for item in _sequence(
                value["authority_delta_paths"], "authority delta paths"
            )
        ),
    )
    _validate_record(record)
    return record


def _validate_record(record: ExecutionClosureRecord) -> None:
    if record.schema != SCHEMA or record.git_object_format not in {"sha1", "sha256"}:
        raise ExecutionClosureError("closure record identity differs")
    oid_length = 40 if record.git_object_format == "sha1" else 64
    _hex(record.implementation_commit, oid_length, "implementation commit")
    if (
        InterpreterIdentity.from_mapping(record.interpreter.to_mapping())
        != record.interpreter
    ):
        raise ExecutionClosureError("interpreter identity is noncanonical")
    if not record.files:
        raise ExecutionClosureError("execution closure has no files")
    if tuple(sorted(record.files, key=lambda item: item.path)) != record.files:
        raise ExecutionClosureError("closure files are not canonical")
    if tuple(sorted(record.imports, key=lambda item: item.module)) != record.imports:
        raise ExecutionClosureError("import pins are not canonical")
    if (
        tuple(sorted(record.namespaces, key=lambda item: item.module))
        != record.namespaces
    ):
        raise ExecutionClosureError("namespace pins are not canonical")
    for values, label in (
        (record.forbidden_module_prefixes, "forbidden module prefixes"),
        (record.forbidden_script_prefixes, "forbidden script prefixes"),
        (record.authority_delta_paths, "authority delta paths"),
    ):
        if tuple(sorted(set(values))) != values:
            raise ExecutionClosureError(f"{label} are duplicated or noncanonical")
    for prefix in record.forbidden_module_prefixes:
        _module(prefix, "forbidden module prefix")
    for prefix in record.forbidden_script_prefixes:
        _prefix_path(prefix, "forbidden script prefix")
    file_map = {item.path: item for item in record.files}
    if len(file_map) != len(record.files):
        raise ExecutionClosureError("closure file path is duplicated")
    if set(record.authority_delta_paths) & set(file_map):
        raise ExecutionClosureError("authority delta overlaps an implementation file")
    import_map = {item.module: item for item in record.imports}
    if len(import_map) != len(record.imports) or len(
        {item.path for item in record.imports}
    ) != len(record.imports):
        raise ExecutionClosureError("import name or origin is duplicated")
    namespace_map = {item.module: item for item in record.namespaces}
    if len(namespace_map) != len(record.namespaces) or set(namespace_map) & set(
        import_map
    ):
        raise ExecutionClosureError("namespace name is duplicated or executable")
    for item in record.files:
        _relative(item.path, "closure file")
        if item.kind not in FILE_KINDS:
            raise ExecutionClosureError("closure file kind differs")
        _hex(item.git_blob_oid, oid_length, "closure Git blob")
        _hex(item.sha256, 64, "closure SHA-256")
        if any(
            _matches_path(item.path, prefix)
            for prefix in record.forbidden_script_prefixes
        ):
            raise ExecutionClosureError("forbidden script path is present in closure")
    for item in record.imports:
        _module(item.module, "import module")
        pin = file_map.get(item.path)
        if (
            pin is None
            or pin.kind != "python"
            or pin.sha256 != item.sha256
            or not item.path.endswith(".py")
        ):
            raise ExecutionClosureError(
                "import origin is not its exact Python file pin"
            )
        if any(
            _matches_module(item.module, prefix)
            for prefix in record.forbidden_module_prefixes
        ):
            raise ExecutionClosureError("forbidden module is present in import map")
        parts = item.module.split(".")
        for length in range(1, len(parts)):
            parent = ".".join(parts[:length])
            parent_import = import_map.get(parent)
            if parent not in namespace_map and (
                parent_import is None or not parent_import.path.endswith("/__init__.py")
            ):
                raise ExecutionClosureError(
                    "import parent lacks a pinned package or namespace"
                )
    for item in record.namespaces:
        _module(item.module, "namespace module")
        expected_suffix = item.module.replace(".", "/")
        if item.path not in {expected_suffix, f"src/{expected_suffix}"}:
            raise ExecutionClosureError("namespace path does not match its module")
        if any(
            _matches_module(item.module, prefix)
            for prefix in record.forbidden_module_prefixes
        ):
            raise ExecutionClosureError("forbidden namespace is present")
    for path in record.authority_delta_paths:
        _relative(path, "authority delta path")
        if _is_importable_source(path):
            raise ExecutionClosureError(
                "authority delta may not contain importable source"
            )
        if any(
            _matches_path(path, prefix) for prefix in record.forbidden_script_prefixes
        ):
            raise ExecutionClosureError(
                "authority delta contains a forbidden script path"
            )


def _repository_root(value: Path) -> Path:
    supplied = Path(os.path.abspath(os.fspath(value)))
    try:
        metadata = supplied.lstat()
    except OSError as exc:
        raise ExecutionClosureError("repository root is absent") from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise ExecutionClosureError("repository root is not a direct directory")
    return Path(os.path.realpath(supplied))


def _git(
    root: Path, *arguments: str, allowed: frozenset[int] = frozenset({0})
) -> subprocess.CompletedProcess[bytes]:
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"})
    try:
        result = subprocess.run(
            (
                "git",
                "-C",
                str(root),
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                *arguments,
            ),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            env=environment,
        )
    except OSError as exc:
        raise ExecutionClosureError("Git authority operation could not start") from exc
    if result.returncode not in allowed:
        raise ExecutionClosureError(
            f"Git authority operation failed: {' '.join(arguments)}"
        )
    return result


def _git_text(root: Path, *arguments: str) -> str:
    try:
        return _git(root, *arguments).stdout.decode("utf-8").strip()
    except UnicodeDecodeError as exc:
        raise ExecutionClosureError("Git authority output is not UTF-8") from exc


def _require_isolated_flags() -> None:
    required = {
        "isolated": 1,
        "ignore_environment": 1,
        "no_user_site": 1,
        "dont_write_bytecode": 1,
        "safe_path": True,
    }
    if any(
        getattr(sys.flags, name, None) != expected
        for name, expected in required.items()
    ):
        raise ExecutionClosureError("bootstrap requires python -I -B")
    if any(not item for item in sys.path):
        raise ExecutionClosureError(
            "isolated base sys.path contains the current directory"
        )


def _require_nofollow_support() -> None:
    if (
        not hasattr(os, "O_NOFOLLOW")
        or not hasattr(os, "O_DIRECTORY")
        or os.open not in os.supports_dir_fd
        or os.stat not in os.supports_dir_fd
        or os.stat not in os.supports_follow_symlinks
    ):
        raise ExecutionClosureError(
            "host lacks descriptor-relative no-follow filesystem support"
        )


def _absolute_regular_sha256(path: str) -> str:
    descriptor = -1
    try:
        descriptor = os.open(
            path,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ExecutionClosureError("interpreter executable is not a regular file")
        digest = sha256()
        while block := os.read(descriptor, 1 << 20):
            digest.update(block)
        return digest.hexdigest()
    except OSError as exc:
        raise ExecutionClosureError(
            "interpreter executable cannot be hashed safely"
        ) from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)


def current_interpreter_identity() -> InterpreterIdentity:
    """Capture the isolated interpreter before repository paths are installed."""
    _require_isolated_flags()
    _require_nofollow_support()
    executable = os.path.abspath(sys.executable)
    executable_realpath = os.path.realpath(executable)
    flag_values = tuple(
        sorted(
            (name, value)
            for name in dir(sys.flags)
            if not name.startswith("_")
            and isinstance((value := getattr(sys.flags, name)), (int, bool))
        )
    )
    return InterpreterIdentity(
        implementation=sys.implementation.name,
        implementation_version=platform.python_version(),
        cache_tag=sys.implementation.cache_tag or "",
        hexversion=sys.hexversion,
        executable=executable,
        executable_realpath=executable_realpath,
        executable_sha256=_absolute_regular_sha256(executable_realpath),
        prefix=os.path.realpath(sys.prefix),
        base_prefix=os.path.realpath(sys.base_prefix),
        exec_prefix=os.path.realpath(sys.exec_prefix),
        base_exec_prefix=os.path.realpath(sys.base_exec_prefix),
        platform=sys.platform,
        platform_release=platform.release(),
        machine=platform.machine(),
        byteorder=sys.byteorder,
        filesystem_encoding=sys.getfilesystemencoding(),
        filesystem_errors=sys.getfilesystemencodeerrors(),
        sysconfig_platform=sysconfig.get_platform(),
        soabi=str(sysconfig.get_config_var("SOABI") or ""),
        flags=flag_values,
        base_sys_path=tuple(sys.path),
    )


def _nofollow_stat(root: Path, relative: str) -> os.stat_result | None:
    parts = PurePosixPath(_relative(relative, "repository path")).parts
    flags_dir = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_DIRECTORY", 0)
    )
    opened: list[int] = []
    try:
        parent = os.open(root, flags_dir)
        opened.append(parent)
        for part in parts[:-1]:
            try:
                before = os.stat(part, dir_fd=parent, follow_symlinks=False)
            except FileNotFoundError:
                return None
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise ExecutionClosureError("repository path has an unsafe ancestor")
            child = os.open(part, flags_dir, dir_fd=parent)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                raise ExecutionClosureError("repository path changed during traversal")
            opened.append(child)
            parent = child
        try:
            return os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            return None
    except ExecutionClosureError:
        raise
    except OSError as exc:
        raise ExecutionClosureError(
            "repository path cannot be inspected safely"
        ) from exc
    finally:
        for descriptor in reversed(opened):
            os.close(descriptor)


def _nofollow_regular_bytes(root: Path, relative: str) -> tuple[bytes, os.stat_result]:
    parts = PurePosixPath(_relative(relative, "closure path")).parts
    flags_dir = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_DIRECTORY", 0)
    )
    flags_file = (
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    )
    opened: list[int] = []
    leaf = -1
    try:
        parent = os.open(root, flags_dir)
        opened.append(parent)
        for part in parts[:-1]:
            before = os.stat(part, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise ExecutionClosureError("closure path has an unsafe ancestor")
            child = os.open(part, flags_dir, dir_fd=parent)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                raise ExecutionClosureError("closure path changed during traversal")
            opened.append(child)
            parent = child
        before = os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
            raise ExecutionClosureError("closure path is not a regular file")
        leaf = os.open(parts[-1], flags_file, dir_fd=parent)
        actual = os.fstat(leaf)
        if (before.st_dev, before.st_ino) != (actual.st_dev, actual.st_ino):
            raise ExecutionClosureError("closure file changed before read")
        if actual.st_nlink != 1:
            raise ExecutionClosureError(
                "closure file has non-unique hard-link ownership"
            )
        chunks: list[bytes] = []
        while block := os.read(leaf, 1 << 20):
            chunks.append(block)
        content = b"".join(chunks)
        after = os.fstat(leaf)
        stable = (before.st_size, before.st_mtime_ns, before.st_ctime_ns)
        if len(content) != before.st_size or stable != (
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise ExecutionClosureError("closure file changed during read")
        return content, after
    except ExecutionClosureError:
        raise
    except (FileNotFoundError, OSError) as exc:
        raise ExecutionClosureError("closure path cannot be read safely") from exc
    finally:
        if leaf != -1:
            os.close(leaf)
        for descriptor in reversed(opened):
            os.close(descriptor)


def _tree_blob(root: Path, commit: str, relative: str) -> tuple[str, bytes]:
    entries = [
        entry
        for entry in _git(root, "ls-tree", "-z", commit, "--", relative).stdout.split(
            b"\0"
        )
        if entry
    ]
    if len(entries) != 1:
        raise ExecutionClosureError(
            f"closure path is not one Git tree entry: {relative}"
        )
    try:
        identity, raw_path = entries[0].split(b"\t", 1)
        mode, kind, oid = identity.split(b" ", 2)
        observed_path = raw_path.decode("utf-8")
        oid_text = oid.decode("ascii")
    except (ValueError, UnicodeDecodeError) as exc:
        raise ExecutionClosureError("Git tree entry is malformed") from exc
    if (
        mode not in {b"100644", b"100755"}
        or kind != b"blob"
        or observed_path != relative
    ):
        raise ExecutionClosureError(
            f"closure path is not a regular Git blob: {relative}"
        )
    content = _git(root, "cat-file", "blob", oid_text).stdout
    return oid_text, content


def _require_clean(root: Path) -> None:
    for arguments in (
        ("diff", "--quiet", "--no-ext-diff", "--ignore-submodules=none", "--"),
        (
            "diff",
            "--cached",
            "--quiet",
            "--no-ext-diff",
            "--ignore-submodules=none",
            "--",
        ),
    ):
        result = _git(root, *arguments, allowed=frozenset({0, 1}))
        if result.returncode == 1:
            raise ExecutionClosureError(
                "tracked working tree or index drift is present"
            )


def _resolve_commit(root: Path, commit: str, object_format: str, label: str) -> str:
    length = 40 if object_format == "sha1" else 64
    commit = _hex(commit, length, label)
    resolved = _git_text(root, "rev-parse", "--verify", f"{commit}^{{commit}}")
    if resolved != commit:
        raise ExecutionClosureError(f"{label} does not resolve exactly")
    return commit


def _changed_paths(root: Path, older: str, newer: str) -> tuple[str, ...]:
    raw = _git(
        root, "diff", "--no-renames", "--name-only", "-z", older, newer, "--"
    ).stdout
    paths: list[str] = []
    for item in raw.split(b"\0"):
        if not item:
            continue
        try:
            paths.append(_relative(item.decode("utf-8"), "authority changed path"))
        except UnicodeDecodeError as exc:
            raise ExecutionClosureError(
                "authority delta contains a non-UTF-8 path"
            ) from exc
    return tuple(sorted(set(paths)))


def _require_repository(
    root: Path,
    record: ExecutionClosureRecord,
    authority_commit: str,
    *,
    require_detached: bool,
) -> None:
    top = Path(os.path.realpath(_git_text(root, "rev-parse", "--show-toplevel")))
    if top != root:
        raise ExecutionClosureError("supplied path is not the Git repository root")
    observed_format = _git_text(root, "rev-parse", "--show-object-format")
    if observed_format != record.git_object_format:
        raise ExecutionClosureError("live Git object format differs")
    implementation = _resolve_commit(
        root, record.implementation_commit, observed_format, "implementation commit"
    )
    authority = _resolve_commit(
        root, authority_commit, observed_format, "authority commit"
    )
    if _git_text(root, "rev-parse", "HEAD") != authority:
        raise ExecutionClosureError(
            "current HEAD differs from the external authority commit"
        )
    ancestor = _git(
        root,
        "merge-base",
        "--is-ancestor",
        implementation,
        authority,
        allowed=frozenset({0, 1}),
    )
    if ancestor.returncode != 0:
        raise ExecutionClosureError(
            "implementation commit is not an authority ancestor"
        )
    if require_detached:
        symbolic = _git(root, "symbolic-ref", "-q", "HEAD", allowed=frozenset({0, 1}))
        if symbolic.returncode == 0:
            raise ExecutionClosureError("detached authority image was required")
    _require_clean(root)
    changed = _changed_paths(root, implementation, authority)
    if changed != record.authority_delta_paths:
        raise ExecutionClosureError("implementation-to-authority changed paths differ")
    if any(_is_importable_source(path) for path in changed):
        raise ExecutionClosureError("authority descendant contains later source drift")


def _capture_files(
    root: Path, record: ExecutionClosureRecord, *, authority_commit: str | None
) -> tuple[CapturedFile, ...]:
    output: list[CapturedFile] = []
    for pin in record.files:
        implementation_oid, committed = _tree_blob(
            root, record.implementation_commit, pin.path
        )
        if (
            implementation_oid != pin.git_blob_oid
            or sha256(committed).hexdigest() != pin.sha256
        ):
            raise ExecutionClosureError(f"implementation blob pin differs: {pin.path}")
        if authority_commit is not None:
            authority_oid, authority_bytes = _tree_blob(
                root, authority_commit, pin.path
            )
            if authority_oid != implementation_oid or authority_bytes != committed:
                raise ExecutionClosureError(
                    f"authority commit changed a closure blob: {pin.path}"
                )
        live, metadata = _nofollow_regular_bytes(root, pin.path)
        if live != committed or sha256(live).hexdigest() != pin.sha256:
            raise ExecutionClosureError(f"live closure bytes differ: {pin.path}")
        output.append(
            CapturedFile(
                pin=pin, content=live, device=metadata.st_dev, inode=metadata.st_ino
            )
        )
    identities = {(item.device, item.inode) for item in output}
    if len(identities) != len(output):
        raise ExecutionClosureError("closure paths alias one physical inode")
    return tuple(output)


def _module_candidates(name: str) -> tuple[str, ...]:
    stem = name.replace(".", "/")
    values: list[str] = []
    for root in IMPORT_ROOTS:
        prefix = f"{root}/" if root else ""
        values.extend((f"{prefix}{stem}.py", f"{prefix}{stem}/__init__.py"))
    return tuple(values)


def _namespace_candidates(name: str) -> tuple[str, ...]:
    stem = name.replace(".", "/")
    return tuple(f"{root}/{stem}" if root else stem for root in IMPORT_ROOTS)


def _candidate_exists(root: Path, relative: str) -> bool:
    metadata = _nofollow_stat(root, relative)
    if metadata is None:
        return False
    if stat.S_ISLNK(metadata.st_mode):
        raise ExecutionClosureError(f"import candidate is a symlink: {relative}")
    if not (stat.S_ISREG(metadata.st_mode) or stat.S_ISDIR(metadata.st_mode)):
        raise ExecutionClosureError(f"import candidate is a special file: {relative}")
    return True


def _repository_candidate_paths(root: Path, name: str) -> tuple[str, ...]:
    candidates = (*_module_candidates(name), *_namespace_candidates(name))
    return tuple(path for path in candidates if _candidate_exists(root, path))


def _git_inventory(root: Path, *, ignored: bool) -> tuple[str, ...]:
    arguments = ["ls-files", "--others", "-z", "--exclude-standard"]
    if ignored:
        arguments.insert(2, "--ignored")
    raw = _git(root, *arguments, "--").stdout
    values: list[str] = []
    for item in raw.split(b"\0"):
        if not item:
            continue
        try:
            values.append(_relative(item.decode("utf-8"), "untracked import inventory"))
        except UnicodeDecodeError as exc:
            raise ExecutionClosureError(
                "untracked import inventory is not UTF-8"
            ) from exc
    return tuple(sorted(set(values)))


def _path_shadows_module(path: str, prefix: str) -> bool:
    if "__pycache__" in PurePosixPath(path).parts or not _is_importable_source(path):
        return False
    stem = prefix.replace(".", "/")
    for root in IMPORT_ROOTS:
        base = f"{root}/{stem}" if root else stem
        if path == base or path.startswith(base + ".") or path.startswith(base + "/"):
            return True
    return False


def _require_no_untracked_shadows(root: Path, record: ExecutionClosureRecord) -> None:
    protected = {
        *(item.module.split(".", 1)[0] for item in record.imports),
        *(item.module.split(".", 1)[0] for item in record.namespaces),
        *record.forbidden_module_prefixes,
    }
    for label, paths in (
        ("untracked", _git_inventory(root, ignored=False)),
        ("ignored", _git_inventory(root, ignored=True)),
    ):
        for path in paths:
            if any(
                _matches_path(path, prefix)
                for prefix in record.forbidden_script_prefixes
            ):
                raise ExecutionClosureError(
                    f"{label} forbidden script shadow is present: {path}"
                )
            if any(_path_shadows_module(path, prefix) for prefix in protected):
                raise ExecutionClosureError(
                    f"{label} repository import shadow is present: {path}"
                )


def _require_import_layout(root: Path, record: ExecutionClosureRecord) -> None:
    for item in record.imports:
        conventional = _module_candidates(item.module)
        if "." in item.module and item.path not in conventional:
            raise ExecutionClosureError(
                "dotted import origin is not under src or repository root"
            )
        observed = set(path for path in conventional if _candidate_exists(root, path))
        observed.add(item.path)
        existing = {path for path in observed if _candidate_exists(root, path)}
        alternate_directories: set[str] = set()
        for path in _namespace_candidates(item.module):
            candidate = _nofollow_stat(root, path)
            if candidate is None:
                continue
            if stat.S_ISLNK(candidate.st_mode):
                raise ExecutionClosureError(f"import candidate is a symlink: {path}")
            if not stat.S_ISDIR(candidate.st_mode):
                raise ExecutionClosureError(
                    f"import namespace candidate is not a directory: {path}"
                )
            alternate_directories.add(path)
        if item.path.endswith("/__init__.py"):
            alternate_directories.discard(item.path.removesuffix("/__init__.py"))
        if existing != {item.path} or alternate_directories:
            raise ExecutionClosureError(
                f"import has duplicate or missing origins: {item.module}"
            )
    for item in record.namespaces:
        metadata = _nofollow_stat(root, item.path)
        if (
            metadata is None
            or stat.S_ISLNK(metadata.st_mode)
            or not stat.S_ISDIR(metadata.st_mode)
        ):
            raise ExecutionClosureError(
                f"namespace path is not a safe directory: {item.module}"
            )
        existing_directories: set[str] = set()
        for path in _namespace_candidates(item.module):
            candidate = _nofollow_stat(root, path)
            if candidate is None:
                continue
            if stat.S_ISLNK(candidate.st_mode):
                raise ExecutionClosureError(f"namespace candidate is a symlink: {path}")
            if not stat.S_ISDIR(candidate.st_mode):
                raise ExecutionClosureError(
                    f"namespace candidate is not a directory: {path}"
                )
            existing_directories.add(path)
        alternate_files = {
            path
            for path in _module_candidates(item.module)
            if _candidate_exists(root, path)
        }
        # An __init__.py in the declared directory is intentionally bypassed by
        # the explicit synthetic namespace, matching the HLT15 bootstrap
        # precedent.  A package/module at any other search root is a duplicate.
        declared_init = f"{item.path}/__init__.py"
        alternate_files.discard(declared_init)
        if existing_directories != {item.path} or alternate_files:
            raise ExecutionClosureError(
                f"namespace has duplicate origins: {item.module}"
            )
    _require_no_untracked_shadows(root, record)


def _preloaded_repository_modules(root: Path, record: ExecutionClosureRecord) -> None:
    protected = {
        *(item.module.split(".", 1)[0] for item in record.imports),
        *(item.module.split(".", 1)[0] for item in record.namespaces),
    }
    expected = {item.module for item in record.imports} | {
        item.module for item in record.namespaces
    }
    for name, module in tuple(sys.modules.items()):
        if not isinstance(name, str):
            raise ExecutionClosureError("sys.modules contains a non-text key")
        if module is None:
            continue
        if name in expected or any(
            _matches_module(name, prefix) for prefix in protected
        ):
            raise ExecutionClosureError(
                f"repository module was loaded before the guard: {name}"
            )
        origin = getattr(module, "__file__", None)
        if isinstance(origin, str) and _relative_origin(root, origin) is not None:
            raise ExecutionClosureError(
                f"repository file was executed before the guard: {name}"
            )
        search = getattr(module, "__path__", None)
        if search is not None:
            try:
                if any(
                    _relative_origin(root, str(path)) is not None for path in search
                ):
                    raise ExecutionClosureError(
                        f"repository namespace was loaded before the guard: {name}"
                    )
            except TypeError:
                raise ExecutionClosureError(
                    "loaded module has an invalid namespace path"
                )
        # A repository shadow of an already loaded stdlib/site module is still
        # ambiguous even though sys.modules would otherwise mask it.
        if name and all(part.isidentifier() for part in name.split(".")):
            candidates = _repository_candidate_paths(root, name)
            if candidates:
                raise ExecutionClosureError(
                    f"repository shadows an already loaded module: {name}"
                )


def _relative_origin(root: Path, origin: str) -> str | None:
    if origin in {"", "built-in", "frozen", "namespace", "proto19-verified-namespace"}:
        return None
    absolute = os.path.abspath(origin)
    real = os.path.realpath(absolute)
    contained: str | None = None
    for candidate in (absolute, real):
        try:
            if os.path.commonpath((str(root), candidate)) == str(root):
                contained = candidate
                break
        except ValueError:
            continue
    if contained is None:
        return None
    relative = os.path.relpath(contained, root).replace(os.sep, "/")
    try:
        return _relative(relative, "module origin")
    except ExecutionClosureError:
        return None


def build_execution_closure_record(
    root: Path,
    *,
    implementation_commit: str,
    files: Mapping[str, str],
    imports: Mapping[str, str],
    namespaces: Mapping[str, str] | None = None,
    forbidden_module_prefixes: Iterable[str] = (),
    forbidden_script_prefixes: Iterable[str] = (),
    authority_delta_paths: Iterable[str] = (),
) -> ExecutionClosureRecord:
    """Construct a prospective record at clean implementation commit A.

    The return value is content only.  It is suitable for later serialization
    into an authority artifact, but it does not authenticate that artifact or
    grant any execution permission.
    """
    _require_isolated_flags()
    _require_nofollow_support()
    repository = _repository_root(root)
    object_format = _git_text(repository, "rev-parse", "--show-object-format")
    if object_format not in {"sha1", "sha256"}:
        raise ExecutionClosureError("unsupported Git object format")
    implementation = _resolve_commit(
        repository, implementation_commit, object_format, "implementation commit"
    )
    if _git_text(repository, "rev-parse", "HEAD") != implementation:
        raise ExecutionClosureError(
            "prospective closure must be built at implementation HEAD"
        )
    _require_clean(repository)
    top = Path(os.path.realpath(_git_text(repository, "rev-parse", "--show-toplevel")))
    if top != repository:
        raise ExecutionClosureError("supplied path is not the Git repository root")

    pins: list[FilePin] = []
    for raw_path, kind in sorted(files.items()):
        path = _relative(raw_path, "closure file")
        if kind not in FILE_KINDS:
            raise ExecutionClosureError("closure file kind differs")
        oid, committed = _tree_blob(repository, implementation, path)
        live, _metadata = _nofollow_regular_bytes(repository, path)
        if live != committed:
            raise ExecutionClosureError(f"prospective live bytes differ: {path}")
        pins.append(FilePin(path, kind, oid, sha256(live).hexdigest()))
    pin_map = {item.path: item for item in pins}
    import_pins: list[ImportPin] = []
    for raw_module, raw_path in sorted(imports.items()):
        name = _module(raw_module, "import module")
        path = _relative(raw_path, "import origin")
        pin = pin_map.get(path)
        if pin is None:
            raise ExecutionClosureError("import origin is absent from file closure")
        import_pins.append(ImportPin(name, path, pin.sha256))
    namespace_pins = tuple(
        NamespacePin(
            _module(name, "namespace module"), _relative(path, "namespace path")
        )
        for name, path in sorted((namespaces or {}).items())
    )
    record = ExecutionClosureRecord(
        implementation_commit=implementation,
        git_object_format=object_format,
        files=tuple(pins),
        interpreter=current_interpreter_identity(),
        imports=tuple(import_pins),
        namespaces=namespace_pins,
        forbidden_module_prefixes=tuple(
            sorted(
                {
                    _module(item, "forbidden module prefix")
                    for item in forbidden_module_prefixes
                }
            )
        ),
        forbidden_script_prefixes=tuple(
            sorted(
                {
                    _prefix_path(item, "forbidden script prefix")
                    for item in forbidden_script_prefixes
                }
            )
        ),
        authority_delta_paths=tuple(
            sorted(
                {
                    _relative(item, "authority delta path")
                    for item in authority_delta_paths
                }
            )
        ),
    )
    _validate_record(record)
    _capture_files(repository, record, authority_commit=None)
    _require_import_layout(repository, record)
    return record


def verify_execution_closure(
    root: Path,
    record: ExecutionClosureRecord | Mapping[str, object],
    *,
    authority_commit: str,
    require_detached: bool = False,
) -> VerifiedExecutionClosure:
    """Verify clean authority commit C and capture every A-owned file byte."""
    _require_isolated_flags()
    _require_nofollow_support()
    parsed = (
        record
        if isinstance(record, ExecutionClosureRecord)
        else parse_execution_closure_record(record)
    )
    _validate_record(parsed)
    observed_identity = current_interpreter_identity()
    if observed_identity != parsed.interpreter:
        raise ExecutionClosureError("isolated interpreter/environment identity differs")
    if tuple(sys.path) != parsed.interpreter.base_sys_path:
        raise ExecutionClosureError(
            "repository paths were installed before closure verification"
        )
    repository = _repository_root(root)
    _require_repository(
        repository, parsed, authority_commit, require_detached=require_detached
    )
    captured = _capture_files(repository, parsed, authority_commit=authority_commit)
    _require_import_layout(repository, parsed)
    _preloaded_repository_modules(repository, parsed)
    return VerifiedExecutionClosure(
        root=repository,
        authority_commit=authority_commit,
        record=parsed,
        files=captured,
        import_path=(str(repository / "src"), str(repository)),
        require_detached=require_detached,
    )


class _VerifiedSourceLoader(importlib.abc.Loader):
    def __init__(self, guard: "ImportOriginGuard", pin: ImportPin) -> None:
        self.guard = guard
        self.pin = pin

    def create_module(self, spec: importlib.machinery.ModuleSpec) -> ModuleType | None:
        return None

    def exec_module(self, module: ModuleType) -> None:
        self.guard._require_active()
        self.guard._require_pin_live(self.pin)
        content = self.guard.verified.bytes_for_path(self.pin.path)
        filename = str(self.guard.verified.root / self.pin.path)
        code = compile(content, filename, "exec", dont_inherit=True)
        exec(code, module.__dict__)
        if getattr(module, "__file__", None) != filename:
            raise ExecutionClosureError(
                f"loaded module changed its origin: {self.pin.module}"
            )


class _SyntheticNamespaceLoader(importlib.abc.Loader):
    def __init__(self, guard: "ImportOriginGuard", pin: NamespacePin) -> None:
        self.guard = guard
        self.pin = pin

    def create_module(self, spec: importlib.machinery.ModuleSpec) -> ModuleType | None:
        return None

    def exec_module(self, module: ModuleType) -> None:
        self.guard._require_active()
        expected = [str(self.guard.verified.root / self.pin.path)]
        if list(getattr(module, "__path__", ())) != expected:
            raise ExecutionClosureError(
                f"synthetic namespace path differs: {self.pin.module}"
            )


class _ClosureFinder(importlib.abc.MetaPathFinder):
    def __init__(self, guard: "ImportOriginGuard") -> None:
        self.guard = guard

    def find_spec(
        self,
        fullname: str,
        path: Sequence[str] | None = None,
        target: ModuleType | None = None,
    ) -> importlib.machinery.ModuleSpec | None:
        del path, target
        self.guard._require_active()
        record = self.guard.verified.record
        if any(
            _matches_module(fullname, prefix)
            for prefix in record.forbidden_module_prefixes
        ):
            raise ExecutionClosureError(f"forbidden module import refused: {fullname}")
        namespace = self.guard.namespace_map.get(fullname)
        if namespace is not None:
            loader = _SyntheticNamespaceLoader(self.guard, namespace)
            spec = importlib.machinery.ModuleSpec(
                fullname,
                loader,
                origin="proto19-verified-namespace",
                is_package=True,
            )
            spec.submodule_search_locations = [
                str(self.guard.verified.root / namespace.path)
            ]
            spec.loader_state = self.guard.token
            return spec
        pin = self.guard.import_map.get(fullname)
        if pin is not None:
            loader = _VerifiedSourceLoader(self.guard, pin)
            package = pin.path.endswith("/__init__.py")
            origin = str(self.guard.verified.root / pin.path)
            spec = importlib.machinery.ModuleSpec(
                fullname, loader, origin=origin, is_package=package
            )
            spec.has_location = True
            spec.cached = None
            spec.loader_state = self.guard.token
            if package:
                spec.submodule_search_locations = [
                    str((self.guard.verified.root / pin.path).parent)
                ]
            return spec
        protected = {
            *(item.module.split(".", 1)[0] for item in record.imports),
            *(item.module.split(".", 1)[0] for item in record.namespaces),
        }
        if any(_matches_module(fullname, prefix) for prefix in protected):
            raise ExecutionClosureError(
                f"module is outside the permitted closure: {fullname}"
            )
        candidates = _repository_candidate_paths(self.guard.verified.root, fullname)
        if candidates:
            raise ExecutionClosureError(
                f"unlisted repository module import refused: {fullname}"
            )
        return None


class ImportOriginGuard:
    """Context-scoped importer for captured repository bytes.

    Keep the guard installed for the complete downstream operation.  Exiting
    the context removes modules loaded by this guard so a later unguarded
    repository import cannot silently reuse them.
    """

    def __init__(self, verified: VerifiedExecutionClosure) -> None:
        self.verified = verified
        self.import_map = {item.module: item for item in verified.record.imports}
        self.namespace_map = {item.module: item for item in verified.record.namespaces}
        self.token = object()
        self.finder = _ClosureFinder(self)
        self._active = False
        self._old_path: list[str] | None = None
        self._old_meta_path: list[object] | None = None
        self._old_path_hooks: list[object] | None = None

    def _require_active(self) -> None:
        if not self._active:
            raise ExecutionClosureError("import-origin guard is not active")

    def _require_pin_live(self, pin: ImportPin) -> None:
        live, _metadata = _nofollow_regular_bytes(self.verified.root, pin.path)
        captured = self.verified.bytes_for_path(pin.path)
        if live != captured or sha256(live).hexdigest() != pin.sha256:
            raise ExecutionClosureError(
                f"import source changed after verification: {pin.module}"
            )

    def install(self) -> "ImportOriginGuard":
        if self._active:
            raise ExecutionClosureError("import-origin guard is already active")
        if tuple(sys.path) != self.verified.record.interpreter.base_sys_path:
            raise ExecutionClosureError(
                "base sys.path changed before guard installation"
            )
        src_metadata = _nofollow_stat(self.verified.root, "src")
        if (
            src_metadata is None
            or not stat.S_ISDIR(src_metadata.st_mode)
            or stat.S_ISLNK(src_metadata.st_mode)
        ):
            raise ExecutionClosureError("src import root is not a safe directory")
        self._old_path = list(sys.path)
        self._old_meta_path = list(sys.meta_path)
        self._old_path_hooks = list(sys.path_hooks)
        sys.path[:] = [*self.verified.import_path, *self._old_path]
        sys.path_importer_cache.pop(self.verified.import_path[0], None)
        sys.path_importer_cache.pop(self.verified.import_path[1], None)
        sys.meta_path.insert(0, self.finder)
        self._active = True
        try:
            self.audit()
        except BaseException:
            self.close(audit=False)
            raise
        return self

    def close(self, *, audit: bool = True) -> None:
        if not self._active:
            return
        pending: BaseException | None = None
        if audit:
            try:
                self.audit()
            except BaseException as exc:  # Preserve cleanup on a failed final audit.
                pending = exc
        protected = {
            *(item.module.split(".", 1)[0] for item in self.verified.record.imports),
            *(item.module.split(".", 1)[0] for item in self.verified.record.namespaces),
        }
        for name, module in tuple(sys.modules.items()):
            if not isinstance(name, str):
                if pending is None:
                    pending = ExecutionClosureError(
                        "sys.modules contains a non-text key"
                    )
                sys.modules.pop(name, None)
                continue
            spec = getattr(module, "__spec__", None) if module is not None else None
            origin = getattr(module, "__file__", None) if module is not None else None
            search = getattr(module, "__path__", None) if module is not None else None
            repository_namespace = search is not None and any(
                _relative_origin(self.verified.root, str(path)) is not None
                for path in search
            )
            if (
                (spec is not None and getattr(spec, "loader_state", None) is self.token)
                or any(_matches_module(name, prefix) for prefix in protected)
                or (
                    isinstance(origin, str)
                    and _relative_origin(self.verified.root, origin) is not None
                )
                or repository_namespace
            ):
                sys.modules.pop(name, None)
        if self._old_meta_path is not None:
            sys.meta_path[:] = self._old_meta_path
        if self._old_path_hooks is not None:
            sys.path_hooks[:] = self._old_path_hooks
        if self._old_path is not None:
            sys.path[:] = self._old_path
        sys.path_importer_cache.pop(self.verified.import_path[0], None)
        sys.path_importer_cache.pop(self.verified.import_path[1], None)
        self._active = False
        if pending is not None:
            raise pending

    def __enter__(self) -> "ImportOriginGuard":
        return self.install()

    def __exit__(
        self, exc_type: object, exc: BaseException | None, traceback: object
    ) -> bool:
        self.close(audit=exc is None)
        return False

    def import_module(self, name: str) -> ModuleType:
        self._require_active()
        module = importlib.import_module(_module(name, "requested module"))
        self.audit()
        return module

    def assert_script_allowed(self, script: Path) -> str:
        """Validate a prospective script path; this function does not execute it."""
        self._require_active()
        absolute = os.path.abspath(os.fspath(script))
        relative = _relative_origin(self.verified.root, absolute)
        if relative is None:
            raise ExecutionClosureError("script is outside the repository")
        if any(
            _matches_path(relative, prefix)
            for prefix in self.verified.record.forbidden_script_prefixes
        ):
            raise ExecutionClosureError(f"forbidden script refused: {relative}")
        pin = next(
            (item.pin for item in self.verified.files if item.pin.path == relative),
            None,
        )
        if pin is None or pin.kind != "python":
            raise ExecutionClosureError(
                "script is outside the permitted Python closure"
            )
        live, _metadata = _nofollow_regular_bytes(self.verified.root, relative)
        if live != self.verified.bytes_for_path(relative):
            raise ExecutionClosureError("script bytes changed after verification")
        return relative

    def audit(self) -> None:
        """Recheck Git/live bytes and every currently observed local module."""
        self._require_active()
        expected_path = (
            *self.verified.import_path,
            *self.verified.record.interpreter.base_sys_path,
        )
        expected_meta_path = [self.finder, *(self._old_meta_path or ())]
        if (
            tuple(sys.path) != expected_path
            or sys.meta_path != expected_meta_path
            or sys.path_hooks != (self._old_path_hooks or [])
        ):
            raise ExecutionClosureError("guarded import configuration drifted")
        _require_repository(
            self.verified.root,
            self.verified.record,
            self.verified.authority_commit,
            require_detached=self.verified.require_detached,
        )
        _capture_files(
            self.verified.root,
            self.verified.record,
            authority_commit=self.verified.authority_commit,
        )
        _require_import_layout(self.verified.root, self.verified.record)
        for name, module in tuple(sys.modules.items()):
            if not isinstance(name, str):
                raise ExecutionClosureError("sys.modules contains a non-text key")
            if module is None:
                continue
            pin = self.import_map.get(name)
            namespace = self.namespace_map.get(name)
            spec = getattr(module, "__spec__", None)
            if pin is not None:
                expected = str(self.verified.root / pin.path)
                if (
                    not isinstance(
                        getattr(module, "__loader__", None), _VerifiedSourceLoader
                    )
                    or module.__loader__.guard is not self  # type: ignore[union-attr]
                    or getattr(module, "__file__", None) != expected
                    or spec is None
                    or getattr(spec, "loader_state", None) is not self.token
                ):
                    raise ExecutionClosureError(
                        f"observed import origin/loader differs: {name}"
                    )
                self._require_pin_live(pin)
                continue
            if namespace is not None:
                expected = [str(self.verified.root / namespace.path)]
                if (
                    not isinstance(
                        getattr(module, "__loader__", None), _SyntheticNamespaceLoader
                    )
                    or module.__loader__.guard is not self  # type: ignore[union-attr]
                    or list(getattr(module, "__path__", ())) != expected
                    or spec is None
                    or getattr(spec, "loader_state", None) is not self.token
                ):
                    raise ExecutionClosureError(
                        f"observed namespace origin/loader differs: {name}"
                    )
                continue
            protected = {
                *(
                    item.module.split(".", 1)[0]
                    for item in self.verified.record.imports
                ),
                *(
                    item.module.split(".", 1)[0]
                    for item in self.verified.record.namespaces
                ),
            }
            if any(_matches_module(name, prefix) for prefix in protected):
                raise ExecutionClosureError(
                    f"observed module is outside permitted closure: {name}"
                )
            origin = getattr(module, "__file__", None)
            if (
                isinstance(origin, str)
                and _relative_origin(self.verified.root, origin) is not None
            ):
                raise ExecutionClosureError(
                    f"observed repository module is outside closure: {name}"
                )
            search = getattr(module, "__path__", None)
            if search is not None and any(
                _relative_origin(self.verified.root, str(path)) is not None
                for path in search
            ):
                raise ExecutionClosureError(
                    f"observed repository namespace is outside closure: {name}"
                )


@contextmanager
def isolated_import_bootstrap(
    root: Path,
    record: ExecutionClosureRecord | Mapping[str, object],
    *,
    authority_commit: str,
    require_detached: bool = False,
) -> Iterator[ImportOriginGuard]:
    """Verify the image and install the guard for one bounded operation."""
    verified = verify_execution_closure(
        root,
        record,
        authority_commit=authority_commit,
        require_detached=require_detached,
    )
    with ImportOriginGuard(verified) as guard:
        yield guard


__all__ = [
    "CapturedFile",
    "ExecutionClosureError",
    "ExecutionClosureRecord",
    "FILE_KINDS",
    "FilePin",
    "ImportOriginGuard",
    "ImportPin",
    "InterpreterIdentity",
    "NamespacePin",
    "SCHEMA",
    "VerifiedExecutionClosure",
    "build_execution_closure_record",
    "current_interpreter_identity",
    "isolated_import_bootstrap",
    "parse_execution_closure_record",
    "verify_execution_closure",
]
