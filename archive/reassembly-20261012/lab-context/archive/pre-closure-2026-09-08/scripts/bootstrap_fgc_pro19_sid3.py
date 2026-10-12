#!/usr/bin/env python3
"""Isolated committed-byte bootstrap for SID3 preflight or resume.

Intended invocation (the script body itself comes from implementation commit
``A`` rather than the live worktree)::

    git show A:scripts/bootstrap_fgc_pro19_sid3.py \
      | python -I -B - --root REPO --implementation-commit A \
          --authority-commit C --preflight

Everything imported before the execution-origin guard is standard-library
code.  This bootstrap performs no store mutation itself.  In ``--resume``
mode it delegates mutation to the guarded runner, which must repeat the
read-only SID3/store preflight immediately before acquiring its writer lease.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
import sys
from types import ModuleType
from typing import Any, Iterable, Mapping, Sequence


BOOTSTRAP_PATH = "scripts/bootstrap_fgc_pro19_sid3.py"
CLOSURE_MODULE_PATH = (
    "src/recursive_horizons/fgc/evolution/proto19_execution_closure.py"
)
CONFIG_PATH = "configs/fgc/fgc-1-pro19-sid3-auth1.toml"
CLOSURE_PATH = "configs/fgc/fgc-1-pro19-sid3-execution-closure.json"
RESULT_PATH = "results/fgc-1-pro19-sid3-auth1.json"
AUTHORITY_MODULE = "recursive_horizons.fgc.evolution.proto19_resume_authority"
RUNNER_MODULE = "scripts.run_fgc_pro19_event1"
DEFAULT_STORE_PATH = "runs/fgc-2-sf1/proto17/calibration"
PROCESS_MARKER = "FGC-PRO19-SID3-EVENT1"


class SID3BootstrapError(RuntimeError):
    """The isolated committed image cannot be authenticated or dispatched."""


def _canonical(value: object) -> bytes:
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
    except (TypeError, ValueError) as exc:
        raise SID3BootstrapError("bootstrap output is not canonicalizable") from exc


def _reject_duplicates(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in items:
        if key in value:
            raise ValueError(key)
        value[key] = item
    return value


def _json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=_reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise SID3BootstrapError(f"{label} is malformed") from exc
    if not isinstance(value, dict) or _canonical(value) != raw:
        raise SID3BootstrapError(f"{label} is noncanonical")
    return value


def _relative(value: str, label: str) -> tuple[str, ...]:
    if not value or "\\" in value or "\x00" in value:
        raise SID3BootstrapError(f"{label} is not a repository path")
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or path.as_posix() != value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise SID3BootstrapError(f"{label} escapes or aliases the repository")
    return path.parts


def _repository_root(value: Path) -> Path:
    root = Path(os.path.abspath(os.fspath(value)))
    try:
        metadata = root.lstat()
    except OSError as exc:
        raise SID3BootstrapError("repository root is absent") from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise SID3BootstrapError("repository root is unsafe")
    return root


def _read_nofollow(root: Path, relative: str, label: str) -> bytes:
    parts = _relative(relative, label)
    directory_flags = (
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    )
    leaf_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    directories: list[int] = []
    leaf = -1
    try:
        parent = os.open(root, directory_flags)
        directories.append(parent)
        for component in parts[:-1]:
            before = os.stat(component, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise SID3BootstrapError(f"{label} ancestor is unsafe")
            child = os.open(component, directory_flags, dir_fd=parent)
            active = os.fstat(child)
            if (before.st_dev, before.st_ino) != (active.st_dev, active.st_ino):
                os.close(child)
                raise SID3BootstrapError(f"{label} ancestor raced")
            directories.append(child)
            parent = child
        before = os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
        if (
            stat.S_ISLNK(before.st_mode)
            or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
        ):
            raise SID3BootstrapError(f"{label} is not a single-link regular file")
        leaf = os.open(parts[-1], leaf_flags, dir_fd=parent)
        active = os.fstat(leaf)
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
            raise SID3BootstrapError(f"{label} raced before read")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(leaf, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        after = os.fstat(leaf)
        if identity != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
            after.st_nlink,
        ):
            raise SID3BootstrapError(f"{label} changed during read")
        return b"".join(chunks)
    except SID3BootstrapError:
        raise
    except OSError as exc:
        raise SID3BootstrapError(f"{label} cannot be read safely") from exc
    finally:
        if leaf != -1:
            os.close(leaf)
        for descriptor in reversed(directories):
            os.close(descriptor)


def _git(root: Path, arguments: Sequence[str], label: str) -> bytes:
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment["GIT_OPTIONAL_LOCKS"] = "0"
    environment["LC_ALL"] = "C"
    process = subprocess.run(
        (
            "git",
            "--no-optional-locks",
            "-c",
            "core.fsmonitor=false",
            "-c",
            "core.untrackedCache=false",
            "-C",
            str(root),
            *arguments,
        ),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=environment,
    )
    if process.returncode:
        raise SID3BootstrapError(f"{label} cannot be resolved")
    return process.stdout


def _resolve_commit(root: Path, value: str, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise SID3BootstrapError(f"{label} is absent")
    raw = _git(root, ("rev-parse", "--verify", f"{value}^{{commit}}"), label)
    try:
        commit = raw.decode("ascii").strip()
    except UnicodeDecodeError as exc:
        raise SID3BootstrapError(f"{label} is malformed") from exc
    if len(commit) not in {40, 64} or any(c not in "0123456789abcdef" for c in commit):
        raise SID3BootstrapError(f"{label} is malformed")
    if value != commit:
        raise SID3BootstrapError(f"{label} must be the full commit id")
    return commit


def _git_show(root: Path, commit: str, relative: str, label: str) -> bytes:
    _relative(relative, label)
    return _git(root, ("show", f"{commit}:{relative}"), label)


def _authority_bytes(
    root: Path, authority_commit: str, relative: str, label: str
) -> bytes:
    live = _read_nofollow(root, relative, label)
    committed = _git_show(root, authority_commit, relative, label)
    if live != committed:
        raise SID3BootstrapError(f"{label} differs from authority commit")
    return live


def _load_closure_module(
    root: Path, implementation_commit: str
) -> tuple[ModuleType, bytes]:
    source = _git_show(
        root,
        implementation_commit,
        CLOSURE_MODULE_PATH,
        "execution-closure implementation",
    )
    name = f"_fgc_pro19_execution_closure_{sha256(source).hexdigest()}"
    origin = f"/__fgc_committed__/{implementation_commit}/{CLOSURE_MODULE_PATH}"
    module = ModuleType(name)
    module.__file__ = origin
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = importlib.util.spec_from_loader(name, loader=None, origin=origin)
    sys.modules[name] = module
    try:
        code = compile(source, origin, "exec", dont_inherit=True)
        exec(code, module.__dict__)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    required = {
        "parse_execution_closure_record",
        "verify_execution_closure",
        "ImportOriginGuard",
    }
    if any(not hasattr(module, item) for item in required):
        sys.modules.pop(name, None)
        raise SID3BootstrapError("execution-closure implementation API differs")
    return module, source


def _require_isolated_interpreter() -> None:
    if not sys.flags.isolated or not sys.dont_write_bytecode:
        raise SID3BootstrapError("bootstrap requires python -I -B")
    if not sys.argv or sys.argv[0] != "-":
        raise SID3BootstrapError(
            "bootstrap must be executed from committed stdin bytes"
        )


def _mark_committed_stdin_origin(implementation_commit: str) -> None:
    """Keep the ``<stdin>`` pseudo-origin outside repository-origin audits."""

    main_module = sys.modules.get("__main__")
    if main_module is None:
        raise SID3BootstrapError("bootstrap __main__ module is absent")
    origin = getattr(main_module, "__file__", None)
    if origin == "<stdin>":
        main_module.__file__ = (
            f"/__fgc_committed__/{implementation_commit}/{BOOTSTRAP_PATH}"
        )


def execute(
    *,
    root: Path,
    implementation_commit: str,
    authority_commit: str,
    process_marker: str,
    mode: str,
    store_root: Path | None = None,
) -> Mapping[str, Any]:
    """Authenticate then dispatch exactly one guarded operation."""

    _require_isolated_interpreter()
    if process_marker != PROCESS_MARKER:
        raise SID3BootstrapError("SID3 process marker differs")
    repository = _repository_root(root)
    implementation = _resolve_commit(
        repository, implementation_commit, "implementation commit"
    )
    _mark_committed_stdin_origin(implementation)
    authority = _resolve_commit(repository, authority_commit, "authority commit")
    closure_module, bootstrap_closure_source = _load_closure_module(
        repository, implementation
    )
    config_raw = _authority_bytes(repository, authority, CONFIG_PATH, "SID3 config")
    closure_raw = _authority_bytes(
        repository, authority, CLOSURE_PATH, "SID3 execution closure"
    )
    result_raw = _authority_bytes(repository, authority, RESULT_PATH, "SID3 result")
    closure_mapping = _json(closure_raw, "SID3 execution closure")
    try:
        record = closure_module.parse_execution_closure_record(closure_mapping)
    except Exception as exc:
        raise SID3BootstrapError("SID3 execution closure cannot be parsed") from exc
    if record.implementation_commit != implementation:
        raise SID3BootstrapError("bootstrap implementation commit differs")
    pins = {item.path: item for item in record.files}
    for path, expected_source in (
        (CLOSURE_MODULE_PATH, bootstrap_closure_source),
        (
            BOOTSTRAP_PATH,
            _git_show(repository, implementation, BOOTSTRAP_PATH, "bootstrap image"),
        ),
    ):
        pin = pins.get(path)
        if (
            pin is None
            or pin.kind != "python"
            or pin.sha256 != sha256(expected_source).hexdigest()
        ):
            raise SID3BootstrapError(f"committed bootstrap pin differs: {path}")
    try:
        verified = closure_module.verify_execution_closure(
            repository,
            record,
            authority_commit=authority,
            require_detached=False,
        )
    except Exception as exc:
        raise SID3BootstrapError("live committed execution closure differs") from exc
    captured = {item.pin.path: item.content for item in verified.files}
    guard = closure_module.ImportOriginGuard(verified)
    with guard:
        try:
            authority_api = guard.import_module(AUTHORITY_MODULE)
            runner = guard.import_module(RUNNER_MODULE)
            receipt = authority_api.authorize_resume_image(
                config_raw,
                closure_raw,
                result_raw,
                captured_files=captured,
                authority_commit=authority,
            )
        except Exception as exc:
            raise SID3BootstrapError("guarded SID3 authority differs") from exc
        if (
            receipt.authority_commit != authority
            or receipt.implementation_commit != implementation
            or receipt.store_path != DEFAULT_STORE_PATH
        ):
            raise SID3BootstrapError("guarded SID3 receipt identity differs")
        expected_store = repository / receipt.store_path
        selected_store = (
            expected_store
            if store_root is None
            else Path(os.path.abspath(os.fspath(store_root)))
        )
        if selected_store != expected_store:
            raise SID3BootstrapError("store path differs from SID3 authority")
        guard.audit()
        if mode == "preflight":
            operation = getattr(runner, "sid3_resume_preflight", None)
        elif mode == "resume":
            operation = getattr(runner, "sid3_resume", None)
        else:
            raise SID3BootstrapError("unknown SID3 bootstrap mode")
        if not callable(operation):
            raise SID3BootstrapError("guarded runner SID3 API is absent")
        try:
            outcome = operation(
                repository,
                authority=receipt,
                store_root=selected_store,
            )
        except Exception as exc:
            raise SID3BootstrapError(f"guarded SID3 {mode} failed") from exc
        guard.audit()
        if not isinstance(outcome, Mapping):
            raise SID3BootstrapError("guarded SID3 outcome is not a mapping")
        # Ensure the printed representation is finite/canonical before leaving
        # the guard; this does not create any tracked or run-store output.
        _canonical(dict(outcome))
        return dict(outcome)


def _arguments(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--implementation-commit", required=True)
    parser.add_argument("--authority-commit", required=True)
    parser.add_argument("--process-marker", required=True)
    parser.add_argument("--store-root", type=Path)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--preflight", action="store_true")
    modes.add_argument("--resume", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _arguments(argv)
    mode = "preflight" if arguments.preflight else "resume"
    try:
        outcome = execute(
            root=arguments.root,
            implementation_commit=arguments.implementation_commit,
            authority_commit=arguments.authority_commit,
            process_marker=arguments.process_marker,
            mode=mode,
            store_root=arguments.store_root,
        )
    except SID3BootstrapError as exc:
        print(f"SID3 bootstrap stopped: {exc}", file=sys.stderr)
        return 2
    sys.stdout.buffer.write(_canonical(dict(outcome)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
