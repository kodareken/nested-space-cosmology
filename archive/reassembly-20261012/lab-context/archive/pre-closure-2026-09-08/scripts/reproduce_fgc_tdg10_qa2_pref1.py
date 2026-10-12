#!/usr/bin/env python3
"""Construct once or compactly verify FGC-1-TDG10-QA2-PREF1.

Default mode validates an existing compact result and never reconstructs a
shadow or reads the raw namespace/store.  If the result is absent, default
mode fails clearly.  Explicit ``--live`` is the sole path that imports the
independent binder's live entry point; it binds an already-existing QA2 raw
terminal and writes the compact result with atomic no-clobber publication.
It never launches QA2.
"""

from __future__ import annotations

import argparse
import importlib
import os
from pathlib import Path
import secrets
import stat
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

CONFIG_PATH = "configs/fgc/fgc-1-tdg10-qa2-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg10-qa2-pref1.json"
BINDER_MODULE = "recursive_horizons.fgc.evolution.tdg10_qa2_pref1_binder"
ABSENT_RESULT_MESSAGE = (
    "FGC-1-TDG10-QA2-PREF1 compact result is absent at "
    f"{RESULT_PATH}. Default mode validates a compact result only when that "
    "file exists; it never runs live binding. Pass --live only for the "
    "explicit one-time independent binder of the existing QA2 terminal."
)
MISSING_BINDER_MESSAGE = (
    "FGC-1-TDG10-QA2-PREF1 independent binder module is not present "
    f"({BINDER_MODULE}). This reproducer never invents a compact result."
)


def load_binder():
    try:
        return importlib.import_module(BINDER_MODULE)
    except ImportError as exc:
        raise SystemExit(MISSING_BINDER_MESSAGE) from exc


def _identity(metadata: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
        metadata.st_nlink,
    )


def _open_parent(root: Path, relative: str, *, create: bool) -> tuple[int, str]:
    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        raise SystemExit(f"unsafe PREF1 repository-relative path: {relative}")
    root_before = root.lstat()
    if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
        raise SystemExit("unsafe PREF1 repository root")
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        if _identity(os.fstat(descriptor)) != _identity(root_before):
            raise SystemExit("PREF1 repository root changed during access")
        for part in candidate.parent.parts:
            try:
                before = os.stat(part, dir_fd=descriptor, follow_symlinks=False)
            except FileNotFoundError:
                if not create:
                    raise
                os.mkdir(part, mode=0o755, dir_fd=descriptor)
                before = os.stat(part, dir_fd=descriptor, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise SystemExit(f"unsafe PREF1 parent directory: {part}")
            child = os.open(
                part,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=descriptor,
            )
            active = os.fstat(child)
            after = os.stat(part, dir_fd=descriptor, follow_symlinks=False)
            if _identity(before) != _identity(active) or _identity(before) != _identity(
                after
            ):
                os.close(child)
                raise SystemExit(f"PREF1 parent directory changed: {part}")
            os.close(descriptor)
            descriptor = child
        return descriptor, candidate.name
    except Exception:
        os.close(descriptor)
        raise


def read_repository_leaf(
    root: Path, relative: str, *, missing_ok: bool = False
) -> bytes | None:
    try:
        parent_fd, name = _open_parent(root, relative, create=False)
    except FileNotFoundError:
        if missing_ok:
            return None
        raise SystemExit(f"required PREF1 leaf is absent: {relative}") from None
    descriptor = -1
    try:
        try:
            before = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            if missing_ok:
                return None
            raise SystemExit(f"required PREF1 leaf is absent: {relative}") from None
        if (
            stat.S_ISLNK(before.st_mode)
            or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_size < 0
            or before.st_size > 512 * 1024 * 1024
        ):
            raise SystemExit(f"unsafe PREF1 leaf: {relative}")
        descriptor = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        active = os.fstat(descriptor)
        if _identity(active) != _identity(before):
            raise SystemExit(f"PREF1 leaf changed before read: {relative}")
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                raise SystemExit(f"short PREF1 leaf read: {relative}")
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            raise SystemExit(f"PREF1 leaf grew during read: {relative}")
        after = os.fstat(descriptor)
        final = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if _identity(before) != _identity(after) or _identity(before) != _identity(
            final
        ):
            raise SystemExit(f"PREF1 leaf changed during read: {relative}")
        return b"".join(chunks)
    finally:
        if descriptor != -1:
            os.close(descriptor)
        os.close(parent_fd)


def write_canonical_result(root: Path, relative: str, payload: bytes) -> None:
    """Durably publish one result without replacing any existing path."""

    parent_fd, destination_name = _open_parent(root, relative, create=True)
    try:
        try:
            os.stat(destination_name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise SystemExit(f"refusing to overwrite existing PREF1 result: {relative}")
    except Exception:
        os.close(parent_fd)
        raise
    temporary_name = f".{destination_name}.tmp-{os.getpid()}-{secrets.token_hex(8)}"
    descriptor = -1
    linked = False
    try:
        descriptor = os.open(
            temporary_name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o644,
            dir_fd=parent_fd,
        )
        pending = memoryview(payload)
        while pending:
            written = os.write(descriptor, pending)
            if written <= 0:
                raise OSError("short PREF1 compact-result write")
            pending = pending[written:]
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = -1
        try:
            os.link(
                temporary_name,
                destination_name,
                src_dir_fd=parent_fd,
                dst_dir_fd=parent_fd,
                follow_symlinks=False,
            )
        except FileExistsError as exc:
            raise SystemExit(
                f"refusing to overwrite existing PREF1 result: {relative}"
            ) from exc
        linked = True
        temporary_stat = os.stat(
            temporary_name, dir_fd=parent_fd, follow_symlinks=False
        )
        destination_stat = os.stat(
            destination_name, dir_fd=parent_fd, follow_symlinks=False
        )
        if (
            not stat.S_ISREG(destination_stat.st_mode)
            or destination_stat.st_nlink != 2
            or _identity(temporary_stat) != _identity(destination_stat)
        ):
            raise RuntimeError("PREF1 destination identity differs after link")
        destination_fd = os.open(
            destination_name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        try:
            if _identity(os.fstat(destination_fd)) != _identity(destination_stat):
                raise RuntimeError("PREF1 destination changed before fsync")
            os.fsync(destination_fd)
            if _identity(
                os.stat(destination_name, dir_fd=parent_fd, follow_symlinks=False)
            ) != _identity(destination_stat):
                raise RuntimeError("PREF1 destination changed before temp unlink")
        finally:
            os.close(destination_fd)
        os.unlink(temporary_name, dir_fd=parent_fd)
        os.fsync(parent_fd)
    except Exception:
        if descriptor != -1:
            os.close(descriptor)
        try:
            os.unlink(temporary_name, dir_fd=parent_fd)
        except FileNotFoundError:
            pass
        if linked:
            try:
                os.stat(destination_name, dir_fd=parent_fd, follow_symlinks=False)
            except FileNotFoundError as exc:
                raise RuntimeError(
                    "PREF1 result link disappeared during publication"
                ) from exc
        raise
    finally:
        os.close(parent_fd)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--output", default=RESULT_PATH)
    arguments = parser.parse_args(argv)
    if arguments.live:
        binder = load_binder()
        config = read_repository_leaf(ROOT, binder.CONFIG_PATH)
        assert config is not None
        result = binder.build_pref1_result(config, ROOT, live=True)
        write_canonical_result(ROOT, arguments.output, binder.canonical_result(result))
        print(
            "FGC-1-TDG10-QA2-PREF1 live binder wrote the compact result; "
            "the tested SSPRK3-on-SBP4 exact-C remedy is nonrobust on retries "
            "4/5, and no state, production method, calibration, or physics "
            "claim was opened"
        )
        return 0
    result_raw = read_repository_leaf(ROOT, arguments.output, missing_ok=True)
    if result_raw is None:
        raise SystemExit(ABSENT_RESULT_MESSAGE)
    binder = load_binder()
    config = read_repository_leaf(ROOT, binder.CONFIG_PATH)
    assert config is not None
    binder.validate_compact_result(config, result_raw)
    print(
        "FGC-1-TDG10-QA2-PREF1 compact result verified; this rejects only "
        "the tested remedy on the retry-4/5 neighborhood and opens no state, "
        "candidate, mechanism, or physics claim"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
