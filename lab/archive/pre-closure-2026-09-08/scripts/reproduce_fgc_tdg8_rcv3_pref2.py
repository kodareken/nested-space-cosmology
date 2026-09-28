#!/usr/bin/env python3
"""Build or verify the outcome-neutral FGC-1-TDG8-RCV3-PREF2 binder."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import stat
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg8_rcv3_pref2_binder as binder  # noqa: E402


_MAX_COMPACT_BYTES = 8 * 1024 * 1024


def _parts(relative: str, label: str) -> tuple[str, ...]:
    path = Path(relative)
    if path.is_absolute() or path.as_posix() != relative or any(
        part in {"", ".", ".."} for part in path.parts
    ):
        raise ValueError(f"{label} path is unsafe")
    return path.parts


def _flags() -> int:
    return os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)


def _read(relative: str, label: str) -> bytes:
    parts = _parts(relative, label)
    parent = os.open(ROOT, _flags())
    leaf = -1
    try:
        for component in parts[:-1]:
            before = os.stat(component, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise ValueError(f"{label} ancestor is unsafe")
            child = os.open(component, _flags(), dir_fd=parent)
            active = os.fstat(child)
            if (before.st_dev, before.st_ino) != (active.st_dev, active.st_ino):
                os.close(child)
                raise ValueError(f"{label} ancestor raced")
            os.close(parent)
            parent = child
        before = os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
        if (
            stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1 or before.st_size < 0
            or before.st_size > _MAX_COMPACT_BYTES
        ):
            raise ValueError(f"{label} is not a bounded single-linked file")
        leaf = os.open(
            parts[-1], os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0), dir_fd=parent,
        )
        opened = os.fstat(leaf)
        identity = (
            before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
            before.st_ctime_ns, before.st_nlink,
        )
        if identity != (
            opened.st_dev, opened.st_ino, opened.st_size, opened.st_mtime_ns,
            opened.st_ctime_ns, opened.st_nlink,
        ):
            raise ValueError(f"{label} changed before read")
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(leaf, min(1 << 20, remaining))
            if not block:
                raise ValueError(f"{label} ended early")
            chunks.append(block)
            remaining -= len(block)
        if os.read(leaf, 1):
            raise ValueError(f"{label} grew during read")
        after = os.fstat(leaf)
        path_after = os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
        if identity != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
            after.st_ctime_ns, after.st_nlink,
        ) or (path_after.st_dev, path_after.st_ino) != (before.st_dev, before.st_ino):
            raise ValueError(f"{label} changed during read")
        return b"".join(chunks)
    finally:
        if leaf != -1:
            os.close(leaf)
        os.close(parent)


def _live() -> bytes:
    config = _read(binder.CONFIG_PATH, "PREF2 config")
    return binder.canonical_result(binder.build_pref2_result(config, ROOT))


def _publish_new(raw: bytes) -> None:
    destination = ROOT / binder.RESULT_PATH
    directory = os.open(destination.parent, _flags())
    leaf = -1
    created = False
    try:
        leaf = os.open(
            destination.name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o644,
            dir_fd=directory,
        )
        created = True
        view = memoryview(raw)
        while view:
            written = os.write(leaf, view)
            if written <= 0:
                raise OSError("short PREF2 result write")
            view = view[written:]
        os.fsync(leaf)
        os.close(leaf)
        leaf = -1
        os.fsync(directory)
    except Exception:
        if leaf != -1:
            os.close(leaf)
        if created:
            try:
                os.unlink(destination.name, dir_fd=directory)
            except OSError:
                pass
        raise
    finally:
        os.close(directory)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--verify", action="store_true")
    modes.add_argument("--verify-compact", action="store_true")
    args = parser.parse_args(argv)

    config = _read(binder.CONFIG_PATH, "PREF2 config")
    if args.write:
        raw = _live()
        _publish_new(raw)
        if _read(binder.RESULT_PATH, "PREF2 result") != raw:
            raise SystemExit("published PREF2 result bytes differ")
        live = True
    elif args.verify:
        raw = _read(binder.RESULT_PATH, "PREF2 result")
        if raw != _live():
            raise SystemExit("PREF2 result differs from the live terminal store")
        live = True
    else:
        raw = _read(binder.RESULT_PATH, "PREF2 result")
        binder.validate_compact_result(config, raw)
        live = False
    print(
        f"{binder.ARTIFACT_ID}: pass; "
        f"terminal_store_read={str(live).lower()}; candidate_or_physics=false"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
