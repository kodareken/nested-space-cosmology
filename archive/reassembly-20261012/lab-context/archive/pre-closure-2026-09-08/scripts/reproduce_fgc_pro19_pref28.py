#!/usr/bin/env python3
"""Build or verify the outcome-neutral PRO19 generation-nine terminal binder."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import stat
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto19_pref28_binder import (  # noqa: E402
    ARTIFACT_ID,
    CONFIG_PATH,
    RESULT_PATH,
    build_pref28_result,
    canonical_result,
    validate_compact_result,
)


_MAX_COMPACT_LEAF_BYTES = 8 * 1024 * 1024


def _relative_parts(relative: str, label: str) -> tuple[str, ...]:
    path = Path(relative)
    if (
        path.is_absolute()
        or not path.parts
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ValueError(f"{label} has an unsafe repository-relative path")
    return tuple(path.parts)


def _directory_flags() -> int:
    return os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)


def _open_parent(repository: Path, relative: str, label: str) -> tuple[int, str]:
    """Open a repository-relative parent chain without following any link."""

    repository = Path(repository)
    parts = _relative_parts(relative, label)
    root_before = repository.lstat()
    if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
        raise ValueError(f"{label} repository root is unsafe")
    current = os.open(repository, _directory_flags())
    try:
        root_opened = os.fstat(current)
        if not stat.S_ISDIR(root_opened.st_mode) or (
            root_before.st_dev,
            root_before.st_ino,
        ) != (root_opened.st_dev, root_opened.st_ino):
            raise ValueError(f"{label} repository root changed before access")
        for component in parts[:-1]:
            before = os.stat(component, dir_fd=current, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise ValueError(f"{label} ancestor is unsafe: {component}")
            child = os.open(component, _directory_flags(), dir_fd=current)
            opened = os.fstat(child)
            if not stat.S_ISDIR(opened.st_mode) or (before.st_dev, before.st_ino) != (
                opened.st_dev,
                opened.st_ino,
            ):
                os.close(child)
                raise ValueError(f"{label} ancestor changed: {component}")
            os.close(current)
            current = child
        return current, parts[-1]
    except BaseException:
        os.close(current)
        raise


def _identity(metadata: os.stat_result) -> tuple[int, ...]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_mode,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
        metadata.st_nlink,
    )


def _read_unique_regular(repository: Path, relative: str, label: str) -> bytes:
    """Read one bounded leaf through descriptor-relative no-follow traversal."""

    parent, leaf = _open_parent(repository, relative, label)
    descriptor = -1
    try:
        before = os.stat(leaf, dir_fd=parent, follow_symlinks=False)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_size < 0
            or before.st_size > _MAX_COMPACT_LEAF_BYTES
        ):
            raise ValueError(f"{label} is not a unique bounded regular file")
        descriptor = os.open(
            leaf,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent,
        )
        opened = os.fstat(descriptor)
        identity = _identity(before)
        if identity != _identity(opened):
            raise ValueError(f"{label} changed before read")
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1024 * 1024, remaining))
            if not block:
                raise ValueError(f"{label} ended early")
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            raise ValueError(f"{label} grew during read")
        if identity != _identity(os.fstat(descriptor)):
            raise ValueError(f"{label} changed during read")
        if identity != _identity(os.stat(leaf, dir_fd=parent, follow_symlinks=False)):
            raise ValueError(f"{label} path changed during read")
        return b"".join(chunks)
    finally:
        if descriptor != -1:
            os.close(descriptor)
        os.close(parent)


def live_result_bytes(repository: Path) -> bytes:
    config_raw = _read_unique_regular(repository, CONFIG_PATH, "PREF28 configuration")
    return canonical_result(build_pref28_result(config_raw, repository))


def write_result(repository: Path) -> None:
    destination = repository / RESULT_PATH
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    if temporary.exists():
        raise SystemExit("stale PREF28 temporary result exists")
    with temporary.open("xb") as handle:
        handle.write(live_result_bytes(repository))
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, destination)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--verify-compact", action="store_true")
    args = parser.parse_args(argv)

    config_raw = _read_unique_regular(ROOT, CONFIG_PATH, "PREF28 configuration")
    if args.write:
        write_result(ROOT)
        live_verified = True
    elif args.verify:
        result_raw = _read_unique_regular(ROOT, RESULT_PATH, "PREF28 result")
        if result_raw != live_result_bytes(ROOT):
            raise SystemExit("PRO19 PREF28 result differs from live terminal boundary")
        live_verified = True
    else:
        result_raw = _read_unique_regular(ROOT, RESULT_PATH, "PREF28 result")
        validate_compact_result(config_raw, result_raw)
        live_verified = False

    print(
        json.dumps(
            {
                "artifact_id": ARTIFACT_ID,
                "accepted_state_advanced": False,
                "candidate_branch_opened": False,
                "common_event_completed": False,
                "live_boundary_verified": live_verified,
                "physical_result_earned": False,
                "root_cause_localized": False,
                "store_mutated": False,
                "verified": bool(args.verify or args.verify_compact),
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
