#!/usr/bin/env python3
"""Write once or compactly verify the TDG8 successor prelaunch authority.

``--write`` performs the one allowed destination-absence observation before
publishing the compact result.  ``--verify-compact`` reads only the tracked
config/result pair; it never reopens either the terminal source campaign or a
now-created successor destination.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import stat
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.tdg8_successor_authority import (  # noqa: E402
    ARTIFACT_ID,
    CONFIG_PATH,
    RESULT_PATH,
    _read_repository_leaf,
    build_prelaunch,
    canonical_result,
    validate_compact,
)


def _read_config(root: Path = ROOT) -> bytes:
    return _read_repository_leaf(
        root,
        CONFIG_PATH,
        label="TDG8 successor config",
        maximum_bytes=1024 * 1024,
    )


def _read_result(root: Path = ROOT) -> bytes:
    return _read_repository_leaf(
        root,
        RESULT_PATH,
        label="TDG8 successor result",
        maximum_bytes=4 * 1024 * 1024,
    )


def _open_parent(root: Path, relative: str) -> tuple[int, str]:
    """Open a fixed result parent through a root-relative no-follow walk."""

    components = Path(relative).parts
    directory_fd = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        if not stat.S_ISDIR(os.fstat(directory_fd).st_mode):
            raise ValueError("TDG8 successor repository root is unsafe")
        for component in components[:-1]:
            child_fd = os.open(
                component,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=directory_fd,
            )
            if not stat.S_ISDIR(os.fstat(child_fd).st_mode):
                os.close(child_fd)
                raise ValueError("TDG8 successor result parent is unsafe")
            os.close(directory_fd)
            directory_fd = child_fd
        return directory_fd, components[-1]
    except BaseException:
        os.close(directory_fd)
        raise


def _publish_new_result(root: Path, raw: bytes) -> None:
    """Create the compact result once; never replace an existing authority."""

    directory_fd, leaf = _open_parent(root, RESULT_PATH)
    result_fd = -1
    created = False
    try:
        result_fd = os.open(
            leaf,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o644,
            dir_fd=directory_fd,
        )
        created = True
        view = memoryview(raw)
        while view:
            written = os.write(result_fd, view)
            if written <= 0:
                raise OSError("TDG8 successor result write made no progress")
            view = view[written:]
        os.fsync(result_fd)
        os.close(result_fd)
        result_fd = -1
        os.fsync(directory_fd)
    except BaseException:
        if result_fd != -1:
            os.close(result_fd)
        if created:
            try:
                os.unlink(leaf, dir_fd=directory_fd)
            except OSError:
                pass
        raise
    finally:
        os.close(directory_fd)


def _write(root: Path = ROOT) -> bytes:
    config_raw = _read_config(root)
    result_raw = canonical_result(
        build_prelaunch(config_raw, root, observe_destination=True)
    )
    _publish_new_result(root, result_raw)
    # The publication check is compact-only.  It does not repeat the one-time
    # destination observation made by build_prelaunch.
    published = _read_result(root)
    validate_compact(config_raw, published)
    if published != result_raw:
        raise ValueError("published TDG8 successor result bytes differ")
    return published


def _verify_compact(root: Path = ROOT) -> bytes:
    config_raw = _read_config(root)
    result_raw = _read_result(root)
    validate_compact(config_raw, result_raw)
    return result_raw


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify-compact", action="store_true")
    args = parser.parse_args(argv)

    if args.write:
        result_raw = _write(ROOT)
    else:
        result_raw = _verify_compact(ROOT)

    print(
        json.dumps(
            {
                "artifact_id": ARTIFACT_ID,
                "candidate_branch_opened": False,
                "compact_result_sha256": sha256(result_raw).hexdigest(),
                "destination_absence_observed": bool(args.write),
                "destination_absence_reobserved": False,
                "source_campaign_opened": False,
                "verified_compact": bool(args.verify_compact),
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
