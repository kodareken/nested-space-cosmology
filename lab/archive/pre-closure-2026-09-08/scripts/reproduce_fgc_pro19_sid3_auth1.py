#!/usr/bin/env python3
"""Build or compactly verify the PRO19 SID3 committed-image authority.

This reproducer consumes only the tracked SID3 TOML contract, execution-
closure JSON, and compact result.  It has no live-preflight mode and never
opens the authenticated run store.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import stat
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto19_resume_authority import (  # noqa: E402
    ARTIFACT_ID,
    CLOSURE_PATH,
    CONFIG_PATH,
    RESULT_PATH,
    build_sid3_result,
    canonical_result,
    validate_sid3_result,
)


def _read_unique_regular(relative: str) -> bytes:
    """Read one repository leaf without following links or hard links."""

    path = ROOT / relative
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
        raise ValueError(f"unsafe SID3 compact leaf: {relative}")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        observed = os.fstat(descriptor)
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_dev != before.st_dev
            or observed.st_ino != before.st_ino
            or observed.st_size != before.st_size
        ):
            raise ValueError(f"SID3 compact leaf changed during read: {relative}")
        chunks: list[bytes] = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        after = os.fstat(descriptor)
        if (
            after.st_dev != observed.st_dev
            or after.st_ino != observed.st_ino
            or after.st_size != observed.st_size
        ):
            raise ValueError(f"SID3 compact leaf changed during read: {relative}")
        return b"".join(chunks)
    finally:
        os.close(descriptor)


def expected_result() -> bytes:
    config_raw = _read_unique_regular(CONFIG_PATH)
    closure_raw = _read_unique_regular(CLOSURE_PATH)
    return canonical_result(build_sid3_result(config_raw, closure_raw))


def _write_result() -> None:
    destination = ROOT / RESULT_PATH
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    if temporary.exists():
        raise ValueError("stale SID3 temporary result exists")
    try:
        with temporary.open("xb") as handle:
            handle.write(expected_result())
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        directory = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def _verify_result() -> None:
    config_raw = _read_unique_regular(CONFIG_PATH)
    closure_raw = _read_unique_regular(CLOSURE_PATH)
    result_raw = _read_unique_regular(RESULT_PATH)
    validate_sid3_result(config_raw, closure_raw, result_raw)
    if result_raw != canonical_result(build_sid3_result(config_raw, closure_raw)):
        raise ValueError("PRO19 SID3 compact result differs")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--build", action="store_true")
    mode.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)

    if args.build:
        _write_result()
    else:
        _verify_result()

    print(
        json.dumps(
            {
                "artifact_id": ARTIFACT_ID,
                "candidate_branch_opened": False,
                "compact_only": True,
                "physical_state_advanced": False,
                "store_mutated": False,
                "trajectory_resumed": False,
                "verified": args.verify,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
