#!/usr/bin/env python3
"""Build or verify the read-only PRO19 SID1 authority certificate."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto19_sid1_authority import (  # noqa: E402
    ARTIFACT_ID, CONFIG_PATH, RESULT_PATH, build_sid1_result, derive_live_anchor,
    expected_anchor,
)


def canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def expected_result(repository: Path) -> bytes:
    return canonical(build_sid1_result((repository / CONFIG_PATH).read_bytes(), derive_live_anchor(repository)))


def expected_compact_result(repository: Path) -> bytes:
    """Rebuild tracked bytes without reopening the one-time live anchor."""
    return canonical(
        build_sid1_result(
            (repository / CONFIG_PATH).read_bytes(),
            expected_anchor(),
        )
    )


def write_result(repository: Path) -> None:
    destination = repository / RESULT_PATH
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    temporary.write_bytes(expected_result(repository))
    os.replace(temporary, destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--verify-compact", action="store_true")
    args = parser.parse_args()
    if args.write:
        expected = expected_result(ROOT)
        write_result(ROOT)
    else:
        expected = (
            expected_result(ROOT)
            if args.verify
            else expected_compact_result(ROOT)
        )
        if (ROOT / RESULT_PATH).read_bytes() != expected:
            raise SystemExit("PRO19 SID1 result differs")
    print(json.dumps({"artifact_id": ARTIFACT_ID, "verified": bool(args.verify or args.verify_compact), "live_anchor_verified": bool(args.write or args.verify), "state_advanced": False, "store_mutated": False, "candidate_branch_opened": False}, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
