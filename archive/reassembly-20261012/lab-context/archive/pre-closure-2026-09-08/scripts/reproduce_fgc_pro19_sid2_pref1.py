#!/usr/bin/env python3
"""Build or verify the independent PRO19 SID2 post-recovery binder."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto19_sid2_binder import (  # noqa: E402
    ARTIFACT_ID,
    CONFIG_PATH,
    RESULT_PATH,
    build_sid2_result,
    canonical_result,
    validate_compact_result,
)


def live_result_bytes(repository: Path) -> bytes:
    config_raw = (repository / CONFIG_PATH).read_bytes()
    return canonical_result(build_sid2_result(config_raw, repository))


def write_result(repository: Path) -> None:
    destination = repository / RESULT_PATH
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    if temporary.exists():
        raise SystemExit("stale SID2 temporary result exists")
    with temporary.open("xb") as handle:
        handle.write(live_result_bytes(repository))
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--verify-compact", action="store_true")
    args = parser.parse_args()

    config_raw = (ROOT / CONFIG_PATH).read_bytes()
    result_path = ROOT / RESULT_PATH
    if args.write:
        write_result(ROOT)
        live_verified = True
    elif args.verify:
        if result_path.read_bytes() != live_result_bytes(ROOT):
            raise SystemExit("PRO19 SID2 result differs from live recovered boundary")
        live_verified = True
    else:
        validate_compact_result(config_raw, result_path.read_bytes())
        live_verified = False

    print(
        json.dumps(
            {
                "artifact_id": ARTIFACT_ID,
                "candidate_branch_opened": False,
                "live_boundary_verified": live_verified,
                "physical_state_advanced": False,
                "resume_authorized": False,
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
