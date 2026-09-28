#!/usr/bin/env python3
"""Write once or compactly verify FGC-1-TDG10-QA1-FRZ1."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons.fgc.evolution import tdg10_qa1_authority as authority  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--verify-compact", action="store_true")
    arguments = parser.parse_args(argv)
    config = authority.read_leaf(ROOT, authority.CONFIG_PATH)
    path = ROOT / authority.RESULT_PATH
    if arguments.write:
        snapshot = authority._snapshot_store(ROOT)
        raw = authority.canonical_pretty(
            authority.build_prelaunch(config, ROOT, store_snapshot=snapshot)
        )
        path.write_bytes(raw)
        mode = "write"
    else:
        raw = authority.read_leaf(ROOT, authority.RESULT_PATH)
        authority.validate_compact(config, raw)
        mode = "verify_compact"
    print(
        json.dumps(
            {
                "artifact_id": authority.ARTIFACT_ID,
                "mode": mode,
                "result_sha256": sha256(raw).hexdigest(),
                "restored_predecessor_fingerprint_count": 1,
                "shadow_executed": False,
                "shadow_proposals_constructed": 0,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
