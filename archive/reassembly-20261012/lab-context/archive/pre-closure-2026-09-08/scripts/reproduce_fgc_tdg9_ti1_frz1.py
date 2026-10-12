#!/usr/bin/env python3
"""Write once or compactly verify FGC-1-TDG9-TI1-FRZ1."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts import run_fgc_tdg9_loc1 as loc1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ti1_authority as authority  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--verify-compact", action="store_true")
    arguments = parser.parse_args(argv)
    config = authority.read_leaf(ROOT, authority.CONFIG_PATH)
    path = ROOT / authority.RESULT_PATH
    if arguments.write:
        snapshot = loc1._snapshot_store(ROOT)
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
                "shadow_executed": False,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
