#!/usr/bin/env python3
"""Write once or compactly verify FGC-1-TDG9-LOC2-FRZ1."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg9_loc2_authority as authority  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--verify-compact", action="store_true")
    arguments = parser.parse_args(argv)
    config = (ROOT / authority.CONFIG_PATH).read_bytes()
    path = ROOT / authority.RESULT_PATH
    if arguments.write:
        sys.path.insert(0, str(ROOT))
        from scripts.run_fgc_tdg9_loc2 import prospective_coefficient_depth_audit

        audit = prospective_coefficient_depth_audit(ROOT)
        raw = authority.canonical_pretty(authority.build_prelaunch(config, ROOT, audit))
        path.write_bytes(raw)
        mode = "write"
    else:
        raw = path.read_bytes()
        authority.validate_compact(config, raw)
        mode = "verify_compact"
    print(json.dumps({"artifact_id": authority.ARTIFACT_ID, "mode": mode, "result_sha256": sha256(raw).hexdigest(), "diagnostic_executed": False}, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
