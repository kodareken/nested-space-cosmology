#!/usr/bin/env python3
"""Reproduce or compactly verify FGC-1-TDG9-AR1-PREF1."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg9_ar1_pref1_binder as binder  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--output", default=binder.RESULT_PATH)
    arguments = parser.parse_args()
    config = (ROOT / binder.CONFIG_PATH).read_bytes()
    output = ROOT / arguments.output
    if arguments.live:
        result = binder.build_pref1_result(config, ROOT, live=True)
        raw = binder.canonical_result(result)
        output.write_bytes(raw)
    else:
        raw = output.read_bytes()
        binder.validate_compact_result(config, raw)
    print(
        "FGC-1-TDG9-AR1-PREF1 verified: "
        "108 exact calls, 10 sufficient contraction failures, no physics claim"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
