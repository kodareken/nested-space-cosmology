#!/usr/bin/env python3
"""Construct once or compactly verify FGC-1-TDG9-UR1-PREF1."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from recursive_horizons.fgc.evolution import tdg9_ur1_pref1_binder as binder  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--output", default=binder.RESULT_PATH)
    arguments = parser.parse_args()
    config = (ROOT / binder.CONFIG_PATH).read_bytes()
    output = ROOT / arguments.output
    if arguments.live:
        result = binder.build_pref1_result(config, ROOT, live=True)
        output.write_bytes(binder.canonical_result(result))
    else:
        binder.validate_compact_result(config, output.read_bytes())
    print(
        "FGC-1-TDG9-UR1-PREF1 verified: u:R rows match TI2, D01/D12 zero "
        "lowers owned by coefficient_plus_arithmetic_debit_clips_positive_"
        "raw_maximum, classification remains order_inconclusive, no "
        "production method or physics claim"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
