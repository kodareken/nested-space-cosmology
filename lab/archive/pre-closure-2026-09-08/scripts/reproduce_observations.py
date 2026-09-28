#!/usr/bin/env python3
"""Write the compact GW170817/GRB 170817A timing-input reconstruction."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.observations import reproduce_gw170817_timing  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY / "results" / "gw170817-timing.json",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    record = reproduce_gw170817_timing()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    validation = record["validation"]
    print(f"wrote {output}")
    print(
        "reconstructed geocentric delay: "
        f"{validation['reconstructed_geocentric_delay_seconds']:.6f} s"
    )
    print(
        "Scope: published timing-input reconstruction only; not raw strain/light-curve "
        "reanalysis or a measurement of an inaccessible parent speed."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
