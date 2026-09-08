#!/usr/bin/env python3
"""Compute the seam-free geometric gap and required 4D Einstein stress.

Exclusive --output writes results/nsc-4-smooth-geometry.json. --check reproduces
every recorded field. This runner does not edit older NSC modules or records.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_smooth_geometry import (  # noqa: E402
    DEFAULT_RESULT,
    build_record,
    jsonable,
)
from check_nsc_scale_closure import compare  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = jsonable(build_record())
    if args.check:
        compare(json.loads(DEFAULT_RESULT.read_text()), record)
        print(
            "Seam-free geometric gap and required Einstein stress reproduced; "
            "every field checked."
        )
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(record, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
