#!/usr/bin/env python3
"""Reproduce every field of the domain-aware full-spinor boundary experiment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_chiral_boundary import DEFAULT_RESULT, build_record, compare, jsonable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--output", type=Path)
    options.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.output is not None and args.output.exists():
        raise FileExistsError(f"refusing to overwrite existing result: {args.output}")
    record = jsonable(build_record())
    if args.check:
        compare(json.loads(DEFAULT_RESULT.read_bytes()), record)
        print("Full-spinor boundary response reproduced; every numeric, exact, convention, and source-hash field checked.")
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(record, stream, sort_keys=True, separators=(",", ":"), allow_nan=False)
        print(json.dumps({"output": str(args.output), "schema": record["schema"], "gate": record["gate"]}, indent=2))
    else:
        print(json.dumps(record, sort_keys=True, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
