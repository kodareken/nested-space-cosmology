#!/usr/bin/env python3
"""Reproduce the bounded free-Dirac normalization-scale flow."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_normalization_flow import (  # noqa: E402
    build_record, compare_record,
)

OUTPUT = ROOT / "results/nsc-12-normalization-flow.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is not None and args.output.exists():
        raise FileExistsError("refusing to overwrite record")
    record = build_record()
    if args.check:
        compare_record(json.loads(OUTPUT.read_text()), record)
        print("Normalization-scale flow, matched coefficient shifts and all-field hashes reproduced.")
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(record, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
