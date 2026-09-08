#!/usr/bin/env python3
"""Reproduce the common-action normalization-profile comparison."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_measure_normalization import build_record

OUTPUT = ROOT / "results/nsc-10-measure-normalization.json"
ABS_TOL, REL_TOL = 2e-7, 2e-8


def compare_record(expected, actual, path="$"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise RuntimeError(f"keys differ at {path}")
        for key in expected:
            compare_record(expected[key], actual[key], path + "/" + str(key))
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise RuntimeError(f"list differs at {path}")
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare_record(left, right, f"{path}/{i}")
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (float, int)) or not np.isclose(
            expected, actual, atol=ABS_TOL, rtol=REL_TOL
        ):
            raise RuntimeError(f"numeric field differs at {path}: {expected!r} != {actual!r}")
    elif type(expected) is not type(actual) or expected != actual:
        raise RuntimeError(f"exact field differs at {path}: {expected!r} != {actual!r}")


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
        print("Normalization-profile identity, metric response and all-field hashes reproduced.")
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
