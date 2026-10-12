#!/usr/bin/env python3
"""Replay corrected region contours from original station chunks.

``--check`` hashes the current source closure and the sealed nf256/nf512 v1
records. It does not parse those records or any trajectory, and it writes
nothing. ``--read`` analyzes the given station chunks with the nodal peak
rate. Before a file is created, each historical case and ordinal must match
the chunk hash stored in the sealed v1 record. The write is exclusive.
A file above 64 MiB is refused. Sealed v1 speeds and moving ledgers are
not rewritten.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

from recursive_horizons.nsc_discovery_regions import (
    MAX_OUTPUT_BYTES,
    assert_creation_bindings,
    jsonable,
    successor_check,
    successor_from_stations,
)


def _create_only(path, text):
    encoded = text.encode("utf-8")
    if len(encoded) > MAX_OUTPUT_BYTES:
        raise ValueError("refusing a region successor larger than 64 MiB")
    path = Path(path)
    if path.exists():
        raise FileExistsError("refusing to replace an existing region successor: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    with os.fdopen(descriptor, "wb") as handle:
        handle.write(encoded)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Corrected discovery-region successor")
    parser.add_argument("--check", action="store_true", help="hash the source closure and sealed v1 records; write nothing")
    parser.add_argument("--read", type=Path, help="episode directory containing the original station chunks")
    parser.add_argument("--case", help="case id inside that directory")
    parser.add_argument("--ordinal", type=int, help="one station ordinal; default is every station chunk")
    parser.add_argument("--output", type=Path, help="creation-only JSON path")
    args = parser.parse_args(argv)
    if args.check == bool(args.read):
        parser.error("choose exactly one of --check or --read")
    if args.output is not None and args.read is None:
        parser.error("--output belongs to --read")
    if args.read is not None and not args.case:
        parser.error("--read requires --case")
    if args.check:
        record = successor_check()
    else:
        if args.output is not None and Path(args.output).resolve().is_relative_to(Path(args.read).resolve()):
            parser.error("--output must stay outside the episode directory")
        record = successor_from_stations(args.read, args.case, ordinal=args.ordinal)
        record["output_written"] = args.output is not None
        assert_creation_bindings(record)
        record["frozen_by_root"] = False
        record["production_campaign"] = False
    text = json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n"
    if args.output is not None:
        _create_only(args.output, text)
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
