#!/usr/bin/env python3
"""Read stage-4 region rows from saved station chunks. No evolution.

``--check`` hashes producers and sealed inputs and writes nothing.
``--read`` analyzes station chunks from an episode directory. ``--output``
creates one new JSON file and refuses to replace an existing file.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

from recursive_horizons.nsc_discovery_regions import check_record, jsonable, read_station_regions


def _create_only(path, text):
    path = Path(path)
    if path.exists():
        raise FileExistsError("refusing to replace an existing region record: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
    descriptor = os.open(path, flags)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        handle.write(text)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Stage-4 discovery region consumer")
    parser.add_argument("--check", action="store_true", help="hash inputs and producers; write nothing")
    parser.add_argument("--read", type=Path, help="episode directory containing station chunks")
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
        record = check_record()
    else:
        if args.output is not None and Path(args.output).resolve().is_relative_to(Path(args.read).resolve()):
            parser.error("--output must stay outside the episode directory")
        record = read_station_regions(args.read, args.case, ordinal=args.ordinal)
        record["output_written"] = args.output is not None
    text = json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n"
    if args.output is not None:
        _create_only(args.output, text)
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
