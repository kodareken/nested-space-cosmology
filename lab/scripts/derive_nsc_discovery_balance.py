#!/usr/bin/env python3
"""Run or check the same-RK4 disjoint ledger; create new records exclusively."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

from recursive_horizons import nsc_discovery_balance as balance

EPISODE_DIR = balance.LAB / "results/development/nsc-discovery-episode-v1"
OUTPUT_DIR = balance.LAB / "results/development/nsc-discovery-balance-v1"


def _destination(path):
    path = Path(path).resolve()
    if OUTPUT_DIR.resolve() not in path.parents or path.suffix != ".json":
        raise PermissionError("new balance JSON records belong under " + str(OUTPUT_DIR))
    if path.exists():
        raise FileExistsError("refusing to replace an existing balance record: " + str(path))
    return path


def _create_only(path, text):
    path = _destination(path)
    if len(text.encode("utf-8")) > balance.CHUNK_LIMIT_BYTES:
        raise RuntimeError("balance record exceeds the 64 MiB chunk cap")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o444)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        handle.write(text)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="forecast without evolution, or validate --record")
    parser.add_argument("--run", "--evolve", dest="run", action="store_true", help="advance the sealed T=0.3 handoff")
    parser.add_argument("--target-time", type=float, default=balance.TARGET_TIME, help="coordinate endpoint (default 1; use 3 for the later station)")
    parser.add_argument("--spatial", action="store_true", help="also compare nf512 at dt=0.0005")
    parser.add_argument("--control", default="coupled", choices=("coupled", "frozen_geometry"))
    parser.add_argument("--write", "--output", dest="write", type=Path, help="exclusive new JSON path under the balance-v1 directory")
    parser.add_argument("--record", type=Path, help="with --check, authenticate and replay recorded scalar accounting")
    args = parser.parse_args(argv)
    if args.check and args.run:
        parser.error("choose --check or --run")
    if not args.check and not args.run and args.write is None:
        parser.error("choose --check, --run, or --write PATH")
    if args.record is not None and not args.check:
        parser.error("--record belongs to --check")
    if args.spatial and args.check:
        parser.error("--spatial belongs to a run")
    balance.acceptance_steps(balance.HANDOFF_TIME, args.target_time, balance.STEP_CAPS[0])
    destination = None if args.write is None else _destination(args.write)
    if args.record is not None:
        record = balance.verify_record(args.record)
    elif args.check:
        record = balance.check_record(target_time=args.target_time)
    else:
        record = balance.evolve_frozen_segment(spatial=args.spatial, control_mode=args.control, target_time=args.target_time)
    record["output_written"] = destination is not None
    text = json.dumps(balance.jsonable(record), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if destination is not None:
        if record.get("ok") is False:
            raise RuntimeError("refusing to write a failed balance check")
        _create_only(destination, text)
    sys.stdout.write(text)
    return 1 if record.get("ok") is False else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, PermissionError, FileExistsError, ValueError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
