#!/usr/bin/env python3
"""Prepare, seal a prediction, then independently measure physical width response."""
from __future__ import annotations
import argparse
import os
import sys

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from recursive_horizons import nsc_discovery_width_response as width


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    for mode in ("prepare", "predict", "measure", "run", "check"):
        modes.add_argument("--"+mode, action="store_true")
    parser.add_argument("--output", required=True)
    parser.add_argument("--nf", type=int, default=128)
    parser.add_argument("--duration", type=float, default=0.01)
    parser.add_argument("--step-cap", type=float, default=width.episode.MATCHED_STEP_CAP)
    parser.add_argument("--production", action="store_true")
    parser.add_argument("--cpu-budget", type=float, default=width.CPU_CAP)
    args = parser.parse_args(argv)
    try:
        if args.prepare or args.run:
            result = width.prepare(args.output, nf=args.nf, duration=args.duration, step_cap=args.step_cap,
                                   production=args.production, cpu_budget=args.cpu_budget)
        elif args.predict:
            result = width.predict(args.output)
        elif args.measure:
            result = width.measure(args.output)
        else:
            result = width.check(args.output)
        if args.run and result["status"] == "prepared":
            result = width.predict(args.output)
            if result["status"] == "predicted":
                result = width.measure(args.output)
        print(width.SCHEMA, result.get("status", "checked"), "cpu", result.get("aggregate_cpu_seconds"), flush=True)
        return 0
    except (ValueError, RuntimeError, PermissionError, FileExistsError, FileNotFoundError) as error:
        print(type(error).__name__+": "+str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
