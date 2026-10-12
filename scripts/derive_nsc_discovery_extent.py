#!/usr/bin/env python3
"""Preview or create a bounded period8/12 physical-domain intervention."""
from __future__ import annotations
import argparse
import json
import os
import sys
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")
from recursive_horizons import nsc_discovery_extent as extent


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    for name in ("prepare", "materialize", "run", "check", "prepare-resume"):
        mode.add_argument("--"+name, action="store_true")
    parser.add_argument("--output")
    parser.add_argument("--reference", help="immutable extent-v1 directory, only for --prepare-resume")
    parser.add_argument("--production", action="store_true")
    parser.add_argument("--periods", type=float, nargs="+", default=[8., 12.])
    parser.add_argument("--points-per-unit", type=int, default=32)
    parser.add_argument("--prefix", type=float, default=0.3)
    parser.add_argument("--stations", type=float, nargs="+", default=[8., 12.])
    parser.add_argument("--step-cap", type=float, default=0.0005)
    parser.add_argument("--cadence", type=float, default=0.1)
    parser.add_argument("--cpu-budget", type=float, default=21600.)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--period16-admission", help="root-reviewed JSON stating extent_discriminates:true")
    args = parser.parse_args(argv)
    try:
        admission = json.loads(open(args.period16_admission).read()) if args.period16_admission else None
        config = dict(periods=args.periods, points_per_unit=args.points_per_unit, prefix=args.prefix,
                      stations=args.stations, step_cap=args.step_cap, cadence=args.cadence,
                      cpu_budget=args.cpu_budget, period16_admission=admission)
        if not any((args.prepare, args.materialize, args.run, args.check, args.prepare_resume)):
            result = extent.settings(**config)
            print(extent.SCHEMA, "preview", [(L, int(L*args.points_per_unit)) for L in result["periods"]], flush=True)
            return 0
        if not args.output:
            parser.error("--output is required for a stage")
        if args.prepare_resume:
            if not args.reference:
                parser.error("--reference is required for --prepare-resume")
            result = extent.prepare_resume(args.output, args.reference)
        elif args.prepare:
            result = extent.prepare(args.output, **config)
        elif args.materialize:
            result = extent.materialize(args.output, production=args.production)
        elif args.run:
            result = extent.run(args.output, production=args.production, workers=args.workers, max_steps=args.max_steps)
        else:
            result = extent.check(args.output)
        print(extent.SCHEMA, result.get("stage", "checked"), "cpu", result.get("aggregate_cpu_seconds"), flush=True)
        return 0
    except (ValueError, RuntimeError, PermissionError, FileNotFoundError, FileExistsError) as error:
        print(type(error).__name__+": "+str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
