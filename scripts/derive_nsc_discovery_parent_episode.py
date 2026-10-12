#!/usr/bin/env python3
"""Pure parent-episode plan, or an explicit frozen prepare / run / check."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

from recursive_horizons import nsc_discovery_parent_episode as parent


def station_list(value):
    try:
        return parent.normalize_stations(value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--prepare", action="store_true", help="write a new campaign from OWNER1 output")
    modes.add_argument("--run", action="store_true", help="root launch: one pool, native FFT, scaled steps")
    modes.add_argument("--check", action="store_true", help="read-only reconstruction and binding check")
    modes.add_argument("--continue-from", type=Path, help="exact completed parent campaign successor; requires explicit remaining CPU")
    parser.add_argument("--source", type=Path, help="one parent.json directory, or a directory of those records")
    parser.add_argument("--output", type=Path, default=parent.OUTPUT)
    parser.add_argument("--producer-commit")
    parser.add_argument("--workers", type=int, default=parent.MAX_WORKERS)
    parser.add_argument("--cpu-budget", type=float)
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--confirm", action="store_true", help="prepare the nf256 dt-cap 0.0005 cases")
    parser.add_argument("--stations", type=station_list, help="prepare/preview ordered positive times, e.g.1,3; run uses frozen stored targets")
    parser.add_argument("--control-mode", choices=("coupled","frozen_geometry","source_free"),
                        help="prepare/preview control; frozen geometry is prescribed")
    parser.add_argument("--step-cap", type=float, help="prepare/preview positive finite cap override")
    args = parser.parse_args(argv)
    if args.continue_from is not None and (args.cpu_budget is None or args.producer_commit is None):
        parser.error("--continue-from requires explicit --cpu-budget remaining and --producer-commit")
    if args.continue_from is not None and (args.control_mode is not None or args.step_cap is not None or args.source is not None):
        parser.error("continuation preserves its source/control/cap; these cannot be overridden")
    budget = parent.CPU_BUDGET_SECONDS if args.cpu_budget is None else args.cpu_budget
    if (args.stations is not None or args.control_mode is not None or args.step_cap is not None) and (args.run or args.check):
        parser.error("--stations/--control-mode/--step-cap belong to --prepare or preview; run/check use stored configuration")
    if args.continue_from is not None:
        result = parent.continue_parent(args.continue_from,args.output,
            stations=args.stations or (1.25,2.25,3.),cpu_budget=budget,producer_commit=args.producer_commit)
    elif args.prepare:
        result = parent.prepare(args.source, args.output, execute=True, producer_commit=args.producer_commit,
                                cpu_budget=budget, confirm=args.confirm, stations=args.stations,
                                control_mode=args.control_mode or "coupled",step_cap=args.step_cap)
    elif args.run:
        result = parent.run(args.output, workers=args.workers, cpu_budget=budget, max_steps=args.max_steps)
    elif args.check:
        result = parent.check(args.output)
    else:
        result = parent.plan(stations=args.stations,control_mode=args.control_mode or "coupled",step_cap=args.step_cap)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
