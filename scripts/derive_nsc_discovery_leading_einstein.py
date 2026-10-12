#!/usr/bin/env python3
"""Pure preflight or explicit leading-Einstein preparation/run/read-only check."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

from recursive_horizons import nsc_discovery_leading_einstein as leading


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--prepare", action="store_true", help="explicitly prepare eight matched coupled/control checkpoints")
    modes.add_argument("--run", action="store_true", help="execute the one existing runner/pool")
    modes.add_argument("--check", action="store_true", help="authenticate and replay constraints without stepping")
    parser.add_argument("--initial", type=Path, default=leading.INITIAL)
    parser.add_argument("--output", type=Path, default=leading.OUTPUT)
    parser.add_argument("--producer-commit", help="frozen source revision required for preparation")
    parser.add_argument("--cpu-budget", type=float, default=600.)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-steps", type=int, help="explicit bounded diagnostic step count")
    args = parser.parse_args(argv)
    if args.prepare:
        result = leading.prepare(args.initial, args.output, execute=True,
                                 producer_commit=args.producer_commit, cpu_budget=args.cpu_budget)
    elif args.run:
        result = leading.run(args.output, workers=args.workers, cpu_budget=args.cpu_budget, max_steps=args.max_steps)
    elif args.check:
        result = leading.check(args.output)
    else:
        result = leading.prepare(args.initial, args.output, cpu_budget=args.cpu_budget)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
