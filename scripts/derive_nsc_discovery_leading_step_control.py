#!/usr/bin/env python3
"""Pure preflight or explicit exact-handoff step-control successor/check/run."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

from recursive_horizons import nsc_discovery_leading_step_control as step


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--prepare", action="store_true", help="copy exact parent handoffs into a new successor")
    modes.add_argument("--run", action="store_true", help="reuse the existing single leading pool with scaled step control")
    modes.add_argument("--check", action="store_true", help="read-only sources/prefixes/handoffs/constraints check")
    modes.add_argument("--paired-check", action="store_true", help="bounded full-step/two-half-step endpoint check")
    parser.add_argument("--source", type=Path, default=step.SOURCE)
    parser.add_argument("--output", type=Path, default=step.OUTPUT)
    parser.add_argument("--producer-commit")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cpu-budget", type=float, default=165.)
    parser.add_argument("--max-steps", type=int)
    args = parser.parse_args(argv)
    if args.prepare:
        result = step.prepare(args.source, args.output, execute=True,
                              producer_commit=args.producer_commit, cpu_budget=args.cpu_budget)
    elif args.run:
        result = step.run(args.output, workers=args.workers, cpu_budget=args.cpu_budget, max_steps=args.max_steps)
    elif args.check:
        result = step.check(args.output)
    elif args.paired_check:
        result = step.paired_step_check(args.source)
    else:
        result = step.prepare(args.source, args.output, cpu_budget=args.cpu_budget)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
