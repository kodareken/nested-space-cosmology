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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--prepare", action="store_true", help="write a new campaign from OWNER1 output")
    modes.add_argument("--run", action="store_true", help="root launch: one pool, native FFT, scaled steps")
    modes.add_argument("--check", action="store_true", help="read-only reconstruction and binding check")
    parser.add_argument("--source", type=Path, help="one parent.json directory, or a directory of those records")
    parser.add_argument("--output", type=Path, default=parent.OUTPUT)
    parser.add_argument("--producer-commit")
    parser.add_argument("--workers", type=int, default=parent.MAX_WORKERS)
    parser.add_argument("--cpu-budget", type=float, default=parent.CPU_BUDGET_SECONDS)
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--confirm", action="store_true", help="prepare the nf256 dt-cap 0.0005 cases")
    args = parser.parse_args(argv)
    if args.prepare:
        result = parent.prepare(args.source, args.output, execute=True, producer_commit=args.producer_commit,
                                cpu_budget=args.cpu_budget, confirm=args.confirm)
    elif args.run:
        result = parent.run(args.output, workers=args.workers, cpu_budget=args.cpu_budget, max_steps=args.max_steps)
    elif args.check:
        result = parent.check(args.output)
    else:
        result = parent.plan()
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
