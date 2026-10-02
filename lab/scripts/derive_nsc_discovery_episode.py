#!/usr/bin/env python3
"""Prepare, check, or run the first separated-pair continuation batch.

Stage 0 writes six exact T=0.3 handoff chunks. New successor files may be
created under ``results/development/nsc-discovery-...``. Existing sealed
records are never overwritten. ``--check`` reloads chunks and does not
evolve. ``--run`` uses one pool, backend ``auto``, 6 workers, and 21600
seconds of summed child-process CPU. This import does not start the batch.
"""
from __future__ import annotations

import argparse
import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

from recursive_horizons.nsc_discovery_episode import (
    DEFAULT_CPU_BUDGET_SECONDS,
    DEFAULT_FORECAST_FACTOR,
    DEFAULT_MEMORY_BYTES,
    DEFAULT_STATIONS,
    DEFAULT_WORKERS,
    check,
    prepare,
    run,
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--output", required=True,
                        help="tmp directory, or lab/results/development/nsc-discovery-episode-v1")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--cpu-budget", type=float, default=DEFAULT_CPU_BUDGET_SECONDS)
    parser.add_argument("--forecast-factor", type=float, default=DEFAULT_FORECAST_FACTOR)
    parser.add_argument("--memory-gib", type=float, default=DEFAULT_MEMORY_BYTES / 1024 ** 3)
    parser.add_argument("--stations", type=float, nargs="*", default=None)
    parser.add_argument("--backend", choices=("auto", "dense", "fft"), default="auto")
    parser.add_argument("--fft", action="store_true", help="force the FFT carrier")
    parser.add_argument("--dense", action="store_true", help="force the dense evolution carrier")
    args = parser.parse_args(argv)
    if args.fft and args.dense:
        parser.error("--fft and --dense cannot both be set")
    stations = None if args.stations is None else tuple(args.stations)
    if args.prepare and stations is None:
        stations = DEFAULT_STATIONS
    if args.check:
        report = check(args.output)
        print(report["schema"], "ok", report["ok"], "commits", len(report["checked"]), flush=True)
        return 0
    if args.prepare:
        manifest = prepare(args.output, stations=stations)
        print(manifest["status"], "cases", len(manifest["cases"]), "evolved", manifest["evolved"], flush=True)
        return 0
    backend = "fft" if args.fft else "dense" if args.dense else args.backend
    manifest = run(args.output, workers=args.workers, cpu_budget_seconds=args.cpu_budget,
                   forecast_factor=args.forecast_factor,
                   memory_limit_bytes=int(args.memory_gib * 1024 ** 3), backend=backend)
    print(manifest["status"], "child_cpu", manifest["child_cpu_seconds"],
          "pools", manifest["executor_pools"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
