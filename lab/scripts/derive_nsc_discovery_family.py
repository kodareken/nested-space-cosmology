#!/usr/bin/env python3
"""Source-family specification, short or production handoff, and one-pool continuation.

``--prepare`` writes the admissible family and the six exploration cases.
It does not solve. ``--materialize`` runs the dense initial solve and
``propagate_baseline``. Duration ``0.3`` requires ``--production`` and is
the root executor's call. ``--run`` opens one episode pool on the
propagated handoff. ``--confirm`` writes a separate directory and does
not regenerate members whose physical signs were not distinct.
"""
from __future__ import annotations

import argparse
import os
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

from recursive_horizons.nsc_discovery_episode import (
    DEFAULT_CPU_BUDGET_SECONDS,
    DEFAULT_FORECAST_FACTOR,
    DEFAULT_MEMORY_BYTES,
    DEFAULT_WORKERS,
    MATCHED_STEP_CAP,
)
from recursive_horizons.nsc_discovery_family import (
    AGGREGATE_CPU_BUDGET_SECONDS,
    CASE_COUNT,
    CHUNK_LIMIT_BYTES,
    EXPLORATION_NF,
    HANDOFF_TIME,
    OCCUPATION_TRACE,
    NONPRODUCTION_DURATION,
    campaign_forecast,
    check,
    family_catalog,
    materialize_directory,
    prepare_confirmation,
    prepare_specification,
    run,
    specification,
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--catalog", action="store_true")
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--confirm", action="store_true")
    mode.add_argument("--materialize", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--output", default=None)
    parser.add_argument("--all-members", action="store_true", help="prepare all 15 members in the same six-worker pool")
    parser.add_argument("--exploration", default=None, help="read-only exploration directory for --confirm")
    parser.add_argument("--nf", type=int, default=EXPLORATION_NF)
    parser.add_argument("--duration", type=float, default=NONPRODUCTION_DURATION)
    parser.add_argument("--production", action="store_true")
    parser.add_argument("--step-cap", type=float, default=MATCHED_STEP_CAP)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--cpu-budget", type=float, default=DEFAULT_CPU_BUDGET_SECONDS)
    parser.add_argument("--forecast-factor", type=float, default=DEFAULT_FORECAST_FACTOR)
    parser.add_argument("--memory-gib", type=float, default=DEFAULT_MEMORY_BYTES / 1024 ** 3)
    parser.add_argument("--max-steps", type=int, default=None)
    parser.add_argument("--backend", choices=("auto", "dense", "fft"), default="auto")
    args = parser.parse_args(argv)
    try:
        return _dispatch(parser, args)
    except (FileExistsError, FileNotFoundError, PermissionError, RuntimeError, ValueError) as error:
        print(type(error).__name__ + ": " + str(error), file=sys.stderr)
        return 2


def _dispatch(parser, args):
    if args.catalog:
        spec = specification()
        members, exploration = family_catalog(args.nf)
        forecast = campaign_forecast(CASE_COUNT, args.nf, args.step_cap)
        print(
            spec["schema"],
            "members", len(members),
            "cases", len(exploration),
            "trace", OCCUPATION_TRACE,
            "budget", AGGREGATE_CPU_BUDGET_SECONDS,
            "chunk", CHUNK_LIMIT_BYTES,
            "pools", spec["coordinator_pools"],
            "forecast_without_inits", forecast["forecast_cpu_seconds_without_inits"],
            flush=True,
        )
        return 0
    if not args.output:
        parser.error("--output is required")
    if args.prepare:
        record = prepare_specification(args.output, nf=args.nf, all_members=args.all_members)
        print(record["stage"], "cases", len(record["cases"]), "materialized", record["materialized"], flush=True)
        return 0
    if args.confirm:
        if not args.exploration:
            parser.error("--confirm requires --exploration")
        record = prepare_confirmation(args.output, args.exploration)
        print(
            record["stage"],
            "selected", len(record["selected"]["cases"]),
            "forced", record["forced_regeneration"],
            flush=True,
        )
        return 0
    if args.materialize:
        manifest = materialize_directory(
            args.output,
            duration=args.duration,
            step_cap=args.step_cap,
            production=args.production,
        )
        print(
            manifest["status"],
            "cases", len(manifest["cases"]),
            "initial_cpu", manifest["initial_cpu_seconds"],
            "continuation_budget", manifest["cpu_budget_seconds"],
            flush=True,
        )
        return 0
    if args.check:
        report = check(args.output)
        print(report["schema"], "ok", report["ok"], "materialized", report["materialized"], flush=True)
        return 0 if report["ok"] else 2
    result = run(
        args.output,
        workers=args.workers,
        cpu_budget_seconds=args.cpu_budget,
        forecast_factor=args.forecast_factor,
        memory_limit_bytes=int(args.memory_gib * 1024 ** 3),
        max_steps=args.max_steps,
        backend=args.backend,
    )
    print(
        result["status"],
        "pools", result["coordinator_pools"],
        "aggregate_cpu", result["aggregate_cpu_seconds"],
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
