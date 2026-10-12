"""Coupled-memory segment writer and runner.

``--write`` creates the opening JSON+NPZ, producer hashes, and CPU forecast
for a segment that starts at the stored T=0.3 state. The full segment is
duration 0.7, ending at T=1. ``--run`` evolves a prepared campaign inside
the CPU budget. ``--check`` reloads those bytes and replays forces without
taking a step. ``--short`` remains the bounded 0.01 comparison.

Root is the executor for the full segment. This file does not start that
segment unless ``--run`` is passed.
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_discovery_coupled_memory as memory


def _fermions(text):
    values = tuple(int(piece) for piece in str(text).split(",") if piece)
    if not values or any(value not in memory.MATCHED_CASES for value in values):
        raise ValueError("fermions must be chosen from 128 and 256")
    return values


def _print(record):
    sys.stdout.write(memory.closure_bytes(record))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Causal coupled-memory segment from the stored T=0.3 state.")
    parser.add_argument("--check", action="store_true", help="reload a campaign, or print the pin closure")
    parser.add_argument("--short", action="store_true", help="bounded 0.01 full-versus-split comparison")
    parser.add_argument("--write", action="store_true", help="write opening JSON+NPZ and the CPU forecast")
    parser.add_argument("--run", action="store_true", help="evolve a prepared campaign inside the CPU budget")
    parser.add_argument("--duration", type=float, default=None)
    parser.add_argument("--step-cap", type=float, default=memory.DEFAULT_STEP_CAP)
    parser.add_argument("--output", default=None)
    parser.add_argument("--cpu-budget", type=float, default=memory.DEFAULT_CPU_BUDGET_SECONDS)
    parser.add_argument("--forecast-factor", type=float, default=memory.DEFAULT_FORECAST_FACTOR)
    parser.add_argument("--fermions", default="128,256")
    args = parser.parse_args(argv)
    fermions = _fermions(args.fermions)
    if args.short and (args.write or args.run):
        raise ValueError("the short check does not write or run the segment")
    if args.short:
        duration = memory.SHORT_DURATION if args.duration is None else float(args.duration)
        if duration > memory.MAX_SHORT_DURATION:
            raise ValueError("the short check stops at 0.02; use --write/--run for the 0.7 segment")
        record = memory.dependency_closure()
        record["short_check_ran"] = True
        record["short_check_duration"] = duration
        with threadpool_limits(limits=1), backend.fft_thread_limit(1):
            record["short_checks"] = [memory.short_coupled_check(nf, duration) for nf in fermions]
        _print(record)
        return 0
    if args.write or args.run or (args.check and args.output):
        if args.output is None:
            raise ValueError("a campaign needs --output outside the repository")
        if args.write:
            duration = memory.SEGMENT_DURATION if args.duration is None else float(args.duration)
            manifest = memory.prepare_campaign(
                args.output, fermions=fermions, duration=duration, step_cap=args.step_cap,
                cpu_budget_seconds=args.cpu_budget, forecast_factor=args.forecast_factor,
            )
            _print(manifest)
        if args.run:
            if not args.write and memory._manifest_path(args.output).is_file() is False:
                duration = memory.SEGMENT_DURATION if args.duration is None else float(args.duration)
                memory.prepare_campaign(
                    args.output, fermions=fermions, duration=duration, step_cap=args.step_cap,
                    cpu_budget_seconds=args.cpu_budget, forecast_factor=args.forecast_factor,
                )
            manifest = memory.run_campaign(args.output, cpu_budget_seconds=args.cpu_budget)
            if not args.write:
                _print(manifest)
        if args.check:
            report = memory.check_campaign(args.output)
            if not args.write and not args.run:
                _print(report)
            elif report["stepped"]:
                raise RuntimeError("check took a time step")
        return 0
    record = memory.dependency_closure()
    text_path = args.output
    _print(record)
    if text_path is not None:
        path = memory.assert_output_outside_repository(text_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(memory.closure_bytes(record))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, PermissionError, FileExistsError, FileNotFoundError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
