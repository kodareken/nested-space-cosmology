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
    for mode in ("prepare", "predict", "measure", "run", "check", "prepare-curvature", "predict-curvature", "measure-curvature", "preflight", "confirm-curvature"):
        modes.add_argument("--"+mode, action="store_true")
    parser.add_argument("--output")
    parser.add_argument("--nonlinear", action="store_true")
    parser.add_argument("--reference", help="immutable nonlinear-v2 forecast and measurement for one numerical confirmation")
    parser.add_argument("--baseline", help="read-only nominal linear preparation/prediction; held-out data are not training")
    parser.add_argument("--redo-baseline", action="store_true")
    parser.add_argument("--nf", type=int, default=128)
    parser.add_argument("--duration", type=float, default=0.01)
    parser.add_argument("--step-cap", type=float, default=width.episode.MATCHED_STEP_CAP)
    parser.add_argument("--production", action="store_true")
    parser.add_argument("--cpu-budget", type=float, default=width.CPU_CAP)
    args = parser.parse_args(argv)
    try:
        if args.preflight:
            if not args.baseline:
                parser.error("--preflight requires --baseline")
            result = width.nonlinear_preflight(args.baseline)
            print(width.NONLINEAR_SCHEMA, "readonly_preflight", result["fresh_batch_forecast_cpu_indicator"], flush=True)
            return 0
        if not args.output:
            parser.error("--output is required for an execution stage")
        if args.confirm_curvature:
            if not args.reference:
                parser.error("--confirm-curvature requires --reference")
            result = width.confirm_curvature(args.output, args.reference, production=args.production,
                                             cpu_budget=min(args.cpu_budget, width.CONFIRMATION_CPU_CAP))
            print(width.CONFIRMATION_SCHEMA, result["status"], "cpu", result["aggregate_cpu_seconds"], flush=True)
            return 0
        nonlinear = args.nonlinear or args.prepare_curvature or args.predict_curvature or args.measure_curvature
        if nonlinear:
            if args.prepare_curvature or args.prepare or args.run:
                result = width.prepare_curvature(args.output, baseline=args.baseline, nf=args.nf,
                    duration=args.duration, step_cap=args.step_cap, production=args.production,
                    cpu_budget=args.cpu_budget, redo_baseline=args.redo_baseline)
            elif args.predict_curvature or args.predict:
                result = width.predict_curvature(args.output)
            elif args.measure_curvature or args.measure:
                result = width.measure_curvature(args.output)
            else:
                result = width.check_curvature(args.output)
            if args.run and result["status"] == "prepared":
                result = width.predict_curvature(args.output)
                if result["status"] == "predicted":
                    result = width.measure_curvature(args.output)
            print(width.NONLINEAR_SCHEMA, result.get("status", "checked"), "cpu", result.get("aggregate_cpu_seconds"), flush=True)
            return 0
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
