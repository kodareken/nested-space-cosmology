#!/usr/bin/env python3
"""Create, run, or check the prediction comparison.

``--create`` writes the declaration only. ``--run`` writes
``prediction-run.json`` and ``prediction-run.npz`` for the selected
``nf`` values. The default duration is ``0.3``; ``1`` is optional.
``--check`` reloads a declaration, or replays observables from a run
record without stepping. One process, FFT evolution, no acceptance gate.
"""
from __future__ import annotations

import argparse
import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

from recursive_horizons.nsc_discovery_prediction import (
    check_successor,
    create_successor,
    run_saved_comparison,
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--create", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--output", required=True,
                        help="tmp directory, or lab/results/development/nsc-discovery-prediction-v1")
    parser.add_argument("--duration", type=float, default=0.3,
                        help="coordinate duration: short test <=0.01, first target 0.3, optional 1")
    parser.add_argument("--nf", type=int, action="append", choices=(128, 256),
                        help="selected resolution; repeat for both. Default is 128 and 256")
    parser.add_argument("--step-cap", type=float, default=0.0005)
    args = parser.parse_args(argv)
    if args.create:
        manifest = create_successor(args.output)
        print(manifest["status"], "frames", len(manifest["frames"]),
              "evolved", manifest["evolved"], "later", manifest["later_outputs_written"], flush=True)
        return 0
    if args.run:
        nfs = (128, 256) if not args.nf else tuple(args.nf)
        report = run_saved_comparison(nfs, args.duration, args.output, step_cap=args.step_cap)
        print(report["status"], "cases", len(report["cases"]), "evolved", report["evolved"],
              "forecast_cpu", report["cpu_forecast"]["forecast_cpu_seconds"], flush=True)
        return 0
    report = check_successor(args.output)
    print(report["schema"], "ok", report["ok"], "replayed", report.get("replayed", False),
          "evolved", report["evolved"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
