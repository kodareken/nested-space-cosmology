#!/usr/bin/env python3
"""Preview exact saved-field free transport; explicitly create/check its new record."""
from __future__ import annotations

import argparse
import json
import math
import resource
import time

from recursive_horizons import nsc_discovery_free_return as free


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true", help="create the fixed free-return-v1 JSON owner once")
    modes.add_argument("--check", action="store_true", help="authenticate and replay exact matrix transport, no ODE")
    parser.add_argument("--controls", action=argparse.BooleanOptionalAction, default=True,
                        help="include saved nf128/timestep/frozen controls (default selected)")
    args = parser.parse_args(argv)
    if args.write and free.OUTPUT.exists():
        raise FileExistsError("refusing to overwrite "+str(free.OUTPUT))
    soft, hard = resource.getrlimit(resource.RLIMIT_CPU)
    limit = math.ceil(time.process_time()+free.CPU_BUDGET)
    if hard != resource.RLIM_INFINITY:
        limit = min(limit, hard)
    if soft == resource.RLIM_INFINITY or soft > limit:
        resource.setrlimit(resource.RLIMIT_CPU, (limit, hard))
    if args.check:
        print(json.dumps(free.check_record(), indent=2))
        return 0
    report = free.measure(include_controls=args.controls)
    if args.write:
        free.write_record(report)
    print(json.dumps({"schema": report["schema"], "cpu_seconds": report["cpu_seconds"],
                      "interpretation": report["interpretation"], "cases": report["cases"]},
                     indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
