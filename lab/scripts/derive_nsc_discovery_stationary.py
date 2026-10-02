#!/usr/bin/env python3
"""Preview a radial seed or execute one bounded NF32 stationary query."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

from recursive_horizons import nsc_discovery_stationary as stationary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve", action="store_true", help="explicit one bounded stationary-critical attempt")
    parser.add_argument("--write", type=Path, nargs="?", const=stationary.OUTPUT,
                        help="create new JSON+NPZ at the v1 stem or a /tmp stem")
    parser.add_argument("--check", action="store_true", help="read-only hashes and saved-array consistency")
    parser.add_argument("--record", type=Path, default=stationary.OUTPUT, help="existing record stem for --check")
    parser.add_argument("--cpu-limit", type=float, default=stationary.CPU_LIMIT, help="CPU seconds, at most 120")
    parser.add_argument("--max-nfev", type=int, default=1000)
    parser.add_argument("--amplitude", type=float, default=.35, help="nonconstant cosine radial seed amplitude")
    parser.add_argument("--producer-commit", help="optional commit, accepted only for exact current producer bytes")
    args = parser.parse_args(argv)
    if args.check:
        if args.solve or args.write is not None or args.producer_commit:
            parser.error("--check cannot solve, write, or change a producer pin")
        report = stationary.check_record(args.record)
    else:
        if args.write is not None:
            stationary.validate_new_output(args.write)
        report, arrays = stationary.run_stationary(solve=args.solve, cpu_limit=args.cpu_limit,
                                                   max_nfev=args.max_nfev, amplitude=args.amplitude,
                                                   producer_commit=args.producer_commit)
        if args.write is not None:
            report = stationary.write_record(report, arrays, args.write)
    sys.stdout.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    # An unsatisfied finite query is a measured scientific outcome, not a CLI
    # failure or a nonexistence claim. Errors alone return a nonzero code.
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
