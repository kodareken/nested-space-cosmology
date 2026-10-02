#!/usr/bin/env python3
"""Preflight or explicitly construct general-Q chi=0 Cauchy data; never evolve."""
import argparse
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[name]="1"

from recursive_horizons import nsc_discovery_dynamic_preparation as preparation


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare",action="store_true",help="explicit bounded initial-data construction")
    parser.add_argument("--write",type=Path,help="create a new JSON+NPZ checkpoint after preparation")
    parser.add_argument("--check",type=Path,help="read-only authentication/replay of a saved checkpoint")
    parser.add_argument("--cpu-limit",type=float,default=30.)
    args=parser.parse_args(argv)
    if args.check is not None:
        if args.prepare or args.write is not None:
            parser.error("check cannot prepare or write")
        report=preparation.check_record(args.check)
    else:
        if args.write is not None and not args.prepare:
            parser.error("write requires explicit prepare")
        if args.write is not None:
            if any(p.exists() for p in preparation.output_paths(args.write)):
                raise FileExistsError("checkpoint already exists")
        report,arrays=preparation.build_record(execute=args.prepare,cpu_limit=args.cpu_limit)
        if args.write is not None:
            report=preparation.write_record(report,arrays,args.write)
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
