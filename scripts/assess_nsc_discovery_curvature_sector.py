#!/usr/bin/env python3
"""Read-only symbolic flat-sector assessment, with explicit immutable creation."""
import argparse
import json
from pathlib import Path

from recursive_horizons import nsc_discovery_curvature_sector as sector


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer-commit",help="frozen source revision required for record creation")
    parser.add_argument("--write",type=Path,nargs="?",const=sector.OUTPUT,help="explicit creation-only record")
    parser.add_argument("--check",type=Path,nargs="?",const=sector.OUTPUT,help="read-only source/input/symbolic replay")
    args=parser.parse_args(argv)
    if args.check is not None:
        if args.write is not None or args.producer_commit:
            parser.error("check cannot create or change a producer pin")
        result=sector.check_record(args.check)
    else:
        if args.write is not None and not args.producer_commit:
            parser.error("creation requires a frozen producing commit")
        if args.write is not None and args.write.exists():
            raise FileExistsError("record already exists")
        result=sector.assess(producer_commit=args.producer_commit)
        if args.write is not None:
            result=sector.write_record(result,args.write)
    print(json.dumps(result,indent=2,allow_nan=False))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
