#!/usr/bin/env python3
"""Compute finite saved-slice geometric vertices; writing is explicit/exclusive."""
import argparse
import json
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

from recursive_horizons import nsc_discovery_geometric_vertices as vertices


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument("--write", action="store_true")
    actions.add_argument("--check", action="store_true")
    parser.add_argument("--output", default=str(vertices.DEFAULT_OUTPUT))
    parser.add_argument("--science-commit", help="full frozen Git commit, required for an exclusive write")
    args = parser.parse_args(argv)
    if args.write and args.science_commit is None:
        parser.error("--write requires --science-commit after the source checkpoint")
    if args.check:
        record = vertices.check_record(args.output)
    elif args.write:
        record = vertices.write_record(args.output, args.science_commit)
    else:
        record = vertices.calculate()
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
