#!/usr/bin/env python3
"""Preview the bounded exact BR control; explicitly create or check a new record."""
from __future__ import annotations

import argparse
import json

from recursive_horizons import nsc_discovery_vacuum_control as control


def digest(report):
    return {
        "schema": report["schema"], "cpu_seconds": report["cpu_seconds"],
        "parameters": report["parameters"], "scope": report["scope"],
        "exact_T8": {nf: rows[-1] for nf, rows in report["exact"].items()},
        "numerical_T3": {nf: case["stations"][-1] for nf, case in report["numerical"].items()},
        "source_hashes_unchanged": report["source_hashes_unchanged"],
        "input_hashes_unchanged": report["input_hashes_unchanged"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", nargs="?", const=str(control.OUTPUT), metavar="DIRECTORY",
                      help="create immutable JSON+NPZ below the owned vacuum-control-v1 directory")
    mode.add_argument("--check", nargs="?", const=str(control.OUTPUT), metavar="DIRECTORY",
                      help="read-only check; does not run rates or RK4")
    args = parser.parse_args(argv)
    if args.check is not None:
        print(json.dumps(control.check_record(args.check), indent=2, sort_keys=True))
        return 0
    if args.write is not None:
        destination = control.output_directory(args.write)
        if destination.exists():
            raise FileExistsError("refusing to overwrite " + str(destination))
    report, arrays = control.run_control()
    if args.write is not None:
        control.write_record(report, arrays, destination)
        control.check_record(destination)
    print(json.dumps(digest(report), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
