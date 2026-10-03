#!/usr/bin/env python3
"""Preview or run an explicit parent-response stage.

The default prints the contract and writes nothing. predict, measure and
check read an existing stage directory. prepare is not invented here: root
calls ``prepare(directory, binding, producer_commit=...)`` with a frozen
physical binding. A missing seal leaves measure open. This driver does not
launch the physical held campaign by itself.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

from recursive_horizons import nsc_discovery_parent_response as parent


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    for mode in ("predict", "measure", "check"):
        modes.add_argument("--" + mode, action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    selected = next((mode for mode in ("predict", "measure", "check") if getattr(args, mode)), None)
    if selected is None:
        print(json.dumps({
            "specification": parent.specification(),
            "dependencies": parent.dependencies(),
            "callbacks": {
                "owner1": "nsc_discovery_parent.prepare_parent",
                "owner2": "nsc_discovery_parent_episode.reconstruct_parent_pair",
                "called": False,
            },
        }, indent=2, sort_keys=True))
        return 0
    if args.output is None:
        raise SystemExit(selected + " requires --output")
    report = getattr(parent, selected)(args.output.expanduser().resolve())
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
