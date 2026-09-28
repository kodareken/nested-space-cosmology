#!/usr/bin/env python3
"""Reproduce the first finite curved Dirac boundary maps without overwriting them.

This runner compares every scientific numeric field, exact flag, convention
entry, and source hash of results/nsc-3-boundary-response.json. It does not
launch a campaign, modify historical outputs, or select a physical scale.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from recursive_horizons.nsc_boundary import (  # noqa: E402
    DEFAULT_RESULT,
    build_record,
    compare,
    jsonable,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check and args.output:
        parser.error("choose --output or --check")
    if args.output is not None and args.output.exists():
        raise FileExistsError(
            f"refusing to overwrite existing result: {args.output}"
        )
    record = jsonable(build_record())
    if args.check:
        expected = json.loads(DEFAULT_RESULT.read_bytes())
        compare(expected, record)
        print(
            "Boundary-response calculation reproduced; "
            "every recorded numeric, exact, convention, and hash field checked."
        )
    elif args.output:
        destination = args.output
        payload = json.dumps(
            record, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
        with destination.open("xb") as stream:
            stream.write(payload)
        print(
            json.dumps(
                {
                    "output": str(destination),
                    "schema": record["schema"],
                    "classification": record["classification"],
                    "gate": record["gate"],
                    "nonclaims": record["nonclaims"],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        print(json.dumps(record, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
