#!/usr/bin/env python3
"""Deterministic PRO20-EV1 freeze-config emitter.

``--emit-config`` prints config bytes and does not write a file, authorize a
run, or construct a source. ``--check`` compares a later tracked config to a
fresh emission and still does not authorize or execute.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = str(ROOT / "src")
if _SOURCE_ROOT not in sys.path:
    sys.path.insert(0, _SOURCE_ROOT)

from recursive_horizons.fgc.evolution import (  # noqa: E402
    pro20_ev1_authority as authority,
)


ARTIFACT_ID = "FGC-1-PRO20-EV1-FRZ1"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--emit-config",
        action="store_true",
        help="print deterministic freeze-config TOML (default)",
    )
    modes.add_argument(
        "--check",
        action="store_true",
        help="compare the tracked freeze config to a fresh emission",
    )
    arguments = parser.parse_args(argv)
    if arguments.check:
        record = authority.check_tracked_config(ROOT)
        sys.stdout.write(
            json.dumps(record, sort_keys=True, indent=2, ensure_ascii=True) + "\n"
        )
        return 0
    sys.stdout.buffer.write(authority.emit_config_bytes(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
