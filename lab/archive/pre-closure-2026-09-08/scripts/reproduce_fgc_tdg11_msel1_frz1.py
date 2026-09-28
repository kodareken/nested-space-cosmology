#!/usr/bin/env python3
"""Outcome-blind compact freeze reproduction for FGC-1-TDG11-MSEL1-FRZ1.

Default ``--check`` verifies the tracked compact bundle only. Emit modes print
deterministic config or compact-result bytes and do not write files, inspect
the raw store, launch a shadow, or authorize execution.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_msel1_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_msel1_contract as contract  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--check",
        action="store_true",
        help="validate the tracked compact freeze (default)",
    )
    modes.add_argument(
        "--emit-config",
        action="store_true",
        help="print deterministic config TOML from the live environment",
    )
    modes.add_argument(
        "--emit-result",
        action="store_true",
        help="print compact_freeze JSON for the existing tracked config",
    )
    arguments = parser.parse_args(argv)
    if arguments.emit_config:
        sys.stdout.buffer.write(authority.emit_config_bytes(ROOT))
        return 0
    if arguments.emit_result:
        sys.stdout.buffer.write(
            canonical_json_bytes(authority.emit_result_object(ROOT)) + b"\n"
        )
        return 0
    bundle = authority.validate_compact_bundle(ROOT)
    sys.stdout.write(
        json.dumps(
            {
                "artifact_id": contract.ARTIFACT_ID,
                "mode": "check",
                "config_sha256": bundle["config_sha256"],
                "diagnostic_executed": bundle["diagnostic_executed"],
                "selection_result_earned": bundle["selection_result_earned"],
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
