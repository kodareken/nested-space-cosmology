#!/usr/bin/env python3
"""Verify SOL1-PREF1 compact bytes or explicitly reconstruct them."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.sgb1_ctl1_sol1_pref1_binder import (  # noqa: E402
    ARTIFACT_ID,
    CONFIG_PATH,
    RESULT_PATH,
    verify_compact,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--emit-config", action="store_true")
    modes.add_argument("--emit-result", action="store_true")
    modes.add_argument("--write-result", action="store_true")
    args = parser.parse_args(argv)
    if args.emit_config or args.emit_result or args.write_result:
        from recursive_horizons.evidence_io import publish_exclusive_file
        from recursive_horizons.fgc.sgb1_ctl1_sol1_pref1_binder import (
            compose_canonical_artifacts,
        )

        config_raw, result_raw = compose_canonical_artifacts(ROOT)
        if args.emit_config:
            sys.stdout.buffer.write(config_raw)
            return 0
        if args.emit_result:
            sys.stdout.buffer.write(result_raw)
            return 0
        if os.path.lexists(ROOT / CONFIG_PATH) or os.path.lexists(ROOT / RESULT_PATH):
            parser.error("compact destination already exists; nothing was replaced")
        publish_exclusive_file(ROOT, CONFIG_PATH, config_raw)
        publish_exclusive_file(ROOT, RESULT_PATH, result_raw)
        result = verify_compact(ROOT)
        print(json.dumps({
            "artifact_id": ARTIFACT_ID,
            "classification": result["classification"],
            "mode": "write_result",
        }, sort_keys=True, separators=(",", ":")))
        return 0
    result = verify_compact(ROOT)
    print(json.dumps({
        "artifact_id": ARTIFACT_ID,
        "classification": result["classification"],
        "mode": "check",
        "nominal_classification": result["nominal"]["classification"],
        "nominal_obstruction": result["nominal"]["obstruction"],
    }, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
