#!/usr/bin/env python3
"""Compact freeze reproduction for FGC-1-DEF1-STAB1-FRZ1.

Default ``--check`` validates tracked compact bytes only.  It is Git/raw/
runs/source/trajectory blind and never executes a runner.  ``--emit-*`` and
``--write-result`` are explicit live reconstruction from DEF1 owner modules.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.def1_stab1_frz1_certificate import (  # noqa: E402
    ARTIFACT_ID,
    CONFIG_PATH,
    RESULT_PATH,
    verify_compact,
)


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
        help="print deterministic config TOML from live DEF1 owners",
    )
    modes.add_argument(
        "--emit-result",
        action="store_true",
        help="print compact freeze JSON from live DEF1 owners",
    )
    modes.add_argument(
        "--write-result",
        action="store_true",
        help="explicit initial construction; refuses an existing result",
    )
    arguments = parser.parse_args(argv)
    if arguments.emit_config or arguments.emit_result or arguments.write_result:
        from recursive_horizons.evidence_io import publish_exclusive_file
        from recursive_horizons.fgc.def1_stab1_frz1_certificate import (
            compose_canonical_artifacts,
        )

        if arguments.emit_config:
            config_raw, _result_raw = compose_canonical_artifacts(ROOT)
            sys.stdout.buffer.write(config_raw)
            return 0
        if arguments.emit_result:
            _config_raw, result_raw = compose_canonical_artifacts(ROOT)
            sys.stdout.buffer.write(result_raw)
            return 0
        if os.path.lexists(ROOT / RESULT_PATH) or os.path.lexists(ROOT / CONFIG_PATH):
            parser.error("compact destination already exists; nothing was replaced")
        config_raw, result_raw = compose_canonical_artifacts(ROOT)
        publish_exclusive_file(ROOT, CONFIG_PATH, config_raw)
        publish_exclusive_file(ROOT, RESULT_PATH, result_raw)
        result = verify_compact(ROOT)
        sys.stdout.write(
            json.dumps(
                {
                    "artifact_id": ARTIFACT_ID,
                    "mode": "write_result",
                    "classification": result["classification"],
                    "def1_error_map_passed": result["def1_error_map_passed"],
                },
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n"
        )
        return 0
    result = verify_compact(ROOT)
    sys.stdout.write(
        json.dumps(
            {
                "artifact_id": ARTIFACT_ID,
                "mode": "check",
                "classification": result["classification"],
                "conversion_instrument_contract_complete": result["claims"][
                    "conversion_instrument_contract_complete"
                ],
                "def1_error_map_passed": result["def1_error_map_passed"],
                "licenses_separate_independent_pref1_binder": result[
                    "licenses_separate_independent_pref1_binder"
                ],
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
