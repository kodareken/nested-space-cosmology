#!/usr/bin/env python3
"""Verify compact PREF1 by default; explicitly bind raw evidence only on request."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons import evidence_io as io  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_msel1_pref1_binder as binder  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--check",
        action="store_true",
        help="compact, raw/store/shadow/Git-blind check (default)",
    )
    action.add_argument(
        "--emit-config",
        action="store_true",
        help="print config for the fixed published raw pins",
    )
    action.add_argument(
        "--bind", action="store_true", help="explicit independent live reconstruction"
    )
    parser.add_argument(
        "--write-result",
        action="store_true",
        help="with --bind, exclusively create the declared compact result file",
    )
    args = parser.parse_args(argv)
    if args.write_result and not args.bind:
        parser.error("--write-result is available only with --bind")
    if args.emit_config:
        sys.stdout.buffer.write(binder.emit_config_bytes(ROOT))
        return 0
    if args.bind:
        if args.write_result and os.path.lexists(ROOT / binder.RESULT_PATH):
            raise binder.TDG11PREF1Error(
                "compact destination already exists; nothing was replayed or replaced"
            )
        config_raw = io.read_regular_file(ROOT, binder.CONFIG_PATH)
        result = binder.bind_raw_result(
            config_raw,
            ROOT,
            progress=lambda event: print(
                io.canonical_json_bytes(event).decode(), file=sys.stderr, flush=True
            ),
        )
        result_raw = io.canonical_json_bytes(result) + b"\n"
        if args.write_result:
            binder.require_live_config(config_raw, ROOT)
            binder.validate_compact_result(config_raw, result_raw, ROOT)
            io.publish_exclusive_file(ROOT, binder.RESULT_PATH, result_raw)
            # A postpublication drift is an error, never a successful binding.
            # The exclusive output remains forensic evidence and is not removed.
            binder.require_live_config(config_raw, ROOT)
            binder.verify_compact(ROOT)
        else:
            sys.stdout.buffer.write(result_raw)
            return 0
    else:
        result = binder.verify_compact(ROOT)
    print(
        json.dumps(
            {
                "artifact_id": binder.ARTIFACT_ID,
                "mode": "independent_live_bind" if args.bind else "compact_check",
                "classification": result["classification"],
                "selected_candidate": result["selected_candidate"],
                "replay_accounting": {
                    key: value
                    for key, value in result["independent_replay"].items()
                    if key != "widths"
                },
                "conclusion": result["conclusion"],
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
