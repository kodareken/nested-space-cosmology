#!/usr/bin/env python3
"""Construct once or compactly verify FGC-1-TDG10-QA1-PREF1.

Default mode validates an existing compact result and never runs live binding.
If the compact result is absent, default mode fails clearly. Explicit --live
is the only path that may import the independent binder, call build, and write
the compact result atomically under the canonical pretty-JSON convention.
"""

from __future__ import annotations

import argparse
import importlib
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

CONFIG_PATH = "configs/fgc/fgc-1-tdg10-qa1-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg10-qa1-pref1.json"
BINDER_MODULE = "recursive_horizons.fgc.evolution.tdg10_qa1_pref1_binder"
ABSENT_RESULT_MESSAGE = (
    "FGC-1-TDG10-QA1-PREF1 compact result is absent at "
    f"{RESULT_PATH}. Default mode validates a compact result only when that "
    "file exists; it never runs live binding. Pass --live only for the "
    "explicit one-time independent binder."
)
MISSING_BINDER_MESSAGE = (
    "FGC-1-TDG10-QA1-PREF1 independent binder module is not present "
    f"({BINDER_MODULE}). This scaffold never invents a compact result."
)


def load_binder():
    try:
        return importlib.import_module(BINDER_MODULE)
    except ImportError as exc:
        raise SystemExit(MISSING_BINDER_MESSAGE) from exc


def write_canonical_result(destination: Path, payload: bytes) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    if temporary.exists():
        raise SystemExit("stale QA1 PREF1 temporary result exists")
    try:
        with temporary.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--output", default=RESULT_PATH)
    arguments = parser.parse_args()
    output = ROOT / arguments.output
    if arguments.live:
        binder = load_binder()
        config = (ROOT / binder.CONFIG_PATH).read_bytes()
        result = binder.build_pref1_result(config, ROOT, live=True)
        write_canonical_result(output, binder.canonical_result(result))
        print(
            "FGC-1-TDG10-QA1-PREF1 live binder wrote the compact result; "
            "this is not RA1 authorization, a production comparator, GR-0 "
            "calibration, or physics"
        )
        return 0
    if not output.is_file():
        raise SystemExit(ABSENT_RESULT_MESSAGE)
    binder = load_binder()
    config = (ROOT / binder.CONFIG_PATH).read_bytes()
    binder.validate_compact_result(config, output.read_bytes())
    print(
        "FGC-1-TDG10-QA1-PREF1 compact result verified; this is not RA1 "
        "authorization, a production comparator, GR-0 calibration, or physics"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
