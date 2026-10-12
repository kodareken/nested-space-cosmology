"""Timing correction: the poll count was poll-age, not march-age.

The round-100 annals quoted the historical poll as proof that the
uncontended march reached station-20 in about ten minutes. Re-examining
the poll loop: iteration 10 means the poll's own start-to-landing watch
was about ten minutes, while the march it watched had been launched
before the polling began. The march's total duration at that landing is
therefore unrecorded, and the ten-minute figure is retracted as a
march-duration statement. The t30 run's pace is judged by its own
monotone CPU accrual, which continues to hold.

    python scripts/lab.py scripts/derive_nsc_t30_timing_correction.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-timing-correction-v1.json"
SCHEMA = "NSC-T30-TIMING-CORRECTION-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "corrects": ["nsc-centenary-v1.json", "nsc-t30-launch-log-v1.json"],
        "retracted": (
            "the statement that the historical march reached station-20 "
            "in about ten minutes"
        ),
        "evidence": (
            "the polling loop counts its own iterations: each is a check "
            "then a sixty-second sleep, and iteration ten is about ten "
            "minutes of polling. The poll began after the march had "
            "already been launched, so the march's age at landing is not "
            "measured by that count at all."
        ),
        "correct_reading": (
            "the prior march's total duration is unrecorded; the current "
            "t30 run is judged by its own monotone CPU accrual, which "
            "holds, and by the absence of any failure trail"
        ),
        "notes": [
            "discipline note: a benchmark record repeated a lower bound as "
            "a measurement; re-reading the poll revealed the counting "
            "semantics and this correction is the fix",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_timing_correction.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_timing_correction.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  retracted:", record["retracted"][:60] + "...")
    print("  correct reading:", record["correct_reading"][:70] + "...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())