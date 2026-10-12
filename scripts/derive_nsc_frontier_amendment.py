"""Frontier amendment: the t>20 cell is closed by the recorded bytes.

The round-63 frontier assessment listed 'wall morphology beyond t=20'
as an open cell whose blocker was the lock-listed continuation. Rounds
73-77 answered the cell from the relaxation and march payloads, so the
assessment is amended: the cell is closed, and the open production list
shrinks accordingly. The amendment quotes the closing records and does
not touch the sealed originals.

    python scripts/lab.py scripts/derive_nsc_frontier_amendment.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

ASSESSMENT = LAB / "results" / "development" / "nsc-frontier-assessment-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-frontier-amendment-v1.json"
SCHEMA = "NSC-FRONTIER-AMENDMENT-v1"


def main():
    if not ASSESSMENT.is_file():
        raise SystemExit("missing frontier assessment record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    assessment = json.loads(ASSESSMENT.read_text(encoding="utf-8"))
    closed = [
        row for row in assessment["open"]
        if row["cell"].startswith("wall morphology beyond t=20")
        or row["cell"].startswith("derivation of p_Q")
    ]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "amends": "nsc-frontier-assessment-v1.json",
        "closed_cells": closed,
        "closed_by_records": {
            "wall_morphology": [
                "nsc-relaxation-wall-v1.json",
                "nsc-wall-afterlife-v1.json",
                "nsc-march-balance-trajectory-v1.json",
                "nsc-wall-extinction-v1.json",
            ],
            "p_q_derivation": [
                "nsc-q-decomposition-v1.json",
                "nsc-q-analytic-envelope-v1.json",
            ],
        },
        "remaining_open_cells": [
            row for row in assessment["open"]
            if not row["cell"].startswith("wall morphology beyond t=20")
            and not row["cell"].startswith("derivation of p_Q")
        ],
        "reading": (
            "two round-63 cells are closed by later records: the wall's "
            "t>20 morphology (rounds 73-77) and the p_Q derivation (rounds "
            "67-68). The t=40 campaign is no longer required for either."
        ),
        "notes": [
            "the round-63 assessment itself remains creation-only and "
            "untouched; this amendment stands beside it",
        ],
        "input_sha256": {
            "lab/results/development/nsc-frontier-assessment-v1.json":
                chain.path_sha256(ASSESSMENT),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_frontier_amendment.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_frontier_amendment.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  closed cells:", [row["cell"] for row in record["closed_cells"]])
    print("  remaining open cells:", [row["cell"] for row in record["remaining_open_cells"]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())