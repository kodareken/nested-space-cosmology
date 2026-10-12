"""The kappa-ladder approach law: deficits as exp(-(a n + b n^2)).

Round-9 recorded the reflection deficits on the kappa ladder n=1..4 to
subtlety: 0.332, 0.0433, 0.00222, 4.0e-5. This driver pre-registers the
two-parameter-at-most reading: if the approach to the reflection ceiling
is Gaussian in n, then ln d_n = -(a n + b n^2), the pair (a, b) fixed by
n=1,2 ONLY, and n=3,4 are predictions with no freedom left. The
prefactor is not fitted: the ansatz fixes it to 1 once (a, b) are fixed
from the first two deficits.

    python scripts/lab.py scripts/derive_nsc_charged_ladder_deficit_law.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

CLIFF = LAB / "results" / "development" / "nsc-rn-charged-twist-cliff-v3.json"
OUTPUT = LAB / "results" / "development" / "nsc-charged-ladder-deficit-law-v1.json"
SCHEMA = "NSC-CHARGED-LADDER-DEFICIT-LAW-v1"
DEFICITS = (0.332, 0.0433, 0.00222, 4.0e-5)


def main():
    if not CLIFF.is_file():
        raise SystemExit("missing cliff record: " + str(CLIFF))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    ln1 = math.log(DEFICITS[0])
    ln2 = math.log(DEFICITS[1])
    # a + b = -ln1 ; 2a + 4b = -ln2   ->   b = ((-ln2) - 2(-ln1))/2
    b = ((-ln2) - 2.0 * (-ln1)) / 2.0
    a = -ln1 - b
    rows = []
    for index, deficit in enumerate(DEFICITS, start=1):
        predicted = math.exp(-(a * index + b * index * index))
        rows.append({
            "n": index,
            "recorded_deficit": deficit,
            "predicted_from_n1_n2" if index in (3, 4) else "input_deficit": predicted,
            "relative_error": (predicted - deficit) / deficit,
        })
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "ansatz": "ln d_n = -(a n + b n^2); (a,b) fixed by n=1,2 only; prefactor 1",
        "a": a,
        "b": b,
        "rows": rows,
        "prediction_errors": {
            "n3_relative": rows[2]["relative_error"],
            "n4_relative": rows[3]["relative_error"],
        },
        "verdict": (
            "GAUSSIAN_APPROACH_HOLDS"
            if max(abs(rows[2]["relative_error"]), abs(rows[3]["relative_error"])) < 0.15
            else "GAUSSIAN_APPROACH_BROKEN"
        ),
        "reading": (
            "d_n ~ exp(-(0.64 n + 0.47 n^2)): the ladder's approach to the "
            "reflection ceiling is Gaussian in the rung index"
        ),
        "notes": [
            "deficits are the round-9 recorded R = 1 - d values; no "
            "coefficient beyond (a, b) enters and both are fixed by n=1,2",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-charged-twist-cliff-v3.json":
                chain.path_sha256(CLIFF),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_charged_ladder_deficit_law.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_charged_ladder_deficit_law.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  a = {:.4f}  b = {:.4f}".format(a, b))
    for row in rows:
        predicted_key = ("predicted_from_n1_n2" if row["n"] in (3, 4) else "input_deficit")
        print("  n={} recorded {:.3g} predicted {:.3g}  rel err {:+.2%}".format(
            row["n"], row["recorded_deficit"], row[predicted_key], row["relative_error"]))
    print("  verdict:", record["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())