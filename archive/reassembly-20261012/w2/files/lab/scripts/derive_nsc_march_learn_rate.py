"""March learn rates: the interval shrink as a per-case decay factor.

The interval lists are short (3 points after the first); this driver
fits the geometric decay I(j) = I2 * s^(j-1) on the recorded intervals
per case and records the learn-rate factors, the compact companion to
the round-106 partial join. No further claim than the data carries.

    python scripts/lab.py scripts/derive_nsc_march_learn_rate.py
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

DURATION = LAB / "results" / "development" / "nsc-march-duration-table-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-march-learn-rate-v1.json"
SCHEMA = "NSC-MARCH-LEARN-RATE-v1"


def main():
    if not DURATION.is_file():
        raise SystemExit("missing duration table")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    duration = json.loads(DURATION.read_text(encoding="utf-8"))
    rows = []
    for row in duration["rows"]:
        intervals = row["interval_minutes"][1:]
        factors = [intervals[i + 1] / intervals[i] for i in range(len(intervals) - 1)]
        log_mean = sum(math.log(f) for f in factors) / len(factors)
        rows.append({
            "case": row["case"],
            "after_first_intervals": intervals,
            "pairwise_factors": factors,
            "geometric_factor_s": math.exp(log_mean),
        })
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "reading": (
            "the per-case exponential-in-index decay factor of the march's "
            "station intervals; the heaviest case decays fastest"
        ),
        "notes": [
            "three intervals per case; the factor is a fit, not a law",
        ],
        "input_sha256": {
            "lab/results/development/nsc-march-duration-table-v1.json":
                chain.path_sha256(DURATION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_march_learn_rate.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_march_learn_rate.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} intervals {} factors {} s {:.4f}".format(
            row["case"], row["after_first_intervals"], row["pairwise_factors"],
            row["geometric_factor_s"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())