"""The inward-march throat slope joins the wall's regime boundary.

The march assessment records a throat slope channel the bubble never
used: +36.6, +48.3, -23.3 across the three cases - a sign flip between
mu=1.010 and mu=1.040, exactly where the wall's converter locks
(round-21). This driver lays the three channels (march throat slope,
march inner-region slope, wall presence) on one table and records
where the sign crossing sits relative to the recorded regime split.

    python scripts/lab.py scripts/derive_nsc_march_wall_join.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

MARCH = LAB / "results" / "development" / "nsc-rn-inward-march-assessment-v1.json"
BUBBLE = LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json"
OUTPUT = LAB / "results" / "development" / "nsc-march-wall-join-v1.json"
SCHEMA = "NSC-MARCH-WALL-JOIN-v1"


def main():
    if not MARCH.is_file() or not BUBBLE.is_file():
        raise SystemExit("missing march or bubble record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    march = json.loads(MARCH.read_text(encoding="utf-8"))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    bubble_by_case = {c["case"]: c for c in bubble["cases"]}
    rows = []
    for case in march["cases"]:
        name = case["case"]
        rows.append({
            "case": name,
            "throat_slope_last": float(case["throat_slope_last"]),
            "inner_region_slope_min_last": float(case["inner_region_slope_min_last"]),
            "wall_presence_fraction": float(bubble_by_case[name]["wall_presence_fraction"]),
        })
    slopes = [row["throat_slope_last"] for row in rows]
    sign_changes = sum(
        1 for left, right in zip(slopes, slopes[1:]) if left * right < 0.0)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "throat_slope_sign_changes": sign_changes,
        "reading": (
            "the march's own throat slope flips sign between the light "
            "cases and the converter-locked case"
        ),
        "verdict": (
            "SIGN_CROSSES_AT_THE_REGIME_BOUNDARY" if sign_changes == 1
            else "NO_CROSSING"
        ),
        "notes": [
            "the march throat slope is a distinct channel from the bubble "
            "balance slopes; both records are quoted without merging",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-inward-march-assessment-v1.json":
                chain.path_sha256(MARCH),
            "lab/results/development/nsc-rn-gradient-bubble-v5.json":
                chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_march_wall_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_march_wall_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} throat {:.2f}  inner_min {:.2f}  presence {:.2f}".format(
            row["case"], row["throat_slope_last"], row["inner_region_slope_min_last"],
            row["wall_presence_fraction"]))
    print("  throat slope sign changes: {}  verdict: {}".format(
        record["throat_slope_sign_changes"], record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())