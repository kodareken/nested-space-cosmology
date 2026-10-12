"""The composed covariance constant against the infall clock.

The identity-implied composed values (round 60) sit near 15.1-15.9;
the infall clock T_total (round 25) sits near 16.3-16.5 and is
mu-blind. This driver lays the two recorded series case by case and
records the ratios. The proximity is observation-level and flagged as
scale-proximity, not identity.

    python scripts/lab.py scripts/derive_nsc_composed_vs_clock.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

LADDER_TEST = LAB / "results" / "development" / "nsc-ladder-covariance-test-v1.json"
TIMESCALE = LAB / "results" / "development" / "nsc-wall-breathing-timescale-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-composed-vs-clock-v1.json"
SCHEMA = "NSC-COMPOSED-VS-CLOCK-v1"


def main():
    if not LADDER_TEST.is_file() or not TIMESCALE.is_file():
        raise SystemExit("missing composed or timescale record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    ladder_test = json.loads(LADDER_TEST.read_text(encoding="utf-8"))
    timescale = json.loads(TIMESCALE.read_text(encoding="utf-8"))
    t_by_case = {
        name: entry["T_total"] for name, entry in timescale["cases"].items()
    }
    rows = []
    for entry in ladder_test["rows"]:
        composed = entry["composed"]
        t_total = t_by_case[entry["case"]]
        rows.append({
            "case": entry["case"],
            "composed": composed,
            "T_total": t_total,
            "ratio_composed_over_clock": composed / t_total,
        })
    ratios = [row["ratio_composed_over_clock"] for row in rows]
    mean_ratio = sum(ratios) / len(ratios)
    spread = (sum((v - mean_ratio) ** 2 for v in ratios) / len(ratios)) ** 0.5 / abs(mean_ratio)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "mean_ratio": mean_ratio,
        "ratio_spread": spread,
        "status": (
            "observation-level scale proximity; flagged, not claimed"
        ),
        "notes": [
            "T_total is the mu-blind cap clock of round 25; the composed "
            "values are the ladder-implied covariance constants",
        ],
        "input_sha256": {
            "lab/results/development/nsc-ladder-covariance-test-v1.json":
                chain.path_sha256(LADDER_TEST),
            "lab/results/development/nsc-wall-breathing-timescale-v1.json":
                chain.path_sha256(TIMESCALE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_composed_vs_clock.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_composed_vs_clock.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} composed {:.3f}  T_total {:.3f}  ratio {:.3f}".format(
            row["case"], row["composed"], row["T_total"], row["ratio_composed_over_clock"]))
    print("  mean ratio {:.3f}  spread {:.2%}  status: {}".format(
        record["mean_ratio"], record["ratio_spread"], record["status"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())