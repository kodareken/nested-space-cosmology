"""Closure of the two-ends join: both Q and rise as gap powers.

Round-51 recorded corr(Q, rise) = -0.93 (n=3). The natural closure
reading is covariance: both series are powers of the horizon gap. This
driver measures the gap exponent of Q from the recorded Q values
(p_Q, post-hoc flagged), uses the recorded rise power 0.5, and forms
the composed quantity Q * rise^(p_Q / p_rise). Its spread against the
raw Q*rise spread is the measurement: a strong collapse means the
correlation is absorbed by the two gap powers.

    python scripts/lab.py scripts/derive_nsc_two_ends_covariance.py
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
from recursive_horizons import nsc_rn_source as source  # noqa: E402

JOIN = LAB / "results" / "development" / "nsc-two-ends-join-v1.json"
BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-two-ends-covariance-v1.json"
SCHEMA = "NSC-TWO-ENDS-COVARIANCE-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}
RISE_POWER = 0.5


def spread(values):
    mean = sum(values) / len(values)
    return (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5 / abs(mean)


def main():
    if not JOIN.is_file() or not BREATHING.is_file():
        raise SystemExit("missing join or breathing record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    join = json.loads(JOIN.read_text(encoding="utf-8"))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    rows = []
    for entry in join["table"]:
        mass = MASSES[entry["case"]]
        model = source.action_reference(mass)
        gap = float(model.horizons["r_plus"]) - float(model.horizons["r_minus"])
        rows.append({
            "case": entry["case"],
            "gap": gap,
            "q": entry["room_quality_ratio"],
            "rise": breathing["cases"][entry["case"]]["mean_rise_rate"],
        })
    largest = max(rows, key=lambda row: row["gap"])
    smallest = min(rows, key=lambda row: row["gap"])
    p_q = math.log(largest["q"] / smallest["q"]) / math.log(largest["gap"] / smallest["gap"])
    composed_exponent = p_q / RISE_POWER
    raw_products = [row["q"] * row["rise"] for row in rows]
    composed = [row["q"] * (row["rise"] ** composed_exponent) for row in rows]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "q_gap_power_post_hoc": p_q,
        "rise_gap_power_recorded": RISE_POWER,
        "composed_exponent": composed_exponent,
        "rows": rows,
        "raw_product_spread": spread(raw_products),
        "composed_spread": spread(composed),
        "collapse_factor": spread(raw_products) / spread(composed),
        "verdict": (
            "COVARIANCE_CLOSED_BY_GAP"
            if spread(composed) < 0.06 else "COVARIANCE_NOT_CLOSED"
        ),
        "reading": (
            "the composed law Q * rise^0.234 absorbs the round-51 "
            "correlation when both series are gap powers; the Q exponent "
            "is post-hoc flagged as measured from these three points"
        ),
        "notes": [
            "the collapse factor is raw spread over composed spread",
            "no provenance for p_Q is claimed beyond the measurement",
        ],
        "input_sha256": {
            "lab/results/development/nsc-two-ends-join-v1.json": chain.path_sha256(JOIN),
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_two_ends_covariance.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_two_ends_covariance.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  p_Q = {:.4f}  composed exponent = {:.4f}".format(p_q, composed_exponent))
    for row in rows:
        print("  {:<20} Q {:.2f} rise {:.4f} composed {:.3f}".format(
            row["case"], row["q"], row["rise"],
            row["q"] * (row["rise"] ** composed_exponent)))
    print("  raw product spread {:.1%} -> composed {:.1%}  collapse {:.1f}x  verdict: {}".format(
        record["raw_product_spread"], record["composed_spread"],
        record["collapse_factor"], record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())