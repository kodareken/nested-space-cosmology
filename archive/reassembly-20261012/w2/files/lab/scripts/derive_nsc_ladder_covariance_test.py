"""Predictive test of the 4*p_Q = b identity on the covariance.

Round-58 flagged the sub-percent proximity between 4*p_Q and the
kappa-ladder's quadratic coefficient b. The pre-registered test: take
the identity as true (p_Q := b/4, no freedom), close the round-52
covariance with that exponent, and measure the resulting spread of
Q * rise^(2*p_Q_predicted) against the recorded 2.1% of the measured
exponent. A spread below 5% is the weak predictive pass; the force of
the test is stated honestly (three gap points).

    python scripts/lab.py scripts/derive_nsc_ladder_covariance_test.py
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

LADDER = LAB / "results" / "development" / "nsc-charged-ladder-deficit-law-v1.json"
BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
COVARIANCE = LAB / "results" / "development" / "nsc-two-ends-covariance-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-ladder-covariance-test-v1.json"
SCHEMA = "NSC-LADDER-COVARIANCE-TEST-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}


def spread(values):
    mean = sum(values) / len(values)
    return (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5 / abs(mean)


def main():
    if any(not p.is_file() for p in (LADDER, BREATHING, COVARIANCE)):
        raise SystemExit("missing input record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    ladder = json.loads(LADDER.read_text(encoding="utf-8"))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    covariance = json.loads(COVARIANCE.read_text(encoding="utf-8"))
    b = ladder["b"]
    p_q_predicted = b / 4.0
    composed_exponent = 2.0 * p_q_predicted
    rows = []
    for entry in covariance["rows"]:
        rise = entry["rise"]
        rows.append({
            "case": entry["case"],
            "q": entry["q"],
            "rise": rise,
            "composed": entry["q"] * (rise ** composed_exponent),
        })
    values = [row["composed"] for row in rows]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "hypothesis_pre_registered": (
            "4*p_Q = b taken as identity; the covariance closed with the "
            "ladder-implied exponent"
        ),
        "ladder_b": b,
        "p_Q_predicted": p_q_predicted,
        "composed_exponent": composed_exponent,
        "rows": rows,
        "predicted_composed_spread": spread(values),
        "recorded_composed_spread_reference": covariance["composed_spread"],
        "verdict": (
            "PREDICTIVE_PASS" if spread(values) < 0.05 else "PREDICTIVE_FAIL"
        ),
        "notes": [
            "three gap points: the pass is weak by construction and the "
            "record says so; nothing beyond the identity is claimed",
        ],
        "input_sha256": {
            "lab/results/development/nsc-charged-ladder-deficit-law-v1.json":
                chain.path_sha256(LADDER),
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
            "lab/results/development/nsc-two-ends-covariance-v1.json":
                chain.path_sha256(COVARIANCE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_ladder_covariance_test.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_ladder_covariance_test.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  p_Q predicted from ladder = b/4 = {:.6f}".format(p_q_predicted))
    for row in rows:
        print("  {:<20} composed {:.3f}".format(row["case"], row["composed"]))
    print("  predicted spread {:.2%} vs recorded {:.2%}  verdict: {}".format(
        record["predicted_composed_spread"], record["recorded_composed_spread_reference"],
        record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())