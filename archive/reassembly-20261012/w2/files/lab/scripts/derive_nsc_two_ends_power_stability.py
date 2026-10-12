"""Temporal stability of the Q gap power across recorded stations.

Round-52 measured p_Q = 0.1173 from the last-station ratios. This driver
repeats the power at every recorded station time of the three cases
(only the earned ingredients: per-station quality ratios and the case
gaps). A stable p_Q across station times makes the round-52 covariance a
property of the run; a drifting power would flag it snapshot-fluked.

    python scripts/lab.py scripts/derive_nsc_two_ends_power_stability.py
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

TWO_ENDS = LAB / "results" / "development" / "nsc-rn-two-ends-v2.json"
OUTPUT = LAB / "results" / "development" / "nsc-two-ends-power-stability-v1.json"
SCHEMA = "NSC-TWO-ENDS-POWER-STABILITY-v1"
MASSES = {"mu-1.002-eps-0.001": 1.002, "mu-1.010-eps-0.001": 1.010,
          "mu-1.040-eps-0.001": 1.040}
LN_GAP_RATIO = 1.0  # placeholder, replaced below


def main():
    if not TWO_ENDS.is_file():
        raise SystemExit("missing two-ends record: " + str(TWO_ENDS))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = json.loads(TWO_ENDS.read_text(encoding="utf-8"))
    gaps = {}
    for name, mass in MASSES.items():
        model = source.action_reference(mass)
        gaps[name] = float(model.horizons["r_plus"]) - float(model.horizons["r_minus"])
    ln_gap_ratio = math.log(gaps["mu-1.040-eps-0.001"] / gaps["mu-1.002-eps-0.001"])
    times = sorted({
        round(float(station["t_over_rm"]), 6)
        for case in record["cases"] for station in case["stations"]
    })
    rows = []
    for t in times:
        values = {}
        for case in record["cases"]:
            station = next(
                (s for s in case["stations"] if abs(float(s["t_over_rm"]) - t) < 1e-6), None)
            if station is not None:
                values[case["case"]] = float(station["quality_ratio_inner_over_outer"])
        light = values.get("mu-1.002-eps-0.001")
        heavy = values.get("mu-1.040-eps-0.001")
        if light and heavy:
            power = math.log(heavy / light) / ln_gap_ratio
            rows.append({"t_over_rm": t, "p_Q": power,
                         "ratio_light": light, "ratio_heavy": heavy})
    powers = [row["p_Q"] for row in rows]
    mean = sum(powers) / len(powers)
    spread = (sum((v - mean) ** 2 for v in powers) / len(powers)) ** 0.5 / abs(mean)
    record_out = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "p_Q_mean_across_stations": mean,
        "p_Q_spread": spread,
        "verdict": (
            "POWER_STABLE" if spread < 0.15 else "POWER_DRIFTS"
        ),
        "notes": [
            "p_Q(t) built from the recorded per-station inner/outer ratios "
            "of the light and heavy cases at matched station times",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-two-ends-v2.json": chain.path_sha256(TWO_ENDS),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_two_ends_power_stability.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_two_ends_power_stability.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record_out, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  t={:<5} p_Q {:.4f}  light {:.3f} heavy {:.3f}".format(
            row["t_over_rm"], row["p_Q"], row["ratio_light"], row["ratio_heavy"]))
    print("  mean {:.4f}  spread {:.1%}  verdict: {}".format(
        record_out["p_Q_mean_across_stations"], record_out["p_Q_spread"],
        record_out["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())