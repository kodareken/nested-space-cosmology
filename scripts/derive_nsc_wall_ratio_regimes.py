"""Ratio regimes: is the scattered mass-energy exchange rate split by the wall?

Follows the round-20 honest negative. The station ratio
R_i = -dE_normal/dM_inner is recomputed per family case and each interval
is assigned to the regime of its END station's bubble wall flag. If the
converter has two settings, the within-regime scatter collapses relative
to the pooled values, and the two settings become two numbers.

    python scripts/lab.py scripts/derive_nsc_wall_ratio_regimes.py
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

FAMILY = LAB / "results" / "development" / "nsc-rn-neutral-family-v1"
BUBBLE = LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-ratio-regimes-v1.json"
SCHEMA = "NSC-WALL-RATIO-REGIMES-v1"
CASES = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")


def stats(values):
    if not values:
        return None
    mean = sum(values) / len(values)
    variance = sum((v - mean) ** 2 for v in values) / len(values)
    std = variance ** 0.5
    return {
        "count": len(values),
        "mean": mean,
        "std": std,
        "relative_scatter": (std / abs(mean)) if mean else None,
    }


def main():
    if not BUBBLE.is_file():
        raise SystemExit("missing bubble record: " + str(BUBBLE))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    cases = {}
    for name in CASES:
        wall_case = next(c for c in bubble["cases"] if c["case"] == name)
        flags = {}
        for wall_row in wall_case["stations"]:
            flags[round(float(wall_row["t_over_rm"]), 6)] = int(wall_row["wall_count"])
        files = sorted(glob.glob(str(FAMILY / name / "station-*.json")))
        rows = []
        for path in files:
            if Path(path).stem == "station-failure":
                continue
            record = json.loads(Path(path).read_text(encoding="utf-8"))
            sample = record["sample"]
            rows.append({
                "t_over_rm": float(sample["t_over_rm"]),
                "inner_mass": float(sample["coupled_inner_mass"]),
                "E_normal": float(sample["coupled_E_normal"]),
            })
        rows.sort(key=lambda row: row["t_over_rm"])
        wall_ratios = []
        quiet_ratios = []
        pooled = []
        for previous, current in zip(rows, rows[1:]):
            dm = current["inner_mass"] - previous["inner_mass"]
            de = current["E_normal"] - previous["E_normal"]
            if abs(dm) < 1e-18:
                continue
            ratio = -de / dm
            regime = flags.get(round(current["t_over_rm"], 6), 0)
            pooled.append(ratio)
            (wall_ratios if regime > 0 else quiet_ratios).append(ratio)
        cases[name] = {
            "pooled": stats(pooled),
            "wall_destined": stats(wall_ratios),
            "quiet_destined": stats(quiet_ratios),
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "cases": cases,
        "reading": (
            "regime = wall flag of the interval's end station; a collapse "
            "of within-regime scatter below the pooled value means the "
            "converter has two settings"
        ),
        "notes": [
            "no parameter is moved; the split is the recorded bubble wall",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_ratio_regimes.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_ratio_regimes.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in cases.items():
        pooled = entry["pooled"]
        wall = entry["wall_destined"]
        quiet = entry["quiet_destined"]
        def fmt(item):
            if item is None:
                return "none"
            return "n={} mean={: .3f} scatter={!s}".format(
                item["count"], item["mean"],
                "none" if item["relative_scatter"] is None
                else "{:.0%}".format(item["relative_scatter"]))
        print("  {:<20} pooled[{}]".format(name, fmt(pooled)))
        print("  {:<20} wall  [{}]  quiet[{}]".format("", fmt(wall), fmt(quiet)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())