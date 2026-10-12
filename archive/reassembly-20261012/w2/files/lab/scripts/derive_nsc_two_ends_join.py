"""Three-way join: room sharpness, wall concentration, breathing rate.

The two-ends record measured the inner room 21.1-25.2x sharper than the
outer (Q); the wall programme measured 10.9-22.4x inflow concentration
and the gap-clocked rise rates. They were never laid on one table. This
driver joins the three recorded families case by case and evaluates the
two pre-registered directional signatures: (a) a sharper inner room
corresponds to slower breathing rise; (b) whether Q and the wall
concentration move together. n=3 is stated everywhere.

    python scripts/lab.py scripts/derive_nsc_two_ends_join.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

TWO_ENDS = LAB / "results" / "development" / "nsc-rn-two-ends-v2.json"
WALL = LAB / "results" / "development" / "nsc-wall-junction-family-v1.json"
BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-two-ends-join-v1.json"
SCHEMA = "NSC-TWO-ENDS-JOIN-v1"


def pearson(x, y):
    n = len(x)
    xm = sum(x) / n
    ym = sum(y) / n
    sxy = sum((a - xm) * (b - ym) for a, b in zip(x, y))
    sxx = sum((a - xm) ** 2 for a in x)
    syy = sum((b - ym) ** 2 for b in y)
    return None if (sxx * syy) == 0.0 else sxy / (sxx * syy) ** 0.5


def main():
    if any(not p.is_file() for p in (TWO_ENDS, WALL, BREATHING)):
        raise SystemExit("missing join input record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    two_ends = json.loads(TWO_ENDS.read_text(encoding="utf-8"))
    wall = json.loads(WALL.read_text(encoding="utf-8"))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    q_summary = two_ends["summary"]["quality_ratio_inner_over_outer_last"]
    cases = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")
    table = []
    for name in cases:
        table.append({
            "case": name,
            "room_quality_ratio": q_summary[name],
            "wall_inflow_ratio": wall["cases"][name]["wall_quiet_mass_inflow_ratio"],
            "rise_rate": breathing["cases"][name]["mean_rise_rate"],
        })
    q_values = [row["room_quality_ratio"] for row in table]
    wall_values = [row["wall_inflow_ratio"] for row in table]
    rise_values = [row["rise_rate"] for row in table]
    corr_q_wall = pearson(q_values, wall_values)
    corr_q_rise = pearson(q_values, rise_values)
    corr_wall_rise = pearson(wall_values, rise_values)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "table": table,
        "correlations": {
            "q_vs_wall_inflow": corr_q_wall,
            "q_vs_rise": corr_q_rise,
            "wall_inflow_vs_rise": corr_wall_rise,
        },
        "pre_registered_signs": {
            "q_vs_rise_negative": corr_q_rise < 0.0,
            "q_vs_wall_uncommitted": True,
        },
        "verdict": (
            "SHARPER_ROOM_SLOWER_RISE"
            if corr_q_rise is not None and corr_q_rise < -0.5 else "NOT_THIS_JOIN"
        ),
        "notes": [
            "n = 3 cases throughout; correlations are indicators",
            "the two ratios are different instruments (geometric sharpness "
            "vs dynamic concentration) and are not identified",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-two-ends-v2.json": chain.path_sha256(TWO_ENDS),
            "lab/results/development/nsc-wall-junction-family-v1.json": chain.path_sha256(WALL),
            "lab/results/development/nsc-wall-breathing-family-v1.json": chain.path_sha256(BREATHING),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_two_ends_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_two_ends_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in table:
        print("  {:<20} Q {:.2f}  wall_ratio {:.2f}  rise {:.4f}".format(
            row["case"], row["room_quality_ratio"], row["wall_inflow_ratio"],
            row["rise_rate"]))
    print("  corr(Q,wall) {!s}  corr(Q,rise) {!s}  corr(wall,rise) {!s}".format(
        corr_q_wall, corr_q_rise, corr_wall_rise))
    print("  verdict:", record["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())