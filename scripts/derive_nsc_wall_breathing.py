"""Wall breathing: rise/fall radial rates and their inflow correlation.

Round-22 recorded the wall radii breathing inside the window shell
(0.97..1.53 r_+). This driver turns the 12 wall epochs into a radial
motion: per adjacent wall-station pair the migration rate dr_wall/dt,
split into rising and falling phases, and the same-station correlation
of the wall radius with the geometric mass inflow of that interval.

    python scripts/lab.py scripts/derive_nsc_wall_breathing.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

JOIN = LAB / "results" / "development" / "nsc-wall-window-join-v1.json"
JUNCTION = LAB / "results" / "development" / "nsc-interior-station-balance-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-breathing-v1.json"
SCHEMA = "NSC-WALL-BREATHING-v1"


def pearson(x, y):
    n = len(x)
    if n < 3:
        return None, n
    xm = sum(x) / n
    ym = sum(y) / n
    sxy = sum((a - xm) * (b - ym) for a, b in zip(x, y))
    sxx = sum((a - xm) ** 2 for a in x)
    syy = sum((b - ym) ** 2 for b in y)
    if sxx <= 0.0 or syy <= 0.0:
        return None, n
    return sxy / (sxx * syy) ** 0.5, n


def main():
    if not JOIN.is_file() or not JUNCTION.is_file():
        raise SystemExit("missing join or junction record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    join = json.loads(JOIN.read_text(encoding="utf-8"))
    junction = json.loads(JUNCTION.read_text(encoding="utf-8"))
    wall_rows = [row for row in join["rows"] if row["wall_ratio_to_horizon"] is not None]
    wall_rows.sort(key=lambda row: row["t_over_rm"])
    inflow_by_start = {}
    for row in junction["junction"]["rows"]:
        if row["interval_geometric_mass_inflow"] is not None:
            inflow_by_start[round(row["t_over_rm"], 6)] = \
                abs(row["interval_geometric_mass_inflow"])
    velocities = []
    rising = []
    falling = []
    for previous, current in zip(wall_rows, wall_rows[1:]):
        dt = current["t_over_rm"] - previous["t_over_rm"]
        if dt <= 0.0:
            continue
        velocity = (current["wall_ratio_to_horizon"] - previous["wall_ratio_to_horizon"]) / dt
        velocities.append({
            "from_t": previous["t_over_rm"],
            "to_t": current["t_over_rm"],
            "velocity_in_horizon_units_per_t": velocity,
        })
        (rising if velocity > 0.0 else falling).append(velocity)
    same_radius = [row["wall_ratio_to_horizon"] for row in wall_rows]
    same_inflow = [inflow_by_start.get(round(row["t_over_rm"], 6), 0.0) for row in wall_rows]
    lag_inflow = [inflow_by_start.get(round(row["t_over_rm"] - 0.5, 6), 0.0)
                  for row in wall_rows]
    same_corr = pearson(same_radius, same_inflow)
    lag_corr = pearson(same_radius, lag_inflow)
    mean = lambda items: (sum(items) / len(items)) if items else None
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "wall_epochs": len(wall_rows),
        "radial_range": [
            min(row["wall_ratio_to_horizon"] for row in wall_rows),
            max(row["wall_ratio_to_horizon"] for row in wall_rows),
        ],
        "rise_epochs": len(rising),
        "fall_epochs": len(falling),
        "mean_rise_rate": mean(rising),
        "mean_fall_rate": mean(falling),
        "radius_inflow_same_corr": same_corr,
        "radius_inflow_lagged_corr": lag_corr,
        "velocities": velocities,
        "reading": (
            "the wall peaks near t=7.5-8.0 at ~1.53 r_+ and sinks toward "
            "0.97 r_+ by t=10; rates are in horizon units per coordinate t"
        ),
        "notes": [
            "rates use adjacent WALL epochs only; the t=2.5..5.5 gap has no "
            "wall and is not interpolated",
        ],
        "input_sha256": {
            "lab/results/development/nsc-wall-window-join-v1.json": chain.path_sha256(JOIN),
            "lab/results/development/nsc-interior-station-balance-v1.json":
                chain.path_sha256(JUNCTION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_breathing.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_breathing.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  epochs {}  radial range [{:.3f}, {:.3f}] r_+".format(
        record["wall_epochs"], record["radial_range"][0], record["radial_range"][1]))
    print("  rise epochs {} mean rate {: .4f}   fall epochs {} mean rate {: .4f}".format(
        record["rise_epochs"],
        0.0 if record["mean_rise_rate"] is None else record["mean_rise_rate"],
        record["fall_epochs"],
        0.0 if record["mean_fall_rate"] is None else record["mean_fall_rate"],
    ))
    print("  r_wall vs inflow same-station corr: {!s}".format(
        "none" if same_corr is None else "{:.3f} (n={})".format(same_corr[0], same_corr[1])))
    print("  r_wall vs inflow lagged(-0.5) corr: {!s}".format(
        "none" if lag_corr is None else "{:.3f} (n={})".format(lag_corr[0], lag_corr[1])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())