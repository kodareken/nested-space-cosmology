"""Station-level interior closure and its junction with the force wall.

Extends the run-level identity (round 14) to every station of the
mu-1.040 family (39 stations, t = 0.5..19.5) and joins the overlap
t <= 8.0 with the gradient-bubble v5 wall flags, asking the pre-
registered junction question:

    does the throat force wall live where the excision inflow arrives?

    python scripts/lab.py scripts/derive_nsc_interior_station_balance.py
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

FAMILY_DIR = LAB / "results" / "development" / "nsc-rn-neutral-family-v1" / "mu-1.040-eps-0.001"
BUBBLE = LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json"
OUTPUT = LAB / "results" / "development" / "nsc-interior-station-balance-v1.json"
SCHEMA = "NSC-INTERIOR-STATION-BALANCE-v1"
CLOSURE_TOLERANCE = 1.0e-9


def main():
    if not BUBBLE.is_file():
        raise SystemExit("missing bubble record: " + str(BUBBLE))
    files = sorted(glob.glob(str(FAMILY_DIR / "station-*.json")))
    if not files:
        raise SystemExit("no family station files under " + str(FAMILY_DIR))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    case = next(c for c in bubble["cases"] if c["case"] == "mu-1.040-eps-0.001")
    rows = []
    for path in files:
        if Path(path).stem == "station-failure":
            continue
        record = json.loads(Path(path).read_text(encoding="utf-8"))
        sample = record["sample"]
        rows.append({
            "t_over_rm": float(sample["t_over_rm"]),
            "clock_t": float(sample["clock_t"]),
            "inner_mass": float(sample["coupled_inner_mass"]),
            "excision_probability": float(sample["coupled_excision"]),
            "outer_probability": float(sample["coupled_outer"]),
            "remaining_probability": float(sample["coupled_remaining"]),
            "backreaction_delta_inner_mass": float(sample["delta_inner_mass"]),
            "backreaction_delta_excision": float(sample["delta_excision"]),
            "recorded_closure_residual": float(sample["closure_residual"]),
        })
    rows.sort(key=lambda row: row["t_over_rm"])
    intervals = []
    for previous, current in zip(rows, rows[1:]):
        dm = current["inner_mass"] - previous["inner_mass"]
        intervals.append({
            "interval": [previous["t_over_rm"], current["t_over_rm"]],
            "delta_geometric_mass": dm,
        })
    # Book check: remaining + excision + outer should partition probability.
    partition_gaps = [
        1.0 - row["remaining_probability"] - row["excision_probability"] - row["outer_probability"]
        for row in rows
    ]
    max_partition_gap = max(abs(value) for value in partition_gaps)
    junction = []
    rows_by_t = {round(row["t_over_rm"], 6): row for row in rows}
    for wall_row in case["stations"]:
        t = round(float(wall_row["t_over_rm"]), 6)
        row = rows_by_t.get(t)
        if row is None:
            continue
        interval_row = next(
            (item for item in intervals if abs(item["interval"][0] - t) < 1e-9), None)
        inflow = None if interval_row is None else interval_row["delta_geometric_mass"]
        junction.append({
            "t_over_rm": t,
            "wall_count": int(wall_row["wall_count"]),
            "balance_at_throat": float(wall_row["balance_at_throat"]),
            "throat_slope_inner": float(wall_row["throat_slope_inner"]),
            "interval_geometric_mass_inflow": inflow,
        })
    wall_rows = [j for j in junction if j["wall_count"] > 0]
    quiet_rows = [j for j in junction if j["wall_count"] == 0]
    wall_inflow = [abs(j["interval_geometric_mass_inflow"]) for j in wall_rows
                   if j["interval_geometric_mass_inflow"] is not None]
    quiet_inflow = [abs(j["interval_geometric_mass_inflow"]) for j in quiet_rows
                    if j["interval_geometric_mass_inflow"] is not None]
    mean = lambda items: (sum(items) / len(items)) if items else None
    wall_ratio = (mean(wall_inflow) / mean(quiet_inflow)) if quiet_rows else None
    import math
    pairs = [(j["throat_slope_inner"], j["interval_geometric_mass_inflow"])
             for j in junction if j["interval_geometric_mass_inflow"] is not None]
    if len(pairs) > 2:
        xm = sum(a for a, _ in pairs) / len(pairs)
        ym = sum(b for _, b in pairs) / len(pairs)
        num = sum((a - xm) * (b - ym) for a, b in pairs)
        den_x = sum((a - xm) ** 2 for a, _ in pairs)
        den_y = sum((b - ym) ** 2 for _, b in pairs)
        correlation = None if (den_x * den_y) == 0.0 else num / math.sqrt(den_x * den_y)
    else:
        correlation = None

    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "station_files": len(rows),
        "closure_intervals": len(intervals),
        "first_t": rows[0]["t_over_rm"] if rows else None,
        "last_t": rows[-1]["t_over_rm"] if rows else None,
        "book_reading": (
            "the family runner's coupled_excision/coupled_outer partition "
            "probability (with coupled_remaining); the geometric mass book "
            "is coupled_inner_mass, whose intervals feed this junction"
        ),
        "probability_partition_max_gap": max_partition_gap,
        "junction": {
            "overlap_stations": len(junction),
            "wall_stations": len(wall_rows),
            "quiet_stations": len(quiet_rows),
            "wall_quiet_mass_inflow_ratio": wall_ratio,
            "slope_mass_inflow_correlation": correlation,
            "rows": junction,
        },
        "notes": [
            "the family's own delta_* fields are coupled-minus-fixed "
            "backreaction and are recorded only for provenance",
            "the junction covers the family/bubble overlap stations; the "
            "correlation is an indicator over these stations, not a fit",
            "wall_count comes from the recorded bubble v5 flags",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_interior_station_balance.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_interior_station_balance.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("stations: {} ({} mass intervals) first {:.1f} last {:.1f}".format(
        len(rows), len(intervals), rows[0]["t_over_rm"], rows[-1]["t_over_rm"]))
    print("probability partition max gap: {:.3e}".format(max_partition_gap))
    correlation_text = "None" if correlation is None else "{:.4f}".format(correlation)
    print("junction: overlap {} wall {} quiet {} mass inflow ratio {} slope-mass corr {}".format(
        len(junction), len(wall_rows), len(quiet_rows), wall_ratio, correlation_text))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())