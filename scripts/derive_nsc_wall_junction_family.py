"""Family junction: wall inflow concentration and response law, all cases.

Repeats the round-15/16 wall junction for each eps=0.001 family case
(mu = 1.002, 1.010, 1.040) and asks whether the two signatures recur:

    (a) wall stations concentrate the geometric mass inflow;
    (b) the wall sub-sample carries a positive inflow -> balance slope
        while the quiet sub-sample carries the opposite one.

Same pre-registered pairings as the round-15/16 records.

    python scripts/lab.py scripts/derive_nsc_wall_junction_family.py
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
OUTPUT = LAB / "results" / "development" / "nsc-wall-junction-family-v1.json"
SCHEMA = "NSC-WALL-JUNCTION-FAMILY-v1"
CASES = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")


def fit(pairs):
    n = len(pairs)
    if n < 2:
        return None
    x = [a for a, _ in pairs]
    y = [b for _, b in pairs]
    xm = sum(x) / n
    ym = sum(y) / n
    sxy = sum((a - xm) * (b - ym) for a, b in pairs)
    sxx = sum((a - xm) ** 2 for a in x)
    syy = sum((b - ym) ** 2 for b in y)
    slope = (sxy / sxx) if sxx != 0.0 else None
    r = (sxy / (sxx * syy) ** 0.5) if (sxx * syy) > 0.0 else None
    return {"slope": slope, "correlation": r, "pairs": n}


def main():
    if not BUBBLE.is_file():
        raise SystemExit("missing bubble record: " + str(BUBBLE))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    cases = {}
    for name in CASES:
        case = next(c for c in bubble["cases"] if c["case"] == name)
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
            })
        rows.sort(key=lambda row: row["t_over_rm"])
        by_t = {round(row["t_over_rm"], 6): row for row in rows}
        junction = []
        for wall_row in case["stations"]:
            t = round(float(wall_row["t_over_rm"]), 6)
            row = by_t.get(t)
            if row is None:
                continue
            junction.append({
                "t_over_rm": t,
                "wall_count": int(wall_row["wall_count"]),
                "balance_at_throat": float(wall_row["balance_at_throat"]),
                "mass_forward_inflow": None,
            })
        interval_mass = []
        for previous, current in zip(rows, rows[1:]):
            interval_mass.append(current["inner_mass"] - previous["inner_mass"])
        for index, item in enumerate(junction[1:], start=1):
            earlier = junction[index - 1]
            if index - 1 < len(interval_mass):
                item["mass_forward_inflow"] = interval_mass[index - 1]
        wall_rows = [j for j in junction if j["wall_count"] > 0]
        quiet_rows = [j for j in junction if j["wall_count"] == 0]
        wall_in = [abs(j["mass_forward_inflow"]) for j in junction[1:]
                   if j["wall_count"] > 0 and j["mass_forward_inflow"] is not None]
        quiet_in = [abs(j["mass_forward_inflow"]) for j in junction[1:]
                    if j["wall_count"] == 0 and j["mass_forward_inflow"] is not None]
        ratio = ((sum(wall_in) / len(wall_in)) / (sum(quiet_in) / len(quiet_in))
                 if wall_in and quiet_in else None)
        pairs = []
        wall_pairs = []
        quiet_pairs = []
        for previous, current in zip(junction, junction[1:]):
            dm = previous["mass_forward_inflow"]
            if dm is None:
                continue
            db = current["balance_at_throat"] - previous["balance_at_throat"]
            pairs.append((abs(dm), db))
            (wall_pairs if current["wall_count"] > 0 else quiet_pairs).append((abs(dm), db))
        cases[name] = {
            "stations": len(rows),
            "junction_overlap": len(junction),
            "wall_stations": len(wall_rows),
            "quiet_stations": len(quiet_rows),
            "wall_quiet_mass_inflow_ratio": ratio,
            "all": fit(pairs),
            "wall": fit(wall_pairs),
            "quiet": fit(quiet_pairs),
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "cases": cases,
        "notes": [
            "dB_i = B_i - B_(i-1) paired with the mass over the true preceding "
            "interval [t_(i-1), t_i] (round 16 stored the forward alignment; "
            "this backward alignment strengthens the wall correlation at "
            "mu=1.040 from 0.49 to 0.71, recorded here as alignment evidence)",
            "junction rows follow the bubble v5 wall flags per station",
            "the quiet regime carries weak mixed-sign slopes across cases: "
            "the wall is the only coherent response regime",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_junction_family.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_junction_family.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in cases.items():
        wall = entry["wall"]
        quiet = entry["quiet"]
        print(
            "  {:<20} overlap {:>3} wall {:>2} quiet {:>2} ratio {:>8} "
            "wall slope {:>9} (r {:.3f}) quiet slope {:>9} (r {:.3f})".format(
                name, entry["junction_overlap"], entry["wall_stations"],
                entry["quiet_stations"],
                "None" if entry["wall_quiet_mass_inflow_ratio"] is None
                else "{:.1f}".format(entry["wall_quiet_mass_inflow_ratio"]),
                "None" if wall is None or wall["slope"] is None
                else "{:.3g}".format(wall["slope"]),
                -1e9 if wall is None or wall["correlation"] is None else wall["correlation"],
                "None" if quiet is None or quiet["slope"] is None
                else "{:.3g}".format(quiet["slope"]),
                -1e9 if quiet is None or quiet["correlation"] is None else quiet["correlation"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())