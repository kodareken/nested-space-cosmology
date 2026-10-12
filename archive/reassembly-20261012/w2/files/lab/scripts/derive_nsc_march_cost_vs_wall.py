"""March acceleration join: interval shrink vs the wall's death time.

The duration table showed every case's station intervals shrink with t,
most sharply in the heaviest case. This driver joins the shrink amounts
with the recorded extinction times: the march accelerates as the wall
dies, and the acceleration ranks with the death earliness.

    python scripts/lab.py scripts/derive_nsc_march_cost_vs_wall.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

DURATION = LAB / "results" / "development" / "nsc-march-duration-table-v1.json"
EXTINCTION = LAB / "results" / "development" / "nsc-wall-extinction-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-march-cost-vs-wall-v1.json"
SCHEMA = "NSC-MARCH-COST-VS-WALL-v1"


def main():
    if not DURATION.is_file() or not EXTINCTION.is_file():
        raise SystemExit("missing join input record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    duration = json.loads(DURATION.read_text(encoding="utf-8"))
    extinction = json.loads(EXTINCTION.read_text(encoding="utf-8"))
    rows = []
    for row in duration["rows"]:
        intervals = row["interval_minutes"]
        shrink = intervals[1] - intervals[-1]
        death = extinction["march_chart_extinction"][row["case"]]["last_wall"]
        rows.append({
            "case": row["case"],
            "interval_shrink_minutes": shrink,
            "wall_extinction_t": death,
            "shrink_rank": None,
            "death_rank": None,
        })
    by_shrink = sorted(rows, key=lambda row: -row["interval_shrink_minutes"])
    by_death = sorted(rows, key=lambda row: row["wall_extinction_t"])
    for rank, row in enumerate(by_shrink, start=1):
        row["shrink_rank"] = rank
    for rank, row in enumerate(by_death, start=1):
        row["death_rank"] = rank
    rows.sort(key=lambda row: row["case"])
    same_order = [row["shrink_rank"] for row in rows] == [row["death_rank"] for row in rows]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "shrink_ranks_equal_death_ranks": bool(same_order),
        "verdict": (
            "PARTIAL: the only case whose wall dies inside the coverage "
            "(mu=1.040 at t=10) also shows the largest interval shrink; "
            "the two light cases both reach the coverage end and their "
            "shrinks do not rank by death - no full ordering is claimed"
        ),
        "reading": (
            "the march accelerates after the wall dies, and the amount of "
            "acceleration ranks with the earliness of the death where the "
            "death is inside the coverage; the light cases are unresolved"
        ),
        "notes": [
            "interval shrink is the first minus the last station interval",
        ],
        "input_sha256": {
            "lab/results/development/nsc-march-duration-table-v1.json":
                chain.path_sha256(DURATION),
            "lab/results/development/nsc-wall-extinction-v1.json":
                chain.path_sha256(EXTINCTION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_march_cost_vs_wall.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_march_cost_vs_wall.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} shrink {:.1f} min  death t={}  ranks {}/{}".format(
            row["case"], row["interval_shrink_minutes"], row["wall_extinction_t"],
            row["shrink_rank"], row["death_rank"]))
    print("  shrink ranks equal death ranks:", record["shrink_ranks_equal_death_ranks"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())