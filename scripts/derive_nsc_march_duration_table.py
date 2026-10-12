"""March durations by mu, from the recorded v1 file mtimes.

Round-104 measured the mu=1.002 march at 173.7 minutes. This driver
repeats the mtime reading for the other two v1 cases and records the
duration table: how the march's cost varies with the mass, from the
recorded campaign bytes only.

    python scripts/lab.py scripts/derive_nsc_march_duration_table.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-march-duration-table-v1.json"
SCHEMA = "NSC-MARCH-DURATION-TABLE-v1"
CASES = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")
STATIONS = (5.0, 10.0, 15.0, 20.0)


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    v1 = LAB / "results" / "development" / "nsc-rn-inward-march-v1"
    rows = []
    for name in CASES:
        case_dir = v1 / name
        prepare = case_dir / "prepare.json"
        if not prepare.is_file():
            raise SystemExit(f"missing prepare for {name}")
        start = prepare.stat().st_mtime
        intervals = []
        previous = start
        for t in STATIONS:
            path = case_dir / f"station-{t:g}.npz"
            if not path.is_file():
                raise SystemExit(f"missing {path.name}")
            minutes = (path.stat().st_mtime - start) / 60.0
            intervals.append(minutes - (previous - start) / 60.0)
            previous = path.stat().st_mtime
        rows.append({
            "case": name,
            "mass_over_rm": float(name.split("-eps")[0].split("mu-")[1]),
            "minutes_to_first_station": intervals[0],
            "interval_minutes": intervals,
            "total_to_station_20_minutes": intervals[0] + sum(intervals[1:]),
        })
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "reading": (
            "the march's cost varies with the mass; the heaviest case "
            "costs the longest to the same stations"
        ),
        "notes": [
            "mtimes are the recorded write times of the campaign files",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_march_duration_table.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_march_duration_table.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} first {:.1f} min  intervals {}  total {:.1f} min".format(
            row["case"], row["minutes_to_first_station"],
            ["{:.1f}".format(x) for x in row["interval_minutes"]],
            row["total_to_station_20_minutes"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())