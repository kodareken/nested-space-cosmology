"""Relaxation campaign durations from recorded mtimes.

The v1 duration table came from the march campaign. The relaxation
campaign (six stations through t=80 for mu=1.040) carries its own
mtimes: a second, longer path through the cost law. This driver
extracts the per-station intervals and the t80 era's continuation of
the interval decay.

    python scripts/lab.py scripts/derive_nsc_relaxation_duration.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-relaxation-duration-v1.json"
SCHEMA = "NSC-RELAXATION-DURATION-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    case_dir = (LAB / "results" / "development" / "nsc-rn-relaxation-v1"
                / "mu-1.040-eps-0.001")
    prepare = case_dir / "prepare.json"
    if not prepare.is_file():
        raise SystemExit("missing relaxation prepare.json")
    start = prepare.stat().st_mtime
    stations = [5.0, 10.0, 15.0, 20.0, 40.0, 80.0]
    rows = []
    previous = start
    for t in stations:
        path = case_dir / f"station-{t:g}.npz"
        if not path.is_file():
            raise SystemExit(f"missing relaxation station {t}")
        minutes_total = (path.stat().st_mtime - start) / 60.0
        minutes_interval = (path.stat().st_mtime - previous) / 60.0
        rows.append({
            "t_over_rm": t,
            "minutes_from_prepare": minutes_total,
            "interval_minutes_plus": minutes_interval,
        })
        previous = path.stat().st_mtime
    per_unit = []
    for left, right in zip(rows, rows[1:]):
        dt = right["t_over_rm"] - left["t_over_rm"]
        per_unit.append({
            "window": [left["t_over_rm"], right["t_over_rm"]],
            "minutes": right["interval_minutes_plus"],
            "minutes_per_t": right["interval_minutes_plus"] / dt,
        })
    late = per_unit[-2:]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "per_unit": per_unit,
        "late_windows_t40_t80": late,
        "total_prepare_to_t80_minutes": rows[-1]["minutes_from_prepare"],
        "reading": (
            "the t80-era march cost per unit of t, measured from the "
            "relaxation campaign's own mtimes"
        ),
        "notes": [
            "mtimes from the recorded campaign directory",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_relaxation_duration.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_relaxation_duration.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  t={:<4} total {:7.1f} min  interval {:6.1f}".format(
            row["t_over_rm"], row["minutes_from_prepare"], row["interval_minutes_plus"]))
    for entry in per_unit:
        print("  window {} -> {}: {:.1f} min ({:.1f} min per t)".format(
            *entry["window"], entry["minutes"], entry["minutes_per_t"]))
    print("  t80 total: {:.1f} min".format(record["total_prepare_to_t80_minutes"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())