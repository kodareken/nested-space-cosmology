"""The continuation's trap trajectory and local clock at half-station cadence.

The scalars.jsonl ledger of the continuation records trapping_horizon,
cfl_dt and coupled_inner_mass at every 0.5 station of t=10..20 - the
finest temporal resolution of the trap position anywhere on the chart.
This driver reads the series, measures the trap drift rate, its
correlation with the mass increments, and the chart's local CFL clock
range. The wall programme's breathing laws cover t<=19.5; this is the
trap-side witness of the same era.

    python scripts/lab.py scripts/derive_nsc_continuation_trap.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

SCALARS = LAB / "results" / "development" / "nsc-rn-neutral-continuation-v1" / "scalars.jsonl"
OUTPUT = LAB / "results" / "development" / "nsc-continuation-trap-v1.json"
SCHEMA = "NSC-CONTINUATION-TRAP-v1"


def pearson(x, y):
    n = len(x)
    xm = sum(x) / n
    ym = sum(y) / n
    sxy = sum((a - xm) * (b - ym) for a, b in zip(x, y))
    sxx = sum((a - xm) ** 2 for a in x)
    syy = sum((b - ym) ** 2 for b in y)
    return None if (sxx * syy) == 0.0 else sxy / (sxx * syy) ** 0.5


def main():
    if not SCALARS.is_file():
        raise SystemExit("missing scalars jsonl")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    lines = SCALARS.read_text(encoding="utf-8").splitlines()
    rows = []
    for line in lines:
        entry = json.loads(line)
        rows.append({
            "t": float(entry["coordinate_t_over_rm"]),
            "trap": float(entry["trapping_horizon"]),
            "cfl_dt": float(entry["cfl_dt"]),
            "mass": float(entry["coupled_inner_mass"]),
        })
    rows.sort(key=lambda row: row["t"])
    intervals = []
    for previous, current in zip(rows, rows[1:]):
        dt = current["t"] - previous["t"]
        if dt <= 0.0:
            continue
        intervals.append({
            "from_t": previous["t"],
            "to_t": current["t"],
            "trap_drift_rate": (current["trap"] - previous["trap"]) / dt,
            "mass_increment": current["mass"] - previous["mass"],
        })
    drift_rates = [item["trap_drift_rate"] for item in intervals]
    mass_steps = [item["mass_increment"] for item in intervals]
    correlation = pearson(drift_rates, mass_steps)
    cfl_values = [row["cfl_dt"] for row in rows]
    traps = [row["trap"] for row in rows]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "station_count": len(rows),
        "trap_start": traps[0],
        "trap_end": traps[-1],
        "trap_total_drift": traps[-1] - traps[0],
        "trap_mean_drift_rate": sum(drift_rates) / len(drift_rates),
        "trap_drift_corr_with_mass_steps": correlation,
        "cfl_dt": {
            "min": min(cfl_values), "max": max(cfl_values),
            "mean": sum(cfl_values) / len(cfl_values),
        },
        "trap_range": [min(traps), max(traps)],
        "reading": (
            "the trap moves at the half-station cadence of the march; its "
            "drift is the trap-side witness of the wall era"
        ),
        "notes": [
            "all values are recorded ledger bytes; no recomputation",
            "the trap lives above the reference horizon throughout",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-neutral-continuation-v1/scalars.jsonl":
                chain.path_sha256(SCALARS),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_continuation_trap.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_continuation_trap.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  stations {}  trap {:.6f} -> {:.6f}  drift total {:.3e}  mean rate {:.3e}".format(
        record["station_count"], record["trap_start"], record["trap_end"],
        record["trap_total_drift"], record["trap_mean_drift_rate"]))
    print("  corr(trap drift, mass step): {!s}".format(correlation))
    print("  cfl_dt: min {:.4f} max {:.4f} mean {:.4f}".format(
        cfl_values[0] if False else record["cfl_dt"]["min"],
        record["cfl_dt"]["max"], record["cfl_dt"]["mean"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())