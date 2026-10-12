"""Inflow -> balance response of the wall stations.

Reads the round-15 station junction and pairs, for consecutive stations,
the balance displacement dB_i = B_i - B_{i-1} with the geometric mass
delivered over the preceding interval dM_{i-1}. Least-squares slopes are
reported for all pairs, the wall stations, and the quiet stations. The
slope is a response law; causality is not claimed by construction.

    python scripts/lab.py scripts/derive_nsc_wall_inflow_response.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

JUNCTION = LAB / "results" / "development" / "nsc-interior-station-balance-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-inflow-response-v1.json"
SCHEMA = "NSC-WALL-INFLOW-RESPONSE-v1"


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
    return {
        "slope": slope,
        "correlation": r,
        "pairs": n,
        "mean_dM": xm,
        "mean_dB": ym,
    }


def main():
    if not JUNCTION.is_file():
        raise SystemExit("missing junction record: " + str(JUNCTION))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    junction = json.loads(JUNCTION.read_text(encoding="utf-8"))
    rows = junction["junction"]["rows"]
    rows.sort(key=lambda row: row["t_over_rm"])
    pairs = []
    wall_pairs = []
    quiet_pairs = []
    for previous, current in zip(rows, rows[1:]):
        dm = previous["interval_geometric_mass_inflow"]
        db = current["balance_at_throat"] - previous["balance_at_throat"]
        if dm is None:
            continue
        pairs.append((abs(dm), db))
        (wall_pairs if current["wall_count"] > 0 else quiet_pairs).append((abs(dm), db))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "pairing": "dB_i = B_i - B_(i-1) with the mass over [t_(i-1), t_i]",
        "all": fit(pairs),
        "wall": fit(wall_pairs),
        "quiet": fit(quiet_pairs),
        "residual_medians": None,
        "notes": [
            "response law, not a conversion-efficiency claim; the pairing is "
            "pre-registered in this file's schema before the numbers",
            "abs(dM) enters the x axis: the response measures magnitude",
        ],
        "input_sha256": {
            "lab/results/development/nsc-interior-station-balance-v1.json":
                chain.path_sha256(JUNCTION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_inflow_response.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_inflow_response.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name in ("all", "wall", "quiet"):
        entry = record[name]
        if entry is None:
            print("  {:<6} insufficient pairs".format(name))
            continue
        print(
            "  {:<6} slope={!s} corr={!s} pairs={} mean_dM={:.3e} mean_dB={!s}".format(
                name,
                None if entry["slope"] is None else "{:.4e}".format(entry["slope"]),
                None if entry["correlation"] is None else "{:.4f}".format(entry["correlation"]),
                entry["pairs"], entry["mean_dM"],
                None if entry["mean_dB"] is None else "{:.4f}".format(entry["mean_dB"]),
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())