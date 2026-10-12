"""Lead-lag profile between mass inflow and throat balance.

For each family case the two station series are extracted from the
round-17 junction payload: m'(t_i) = geometric mass delivered over
[t_(i-1), t_i], and B(t_i) = balance at the throat. The lag profile is
the Pearson correlation of m'(t_i) with B(t_(i+k)) for k = -2..+2. The
naming is pre-registered: the profile PEAKS at the lag where the two
series bind; a peak at k = +1 says inflow leads balance by one station
interval. Correlation is not causality, and the record says so.

    python scripts/lab.py scripts/derive_nsc_wall_lead_lag.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

JUNCTION = LAB / "results" / "development" / "nsc-wall-junction-family-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-lead-lag-v1.json"
SCHEMA = "NSC-WALL-LEAD-LAG-v1"
LAGS = (-2, -1, 0, 1, 2)


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
    if not JUNCTION.is_file():
        raise SystemExit("missing family junction record: " + str(JUNCTION))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    junction = json.loads(JUNCTION.read_text(encoding="utf-8"))
    result = {}
    for name in junction["cases"]:
        rows = junction_local(name)
        inflow = [row["mass_forward_inflow"] for row in rows[1:]
                  if row["mass_forward_inflow"] is not None]
        balance = [row["balance_at_throat"] for row in rows]
        profile = {}
        for lag in LAGS:
            x = []
            y = []
            for index in range(len(inflow)):
                bound = index + lag
                if 0 <= bound < len(balance):
                    x.append(inflow[index])
                    y.append(balance[bound])
            r, n = pearson(x, y)
            profile[str(lag)] = {"correlation": r, "pairs": n}
        peak = max(
            (lag for lag in LAGS if profile[str(lag)]["correlation"] is not None),
            key=lambda lag: abs(profile[str(lag)]["correlation"]),
            default=None,
        )
        result[name] = {
            "stations": len(rows),
            "lag_profile": profile,
            "peak_lag": peak,
            "reading": (
                "peak at +1: inflow leads balance by one station interval; "
                "peak at 0: simultaneous; correlation is not causality"
            ),
            "caveat": (
                "inflow series is the forward interval of the junction "
                "payload; balance sampled at station times"
            ),
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "cases": result,
        "notes": [
            "correlation alignment, not a causal proof; the round-17 "
            "alignment effect predicted peak at +1 if inflow leads",
        ],
        "input_sha256": {
            "lab/results/development/nsc-wall-junction-family-v1.json":
                chain.path_sha256(JUNCTION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_lead_lag.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_lead_lag.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in result.items():
        profile = entry["lag_profile"]
        text = ", ".join(
            "k{:+d}:{!s}".format(lag,
                                 "None" if profile[str(lag)]["correlation"] is None
                                 else "{:.3f}".format(profile[str(lag)]["correlation"]))
            for lag in LAGS
        )
        print("  {:<20} peak k={!s}  {}".format(name, entry["peak_lag"], text))
    return 0


def junction_local(name):
    # Reread the source record for the per-case rows (not stored flat).
    import glob as _glob
    from pathlib import Path as _Path
    bubble = json.loads(
        (LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json").read_text()
    )
    case = next(c for c in bubble["cases"] if c["case"] == name)
    files = sorted(_glob.glob(
        str(LAB / "results" / "development" / "nsc-rn-neutral-family-v1" / name
            / "station-*.json")))
    rows = []
    for path in files:
        if _Path(path).stem == "station-failure":
            continue
        record = json.loads(_Path(path).read_text(encoding="utf-8"))
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
    interval_mass = [
        current["inner_mass"] - previous["inner_mass"]
        for previous, current in zip(rows, rows[1:])
    ]
    for index, item in enumerate(junction[1:], start=1):
        if index - 1 < len(interval_mass):
            item["mass_forward_inflow"] = interval_mass[index - 1]
    return junction


if __name__ == "__main__":
    raise SystemExit(main())