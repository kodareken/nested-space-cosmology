"""Energy side-channel of the wall: does energy arrive with the mass?

Pre-registered from the round-18 lag profile: the wall balance drops one
station after the mass arrives. This driver pulls the energy ledger from
the same family station files and asks

    (a) does the normal-energy increment e'(t_i) correlate with the mass
        increment m'(t_i) at the same station (energy rides the mass)?
    (b) where does the e' vs balance lag profile peak?
    (c) does energy lead or trail the mass (cross-lag m' vs e')?

Correlation is alignment evidence, not a causal proof; the record says so.

    python scripts/lab.py scripts/derive_nsc_wall_energy_channel.py
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
OUTPUT = LAB / "results" / "development" / "nsc-wall-energy-channel-v1.json"
SCHEMA = "NSC-WALL-ENERGY-CHANNEL-v1"
CASES = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")
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
    if not BUBBLE.is_file():
        raise SystemExit("missing bubble record: " + str(BUBBLE))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    cases = {}
    for name in CASES:
        wall_case = next(c for c in bubble["cases"] if c["case"] == name)
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
                "E_normal": float(sample["coupled_E_normal"]),
            })
        rows.sort(key=lambda row: row["t_over_rm"])
        by_t = {round(row["t_over_rm"], 6): row for row in rows}
        wall_flags = {}
        for wall_row in wall_case["stations"]:
            t = round(float(wall_row["t_over_rm"]), 6)
            if t in by_t:
                wall_flags[t] = int(wall_row["wall_count"])
        mass_inc = [current["inner_mass"] - previous["inner_mass"]
                    for previous, current in zip(rows, rows[1:])]
        energy_inc = [current["E_normal"] - previous["E_normal"]
                      for previous, current in zip(rows, rows[1:])]
        balance = [row["t_over_rm"] for row in rows]  # placeholder; balance per rows below
        flags = [wall_flags.get(round(row["t_over_rm"], 6), 0) for row in rows]
        # balance series from the bubble case at matching station times
        balance_series = []
        for row in rows:
            entry = next(
                (w for w in wall_case["stations"]
                 if round(float(w["t_over_rm"]), 6) == round(row["t_over_rm"], 6)),
                None,
            )
            balance_series.append(None if entry is None
                                  else float(entry["balance_at_throat"]))
        same_r, same_n = pearson(mass_inc, energy_inc)
        wall_mass = [m for m, flag in zip(mass_inc, flags[1:]) if flag > 0]
        wall_energy = [e for e, flag in zip(energy_inc, flags[1:]) if flag > 0]
        wall_same = pearson(wall_mass, wall_energy)
        energy_balance_lag = {}
        for lag in LAGS:
            x = []
            y = []
            for index in range(len(energy_inc)):
                bound = index + lag + 0
                balance_index = bound + 1
                if 0 <= balance_index < len(balance_series) and balance_series[balance_index] is not None:
                    x.append(energy_inc[index])
                    y.append(balance_series[balance_index])
            r, n = pearson(x, y)
            energy_balance_lag[str(lag)] = {"correlation": r, "pairs": n}
        mass_energy_lag = {}
        for lag in (-1, 0, 1):
            x = []
            y = []
            for index in range(len(mass_inc)):
                bound = index + lag
                if 0 <= bound < len(energy_inc):
                    x.append(mass_inc[index])
                    y.append(energy_inc[bound])
            r, n = pearson(x, y)
            mass_energy_lag[str(lag)] = {"correlation": r, "pairs": n}
        cases[name] = {
            "stations": len(rows),
            "mass_energy_same_station": {"correlation": same_r, "pairs": same_n},
            "mass_energy_at_wall_stations": {
                "correlation": None if wall_same is None else wall_same[0],
                "pairs": 0 if wall_same is None else wall_same[1],
            },
            "energy_balance_lag": energy_balance_lag,
            "mass_energy_cross_lag": mass_energy_lag,
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "cases": cases,
        "notes": [
            "e' = consecutive difference of the recorded coupled_E_normal",
            "the energy_balance lag index reuses the round-18 alignment "
            "convention; correlation is alignment evidence, not causality",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_energy_channel.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_energy_channel.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in cases.items():
        same = entry["mass_energy_same_station"]["correlation"]
        same_text = "None" if same is None else "{:.3f}".format(same)
        wall_text = "None" if entry["mass_energy_at_wall_stations"]["correlation"] is None \
            else "{:.3f}".format(entry["mass_energy_at_wall_stations"]["correlation"])
        lag_text = ", ".join(
            "k{:+d}:{!s}".format(lag,
                                 "None" if entry["energy_balance_lag"][str(lag)]["correlation"] is None
                                 else "{:.3f}".format(entry["energy_balance_lag"][str(lag)]["correlation"]))
            for lag in LAGS
        )
        cross_text = ", ".join(
            "k{:+d}:{!s}".format(lag,
                                 "None" if entry["mass_energy_cross_lag"][str(lag)]["correlation"] is None
                                 else "{:.3f}".format(entry["mass_energy_cross_lag"][str(lag)]["correlation"]))
            for lag in (-1, 0, 1)
        )
        print("  {:<20} same {:<6} wall {:<6} e'vsB [{}]  m'vse' [{}]".format(
            name, same_text, wall_text, lag_text, cross_text))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())