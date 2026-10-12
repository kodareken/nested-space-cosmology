"""Mass-energy ratio of the wall channel versus 1/G_N.

Round-19's anti-correlated pair invites the ratio test. The action's mass
function updates as G_N times the integrated Hamiltonian density, so when
the channel is one conserved transport the station ratio

    R_i = - dE_normal / dM_inner

should be close to 1/G_N up to lapse-profile structure factors, and close
to a constant across stations. This driver measures the mean and scatter
of R_i per family case and holds it against the audited 1/G_N.

    python scripts/lab.py scripts/derive_nsc_wall_energy_mass_ratio.py
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
OUTPUT = LAB / "results" / "development" / "nsc-wall-energy-mass-ratio-v1.json"
SCHEMA = "NSC-WALL-ENERGY-MASS-RATIO-v1"
A = 0.04501936182826115
CASES = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    g_n = 1.0 / (16.0 * 3.141592653589793 * A)
    cases = {}
    for name in CASES:
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
        ratios = []
        for previous, current in zip(rows, rows[1:]):
            dm = current["inner_mass"] - previous["inner_mass"]
            de = current["E_normal"] - previous["E_normal"]
            if abs(dm) < 1e-18:
                continue
            ratios.append({
                "t_from": previous["t_over_rm"],
                "t_to": current["t_over_rm"],
                "dM": dm,
                "dE": de,
                "ratio": -de / dm,
            })
        values = [item["ratio"] for item in ratios]
        mean = sum(values) / len(values)
        variance = sum((v - mean) ** 2 for v in values) / len(values)
        std = variance ** 0.5
        cases[name] = {
            "intervals": len(ratios),
            "ratio_mean": mean,
            "ratio_std": std,
            "ratio_relative_scatter": (std / abs(mean)) if mean else None,
            "over_g_n_inverse": mean / (1.0 / g_n),
            "g_n": g_n,
            "g_n_inverse": 1.0 / g_n,
            "rows": ratios,
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "cases": cases,
        "audited_A": A,
        "notes": [
            "R_i = -dE_normal / dM_inner per station interval; the mass "
            "function update is G_N times the integrated Hamiltonian "
            "density, so R ~ 1/G_N is the one-transport reading",
            "mean/relative scatter measured; closeness to 1/G_N includes "
            "lapse-profile structure factors and is reported, not asserted",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_energy_mass_ratio.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_energy_mass_ratio.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("1/G_N = {:.6f}  (G_N = {:.6f})".format(1.0 / g_n, g_n))
    for name, entry in cases.items():
        print(
            "  {:<20} intervals {:>3} R_mean {: .4f} std {:.4f} "
            "scatter {:.2%} R/(1/G_N) {:.4f}".format(
                name, entry["intervals"], entry["ratio_mean"], entry["ratio_std"],
                entry["ratio_relative_scatter"], entry["over_g_n_inverse"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())