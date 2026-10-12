"""Angular-ladder cliff search: does reflection ever cross unity?

Finer omega grid for kappa = sqrt(6) (n = 2) than the 0.005 ladder scan,
and the first scan of kappa = sqrt(12) (n = 3). The crossing question is
pre-registered: a channel with reflection above unity inside the electric
window would be a radiative (superradiant) window.

    python scripts/lab.py scripts/derive_nsc_rn_charged_twist_cliff.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_rn_charged_lll as lll  # noqa: E402

A = 0.04501936182826115
C_F = 0.01125484045706527
Q_E = 4.0
MASS_RATIO = 1.01
SIGN = 1
OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-twist-cliff-v1.json"
REFINED_OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-twist-cliff-v2.json"
RUN3_OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-twist-cliff-v3.json"
SCHEMA = "NSC-RN-CHARGED-TWIST-CLIFF-v1"


def main():
    if OUTPUT.exists() and REFINED_OUTPUT.exists() and RUN3_OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing outputs")
    ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=1.0, Q_E=Q_E, mass_over_total_radius=MASS_RATIO,
    )
    grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
    chart = lll.charged_chart(ref["model"], grid["radius"])

    def sweep(kappa, omegas):
        rows = []
        radiative = False
        for omega in omegas:
            result = lll.charged_twisted_column(
                model=ref["model"], radius=grid["radius"], omega=omega, sign=SIGN,
                Q_E=Q_E, C_F=C_F, kappa=kappa,
            )
            value = lll.reflection_ratio(result["mode"], chart=chart)["reflection_ratio"]
            rows.append({"omega": float(omega), "reflection_ratio": value})
            radiative = radiative or bool(value > 1.0 + 1e-12)
        peak = max(rows, key=lambda row: row["reflection_ratio"])
        return {
            "kappa": float(kappa),
            "rows": rows,
            "peak": peak,
            "radiative": bool(radiative),
        }

    fine_6 = tuple(round(0.095 + 0.001 * index, 4) for index in range(31))  # 0.095..0.125
    coarse_12 = tuple(round(0.05 + 0.005 * index, 4) for index in range(21))  # 0.05..0.15
    fine_12 = tuple(round(0.104 + 0.0005 * index, 5) for index in range(41))  # 0.104..0.124
    entries = [
        sweep(6.0 ** 0.5, fine_6),
        sweep(12.0 ** 0.5, coarse_12),
    ]
    refined = {
        "kappa": float(12.0 ** 0.5),
        "grid": "0.0005 steps over 0.104..0.124",
        "entry": sweep(12.0 ** 0.5, fine_12),
    }
    record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "Phi_H": ref["Phi_H"],
        "entries": entries,
        "any_radiative": bool(any(entry["radiative"] for entry in entries)),
        "notes": [
            "crossing pre-registered: reflection above unity means the channel "
            "emits more through the plus characteristic than it ingests",
            "n=1 never crossed (rounds 5-7 records)",
        ],
    }
    import hashlib
    sources = {}
    for relative, target in (
        ("lab/src/recursive_horizons/nsc_rn_charged_lll.py",
         LAB / "src" / "recursive_horizons" / "nsc_rn_charged_lll.py"),
        ("lab/scripts/derive_nsc_rn_charged_twist_cliff.py",
         LAB / "scripts" / "derive_nsc_rn_charged_twist_cliff.py"),
    ):
        digest = hashlib.sha256()
        with target.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        sources[relative] = digest.hexdigest()
    record["source_sha256"] = sources
    if not OUTPUT.exists():
        with OUTPUT.open("w", encoding="utf-8") as handle:
            handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
        print("wrote", OUTPUT.relative_to(LAB))
    for entry in entries:
        print(
            "  kappa={: .4f} peak omega={:.4f} R={:.6f} radiative={} points={}".format(
                entry["kappa"], entry["peak"]["omega"], entry["peak"]["reflection_ratio"],
                entry["radiative"], len(entry["rows"]),
            )
        )
    print("  any_radiative:", record["any_radiative"])

    if not REFINED_OUTPUT.exists():
        refined_record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "supersedes_grid_only": "cliff-v1 kappa=sqrt(12) 0.005 grid",
        "Phi_H": ref["Phi_H"],
        "kappa": refined["kappa"],
        "result": refined["entry"],
        "notes": [
            "finer grid for the n=3 multiplet to rule out an undersampled crossing",
        ],
        "source_sha256": sources,
        }
        with REFINED_OUTPUT.open("w", encoding="utf-8") as handle:
            handle.write(json.dumps(refined_record, indent=2, sort_keys=True, allow_nan=False) + "\n")
        print("wrote", REFINED_OUTPUT.relative_to(LAB))
    entry = refined["entry"]
    print(
        "  kappa={: .4f} fine peak omega={:.5f} R={:.6f} radiative={}".format(
            entry["kappa"], entry["peak"]["omega"], entry["peak"]["reflection_ratio"],
            entry["radiative"],
        )
    )

    if RUN3_OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(RUN3_OUTPUT))
    fine_20 = tuple(round(0.110 + 0.0005 * index, 5) for index in range(33))  # 0.110..0.126
    entry20 = sweep(20.0 ** 0.5, fine_20)
    run3_record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "Phi_H": ref["Phi_H"],
        "kappa": float(20.0 ** 0.5),
        "grid": "0.0005 steps over 0.110..0.126",
        "result": entry20,
        "deficit_sequence": {
            "n1_deficit": 1.0 - 0.667912,
            "n2_deficit": 1.0 - 0.956659,
            "n3_deficit": 1.0 - 0.997778,
            "n4_deficit": 1.0 - entry20["peak"]["reflection_ratio"],
        },
        "notes": [
            "n=4 multiplet; if the deficit collapses to ~1e-4 the ladder "
            "approaches a reflective asymptote from below",
        ],
        "source_sha256": sources,
    }
    with RUN3_OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(run3_record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", RUN3_OUTPUT.relative_to(LAB))
    print(
        "  kappa={: .4f} peak omega={:.5f} R={:.6f} radiative={}".format(
            entry20["kappa"], entry20["peak"]["omega"],
            entry20["peak"]["reflection_ratio"], entry20["radiative"],
        )
    )
    print(
        "  deficits n1..n4: {:.4e} {:.4e} {:.4e} {:.4e}".format(
            1.0 - 0.667912, 1.0 - 0.956659, 1.0 - 0.997778,
            1.0 - entry20["peak"]["reflection_ratio"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())