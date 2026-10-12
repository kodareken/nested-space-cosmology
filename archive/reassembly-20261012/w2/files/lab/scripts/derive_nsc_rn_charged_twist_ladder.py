"""Angular ladder and omega_star relation for the twisted reflection.

Answers two discovery questions on frozen machinery:

    1. Does the reflection resonance track omega_star across slices?
       Scan the peak for eta = 1.01 and eta = 1.5 and compare with the
       phase-averaged potential of each slice.
    2. Does the next multiplet (kappa = sqrt(6), n = 2) keep its
       reflection sub-unity, and where does its peak sit?

    python scripts/lab.py scripts/derive_nsc_rn_charged_twist_ladder.py
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
SIGN = 1
SCAN = tuple(round(0.05 + 0.005 * index, 4) for index in range(21))  # 0.05..0.15
OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-twist-ladder-v1.json"
SCHEMA = "NSC-RN-CHARGED-TWIST-LADDER-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))

    def scan_peak(eta, kappa):
        ref = lll.dyonic_reference(
            A=A, C_F=C_F, q_mag=1.0, Q_E=Q_E, mass_over_total_radius=eta,
        )
        grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
        chart = lll.charged_chart(ref["model"], grid["radius"])
        winding = lll.infall_winding(
            model=ref["model"], omega=0.0, sign=SIGN, Q_E=Q_E, C_F=C_F,
            r_outer=float(grid["outer"]), points=2400,
        )
        best = None
        rows = []
        for omega in SCAN:
            result = lll.charged_twisted_column(
                model=ref["model"], radius=grid["radius"], omega=omega, sign=SIGN,
                Q_E=Q_E, C_F=C_F, kappa=kappa,
            )
            ratio = lll.reflection_ratio(result["mode"], chart=chart)["reflection_ratio"]
            rows.append({"omega": float(omega), "reflection_ratio": ratio})
            if best is None or ratio > best["reflection_ratio"]:
                best = {"omega": float(omega), "reflection_ratio": ratio}
        return {
            "eta": float(eta),
            "kappa": float(kappa),
            "omega_star": float(winding["omega_star"]),
            "omega_star_half": float(winding["omega_star"] / 2.0),
            "Phi_H": float(ref["Phi_H"]),
            "peak": best,
            "peak_over_omega_star": (best["omega"] / winding["omega_star"]) if winding["omega_star"] else None,
            "rows": rows,
        }

    entries = [
        scan_peak(1.01, 2.0 ** 0.5),
        scan_peak(1.5, 2.0 ** 0.5),
        scan_peak(1.01, 6.0 ** 0.5),
    ]
    record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "scan_omega_min": float(SCAN[0]),
        "scan_omega_max": float(SCAN[-1]),
        "entries": entries,
        "notes": [
            "the omega_star/2 comparison is a recorded relation test, not a fit",
            "kappa = sqrt(6) is the n = 2 multiplet of the same lowest shell",
        ],
    }
    import hashlib
    sources = {}
    for relative, target in (
        ("lab/src/recursive_horizons/nsc_rn_charged_lll.py",
         LAB / "src" / "recursive_horizons" / "nsc_rn_charged_lll.py"),
        ("lab/scripts/derive_nsc_rn_charged_twist_ladder.py",
         LAB / "scripts" / "derive_nsc_rn_charged_twist_ladder.py"),
    ):
        digest = hashlib.sha256()
        with target.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        sources[relative] = digest.hexdigest()
    record["source_sha256"] = sources
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for entry in entries:
        print(
            "  eta={:g} kappa={: .4f}: omega_star={:.6f} star/2={:.6f} "
            "peak omega={:.4f} R={:.6f} peak/star={:.4f}".format(
                entry["eta"], entry["kappa"], entry["omega_star"],
                entry["omega_star_half"], entry["peak"]["omega"],
                entry["peak"]["reflection_ratio"], entry["peak_over_omega_star"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())