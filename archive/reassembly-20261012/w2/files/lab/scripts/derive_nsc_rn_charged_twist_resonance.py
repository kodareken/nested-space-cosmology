"""Fine resonance scan of the twisted reflection inside the electric window.

Locates the reflection resonance between omega = 0.06 and 0.14 on the
admitted slice (round-5/6 structure), tests the kappa -> -kappa chirality
symmetry at the resonance, and takes one matching point on the eta = 1.5
slice to see how the resonance moves away from extremality.

    python scripts/lab.py scripts/derive_nsc_rn_charged_twist_resonance.py
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
KAPPA = 2.0 ** 0.5
SCAN = tuple(round(0.06 + 0.005 * index, 4) for index in range(17))  # 0.06..0.14
OUTPUT = LAB / "results" / "development" / "nsc-rn-charged-twist-resonance-v1.json"
SCHEMA = "NSC-RN-CHARGED-TWIST-RESONANCE-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))

    def run_slice(eta):
        ref = lll.dyonic_reference(
            A=A, C_F=C_F, q_mag=1.0, Q_E=Q_E, mass_over_total_radius=eta,
        )
        grid = lll.charged_grid(model=ref["model"], r_m_tot=ref["r_m_tot"], points=876)
        chart = lll.charged_chart(ref["model"], grid["radius"])
        return ref, grid, chart

    ref, grid, chart = run_slice(1.01)
    rows = []
    best = None
    for omega in SCAN:
        result = lll.charged_twisted_column(
            model=ref["model"], radius=grid["radius"], omega=omega, sign=SIGN,
            Q_E=Q_E, C_F=C_F, kappa=KAPPA,
        )
        ratio = lll.reflection_ratio(result["mode"], chart=chart)
        rows.append({
            "omega": float(omega),
            "omega_over_Phi_H": float(omega / ref["Phi_H"]),
            "reflection_ratio": ratio["reflection_ratio"],
        })
        if best is None or ratio["reflection_ratio"] > best["reflection_ratio"]:
            best = {
                "omega": float(omega),
                "omega_over_Phi_H": float(omega / ref["Phi_H"]),
                "reflection_ratio": ratio["reflection_ratio"],
                "position": "grid",
            }
    # Parabolic peak location around the grid maximum.
    peaks = [row for row in rows if row["omega"] in
             (best["omega"] - 0.005, best["omega"], best["omega"] + 0.005)]
    if len(peaks) == 3 and peaks[0]["omega"] + 0.005 == peaks[1]["omega"] \
            and peaks[1]["omega"] + 0.005 == peaks[2]["omega"] \
            and peaks[0]["reflection_ratio"] < peaks[2]["reflection_ratio"] * 1.0:
        y0, y1, y2 = (p["reflection_ratio"] for p in peaks)
        denominator = y0 - 2.0 * y1 + y2
        if denominator != 0.0:
            shift = 0.5 * (y0 - y2) / denominator * 0.005
            best["omega_parabolic"] = peaks[1]["omega"] + shift
            best["reflection_ratio_parabolic"] = y1 - 0.25 * (y0 - y2) * shift / 0.005

    # Chirality symmetry at the resonance and one partial point.
    symmetry = {}
    for kappa in (KAPPA, -KAPPA):
        result = lll.charged_twisted_column(
            model=ref["model"], radius=grid["radius"], omega=best["omega"], sign=SIGN,
            Q_E=Q_E, C_F=C_F, kappa=kappa,
        )
        value = lll.reflection_ratio(result["mode"], chart=chart)["reflection_ratio"]
        symmetry["kappa_plus_sqrt2" if kappa > 0 else "kappa_minus_sqrt2"] = {
            "omega": best["omega"],
            "reflection_ratio": value,
        }

    ref15, grid15, chart15 = run_slice(1.5)
    eta15_row = []
    for omega in (0.1, best["omega_over_Phi_H"] * ref15["Phi_H"]):
        result = lll.charged_twisted_column(
            model=ref15["model"], radius=grid15["radius"], omega=omega, sign=SIGN,
            Q_E=Q_E, C_F=C_F, kappa=KAPPA,
        )
        ratio = lll.reflection_ratio(result["mode"], chart=chart15)
        eta15_row.append({
            "omega": float(omega),
            "omega_over_Phi_H_15": float(omega / ref15["Phi_H"]),
            "reflection_ratio": ratio["reflection_ratio"],
            "eta": 1.5,
        })

    record = {
        "schema": SCHEMA,
        "api_version": lll.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "slice": {
            "A": A, "C_F": C_F, "q_mag": 1.0, "Q_E": Q_E, "fermion_sign": SIGN,
            "kappa": KAPPA, "Phi_H": ref["Phi_H"],
            "eta_admitted": 1.01, "eta_far": 1.5,
        },
        "scan": rows,
        "resonance": best,
        "chirality_symmetry": symmetry,
        "symmetry_offset": abs(symmetry["kappa_plus_sqrt2"]["reflection_ratio"]
                              - symmetry["kappa_minus_sqrt2"]["reflection_ratio"]),
        "eta_1_5_points": eta15_row,
        "notes": [
            "kappa and -kappa are the two signs of the same multiplet; "
            "reflection is a flux ratio and must agree to integration error",
            "the eta = 1.5 points probe resonance migration away from "
            "extremality",
        ],
    }
    import hashlib
    sources = {}
    for relative, target in (
        ("lab/src/recursive_horizons/nsc_rn_charged_lll.py",
         LAB / "src" / "recursive_horizons" / "nsc_rn_charged_lll.py"),
        ("lab/scripts/derive_nsc_rn_charged_twist_resonance.py",
         LAB / "scripts" / "derive_nsc_rn_charged_twist_resonance.py"),
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
    for row in rows:
        marker = " <-- grid max" if row["omega"] == best["omega"] else ""
        print("  omega={:.4f}  reflection={:.6e}{}".format(
            row["omega"], row["reflection_ratio"], marker))
    if "omega_parabolic" in best:
        print("  parabolic peak: omega={:.5f} (Phi_H ratio {:.5f}) R={:.6e}".format(
            best["omega_parabolic"], best["omega_over_Phi_H"], best["reflection_ratio_parabolic"]))
    print("  symmetry offset: {:.3e}".format(record["symmetry_offset"]))
    for row in eta15_row:
        print("  eta=1.5 omega={:.4f} (of Phi_H15: {:.4f}) reflection={:.6e}".format(
            row["omega"], row["omega_over_Phi_H_15"], row["reflection_ratio"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())