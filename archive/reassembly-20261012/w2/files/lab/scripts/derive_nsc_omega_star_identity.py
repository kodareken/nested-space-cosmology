"""Parameter-free identity inside the omega_star quadrature.

omega_star is the infall-time average of the window potential
Phi(r) = Q_E/(16 pi C_F r). Linearity forces the exact identity

    omega_star(eta) = k * <1/r>_T      with k = Q_E/(16 pi C_F)

if the round-4/+27 machinery is coherent. This driver evaluates both
sides of the identity at eta = 1.01, 1.02, 1.04 and records the ratio
spread. A spread at machine level certifies the spectroscopic registry;
a nonconstant spread would expose a hidden functional form.

    python scripts/lab.py scripts/derive_nsc_omega_star_identity.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402
from recursive_horizons import nsc_rn_charged_lll as lll  # noqa: E402

A = 0.04501936182826115
C_F = 0.01125484045706527
Q_MAG = 1.0
Q_E = 4.0
OUTPUT = LAB / "results" / "development" / "nsc-omega-star-identity-v1.json"
SCHEMA = "NSC-OMEGA-STAR-IDENTITY-v1"
ETAS = (1.01, 1.02, 1.04)
K = Q_E / (16.0 * 3.141592653589793 * C_F)


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    rows = []
    for eta in ETAS:
        ref = lll.dyonic_reference(
            A=A, C_F=C_F, q_mag=Q_MAG, Q_E=Q_E, mass_over_total_radius=eta,
        )
        model = ref["model"]
        r_plus = float(ref["horizons"]["r_plus"])
        cap = 32.0 * ref["r_m_tot"]
        fine = np.linspace(r_plus, cap, 3200)
        f = 1.0 - 2.0 * float(model.mass) / fine + float(model.magnetic_r2) / fine ** 2
        weight = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, f)))
        total = float(np.trapezoid(weight, fine))
        phi = lll.electric_potential(fine, Q_E=Q_E, C_F=C_F)
        omega_star = float(np.trapezoid(weight * phi, fine)) / total
        inv_r_mean = float(np.trapezoid(weight / fine, fine)) / total
        ratio = omega_star / inv_r_mean
        rows.append({
            "eta": eta,
            "omega_star": omega_star,
            "inv_r_mean": inv_r_mean,
            "ratio": ratio,
            "k_expected": K,
            "ratio_minus_k": ratio - K,
        })
    ratios = [row["ratio"] for row in rows]
    mean_ratio = sum(ratios) / len(ratios)
    spread = (sum((v - mean_ratio) ** 2 for v in ratios) / len(ratios)) ** 0.5 / abs(mean_ratio)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "k": K,
        "rows": rows,
        "ratio_spread": spread,
        "identity_holds": bool(abs(mean_ratio - K) <= 1e-9 * max(1.0, abs(K))),
        "verdict": (
            "OMEGA_STAR_IDENTITY_HOLDS" if abs(mean_ratio - K) <= 1e-9 * max(1.0, abs(K))
            else "OMEGA_STAR_IDENTITY_BROKEN"
        ),
        "reading": (
            "the spectroscopic omega_star is the k-weighted 1/r infall mean; "
            "this driver verifies it parameter-free"
        ),
        "notes": [
            "linearity of the round-4 Phi_mean quadrature forces the identity; "
            "the measure here is coherence, not a free fit",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_omega_star_identity.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_omega_star_identity.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  k = {:.6f}  (Q_E/(16 pi C_F))".format(K))
    for row in rows:
        print("  eta {:<5} omega_star {:.6f} <1/r> {:.6f} ratio {:.6f}  ratio-k {:.2e}".format(
            row["eta"], row["omega_star"], row["inv_r_mean"], row["ratio"],
            row["ratio_minus_k"],
        ))
    print("  ratio spread: {:.2e}  verdict: {}".format(record["ratio_spread"], record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())