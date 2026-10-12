"""Is the charged spectroscopy clocked by the same sqrt(gap) law?

Round-26 found the neutral wall's inflation phase obeys rise ~ Delta^-1/2
(family spread collapsed 34.5% -> 4.0%). This driver asks whether the
charged programme's phase-averaged window frequency omega_star follows
the same horizon-gap law across the mass ratio eta. omega_star is the
round-4 integral (Phi_mean of the electric window weighted by infall
time), recomputed for eta = 1.01, 1.02, 1.04 on the dyonic chart; the
gap is 2 sqrt(M^2 - P_tot^2) of the same chart. If the spectroscopy is
gap-clocked, omega_star * sqrt(Delta) is family-constant.

    python scripts/lab.py scripts/derive_nsc_charged_clock_horizon.py
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
OUTPUT = LAB / "results" / "development" / "nsc-charged-clock-horizon-v1.json"
SCHEMA = "NSC-CHARGED-CLOCK-HORIZON-v1"
ETAS = (1.01, 1.02, 1.04)


def phi_mean(model, r_m_tot, r_plus_h):
    cap = 32.0 * r_m_tot
    fine = np.linspace(r_plus_h, cap, 3200)
    f = 1.0 - 2.0 * float(model.mass) / fine + float(model.magnetic_r2) / fine ** 2
    weight = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, f)))
    phi = lll.electric_potential(fine, Q_E=Q_E, C_F=C_F)
    total = float(np.trapezoid(weight, fine))
    return float(np.trapezoid(weight * phi, fine)) / total, total


def spread(values):
    mean = sum(values) / len(values)
    return (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5 / abs(mean)


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
        r_minus = float(ref["horizons"]["r_minus"])
        gap = r_plus - r_minus
        omega_star, t_total = phi_mean(model, ref["r_m_tot"], r_plus)
        rows.append({
            "eta": eta,
            "r_plus": r_plus,
            "r_minus": r_minus,
            "horizon_gap": gap,
            "omega_star": omega_star,
            "T_total": t_total,
            "omega_star_times_sqrt_gap": omega_star * (gap ** 0.5),
        })
    omegas = [row["omega_star"] for row in rows]
    gapped = [row["omega_star_times_sqrt_gap"] for row in rows]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "omega_star_raw_spread": spread(omegas),
        "omega_star_sqrt_gap_spread": spread(gapped),
        "reading": (
            "if the charged spectroscopy is gap-clocked like the wall rise "
            "phase, the gapped spread collapses below the raw spread"
        ),
        "notes": [
            "omega_star = phase-averaged window potential at sign +1, the "
            "round-4 definition with the same quadrature",
            "the dyonic gap is 2 sqrt(M^2 - P_tot^2) with fixed A, C_F, "
            "q_mag=1, Q_E=4 as eta varies",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_charged_clock_horizon.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_charged_clock_horizon.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  eta {:<5} gap {:.4f} omega_star {:.4f} omega*sqrt {: .4f}".format(
            row["eta"], row["horizon_gap"], row["omega_star"],
            row["omega_star_times_sqrt_gap"],
        ))
    print("  spreads: raw {:.1%} -> gapped {:.1%}".format(
        record["omega_star_raw_spread"], record["omega_star_sqrt_gap_spread"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())