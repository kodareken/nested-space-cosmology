"""Does the fringe quantum ladder thin like the awareness profile?

Round-47 counted the fringe delivery {3, 3, 2, 1} quanta at the four
charged cap-fraction seams; round-37 recorded the neutral awareness
beta(s). This driver builds the CHARGED beta at the same four fractions
(lapse fee from round-43, clock from the dyonic quadrature) and lays the
counts against it. No law is asserted in advance: the table and the
spread of count/beta are the measurement.

    python scripts/lab.py scripts/derive_nsc_fringe_beta_join.py
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

SEAM_INVARIANCE = LAB / "results" / "development" / "nsc-charged-seam-invariance-v1.json"
LADDER = LAB / "results" / "development" / "nsc-tilt-quantum-ladder-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-fringe-beta-join-v1.json"
SCHEMA = "NSC-FRINGE-BETA-JOIN-v1"
A = 0.04501936182826115
C_F = 0.01125484045706527
FRACTIONS = (0.125, 0.25, 0.5, 0.75)


def main():
    if not SEAM_INVARIANCE.is_file() or not LADDER.is_file():
        raise SystemExit("missing seam-invariance or ladder record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    invariance = json.loads(SEAM_INVARIANCE.read_text(encoding="utf-8"))
    ladder = json.loads(LADDER.read_text(encoding="utf-8"))
    ref = lll.dyonic_reference(
        A=A, C_F=C_F, q_mag=1.0, Q_E=4.0, mass_over_total_radius=1.01,
    )
    cap = 32.0 * ref["r_m_tot"]
    r_plus = float(ref["horizons"]["r_plus"])
    fine = np.linspace(r_plus, cap, 3200)
    f = 1.0 - 2.0 * float(ref["model"].mass) / fine + float(ref["model"].magnetic_r2) / fine ** 2
    weight = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, f)))
    total = float(np.trapezoid(weight, fine))
    fee_eta_101 = invariance["profiles"]["eta-1.010"]
    rows = []
    for fraction in FRACTIONS:
        seam = fraction * cap
        inside = fine < seam
        t_below = float(np.trapezoid(weight[inside], fine[inside])) / total
        t_above = 1.0 - t_below
        f_seam = float(fee_eta_101[str(fraction)])
        beta = (1.0 - f_seam) / t_above if t_above > 1e-15 else None
        count = ladder["seams"][str(fraction)]["nearest_integer"]
        rows.append({
            "cap_fraction": fraction,
            "seam_rm_tot": seam,
            "lapse_fee": f_seam,
            "t_above": t_above,
            "beta": beta,
            "quantum_count": count,
            "count_over_beta": None if beta is None else count / beta,
        })
    ratios = [row["count_over_beta"] for row in rows if row["count_over_beta"] is not None]
    mean_ratio = sum(ratios) / len(ratios)
    spread = (sum((v - mean_ratio) ** 2 for v in ratios) / len(ratios)) ** 0.5 / abs(mean_ratio)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "count_over_beta_mean": mean_ratio,
        "count_over_beta_spread": spread,
        "reading": (
            "the ladder's per-seam quanta laid against the charged chart's "
            "own beta; no law pre-asserted"
        ),
        "notes": [
            "lapse fees are the round-43 eta=1.01 baseline; the clock is "
            "the dyonic vacuum quadrature at eta=1.01",
        ],
        "input_sha256": {
            "lab/results/development/nsc-charged-seam-invariance-v1.json":
                chain.path_sha256(SEAM_INVARIANCE),
            "lab/results/development/nsc-tilt-quantum-ladder-v1.json":
                chain.path_sha256(LADDER),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_fringe_beta_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_fringe_beta_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  fraction {:<5} beta {:.4f} count {} count/beta {!s}".format(
            row["cap_fraction"], row["beta"], row["quantum_count"],
            None if row["count_over_beta"] is None else "{:.2f}".format(row["count_over_beta"])))
    print("  count/beta mean {:.2f}  spread {:.1%}".format(mean_ratio, spread))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())