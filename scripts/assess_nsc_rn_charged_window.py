"""Decide whether the charged electric window fits the declared parent cap.

The [charged source](../docs/nsc-rn-charged-source.md) note declares the
dyonic window B_< and the parent cap 32 r_m,tot. Before any mode solver or
march is written for that branch, this driver asks the one question the
declaration must answer: does B_< contain a frequency the cap can carry?

No time step is taken. No mode is built. The scan is over a declared grid of
black-hole Gauss charges and extremality ratios at the audited coefficients,
for both fermion signs, with the pure-magnetic arm as the control.

    python scripts/lab.py scripts/assess_nsc_rn_charged_window.py \\
        --record results/development/nsc-rn-charged-window-assessment-v1.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_rn_charged_window as window  # noqa: E402
from recursive_horizons import nsc_rn_reference as reference  # noqa: E402

RECORD = LAB / "results/development/nsc-rn-charged-window-assessment-v1.json"

#: Declared scan, not a fit. Q_E = 0 is the pure-magnetic control; the
#: extremality ratios bracket a safely subextremal hole (ratio 2) down to a
#: nearly extremal one (ratio 1.01).
GAUSS_CHARGES = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)
EXTREMALITY_RATIOS = (1.01, 1.1, 1.5, 2.0, 4.0)
SIGNS = (1, -1)
FLUX = 1.0


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _producing_commit() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _row(Q_E: float, eta: float, sign: int) -> Dict[str, Any]:
    result = window.window_gate(
        A=reference.AUDITED_A, C_F=reference.AUDITED_C_F, q_mag=FLUX,
        Q_E=float(Q_E), mass_over_total_radius=float(eta), sign=sign,
    )
    return {
        "Q_E": float(Q_E),
        "mass_over_total_radius": float(eta),
        "fermion_sign": int(sign),
        "arm": "pure_magnetic_control" if Q_E == 0.0 else "dyonic",
        "admitted": bool(result["admitted"]),
        "window_open": bool(result["window_open"]),
        "reason": result["reason"],
        "Phi_H": result["Phi_H"],
        "omega_min": result["omega_min"],
        "window_upper": result["window_upper"],
        "omega_central": result["omega_central"],
        "band_width": result["band_width"],
        "r_plus": result["r_plus"],
        "r_minus": result["r_minus"],
        "cap_room": result["cap_room"],
        "total_charge_radius": result["charge_squares"]["total_charge_radius"],
        "magnetic_radius": result["charge_squares"]["magnetic_radius"],
        "P_total_squared": result["charge_squares"]["P_total_squared"],
        "neutral_frequency": result["neutral_frequency"],
        "neutral_frequency_band": result["neutral_frequency_band"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", default=str(RECORD.relative_to(LAB)))
    args = parser.parse_args()
    record_path = (LAB / args.record).resolve()
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    rows: List[Dict[str, Any]] = [
        _row(Q_E, eta, sign)
        for Q_E in GAUSS_CHARGES
        for eta in EXTREMALITY_RATIOS
        for sign in SIGNS
    ]
    dyonic = [r for r in rows if r["arm"] == "dyonic"]
    admitted = [r for r in dyonic if r["admitted"]]
    controls = [r for r in rows if r["arm"] == "pure_magnetic_control"]
    controls_closed = all(not r["window_open"] for r in controls)

    best = None
    if admitted:
        best = max(admitted, key=lambda r: (r["band_width"], -r["omega_min"]))

    record = {
        "schema": "NSC-RN-CHARGED-WINDOW-v1",
        "api_version": window.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "question": (
            "does the electric window B_< contain a frequency the declared "
            "parent cap 32 r_m,tot can carry?"
        ),
        "coefficients": {
            "A": reference.AUDITED_A,
            "C_F": reference.AUDITED_C_F,
            "flux": reference.AUDITED_FLUX,
            "scan_flux": FLUX,
            "V4": reference.AUDITED_V4,
            "source": "nsc_rn_reference audited coefficients; flux scanned at 1",
        },
        "declaration": {
            "parent_cap_multiple": window.PARENT_CAP_MULTIPLE,
            "budget_numerator": window.BUDGET_NUMERATOR,
            "gate": "admitted iff omega_min < s Phi_H, with omega_min = 4/(32 r_m,tot - r_+)",
            "scan_gauss_charges": list(GAUSS_CHARGES),
            "scan_extremality_ratios": list(EXTREMALITY_RATIOS),
            "fermion_signs": list(SIGNS),
        },
        "sources": [
            {"path": "lab/src/recursive_horizons/nsc_rn_charged_window.py",
             "sha256": _sha256_file(LAB / "src/recursive_horizons/nsc_rn_charged_window.py")},
            {"path": "lab/src/recursive_horizons/nsc_rn_reference.py",
             "sha256": _sha256_file(LAB / "src/recursive_horizons/nsc_rn_reference.py")},
            {"path": "lab/docs/nsc-rn-charged-source.md",
             "sha256": _sha256_file(LAB / "docs/nsc-rn-charged-source.md")},
        ],
        "limits": [
            "source-free algebra on a declared slice: no mode, no packet, no march",
            "the gate is the note's own declared cap; it is not a statement about the actual universe",
            "admitted means the declaration is self-consistent, not that a response exists",
        ],
        "coverage": {"rows": len(rows), "dyonic_rows": len(dyonic),
                     "admitted_rows": len(admitted), "control_rows": len(controls)},
        "controls": {
            "pure_magnetic_window_closed": bool(controls_closed),
            "statement": "Q_E = 0 gives Phi_H = 0, so the control cannot enter B_<",
        },
        "verdict": {
            "gate": "open" if admitted else "closed",
            "best_admitted": best,
            "reading": (
                "the dyonic declaration admits a packet inside the cap"
                if admitted else
                "no declared dyonic slice admits a packet inside the cap"
            ),
        },
        "rows": rows,
    }
    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  rows {len(rows)}  dyonic {len(dyonic)}  admitted {len(admitted)}")
    print(f"  pure-magnetic window closed: {controls_closed}")
    print(f"  gate: {record['verdict']['gate']}")
    if best is not None:
        print(
            f"  best: Q_E={best['Q_E']:g} eta={best['mass_over_total_radius']:g} "
            f"sign={best['fermion_sign']:+d} omega_min={best['omega_min']:.6g} "
            f"Phi_H={best['Phi_H']:.6g} width={best['band_width']:.6g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
