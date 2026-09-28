#!/usr/bin/env python3
"""Regenerate the GMF-1B-ECD-TOR1 spin-tensor/torsion/contact-stress artifact."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_ecd_torsion import (  # noqa: E402
    contact_tetrad_stress_interpretation,
    effective_contact_lagrangian_crosscheck,
    required_nonclaims,
    spin_tensor_components,
    torsion_quadratic_invariant,
)


def _json(value: object) -> object:
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    raise TypeError(f"cannot serialize {type(value).__name__}")


def record() -> dict[str, object]:
    radius = 1.0 / (2.0 * 3.141592653589793**0.5)
    F, G = 0.75 + 0j, 0.5 + 0j
    kappa = 1.0
    theta, phi = 0.7, 0.2
    metric = (
        (-1.0, 0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    )
    spin = spin_tensor_components(F, G, radius, 1.0, theta, phi)
    torsion = torsion_quadratic_invariant(F, G, radius, 1.0, kappa, theta, phi)
    contact = effective_contact_lagrangian_crosscheck(
        F, G, radius, 1.0, kappa, theta, phi
    )
    stress = contact_tetrad_stress_interpretation(kappa, spin["axial_squared"], metric)
    return {
        "schema_version": 2,
        "project_version": "0.11.0",
        "model_id": "GMF-1B-ECD-TOR1",
        "artifact": "ecd_spin_tensor_torsion_and_effective_contact_crosscheck",
        "classification": (
            "spin_tensor_torsion_invariant_and_reduced_contact_tetrad_stress_gate"
        ),
        "conventions": {
            "metric_signature": "-+++",
            "orientation_0123": 1,
            "kappa": kappa,
        },
        "fixture": {
            "F": F,
            "G": G,
            "areal_radius": radius,
            "radial_metric": 1.0,
            "theta": theta,
            "phi": phi,
        },
        "axial_current": {
            "A_hat_0": spin["A_hat_0"],
            "A_hat_r": spin["A_hat_r"],
            "axial_squared": spin["axial_squared"],
        },
        "spin_tensor": {
            "hodge_dual_of_lowered_axial_current": True,
            "components": spin["spin_tensor_components"],
            "hodge_dual_residual": spin["hodge_dual_residual"],
            "hodge_dual_verified": spin["hodge_dual_verified"],
            "totally_antisymmetric": spin["totally_antisymmetric"],
        },
        "torsion": torsion,
        "contact_derivation": contact,
        "contact_tetrad_stress": stress,
        "nonclaims": required_nonclaims(),
        "provenance": {
            "full_torsion_eliminated_contact": {
                "reference": (
                    "Choudhury-Maity-Lahiri 2024 Eq. 7 and Cabral-Lobo-"
                    "Rubiera-Garcia 2019 Eqs. 20-24"
                ),
                "use": (
                    "cross-check the full reduced contact after the independent "
                    "connection is eliminated everywhere; do not double count "
                    "intermediate gravitational and Dirac connection pieces"
                ),
            }
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=REPOSITORY / "results" / "gmf-1b-ecd-torsion.json"
    )
    output = parser.parse_args().output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(record(), default=_json, indent=2, sort_keys=True, allow_nan=False) + "\n"
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
