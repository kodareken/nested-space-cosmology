#!/usr/bin/env python3
"""Regenerate the corrected GMF-1B-ECD-KIN1 tetrad-variation artifact."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_ecd_kinetic import (  # noqa: E402
    anticommutator_identity,
    contact_stress_null_energy_gate,
    effective_action_bookkeeping,
    free_dirac_positive_frequency_plane_wave_control,
    required_nonclaims,
    tetrad_current_invariant_cancellation,
)


def _json(value: object) -> object:
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    raise TypeError(f"cannot serialize {type(value).__name__}")


def record() -> dict[str, object]:
    radius = 1.0 / (2.0 * 3.141592653589793**0.5)
    F, G = 0.75 + 0j, 0.5 + 0j
    kappa = 1.0
    metric = (
        (-1.0, 0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    )
    return {
        "schema_version": 2,
        "project_version": "0.11.0",
        "model_id": "GMF-1B-ECD-KIN1",
        "artifact": "ecd_contact_tetrad_variation_and_free_dirac_control",
        "classification": (
            "corrected_contact_stress_and_bounded_torsion_free_dirac_preflight"
        ),
        "conventions": {
            "canonical_signature": "+---",
            "vc_signature": "-+++",
            "gamma_conversion": "Gamma^a=-i*gamma_VC^a",
            "orientation_0123": 1,
            "kappa": kappa,
        },
        "fixture": {
            "F": F,
            "G": G,
            "areal_radius": radius,
            "radial_metric": 1.0,
        },
        "anticommutator_identity": anticommutator_identity(),
        "effective_action_bookkeeping": effective_action_bookkeeping(
            F, G, radius, 1.0, kappa
        ),
        "tetrad_current_invariant_cancellation": tetrad_current_invariant_cancellation(),
        "contact_stress_null_energy": contact_stress_null_energy_gate(
            kappa,
            -9.0 / 4.0,
            metric,
            ((1.0, 1.0, 0.0, 0.0), (1.0, 0.0, 0.0, 1.0)),
        ),
        "free_dirac_plane_wave_control": free_dirac_positive_frequency_plane_wave_control(),
        "correction": {
            "retracted_claim": (
                "A fixed-coordinate-current variation of J^2 supplies a physical "
                "anisotropic contact stress and local NEC violation."
            ),
            "reason": (
                "J_mu=e^I_mu J_I varies with the coframe; its response cancels "
                "the fixed-current metric partial because J^I J_I is tetrad independent."
            ),
            "correct_contact_result": "T_ab=L_4 g_ab and T_ab k^a k^b=0",
            "homogeneous_spin_fluid_bounce_not_derived": True,
        },
        "nonclaims": required_nonclaims(),
        "provenance": {
            "torsion_eliminated_action_and_stress": {
                "reference": (
                    "Choudhury-Maity-Lahiri 2024 Eqs. 7-10; "
                    "Cabral-Lobo-Rubiera-Garcia 2019 Eqs. 20-24"
                ),
                "use": (
                    "full reduced contact sign, tetrad stress, and separation "
                    "from the torsion-free Dirac kinetic stress"
                ),
            },
            "free_dirac_scope": {
                "reference": "direct canonical gamma-matrix plane-wave fixture",
                "use": "positive-frequency control only; no general Dirac NEC theorem",
            },
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=REPOSITORY / "results" / "gmf-1b-ecd-kinetic.json"
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
