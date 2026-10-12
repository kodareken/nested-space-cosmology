#!/usr/bin/env python3
"""Regenerate the GMF-1B-ECD-INT1 interaction-identity artifact."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_ecd_identity import (  # noqa: E402
    contact_tetrad_stress,
    metric_null_contact_contraction,
    required_nonclaims,
    vc_action_projection_identity,
    vc_contact_action_sources,
    vc_reduced_contact_density,
    vc_signature_bridge,
)
from recursive_horizons.gmf1b_ecd_symmetry import ventrella_chiral_axial_bilinears  # noqa: E402


def _json(value: object) -> object:
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    raise TypeError(f"cannot serialize {type(value).__name__}")


def record() -> dict[str, object]:
    radius = 1.0 / (2.0 * 3.141592653589793**0.5)
    F, G = 0.75 + 0j, 0.5 + 0j
    kappa = 1.0
    metric = ((-1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0), (0.0, 0.0, 0.0, 1.0))
    axial = ventrella_chiral_axial_bilinears(F, G, radius, 1.0)
    one_amplitude_axial = ventrella_chiral_axial_bilinears(1.0 + 0j, 0.0 + 0j, radius, 1.0)
    return {
        "schema_version": 2,
        "project_version": "0.11.0",
        "model_id": "GMF-1B-ECD-INT1",
        "artifact": "ecd_interaction_action_gamma_identity_gate",
        "classification": "bounded_ecd_contact_identity_and_tested_local_reconstruction_preflight",
        "exact_fixture": {
            "F": F,
            "G": G,
            "areal_radius": radius,
            "radial_metric": 1.0,
            "lapse": 1.0,
            "kappa": kappa,
            "axial": axial,
            "signature_bridge": vc_signature_bridge(F, G, 1.0, radius, kappa, 0.7, 0.2),
            "reduced_contact_density": vc_reduced_contact_density(F, G, 1.0, 1.0, radius, kappa),
            "action_sources": vc_contact_action_sources(F, G, 1.0, 1.0, radius, kappa),
            "identity": vc_action_projection_identity(F, G, 1.0, 1.0, radius, kappa, ((0.7, 0.2), (1.1, -0.9), (2.0, 1.4))),
        },
        "contact_tetrad_stress": contact_tetrad_stress(
            kappa, axial["axial_squared"], metric
        ),
        "null_contact_gate": metric_null_contact_contraction(
            kappa, axial["axial_squared"], metric, (1.0, 1.0, 0.0, 0.0)
        ),
        "one_amplitude_semantics": {
            "F": {"real": 1.0, "imag": 0.0},
            "G": {"real": 0.0, "imag": 0.0},
            "axial": one_amplitude_axial,
            "reduced_contact_density": vc_reduced_contact_density(1.0, 0.0, 1.0, 1.0, radius, kappa),
            "interpretation": "A nonzero null axial/torsion current can have A^2=0 and a vanishing projected cubic contact source.",
        },
        "nonclaims": required_nonclaims(),
        "provenance": {
            "cabral_lobo_rubiera_garcia_2019": {
                "arxiv_url": "https://arxiv.org/abs/1902.02222",
                "equations": "8-11, 18, 20, 23, and 83",
                "use": "combined torsion-eliminated contact action and Hehl-Datta equation in the canonical layer, then explicit signature translation",
            }
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPOSITORY / "results" / "gmf-1b-ecd-identity.json")
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
