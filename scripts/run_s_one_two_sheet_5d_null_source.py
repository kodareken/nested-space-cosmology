#!/usr/bin/env python3
"""Derive the null source required by the explicit two-sheet 5D carrier."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-two-sheet-5d-null-source.json"
GEOMETRY = ROOT / "results/nsc-1-s-one-two-sheet-5d-geometry.json"
SPECTRAL = ROOT / "results/nsc-1-s-one-spectral-room-bind.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("two-sheet 5D null-source output already exists")

    authenticated = []
    for path in (GEOMETRY, SPECTRAL):
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item.get("artifact_id"),
            }
        )

    rho, outside = sp.symbols("rho y", real=True)
    warp_log = -sp.Rational(18, 1015) * outside**2
    room_ricci_null = -2 / (1 + rho**2) ** 2
    bulk_ricci_null = sp.simplify(sp.exp(-2 * warp_log) * room_ricci_null)
    sign_removed = sp.simplify(
        bulk_ricci_null / (-2 * sp.exp(sp.Rational(36, 1015) * outside**2))
    )

    samples = []
    for rho_value in (
        Fraction(1, 4),
        Fraction(1, 2),
        Fraction(3, 4),
    ):
        for outside_value in (-1, 0, 1):
            value = bulk_ricci_null.subs(
                {rho: sp.Rational(rho_value.numerator, rho_value.denominator), outside: outside_value}
            )
            samples.append(
                {
                    "rho": str(rho_value),
                    "outside_y": outside_value,
                    "Ricci5_null": str(value),
                    "decimal": float(value.evalf(30)),
                    "strictly_negative": bool(value < 0),
                }
            )

    record = {
        "artifact_id": "NSC-1-S-ONE-TWO-SHEET-5D-NULL-SOURCE",
        "schema": "NSC-1-S-ONE-TWO-SHEET-5D-NULL-SOURCE-v1",
        "classification": "the_explicit_even_5D_carrier_preserves_the_negative_null_curvature_and_requires_a_geometric_spectral_source_beyond_canonical_Einstein_matter",
        "authenticated_inputs": authenticated,
        "conformal_null_identity": {
            "metric": "g5_AB=exp(2W(y))*[g4_mu_nu direct_sum -1]",
            "warp_log_W": str(warp_log),
            "tangent_null_ray": "khat^A=exp(-W)*(k4^mu,0)",
            "reason_cross_terms_vanish": "khat_has_no_y_component_and_g4_mu_nu*k4^mu*k4^nu=0",
            "room_Ricci_null": str(room_ricci_null),
            "bulk_Ricci_null": str(bulk_ricci_null),
            "reduced_exact_factor": str(sign_removed),
        },
        "trapped_interval_samples": samples,
        "source_selection": {
            "pure_5D_Einstein_equation": "R5_AB*k^A*k^B=kappa5^2*Tbulk_AB*k^A*k^B",
            "canonical_minimally_coupled_scalar": "Tkk=(k^A*partial_A*phi)^2>=0",
            "canonical_scalar_alone_can_source_carrier": False,
            "required_one_operator_sector": [
                "full_connection_or_torsion",
                "spectral_Weyl_or_higher_curvature_terms",
                "nonlocal_Schur_complement_from_the_other_sheet",
            ],
            "separate_phantom_field_required": False,
        },
        "known_mathematical_launch_surfaces": [
            {
                "source": "arXiv:1607.07791",
                "result": "Einstein_Cartan_black_universe_with_regular_torsion_and_positive_kinetic_scalar",
                "status": "imported_candidate_connection_reduction_not_yet_the_same_spectral_torsion_operator",
            },
            {
                "source": "arXiv:0911.5074",
                "result": "Standard_Model_spectral_action_with_skew_symmetric_torsion",
                "status": "imported_candidate_spectral_connection_reduction_with_a_different_torsion_class",
            },
        ],
        "next_result": "derive_the_required_null_projection_from_one_admissible_torsion_or_higher_curvature_component_of_mathbb_D_and_check_its_complete_kinetic_spectrum",
        "nonclaims": {
            "Einstein_Cartan_trace_torsion_equals_spectral_skew_torsion": False,
            "bulk_source_action_solved": False,
            "stability_proved": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == 2,
            "conformal_null_factor_exact": sign_removed == 1 / (1 + rho**2) ** 2,
            "all_trapped_samples_negative": all(
                item["strictly_negative"] for item in samples
            ),
            "canonical_scalar_alone_excluded": True,
            "specific_healthy_spectral_source_solved": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "specific_healthy_spectral_source_solved" and value is not True:
            raise RuntimeError(f"two-sheet 5D null source failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
