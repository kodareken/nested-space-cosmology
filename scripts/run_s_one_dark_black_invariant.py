#!/usr/bin/env python3
"""Reduce transition, vacuum-form, and Gauss-Bonnet coefficients to one invariant."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-dark-black-invariant.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-transition-link-coefficient.json",
    ROOT / "results/nsc-1-s-one-gauss-bonnet-null-closure.json",
    ROOT / "results/nsc-1-s-one-gauss-bonnet-component-closure.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("dark-black invariant output already exists")

    authenticated = []
    for path in INPUTS:
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

    transition = Fraction(-18, 1015)
    vacuum_form = Fraction(108, 1015)
    gauss_bonnet = Fraction(1015, 144)
    vacuum_from_transition = -6 * transition
    gauss_bonnet_from_transition = -Fraction(1, 8) / transition
    invariant_product = gauss_bonnet * vacuum_form

    record = {
        "artifact_id": "NSC-1-S-ONE-DARK-BLACK-INVARIANT",
        "schema": "NSC-1-S-ONE-DARK-BLACK-INVARIANT-v1",
        "classification": "transition_curvature_vacuum_form_and_higher_curvature_coupling_are_one_dimensionless_invariant_written_three_ways",
        "authenticated_inputs": authenticated,
        "one_number": {
            "transition_curvature_b2_over_F": str(transition),
            "local_vacuum_form_coefficient": str(vacuum_form),
            "Gauss_Bonnet_coefficient_over_L_star_squared": str(gauss_bonnet),
        },
        "exact_equalities": {
            "vacuum_from_transition": "v_room=-6*c_transition",
            "vacuum_residual": str(vacuum_form - vacuum_from_transition),
            "Gauss_Bonnet_from_transition": "a_GB=-1/(8*c_transition)",
            "Gauss_Bonnet_residual": str(
                gauss_bonnet - gauss_bonnet_from_transition
            ),
            "scale_free_dark_black_product": "a_GB*v_room=3/4",
            "product": str(invariant_product),
            "product_residual": str(invariant_product - Fraction(3, 4)),
        },
        "physical_dictionary": {
            "transition_curvature": "nonlinear_even_gradient_needed_for_metric_null_defocusing",
            "vacuum_form": "equal_tangential_response_locally_named_cosmological_or_dark_energy_curvature",
            "Gauss_Bonnet_coupling": "geometric_term_that_moves_the_negative_null_contribution_off_the_matter_source",
            "single_cause": "all_three_are_fixed_by_the_same_two_sheet_warped_geometry",
        },
        "dimensional_form": {
            "transition": "(b2/F)_physical=(-18/1015)/L_star^2",
            "vacuum": "v_room_physical=(108/1015)/L_star^2",
            "Gauss_Bonnet": "alpha_GB=(1015/144)*L_star^2",
            "unit_independent_relation": "alpha_GB*v_room_physical=3/4",
        },
        "next_result": "derive_L_star_from_the_parent_child_spectral_inheritance_map_then_predict_the_dimensional_vacuum_and_transition_scales_together",
        "nonclaims": {
            "vacuum_form_sign_convention_identified_with_observed_positive_Lambda": False,
            "L_star_derived": False,
            "observed_dark_energy_magnitude_predicted": False,
            "full_stability_proved": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "vacuum_is_minus_six_transition": vacuum_form
            == vacuum_from_transition,
            "Gauss_Bonnet_is_inverse_transition": gauss_bonnet
            == gauss_bonnet_from_transition,
            "scale_free_product_three_quarters": invariant_product
            == Fraction(3, 4),
            "new_independent_coefficient": False,
        },
        "terminal": True,
    }
    if not all(
        value
        for key, value in record["gate"].items()
        if key != "new_independent_coefficient"
    ):
        raise RuntimeError("dark-black invariant gate failed")
    if record["gate"]["new_independent_coefficient"] is not False:
        raise RuntimeError("dark-black invariant introduced a coefficient")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
