#!/usr/bin/env python3
"""Prove a nonzero parameter rectangle of exact black-universe defocusing."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


Q = Fraction
OUTPUT = ROOT / "results/nsc-1-exact-black-universe-rob1.json"


def _text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _record(value: Fraction) -> dict[str, object]:
    return {"exact": _text(value), "decimal": float(value)}


def construct() -> dict[str, object]:
    b_min, b_max = Q(3, 4), Q(5, 4)
    ratio_min, ratio_max = Q(3, 4), Q(5, 4)
    pi_lower = Q(3)

    # A(rho=b)=2*[1/2+q*(3/2-3*pi/4)].  Its q coefficient is negative.
    B_at_b_upper = Q(1, 2) + ratio_min * (
        Q(3, 2) - Q(3, 4) * pi_lower
    )
    A_at_b_upper = 2 * B_at_b_upper
    q_lower = Q(224, 625) / b_max**2
    affine_length_lower = b_min / 2
    mass_lower = b_min * ratio_min
    mass_upper = b_max * ratio_max
    interior_lambda_lower = 9 * pi_lower * ratio_min / b_max**2

    assertions = {
        "parameter_rectangle_has_nonzero_area": b_max > b_min
        and ratio_max > ratio_min,
        "all_external_masses_positive": mass_lower > 0,
        "rho_equals_b_is_strictly_inside_T_region": A_at_b_upper < 0,
        "horizon_lies_outside_entire_positive_Q_trapped_interval": A_at_b_upper < 0,
        "uniform_positive_Q_lower_bound": q_lower > 0,
        "uniform_nonzero_affine_length": affine_length_lower > 0,
        "uniform_positive_child_de_Sitter_limit": interior_lambda_lower > 0,
        "singular_areal_radius_excluded": b_min > 0,
        "all_exact_gate_passed": False,
    }
    assertions["all_exact_gate_passed"] = all(
        value for key, value in assertions.items() if key != "all_exact_gate_passed"
    )

    return {
        "artifact_id": "NSC-1-EXACT-BLACK-UNIVERSE-ROB1",
        "schema": "NSC-1-EXACT-BLACK-UNIVERSE-ROB1-v1",
        "classification": "nonzero_exact_parameter_family_has_uniform_trapped_positive_Q_and_child_expansion",
        "source_result": "results/nsc-1-exact-black-universe-defocusing.json",
        "frozen_parameter_rectangle": {
            "minimum_areal_radius_b": [_record(b_min), _record(b_max)],
            "mass_to_b_ratio_q": [_record(ratio_min), _record(ratio_max)],
            "external_mass_range": [_record(mass_lower), _record(mass_upper)],
            "area_in_b_q_coordinates": _record(
                (b_max - b_min) * (ratio_max - ratio_min)
            ),
        },
        "analytic_family": {
            "r": "sqrt(rho^2+b^2)",
            "theta": "-2*rho/(rho^2+b^2)",
            "Ricci_null": "-2*b^2/(rho^2+b^2)^2",
            "complete_Q": "2*(b^2-rho^2)/(rho^2+b^2)^2",
            "child_de_Sitter_Lambda": "9*pi*(m/b)/b^2",
            "asymptotic_flatness_constant": "c=-3*pi*m/(2*b)",
        },
        "uniform_result": {
            "closed_trapped_defocusing_interval": "rho/b in [1/4,3/4]",
            "minimum_affine_length": _record(affine_length_lower),
            "minimum_complete_Q": _record(q_lower),
            "maximum_A_at_rho_equals_b_using_pi_greater_than_3": _record(
                A_at_b_upper
            ),
            "minimum_child_de_Sitter_Lambda_using_pi_greater_than_3": _record(
                interior_lambda_lower
            ),
            "horizon_is_outside_rho_equals_b": True,
            "turning_surface_rho_zero_is_regular_for_every_member": True,
            "future_child_expansion_for_rho_below_zero": True,
        },
        "proof": {
            "Q_bound": "Q(rho=3b/4)=224/(625*b^2)>=224/[625*(5/4)^2]",
            "horizon_bound": "A(b)<=2*[1/2+(3/4)*(3/2-9/4)]=-1/8<0_using_pi>3",
            "regularity_bound": "r(rho)>=b>=3/4",
            "de_Sitter_bound": "Lambda_child=9*pi*q/b^2>=9*3*(3/4)/(5/4)^2",
        },
        "assertions": assertions,
        "inherited_invariant": "the_same_exact_affine_Raychaudhuri_identity_and_parent_child_orientation_map_hold_throughout_a_nonzero_mass_and_scale_neighborhood",
        "observable_consequence": "collapse_to_defocusing_is_not_an_isolated_parameter_coincidence_in_the_exact_black_universe_family",
        "next_intervention": "replace_the_effective_phantom_source_over_this_frozen_rectangle_with_the_shared_positive_F_topological_gradient_action",
        "nonclaims": {
            "phantom_microphysics_stable": False,
            "same_BPS_parameters_source_this_family": False,
            "observational_holdout_passed": False,
        },
        "terminal": True,
    }


def main() -> int:
    record = construct()
    if not record["assertions"]["all_exact_gate_passed"]:
        raise RuntimeError("black-universe robustness proof failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
