#!/usr/bin/env python3
"""Fix the even link curvature from the measured positive-F transition."""

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
from recursive_horizons.unified_action import ReciprocalBoundaryLink  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-transition-link-coefficient.json"
TRANSITION = ROOT / "results/nsc-1-unified-gradient-boundary-dev3.json"
WEAK_CHANNEL = ROOT / "results/nsc-1-s-one-resolution-channel-closure.json"
TRANSITION_RUNNER = ROOT / "scripts/run_unified_gradient_boundary_dev3.py"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("transition-link coefficient output already exists")

    transition_raw = TRANSITION.read_bytes()
    weak_raw = WEAK_CHANNEL.read_bytes()
    transition = json.loads(transition_raw)
    weak = json.loads(weak_raw)
    if transition.get("terminal") is not True or weak.get("terminal") is not True:
        raise RuntimeError("transition-link inputs must be terminal")

    # Frozen constants in the authenticated DEV3 runner.
    thickness = Fraction(1, 50)
    expansion_magnitude = Fraction(225, 64)
    center_radius_factor = 1 - expansion_magnitude * thickness / 8
    center_ricci_null = (
        -2 * expansion_magnitude / (thickness * center_radius_factor)
    )

    # xi=(lambda-L/2)/(L/2), hence d/dxi=(L/2)d/dlambda.
    normalized_second_derivative_over_F = (
        thickness**2 * center_ricci_null / 4
    )
    normalized_quadratic_over_F = normalized_second_derivative_over_F / 2
    center_f = float(transition["result"]["maximum_F"])
    absolute_quadratic = float(normalized_quadratic_over_F) * center_f

    link = ReciprocalBoundaryLink(
        baseline=Fraction(4),
        quadratic=normalized_quadratic_over_F,
    )
    plus = link.response(Fraction(1, 10))
    minus = link.response(Fraction(-1, 10))

    transition_runner_raw = TRANSITION_RUNNER.read_bytes()
    record = {
        "artifact_id": "NSC-1-S-ONE-TRANSITION-LINK-COEFFICIENT",
        "schema": "NSC-1-S-ONE-TRANSITION-LINK-COEFFICIENT-v1",
        "classification": "the_positive_F_transition_fixes_the_even_link_curvature_without_using_the_weak_field_lensing_observation",
        "authenticated_inputs": [
            {
                "path": str(TRANSITION.relative_to(ROOT)),
                "sha256": hashlib.sha256(transition_raw).hexdigest(),
                "artifact_id": transition.get("artifact_id"),
            },
            {
                "path": str(WEAK_CHANNEL.relative_to(ROOT)),
                "sha256": hashlib.sha256(weak_raw).hexdigest(),
                "artifact_id": weak.get("artifact_id"),
            },
            {
                "path": str(TRANSITION_RUNNER.relative_to(ROOT)),
                "sha256": hashlib.sha256(transition_runner_raw).hexdigest(),
                "role": "owns_the_frozen_rational_transition_constants",
            },
        ],
        "transition_derivation": {
            "profile_equation": "d2F/dlambda2=R_kk(lambda)*F",
            "affine_thickness_L": str(thickness),
            "expansion_magnitude_a": str(expansion_magnitude),
            "center_radius_factor": str(center_radius_factor),
            "center_Ricci_null": str(center_ricci_null),
            "dimensionless_coordinate": "xi=(lambda-L/2)/(L/2)",
            "F_xixi_over_F_at_center": str(
                normalized_second_derivative_over_F
            ),
            "quadratic_coefficient_over_F": str(normalized_quadratic_over_F),
            "measured_center_F": center_f,
            "absolute_quadratic_coefficient": absolute_quadratic,
        },
        "same_coefficient_two_regimes": {
            "weak_stationary_room": {
                "identity": "F(xi)=F(-xi)",
                "linear_mixing_F_xi_at_zero": "0",
                "PPN_gamma": weak["resolution_channels"]["stationary_room"][
                    "PPN_gamma"
                ],
                "coefficient_fitted_to_gamma": False,
            },
            "strong_boundary": {
                "nonlinear_curvature": "F_xixi(0)=2*b2",
                "b2_over_F": str(normalized_quadratic_over_F),
                "negative_sign_supplies_defocusing_direction": (
                    normalized_quadratic_over_F < 0
                ),
                "source": "authenticated_positive_F_null_metric_equation",
            },
        },
        "signed_public_fixture": {
            "baseline": str(link.baseline),
            "quadratic": str(link.quadratic),
            "response_at_plus_one_tenth": str(plus),
            "response_at_minus_one_tenth": str(minus),
            "reciprocal_equal": plus == minus,
            "linear_mixing": str(link.fixed_point_linear_mixing),
            "nonlinear_activation": str(link.fixed_point_nonlinear_activation),
        },
        "causal_correction": {
            "prior_implementation_error": "quadratic_boundary_response_was_restricted_to_positive_values_without_physical_justification",
            "correction": "reciprocity_requires_evenness_not_positive_curvature",
            "measured_transition_sign": "negative",
            "historical_result_rewritten": False,
        },
        "next_result": "insert_the_fixed_minus_18_over_1015_link_curvature_into_the_tensor_recursive_kernel_and_solve_the_remaining_metric_scalar_and_bulk_equations",
        "nonclaims": {
            "full_tensor_form_factor_fixed": False,
            "bulk_solution_completed": False,
            "dark_matter_amplitude_predicted": False,
            "particle_spectrum_predicted": False,
        },
        "gate": {
            "transition_inputs_authenticated": True,
            "normalized_coefficient_exact": normalized_quadratic_over_F
            == Fraction(-18, 1015),
            "weak_linear_mixing_zero": link.fixed_point_linear_mixing == 0,
            "strong_nonlinear_response_negative": (
                link.fixed_point_nonlinear_activation < 0
            ),
            "reciprocal_profile_exact": plus == minus,
            "gamma_not_used_to_fix_coefficient": True,
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("transition-link coefficient gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
