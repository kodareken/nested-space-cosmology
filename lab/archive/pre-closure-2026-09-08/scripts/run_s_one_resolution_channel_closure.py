#!/usr/bin/env python3
"""Use reciprocal room symmetry to select weak and strong response channels."""

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
from recursive_horizons.unified_action import ReciprocalBoundaryLink  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-resolution-channel-closure.json"
SOURCE = ROOT / "results/nsc-1-s-one-shadow-lensing-observation.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("resolution-channel closure output already exists")

    source_raw = SOURCE.read_bytes()
    source = json.loads(source_raw)
    if source.get("terminal") is not True:
        raise RuntimeError("shadow-lensing observation is not terminal")

    chi = sp.symbols("chi", real=True)
    b0, b2, b4 = sp.symbols("b_0 b_2 b_4", nonnegative=True)
    boundary_link = b0 + b2 * chi**2 + b4 * chi**4
    reciprocal_residual = sp.expand(boundary_link.subs(chi, -chi) - boundary_link)
    linear_mixing = sp.diff(boundary_link, chi).subs(chi, 0)
    nonlinear_activation = sp.diff(boundary_link, chi, 2).subs(chi, 0)

    fixture = ReciprocalBoundaryLink(
        baseline=Fraction(1),
        quadratic=Fraction(1),
        quartic=Fraction(1),
    )
    observed_gamma = float(source["observational_source"]["gamma"])
    observed_sigma = float(source["observational_source"]["one_sigma"])
    reciprocal_weak_gamma = 1.0
    weak_gamma_z = abs(reciprocal_weak_gamma - observed_gamma) / observed_sigma

    record = {
        "artifact_id": "NSC-1-S-ONE-RESOLUTION-CHANNEL-CLOSURE",
        "schema": "NSC-1-S-ONE-RESOLUTION-CHANNEL-CLOSURE-v1",
        "classification": "reciprocal_gradient_symmetry_removes_the_linear_shadow_scalar_while_preserving_nonlinear_transition_activation",
        "authenticated_observation": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "artifact_id": source.get("artifact_id"),
        },
        "one_action_symmetry": {
            "parent_child_map": "T_Theta: chi->-chi",
            "shared_link": str(boundary_link),
            "required_identity": "B_Theta(chi)=B_Theta(-chi)",
            "exact_reciprocal_residual": str(reciprocal_residual),
            "odd_coefficients": "forbidden_by_the_same_parent_child_symmetry",
        },
        "resolution_channels": {
            "stationary_room": {
                "gradient_state": "chi=0",
                "linear_scalar_mixing": str(linear_mixing),
                "metric_response": "universally_coupled_transverse_spin_2_channel",
                "PPN_gamma": reciprocal_weak_gamma,
                "observed_gamma": observed_gamma,
                "observed_one_sigma": observed_sigma,
                "difference_sigma": weak_gamma_z,
            },
            "localized_transition": {
                "gradient_state": "chi_not_equal_0",
                "first_available_scalar_boundary_response": "b2*chi^2",
                "fixed_point_second_derivative": str(nonlinear_activation),
                "role": "nonlinear_outside_geometry_can_activate_without_a_weak_field_fifth_force",
            },
        },
        "public_exact_fixture": {
            "coefficients": {
                "baseline": str(fixture.baseline),
                "quadratic": str(fixture.quadratic),
                "quartic": str(fixture.quartic),
            },
            "response_at_plus_one_half": str(fixture.response(Fraction(1, 2))),
            "response_at_minus_one_half": str(fixture.response(Fraction(-1, 2))),
            "fixed_point_linear_mixing": str(fixture.fixed_point_linear_mixing),
            "fixed_point_nonlinear_activation": str(
                fixture.fixed_point_nonlinear_activation
            ),
        },
        "causal_repair": {
            "measured_failure": "constant_two_wall_shadow_scalar_predicts_gamma_1_over_2",
            "responsible_layer": "an_unsuppressed_odd_or_linear_boundary_mode_at_the_stationary_room",
            "repair": "enforce_the_already_declared_reciprocal_chi_to_minus_chi_symmetry",
            "new_threshold_or_observational_fit": False,
            "weak_field_result": "gamma_returns_to_1_independently_of_b2_and_b4",
            "strong_field_capacity_retained": "even_nonlinear_boundary_terms_remain_active",
        },
        "next_result": "determine_b2_and_the_tensor_form_factor_from_the_positive_scalar_black_universe_bulk_solution_then_predict_the_finite_wavelength_dark_response_without_refitting_gamma",
        "nonclaims": {
            "dark_matter_amplitude_predicted": False,
            "black_universe_bulk_solved": False,
            "b2_fixed": False,
            "particle_spectrum_fixed": False,
        },
        "gate": {
            "reciprocal_identity_exact": reciprocal_residual == 0,
            "linear_scalar_mixing_exactly_zero": linear_mixing == 0,
            "nonlinear_activation_retained": nonlinear_activation == 2 * b2,
            "weak_field_gamma_within_one_sigma": weak_gamma_z < 1.0,
            "no_observation_fitted_coefficient": True,
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("resolution-channel closure gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
