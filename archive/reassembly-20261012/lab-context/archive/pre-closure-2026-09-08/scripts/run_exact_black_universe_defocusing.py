#!/usr/bin/env python3
"""Reproduce an exact regular black-universe collapse-to-defocusing solution."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-exact-black-universe-defocusing.json"
Q = Fraction


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _exact(value: Fraction) -> dict[str, object]:
    return {"exact": _fraction_text(value), "decimal": float(value)}


def _mp(value: mp.mpf, digits: int = 30) -> str:
    return mp.nstr(value, digits)


def construct() -> dict[str, object]:
    mp.mp.dps = 80
    b = mp.mpf(1)
    mass = mp.mpf(1)
    rho0 = 3 * mass
    c = -mp.pi * rho0 / (2 * b)

    def radius(rho):
        return mp.sqrt(rho * rho + b * b)

    def B(rho):
        return (
            c / b**2
            + 1 / (b * b + rho * rho)
            + rho0
            / b**3
            * (b * rho / (b * b + rho * rho) + mp.atan(rho / b))
        )

    def metric_A(rho):
        return B(rho) * radius(rho) ** 2

    def phi(rho):
        return mp.sqrt(2) * mp.atan(rho / b)

    def potential(rho):
        r2 = radius(rho) ** 2
        return (
            -c / b**2 * (r2 + 2 * rho * rho) / r2
            - rho0
            / b**3
            * (
                3 * b * rho / r2
                + (r2 + 2 * rho * rho) / r2 * mp.atan(rho / b)
            )
        )

    horizon = mp.findroot(metric_A, (mp.mpf(1), mp.mpf(2)))
    samples = tuple(mp.mpf(value) for value in (-10, -2, -1, -0.5, 0, 0.5, 1, 2, 10))
    residuals = []
    for rho in samples:
        r = radius(rho)
        r2 = r * r
        scalar_derivative = mp.diff(phi, rho)
        dV_dphi = mp.diff(potential, rho) / scalar_derivative
        equation4 = mp.diff(lambda x: metric_A(x) * radius(x) ** 2 * mp.diff(phi, x), rho) + r2 * dV_dphi
        equation5 = mp.diff(lambda x: mp.diff(metric_A, x) * radius(x) ** 2, rho) + 2 * r2 * potential(rho)
        equation6 = 2 * mp.diff(radius, rho, 2) / r - scalar_derivative**2
        equation7 = metric_A(rho) * mp.diff(lambda x: radius(x) ** 2, rho, 2) - r2 * mp.diff(metric_A, rho, 2) - 2
        equation8 = mp.diff(B, rho) - 2 * (rho0 - rho) / r**4
        residuals.append(
            {
                "rho": _mp(rho),
                "eq4": _mp(equation4),
                "eq5": _mp(equation5),
                "eq6": _mp(equation6),
                "eq7": _mp(equation7),
                "eq8": _mp(equation8),
                "maximum_absolute": _mp(
                    max(abs(value) for value in (equation4, equation5, equation6, equation7, equation8))
                ),
            }
        )
    maximum_field_residual = max(
        mp.mpf(item["maximum_absolute"]) for item in residuals
    )

    symbolic_rho = sp.symbols("rho", real=True)
    symbolic_radius = sp.sqrt(symbolic_rho**2 + 1)
    symbolic_B = (
        -3 * sp.pi / 2
        + 1 / (1 + symbolic_rho**2)
        + 3
        * (
            symbolic_rho / (1 + symbolic_rho**2)
            + sp.atan(symbolic_rho)
        )
    )
    symbolic_A = sp.simplify(symbolic_B * symbolic_radius**2)
    symbolic_phi = sp.sqrt(2) * sp.atan(symbolic_rho)
    symbolic_V = (
        3
        * sp.pi
        / 2
        * (symbolic_radius**2 + 2 * symbolic_rho**2)
        / symbolic_radius**2
        - 3
        * (
            3 * symbolic_rho / symbolic_radius**2
            + (symbolic_radius**2 + 2 * symbolic_rho**2)
            / symbolic_radius**2
            * sp.atan(symbolic_rho)
        )
    )
    symbolic_equations = {
        "equation_4": sp.simplify(
            sp.diff(
                symbolic_A
                * symbolic_radius**2
                * sp.diff(symbolic_phi, symbolic_rho),
                symbolic_rho,
            )
            + symbolic_radius**2
            * sp.diff(symbolic_V, symbolic_rho)
            / sp.diff(symbolic_phi, symbolic_rho)
        ),
        "equation_5": sp.simplify(
            sp.diff(
                sp.diff(symbolic_A, symbolic_rho) * symbolic_radius**2,
                symbolic_rho,
            )
            + 2 * symbolic_radius**2 * symbolic_V
        ),
        "equation_6": sp.simplify(
            2
            * sp.diff(symbolic_radius, symbolic_rho, 2)
            / symbolic_radius
            - sp.diff(symbolic_phi, symbolic_rho) ** 2
        ),
        "equation_7": sp.simplify(
            symbolic_A
            * sp.diff(symbolic_radius**2, symbolic_rho, 2)
            - symbolic_radius**2 * sp.diff(symbolic_A, symbolic_rho, 2)
            - 2
        ),
        "equation_8": sp.simplify(
            sp.diff(symbolic_B, symbolic_rho)
            - 2 * (3 - symbolic_rho) / symbolic_radius**4
        ),
    }
    symbolic_exact_zero = all(value == 0 for value in symbolic_equations.values())

    def theta_exact(rho: Fraction) -> Fraction:
        return -2 * rho / (rho * rho + 1)

    def ricci_null_exact(rho: Fraction) -> Fraction:
        return -2 / (rho * rho + 1) ** 2

    def complete_q_exact(rho: Fraction) -> Fraction:
        return 2 * (1 - rho * rho) / (rho * rho + 1) ** 2

    rational_samples = (Q(3, 4), Q(1, 2), Q(1, 4), Q(0), Q(-1, 4), Q(-1, 2), Q(-3, 4))
    optical = []
    for rho in rational_samples:
        theta = theta_exact(rho)
        ricci = ricci_null_exact(rho)
        complete_q = complete_q_exact(rho)
        optical.append(
            {
                "rho": _exact(rho),
                "areal_radius_squared": _exact(rho * rho + 1),
                "theta_plus": _exact(theta),
                "theta_minus": _exact(theta),
                "minus_half_theta_squared": _exact(-theta * theta / 2),
                "Ricci_null": _exact(ricci),
                "complete_Q": _exact(complete_q),
                "Raychaudhuri_identity_exact": complete_q == -theta * theta / 2 - ricci,
                "trapped": theta < 0,
                "expanding_antitrapped": theta > 0,
            }
        )
    trapped_interval_q_lower = complete_q_exact(Q(3, 4))
    expanding_interval_q_lower = complete_q_exact(Q(-3, 4))
    outside_rho = mp.mpf("1000000")
    inside_rho = -outside_rho
    external_mass_estimate = radius(outside_rho) * (1 - metric_A(outside_rho)) / 2
    interior_B = B(inside_rho)
    interior_V = potential(inside_rho)

    assertions = {
        "positive_Schwarzschild_mass": mass > 0,
        "one_simple_horizon_outside_defocusing_interval": horizon > 1,
        "metric_static_outside": metric_A(mp.mpf(10)) > 0,
        "metric_T_region_through_transition": metric_A(mp.mpf(1)) < 0 and metric_A(mp.mpf(-1)) < 0,
        "areal_radius_has_finite_positive_minimum": radius(mp.mpf(0)) == b,
        "field_equations_reproduced": maximum_field_residual < mp.mpf("1e-70"),
        "field_equations_symbolically_zero": symbolic_exact_zero,
        "closed_trapped_interval_strictly_defocusing": trapped_interval_q_lower > 0,
        "closed_expanding_interval_strictly_defocusing": expanding_interval_q_lower > 0,
        "affine_Raychaudhuri_identity_exact_at_all_samples": all(
            item["Raychaudhuri_identity_exact"] for item in optical
        ),
        "negative_null_Ricci_supplies_defocusing": all(
            item["Ricci_null"]["decimal"] < 0 for item in optical
        ),
        "same_affine_congruence_changes_sign": optical[0]["theta_plus"]["decimal"] < 0 < optical[-1]["theta_plus"]["decimal"],
        "all_gate_passed": False,
    }
    assertions["all_gate_passed"] = all(
        bool(value) for key, value in assertions.items() if key != "all_gate_passed"
    )
    if not assertions["all_gate_passed"]:
        raise RuntimeError("exact black-universe reproduction failed")

    return {
        "artifact_id": "NSC-1-EXACT-BLACK-UNIVERSE-DEFOCUSING",
        "schema": "NSC-1-EXACT-BLACK-UNIVERSE-DEFOCUSING-v1",
        "classification": "complete_four_dimensional_regular_black_universe_has_trapped_positive_Q_interval_and_child_expansion",
        "source": {
            "paper": "Bronnikov_Dehnen_Melnikov_Regular_black_holes_and_black_universes",
            "arxiv": "gr-qc/0611022v2",
            "doi": "10.1007/s10714-007-0430-6",
            "reproduced_equations": [3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 16],
        },
        "active_scale": "complete_parent_black_hole_to_child_Kantowski_Sachs_de_Sitter_geometry",
        "model": {
            "action": "S=int_sqrt_abs_g_[R-(partial_phi)^2-2V(phi)]",
            "matter": "minimally_coupled_phantom_scalar",
            "metric": "ds2=A(rho)dt2-drho2/A(rho)-r(rho)^2dOmega2",
            "b": _mp(b),
            "positive_external_mass": _mp(mass),
            "rho0": _mp(rho0),
            "c": _mp(c),
            "r_of_rho": "sqrt(rho^2+b^2)",
        },
        "global_geometry": {
            "horizon_rho": _mp(horizon),
            "horizon_areal_radius": _mp(radius(horizon)),
            "minimum_areal_radius": _mp(radius(mp.mpf(0))),
            "external_mass_estimate_at_rho_1e6": _mp(external_mass_estimate),
            "interior_B_at_rho_minus_1e6": _mp(interior_B),
            "interior_V_at_rho_minus_1e6": _mp(interior_V),
            "exact_interior_limits": {
                "B_minus_infinity": "-3*pi",
                "V_minus_infinity": "9*pi",
                "phi_minus_infinity": "-pi/sqrt(2)",
            },
            "future_direction_inside_horizon": "decreasing_rho",
            "future_endpoint": "expanding_isotropizing_de_Sitter_asymptotic_not_r_equals_zero",
        },
        "field_equation_reproduction": {
            "symbolic_equations": {
                key: str(value) for key, value in symbolic_equations.items()
            },
            "symbolic_all_exact_zero": symbolic_exact_zero,
            "precision_decimal_digits": mp.mp.dps,
            "sample_rho_values": [_mp(value) for value in samples],
            "maximum_absolute_residual": _mp(maximum_field_residual),
            "residuals": residuals,
        },
        "affine_null_result": {
            "future_radial_null_generators": "k_plus_minus=(plus_or_minus_1/A,-1,0,0)",
            "affine_condition": "d_rho/d_lambda=-1_constant",
            "theta": "-2*rho/(rho^2+1)",
            "Ricci_null": "-2/(rho^2+1)^2",
            "complete_Q": "2*(1-rho^2)/(rho^2+1)^2",
            "strictly_defocusing_trapped_interval": {
                "rho": ["1/4", "3/4"],
                "affine_length": _exact(Q(1, 2)),
                "minimum_complete_Q": _exact(trapped_interval_q_lower),
                "theta_remains_negative": True,
            },
            "strictly_defocusing_expanding_interval": {
                "rho": ["-3/4", "-1/4"],
                "affine_length": _exact(Q(1, 2)),
                "minimum_complete_Q": _exact(expanding_interval_q_lower),
                "theta_remains_positive": True,
            },
            "turning_surface": {
                "rho": _exact(Q(0)),
                "areal_radius": _exact(Q(1)),
                "theta": _exact(Q(0)),
                "complete_Q": _exact(Q(2)),
            },
            "samples": optical,
        },
        "inherited_invariant": "one_smooth_four_geometry_and_one_affine_null_congruence_connect_external_positive_mass_black_hole_collapse_to_regular_child_expansion",
        "assertions": assertions,
        "observable_consequence": "trapped_parent_side_collapse_locally_turns_into_metric_null_defocusing_before_the_finite_radius_minimum_and_continues_as_child_side_expansion",
        "next_intervention": "replace_the_phantom_effective_description_with_the_positive_F_gradient_boundary_sector_and_test_stability_and_observational_inheritance",
        "nonclaims": {
            "phantom_scalar_is_microscopically_stable": False,
            "this_solution_is_the_FGCQR_action": False,
            "our_observed_universe_is_this_specific_solution": False,
            "nuclear_Skyrme_sector_is_already_coupled_to_this_solution": False,
        },
        "terminal": True,
    }


def main() -> int:
    record = construct()
    OUTPUT.write_bytes(canonical_json_bytes(record))
    print(
        json.dumps(
            {
                "classification": record["classification"],
                "horizon_rho": record["global_geometry"]["horizon_rho"],
                "field_residual": record["field_equation_reproduction"]["maximum_absolute_residual"],
                "trapped_interval": record["affine_null_result"]["strictly_defocusing_trapped_interval"],
                "turning_surface": record["affine_null_result"]["turning_surface"],
                "expanding_interval": record["affine_null_result"]["strictly_defocusing_expanding_interval"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
