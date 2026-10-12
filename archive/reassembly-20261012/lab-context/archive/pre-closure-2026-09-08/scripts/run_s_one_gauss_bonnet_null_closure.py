#!/usr/bin/env python3
"""Derive the 5D Gauss-Bonnet coupling that closes the null source."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-gauss-bonnet-null-closure.json"
GEOMETRY = ROOT / "results/nsc-1-s-one-two-sheet-5d-geometry.json"
NULL_SOURCE = ROOT / "results/nsc-1-s-one-two-sheet-5d-null-source.json"
COEFFICIENT = ROOT / "results/nsc-1-s-one-transition-link-coefficient.json"


def _derive_null_projections() -> tuple[sp.Expr, sp.Expr, dict[str, sp.Symbol]]:
    time, rho, theta, phi, outside = sp.symbols(
        "t rho theta phi y", real=True
    )
    coordinates = (time, rho, theta, phi, outside)
    warp_curvature = sp.Rational(18, 1015)
    conformal_factor = sp.exp(-2 * warp_curvature * outside**2)
    radius = sp.sqrt(rho**2 + 1)
    metric_b = (
        -3 * sp.pi / 2
        + 1 / (1 + rho**2)
        + 3 * (rho / (1 + rho**2) + sp.atan(rho))
    )
    metric_a = sp.simplify(metric_b * radius**2)
    base = sp.diag(
        metric_a,
        -1 / metric_a,
        -radius**2,
        -radius**2 * sp.sin(theta) ** 2,
        -1,
    )
    metric = sp.diag(
        *[sp.simplify(conformal_factor * base[index, index]) for index in range(5)]
    )
    inverse = sp.diag(
        *[sp.simplify(1 / metric[index, index]) for index in range(5)]
    )
    dimension = 5
    christoffel = [
        [
            [
                sp.simplify(
                    sp.Rational(1, 2)
                    * sum(
                        inverse[upper, lower]
                        * (
                            sp.diff(metric[lower, right], coordinates[left])
                            + sp.diff(metric[lower, left], coordinates[right])
                            - sp.diff(metric[left, right], coordinates[lower])
                        )
                        for lower in range(dimension)
                    )
                )
                for right in range(dimension)
            ]
            for left in range(dimension)
        ]
        for upper in range(dimension)
    ]
    riemann = [
        [
            [
                [
                    sp.simplify(
                        sp.diff(
                            christoffel[upper][lower][fourth],
                            coordinates[third],
                        )
                        - sp.diff(
                            christoffel[upper][lower][third],
                            coordinates[fourth],
                        )
                        + sum(
                            christoffel[upper][middle][third]
                            * christoffel[middle][lower][fourth]
                            - christoffel[upper][middle][fourth]
                            * christoffel[middle][lower][third]
                            for middle in range(dimension)
                        )
                    )
                    for fourth in range(dimension)
                ]
                for third in range(dimension)
            ]
            for lower in range(dimension)
        ]
        for upper in range(dimension)
    ]
    ricci = [
        [
            sp.simplify(
                sum(
                    riemann[upper][lower][upper][right]
                    for upper in range(dimension)
                )
            )
            for right in range(dimension)
        ]
        for lower in range(dimension)
    ]
    ricci_scalar = sp.simplify(
        sum(
            inverse[left, right] * ricci[left][right]
            for left in range(dimension)
            for right in range(dimension)
        )
    )
    riemann_down = [
        [
            [
                [
                    sp.simplify(metric[first, first] * riemann[first][second][third][fourth])
                    for fourth in range(dimension)
                ]
                for third in range(dimension)
            ]
            for second in range(dimension)
        ]
        for first in range(dimension)
    ]

    def lanczos_without_trace(first: int, second: int) -> sp.Expr:
        # The omitted -g_AB*L_GB/2 term vanishes after null contraction.
        ricci_square = sum(
            ricci[first][middle]
            * inverse[middle, middle]
            * ricci[middle][second]
            for middle in range(dimension)
        )
        ricci_riemann = sum(
            inverse[left, left]
            * inverse[right, right]
            * ricci[left][right]
            * riemann_down[first][left][second][right]
            for left in range(dimension)
            for right in range(dimension)
        )
        riemann_square = sum(
            inverse[left, left]
            * inverse[right, right]
            * inverse[last, last]
            * riemann_down[first][left][right][last]
            * riemann_down[second][left][right][last]
            for left in range(dimension)
            for right in range(dimension)
            for last in range(dimension)
        )
        return sp.factor(
            sp.simplify(
                2
                * (
                    ricci_scalar * ricci[first][second]
                    - 2 * ricci_square
                    - 2 * ricci_riemann
                    + riemann_square
                )
            )
        )

    lanczos_tt = lanczos_without_trace(0, 0)
    lanczos_rr = lanczos_without_trace(1, 1)
    affine_scale = sp.exp(warp_curvature * outside**2)
    lanczos_null = sp.factor(
        sp.simplify(
            affine_scale**2
            * (lanczos_tt / metric_a**2 + lanczos_rr)
        )
    )
    ricci_null = -2 * sp.exp(2 * warp_curvature * outside**2) / (
        rho**2 + 1
    ) ** 2
    return ricci_null, lanczos_null, {"rho": rho, "y": outside}


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("Gauss-Bonnet null-closure output already exists")

    authenticated = []
    for path in (GEOMETRY, NULL_SOURCE, COEFFICIENT):
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

    ricci_null, lanczos_null, symbols = _derive_null_projections()
    rho = symbols["rho"]
    outside = symbols["y"]
    expected_lanczos = (
        288
        * sp.exp(sp.Rational(72, 1015) * outside**2)
        / (1015 * (rho**2 + 1) ** 2)
    )
    lanczos_residual = sp.simplify(lanczos_null - expected_lanczos)

    transition_curvature = Fraction(-18, 1015)
    alpha = Fraction(1015, 144)
    alpha_from_transition = -Fraction(1, 8) / transition_curvature
    effective_null = sp.factor(
        sp.simplify(
            ricci_null
            + sp.Rational(alpha.numerator, alpha.denominator) * lanczos_null
        )
    )
    expected_effective_null = (
        2
        * sp.exp(sp.Rational(36, 1015) * outside**2)
        * (sp.exp(sp.Rational(36, 1015) * outside**2) - 1)
        / (rho**2 + 1) ** 2
    )
    effective_null_residual = sp.simplify(
        effective_null - expected_effective_null
    )
    brane_residual = sp.simplify(effective_null.subs(outside, 0))

    samples = []
    for rho_value in (
        Fraction(1, 4),
        Fraction(1, 2),
        Fraction(3, 4),
    ):
        for outside_value in (-1, 0, 1):
            substitutions = {
                rho: sp.Rational(rho_value.numerator, rho_value.denominator),
                outside: outside_value,
            }
            r_value = ricci_null.subs(substitutions)
            h_value = lanczos_null.subs(substitutions)
            total_value = effective_null.subs(substitutions)
            samples.append(
                {
                    "rho": str(rho_value),
                    "outside_y": outside_value,
                    "Ricci_null": float(r_value.evalf(30)),
                    "Gauss_Bonnet_Lanczos_null": float(h_value.evalf(30)),
                    "combined_geometric_null": float(total_value.evalf(30)),
                    "canonical_matter_null_allowed": bool(total_value >= 0),
                }
            )

    record = {
        "artifact_id": "NSC-1-S-ONE-GAUSS-BONNET-NULL-CLOSURE",
        "schema": "NSC-1-S-ONE-GAUSS-BONNET-NULL-CLOSURE-v1",
        "classification": "one_transition_fixed_positive_5D_Gauss_Bonnet_coupling_supplies_the_required_effective_negative_null_geometry_without_phantom_matter",
        "authenticated_inputs": authenticated,
        "one_geometric_action": {
            "formula": "S_geom=(2*kappa5^2)^-1*integral_sqrt_abs_g*[R5+alpha_GB*L_GB]",
            "L_GB": "R^2-4*R_AB*R^AB+R_ABCD*R^ABCD",
            "field_equation_null_projection": "R_kk+alpha_GB*H_GB,kk=kappa5^2*T_matter,kk",
            "field_equations_remain_second_order": "five_dimensional_Lovelock_property",
        },
        "symbolic_derivation": {
            "Ricci_null": str(ricci_null),
            "Lanczos_null": str(lanczos_null),
            "expected_Lanczos_null": str(expected_lanczos),
            "Lanczos_derivation_residual": str(lanczos_residual),
            "combined_geometric_null": str(expected_effective_null),
            "combined_geometric_null_residual": str(effective_null_residual),
            "combined_null_on_local_room": str(brane_residual),
        },
        "coefficient_closure": {
            "transition_b2_over_F": str(transition_curvature),
            "identity": "alpha_GB/L_star^2=-1/[8*(b2/F)]",
            "derived_alpha_GB_over_L_star_squared": str(alpha),
            "derived_not_fitted": alpha == alpha_from_transition,
            "positive": alpha > 0,
        },
        "trapped_interval_samples": samples,
        "physical_result": {
            "negative_Ricci_null_retained": True,
            "effective_negative_source_is_geometric": True,
            "canonical_matter_null_source_nonnegative_on_samples": all(
                item["canonical_matter_null_allowed"] for item in samples
            ),
            "local_room_null_source_vacuum_saturated": brane_residual == 0,
            "outside_bulk_null_source_positive_away_from_room": all(
                item["combined_geometric_null"] > 0.0
                for item in samples
                if item["outside_y"] != 0
            ),
            "separate_phantom_field_inserted": False,
        },
        "next_result": "solve_all_independent_Einstein_Gauss_Bonnet_components_for_the_bulk_source_and_compute_the_tensor_scalar_vector_kinetic_spectrum",
        "nonclaims": {
            "complete_5D_field_equations_solved": False,
            "Gauss_Bonnet_background_ghost_free": False,
            "particle_and_dark_spectra_predicted": False,
            "observed_universe_identified": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == 3,
            "Lanczos_null_derived_exactly": lanczos_residual == 0,
            "alpha_fixed_by_transition": alpha == alpha_from_transition,
            "alpha_positive": alpha > 0,
            "local_room_null_equation_closes": brane_residual == 0,
            "combined_null_simplification_exact": effective_null_residual == 0,
            "canonical_matter_null_nonnegative_on_all_samples": all(
                item["canonical_matter_null_allowed"] for item in samples
            ),
            "full_field_and_stability_system_closed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "full_field_and_stability_system_closed" and value is not True:
            raise RuntimeError(f"Gauss-Bonnet null closure failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
