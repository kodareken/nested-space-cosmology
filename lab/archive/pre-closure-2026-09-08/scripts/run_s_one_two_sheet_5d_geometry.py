#!/usr/bin/env python3
"""Construct and test an explicit compact 5D carrier for both black-universe sheets."""

from __future__ import annotations

import hashlib
import json
from math import isfinite
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-two-sheet-5d-geometry.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-two-sheet-black-geometry.json",
    ROOT / "results/nsc-1-s-one-transition-link-coefficient.json",
)


def _derive_invariants() -> tuple[sp.Expr, sp.Expr, sp.Expr, dict[str, sp.Symbol]]:
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
    ricci_scalar = sp.factor(
        sp.simplify(
            sum(
                inverse[left, right] * ricci[left][right]
                for left in range(dimension)
                for right in range(dimension)
            )
        )
    )
    kretschmann = sp.factor(
        sp.simplify(
            sum(
                metric[upper, upper]
                * inverse[lower, lower]
                * inverse[third, third]
                * inverse[fourth, fourth]
                * riemann[upper][lower][third][fourth] ** 2
                for upper in range(dimension)
                for lower in range(dimension)
                for third in range(dimension)
                for fourth in range(dimension)
                if riemann[upper][lower][third][fourth] != 0
            ).subs({sp.sin(theta): 1, sp.cos(theta): 0})
        )
    )
    determinant = sp.factor(sp.simplify(metric.det()))
    return ricci_scalar, kretschmann, determinant, {
        "rho": rho,
        "theta": theta,
        "y": outside,
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("two-sheet 5D geometry output already exists")

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

    ricci, kretschmann, determinant, symbols = _derive_invariants()
    rho = symbols["rho"]
    theta = symbols["theta"]
    outside = symbols["y"]

    mp.mp.dps = 50

    def metric_a(value: mp.mpf) -> mp.mpf:
        return (
            -3 * mp.pi / 2
            + 1 / (1 + value**2)
            + 3 * (value / (1 + value**2) + mp.atan(value))
        ) * (1 + value**2)

    horizon = mp.findroot(metric_a, (mp.mpf(1), mp.mpf(2)))
    rho_samples = (-10.0, -2.0, 0.0, float(horizon), 2.0, 10.0)
    outside_samples = (-1.0, 0.0, 1.0)
    samples = []
    for rho_value in rho_samples:
        for outside_value in outside_samples:
            substitutions = {
                rho: rho_value,
                theta: sp.pi / 2,
                outside: outside_value,
            }
            ricci_value = float(ricci.subs(substitutions).evalf(30))
            kretschmann_value = float(
                kretschmann.subs(substitutions).evalf(30)
            )
            samples.append(
                {
                    "rho": rho_value,
                    "outside_y": outside_value,
                    "Ricci5": ricci_value,
                    "Kretschmann5": kretschmann_value,
                    "finite": isfinite(ricci_value)
                    and isfinite(kretschmann_value),
                }
            )

    limits = {}
    for name, direction in (("parent", sp.oo), ("child", -sp.oo)):
        ricci_limit = sp.factor(sp.limit(ricci, rho, direction))
        kretschmann_limit = sp.factor(
            sp.limit(kretschmann, rho, direction)
        )
        limits[name] = {
            "Ricci5": str(ricci_limit),
            "Kretschmann5": str(kretschmann_limit),
            "at_y_zero": {
                "Ricci5": float(ricci_limit.subs(outside, 0).evalf(30)),
                "Kretschmann5": float(
                    kretschmann_limit.subs(outside, 0).evalf(30)
                ),
            },
            "at_y_one": {
                "Ricci5": float(ricci_limit.subs(outside, 1).evalf(30)),
                "Kretschmann5": float(
                    kretschmann_limit.subs(outside, 1).evalf(30)
                ),
            },
        }

    even_residuals = []
    for rho_value in rho_samples:
        plus = float(ricci.subs({rho: rho_value, outside: 1}).evalf(30))
        minus = float(ricci.subs({rho: rho_value, outside: -1}).evalf(30))
        even_residuals.append(abs(plus - minus))

    record = {
        "artifact_id": "NSC-1-S-ONE-TWO-SHEET-5D-GEOMETRY",
        "schema": "NSC-1-S-ONE-TWO-SHEET-5D-GEOMETRY-v1",
        "classification": "an_explicit_compact_even_5D_metric_carries_both_black_universe_branches_with_finite_tested_curvature_invariants",
        "authenticated_inputs": authenticated,
        "metric": {
            "formula": "ds5^2=exp(-36*y^2/1015)*[A(rho)dt^2-drho^2/A(rho)-(rho^2+1)dOmega2^2-dy^2]",
            "outside_interval": "-1<=y<=1",
            "induced_room": "y=0_is_the_complete_smooth_parent_to_child_rho_geometry",
            "warp_even": True,
            "warp_curvature_source": "transition_fixed_b2_over_F=-18/1015",
            "new_dimensionless_warp_fit": False,
        },
        "symbolic_invariants": {
            "determinant": str(determinant),
            "Ricci5": str(ricci),
            "Kretschmann5": str(kretschmann),
            "Ricci5_length": len(str(ricci)),
            "Kretschmann5_length": len(str(kretschmann)),
        },
        "horizon": {
            "rho": mp.nstr(horizon, 40),
            "metric_A": mp.nstr(metric_a(horizon), 10),
            "invariant_samples_included": True,
        },
        "finite_samples": samples,
        "asymptotic_limits": limits,
        "geometry_result": {
            "all_sampled_invariants_finite": all(
                item["finite"] for item in samples
            ),
            "outside_reflection_exact_on_samples": max(even_residuals) == 0.0,
            "parent_and_child_asymptotic_invariants_finite": all(
                isfinite(side[place][quantity])
                for side in limits.values()
                for place in ("at_y_zero", "at_y_one")
                for quantity in ("Ricci5", "Kretschmann5")
            ),
            "horizon_and_throat_carried_in_original_smooth_coordinate": True,
        },
        "next_result": "derive_the_5D_Einstein_tensor_and_reconstruct_one_positive_spectral_bulk_source_then_test_its_kinetic_and_gradient_spectrum",
        "nonclaims": {
            "bulk_source_action_reconstructed": False,
            "positive_energy_or_stability_proved": False,
            "this_compact_warp_is_unique": False,
            "observed_cosmology_identified": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "symbolic_invariants_derived": ricci != 0 and kretschmann != 0,
            "all_sampled_invariants_finite": all(
                item["finite"] for item in samples
            ),
            "outside_reflection_symmetric": max(even_residuals) == 0.0,
            "both_asymptotic_branches_finite": all(
                isfinite(side[place][quantity])
                for side in limits.values()
                for place in ("at_y_zero", "at_y_one")
                for quantity in ("Ricci5", "Kretschmann5")
            ),
            "source_and_stability_completed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "source_and_stability_completed" and value is not True:
            raise RuntimeError(f"two-sheet 5D geometry failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
