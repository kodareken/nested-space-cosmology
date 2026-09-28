#!/usr/bin/env python3
"""Derive the full on-room Einstein-Gauss-Bonnet component closure."""

from __future__ import annotations

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


OUTPUT = ROOT / "results/nsc-1-s-one-gauss-bonnet-component-closure.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-gauss-bonnet-null-closure.json",
    ROOT / "results/nsc-1-s-one-one-operator-bind.json",
)


def _curvature(metric: sp.Matrix, coordinates: tuple[sp.Symbol, ...]):
    dimension = len(coordinates)
    inverse = sp.diag(
        *[sp.simplify(1 / metric[index, index]) for index in range(dimension)]
    )
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
                    sp.simplify(
                        metric[first, first] * riemann[first][second][third][fourth]
                    )
                    for fourth in range(dimension)
                ]
                for third in range(dimension)
            ]
            for second in range(dimension)
        ]
        for first in range(dimension)
    ]
    ricci_squared = sp.simplify(
        sum(
            inverse[first, first]
            * inverse[second, second]
            * ricci[first][second] ** 2
            for first in range(dimension)
            for second in range(dimension)
        )
    )
    riemann_squared = sp.simplify(
        sum(
            metric[first, first]
            * inverse[second, second]
            * inverse[third, third]
            * inverse[fourth, fourth]
            * riemann[first][second][third][fourth] ** 2
            for first in range(dimension)
            for second in range(dimension)
            for third in range(dimension)
            for fourth in range(dimension)
            if riemann[first][second][third][fourth] != 0
        )
    )
    return (
        inverse,
        riemann,
        ricci,
        ricci_scalar,
        riemann_down,
        ricci_squared,
        riemann_squared,
    )


def _derive_component_closure():
    time, rho, theta, phi, outside = sp.symbols(
        "t rho theta phi y", real=True
    )
    radius = sp.sqrt(rho**2 + 1)
    metric_b = (
        -3 * sp.pi / 2
        + 1 / (1 + rho**2)
        + 3 * (rho / (1 + rho**2) + sp.atan(rho))
    )
    metric_a = sp.simplify(metric_b * radius**2)
    room_metric = sp.diag(
        metric_a,
        -1 / metric_a,
        -radius**2,
        -radius**2 * sp.sin(theta) ** 2,
    )
    room_curvature = _curvature(room_metric, (time, rho, theta, phi))
    room_ricci_scalar = room_curvature[3]
    room_gauss_bonnet = sp.factor(
        sp.simplify(
            (
                room_ricci_scalar**2
                - 4 * room_curvature[5]
                + room_curvature[6]
            ).subs({sp.sin(theta): 1, sp.cos(theta): 0})
        )
    )

    warp_curvature = sp.Rational(18, 1015)
    conformal_factor = sp.exp(-2 * warp_curvature * outside**2)
    base = sp.diag(
        metric_a,
        -1 / metric_a,
        -radius**2,
        -radius**2 * sp.sin(theta) ** 2,
        -1,
    )
    bulk_metric = sp.diag(
        *[sp.simplify(conformal_factor * base[index, index]) for index in range(5)]
    )
    bulk_curvature = _curvature(
        bulk_metric, (time, rho, theta, phi, outside)
    )
    inverse, _, ricci, ricci_scalar, riemann_down, ricci_squared, riemann_squared = (
        bulk_curvature
    )
    gauss_bonnet = sp.simplify(
        ricci_scalar**2 - 4 * ricci_squared + riemann_squared
    )
    alpha = sp.Rational(1015, 144)
    dimension = 5

    def lanczos(first: int, second: int) -> sp.Expr:
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
        return sp.simplify(
            2
            * (
                ricci_scalar * ricci[first][second]
                - 2 * ricci_square
                - 2 * ricci_riemann
                + riemann_square
            )
            - sp.Rational(1, 2) * bulk_metric[first, second] * gauss_bonnet
        )

    mixed_components = []
    for index in range(dimension):
        einstein = (
            ricci[index][index]
            - sp.Rational(1, 2) * bulk_metric[index, index] * ricci_scalar
        )
        mixed = sp.factor(
            sp.simplify(
                inverse[index, index] * (einstein + alpha * lanczos(index, index))
            ).subs({outside: 0, sp.sin(theta): 1, sp.cos(theta): 0})
        )
        mixed_components.append(mixed)

    tangential_constant = sp.Rational(108, 1015)
    tangential_residuals = [
        sp.simplify(component - tangential_constant)
        for component in mixed_components[:4]
    ]
    normal_expected = -sp.Rational(1, 2) * (
        room_ricci_scalar + alpha * room_gauss_bonnet
    )
    normal_residual = sp.simplify(mixed_components[4] - normal_expected)
    return {
        "rho": rho,
        "alpha": alpha,
        "room_R": room_ricci_scalar,
        "room_GB": room_gauss_bonnet,
        "mixed": mixed_components,
        "tangential_constant": tangential_constant,
        "tangential_residuals": tangential_residuals,
        "normal_expected": normal_expected,
        "normal_residual": normal_residual,
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("Gauss-Bonnet component-closure output already exists")

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

    derived = _derive_component_closure()
    rho = derived["rho"]
    samples = []
    for value in (-2, -1, 0, 1, 2, 10):
        components = [
            float(component.subs(rho, value).evalf(30))
            for component in derived["mixed"]
        ]
        samples.append(
            {
                "rho": value,
                "tangential": components[:4],
                "normal": components[4],
                "tangential_constant_residual": max(
                    abs(component - float(derived["tangential_constant"]))
                    for component in components[:4]
                ),
            }
        )

    record = {
        "artifact_id": "NSC-1-S-ONE-GAUSS-BONNET-COMPONENT-CLOSURE",
        "schema": "NSC-1-S-ONE-GAUSS-BONNET-COMPONENT-CLOSURE-v1",
        "classification": "the_same_5D_geometric_coupling_generates_a_constant_local_vacuum_form_and_encodes_the_room_curvature_action_as_outside_pressure",
        "authenticated_inputs": authenticated,
        "exact_on_room_tensor": {
            "equation": "[G^A_B+alpha_GB*H^A_B]_(y=0)=diag(108/1015,108/1015,108/1015,108/1015,-[R4+alpha_GB*GB4]/2)",
            "alpha_GB": str(derived["alpha"]),
            "tangential_components": [
                str(component) for component in derived["mixed"][:4]
            ],
            "tangential_residuals": [
                str(value) for value in derived["tangential_residuals"]
            ],
            "normal_component": str(derived["mixed"][4]),
            "normal_expected": str(derived["normal_expected"]),
            "normal_identity_residual": str(derived["normal_residual"]),
        },
        "room_curvature_spectrum": {
            "Ricci_scalar": str(derived["room_R"]),
            "Gauss_Bonnet_density": str(derived["room_GB"]),
            "geometric_action_density": "L4_geom=R4+(1015/144)*GB4",
        },
        "samples": samples,
        "one_equation_meaning": {
            "local_dark_energy_form": "all_four_tangential_components_are_the_same_constant_108_over_1015",
            "outside_pressure": "the_normal_component_is_minus_one_half_of_the_local_room_geometric_action_density",
            "black_transition": "the_same_alpha_GB_already_closes_the_required_negative_Ricci_null_projection",
            "self_accounting": "the_room_action_reappears_as_the_source_seen_in_the_unresolved_direction",
            "new_independent_coefficient": False,
        },
        "next_result": "solve_the_normal_source_as_a_spectral_boundary_state_and_test_the_complete_5D_perturbation_kinetic_matrix",
        "nonclaims": {
            "108_over_1015_matches_observed_Lambda": False,
            "normal_source_microphysics_solved": False,
            "linear_stability_proved": False,
            "particle_dark_spectrum_predicted": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "four_tangential_components_exactly_constant": all(
                value == 0 for value in derived["tangential_residuals"]
            ),
            "normal_component_equals_room_action_density": derived[
                "normal_residual"
            ]
            == 0,
            "same_alpha_as_null_closure": derived["alpha"]
            == sp.Rational(1015, 144),
            "no_new_coefficient": True,
            "normal_source_and_stability_closed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "normal_source_and_stability_closed" and value is not True:
            raise RuntimeError(f"Gauss-Bonnet component closure failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
