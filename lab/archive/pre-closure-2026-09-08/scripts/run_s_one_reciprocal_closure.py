#!/usr/bin/env python3
"""Solve the L0/L2/L4/L6 coefficient fixed point from self-equality."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.unified_action import ReciprocalGradientEnergy  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-reciprocal-closure.json"


def _text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one reciprocal closure output already exists")
    a = sp.symbols("a", positive=True)
    k0, k2, k4, k6 = sp.symbols("kappa_0 kappa_2 kappa_4 kappa_6")
    energy = k0 * a**3 + k2 * a + k4 / a + k6 / a**3
    reciprocal = sp.simplify(energy.subs(a, 1 / a))
    difference = sp.collect(sp.expand(a**3 * (energy - reciprocal)), a)
    polynomial = sp.Poly(difference, a)
    coefficient_equations = [sp.Eq(value, 0) for value in polynomial.all_coeffs()]
    solved = sp.solve(coefficient_equations, (k4, k6), dict=True)
    if solved != [{k4: k2, k6: k0}]:
        raise RuntimeError("reciprocal coefficient closure differs")
    closed_energy = sp.simplify(energy.subs(solved[0]))
    fixed_gradient = sp.simplify(sp.diff(closed_energy, a).subs(a, 1))
    fixed_hessian = sp.simplify(sp.diff(closed_energy, a, 2).subs(a, 1))
    force = -sp.diff(closed_energy, a)
    force_identity = sp.simplify(force.subs(a, 1 / a) + a**2 * force)

    fixture = ReciprocalGradientEnergy(
        kappa0=Fraction(3, 2),
        kappa2=Fraction(5, 4),
        kappa4=Fraction(5, 4),
        kappa6=Fraction(3, 2),
    )
    pairs = []
    for scale in (Fraction(1, 4), Fraction(2, 3), Fraction(3, 2), Fraction(4)):
        inverse = 1 / scale
        pairs.append(
            {
                "scale": _text(scale),
                "reciprocal_scale": _text(inverse),
                "energy": _text(fixture.energy(scale)),
                "reciprocal_energy": _text(fixture.energy(inverse)),
                "force": _text(fixture.force(scale)),
                "reciprocal_force": _text(fixture.force(inverse)),
                "energy_equal": fixture.energy(scale) == fixture.energy(inverse),
                "force_reversal_exact": fixture.force(inverse)
                == -scale**2 * fixture.force(scale),
            }
        )
    result = {
        "artifact_id": "NSC-1-S-ONE-RECIPROCAL-CLOSURE",
        "schema": "NSC-1-S-ONE-RECIPROCAL-CLOSURE-v1",
        "classification": "self_equality_removes_two_independent_coefficients_and_creates_one_stable_room_fixed_point",
        "Derrick_scaled_energy": {
            "formula": str(energy),
            "terms": {
                "L0": "kappa0*a^3_bulk_volume",
                "L2": "kappa2*a_long_gradient",
                "L4": "kappa4/a_short_gradient",
                "L6": "kappa6/a^3_topological_compression",
            },
        },
        "self_equality": {
            "requirement": "E(a)=E(1/a)",
            "polynomial_identity": str(difference),
            "unique_coefficient_relations": {
                "kappa4": "kappa2",
                "kappa6": "kappa0",
            },
            "closed_energy": str(closed_energy),
            "independent_coefficients_before": 4,
            "independent_coefficients_after": 2,
        },
        "fixed_room": {
            "scale": "1",
            "first_derivative": str(fixed_gradient),
            "second_derivative": str(fixed_hessian),
            "stable_for_positive_coefficients": fixed_hessian
            == 18 * k0 + 2 * k2,
        },
        "parent_child_force": {
            "identity": "F(1/a)=-a^2*F(a)",
            "symbolic_residual": str(force_identity),
            "fixture_pairs": pairs,
        },
        "Theta_fixed_point_progress": {
            "determined_relations": ["kappa0=kappa6", "kappa2=kappa4"],
            "remaining_dimensionless_unknowns": [
                "kappa2_over_kappa0",
                "nonminimal_and_degenerate_boundary_relation",
                "gauge_spin_mixing_relation",
            ],
            "mass_parameters_inserted": False,
        },
        "gate": {
            "symbolic_energy_self_duality": sp.simplify(
                closed_energy - closed_energy.subs(a, 1 / a)
            )
            == 0,
            "fixed_point_stationary": fixed_gradient == 0,
            "fixed_point_stable": fixed_hessian == 18 * k0 + 2 * k2,
            "force_reversal": force_identity == 0,
            "all_exact_fixtures": all(
                item["energy_equal"] and item["force_reversal_exact"]
                for item in pairs
            ),
        },
        "next_result": "derive_kappa2_over_kappa0_from_the_shared_particle_spin_charge_boundary_stationarity_instead_of_fitting_a_particle_mass",
        "nonclaims": {
            "all_Theta_components_determined": False,
            "numerical_constants_predicted": False,
            "electron_stationary_solution_obtained": False,
            "cosmological_holdout_opened": False,
        },
        "terminal": True,
    }
    if not all(result["gate"].values()):
        raise RuntimeError("S-one reciprocal closure gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(result))
    return result


def main() -> int:
    result = run()
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
