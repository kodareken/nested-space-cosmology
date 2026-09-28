#!/usr/bin/env python3
"""Apply the same self-dual balancing law to gradients of gradients."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.unified_action import NestedDualPairEnergy  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-nested-pair-closure.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one nested-pair closure output already exists")
    a, s = sp.symbols("a s", positive=True)
    k0, k2 = sp.symbols("kappa_0 kappa_2", positive=True)
    pair_06 = k0 * (a**3 + a**-3)
    pair_24 = k2 * (a + a**-1)
    energy = sp.Rational(1, 2) * (s**2 * pair_06 + s**-2 * pair_24)
    state_equation = sp.factor(sp.diff(energy, s))
    positive_state = sp.simplify((pair_24 / pair_06) ** sp.Rational(1, 4))
    stationary_residual = sp.simplify(state_equation.subs(s, positive_state))
    state_hessian = sp.simplify(sp.diff(energy, s, 2).subs(s, positive_state))
    fixed_equation = sp.simplify(state_equation.subs({a: 1, s: 1}))
    coefficient_solution = sp.solve(sp.Eq(fixed_equation, 0), k2)
    if coefficient_solution != [k0]:
        raise RuntimeError("nested-pair fixed point did not determine kappa2=kappa0")
    fully_closed = sp.simplify(energy.subs(k2, k0))
    room_gradient = sp.simplify(sp.diff(fully_closed, a).subs({a: 1, s: 1}))
    room_hessian = sp.simplify(sp.diff(fully_closed, a, 2).subs({a: 1, s: 1}))
    fixture = NestedDualPairEnergy(pair_06=18.0, pair_24=2.0)
    result = {
        "artifact_id": "NSC-1-S-ONE-NESTED-PAIR-CLOSURE",
        "schema": "NSC-1-S-ONE-NESTED-PAIR-CLOSURE-v1",
        "classification": "recursive_self_equality_fixes_all_L0_L2_L4_L6_coefficients_to_one_normalization",
        "nested_energy": {
            "pair_06": str(pair_06),
            "pair_24": str(pair_24),
            "formula": str(energy),
            "gradient_state_equation": str(state_equation),
            "positive_gradient_state": str(positive_state),
            "stationary_residual": str(stationary_residual),
            "positive_state_hessian": str(state_hessian),
        },
        "room_fixed_point": {
            "a": "1",
            "s": "1",
            "field_equation": str(fixed_equation),
            "forced_relation": "kappa2=kappa0",
            "prior_relations": ["kappa6=kappa0", "kappa4=kappa2"],
            "complete_relation": "kappa0=kappa2=kappa4=kappa6",
            "fully_closed_energy": str(fully_closed),
            "scale_gradient": str(room_gradient),
            "scale_hessian": str(room_hessian),
        },
        "recursive_meaning": {
            "first_level": "each_gradient_term_balances_its_reciprocal_dual",
            "second_level": "the_complete_long_short_pair_balances_the_complete_bulk_topological_pair",
            "same_law_each_level": "x^2*A+x^-2*B_is_stationary_when_the_two_weighted_contributions_are_equal",
            "new_coefficient_added": False,
        },
        "public_fixture": {
            "pair_06": fixture.pair_06,
            "pair_24": fixture.pair_24,
            "derived_scale_state": fixture.scale_state,
            "stationary_energy": fixture.energy(),
            "each_pair_contribution": fixture.pair_contribution,
            "self_consistent": fixture.self_consistent,
        },
        "Theta_fixed_point_progress": {
            "Skyrme_scale_coefficients_before_self_equality": 4,
            "Skyrme_scale_coefficients_after_first_duality": 2,
            "Skyrme_scale_coefficients_after_nested_duality": 1,
            "remaining_value": "one_overall_energy_length_normalization_only",
        },
        "gate": {
            "nested_state_variation_closes": stationary_residual == 0,
            "nested_state_stable": state_hessian.is_positive is True,
            "all_four_coefficients_equal": coefficient_solution == [k0],
            "room_stationary": room_gradient == 0,
            "room_stable": room_hessian == 10 * k0,
            "numeric_fixture_self_consistent": fixture.self_consistent,
        },
        "next_result": "use_the_single_remaining_energy_length_anchor_to_solve_the_shared_spinor_gauge_and_topological_charge_sectors",
        "nonclaims": {
            "overall_dimensional_scale_predicted": False,
            "electron_solution_obtained": False,
            "particle_mass_ratios_predicted": False,
            "gravity_boundary_coefficients_closed": False,
        },
        "terminal": True,
    }
    if not all(result["gate"].values()):
        raise RuntimeError("S-one nested-pair closure gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(result))
    return result


def main() -> int:
    result = run()
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
