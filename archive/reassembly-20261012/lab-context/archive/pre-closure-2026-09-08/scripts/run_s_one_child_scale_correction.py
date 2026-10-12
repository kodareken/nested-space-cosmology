#!/usr/bin/env python3
"""Correct the local-vacuum versus child-asymptotic scale identification.

The on-room 5D EGB tensor contains the constant 108/1015, but the exact
black-universe child asymptotes to V=9*pi and R=-36*pi.  In the doubled-FLRW
convention R=-2*Lambda_effective, so the child coefficient is 18*pi/L_star^2.
It is not the local on-room tensor coefficient.

Matching the doubled spectral vacuum to the actual child gives
mu^2=1-3*pi/(2*zeta), zeta=(Lambda*L_star)^2.  A real positive gap requires
zeta>3*pi/2.  The previously explored zeta=1 gap is therefore not the
child-matched physical branch.  Its exact pole/correlator algebra remains a
controlled conditional calculation, while its cosmological interpretation is
withdrawn.
"""

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


OUTPUT = ROOT / "results/nsc-1-s-one-child-scale-correction.json"
INPUTS = (
    ROOT / "results/nsc-1-exact-black-universe-defocusing.json",
    ROOT / "results/nsc-1-s-one-doubled-flrw-closure.json",
    ROOT / "results/nsc-1-s-one-dark-black-invariant.json",
    ROOT / "results/nsc-1-s-one-gap-transition-closure.json",
    ROOT / "results/nsc-1-s-one-child-orientation.json",
    ROOT / "results/nsc-1-s-one-boundary-retarded-pole.json",
    ROOT / "results/nsc-1-s-one-fermionic-relative-observable.json",
)


def _authenticate() -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    authenticated: list[dict[str, object]] = []
    records: dict[str, dict[str, object]] = {}
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        artifact_id = str(item["artifact_id"])
        records[artifact_id] = item
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": artifact_id,
            }
        )
    return authenticated, records


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("child-scale correction output already exists")
    authenticated, records = _authenticate()
    exact_child = records["NSC-1-EXACT-BLACK-UNIVERSE-DEFOCUSING"]
    dark_black = records["NSC-1-S-ONE-DARK-BLACK-INVARIANT"]

    interior = exact_child["global_geometry"]["exact_interior_limits"]
    potential_limit = sp.sympify(interior["V_minus_infinity"])
    child_ricci = -4 * potential_limit
    child_lambda_effective = -child_ricci / 2
    local_room_vacuum = sp.Rational(
        dark_black["one_number"]["local_vacuum_form_coefficient"]
    )

    zeta, mu_squared = sp.symbols(
        "zeta mu_squared", positive=True, real=True
    )
    matching_equation = sp.Eq(
        child_lambda_effective, 12 * zeta * (1 - mu_squared)
    )
    matching_solution = sp.solve(matching_equation, mu_squared)[0]
    one_scale_value = sp.simplify(matching_solution.subs(zeta, 1))
    reality_threshold = sp.Rational(3, 2) * sp.pi
    local_child_difference = sp.simplify(
        child_lambda_effective - local_room_vacuum
    )

    record = {
        "artifact_id": "NSC-1-S-ONE-CHILD-SCALE-CORRECTION",
        "schema": "NSC-1-S-ONE-CHILD-SCALE-CORRECTION-v1",
        "classification": "the_actual_child_asymptotic_rejects_the_direct_one_scale_gap_and_requires_a_derived_recursive_crossover_above_three_pi_over_two",
        "authenticated_inputs": authenticated,
        "exact_child_asymptotic": {
            "V_minus_infinity": str(potential_limit),
            "R_minus_infinity": str(child_ricci),
            "doubled_FLRW_identity": "R_child=-2*Lambda_effective",
            "Lambda_effective_times_L_star_squared": str(
                child_lambda_effective
            ),
        },
        "corrected_cross_scale_equation": {
            "zeta": "(Lambda*L_star)^2",
            "mu_squared": "abs(Phi)^2/Lambda^2",
            "equation": str(matching_equation),
            "solution": str(matching_solution),
            "expanded_solution": "mu_squared=1-3*pi/(2*zeta)",
            "positive_gap_requirement": "zeta>3*pi/2",
            "threshold_exact": str(reality_threshold),
            "threshold_decimal": float(reality_threshold),
        },
        "measured_correction": {
            "local_on_room_vacuum_form": str(local_room_vacuum),
            "child_asymptotic_vacuum": str(child_lambda_effective),
            "difference": str(local_child_difference),
            "same_magnitude": local_child_difference == 0,
            "zeta_one_gap_squared": str(one_scale_value),
            "zeta_one_gap_real_positive": bool(one_scale_value > 0),
            "responsible_assumption": "Lambda_times_L_star_was_set_to_one_before_the_recursive_inheritance_map_was_solved",
        },
        "preserved_and_superseded": {
            "preserved_exact_results": [
                "same_Phi_mass_gap_and_visible_Schur_self_energy",
                "same_Phi_relative_metric_coupling_and_vacuum_shift",
                "full_exponential_flat_relative_heat_kernel_scan",
                "conditional_two_sheet_retarded_pole_and_Lehmann_weight_algebra",
            ],
            "superseded_physical_interpretations": [
                "1006_over_1015_is_the_child_matched_gap",
                "54_over_503_is_a_child_cosmological_invariant",
                "the_conditional_boundary_pole_has_a_derived_physical_mass_before_zeta_is_solved",
            ],
            "historical_results_rewritten": False,
        },
        "causal_result": {
            "first_owner": "recursive_parent_child_scale_inheritance",
            "required_next_result": "derive_zeta_from_the_same_transition_operator_without_using_the_observed_child_vacuum_then_insert_it_into_mu_squared_equals_one_minus_three_pi_over_two_zeta",
            "parameter_fit_authorized": False,
            "theory_rejection": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "child_potential_limit_exact": potential_limit == 9 * sp.pi,
            "child_Ricci_limit_exact": child_ricci == -36 * sp.pi,
            "child_effective_vacuum_exact": child_lambda_effective
            == 18 * sp.pi,
            "local_room_and_child_vacua_distinct": bool(
                local_child_difference > 0
            ),
            "general_scale_solution_exact": sp.simplify(
                matching_solution
                - (1 - 3 * sp.pi / (2 * zeta))
            )
            == 0,
            "zeta_one_positive_gap_rejected": bool(one_scale_value < 0),
            "recursive_zeta_derived": False,
            "physical_gap_fixed": False,
            "conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "observed_Lambda_used": False,
            "zeta_value_predicted": False,
            "electron_mass_predicted": False,
            "child_matched_boundary_mass_predicted": False,
            "observed_universe_identified": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key not in {
            "recursive_zeta_derived",
            "physical_gap_fixed",
            "conviction_proof_seed_completed",
        } and value is not True:
            raise RuntimeError(f"child-scale correction failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
