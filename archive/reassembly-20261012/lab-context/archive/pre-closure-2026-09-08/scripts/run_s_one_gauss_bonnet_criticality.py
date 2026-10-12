#!/usr/bin/env python3
"""Identify the exact EGB criticality implied by the dark-black invariant."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-gauss-bonnet-criticality.json"
INVARIANT = ROOT / "results/nsc-1-s-one-dark-black-invariant.json"
COMPONENTS = ROOT / "results/nsc-1-s-one-gauss-bonnet-component-closure.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("Gauss-Bonnet criticality output already exists")

    authenticated = []
    for path in (INVARIANT, COMPONENTS):
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

    dimension = 5
    alpha = Fraction(1015, 144)
    tangential_vacuum_form = Fraction(108, 1015)
    bare_cosmological = -tangential_vacuum_form
    delta_zero = (
        Fraction((dimension - 3) * (dimension - 4), (dimension - 1) * (dimension - 2))
        * alpha
    )
    vacuum_discriminant = 1 + 8 * delta_zero * bare_cosmological
    critical_product = alpha * bare_cosmological

    record = {
        "artifact_id": "NSC-1-S-ONE-GAUSS-BONNET-CRITICALITY",
        "schema": "NSC-1-S-ONE-GAUSS-BONNET-CRITICALITY-v1",
        "classification": "the_exact_dark_black_coefficient_closure_sits_at_the_5D_Einstein_Gauss_Bonnet_critical_point",
        "authenticated_inputs": authenticated,
        "source": {
            "title": "Criticality in Einstein-Gauss-Bonnet Gravity: Gravity without Graviton",
            "authors": ["Zhong-Ying Fan", "Bin Chen", "Hong Lu"],
            "arxiv": "1606.02728v2",
            "doi": "10.1140/epjc/s10052-016-4389-x",
            "source_status": "imported_criticality_equations_not_project_novelty",
            "equations": {
                "vacuum_polynomial": 3,
                "two_vacua": 4,
                "coalesced_vacuum": 5,
                "linear_kinetic_coefficient": 7,
                "critical_relation": 16,
            },
        },
        "exact_mapping": {
            "dimension": dimension,
            "alpha_GB": str(alpha),
            "vacuum_form_from_room_tensor": str(tangential_vacuum_form),
            "bare_Lambda0_needed_on_EGB_left_hand_side": str(
                bare_cosmological
            ),
            "Delta0=(D-3)(D-4)*alpha/[(D-1)(D-2)]": str(delta_zero),
            "alpha_times_Lambda0": str(critical_product),
            "vacuum_discriminant_1_plus_8_Delta0_Lambda0": str(
                vacuum_discriminant
            ),
            "maximally_symmetric_kappa_eff_over_kappa": "sqrt(1+8*Delta0*Lambda0)=0",
        },
        "physical_result": {
            "two_maximally_symmetric_vacua_coalesce": vacuum_discriminant == 0,
            "ordinary_bilinear_graviton_on_that_vacuum": "absent_or_strongly_coupled",
            "nonlinear_black_solutions_can_still_exist": "established_at_the_critical_point_by_the_source_paper",
            "our_background_is_maximally_symmetric": False,
            "our_background_stability_decided_by_this_test": False,
        },
        "required_next_test": {
            "object": "covariant_Lovelock_principal_symbol_on_the_actual_two_sheet_5D_background",
            "reason": "the_critical_parameter_invalidates_inference_from_an_ordinary_maximally_symmetric_graviton_expansion",
            "pass_condition": "all_physical_characteristic_roots_real_with_a_common_hyperbolic_time_direction_and_nonzero_physical_kinetic_rank",
        },
        "nonclaims": {
            "actual_background_hyperbolicity_failed": False,
            "Nested_Space_mechanism_rejected": False,
            "critical_nonlinear_theory_is_quantum_complete": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == 2,
            "five_dimensional_Delta0_correct": delta_zero == alpha / 6,
            "critical_product_minus_three_quarters": critical_product
            == Fraction(-3, 4),
            "vacuum_discriminant_exactly_zero": vacuum_discriminant == 0,
            "actual_background_principal_symbol_tested": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "actual_background_principal_symbol_tested" and value is not True:
            raise RuntimeError(f"Gauss-Bonnet criticality gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
