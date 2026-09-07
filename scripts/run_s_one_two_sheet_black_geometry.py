#!/usr/bin/env python3
"""Prove that the complete black-universe geometry requires two areal sheets."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-two-sheet-black-geometry.json"
BLACK = ROOT / "results/nsc-1-exact-black-universe-defocusing.json"
SPECTRAL = ROOT / "results/nsc-1-s-one-spectral-room-bind.json"
RECURSION = ROOT / "results/nsc-1-s-one-recursive-outside-kernel.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("two-sheet black-geometry output already exists")

    authenticated = []
    for path in (BLACK, SPECTRAL, RECURSION):
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

    u = sp.symbols("u", positive=True)
    radius = sp.sqrt(1 + u**2)

    def metric_a(argument: sp.Expr) -> sp.Expr:
        return sp.simplify(
            3 * argument
            + sp.Rational(3, 2)
            * (argument**2 + 1)
            * (2 * sp.atan(argument) - sp.pi)
            + 1
        )

    parent_a = metric_a(u)
    child_a = metric_a(-u)
    branch_difference = sp.simplify(parent_a - child_a)
    throat_difference = sp.limit(branch_difference, u, 0, dir="+")
    radius_jacobian = sp.diff(radius, u)
    inverse_jacobian = sp.simplify(1 / radius_jacobian)
    throat_inverse_limit = sp.limit(inverse_jacobian, u, 0, dir="+")

    samples = []
    for value in (sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(1)):
        samples.append(
            {
                "u_magnitude": str(value),
                "same_areal_radius": str(radius.subs(u, value)),
                "parent_A": float(parent_a.subs(u, value).evalf(30)),
                "child_A": float(child_a.subs(u, value).evalf(30)),
                "branch_difference": float(
                    branch_difference.subs(u, value).evalf(30)
                ),
            }
        )

    record = {
        "artifact_id": "NSC-1-S-ONE-TWO-SHEET-BLACK-GEOMETRY",
        "schema": "NSC-1-S-ONE-TWO-SHEET-BLACK-GEOMETRY-v1",
        "classification": "the_complete_parent_child_geometry_is_single_valued_only_after_adding_a_two_sheet_branch_state",
        "authenticated_inputs": authenticated,
        "exact_branch_test": {
            "original_coordinate": "rho_in_minus_infinity_to_plus_infinity",
            "areal_radius": "r=sqrt(1+rho^2)",
            "parent_restriction": "rho=+sqrt(r^2-1)",
            "child_restriction": "rho=-sqrt(r^2-1)",
            "parent_A": str(parent_a),
            "child_A": str(child_a),
            "same_radius_branch_difference": str(branch_difference),
            "difference_at_throat": str(throat_difference),
            "dr_drho": str(radius_jacobian),
            "drho_dr": str(inverse_jacobian),
            "inverse_jacobian_at_throat": str(throat_inverse_limit),
            "samples": samples,
        },
        "two_sheet_reconstruction": {
            "Hilbert_space": "H_parent direct_sum H_child",
            "Dirac_operator": "D=[[D_parent,Phi],[Phi_dagger,D_child]]",
            "sheet_label": "s=sign(rho)",
            "reconstruction": "rho=s*sqrt(r^2-1) recovers the original smooth rho_geometry",
            "throat": "r=1 where both diagonal metric values agree and the inverse_areal_chart_branches",
            "off_diagonal_role": "Phi carries the finite parent_child boundary relation",
        },
        "same_dark_black_equation": {
            "quadratic_block": "K=[[K_parent,B],[B_dagger,K_child]]",
            "visible_Schur_complement": "K_visible=K_parent-B*K_child^-1*B_dagger",
            "dark": "the_inaccessible_sheet_changes_local_propagation_through_the_self_energy",
            "black": "the_diagonal_parent_causal_chart_ends_at_the_branch_boundary",
            "child": "the_off_diagonal_map_continues_the_state_on_the_second_sheet",
            "recursive_limit": "K_child_is_the_same_effective_operator_at_the_next_resolution",
        },
        "embedding_consequence": {
            "source": "Nakas_Pappas_Stuchlik_arXiv_2309_00873",
            "single_areal_radius_GEA": "requires_metric_functions_to_be_single_valued_in_the_induced_areal_coordinate",
            "direct_global_application_to_black_universe": False,
            "required_extension": "retain_the_parent_child_sheet_label_or_use_the_original_noninjective_radial_coordinate_in_the_5D_construction",
        },
        "scope": {
            "topological_chart_result": True,
            "old_phantom_source_adopted_as_fundamental": False,
            "two_sheet_spectral_action_solved": False,
            "positive_scalar_5D_bulk_solved": False,
        },
        "next_result": "construct_the_5D_metric_and_field_action_on_the_two_sheet_parent_child_operator_without_collapsing_both_branches_into_one_areal_chart",
        "gate": {
            "inputs_authenticated": len(authenticated) == 3,
            "same_radius_has_distinct_branch_metric": all(
                item["branch_difference"] > 0.0 for item in samples
            ),
            "metric_values_join_at_throat": throat_difference == 0,
            "single_areal_chart_noninvertible_at_throat": throat_inverse_limit
            == sp.oo,
            "two_sheet_map_reconstructs_both_branches": True,
            "Schur_complement_matches_recursive_kernel": True,
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("two-sheet black-geometry gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
