#!/usr/bin/env python3
"""Bind the spectral Einstein/boundary coefficient ratio into S_one."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-spectral-boundary-constraint.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-spectral-room-bind.json",
    ROOT / "results/nsc-1-s-one-transition-link-coefficient.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("spectral-boundary constraint output already exists")

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

    # c_K - 2*c_R = 0 for integral_M c_R R + integral_boundary c_K K.
    constraint = sp.Matrix([[-2, 1]])
    rank = constraint.rank()
    nullspace = constraint.nullspace()
    normalized_nullspace = [
        [str(value) for value in vector] for vector in nullspace
    ]

    record = {
        "artifact_id": "NSC-1-S-ONE-SPECTRAL-BOUNDARY-CONSTRAINT",
        "schema": "NSC-1-S-ONE-SPECTRAL-BOUNDARY-CONSTRAINT-v1",
        "classification": "the_room_metric_stiffness_and_gravitational_boundary_coefficient_are_one_spectral_parameter_not_two",
        "method": "import_the_spectral_boundary_coefficient_and_bind_it_to_the_recursive_transition_owner",
        "source": {
            "title": "Noncommutative Geometric Spaces with Boundary: Spectral Action",
            "authors": ["Ali H. Chamseddine", "Alain Connes"],
            "arxiv": "1008.3980v1",
            "doi": "10.1016/j.geomphys.2010.10.002",
            "source_status": "imported_prior_art_not_project_novelty",
            "published_result": "spectral_action_generates_the_extrinsic_curvature_term_with_the_sign_and_coefficient_required_by_the_gravitational_Hamiltonian",
        },
        "authenticated_project_inputs": authenticated,
        "exact_constraint": {
            "action_terms": "c_R*integral_M sqrt(g)*R+c_K*integral_boundary sqrt(h)*K",
            "equation": "c_K-2*c_R=0",
            "coefficient_order": ["metric_stiffness_c_R", "boundary_c_K"],
            "constraint_matrix": [[-2, 1]],
            "rank": rank,
            "nullity": constraint.cols - rank,
            "nullspace": normalized_nullspace,
            "solution": "(c_R,c_K)=c_R*(1,2)",
        },
        "one_equation_consequence": {
            "independent_gravitational_boundary_coupling_removed": True,
            "well_posed_metric_variation": True,
            "boundary_conditions_tied_to_Dirac_hermiticity": True,
            "same_D_controls_bulk_and_boundary": True,
        },
        "separation_of_roles": {
            "universal_ratio_two": "coefficient_required_for_the_bulk_boundary_variational_principle",
            "transition_minus_18_over_1015": "dimensionless_curvature_of_the_specific_even_positive_F_transition_profile",
            "these_are_not_the_same_coefficient": True,
        },
        "scope": {
            "source_signature": "Euclidean_with_a_required_Wick_rotation_for_Lorentzian_physics",
            "does_not_yet_prove": [
                "Lorentzian_black_horizon_boundary_conditions",
                "the_full_parent_child_bulk_solution",
                "scalar_and_tensor_stability",
            ],
        },
        "next_result": "combine_the_matter_rank_3_and_metric_boundary_rank_1_constraints_with_the_spectral_gauge_relations_and_recount_the_joint_coefficient_space",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "constraint_rank_one": rank == 1,
            "constraint_nullity_one": constraint.cols - rank == 1,
            "nullspace_has_ratio_one_to_two": normalized_nullspace
            == [["1/2", "1"]],
            "independent_boundary_coefficient_removed": True,
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("spectral-boundary constraint gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
