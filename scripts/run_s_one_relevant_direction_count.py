#!/usr/bin/env python3
"""Count the unfixed directions left by recursive matter self-equality."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-relevant-direction-count.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-reciprocal-closure.json",
    ROOT / "results/nsc-1-s-one-nested-pair-closure.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("relevant-direction count output already exists")

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

    # Columns are kappa0, kappa2, kappa4, kappa6.
    constraints = sp.Matrix(
        [
            [1, 0, 0, -1],  # kappa0-kappa6=0
            [0, 1, -1, 0],  # kappa2-kappa4=0
            [1, -1, 0, 0],  # kappa0-kappa2=0
        ]
    )
    rank = constraints.rank()
    nullspace = constraints.nullspace()
    nullity = constraints.cols - rank
    normalized_nullspace = [
        [str(value) for value in vector] for vector in nullspace
    ]

    record = {
        "artifact_id": "NSC-1-S-ONE-RELEVANT-DIRECTION-COUNT",
        "schema": "NSC-1-S-ONE-RELEVANT-DIRECTION-COUNT-v1",
        "classification": "recursive_self_equality_leaves_exactly_one_overall_unit_direction_in_the_L0246_matter_sector",
        "authenticated_inputs": authenticated,
        "coefficient_order": ["kappa0", "kappa2", "kappa4", "kappa6"],
        "constraint_system": {
            "matrix": [
                [int(value) for value in row] for row in constraints.tolist()
            ],
            "equations": [
                "kappa0-kappa6=0",
                "kappa2-kappa4=0",
                "kappa0-kappa2=0",
            ],
            "rank": rank,
            "nullity": nullity,
            "nullspace": normalized_nullspace,
            "solution": "(kappa0,kappa2,kappa4,kappa6)=kappa*(1,1,1,1)",
        },
        "physical_meaning": {
            "dimensionless_ratios_remaining": 0,
            "overall_dimensional_normalization_remaining": 1,
            "per_operator_retuning_allowed": False,
            "fixed_ratios": [
                "kappa2/kappa0=1",
                "kappa4/kappa0=1",
                "kappa6/kappa0=1",
            ],
        },
        "one_equation_progress": {
            "complete_Theta_dimension": 10,
            "closed_subspace_dimension": 4,
            "closed_subspace_relevant_directions": nullity,
            "remaining_unclosed_coefficients": [
                "nonminimal_chi",
                "degenerate_boundary",
                "gauge_coupling",
                "spin_coupling",
                "mixing_coupling",
                "boundary_coupling",
            ],
        },
        "next_result": "derive_six_independent_cross_sector_constraints_for_the_remaining_boundary_gauge_spin_and_mixing_coefficients_then_recount_the_full_Theta_fixed_point",
        "nonclaims": {
            "full_Theta_has_one_relevant_direction": False,
            "particle_mass_ratios_predicted": False,
            "dark_amplitude_predicted": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "three_constraints_independent": rank == 3,
            "matter_subspace_nullity_one": nullity == 1,
            "nullspace_is_common_normalization": normalized_nullspace
            == [["1", "1", "1", "1"]],
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("relevant-direction count gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
