#!/usr/bin/env python3
"""Correct the missing owner after anomaly-compensated scale flow."""

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


OUTPUT = ROOT / "results/nsc-2-zeta1-anomaly-owner-correction.json"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-anomaly-decomposition.json",
    ROOT / "results/nsc-1-s-one-invariant-partition-bind.json",
    ROOT / "results/nsc-1-s-one-fermionic-geometry-bind.json",
    ROOT / "results/nsc-1-s-one-recursive-outside-kernel.json",
)


def _authenticate() -> list[dict[str, object]]:
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
                "artifact_id": item["artifact_id"],
            }
        )
    return authenticated


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 anomaly-owner correction output already exists")
    authenticated = _authenticate()

    partition, transformed, anomaly = sp.symbols(
        "Z_D Z_Dphi S_anom", positive=True
    )
    source_identity = sp.Eq(anomaly, sp.log(partition / transformed))
    transformed_solution = sp.solve(source_identity, transformed)[0]
    expected_transformed = partition * sp.exp(-anomaly)
    residual = sp.simplify(transformed_solution - expected_transformed)

    local, link, tail = sp.symbols("K B Gamma_tail", nonzero=True)
    joined = sp.simplify(local - link**2 / tail)
    self_equal = sp.symbols("Gamma", nonzero=True)
    fixed_equation = sp.Eq(self_equal, joined.subs(tail, self_equal))
    fixed_solutions = sp.solve(fixed_equation, self_equal)

    record = {
        "artifact_id": "NSC-2-ZETA1-ANOMALY-OWNER-CORRECTION",
        "schema": "NSC-2-ZETA1-ANOMALY-OWNER-CORRECTION-v1",
        "classification": "the_monotone_single_transition_scale_flow_is_missing_the_recursive_child_tail_not_an_independently_weighted_geometric_action",
        "authenticated_inputs": authenticated,
        "primary_anomaly_authority": {
            "title": "Spectral_action_Weyl_anomaly_and_the_Higgs_Dilaton_potential",
            "arxiv": "1106.3263v1",
            "doi": "10.1007/JHEP10(2011)001",
            "equations": ["caseA", "splitting", "zinvprod", "sanom1", "Sanomal", "sanomcoll"],
            "source_identity": str(source_identity),
            "transformed_partition": str(transformed_solution),
            "identity_residual": str(residual),
            "sign_note": "Z_one=integral_Dphi_Z(D_phi)=Z(D)*integral_Dphi_exp(-S_anom)_for_the_caseA_choice_already_used_by_the_project",
        },
        "double_counting_boundary": {
            "fermionic_owner": "minus_log_det_Lambda_D",
            "induced_geometry_owner": "the_compensating_scale_anomaly_of_that_same_determinant",
            "forbidden_repair": "add_the_existing_Einstein_Gauss_Bonnet_or_heat_action_again_with_an_independent_weight",
            "reason": "the_current_one_equation_defines_geometry_as_the_anomaly_induced_response_of_the_same_spectrum",
        },
        "missing_recursive_term": {
            "one_transition_used_so_far": "joined_parent_to_one_child_minus_disconnected_parent_plus_terminal_child",
            "theory_requires": "the_child_Dirichlet_to_Neumann_map_contains_its_own_next_child_tail_under_T_Theta",
            "quadratic_reduction": str(fixed_equation),
            "solutions": [str(value) for value in fixed_solutions],
            "physical_branch": "Gamma=(K+sqrt(K^2-4*B^2))/2_with_Gamma_to_K_as_B_to_zero",
            "full_boundary_equation": "R_tail(E,zeta)=R_local(E,zeta)-Phi_dagger*R_tail(E,zeta)^-1*Phi",
        },
        "causal_correction": {
            "preserved_result": "the_anomaly_compensated_single_transition_child_link_derivative_is_negative_on_the_scanned_domain",
            "wrong_next_owner": "an_independent_classical_or_anomaly_induced_geometric_boundary_weight",
            "correct_next_owner": "the_self_consistent_recursive_Calderon_or_Dirichlet_to_Neumann_tail",
            "next_result": "solve_the_mode_resolved_boundary_Riccati_fixed_point_then_recompute_the_anomaly_compensated_scale_derivative",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "primary_partition_identity_exact": residual == 0,
            "independent_geometry_double_count_forbidden": True,
            "recursive_tail_absent_from_single_transition_result": True,
            "scalar_recursive_equation_recovered": fixed_equation
            == sp.Eq(self_equal, local - link**2 / self_equal),
            "mode_resolved_recursive_tail_solved": False,
            "zeta_derived": False,
        },
        "nonclaims": {
            "scalar_Riccati_solution_is_the_full_transition_spectrum": False,
            "recursive_tail_is_known_to_restore_a_scale_root": False,
            "measure_and_Weyl_invariant_functional_fully_derived": False,
            "physical_zeta_promoted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key in {"mode_resolved_recursive_tail_solved", "zeta_derived"}:
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 anomaly owner correction failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
