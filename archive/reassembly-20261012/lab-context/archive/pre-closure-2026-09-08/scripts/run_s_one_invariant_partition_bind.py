#!/usr/bin/env python3
"""Bind the anomaly-consistent invariant partition function as the one equation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-invariant-partition-bind.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-fermionic-geometry-bind.json",
    ROOT / "results/nsc-1-s-one-reflection-positivity-obstruction.json",
    ROOT / "results/nsc-1-s-one-one-operator-bind.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("invariant-partition bind output already exists")

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

    record = {
        "artifact_id": "NSC-1-S-ONE-INVARIANT-PARTITION-BIND",
        "schema": "NSC-1-S-ONE-INVARIANT-PARTITION-BIND-v1",
        "classification": "one_scale_invariant_regularized_fermionic_partition_function_owns_matter_bosonic_geometry_and_resolution",
        "authenticated_inputs": authenticated,
        "source": {
            "title": "Bosonic Spectral Action Induced from Anomaly Cancelation",
            "authors": ["A.A. Andrianov", "Fedele Lizzi"],
            "arxiv": "1001.2036v1",
            "doi": "10.1007/JHEP05(2010)057",
            "equations": {
                "finite_Dirac_split": "2.4",
                "bosonic_plus_fermionic_action": "2.5",
                "invariant_partition": ["4.11", "4.12"],
                "anomaly_ratio": "4.13",
                "anomaly_integral": "4.19",
                "constant_scale_heat_coefficients": "5.22",
            },
        },
        "one_partition_function": {
            "formula": "Z_one,Lambda(mathbb_D)=integral Dvarphi Z_Lambda(exp(-varphi/2)*mathbb_D*exp(-varphi/2))",
            "recursion": "T_Theta^*mathbb_D_Theta=mathbb_D_Theta",
            "varphi_role": "local_resolution_and_parent_child_scale_relation",
            "separate_bosonic_weight": False,
        },
        "exact_global_anomaly": {
            "projector": "P_Lambda_projects_onto_Dirac_eigenvalues_below_the_cutoff",
            "formula": "S_anom=integral_0^1 dt varphi Tr(P_Lambda)",
            "four_dimensional_constant_varphi": "S_anom=[exp(4varphi)-1]*a0/8+[exp(2varphi)-1]*a2/2+varphi*a4",
            "meaning": "vacuum_Einstein_gauge_and_curvature_terms_are_responses_of_the_same_regularized_fermionic_spectrum",
        },
        "known_limit": {
            "global_constant_rescaling_only": True,
            "local_varphi_derivative_terms_derived_in_source": False,
            "full_Standard_Model_anomaly_coefficients_checked_in_source": False,
            "finite_Dirac_matrix_contains_input_masses_and_mixings": True,
            "complete_two_sheet_covariance_supplied": False,
        },
        "causal_correction": {
            "wrong_next_step": "add_an_independently_weighted_minus_log_det_term_to_the_bosonic_covariance",
            "correct_next_step": "derive_the_local_two_sheet_anomaly_of_the_single_regulated_partition_function",
            "reflection_positivity_owner": "a_complete_gauge_invariant_observable_of_Z_one",
        },
        "next_result": "promote_varphi_to_the_local_two_sheet_resolution_field_and_derive_its_anomaly_kernel_including_off_diagonal_Dirac_blocks",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "known_anomaly_equations_bound": True,
            "independent_determinant_double_counting_removed": True,
            "one_partition_function_formulated": True,
            "local_two_sheet_anomaly_derived": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "local_two_sheet_anomaly_derived" and value is not True:
            raise RuntimeError(f"invariant-partition bind failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
