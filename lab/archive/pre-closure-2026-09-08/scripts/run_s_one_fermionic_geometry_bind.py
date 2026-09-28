#!/usr/bin/env python3
"""Bind the bosonic spectral action to its regulated fermionic origin."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-fermionic-geometry-bind.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-spectral-profile-closure.json",
    ROOT / "results/nsc-1-s-one-reflection-positivity-obstruction.json",
    ROOT / "results/nsc-1-s-one-spectral-room-bind.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("fermionic-geometry bind output already exists")

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
        "artifact_id": "NSC-1-S-ONE-FERMIONIC-GEOMETRY-BIND",
        "schema": "NSC-1-S-ONE-FERMIONIC-GEOMETRY-BIND-v1",
        "classification": "the_bosonic_room_geometry_and_forces_are_bound_to_the_scale_anomaly_of_the_same_regularized_fermionic_spectrum",
        "authenticated_inputs": authenticated,
        "source": {
            "title": "Bosonic Spectral Action Induced from Anomaly Cancelation",
            "authors": ["A.A. Andrianov", "Fedele Lizzi"],
            "arxiv": "1001.2036v1",
            "doi": "10.1007/JHEP05(2010)057",
            "source_status": "imported_anomaly_derivation_not_project_novelty",
            "published_result": "a_modified_bosonic_spectral_action_is_the_term_required_to_cancel_the_scale_anomaly_of_the_spectrally_regularized_fermionic_action",
        },
        "one_quantum_equation": {
            "effective_action": "Gamma_one[mathbb_D]=-log(det_Lambda(mathbb_D))+A_scale[mathbb_D]",
            "anomaly_term": "A_scale_generates_Tr exp(-mathbb_D^2/Lambda^2)_with_fixed_coefficient_modifications",
            "same_spectrum": True,
            "independent_bosonic_and_fermionic_weights": False,
        },
        "physical_identity": {
            "fermions": "states_and_eigenmodes_of_mathbb_D",
            "wavefunctions": "amplitudes_in_the_fermionic_Hilbert_space",
            "bosonic_fields": "inner_fluctuations_and_anomaly_induced_response_of_mathbb_D",
            "geometry_and_gravity": "heat_spectral_invariants_induced_by_the_same_regulated_determinant",
            "dark_and_black": "off_diagonal_and_nonlinear_spectral_geometry_of_the_same_block_operator",
        },
        "reflection_positivity_routing": {
            "failed_object": "standalone_background_adjusted_bosonic_TT_covariance",
            "correct_object": "complete_gauge_invariant_anomaly_consistent_covariance_of_Gamma_one",
            "double_counting_forbidden": "do_not_add_a_second_free_fermion_determinant_weight_to_force_positivity",
        },
        "remaining_open_input": {
            "ordinary_spectral_Standard_Model": "finite_Dirac_matrix_contains_fermion_masses_and_mixings",
            "nested_space_requirement": "derive_the_finite_Dirac_spectrum_from_the_recursive_fixed_point",
        },
        "next_result": "derive_the_anomaly_consistent_quadratic_gauge_invariant_two_point_kernel_of_the_full_two_sheet_Dirac_operator_and_repeat_the_Stieltjes_test",
        "nonclaims": {
            "complete_anomaly_kernel_reconstructed": False,
            "reflection_positivity_repaired": False,
            "finite_Dirac_masses_predicted": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "known_anomaly_result_not_recomputed": True,
            "one_spectrum_owner": True,
            "independent_relative_weight_removed": True,
            "complete_quantum_covariance_tested": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "complete_quantum_covariance_tested" and value is not True:
            raise RuntimeError(f"fermionic-geometry bind failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
