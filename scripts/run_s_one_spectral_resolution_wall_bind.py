#!/usr/bin/env python3
"""Bind the full spectral high-energy limit to the local truncation boundary."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-spectral-resolution-wall-bind.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-spectral-room-bind.json",
    ROOT / "results/nsc-1-s-one-interface-characteristic-scan.json",
    ROOT / "results/nsc-1-s-one-one-operator-bind.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("spectral-resolution wall bind output already exists")

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
        "artifact_id": "NSC-1-S-ONE-SPECTRAL-RESOLUTION-WALL-BIND",
        "schema": "NSC-1-S-ONE-SPECTRAL-RESOLUTION-WALL-BIND-v1",
        "classification": "the_critical_EGB_complex_principal_roots_mark_a_finite_truncation_boundary_while_the_full_spectral_action_shuts_off_local_high_energy_boson_propagation",
        "authenticated_inputs": authenticated,
        "source": {
            "title": "High energy bosons do not propagate",
            "authors": ["M.A. Kurkov", "Fedele Lizzi", "Dmitri Vassilevich"],
            "arxiv": "1312.2235v1",
            "doi": "10.1016/j.physletb.2014.02.053",
            "source_status": "imported_full_spectral_high_momentum_result_not_project_novelty",
            "setting": "Euclidean_flat_background_nonlocal_heat_kernel_quadratic_in_fields",
            "principal_profile": "f(z)=exp(-z)_with_qualitatively_similar_generic_cutoff_behavior",
        },
        "known_full_spectral_result": {
            "fields": ["scalar", "gauge", "graviton"],
            "high_momentum_quadratic_action": "contains_no_positive_powers_of_derivatives",
            "two_point_behavior": "Green_functions_vanish_for_points_nearer_than_the_inverse_cutoff_scale",
            "interpretation_by_source": "high_energy_bosons_do_not_propagate_and_points_effectively_decouple",
        },
        "binding_to_current_result": {
            "tested_local_truncation": "R5+alpha_GB*GB5+induced_R4",
            "tested_limit": "homogeneous_infinite_frequency_principal_symbol",
            "measured_result": "ten_physical_GR_roots_real_but_diagonal_critical_EGB_interface_retains_complex_roots_for_all_401_normalizations",
            "correct_classification": "finite_heat_kernel_truncation_not_authoritative_at_the_spectral_cutoff",
            "full_spectral_operator_tested_by_that_scan": False,
        },
        "one_equation_resolution_story": {
            "below_cutoff": "local_poles_and_wave_propagation_are_resolved",
            "at_cutoff": "local_bosonic_propagation_shuts_off",
            "across_boundary": "Phi_n_can_remain_finite_and_continue_the_state_on_the_next_sheet",
            "black": "local_causal_propagation_termination",
            "child": "off_diagonal_spectral_continuation_at_a_new_resolution",
        },
        "required_proof": {
            "finite_momentum": "derive_the_complete_curved_Lorentzian_form_factors_of_mathbb_D_Theta",
            "positivity": "require_a_positive_spectral_measure_and_no_unwanted_poles_below_Lambda_n",
            "causality": "require_retarded_support_and_a_common_hyperbolic_low_energy_domain",
            "transition": "require_the_low_momentum_limit_to_reproduce_the_exact_dark_black_background_invariants",
        },
        "nonclaims": {
            "full_curved_Lorentzian_spectral_action_solved": False,
            "complex_EGB_roots_proved_harmless": False,
            "particle_and_dark_spectrum_predicted": False,
            "spectral_cutoff_equals_observed_Planck_length": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "known_high_energy_result_not_recomputed": True,
            "finite_truncation_authority_boundary_explicit": True,
            "full_spectral_finite_momentum_test_completed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "full_spectral_finite_momentum_test_completed" and value is not True:
            raise RuntimeError(f"spectral-resolution wall bind failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
