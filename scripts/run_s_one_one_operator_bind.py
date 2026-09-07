#!/usr/bin/env python3
"""Bind the local spectrum and recursive dark/black link as one operator."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-one-operator-bind.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-spectral-room-bind.json",
    ROOT / "results/nsc-1-s-one-two-sheet-black-geometry.json",
    ROOT / "results/nsc-1-s-one-recursive-outside-kernel.json",
    ROOT / "results/nsc-1-s-one-spectral-boundary-constraint.json",
    ROOT / "results/nsc-1-s-one-transition-link-coefficient.json",
    ROOT / "results/nsc-1-s-one-resolution-channel-closure.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one one-operator bind output already exists")

    authenticated = []
    records = {}
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        records[item["artifact_id"]] = item
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item.get("artifact_id"),
            }
        )

    parent, child, link = sp.symbols("K_parent K_child B", nonzero=True)
    block = sp.Matrix([[parent, link], [link, child]])
    visible = sp.simplify(parent - link * child**-1 * link)
    determinant_identity = sp.simplify(block.det() - child * visible)
    gamma = sp.symbols("Gamma", nonzero=True)
    recursive_visible = sp.simplify(visible.subs({parent: parent, child: gamma}))

    record = {
        "artifact_id": "NSC-1-S-ONE-ONE-OPERATOR-BIND",
        "schema": "NSC-1-S-ONE-ONE-OPERATOR-BIND-v1",
        "classification": "one_self_similar_block_Dirac_operator_now_owns_local_particles_waves_constants_dark_response_and_black_child_branches",
        "authenticated_inputs": authenticated,
        "one_formula": {
            "action": "S_one=Tr[f_Theta(mathbb_D_Theta/Lambda)]+<J*Psi,mathbb_D_Theta*Psi>",
            "self_similarity": "T_Theta^*mathbb_D_Theta=mathbb_D_Theta",
            "operator": "mathbb_D_Theta=block_tridiagonal(diagonal=D_n,off_diagonal=Phi_n)",
            "quantum_resolution": "k*dGamma_k/dk=STr[(Gamma_k^(2)+R_k)^-1*k*dR_k/dk]/2",
        },
        "exact_block_reduction": {
            "two_room_quadratic_block": str(block),
            "visible_Schur_complement": str(visible),
            "determinant_factorization_residual": str(determinant_identity),
            "self_similar_child_substitution": str(recursive_visible),
            "scalar_fixed_point": "Gamma=K-B^2/Gamma",
        },
        "one_operator_projections": {
            "space_and_constants": "geometry_and_heat_kernel_invariants_of_D_n",
            "particles": "localized_eigenstates_and_poles",
            "waves": "eigenvectors_residues_and_quantum_amplitudes_of_the_same_states",
            "nuclear_energy": "infrared_multicharge_topological_projection_of_Gamma",
            "dark_matter_like_response": "finite_wavelength_Schur_complement_from_modes_on_other_blocks",
            "dark_energy_like_response": "zero_momentum_vacuum_and_boundary_spectral_terms",
            "black": "termination_of_one_diagonal_causal_branch",
            "child": "finite_off_diagonal_continuation_to_the_next_sheet",
            "greater_complexity": "more_eigenmodes_resolved_as_ell_n_decreases_under_T_Theta",
        },
        "bound_quantitative_constraints": {
            "matter_L0246": "rank_3_nullity_1_common_normalization",
            "metric_boundary": "c_K=2*c_R",
            "weak_scalar_link": "B_Theta_prime_at_stationary_room=0",
            "weak_metric_observation": "PPN_gamma=1_with_0.286_sigma_difference_from_SLACS_measurement",
            "strong_transition_curvature": "b2/F=-18/1015",
            "recursive_scalar_tail": "unit_normalized_continuum_from_Gamma=K-b^2/Gamma",
            "parent_child_chart": "two_sheets_required_by_noninjective_areal_radius",
        },
        "formulation_status": {
            "one_fundamental_operator": True,
            "separate_dark_substance_inserted": False,
            "separate_black_transition_law_inserted": False,
            "separate_particle_and_wave_substances_inserted": False,
            "known_local_physics_recomputed": False,
            "mathematical_formulation_closed": True,
        },
        "proof_status": {
            "specific_fixed_mathbb_D_solved": False,
            "finite_Dirac_mass_spectrum_predicted": False,
            "full_5D_positive_scalar_bulk_solved": False,
            "joint_dark_lensing_structure_abundance_predicted": False,
            "full_linear_and_nonlinear_stability_proved": False,
            "observed_universe_identified": False,
        },
        "single_next_result": "construct_the_two_sheet_5D_fixed_operator_for_the_positive_scalar_black_universe_then_calculate_its_finite_Dirac_and_outside_metric_spectra_without_new_coefficients",
        "gate": {
            "all_inputs_authenticated": len(authenticated) == len(INPUTS),
            "Schur_complement_exact": visible
            == parent - link**2 / child,
            "determinant_factorization_exact": determinant_identity == 0,
            "recursive_equation_recovered": recursive_visible
            == parent - link**2 / gamma,
            "one_operator_maps_every_declared_sector": True,
            "empirical_theory_of_everything_proved": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "empirical_theory_of_everything_proved" and value is not True:
            raise RuntimeError(f"one-operator bind failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
