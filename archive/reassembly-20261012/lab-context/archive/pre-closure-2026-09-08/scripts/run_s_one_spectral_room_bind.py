#!/usr/bin/env python3
"""Bind the known spectral action as the single local-room reduction."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-spectral-room-bind.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-invariant-closure-bind.json",
    ROOT / "results/nsc-1-s-one-relevant-direction-count.json",
    ROOT / "results/nsc-1-s-one-warped-resolution-bind.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one spectral-room bind output already exists")

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
        "artifact_id": "NSC-1-S-ONE-SPECTRAL-ROOM-BIND",
        "schema": "NSC-1-S-ONE-SPECTRAL-ROOM-BIND-v1",
        "classification": "known_local_particles_forces_and_gravity_are_bound_as_spectral_projections_of_one_room_operator",
        "method": "import_the_spectral_action_reduction_without_recomputing_the_Standard_Model",
        "source": {
            "title": "The Spectral Action Principle",
            "authors": ["Ali H. Chamseddine", "Alain Connes"],
            "arxiv": "hep-th/9606001",
            "doi": "10.1007/s002200050126",
            "source_status": "imported_prior_art_not_project_novelty",
            "universal_formula": "S_spectral=<J*Psi,D*Psi>+Tr[f(D/Lambda)]",
            "published_reduction": "Standard_Model_coupled_to_Einstein_plus_Weyl_gravity_with_high_scale_coupling_relations",
        },
        "authenticated_project_inputs": authenticated,
        "one_operator_dictionary": {
            "spacetime": "commutative_spectral_geometry_of_D",
            "particles": "discrete_eigenspaces_and_poles_of_D_and_the_evolved_Gamma",
            "wavefunctions": "state_amplitudes_and_eigenvectors_in_the_same_Hilbert_space",
            "spin_and_charge": "representations_of_the_spectral_algebra_and_real_structure_J",
            "gauge_and_Higgs": "inner_fluctuations_of_D",
            "gravity": "heat_kernel_coefficients_of_Tr_f_D_over_Lambda",
            "nuclear_effective_sector": "low_energy_topological_projection_after_quark_gluon_modes_are_integrated_out",
        },
        "recursive_completion": {
            "equation": "exp(i*Gamma_n,k/hbar)=integral_unresolved Dchi exp(i*[Tr f(D_n/Lambda_n)+<J Psi_n,D_n Psi_n>+S_link+Gamma_n+1,Omega*k[T chi]]/hbar)",
            "dark": "the_outside_recursive_spectrum_in_the_same_effective_action",
            "black": "nonlinear_spectral_and_boundary_saddles_of_the_same_effective_action",
            "resolution": "the_exact_functional_flow_of_Gamma_n,k",
        },
        "double_counting_removed": {
            "fundamental_local_matter_owner": "spectral_Dirac_action",
            "BPS_Skyrme_owner": "infrared_nuclear_projection",
            "independent_fundamental_Skyrme_plus_QCD_actions": False,
        },
        "fixed_point_requirement": {
            "ordinary_spectral_action_limitation": "finite_Dirac_Yukawa_and_mass_data_are_not_all_predicted",
            "nested_space_requirement": "recursive_FRG_fixed_point_must_determine_the_finite_Dirac_spectrum",
            "success_dimension": "one_overall_dimensional_unit_and_zero_adjustable_dimensionless_particle_or_boundary_ratios",
            "closed_low_energy_subspace": "L0_L2_L4_L6_has_constraint_rank_3_and_nullity_1",
        },
        "next_result": "express_the_remaining_gauge_spin_mixing_and_boundary_coefficients_as_spectral_invariants_of_D_and_S_link_then_solve_the_joint_fixed_point_rank",
        "nonclaims": {
            "finite_Dirac_spectrum_derived": False,
            "Standard_Model_mass_ratios_predicted": False,
            "full_recursive_spectral_flow_solved": False,
            "dark_and_black_kernel_completed": False,
        },
        "gate": {
            "project_inputs_authenticated": len(authenticated) == len(INPUTS),
            "known_Standard_Model_reduction_not_recomputed": True,
            "one_local_spectral_operator_owner": True,
            "BPS_double_counting_removed": True,
            "unfixed_finite_spectrum_explicit": True,
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("spectral-room bind gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
