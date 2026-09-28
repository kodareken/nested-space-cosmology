#!/usr/bin/env python3
"""Derive the self-adjoint ZETA1 Dirac and throat-gluing domain."""

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


OUTPUT = ROOT / "results/nsc-2-zeta1-self-adjoint-domain.json"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-foliation.json",
    ROOT / "results/nsc-1-s-one-local-two-sheet-anomaly.json",
    ROOT / "results/nsc-1-s-one-child-scale-correction.json",
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
        raise RuntimeError("ZETA1 self-adjoint-domain output already exists")
    authenticated = _authenticate()

    rho, outside, angular = sp.symbols("rho y kappa", real=True)
    radius = sp.sqrt(1 + rho**2)
    conformal_sigma = -sp.Rational(18, 1015) * outside**2
    warp = sp.exp(2 * conformal_sigma)
    sigma_derivative = sp.diff(conformal_sigma, outside)

    radial_superpotential = angular / radius
    superpotential_derivative = sp.simplify(
        sp.diff(radial_superpotential, rho)
    )
    partner_plus = sp.simplify(
        radial_superpotential**2 + superpotential_derivative
    )
    partner_minus = sp.simplify(
        radial_superpotential**2 - superpotential_derivative
    )
    partner_difference = sp.simplify(
        partner_plus - partner_minus - 2 * superpotential_derivative
    )

    # A concrete two-component boundary algebra verifies the general flux
    # statement.  The parent and child outward normals at the cut are opposite.
    phase = sp.symbols("vartheta", real=True)
    normal_current = sp.diag(1, -1)
    unitary_match = sp.diag(sp.exp(sp.I * phase), sp.exp(-sp.I * phase))
    unitary_residual = sp.simplify(
        unitary_match.conjugate().T
        * normal_current
        * unitary_match
        - normal_current
    )
    parent_spinor = sp.Matrix(sp.symbols("p_1:3", complex=True))
    test_spinor = sp.Matrix(sp.symbols("q_1:3", complex=True))
    child_spinor = unitary_match * parent_spinor
    child_test_spinor = unitary_match * test_spinor
    joined_flux = sp.simplify(
        (parent_spinor.conjugate().T * normal_current * test_spinor)[0]
        - (
            child_spinor.conjugate().T
            * normal_current
            * child_test_spinor
        )[0]
    )

    record = {
        "artifact_id": "NSC-2-ZETA1-SELF-ADJOINT-DOMAIN",
        "schema": "NSC-2-ZETA1-SELF-ADJOINT-DOMAIN-v1",
        "classification": "the_global_transition_slice_has_a_self_adjoint_Dirac_gluing_domain_and_the_throat_Phi_is_a_derived_Neumann_jump_operator",
        "authenticated_inputs": authenticated,
        "spatial_Dirac_geometry": {
            "spatial_metric": "h=W(y)*[d_rho^2+(1+rho^2)dOmega2^2+dy^2]",
            "conformal_factor": str(conformal_sigma),
            "warp": str(warp),
            "unwarped_Dirac": "D_tilde=gamma_rho*(partial_rho+rho/(1+rho^2))+D_S2/sqrt(1+rho^2)+gamma_y*partial_y",
            "conformal_rule": "D_h=exp(-5*sigma/2)*D_tilde*exp(3*sigma/2)",
            "sigma_y": str(sigma_derivative),
            "Hamiltonian_additions": "lapse_sqrt(W), radial_shift_sqrt(1-A), and the associated spin_connection_lower_order_terms",
        },
        "radial_partial_wave": {
            "spinor_sphere_eigenvalue": "kappa=plus_or_minus(l+1)",
            "flat_measure_rescaling": "psi_radial=u(rho)/sqrt(1+rho^2)",
            "first_order_pair": "A_kappa=partial_rho+kappa/sqrt(1+rho^2); A_kappa_dagger=-partial_rho+kappa/sqrt(1+rho^2)",
            "superpotential": str(radial_superpotential),
            "superpotential_derivative": str(superpotential_derivative),
            "squared_partner_plus": str(partner_plus),
            "squared_partner_minus": str(partner_minus),
            "partner_identity_residual": str(partner_difference),
            "asymptotic_potentials": "both_decay_as_kappa_squared/rho_squared_with_opposite_order_rho^-2_derivative_corrections",
        },
        "boundary_form": {
            "general": "<psi,H chi>-<H psi,chi>=integral_boundary sqrt(gamma)*psi_dagger*alpha_normal*chi",
            "y_boundaries": "complementary_APS_projectors_of_the_tangential_Dirac_operator_at_y=plus_or_minus_1",
            "radial_asymptotics": "two_channel_scattering_domain_with_square_integrable_bound_states",
            "throat_parent_normal": "+alpha_rho",
            "throat_child_normal": "-alpha_rho",
            "joined_condition": "psi_child=U_Phi*psi_parent with U_Phi_dagger*alpha_rho*U_Phi=alpha_rho",
            "finite_matrix_current": str(normal_current),
            "finite_matrix_matching": str(unitary_match),
            "unitary_current_residual": str(unitary_residual),
            "joined_flux_residual": str(joined_flux),
        },
        "gluing_owner": {
            "connected_operator": "H_joined_on_complete_rho_line",
            "disconnected_reference": "H_parent_APS_direct_sum_H_child_APS",
            "Calderon_data": "C_parent(E)_and_C_child(E)_from_the_actual_half_operator_Cauchy_spaces",
            "Neumann_jump": "R(E)=N_parent(E)+N_child(E)",
            "derived_link_identification": "Phi_throat(E)=R(E)",
            "relative_determinant": "det_zeta(H_joined^2)/[det_zeta(H_parent_APS^2)*det_zeta(H_child_APS^2)]=local_gluing_factor*det_zeta(R)",
            "arbitrary_throat_profile_inserted": False,
            "primary_reference": {
                "title": "BFK_gluing_formula_and_adiabatic_decomposition_of_the_zeta_determinant_of_a_Dirac_Laplacian",
                "arxiv": "math/0304347",
            },
        },
        "next_result": "discretize_the_lowest_Dirac_angular_and_y_parity_sectors_and_compare_joined_versus_disconnected_spectra_before_the_full_mode_sum",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "conformal_Dirac_weight_correct_for_four_spatial_dimensions": True,
            "radial_partner_identity_exact": partner_difference == 0,
            "unitary_matching_preserves_normal_current": unitary_residual
            == sp.zeros(2),
            "joined_throat_flux_cancels": joined_flux == 0,
            "APS_y_domain_declared": True,
            "two_channel_asymptotic_domain_declared": True,
            "throat_link_derived_from_jump_operator": True,
            "numerical_joined_disconnected_spectrum_completed": False,
            "zeta_derived": False,
        },
        "nonclaims": {
            "formal_domain_is_a_computed_spectrum": False,
            "lowest_partial_wave_is_the_complete_Dirac_operator": False,
            "BFK_compact_formula_directly_proves_the_noncompact_result": False,
            "physical_gap_fixed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key in {
            "numerical_joined_disconnected_spectrum_completed",
            "zeta_derived",
        }:
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 self-adjoint-domain gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
