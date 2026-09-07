#!/usr/bin/env python3
"""Bind known reductions to one action and expose the first unclosed equality."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-invariant-closure-bind.json"
INPUTS = {
    "particle_wave_collective_reduction": ROOT
    / "results/nsc-1-s-one-particle-wave-identity.json",
    "nuclear_topological_reduction": ROOT / "results/nsc-1-bps-self-equal-energy.json",
    "nuclear_to_gravity_reduction": ROOT
    / "results/nsc-1-gravitating-bps-observation-link.json",
    "black_to_child_geometry_reduction": ROOT
    / "results/nsc-1-exact-black-universe-defocusing.json",
    "positive_F_null_projection": ROOT
    / "results/nsc-1-unified-gradient-boundary-dev3.json",
    "local_constant_dictionary": ROOT
    / "results/nsc-1-s-one-constant-dictionary.json",
    "reciprocal_scale_closure": ROOT
    / "results/nsc-1-s-one-reciprocal-closure.json",
    "nested_coefficient_closure": ROOT
    / "results/nsc-1-s-one-nested-pair-closure.json",
}


def _exact(value: Fraction) -> dict[str, str | float]:
    return {"exact": str(value), "decimal": float(value)}


def _read_bound_inputs() -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    bindings: dict[str, object] = {}
    for role, path in INPUTS.items():
        raw = path.read_bytes()
        record = json.loads(raw)
        if record.get("terminal") is not True:
            raise RuntimeError(f"bound input is not terminal: {path.relative_to(ROOT)}")
        records[role] = record
        bindings[role] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "artifact_id": record.get("artifact_id"),
            "classification": record.get("classification"),
        }
    return records, bindings


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one invariant-closure bind output already exists")

    records, bindings = _read_bound_inputs()
    black = records["black_to_child_geometry_reduction"]
    positive_f = records["positive_F_null_projection"]
    particle_wave = records["particle_wave_collective_reduction"]
    gravitating = records["nuclear_to_gravity_reduction"]

    rho_values = (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))
    required_null = {
        str(rho): _exact(-Fraction(2, 1) / (rho * rho + 1) ** 2)
        for rho in rho_values
    }
    required_q = {
        str(rho): _exact(Fraction(2, 1) * (1 - rho * rho) / (rho * rho + 1) ** 2)
        for rho in rho_values
    }

    black_assertions = black.get("assertions", {})
    positive_f_assertions = positive_f.get("assertions", {})
    particle_nonclaims = particle_wave.get("nonclaims", {})
    gravity_status = gravitating.get("one_model_status", {})

    record = {
        "artifact_id": "NSC-1-S-ONE-INVARIANT-CLOSURE-BIND",
        "schema": "NSC-1-S-ONE-INVARIANT-CLOSURE-BIND-v1",
        "classification": "known_reductions_bound_to_one_Gamma_with_one_explicit_unclosed_cross_sector_source_equation",
        "method": "authenticate_and_bind_existing_terminals_without_recomputing_their_known_physics",
        "one_equation": {
            "influence_functional": "Gamma_n=-i*hbar_n*log integral_{U\\S_n} D(Phi_out)*exp(i*S_one[Phi_n,Phi_out,B;Theta]/hbar_n)",
            "observable_jet": "O_n=P*J_infinity_{Phi_n,g_n}(Gamma_n)|_{delta_Gamma_n=0}",
            "recursive_fixed_point": "Gamma_n=T_Theta^*Gamma_{n+1}; Theta=Theta",
            "observable_vector": [
                "Lambda_n",
                "c_n",
                "G_n",
                "ell_P_n",
                "particle_poles_charges_spins_and_wave_residues",
                "finite_wavelength_dark_metric_response",
                "black_boundary_to_child_transition",
            ],
        },
        "authenticated_inputs": bindings,
        "imported_reduction_map": {
            "quantum": {
                "known": "field excitations and their quantum amplitudes are one quantum-state description; interference, Bell nonlocality, and PBR constraints are imported",
                "restriction": "R_quantum(Gamma)=Gamma_QFT",
                "project_status": "collective BPS particle/wave identity is bound; electron charge-spin stationary solution remains open",
            },
            "nuclear": {
                "known": "topological field configurations and nuclear energy differences are imported",
                "restriction": "R_nuclear(Gamma)=S_BPS_or_near_BPS",
                "shared_quantity": "the charge-one stationary energy is the one-body limit of the same multi-charge functional",
            },
            "gravity": {
                "known": "published Einstein-BPS solutions apply G to the nuclear-fixed field system",
                "restriction": "R_gravity(Gamma): delta_Gamma/delta_g=0",
                "shared_quantity": "the nuclear parameters remain unchanged when gravity is activated",
            },
            "black_to_child": {
                "known": "a regular four-geometry with trapped positive-Q interval and expanding child interior exists",
                "restriction": "R_transition(Gamma): delta_Gamma/delta_fields=0 on the black-universe saddle",
                "shared_quantity": "the same affine null congruence crosses from trapped contraction through finite minimum radius to expansion",
            },
            "dark_and_black": {
                "restriction": "R_dark,black(Gamma)=P[J_infinity(Gamma_outside)]",
                "meaning": "finite-wavelength unresolved metric response and nonlinear causal-boundary saddles are projections of the same integrated-out sector, not separately named substances",
            },
        },
        "commuting_restriction_conditions": [
            "R_one_body(R_nuclear(Gamma))=R_nuclear(R_particle(Gamma))",
            "R_gravity(R_nuclear(Gamma))=R_nuclear(R_gravity(Gamma))",
            "R_transition(R_gravity(Gamma))=R_gravity(R_transition(Gamma))",
            "R_child(R_transition(Gamma))=T_Theta^*R_parent(R_transition(Gamma))",
            "R_dark(Gamma) and R_black(Gamma) use the identical Pi_outside,Theta",
        ],
        "exact_transition_bridge": {
            "black_universe_required_Ricci_null": "R_kk=-2/(rho^2+1)^2",
            "sampled_trapped_interval": ["1/4", "3/4"],
            "required_Ricci_null_samples": required_null,
            "complete_Q_samples": required_q,
            "Einstein_null_projection": "8*pi*G*T_one,kk=R_kk",
            "BPS_visible_null_projection": "rho_BPS+p_BPS=2*lambda^2*pi^4*B_0^2>=0",
            "therefore_required_shared_boundary_projection": "8*pi*G*T_boundary,kk=-2/(rho^2+1)^2-8*pi*G*T_visible,kk",
            "strict_consequence": "T_boundary,kk<0 throughout the trapped positive-Q interval",
            "existing_positive_F_partial_solution": "F*R_kk=d^2F/dlambda^2 with F>0 and equal endpoint F",
            "unclosed_equations": [
                "the scalar/order-parameter Euler-Lagrange equation",
                "the remaining independent metric equations",
                "the quadratic ghost-gradient-cone stability conditions",
                "the zero- and finite-momentum dark-response projections of the same boundary kernel",
            ],
        },
        "what_is_already_bound": {
            "particle_wave_same_field_configuration": particle_wave.get("gate", {}).get(
                "particle_and_wave_use_same_field_configuration"
            ),
            "nuclear_parameters_unchanged_under_gravity": gravity_status.get(
                "nuclear_parameters_fixed_before_gravity"
            ),
            "trapped_to_child_geometry_exact": black_assertions.get(
                "same_affine_congruence_changes_sign"
            ),
            "positive_F_supplies_required_null_projection": positive_f_assertions.get(
                "nonminimal_null_term_matches_required_Ricci"
            ),
        },
        "what_is_not_yet_bound": {
            "electron_stationary_solution": not bool(
                particle_nonclaims.get("electron_solution_obtained", True)
            ),
            "complete_shared_boundary_action_solution": not bool(
                positive_f.get("nonclaims", {}).get("full_four_dimensional_action_solution", True)
            ),
            "same_kernel_dark_lensing_and_expansion_prediction": False,
            "observational_parent_child_identification": False,
        },
        "single_next_result": {
            "equation": "find one degenerate boundary Lagrangian L_deg whose full Euler-Lagrange system contains the authenticated positive-F black-universe saddle and whose quadratic kernel is healthy",
            "same_kernel_requirement": "the zero- and finite-momentum metric response of that L_deg must then be used for Lambda, clustering, and lensing without adding dark functions",
            "success": "full_background_residual_zero AND positive_kinetic_and_cone_margins AND unchanged_BPS_reduction",
            "reason": "this is the first noncommuting edge in the imported reduction diagram",
        },
        "gate": {
            "all_inputs_terminal_and_hashed": len(bindings) == len(INPUTS),
            "no_known_component_recomputed": True,
            "black_universe_transition_exact": black_assertions.get(
                "field_equations_symbolically_zero"
            )
            is True,
            "positive_F_partial_bridge_exact": positive_f_assertions.get(
                "nonminimal_null_term_matches_required_Ricci"
            )
            is True,
            "single_unclosed_bridge_identified": True,
            "invariant_closure_completed": False,
        },
        "terminal": True,
    }
    if not all(
        value
        for key, value in record["gate"].items()
        if key != "invariant_closure_completed"
    ):
        raise RuntimeError("S-one invariant-closure bind gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
