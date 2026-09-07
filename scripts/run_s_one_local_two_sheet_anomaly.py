#!/usr/bin/env python3
"""Derive the first local two-sheet resolution and shared-link spectrum.

The published Higgs--dilaton anomaly already owns a local *common* Weyl
rescaling.  This calculation isolates the genuinely new parent/child degree of
freedom: the relative sheet resolution.  It then proves, on the smallest
Dirac symbol that retains the off-diagonal link, that the same ``Phi`` opens a
fermion mass gap and appears in the visible Schur-complement self-energy.

This is a controlled spectral-mode theorem.  It deliberately does not promote
the result to a complete curved Lorentzian anomaly, an electron spectrum, or a
black-to-child solution.
"""

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


OUTPUT = ROOT / "results/nsc-1-s-one-local-two-sheet-anomaly.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-invariant-partition-bind.json",
    ROOT / "results/nsc-1-s-one-one-operator-bind.json",
    ROOT / "results/nsc-1-s-one-transition-link-coefficient.json",
)


def _authenticate_inputs() -> list[dict[str, object]]:
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
    return authenticated


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("local two-sheet anomaly output already exists")

    authenticated = _authenticate_inputs()

    momentum, link, mean, relative, pole, cutoff = sp.symbols(
        "p Phi varphi_bar delta z Lambda", positive=True, real=True
    )

    # A minimal Euclidean Dirac symbol.  sigma_3 supplies the oriented local
    # propagation and sigma_1 is the finite off-diagonal sheet/link field.
    unscaled = sp.Matrix([[momentum, link], [link, -momentum]])
    scale = sp.diag(
        sp.exp(-(mean + relative) / 2),
        sp.exp(-(mean - relative) / 2),
    )
    scaled = sp.simplify(scale * unscaled * scale)

    expected_scaled = sp.Matrix(
        [
            [
                momentum * sp.exp(-mean - relative),
                link * sp.exp(-mean),
            ],
            [
                link * sp.exp(-mean),
                -momentum * sp.exp(-mean + relative),
            ],
        ]
    )
    scaled_residual = sp.simplify(scaled - expected_scaled)

    determinant = sp.factor(scaled.det())
    trace_square = sp.factor(sp.trace(scaled**2))
    expected_determinant = -(momentum**2 + link**2) * sp.exp(-2 * mean)
    expected_trace_square = sp.exp(-2 * mean) * (
        2 * link**2
        + momentum**2 * sp.exp(2 * relative)
        + momentum**2 * sp.exp(-2 * relative)
    )

    spectral_radius = sp.sqrt(momentum**2 + link**2)
    common_scaled = scaled.subs(relative, 0)
    common_characteristic = sp.factor(
        (pole * sp.eye(2) - common_scaled).det()
    )
    expected_characteristic = sp.factor(
        pole**2 - sp.exp(-2 * mean) * spectral_radius**2
    )

    resolvent = sp.simplify((pole * sp.eye(2) - unscaled).inv())
    visible_resolvent = sp.factor(resolvent[0, 0])
    visible_inverse = sp.factor(1 / visible_resolvent)
    expected_visible_inverse = sp.factor(
        pole - momentum - link**2 / (pole + momentum)
    )
    positive_residue = sp.simplify(
        sp.limit(
            (pole - spectral_radius) * visible_resolvent,
            pole,
            spectral_radius,
        )
    )
    negative_residue = sp.simplify(
        sp.limit(
            (pole + spectral_radius) * visible_resolvent,
            pole,
            -spectral_radius,
        )
    )

    # The exact heat trace is the controlled finite-momentum kernel needed by
    # the anomaly integral.  Its relative-resolution Hessian proves that the
    # relative sheet mode is invisible at p=0 but is dynamically seen by every
    # propagating mode.  The sign changes at the spectral cutoff, so the
    # bosonic heat trace alone cannot certify the physical relative mode.
    sheet_radius = sp.sqrt(
        momentum**2 * sp.cosh(relative) ** 2 + link**2
    )
    eigenvalue_plus = sp.exp(-mean) * (
        -momentum * sp.sinh(relative) + sheet_radius
    )
    eigenvalue_minus = sp.exp(-mean) * (
        -momentum * sp.sinh(relative) - sheet_radius
    )
    heat_trace = sp.exp(-(eigenvalue_plus / cutoff) ** 2) + sp.exp(
        -(eigenvalue_minus / cutoff) ** 2
    )
    relative_heat_hessian = sp.factor(
        sp.diff(heat_trace, relative, 2).subs({mean: 0, relative: 0})
    )
    expected_relative_heat_hessian = sp.factor(
        8
        * momentum**2
        * (momentum**2 + link**2 - cutoff**2)
        * sp.exp(-(momentum**2 + link**2) / cutoff**2)
        / cutoff**4
    )

    below = sp.simplify(
        relative_heat_hessian.subs(
            {momentum: sp.Rational(1, 4), link: sp.Rational(1, 2), cutoff: 1}
        )
    )
    at_cutoff = sp.simplify(
        relative_heat_hessian.subs(
            {momentum: sp.Rational(3, 5), link: sp.Rational(4, 5), cutoff: 1}
        )
    )
    above = sp.simplify(
        relative_heat_hessian.subs({momentum: 1, link: 1, cutoff: 1})
    )

    # In the individual sheet variables, the mixed second heat moment is the
    # direct local witness that the boundary field couples the two resolutions.
    parent_scale, child_scale = sp.symbols(
        "varphi_parent varphi_child", real=True
    )
    sheet_scaled = sp.diag(
        sp.exp(-parent_scale / 2), sp.exp(-child_scale / 2)
    ) * unscaled * sp.diag(
        sp.exp(-parent_scale / 2), sp.exp(-child_scale / 2)
    )
    mixed_second_moment = sp.simplify(
        sp.diff(
            sp.trace(sheet_scaled**2), parent_scale, child_scale
        ).subs({parent_scale: 0, child_scale: 0})
    )

    record = {
        "artifact_id": "NSC-1-S-ONE-LOCAL-TWO-SHEET-ANOMALY",
        "schema": "NSC-1-S-ONE-LOCAL-TWO-SHEET-ANOMALY-v1",
        "classification": "local_two_sheet_resolution_split_proves_one_Phi_is_both_mass_gap_and_visible_outside_self_energy",
        "authenticated_inputs": authenticated,
        "imported_local_common_anomaly": {
            "source": "Kurkov_and_Lizzi_Higgs_Dilaton_Lagrangian_from_Spectral_Regularization",
            "arxiv": "1210.2663v1",
            "doi": "10.1142/S0217732312502033",
            "published_equations": ["Dstruct", "TRANSzero", "Scoll2", "finalansw"],
            "established_scope": "one_common_local_Weyl_field_with_Higgs_gauge_curvature_and_derivative_terms",
            "not_established_by_source": "two_independent_parent_child_resolution_fields_or_their_relative_mode",
            "known_result_recomputed": False,
        },
        "exact_local_operator_identity": {
            "sheet_fields": "varphi_parent=varphi_bar+delta; varphi_child=varphi_bar-delta",
            "local_diagonal_rule": "exp(-varphi_i/2) D_i exp(-varphi_i/2)=exp(-varphi_i)[D_i-c(dvarphi_i)/2]",
            "off_diagonal_rule": "exp(-varphi_parent/2) Phi exp(-varphi_child/2)=exp(-varphi_bar) Phi",
            "meaning": "the_mean_resolution_scales_the_link_while_the_relative_resolution_changes_the_two_diagonal_propagators",
            "controlled_symbol": str(unscaled),
            "scaled_symbol": str(scaled),
            "scaled_symbol_residual": str(scaled_residual),
            "determinant": str(determinant),
            "trace_D_squared": str(trace_square),
        },
        "same_link_particle_and_outside_identity": {
            "common_scale_characteristic": str(common_characteristic),
            "poles": [
                "-exp(-varphi_bar)*sqrt(p^2+Phi^2)",
                "+exp(-varphi_bar)*sqrt(p^2+Phi^2)",
            ],
            "rest_mass_gap": "abs(Phi)*exp(-varphi_bar)",
            "visible_resolvent": str(visible_resolvent),
            "visible_inverse": str(visible_inverse),
            "outside_self_energy": "-Phi^2/(z+p)",
            "positive_pole_residue": str(positive_residue),
            "negative_pole_residue": str(negative_residue),
            "zero_momentum_residues": [
                str(sp.simplify(positive_residue.subs(momentum, 0))),
                str(sp.simplify(negative_residue.subs(momentum, 0))),
            ],
            "zero_link_control": "Phi=0_removes_both_the_mass_gap_and_the_outside_self_energy",
        },
        "local_heat_kernel_result": {
            "eigenvalue_plus": str(eigenvalue_plus),
            "eigenvalue_minus": str(eigenvalue_minus),
            "heat_trace": str(heat_trace),
            "relative_resolution_hessian_at_fixed_room": str(
                relative_heat_hessian
            ),
            "below_cutoff_sample": str(below),
            "at_cutoff_sample": str(at_cutoff),
            "above_cutoff_sample": str(above),
            "mixed_parent_child_second_moment": str(mixed_second_moment),
            "interpretation": "propagating_modes_see_the_relative_resolution_and_the_heat_trace_curvature_changes_sign_exactly_at_p2_plus_Phi2_equals_Lambda2",
        },
        "first_causal_obstruction": {
            "kind": "missing_relative_two_sheet_anomaly_and_stationary_link_amplitude",
            "common_local_anomaly_available": True,
            "relative_local_anomaly_derived_in_full_continuum": False,
            "Phi_over_Lambda_fixed_by_current_equations": False,
            "why_proof_seed_cannot_yet_emit_Xi_NSC": "the_current_transition_result_fixes_a_normalized_link_curvature_b2_over_F_not_the_stationary_spectral_gap_Phi_over_Lambda",
            "next_equation": "derive_the_matrix_valued_local_anomaly_for_varphi_bar_I_plus_delta_sigma3_then_solve_delta_Gamma_one_over_delta(varphi_bar,delta,Phi)=0",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "published_local_common_anomaly_bound_not_recomputed": True,
            "scaled_two_sheet_symbol_exact": scaled_residual == sp.zeros(2),
            "determinant_exact": sp.simplify(determinant - expected_determinant)
            == 0,
            "second_heat_moment_exact": sp.simplify(
                trace_square - expected_trace_square
            )
            == 0,
            "common_scale_dispersion_exact": sp.simplify(
                common_characteristic - expected_characteristic
            )
            == 0,
            "same_Phi_in_mass_gap_and_Schur_self_energy": sp.simplify(
                visible_inverse - expected_visible_inverse
            )
            == 0,
            "positive_energy_residue_positive": bool(positive_residue > 0),
            "zero_momentum_residues_sum_to_one": sp.simplify(
                positive_residue.subs(momentum, 0)
                + negative_residue.subs(momentum, 0)
                - 1
            )
            == 0,
            "relative_heat_hessian_exact": sp.simplify(
                relative_heat_hessian - expected_relative_heat_hessian
            )
            == 0,
            "relative_heat_hessian_changes_at_cutoff": bool(
                below < 0 and at_cutoff == 0 and above > 0
            ),
            "link_couples_sheet_heat_moments": mixed_second_moment
            == 2 * link**2,
            "full_conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "complete_curved_two_sheet_anomaly_derived": False,
            "reflection_positivity_completed": False,
            "electron_mass_predicted": False,
            "stationary_Phi_over_Lambda_predicted": False,
            "black_child_solution_sourced_by_this_Phi": False,
            "Xi_NSC_predicted": False,
            "observed_universe_identified": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "full_conviction_proof_seed_completed" and value is not True:
            raise RuntimeError(f"local two-sheet anomaly gate failed: {key}")

    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
