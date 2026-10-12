#!/usr/bin/env python3
"""Bind the exact doubled-FLRW spectral link and test its leading stability.

The published doubled spectral action uses one off-diagonal ``Phi`` to set the
interaction between two FLRW metrics and to shift their effective vacuum
curvature.  The preceding project result proves that the same ``Phi`` is the
rest mass gap and visible Schur-complement self-energy of the controlled
two-sheet Dirac symbol.

This calculation composes those equations without identifying unrelated
transition-curvature data.  It also tests the generic de Sitter relative-mode
exponents of the published leading (a0+a2) truncation.  The resulting
incompatibility is a measured authority boundary for that truncation and
selects the full exponential spectral kernel as the next calculation.
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


OUTPUT = ROOT / "results/nsc-1-s-one-doubled-flrw-closure.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-local-two-sheet-anomaly.json",
    ROOT / "results/nsc-1-s-one-spectral-profile-closure.json",
    ROOT / "results/nsc-1-s-one-invariant-partition-bind.json",
)


def _authenticate() -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    authenticated: list[dict[str, object]] = []
    records: dict[str, dict[str, object]] = {}
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        records[str(item["artifact_id"])] = item
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item["artifact_id"],
            }
        )
    return authenticated, records


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("doubled-FLRW closure output already exists")

    authenticated, records = _authenticate()
    profile = records["NSC-1-S-ONE-SPECTRAL-PROFILE-CLOSURE"]
    moments = profile["unique_normalized_profile"][
        "first_five_dimensionless_moments"
    ]
    leading_moment_ratio = sp.Rational(moments["1"]) / sp.Rational(
        moments["0"]
    )

    cutoff_squared, link_squared, coefficient = sp.symbols(
        "Lambda_squared Phi_squared c", positive=True, real=True
    )
    kappa = sp.symbols("kappa", nonzero=True, real=True)
    lambda_effective = sp.symbols(
        "Lambda_effective", positive=True, real=True
    )
    alpha = sp.symbols("alpha", real=True)

    published_lambda = sp.factor(
        sp.Rational(12, 1)
        / coefficient
        * (cutoff_squared - coefficient * kappa * link_squared)
    )
    published_alpha = 12 * kappa * link_squared
    shared_identity = sp.factor(published_lambda + published_alpha)
    expected_shared_identity = 12 * cutoff_squared / coefficient

    fixed_profile_lambda = sp.simplify(
        published_lambda.subs(coefficient, leading_moment_ratio)
    )
    fixed_profile_alpha = sp.simplify(published_alpha)
    fixed_profile_identity = sp.simplify(
        fixed_profile_lambda + fixed_profile_alpha
    )
    positive_gap_alpha = sp.simplify(published_alpha.subs(kappa, 1))

    # The doubled-FLRW de Sitter solution writes Lambda_effective=6*hubble^2.
    # Its generic relative perturbations have the published exponents below.
    hubble = sp.symbols("lambda", positive=True, real=True)
    discriminant = sp.sqrt(21 * hubble**2 + 2 * alpha)
    relative_scale_exponents = (
        sp.factor((-hubble + discriminant) / 2),
        sp.factor((-hubble - discriminant) / 2),
    )
    relative_lapse_exponents = (
        sp.factor((-3 * hubble + discriminant) / 2),
        sp.factor((-3 * hubble - discriminant) / 2),
    )
    generic_scale_decay_bound = -10 * hubble**2
    generic_lapse_decay_bound = -6 * hubble**2
    scale_decay_in_lambda_effective = sp.simplify(
        generic_scale_decay_bound.subs(
            hubble**2, lambda_effective / 6
        )
    )
    lapse_decay_in_lambda_effective = sp.simplify(
        generic_lapse_decay_bound.subs(
            hubble**2, lambda_effective / 6
        )
    )

    # For c>0 and Lambda_cutoff^2>0 the published identity enforces
    # alpha > -Lambda_effective.  Generic decay instead needs
    # alpha < -(5/3)Lambda_effective.  At positive Lambda_effective those
    # half-lines do not intersect.
    spectral_lower_offset = sp.simplify(
        alpha + lambda_effective
    )
    stability_upper_offset = sp.simplify(
        alpha + sp.Rational(5, 3) * lambda_effective
    )
    incompatibility_gap = sp.simplify(
        (-lambda_effective)
        - (-sp.Rational(5, 3) * lambda_effective)
    )

    # Direct examples for both allowed grading signs make the contradiction
    # concrete without selecting a physical value of Phi/Lambda.
    plus_example = {
        cutoff_squared: 1,
        link_squared: sp.Rational(1, 4),
        kappa: 1,
        coefficient: leading_moment_ratio,
    }
    minus_example = {
        cutoff_squared: 1,
        link_squared: sp.Rational(1, 4),
        kappa: -1,
        coefficient: leading_moment_ratio,
    }

    record = {
        "artifact_id": "NSC-1-S-ONE-DOUBLED-FLRW-CLOSURE",
        "schema": "NSC-1-S-ONE-DOUBLED-FLRW-CLOSURE-v1",
        "classification": "one_Phi_is_mass_gap_relative_metric_coupling_and_vacuum_shift_but_the_leading_doubled_FLRW_truncation_cannot_generically_stabilize_positive_de_Sitter",
        "authenticated_inputs": authenticated,
        "imported_doubled_geometry": {
            "operator": "D_doubled=[[D_1,gamma*Phi],[gamma*Phi_star,D_2]]",
            "source_action": {
                "title": "On_stability_of_Friedmann_Lemaitre_Robertson_Walker_solutions_in_doubled_geometries",
                "authors": ["Arkadiusz_Bochniak", "Andrzej_Sitarz"],
                "arxiv": "2012.06401v2",
                "doi": "10.1103/PhysRevD.103.044041",
                "equations": ["S-k1", "S-k1L", "EL01", "EL03", "Eq2", "E1"],
            },
            "source_interaction": {
                "title": "Spectral_interaction_between_universes",
                "authors": ["Arkadiusz_Bochniak", "Andrzej_Sitarz"],
                "arxiv": "2201.03839v1",
                "doi": "10.1088/1475-7516/2022/04/055",
                "equations": ["Sh1h2", "eq:bimW"],
            },
            "established_scope": "leading_two_spectral_terms_for_doubled_FLRW_and_perturbative_metric_interaction_about_flat_Euclidean_geometry_with_a_Lorentzian_FLRW_continuation",
            "known_results_recomputed": False,
        },
        "fixed_exponential_profile": {
            "profile": "exp(-D_squared/Lambda_squared)",
            "zeroth_moment": moments["0"],
            "first_moment": moments["1"],
            "leading_Wodzicki_ratio_c": str(leading_moment_ratio),
            "no_profile_fit": True,
        },
        "exact_same_Phi_closure": {
            "mass_gap_from_predecessor": "m_gap_squared=abs(Phi)^2",
            "doubled_Dirac_square_zero_order": "a_0=kappa*F^2",
            "positive_Euclidean_gap_branch": "kappa=+1",
            "relative_metric_coupling": "alpha=12*kappa*abs(Phi)^2",
            "effective_vacuum": "Lambda_effective=(12/c)*(Lambda_squared-c*kappa*abs(Phi)^2)",
            "eliminated_Phi_identity": "Lambda_effective+alpha=12*Lambda_squared/c",
            "symbolic_identity_residual": str(
                sp.simplify(shared_identity - expected_shared_identity)
            ),
            "fixed_profile_effective_vacuum": str(fixed_profile_lambda),
            "fixed_profile_interaction": str(fixed_profile_alpha),
            "fixed_profile_identity": str(fixed_profile_identity),
            "meaning": "the_same_off_diagonal_spectral_gap_shifts_the_vacuum_and_couples_the_relative_parent_child_metric",
        },
        "relative_de_Sitter_stability": {
            "background": "a(t)=a0*exp(lambda*t); Lambda_effective=6*lambda^2",
            "relative_scale_exponents": [
                str(value) for value in relative_scale_exponents
            ],
            "relative_lapse_exponents": [
                str(value) for value in relative_lapse_exponents
            ],
            "generic_scale_decay_requires": "alpha<-10*lambda^2=-(5/3)*Lambda_effective",
            "generic_lapse_decay_requires": "alpha<-6*lambda^2=-Lambda_effective",
            "strongest_required_bound": str(scale_decay_in_lambda_effective),
            "spectral_identity_requires_for_positive_cutoff_and_c": "alpha>-Lambda_effective",
            "positive_gap_branch_additionally_has": "alpha=12*abs(Phi)^2>0",
            "kappa_minus_one_caveat": "it_reverses_the_F_squared_term_in_the_Euclidean_Dirac_square_and_is_not_the_positive_mass_gap_branch",
            "incompatibility_gap": str(incompatibility_gap),
            "special_decaying_initial_subspace_exists": True,
            "special_initial_subspace_is_generic_stability": False,
        },
        "signed_examples": {
            "kappa_plus_one": {
                "Phi_squared_over_Lambda_squared": "1/4",
                "Lambda_effective_over_Lambda_squared": str(
                    sp.simplify(published_lambda.subs(plus_example))
                ),
                "alpha_over_Lambda_squared": str(
                    sp.simplify(published_alpha.subs(plus_example))
                ),
            },
            "kappa_minus_one": {
                "Phi_squared_over_Lambda_squared": "1/4",
                "Lambda_effective_over_Lambda_squared": str(
                    sp.simplify(published_lambda.subs(minus_example))
                ),
                "alpha_over_Lambda_squared": str(
                    sp.simplify(published_alpha.subs(minus_example))
                ),
            },
        },
        "causal_result": {
            "new_positive_connection": "one_Phi_now_has_three_exact_roles_mass_gap_relative_metric_interaction_and_effective_vacuum_shift",
            "measured_first_obstruction": "the_leading_a0_plus_a2_doubled_FLRW_truncation_has_no_parameter_overlap_between_positive_de_Sitter_and_generic_relative_mode_decay",
            "responsible_layer": "finite_leading_spectral_truncation_of_the_relative_metric_kernel",
            "not_responsible": [
                "the_exact_two_sheet_mass_and_Schur_identity",
                "the_full_exponential_spectral_action",
                "the_nonlinear_parent_child_hypothesis",
            ],
            "next_equation": "compute_the_full_exponential_two_sheet_relative_kernel_K_minus(omega,k;Phi/Lambda)_on_the_symmetric_background_and_require_a_positive_retarded_spectrum",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "exponential_profile_fixes_c_to_one": leading_moment_ratio == 1,
            "same_Phi_elimination_exact": sp.simplify(
                shared_identity - expected_shared_identity
            )
            == 0,
            "mass_interaction_ratio_exact": sp.simplify(
                published_alpha / link_squared - 12 * kappa
            )
            == 0,
            "positive_Euclidean_gap_implies_positive_alpha": bool(
                positive_gap_alpha > 0
            ),
            "generic_scale_decay_bound_exact": scale_decay_in_lambda_effective
            == -sp.Rational(5, 3) * lambda_effective,
            "generic_lapse_decay_bound_exact": lapse_decay_in_lambda_effective
            == -lambda_effective,
            "positive_de_Sitter_stability_intersection_empty": bool(
                incompatibility_gap > 0
            ),
            "full_exponential_relative_kernel_tested": False,
            "conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "Phi_is_the_electron_mass": False,
            "Phi_over_Lambda_predicted": False,
            "leading_truncation_is_the_full_spectral_theory": False,
            "Lorentzian_reflection_positivity_proved": False,
            "black_child_transition_solved": False,
            "Xi_NSC_predicted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key not in {
            "full_exponential_relative_kernel_tested",
            "conviction_proof_seed_completed",
        } and value is not True:
            raise RuntimeError(f"doubled-FLRW closure gate failed: {key}")

    # Keep these symbolic facts visible rather than relying on prose alone.
    if spectral_lower_offset != alpha + lambda_effective:
        raise RuntimeError("spectral lower-bound expression changed")
    if stability_upper_offset != alpha + sp.Rational(5, 3) * lambda_effective:
        raise RuntimeError("stability upper-bound expression changed")

    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
