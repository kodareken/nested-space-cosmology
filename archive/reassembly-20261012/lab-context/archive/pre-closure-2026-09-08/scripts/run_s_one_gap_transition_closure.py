#!/usr/bin/env python3
"""Match the spectral gap, doubled vacuum, and black-boundary scale.

The doubled-FLRW action expresses its effective vacuum through the same Phi
that opens the two-sheet mass gap.  The independently derived black geometry
expresses the room vacuum as 108/(1015 L_star^2), while the warped-resolution
result identifies L_star as the spectral crossover.  Requiring those to be
the same projection fixes Phi/Lambda up to the unresolved sign convention for
moving the vacuum-form tensor across the metric equation.

Both sign routes are retained.  They straddle the spectral wall symmetrically,
and their exact separation is the transition-curvature magnitude.  The full
exponential relative kernel is evaluated at both derived values; no observed
cosmological number or particle mass is used.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from math import exp, pi, sqrt
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-gap-transition-closure.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-doubled-flrw-closure.json",
    ROOT / "results/nsc-1-s-one-dark-black-invariant.json",
    ROOT / "results/nsc-1-s-one-warped-resolution-bind.json",
    ROOT / "results/nsc-1-s-one-full-exponential-relative-kernel.json",
    ROOT / "results/nsc-1-s-one-local-two-sheet-anomaly.json",
)
KERNEL_RUNNER = ROOT / "scripts/run_s_one_full_exponential_relative_kernel.py"


def _load_kernel_owner():
    specification = importlib.util.spec_from_file_location(
        "_nsc_full_relative_kernel_owner", KERNEL_RUNNER
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load full relative-kernel owner")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _authenticate() -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    authenticated: list[dict[str, object]] = []
    records: dict[str, dict[str, object]] = {}
    for path in (*INPUTS, KERNEL_RUNNER):
        raw = path.read_bytes()
        if path.suffix == ".json":
            item = json.loads(raw)
            if item.get("terminal") is not True:
                raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
            artifact_id = str(item["artifact_id"])
            records[artifact_id] = item
        else:
            artifact_id = None
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": artifact_id,
            }
        )
    return authenticated, records


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("gap-transition closure output already exists")
    authenticated, records = _authenticate()

    dark_black = records["NSC-1-S-ONE-DARK-BLACK-INVARIANT"]
    transition = Fraction(
        dark_black["one_number"]["transition_curvature_b2_over_F"]
    )
    vacuum_magnitude = Fraction(
        dark_black["one_number"]["local_vacuum_form_coefficient"]
    )
    if vacuum_magnitude != -6 * transition:
        raise RuntimeError("dark/black vacuum-transition identity changed")

    # c=1 is fixed by the exponential profile.  With Lambda=1/L_star, the
    # doubled identity Lambda_effective=12*(Lambda^2-kappa*Phi^2) gives
    # mu^2=(Phi/Lambda)^2=1-s*vacuum_magnitude/12 for kappa=+1.
    source_side_squared = 1 - vacuum_magnitude / 12
    geometric_left_squared = 1 + vacuum_magnitude / 12
    midpoint = (source_side_squared + geometric_left_squared) / 2
    separation = geometric_left_squared - source_side_squared

    if source_side_squared != Fraction(1006, 1015):
        raise RuntimeError("source-side gap candidate changed")
    if geometric_left_squared != Fraction(1024, 1015):
        raise RuntimeError("geometric-left gap candidate changed")

    kernel_owner = _load_kernel_owner()
    external_ratios = (0.0, 0.5, 1.0, 1.5, 2.0)
    candidate_rows: dict[str, list[dict[str, float | int]]] = {}
    for name, squared in (
        ("source_side_positive", source_side_squared),
        ("geometric_left_negative", geometric_left_squared),
    ):
        mass_ratio = sqrt(float(squared))
        candidate_rows[name] = [
            kernel_owner._one_kernel((64, mass_ratio, external_ratio))
            for external_ratio in external_ratios
        ]

    analytic_values = {
        name: 4.0
        * (float(squared) + 2.0)
        * exp(-float(squared))
        / pi**2
        for name, squared in (
            ("source_side_positive", source_side_squared),
            ("geometric_left_negative", geometric_left_squared),
        )
    }
    zero_mode_errors = {
        name: abs(rows[0]["normalized_kernel"] - analytic_values[name])
        / analytic_values[name]
        for name, rows in candidate_rows.items()
    }
    minimum_kernel = min(
        float(row["normalized_kernel"])
        for rows in candidate_rows.values()
        for row in rows
    )
    maximum_error = max(
        float(row["finite_difference_error"])
        for rows in candidate_rows.values()
        for row in rows
    )

    mu_squared, signed_vacuum = sp.symbols(
        "mu_squared v_signed", real=True
    )
    crossover_squared = sp.symbols("zeta", positive=True, real=True)
    general_matching_equation = sp.Eq(
        signed_vacuum, 12 * crossover_squared * (1 - mu_squared)
    )
    general_matching_solution = sp.solve(
        general_matching_equation, mu_squared
    )[0]
    matching_equation = general_matching_equation.subs(crossover_squared, 1)
    matching_solution = sp.solve(matching_equation, mu_squared)[0]

    record = {
        "artifact_id": "NSC-1-S-ONE-GAP-TRANSITION-CLOSURE",
        "schema": "NSC-1-S-ONE-GAP-TRANSITION-CLOSURE-v1",
        "classification": "the_one_scale_cross_projection_branch_fixes_two_boundary_gap_routes_that_straddle_the_spectral_cutoff",
        "authenticated_inputs": authenticated,
        "cross_projection_equation": {
            "spectral_resolution": "Lambda=1/L_star",
            "resolution_source": "the_strongest_one_scale_fixed_point_identifies_the_warped_boundary_width_with_the_inverse_spectral_cutoff",
            "general_dimensionless_crossover": "zeta=(Lambda*L_star)^2",
            "general_matching_equation": str(general_matching_equation),
            "general_matching_solution": str(general_matching_solution),
            "one_scale_branch": "zeta=1",
            "doubled_vacuum": "Lambda_effective*L_star^2=12*(1-(Phi*L_star)^2)",
            "black_geometry_vacuum_magnitude": str(vacuum_magnitude),
            "transition_curvature": str(transition),
            "vacuum_transition_identity": "v_room=-6*c_transition",
            "symbolic_matching_equation": str(matching_equation),
            "symbolic_solution": str(matching_solution),
            "observational_value_used": False,
        },
        "derived_gap_candidates": {
            "source_side_positive": {
                "signed_vacuum": str(vacuum_magnitude),
                "Phi_squared_over_Lambda_squared": str(source_side_squared),
                "Phi_over_Lambda": sqrt(float(source_side_squared)),
                "distance_squared_below_wall": str(1 - source_side_squared),
            },
            "geometric_left_negative": {
                "signed_vacuum": str(-vacuum_magnitude),
                "Phi_squared_over_Lambda_squared": str(geometric_left_squared),
                "Phi_over_Lambda": sqrt(float(geometric_left_squared)),
                "distance_squared_above_wall": str(
                    geometric_left_squared - 1
                ),
            },
            "exact_midpoint": str(midpoint),
            "exact_squared_separation": str(separation),
            "separation_equals_transition_magnitude": separation
            == abs(transition),
            "interpretation": "the_two_field_equation_sign_routes_are_symmetric_about_the_spectral_wall_and_the_transition_curvature_is_their_exact_split",
        },
        "full_exponential_evaluation": {
            "external_momentum_over_cutoff": list(external_ratios),
            "resolution": 64,
            "candidate_rows": candidate_rows,
            "analytic_zero_mode_kernel": analytic_values,
            "zero_mode_relative_errors": zero_mode_errors,
            "minimum_candidate_kernel": minimum_kernel,
            "maximum_finite_difference_error": maximum_error,
            "all_candidate_kernels_positive": minimum_kernel > 0.0,
        },
        "causal_result": {
            "continuous_Phi_over_Lambda_freedom_removed_in_one_scale_branch": True,
            "one_scale_identification_still_requires_full_operator_validation": True,
            "remaining_ambiguity": "one_discrete_field_equation_orientation_sign",
            "selection_rule": "the_curved_retarded_parent_child_equation_must_choose_the_sign_without_observational_fitting",
            "same_Phi_roles": [
                "boundary_scale_mass_gap",
                "visible_outside_self_energy",
                "relative_metric_interaction",
                "vacuum_shift",
            ],
            "next_equation": "evaluate_both_gap_candidates_in_the_curved_retarded_two_sheet_kernel_and_keep_only_a_causal_positive_orientation_then_test_delta_Gamma_over_delta_Phi",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS) + 1,
            "one_scale_fixed_point_branch_explicit": True,
            "vacuum_transition_identity_exact": vacuum_magnitude
            == -6 * transition,
            "source_side_candidate_exact": source_side_squared
            == Fraction(1006, 1015),
            "geometric_left_candidate_exact": geometric_left_squared
            == Fraction(1024, 1015),
            "candidates_center_on_spectral_wall": midpoint == 1,
            "candidate_split_is_transition_curvature": separation
            == abs(transition),
            "both_candidates_full_kernel_positive": minimum_kernel > 0.0,
            "zero_modes_match_analytic_below_5e_5": max(
                zero_mode_errors.values()
            )
            < 5.0e-5,
            "finite_difference_error_below_1e_5": maximum_error < 1.0e-5,
            "Lorentzian_orientation_selected": False,
            "stationary_gap_equation_solved": False,
            "conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "either_sign_selected_from_observation": False,
            "boundary_gap_is_the_electron_mass": False,
            "Lambda_times_L_star_independently_predicted": False,
            "Euler_Lagrange_stationarity_proved": False,
            "curved_Lorentzian_causality_proved": False,
            "observed_cosmological_constant_predicted": False,
            "Xi_NSC_predicted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key not in {
            "Lorentzian_orientation_selected",
            "stationary_gap_equation_solved",
            "conviction_proof_seed_completed",
        } and value is not True:
            raise RuntimeError(f"gap-transition closure failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
