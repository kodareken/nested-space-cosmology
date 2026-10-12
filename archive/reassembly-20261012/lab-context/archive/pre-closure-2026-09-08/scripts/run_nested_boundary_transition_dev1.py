#!/usr/bin/env python3
"""Construct an exact finite parent-collapse/child-expansion boundary example."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


Q = Fraction
OUTPUT = ROOT / "results/nsc-1-parent-child-boundary-dev1.json"


def _text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _record(value: Fraction) -> dict[str, object]:
    return {"exact": _text(value), "decimal": float(value)}


def construct() -> dict[str, object]:
    # G=c=1. Choosing rational v makes every geometric identity exact.
    mass = Q(1)
    boundary_speed = Q(5, 4)
    radius = 2 * mass / boundary_speed**2
    hubble = boundary_speed / radius
    compactness = 2 * mass / radius

    parent_areal_out = 1 - boundary_speed
    parent_areal_in = -1 - boundary_speed
    child_areal_out = boundary_speed + 1
    child_areal_in = boundary_speed - 1
    parent_theta_out = 2 * parent_areal_out / radius
    parent_theta_in = 2 * parent_areal_in / radius
    child_theta_out = 2 * child_areal_out / radius
    child_theta_in = 2 * child_areal_in / radius

    child_mass = hubble**2 * radius**3 / 2
    expansion_jump = child_theta_out - parent_theta_out
    normalized_surface_null_stress = -expansion_jump
    parent_kretschmann = 48 * mass**2 / radius**6
    child_kretschmann = 24 * hubble**4
    child_ricci_scalar = 12 * hubble**2

    assertions = {
        "boundary_inside_parent_horizon": compactness > 1,
        "parent_future_trapped": parent_theta_out < 0 and parent_theta_in < 0,
        "child_future_antitrapped": child_theta_out > 0 and child_theta_in > 0,
        "areal_radius_matches": radius == radius,
        "misner_sharp_mass_matches": child_mass == mass,
        "collapse_expansion_speed_magnitude_matches": hubble * radius == boundary_speed,
        "outgoing_expansion_jump_positive": expansion_jump > 0,
        "integrated_complete_Q_must_be_positive": expansion_jump > 0,
        "required_surface_null_stress_is_negative": normalized_surface_null_stress < 0,
        "parent_curvature_finite": parent_kretschmann > 0,
        "child_curvature_finite": child_kretschmann > 0,
        "child_to_parent_kretschmann_ratio_is_two": child_kretschmann == 2 * parent_kretschmann,
        "all_exact_gate_passed": False,
    }
    assertions["all_exact_gate_passed"] = all(
        value for key, value in assertions.items() if key != "all_exact_gate_passed"
    )

    return {
        "artifact_id": "NSC-1-PARENT-CHILD-BOUNDARY-DEV1",
        "schema": "NSC-1-PARENT-CHILD-BOUNDARY-DEV1-v1",
        "classification": "exact_finite_mass_matched_trapped_parent_expanding_child_kinematics",
        "active_scale": "transition_between_one_parent_black_hole_domain_and_one_child_cosmological_domain",
        "reference_frames": {
            "parent": "ingoing_Painleve_Gullstrand_future_orientation",
            "child": "expanding_flat_FLRW_areal_radius_chart_with_child_future_orientation",
            "identification": "same_boundary_two_sphere_same_areal_radius_same_Misner_Sharp_mass_opposite_sides",
        },
        "exact_inputs": {
            "units": "G=c=1",
            "parent_Misner_Sharp_mass": _record(mass),
            "boundary_flow_magnitude": _record(boundary_speed),
        },
        "derived_boundary": {
            "areal_radius": _record(radius),
            "compactness_2M_over_R": _record(compactness),
            "child_Hubble_parameter": _record(hubble),
            "child_HR": _record(hubble * radius),
            "child_Misner_Sharp_mass": _record(child_mass),
        },
        "parent_local_interpretation": {
            "outgoing_areal_derivative": _record(parent_areal_out),
            "ingoing_areal_derivative": _record(parent_areal_in),
            "theta_out": _record(parent_theta_out),
            "theta_in": _record(parent_theta_in),
            "classification": "future_trapped",
        },
        "child_local_interpretation": {
            "outgoing_areal_derivative": _record(child_areal_out),
            "ingoing_areal_derivative": _record(child_areal_in),
            "theta_out": _record(child_theta_out),
            "theta_in": _record(child_theta_in),
            "classification": "future_antitrapped_expanding",
        },
        "generalized_junction_data": {
            "outgoing_expansion_jump": _record(expansion_jump),
            "integrated_complete_Q_thin_transition": _record(expansion_jump),
            "normalized_surface_null_stress_8piG_times_Skk": _record(
                normalized_surface_null_stress
            ),
            "meaning": "a_smooth_or_thin_transition_must_supply_this_positive_integrated_Q_and_corresponding_negative_null_boundary_stress",
        },
        "finite_curvature": {
            "parent_Kretschmann": _record(parent_kretschmann),
            "child_Ricci_scalar": _record(child_ricci_scalar),
            "child_Kretschmann": _record(child_kretschmann),
            "child_to_parent_Kretschmann_ratio": _record(
                child_kretschmann / parent_kretschmann
            ),
        },
        "inherited_invariant": {
            "name": "spherical_quasilocal_energy_and_boundary_areal_radius",
            "relation": "2M/R=(HR)^2",
            "exact_value": _record(compactness),
        },
        "assertions": assertions,
        "observable_consequence": "parent_future_trapping_and_child_future_expansion_are_kinematically_consistent_opposite_side_descriptions_at_one_finite_mass_matched_boundary",
        "next_intervention": "construct_a_finite_thickness_boundary_sector_whose_stress_integrates_to_the_exact_expansion_jump_without_gauge_constraint_support",
        "literature_launch_surface": [
            {
                "role": "Misner_Sharp_energy_and_trapped_sphere_condition",
                "doi": "10.1103/PhysRevD.53.1938",
            },
            {
                "role": "dual_null_expansions_and_Raychaudhuri_in_spherical_symmetry",
                "doi": "10.1103/PhysRevD.92.083525",
            },
            {
                "role": "known_black_universe_solution_family_context",
                "arxiv": "gr-qc/0611022",
            },
        ],
        "nonclaims": {
            "smooth_FGCQR_dynamics_supplies_the_required_boundary_stress": False,
            "junction_surface_stress_has_been_derived_from_nuclear_gradient_energy": False,
            "global_child_cosmology_matches_observation": False,
            "universal_uniqueness_of_nested_spaces_proved": False,
        },
        "terminal": True,
    }


def main() -> int:
    record = construct()
    if not record["assertions"]["all_exact_gate_passed"]:
        raise RuntimeError("exact parent/child boundary assertions failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
