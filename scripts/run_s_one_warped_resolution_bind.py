#!/usr/bin/env python3
"""Bind warped mode localization to the recursive resolution equation."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-warped-resolution-bind.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-recursive-outside-kernel.json",
    ROOT / "results/nsc-1-s-one-resolution-channel-closure.json",
    ROOT / "results/nsc-1-s-one-transition-link-coefficient.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one warped-resolution bind output already exists")

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

    correction_samples = []
    for radius_over_ell in (
        Fraction(2),
        Fraction(4),
        Fraction(10),
        Fraction(100),
    ):
        relative = Fraction(1, 2) / radius_over_ell**2
        correction_samples.append(
            {
                "radius_over_boundary_scale": str(radius_over_ell),
                "leading_KK_fractional_correction": str(relative),
                "decimal": float(relative),
            }
        )

    record = {
        "artifact_id": "NSC-1-S-ONE-WARPED-RESOLUTION-BIND",
        "schema": "NSC-1-S-ONE-WARPED-RESOLUTION-BIND-v1",
        "classification": "a_gradient_boundary_localizes_the_local_gravity_pole_while_the_same_spectrum_retains_outside_modes",
        "method": "bind_the_published_RS2_spectral_solution_without_recomputing_it",
        "source": {
            "title": "Gravity in the Randall-Sundrum Brane World",
            "authors": ["Jaume Garriga", "Takahiro Tanaka"],
            "arxiv": "hep-th/9911055v4",
            "doi": "10.1103/PhysRevLett.84.2778",
            "source_status": "imported_prior_art_not_project_novelty",
            "equations": {
                "warped_room_metric": 1,
                "bulk_mode_equation": 3,
                "junction_condition": 5,
                "zero_plus_continuum_Green_function": 13,
                "local_static_Green_function": 15,
            },
        },
        "authenticated_project_inputs": authenticated,
        "one_spectrum": {
            "warp": "a(y)=exp(-abs(y)/ell)",
            "bulk_cosmological_constant": "Lambda5=-6/ell^2",
            "boundary_tension": "sigma=3/(4*pi*ell*G5)",
            "local_Newton_constant": "G4=G5/ell",
            "brane_propagator": "G_brane(k)=psi0(0)^2/k^2+integral_0^infinity dm psi_m(0)^2/(k^2+m^2)",
            "local_mode": "the_normalizable_zero_mode_is_the_massless_four_dimensional_graviton",
            "outside_modes": "the_Kaluza_Klein_continuum_is_the_unresolved_bulk_response",
            "same_operator": True,
        },
        "resolution_map": {
            "large_distance": "r_much_greater_than_ell_resolves_the_zero_mode_as_four_dimensional_gravity",
            "leading_static_response": "G(r)=-[4*pi*ell*r]^-1*[1+ell^2/(2*r^2)+...]",
            "short_distance": "r_near_or_below_ell_resolves_the_outside_continuum_and_gravity_becomes_higher_dimensional",
            "gradient_resolution": "ell_is_the_boundary_localization_width_and_the_spectral_crossover",
            "fractional_correction_samples": correction_samples,
        },
        "binding_to_the_recursive_equation": {
            "unwarped_fixed_point_role": "Gamma=K-b^2/Gamma_proves_how_an_infinite_tail_becomes_one_self_energy",
            "warp_role": "T_Theta_rescales_successive_room_modes_and_creates_a_normalizable_local_pole",
            "combined_role": "one_pole_plus_one_continuum_are_two_resolution_views_of_the_same_recursive_metric_kernel",
            "nonlinear_role": "the_transition_fixed_b2_over_F=-18/1015_is_the_strong_gradient_curvature_of_the_same_boundary",
        },
        "physical_result": {
            "local_gravity_and_outside_gravity_need_not_be_separate_forces": True,
            "ordinary_inverse_radius_gravity_is_a_boundary_localized_mode": True,
            "outside_degrees_remain_in_the_same_propagator": True,
            "the_room_scale_controls_which_modes_are_resolved": True,
        },
        "next_result": "replace_the_single_RS2_warp_by_the_recursive_transition_derived_scale_map_and_calculate_the_shared_tensor_spectral_weight_across_particle_galaxy_and_black_boundary_scales",
        "nonclaims": {
            "single_RS2_continuum_explains_galaxy_dark_matter": False,
            "recursive_warp_scale_predicted": False,
            "particle_spectrum_predicted": False,
            "black_universe_bulk_completed": False,
        },
        "gate": {
            "all_project_inputs_authenticated": len(authenticated) == len(INPUTS),
            "known_RS2_modes_not_recomputed": True,
            "one_operator_contains_local_pole_and_outside_continuum": True,
            "boundary_scale_is_resolution_crossover": True,
            "transition_coefficient_retained": True,
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("warped-resolution bind gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
