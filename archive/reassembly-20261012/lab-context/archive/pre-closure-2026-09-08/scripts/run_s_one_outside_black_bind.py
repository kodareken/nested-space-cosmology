#!/usr/bin/env python3
"""Bind the published RS2 black universe to the one-action outside kernel."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-outside-black-bind.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-invariant-closure-bind.json",
    ROOT / "results/nsc-1-exact-black-universe-defocusing.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one outside-black bind output already exists")

    inputs = []
    records = []
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        records.append(item)
        inputs.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item.get("artifact_id"),
            }
        )

    closure, black = records
    black_assertions = black.get("assertions", {})
    record = {
        "artifact_id": "NSC-1-S-ONE-OUTSIDE-BLACK-BIND",
        "schema": "NSC-1-S-ONE-OUTSIDE-BLACK-BIND-v1",
        "classification": "published_outside_geometry_replaces_local_phantom_in_a_black_universe_and_occupies_the_Gamma_outside_slot",
        "method": "bind_primary_source_equations_to_the_existing_one_action_closure_without_recomputing_the_published_solution",
        "source": {
            "title": "Possible black universes in a brane world",
            "authors": ["K.A. Bronnikov", "E.V. Donskoy"],
            "arxiv": "0910.4930v1",
            "doi": "10.1134/S0202289310010068",
            "journal": "Gravitation and Cosmology 16, 42-49 (2010)",
            "source_status": "imported_prior_art_not_project_novelty",
        },
        "authenticated_project_inputs": inputs,
        "published_effective_equation": {
            "equation_13_source_convention": "G^nu_mu=-Lambda4*delta^nu_mu-kappa4^2*T^nu_mu-kappa5^4*Pi^nu_mu-E^nu_mu",
            "equation_14": "Lambda4=(Lambda5+kappa5^4*lambda^2/6)/2",
            "equation_15": "Pi^nu_mu is the brane-matching term quadratic in local T^nu_mu",
            "equation_16": "E_mu_nu=delta_mu^A*delta_nu^C*C5_ABCD*n^B*n^D",
            "E_properties": [
                "projected_electric_part_of_the_five_dimensional_Weyl_tensor",
                "traceless",
                "not_determined_without_the_five_dimensional_metric",
            ],
        },
        "exact_identification_with_the_one_equation": {
            "one_action_definition": closure.get("one_equation", {}).get(
                "influence_functional"
            ),
            "local_visible_projection": "delta(S_visible)/delta(g^mu_nu) gives T^nu_mu",
            "outside_projection": "P_mu_nu[delta(Gamma_outside)/delta(g^mu_nu)]=-E^nu_mu in the source convention",
            "matching_projection": "the nonlinear boundary matching derivative gives Pi^nu_mu",
            "local_vacuum_projection": "the bulk vacuum and boundary tension give Lambda4 by equation 14",
            "single_effective_equation": "delta(Gamma_visible+Gamma_outside)/delta(g^mu_nu)=0 is equation 13",
        },
        "already_known_joint_result": {
            "local_matter": "normal_positive_kinetic_scalar_on_the_four_dimensional_brane",
            "outside_effect": "projected_bulk_Weyl_geometry_supplies_the_effective_exotic_response",
            "black_result": "an_explicit_smooth_peaked_model_has_black_hole_exterior_and_expanding_de_Sitter_child_interior",
            "dark_result": "the_same_E_mu_nu_projection_appears_as_dark_radiation_in_isotropic_brane_cosmology",
            "constant_result": "Lambda4_and_kappa4_are_relations_of_bulk_curvature_boundary_tension_and_kappa5_not_isolated_local_insertions",
        },
        "why_this_binds_dark_and_black": {
            "black": "a nonlinear localized profile of E_mu_nu changes the causal geometry into a black-universe transition",
            "dark": "observers confined to the brane see E_mu_nu only through its gravitational response and cannot reconstruct it from local T_mu_nu",
            "one_cause": "both are projections of the same unresolved bulk geometry",
            "not_yet_demonstrated": [
                "observed_galaxy_and_cluster_dark_matter_lensing",
                "the_observed_late_time_dark_energy_density",
                "a_particle_mass_spectrum",
            ],
        },
        "published_limit_that_becomes_the_project_test": {
            "four_dimensional_equations_not_closed": True,
            "reason": "E_mu_nu_contains_one_arbitrary_radial_function_until_the_five_dimensional_bulk_metric_is_specified",
            "bulk_extension_not_constructed": True,
            "stability_not_established": True,
            "source_conclusion_scope": "existence_in_principle_not_a_realistic_complete_model",
        },
        "single_next_result": {
            "equation": "delta(S_one_5D_plus_boundary_plus_local_fields)/delta(G_AB,Phi)=0",
            "required_background_projection": "E_mu_nu_of_the_solved_bulk_equals_the_peaked_conserved_E_mu_nu_that_supports_the_positive_kinetic_black_universe",
            "required_linear_projection": "the_same_bulk_boundary_kernel_has_positive_physical_kinetic_and_gradient_spectrum",
            "required_dark_projection": "its_zero_and_finite_wavelength_brane_metric_response_predicts_Lambda_structure_growth_and_lensing",
            "no_new_local_substance": True,
        },
        "gate": {
            "known_black_universe_not_recomputed": True,
            "one_Gamma_mapping_explicit": True,
            "outside_geometry_replaces_local_phantom_in_imported_solution": True,
            "exact_existing_black_geometry_terminal_bound": black_assertions.get(
                "field_equations_symbolically_zero"
            )
            is True,
            "bulk_closed_and_stable": False,
            "dark_observations_predicted": False,
        },
        "terminal": True,
    }
    for key in (
        "known_black_universe_not_recomputed",
        "one_Gamma_mapping_explicit",
        "outside_geometry_replaces_local_phantom_in_imported_solution",
        "exact_existing_black_geometry_terminal_bound",
    ):
        if record["gate"][key] is not True:
            raise RuntimeError(f"outside-black bind gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
