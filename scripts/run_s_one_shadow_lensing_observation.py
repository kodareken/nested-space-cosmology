#!/usr/bin/env python3
"""Confront the published shadow-matter lensing ratio with galaxy data."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-shadow-lensing-observation.json"
SOURCE = ROOT / "results/nsc-1-s-one-shadow-matter-bind.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one shadow-lensing observation output already exists")

    source_raw = SOURCE.read_bytes()
    source = json.loads(source_raw)
    if source.get("terminal") is not True:
        raise RuntimeError("shadow-matter bind is not terminal")

    # Bolton, Rappaport, and Burles, arXiv:astro-ph/0607657.
    observed_gamma = 0.98
    observed_sigma = 0.07
    shadow_deflection_ratio = Fraction(3, 4)
    shadow_gamma = 2 * shadow_deflection_ratio - 1
    pure_shadow_z = (observed_gamma - float(shadow_gamma)) / observed_sigma

    # alpha/alpha_GR=(1+gamma)/2=1-f_shadow/4.
    maximum_likelihood_shadow_fraction = 2.0 * (1.0 - observed_gamma)
    shadow_fraction_sigma = 2.0 * observed_sigma
    two_sided_95_upper = max(
        0.0,
        maximum_likelihood_shadow_fraction + 1.96 * shadow_fraction_sigma,
    )

    record = {
        "artifact_id": "NSC-1-S-ONE-SHADOW-LENSING-OBSERVATION",
        "schema": "NSC-1-S-ONE-SHADOW-LENSING-OBSERVATION-v1",
        "classification": "constant_two_wall_shadow_response_cannot_dominate_the_tested_kiloparsec_gravitating_mass",
        "authenticated_prediction": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "artifact_id": source.get("artifact_id"),
            "prediction_opened_before_comparison": "shadow_deflection/ordinary_deflection=3/4",
        },
        "observational_source": {
            "title": "Constraint on the Post-Newtonian Parameter gamma on Galactic Size Scales",
            "authors": ["Adam S. Bolton", "Saul Rappaport", "Scott Burles"],
            "arxiv": "astro-ph/0607657",
            "doi": "10.1103/PhysRevD.74.061501",
            "sample": "15_SLACS_elliptical_lensing_galaxies",
            "method": "compare_Newtonian_dynamical_masses_with_gamma_dependent_strong_lensing_masses",
            "gamma": observed_gamma,
            "one_sigma": observed_sigma,
            "reported_confidence": "68_percent",
        },
        "model_translation": {
            "PPN_deflection_ratio": "alpha/alpha_GR=(1+gamma)/2",
            "shadow_mixture_ratio": "alpha/alpha_GR=1-f_shadow/4",
            "implied_relation": "gamma=1-f_shadow/2",
            "pure_shadow_gamma": str(shadow_gamma),
            "maximum_likelihood_shadow_fraction": maximum_likelihood_shadow_fraction,
            "shadow_fraction_one_sigma": shadow_fraction_sigma,
            "approximate_two_sided_95_percent_upper": two_sided_95_upper,
            "pure_shadow_discrepancy_sigma": pure_shadow_z,
        },
        "physical_result": {
            "pure_shadow_dominance_on_tested_scales": False,
            "reason": "the_frozen_two_wall_tensor_ratio_predicts_gamma_0_5_while_the_measurement_is_centered_near_1",
            "allowed_interpretation": "some_other_room_contribution_can_remain_but_the_constant_three_quarter_lensing_efficiency_cannot_supply_most_of_the_Newtonian_mass_in_this_sample",
            "not_a_theory_rejection": "the_tested_owner_is_the_constant_low_energy_two_wall_tensor_projection_not_recursive_outside_geometry_as_a_whole",
        },
        "constructive_kernel_revision": {
            "required_property": "the_tensor_link_B_Theta(p,ell) must yield gamma_eff_approximately_1_on_kiloparsec_lensing_and_dynamical scales",
            "preserved_property": "the_same_kernel_must_retain_nonlocal_outside_response_and_the_nonlinear_black_to_child_saddle",
            "forbidden_repair": "do_not_add_an_independent_dark_matter_stress_tensor_or_refit_the_observed_gamma",
            "next_test": "derive_the_scale_dependent_tensor_projection_from_the_recursive_kernel_and_compare_its_gamma_eff_across_solar_galaxy_and_transition_scales",
        },
        "scope": {
            "data_model_dependence": [
                "elliptical_density_profile_and_velocity_anisotropy_distributions",
                "weak_field_PPN_interpretation",
                "the_shadow_fraction_is_a_fraction_of_total_Newtonian_mass_within_the_tested_region",
            ],
            "does_not_directly_test": [
                "cosmic_dark_energy_fraction",
                "all_galaxy_or_cluster_environments",
                "the_parent_child_transition",
            ],
        },
        "gate": {
            "prediction_authenticated": source.get("terminal") is True,
            "pure_shadow_gamma_exactly_one_half": shadow_gamma == Fraction(1, 2),
            "pure_shadow_differs_by_more_than_six_sigma": pure_shadow_z > 6.0,
            "constant_two_wall_shadow_dominance_passes": False,
            "recursive_scale_dependent_kernel_tested": False,
        },
        "terminal": True,
    }
    for key in (
        "prediction_authenticated",
        "pure_shadow_gamma_exactly_one_half",
        "pure_shadow_differs_by_more_than_six_sigma",
    ):
        if record["gate"][key] is not True:
            raise RuntimeError(f"shadow-lensing observation gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
