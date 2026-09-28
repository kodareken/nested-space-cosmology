#!/usr/bin/env python3
"""Bind the nuclear-fitted gravitating BPS maximum mass to observations."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-gravitating-bps-observation-link.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("gravitating BPS observation link output already exists")
    flat_prediction = 24.444890346689895
    exact_maximum = 3.34
    exact_compactness = 0.58
    mean_field_maximum = 3.79
    mean_field_compactness = 0.68
    observed_bh = 3.6
    population_transition = 2.4
    population_sigma = 0.5
    inferred_ns_maximum = 2.15
    record = {
        "artifact_id": "NSC-1-GRAVITATING-BPS-OBSERVATION-LINK",
        "schema": "NSC-1-GRAVITATING-BPS-OBSERVATION-LINK-v1",
        "classification": "nuclear_fitted_self_gravitating_BPS_scale_lands_in_observed_compact_object_transition_region",
        "source": {
            "paper": "Neutron_stars_in_the_BPS_Skyrme_model_mean_field_limit_vs_full_field_theory",
            "arxiv": "1503.03095v2",
            "doi": "10.1103/PhysRevC.92.025802",
            "field_equations": ["IV.1", "IV.5", "IV.7", "IV.8", "IV.9", "IV.10"],
            "nuclear_fit": "IV.15",
            "maximum_mass_table": "Table_I",
            "compactness_table": "Table_II",
            "published_result_not_yet_independently_recomputed": True,
        },
        "one_action_nuclear_inputs": {
            "potential": "U_pi=2h=1-cos(xi)",
            "binding_energy_per_nucleon_MeV": 16.3,
            "nuclear_saturation_density_per_fm3": 0.153,
            "nucleon_energy_MeV": 939.6,
            "lambda_squared_MeV_fm3": 26.88,
            "mu_squared_MeV_per_fm3": 88.26,
            "gravity_coupling_fitted": False,
        },
        "self_gravitating_predictions": {
            "full_field_maximum_mass_solar": exact_maximum,
            "full_field_maximum_compactness_2GM_over_Rc2": exact_compactness,
            "mean_field_maximum_mass_solar": mean_field_maximum,
            "mean_field_maximum_compactness_2GM_over_Rc2": mean_field_compactness,
            "stable_branch_ends_before_horizon": exact_compactness < 1.0,
        },
        "repair_of_flat_extrapolation": {
            "invalid_flat_space_horizon_onset_solar": flat_prediction,
            "self_gravitating_full_field_maximum_solar": exact_maximum,
            "mass_reduction_factor": flat_prediction / exact_maximum,
            "causal_owner": "gravity_changes_the_soliton_profile_radius_and_on_shell_equation_of_state",
        },
        "observational_comparison": {
            "GW230529_likely_BH_mass_solar": observed_bh,
            "difference_from_full_field_maximum_solar": observed_bh - exact_maximum,
            "population_transition_mass_solar": population_transition,
            "population_transition_uncertainty_solar": population_sigma,
            "full_field_offset_in_population_sigma": (
                exact_maximum - population_transition
            )
            / population_sigma,
            "current_inferred_NS_maximum_central_solar": inferred_ns_maximum,
            "full_field_prediction_in_3_to_5_solar_mass_region": 3.0
            <= exact_maximum
            <= 5.0,
        },
        "one_model_status": {
            "same_BPS_L0_L6_action_at_nuclear_and_stellar_scales": True,
            "same_potential": True,
            "nuclear_parameters_fixed_before_gravity": True,
            "collapse_mass_fitted": False,
            "parent_child_defocusing_source_derived": False,
        },
        "inherited_invariant": "one_nuclear_fitted_topological_energy_functional_predicts_its_own_self_gravitating_maximum_mass_when_G_is_applied",
        "next_intervention": "continue_the_full_field_branch_past_loss_of_static_support_into_the_shared_boundary_sector_and_test_for_the_exact_black_universe_Q_interval",
        "nonclaims": {
            "3_34_solar_mass_maximum_is_a_precision_observational_fit": False,
            "static_star_calculation_is_a_collapse_trajectory": False,
            "BPS_positive_energy_stress_already_produces_defocusing": False,
            "complete_parent_child_energy_ledger_closed": False,
        },
        "terminal": True,
    }
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
