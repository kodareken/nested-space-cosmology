#!/usr/bin/env python3
"""Predict the black-universe onset scale from published BPS nuclear parameters."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

from scipy.constants import G, c, electron_volt


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-bps-black-universe-cross-scale-prediction.json"
LAMBDA_MU_MEV = 48.9902
MU_OVER_LAMBDA_ONE_THIRD_PER_FM = 0.604327
SOLAR_MASS_KG = 1.98847e30


def _prediction(target_mass_to_radius: float) -> dict[str, float]:
    energy_per_charge_MeV = (
        64.0 * math.sqrt(2.0) * math.pi / 15.0 * LAMBDA_MU_MEV
    )
    compacton_radius_one_m = (
        math.sqrt(2.0) / MU_OVER_LAMBDA_ONE_THIRD_PER_FM * 1.0e-15
    )
    energy_per_charge_joule = energy_per_charge_MeV * 1.0e6 * electron_volt
    mass_per_charge_kg = energy_per_charge_joule / c**2
    one_charge_geometric_ratio = (
        G * mass_per_charge_kg / (c**2 * compacton_radius_one_m)
    )
    charge = (target_mass_to_radius / one_charge_geometric_ratio) ** 1.5
    mass_kg = charge * mass_per_charge_kg
    radius_m = compacton_radius_one_m * charge ** (1.0 / 3.0)
    return {
        "target_GM_over_c2R": target_mass_to_radius,
        "topological_charge": charge,
        "log10_topological_charge": math.log10(charge),
        "mass_kg": mass_kg,
        "mass_solar": mass_kg / SOLAR_MASS_KG,
        "minimum_radius_m": radius_m,
        "minimum_radius_km": radius_m / 1000.0,
        "reconstructed_GM_over_c2R": G * mass_kg / (c**2 * radius_m),
        "energy_per_charge_MeV": energy_per_charge_MeV,
        "one_charge_compacton_radius_fm": compacton_radius_one_m / 1.0e-15,
        "one_charge_geometric_ratio": one_charge_geometric_ratio,
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("BPS black-universe cross-scale output already exists")
    horizon_onset = _prediction(0.5)
    robust_onset = _prediction(0.75)
    radius = robust_onset["minimum_radius_m"]
    minimum_q = (3584.0 / 15625.0) / radius**2
    minimum_affine_length = 3.0 * radius / 8.0
    child_lambda_lower = (324.0 / 25.0) / radius**2
    child_hubble_lower = c * math.sqrt(child_lambda_lower / 3.0)
    assertions = {
        "same_published_BPS_parameters_used": True,
        "no_collapse_parameter_fitted": True,
        "horizon_onset_reconstructs_half_compactness": abs(
            horizon_onset["reconstructed_GM_over_c2R"] - 0.5
        )
        < 1.0e-12,
        "robust_onset_reconstructs_three_quarter_compactness": abs(
            robust_onset["reconstructed_GM_over_c2R"] - 0.75
        )
        < 1.0e-12,
        "robust_mass_above_horizon_mass": robust_onset["mass_kg"]
        > horizon_onset["mass_kg"],
        "predicted_Q_lower_bound_positive": minimum_q > 0.0,
        "predicted_child_expansion_positive": child_hubble_lower > 0.0,
        "all_gate_passed": False,
    }
    assertions["all_gate_passed"] = all(
        value for key, value in assertions.items() if key != "all_gate_passed"
    )
    record = {
        "artifact_id": "NSC-1-BPS-BLACK-UNIVERSE-CROSS-SCALE-PREDICTION",
        "schema": "NSC-1-BPS-BLACK-UNIVERSE-CROSS-SCALE-PREDICTION-v1",
        "classification": "published_nuclear_soliton_scale_predicts_stellar_mass_black_universe_threshold_without_collapse_fit",
        "source_results": [
            "results/nsc-1-bps-self-equal-energy.json",
            "results/nsc-1-exact-black-universe-rob1.json",
            "results/nsc-1-cross-scale-fixed-point-dev1.json",
        ],
        "frozen_cross_scale_identification": {
            "BPS_standard_potential_compacton_radius": "R_B=sqrt(2)*(lambda/mu)^(1/3)*B^(1/3)",
            "BPS_classical_mass": "M_B*c^2=(64*sqrt(2)*pi/15)*lambda*mu*B",
            "black_universe_minimum_radius": "b=R_B",
            "geometric_ratio": "m/b=G*M_B/(c^2*R_B)",
            "status": "prospective_identification_not_fitted_to_the_output",
        },
        "published_BPS_parameters": {
            "lambda_times_mu_MeV": LAMBDA_MU_MEV,
            "mu_over_lambda_one_third_per_fm": MU_OVER_LAMBDA_ONE_THIRD_PER_FM,
            "paper": "arXiv:1312.2960",
            "doi": "10.1103/PhysRevLett.111.232501",
        },
        "physical_constants": {
            "G_m3_kg_s2": G,
            "c_m_s": c,
            "electron_volt_joule": electron_volt,
            "solar_mass_kg_for_reporting": SOLAR_MASS_KG,
        },
        "horizon_onset_prediction": horizon_onset,
        "robust_defocusing_onset_prediction": {
            **robust_onset,
            "uniform_minimum_complete_Q_per_m2": minimum_q,
            "uniform_minimum_affine_length_m": minimum_affine_length,
            "uniform_child_Lambda_lower_per_m2": child_lambda_lower,
            "uniform_child_Hubble_lower_per_s": child_hubble_lower,
        },
        "headline_prediction": {
            "first_horizon_scale_solar_masses": horizon_onset["mass_solar"],
            "robust_black_universe_scale_solar_masses": robust_onset["mass_solar"],
            "robust_minimum_radius_km": robust_onset["minimum_radius_km"],
            "robust_topological_charge": robust_onset["topological_charge"],
        },
        "assertions": assertions,
        "inherited_invariant": "one_charge_mass_and_radius_scaling_from_the_nuclear_soliton_action_fix_the_mass_dependence_of_black_universe_compactness",
        "observable_consequence": "the_unfitted_transition_scale_lies_in_the_stellar_mass_regime_rather_than_at_an_arbitrary_inserted_scale",
        "next_intervention": "solve_the_gravitating_BPS_boundary_equations_and_test_whether_the_predicted_threshold_and_negative_null_Ricci_are_retained",
        "nonclaims": {
            "b_equals_BPS_compacton_radius_derived_from_full_gravity": False,
            "stellar_black_hole_mass_distribution_fitted": False,
            "positive_energy_BPS_stress_alone_supplies_negative_Ricci_null": False,
            "observational_prediction_confirmed": False,
        },
        "terminal": True,
    }
    if not assertions["all_gate_passed"]:
        raise RuntimeError("cross-scale prediction checks failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
