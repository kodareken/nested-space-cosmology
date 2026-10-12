#!/usr/bin/env python3
"""Compare the unfitted flat-BPS collapse scale with compact-object observations."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


SOURCE = ROOT / "results/nsc-1-bps-black-universe-cross-scale-prediction.json"
OUTPUT = ROOT / "results/nsc-1-cross-scale-observation-dev1.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("cross-scale observation DEV1 output already exists")
    prediction = json.loads(SOURCE.read_text())
    predicted = prediction["headline_prediction"]["first_horizon_scale_solar_masses"]
    lvk_low_mass_bh = 3.6
    population_transition = 2.4
    population_uncertainty = 0.5
    inferred_ns_maximum = 2.15
    inferred_ns_upper = 2.15 + 0.14
    record = {
        "artifact_id": "NSC-1-CROSS-SCALE-OBSERVATION-DEV1",
        "schema": "NSC-1-CROSS-SCALE-OBSERVATION-DEV1-v1",
        "classification": "flat_space_BPS_radius_extrapolation_to_collapse_is_observationally_nonpassing",
        "source_prediction": str(SOURCE.relative_to(ROOT)),
        "frozen_prediction": {
            "horizon_onset_solar_masses": predicted,
            "mapping": "b=R1*B^(1/3)_without_self_gravitating_radius_contraction",
        },
        "observational_comparators": [
            {
                "name": "GW230529_primary",
                "mass_solar": lvk_low_mass_bh,
                "classification": "99_percent_probability_black_hole_under_reported_population_assumptions",
                "source": "https://ligo.org/science-summaries/GW230529/",
            },
            {
                "name": "compact_object_population_transition",
                "mass_solar": population_transition,
                "uncertainty_solar": population_uncertainty,
                "source": "https://dcc.ligo.org/LIGO-P2100403/public",
            },
            {
                "name": "NICER_chiral_EFT_maximum_neutron_star_mass",
                "central_mass_solar": inferred_ns_maximum,
                "upper_quoted_mass_solar": inferred_ns_upper,
                "source": "https://arxiv.org/abs/2407.06790",
            },
        ],
        "comparison": {
            "prediction_over_GW230529_primary": predicted / lvk_low_mass_bh,
            "prediction_over_population_transition": predicted / population_transition,
            "prediction_over_NS_maximum_upper": predicted / inferred_ns_upper,
            "observed_black_hole_below_predicted_onset": lvk_low_mass_bh < predicted,
            "nonpass_is_large_not_roundoff": predicted / lvk_low_mass_bh > 5.0,
        },
        "first_causal_obstruction": {
            "owner": "flat_space_compacton_radius_scaling_used_outside_its_domain",
            "excluded_neighbor": "published_one_charge_mass_and_radius_parameters_remain_the_nuclear_input",
            "smallest_fix": "solve_the_self_gravitating_BPS_or_near_BPS_mass_radius_curve_before_applying_the_horizon_condition",
            "new_collapse_fit_forbidden": True,
        },
        "interpretation": "the_universal_self_dual_energy_kernel_survives_but_the_state_map_from_flat_nuclear_soliton_to_stellar_compact_object_does_not",
        "nonclaims": {
            "Nested_Space_Cosmology_rejected": False,
            "exact_black_universe_solution_rejected": False,
            "BPS_nuclear_field_theory_rejected": False,
            "self_gravitating_BPS_prediction_already_known": False,
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
