#!/usr/bin/env python3
"""Test one bulk/boundary gradient-energy law on held-out nuclear Q-values."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


DATA = ROOT / "collected-data/ame2020/mass_1.mas20.txt"
OUTPUT = ROOT / "results/nsc-1-nuclear-gradient-q-dev1.json"
HOLDOUT = {
    (1, 0),
    (1, 1),
    (2, 1),
    (3, 1),
    (3, 2),
    (4, 2),
    (7, 3),
    (8, 4),
    (92, 36),
    (141, 56),
    (235, 92),
}


def _number(field: str) -> float | None:
    text = field.strip()
    if not text or "*" in text:
        return None
    try:
        return float(text.replace("#", ""))
    except ValueError:
        return None


def _load():
    records = {}
    for line in DATA.read_text().splitlines():
        try:
            neutron = int(line[4:9])
            proton = int(line[9:14])
            mass_number = int(line[14:19])
        except (ValueError, IndexError):
            continue
        if neutron + proton != mass_number:
            continue
        mass_excess = _number(line[28:42])
        binding_per_nucleon = _number(line[54:67])
        if mass_excess is None or binding_per_nucleon is None:
            continue
        records[(mass_number, proton)] = {
            "A": mass_number,
            "Z": proton,
            "N": neutron,
            "element": line[20:23].strip(),
            "mass_excess_keV": mass_excess,
            "binding_energy_MeV": binding_per_nucleon * mass_number / 1000.0,
            "estimated": "#" in line[28:67],
        }
    return records


def _features(A: int, Z: int) -> np.ndarray:
    if A <= 1:
        return np.zeros(5, dtype=np.float64)
    if A % 2 == 1:
        pairing = 0.0
    elif Z % 2 == 0 and (A - Z) % 2 == 0:
        pairing = 1.0 / np.sqrt(A)
    else:
        pairing = -1.0 / np.sqrt(A)
    return np.asarray(
        (
            A,
            -(A ** (2.0 / 3.0)),
            -Z * (Z - 1) / (A ** (1.0 / 3.0)),
            -((A - 2 * Z) ** 2) / A,
            pairing,
        ),
        dtype=np.float64,
    )


def _binding(coefficients: np.ndarray, A: int, Z: int) -> float:
    return float(_features(A, Z) @ coefficients)


def _species(A: int, Z: int, count: int = 1):
    return (A, Z, count)


REACTIONS = (
    (
        "D_plus_D_to_He3_plus_n",
        (_species(2, 1, 2),),
        (_species(3, 2), _species(1, 0)),
    ),
    (
        "D_plus_D_to_T_plus_H1",
        (_species(2, 1, 2),),
        (_species(3, 1), _species(1, 1)),
    ),
    (
        "T_plus_D_to_He4_plus_n",
        (_species(3, 1), _species(2, 1)),
        (_species(4, 2), _species(1, 0)),
    ),
    (
        "Li7_plus_H1_to_two_He4",
        (_species(7, 3), _species(1, 1)),
        (_species(4, 2, 2),),
    ),
    (
        "Be8_to_two_He4",
        (_species(8, 4),),
        (_species(4, 2, 2),),
    ),
    (
        "U235_plus_n_to_Ba141_plus_Kr92_plus_three_n",
        (_species(235, 92), _species(1, 0)),
        (_species(141, 56), _species(92, 36), _species(1, 0, 3)),
    ),
)


def _sum_mass_excess(records, side) -> float:
    return sum(count * records[(A, Z)]["mass_excess_keV"] for A, Z, count in side)


def _sum_predicted_binding(coefficients, side) -> float:
    return sum(count * _binding(coefficients, A, Z) for A, Z, count in side)


def run() -> dict[str, object]:
    records = _load()
    missing = sorted(HOLDOUT - set(records))
    if missing:
        raise RuntimeError(f"AME2020 is missing held-out nuclides: {missing}")
    training = [
        record
        for key, record in records.items()
        if 16 <= record["A"] <= 240
        and not record["estimated"]
        and key not in HOLDOUT
    ]
    design = np.asarray([_features(item["A"], item["Z"]) for item in training])
    target = np.asarray([item["binding_energy_MeV"] for item in training])
    coefficients, _residuals, rank, singular = np.linalg.lstsq(design, target, rcond=None)
    if rank != 5:
        raise RuntimeError("gradient droplet fit is rank deficient")
    fitted = design @ coefficients
    reactions = []
    for name, initial, final in REACTIONS:
        measured = (_sum_mass_excess(records, initial) - _sum_mass_excess(records, final)) / 1000.0
        predicted = _sum_predicted_binding(coefficients, final) - _sum_predicted_binding(
            coefficients, initial
        )
        reactions.append(
            {
                "reaction": name,
                "measured_Q_MeV_from_Delta_m_c2": measured,
                "predicted_Q_MeV_from_gradient_energy_difference": predicted,
                "signed_error_MeV": predicted - measured,
                "absolute_error_MeV": abs(predicted - measured),
                "measured_release_positive": measured > 0.0,
                "predicted_release_positive": predicted > 0.0,
            }
        )
    errors = np.asarray([item["absolute_error_MeV"] for item in reactions])
    signs = [item["measured_release_positive"] == item["predicted_release_positive"] for item in reactions]
    record = {
        "artifact_id": "NSC-1-NUCLEAR-GRADIENT-Q-DEV1",
        "schema": "NSC-1-NUCLEAR-GRADIENT-Q-DEV1-v1",
        "classification": "single_bulk_boundary_gradient_functional_tested_on_held_out_fusion_and_fission_Q_values",
        "active_scale": "nuclear_stable_gradient_configurations",
        "data": {
            "source": "AME2020_mass_1.mas20",
            "sha256": "e8599c6d7f724fac91934e59f1b9de8fb8f63e820f4b39456b790665ed2a3307",
            "doi_part_I": "10.1088/1674-1137/abddb0",
            "doi_part_II": "10.1088/1674-1137/abddaf",
            "training_nuclide_count": len(training),
            "all_reaction_nuclides_held_out": True,
        },
        "coarse_grained_gradient_energy": {
            "binding_formula": "B=a_v*A-a_s*A^(2/3)-a_c*Z*(Z-1)/A^(1/3)-a_a*(A-2Z)^2/A+delta_pair",
            "interpretation": {
                "volume": "stable_interior_gradient_energy",
                "surface": "finite_boundary_gradient_cost",
                "coulomb": "outward_electromagnetic_push",
                "asymmetry": "neutron_proton_state_opposition_cost",
                "pairing": "local_correlation_stability",
            },
            "coefficients_MeV": {
                "a_volume": float(coefficients[0]),
                "a_surface": float(coefficients[1]),
                "a_coulomb": float(coefficients[2]),
                "a_asymmetry": float(coefficients[3]),
                "a_pairing": float(coefficients[4]),
            },
            "training_binding_RMSE_MeV": float(np.sqrt(np.mean((fitted - target) ** 2))),
            "design_condition": float(singular[0] / singular[-1]),
        },
        "energy_identity": "Q=(sum_M_initial-sum_M_final)*c^2=sum_B_final-sum_B_initial_for_A_Z_conserving_reactions",
        "held_out_reactions": reactions,
        "summary": {
            "reaction_count": len(reactions),
            "release_signs_correct": int(sum(signs)),
            "mean_absolute_Q_error_MeV": float(np.mean(errors)),
            "maximum_absolute_Q_error_MeV": float(np.max(errors)),
            "fusion_mean_absolute_Q_error_MeV": float(np.mean(errors[:-1])),
            "fission_absolute_Q_error_MeV": float(errors[-1]),
        },
        "inherited_invariant": "energy_release_is_the_difference_between_two_stable_bounded_configurations_not_a_separate_energy_source",
        "next_intervention": "replace_the_coarse_droplet_terms_with_one_gravitating_near_BPS_Skyrme_field_functional_and_require_both_nuclear_Q_predictions_and_boundary_negative_Ricci_null",
        "nonclaims": {
            "semi_empirical_formula_is_a_microscopic_particle_wave_theory": False,
            "nuclear_fit_derives_the_parent_child_boundary_stress": False,
            "held_out_accuracy_is_known_before_execution": False,
            "Nested_Space_Cosmology_proved": False,
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
