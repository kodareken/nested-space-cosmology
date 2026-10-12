#!/usr/bin/env python3
"""Evaluate the published BPS-Skyrme nuclear mass formula on held-out Q-values."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from scripts.run_nuclear_gradient_q_dev1 import _load  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-published-bps-nuclear-q.json"
LAMBDA_MU_MEV = 48.9902
MU_OVER_LAMBDA_ONE_THIRD_PER_FM = 0.604327
ISOSPIN_BREAKING_MEV = -1.68593
HBAR_MEV_FM = 197.327
ELEMENTARY_CHARGE_C = 1.60218e-19
EPSILON0_INVERSE_MEV_FM = 8.8542e-21 / ELEMENTARY_CHARGE_C


def _soliton_energy(A: int) -> float:
    return 64.0 * math.sqrt(2.0) * math.pi / 15.0 * LAMBDA_MU_MEV * A


def _rotation_energy(A: int, Z: int, spin: float = 0.0) -> float:
    i3 = 0.5 * (2 * Z - A)
    prefactor = (
        105.0
        / (512.0 * math.sqrt(2.0) * math.pi)
        * HBAR_MEV_FM**2
        * MU_OVER_LAMBDA_ONE_THIRD_PER_FM**2
        * A ** (1.0 / 3.0)
        / LAMBDA_MU_MEV
    )
    return prefactor * (
        spin * (spin + 1.0) / A**2
        + 4.0 * abs(i3) * (abs(i3) + 1.0) / (3.0 * A**2 + 1.0)
    )


def _coulomb_energy(A: int, Z: int) -> float:
    i3 = 0.5 * (2 * Z - A)
    bracket = (
        128.0 * A**2 / (315.0 * math.pi**2)
        + 245.0 * A * i3 / 1536.0
        + 805.0 * i3**2 / 5148.0
        + 7.0 * i3**2 / (429.0 * (1.0 + 3.0 * A**2) ** 2)
    )
    return (
        1.0
        / (math.sqrt(2.0) * math.pi * EPSILON0_INVERSE_MEV_FM)
        * (MU_OVER_LAMBDA_ONE_THIRD_PER_FM / A ** (1.0 / 3.0))
        * bracket
    )


def _isospin_energy(A: int, Z: int) -> float:
    return ISOSPIN_BREAKING_MEV * 0.5 * (2 * Z - A)


def _energy(A: int, Z: int, spin: float = 0.0) -> dict[str, float]:
    if A <= 1 or A % 2:
        raise ValueError("published axial BPS formula here is restricted to even A>1")
    parts = {
        "classical_soliton": _soliton_energy(A),
        "spin_isospin_rotation": _rotation_energy(A, Z, spin),
        "Coulomb": _coulomb_energy(A, Z),
        "isospin_breaking": _isospin_energy(A, Z),
    }
    parts["total"] = sum(parts.values())
    return parts


REACTIONS = (
    (
        "He4_plus_He4_to_Be8",
        ((4, 2, 2),),
        ((8, 4, 1),),
        "fusion_transient_resonance",
    ),
    (
        "C12_plus_He4_to_O16",
        ((12, 6, 1), (4, 2, 1)),
        ((16, 8, 1),),
        "fusion_alpha_capture",
    ),
    (
        "O16_plus_He4_to_Ne20",
        ((16, 8, 1), (4, 2, 1)),
        ((20, 10, 1),),
        "fusion_alpha_capture",
    ),
    (
        "C12_plus_C12_to_Mg24",
        ((12, 6, 2),),
        ((24, 12, 1),),
        "fusion_carbon",
    ),
    (
        "U236_to_Ba144_plus_Kr92",
        ((236, 92, 1),),
        ((144, 56, 1), (92, 36, 1)),
        "fission_binary_channel",
    ),
    (
        "U238_to_Xe140_plus_Sr98",
        ((238, 92, 1),),
        ((140, 54, 1), (98, 38, 1)),
        "fission_binary_channel",
    ),
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("published BPS Q output already exists")
    masses = _load()
    evaluated = []
    for name, initial, final, kind in REACTIONS:
        for A, Z, _count in (*initial, *final):
            if (A, Z) not in masses:
                raise RuntimeError(f"AME2020 missing {(A, Z)}")
        initial_model = sum(
            count * _energy(A, Z)["total"] for A, Z, count in initial
        )
        final_model = sum(
            count * _energy(A, Z)["total"] for A, Z, count in final
        )
        predicted = initial_model - final_model
        measured = (
            sum(count * masses[(A, Z)]["mass_excess_keV"] for A, Z, count in initial)
            - sum(count * masses[(A, Z)]["mass_excess_keV"] for A, Z, count in final)
        ) / 1000.0
        evaluated.append(
            {
                "reaction": name,
                "kind": kind,
                "all_nuclei_even_even_spin_zero": True,
                "predicted_Q_MeV": predicted,
                "measured_Q_MeV": measured,
                "signed_error_MeV": predicted - measured,
                "absolute_error_MeV": abs(predicted - measured),
                "release_sign_correct": (predicted > 0.0) == (measured > 0.0),
                "initial_model_terms": [
                    {
                        "A": A,
                        "Z": Z,
                        "count": count,
                        "energy": _energy(A, Z),
                    }
                    for A, Z, count in initial
                ],
                "final_model_terms": [
                    {
                        "A": A,
                        "Z": Z,
                        "count": count,
                        "energy": _energy(A, Z),
                    }
                    for A, Z, count in final
                ],
            }
        )
    fusion = [item for item in evaluated if item["kind"].startswith("fusion")]
    fission = [item for item in evaluated if item["kind"].startswith("fission")]
    record = {
        "artifact_id": "NSC-1-PUBLISHED-BPS-NUCLEAR-Q",
        "schema": "NSC-1-PUBLISHED-BPS-NUCLEAR-Q-v1",
        "classification": "published_one_field_BPS_mass_formula_recovers_heavy_fission_but_not_light_fusion",
        "source": {
            "paper": "BPS_Skyrme_model_and_nuclear_binding_energies",
            "arxiv": "1312.2960",
            "doi": "10.1103/PhysRevLett.111.232501",
            "implemented_equations": [9, 18, 26, 27],
            "published_fit_inputs": ["proton_mass", "neutron_mass", "Ba138_mass"],
            "none_of_the_evaluated_reaction_nuclei_was_a_published_fit_input": True,
        },
        "published_parameters": {
            "lambda_times_mu_MeV": LAMBDA_MU_MEV,
            "mu_over_lambda_one_third_per_fm": MU_OVER_LAMBDA_ONE_THIRD_PER_FM,
            "a_I_MeV": ISOSPIN_BREAKING_MEV,
            "hbar_MeV_fm": HBAR_MEV_FM,
            "epsilon0_MeV_inverse_fm_inverse": EPSILON0_INVERSE_MEV_FM,
        },
        "AME2020": {
            "sha256": "e8599c6d7f724fac91934e59f1b9de8fb8f63e820f4b39456b790665ed2a3307",
            "mass_difference_Q_identity_used": True,
        },
        "reactions": evaluated,
        "summary": {
            "reaction_count": len(evaluated),
            "release_signs_correct": sum(item["release_sign_correct"] for item in evaluated),
            "fusion_signs_correct": sum(item["release_sign_correct"] for item in fusion),
            "fission_signs_correct": sum(item["release_sign_correct"] for item in fission),
            "fusion_mean_absolute_error_MeV": float(
                np.mean([item["absolute_error_MeV"] for item in fusion])
            ),
            "fission_mean_absolute_error_MeV": float(
                np.mean([item["absolute_error_MeV"] for item in fission])
            ),
        },
        "measured_owner": {
            "heavy_collective_BPS_sector": "captures_exothermic_fission_scale",
            "light_nuclear_failure": "missing_L2_L4_pion_single_particle_and_cluster_gradient_structure",
            "parameter_retuning_authorized": False,
        },
        "inherited_invariant": "all_nuclear_energies_are_evaluations_of_one_topological_field_energy_and_Q_is_the_difference_between_stationary_configurations",
        "next_intervention": "solve_the_near_BPS_L2_L4_light_cluster_sector_with_parameters_fixed_by_one_particle_solutions_then_repeat_the_same_held_out_reactions",
        "nonclaims": {
            "published_BPS_formula_is_accurate_for_light_nuclei": False,
            "same_parameters_generate_black_universe_defocusing": False,
            "cross_scale_invariant_closure_completed": False,
        },
        "terminal": True,
    }
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(
        json.dumps(
            {
                "classification": record["classification"],
                "summary": record["summary"],
                "reactions": [
                    {
                        key: item[key]
                        for key in (
                            "reaction",
                            "predicted_Q_MeV",
                            "measured_Q_MeV",
                            "signed_error_MeV",
                            "release_sign_correct",
                        )
                    }
                    for item in record["reactions"]
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
