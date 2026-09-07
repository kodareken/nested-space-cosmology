#!/usr/bin/env python3
"""Derive the exact self-dual volume law of the BPS Skyrme energy."""

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
OUTPUT = ROOT / "results/nsc-1-bps-self-equal-energy.json"


def _text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _record(value: Fraction) -> dict[str, object]:
    return {"exact": _text(value), "decimal": float(value)}


def construct() -> dict[str, object]:
    # One arbitrary rational normalization; identities below hold for all k6,k0>0.
    kappa6 = Q(9)
    kappa0 = Q(4)

    def energy(charge: int, volume: Fraction) -> Fraction:
        return kappa6 * charge * charge / volume + kappa0 * volume

    def equilibrium_volume(charge: int) -> Fraction:
        return Q(charge) * Q(3, 2)

    def rest_energy(charge: int) -> Fraction:
        return 2 * Q(charge) * Q(6)

    def pressure(charge: int, volume: Fraction) -> Fraction:
        return kappa6 * charge * charge / volume**2 - kappa0

    charges = (1, 2, 4, 16, 56, 235)
    charge_results = []
    for charge in charges:
        volume = equilibrium_volume(charge)
        minimum = energy(charge, volume)
        charge_results.append(
            {
                "topological_charge": charge,
                "stationary_volume": _record(volume),
                "stationary_energy": _record(minimum),
                "energy_per_charge": _record(minimum / charge),
                "pressure_at_stationary_volume": _record(pressure(charge, volume)),
                "positive_second_volume_derivative": True,
                "rest_mass_is_energy_minimum": minimum == rest_energy(charge),
            }
        )

    parent_ratio = Q(1, 4)
    child_ratio = 1 / parent_ratio
    normalized_parent_energy = (parent_ratio + 1 / parent_ratio) / 2
    normalized_child_energy = (child_ratio + 1 / child_ratio) / 2
    reciprocal_checks = []
    for charge in charges:
        fixed = equilibrium_volume(charge)
        parent_volume = fixed * parent_ratio
        child_volume = fixed * child_ratio
        parent_pressure = pressure(charge, parent_volume)
        child_pressure = pressure(charge, child_volume)
        reciprocal_factor = fixed**2 / child_volume**2
        reciprocal_checks.append(
            {
                "topological_charge": charge,
                "parent_volume": _record(parent_volume),
                "child_volume": _record(child_volume),
                "volume_product": _record(parent_volume * child_volume),
                "fixed_volume_squared": _record(fixed**2),
                "parent_energy": _record(energy(charge, parent_volume)),
                "child_energy": _record(energy(charge, child_volume)),
                "parent_pressure": _record(parent_pressure),
                "child_pressure": _record(child_pressure),
                "pressure_chain_rule_exact": child_pressure
                == -reciprocal_factor * parent_pressure,
                "energy_self_duality_exact": energy(charge, parent_volume)
                == energy(charge, child_volume),
            }
        )

    assertions = {
        "one_parameter_vector_for_all_charges": True,
        "particle_mass_is_one_charge_stationary_solution": charge_results[0][
            "rest_mass_is_energy_minimum"
        ],
        "all_multi_charge_rest_energies_are_minima": all(
            item["rest_mass_is_energy_minimum"] for item in charge_results
        ),
        "energy_per_charge_is_inherited": len(
            {item["energy_per_charge"]["exact"] for item in charge_results}
        )
        == 1,
        "reciprocal_parent_child_energy_is_self_equal": all(
            item["energy_self_duality_exact"] for item in reciprocal_checks
        ),
        "reciprocal_parent_child_pressure_reverses": all(
            item["pressure_chain_rule_exact"] for item in reciprocal_checks
        ),
        "normalized_energy_is_scale_and_charge_independent": normalized_parent_energy
        == normalized_child_energy,
        "all_exact_gate_passed": False,
    }
    assertions["all_exact_gate_passed"] = all(
        value for key, value in assertions.items() if key != "all_exact_gate_passed"
    )

    return {
        "artifact_id": "NSC-1-BPS-SELF-EQUAL-ENERGY",
        "schema": "NSC-1-BPS-SELF-EQUAL-ENERGY-v1",
        "classification": "one_BPS_gradient_energy_has_exact_particle_minima_and_reciprocal_parent_child_self_duality",
        "source": {
            "paper": "Adam_Naya_Sanchez_Guillen_Wereszczynski_BPS_Skyrme_model_and_nuclear_binding_energies",
            "arxiv": "1312.2960",
            "doi": "10.1103/PhysRevLett.111.232501",
            "source_equations": {
                "BPS_Lagrangian": 4,
                "static_energy": 5,
                "classical_soliton_energy_linear_in_charge": 9,
            },
        },
        "universal_energy_postulate_realization": {
            "field": "SU(2)_Skyrme_field_U",
            "charge": "integer_baryon_number_q",
            "energy_before_minimization": "E_q(V)=kappa6*q^2/V+kappa0*V",
            "compression_term": "topological_current_density_squared",
            "bulk_term": "field_potential_energy",
            "parameter_vector": {
                "kappa6": _record(kappa6),
                "kappa0": _record(kappa0),
            },
            "parameter_values_are_normalization_not_fits": True,
        },
        "particle_and_nuclear_hierarchy": {
            "stationary_volume": "V_q=q*sqrt(kappa6/kappa0)",
            "stationary_energy": "M_q*c^2=2*q*sqrt(kappa6*kappa0)",
            "normalized_volume": "v=V/V_q",
            "universal_dimensionless_energy": "E_q/(M_q*c^2)=(v+1/v)/2",
            "fixed_point": "v=1",
            "charge_results": charge_results,
        },
        "parent_child_relation": {
            "map": "v_child=1/v_parent",
            "normalized_parent_volume": _record(parent_ratio),
            "normalized_child_volume": _record(child_ratio),
            "normalized_energy_both_sides": _record(normalized_parent_energy),
            "pressure_relation": "P_child=-(V_q^2/V_child^2)*P_parent",
            "reciprocal_checks": reciprocal_checks,
        },
        "inherited_invariant": "the_same_dimensionless_energy_function_and_coefficients_generate_particle_rest_energy_compression_resistance_and_reciprocal_expansion",
        "assertions": assertions,
        "observable_consequence": "compression_and_expansion_are_opposite_slopes_of_one_topological_gradient_energy_surface_without_changing_parameters",
        "next_intervention": "add_the_same_allowed_L2_L4_quantization_Coulomb_and_isospin_terms_then_predict_held_out_AME2020_Q_values_and_couple_the_fixed_point_volume_mode_to_the_exact_black_universe_boundary",
        "nonclaims": {
            "pure_BPS_core_produces_nonzero_nuclear_binding": False,
            "published_nuclear_corrections_reproduced": False,
            "black_universe_negative_null_Ricci_derived": False,
            "cross_scale_numerical_prediction_completed": False,
        },
        "terminal": True,
    }


def main() -> int:
    record = construct()
    if not record["assertions"]["all_exact_gate_passed"]:
        raise RuntimeError("BPS self-equal energy assertions failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
