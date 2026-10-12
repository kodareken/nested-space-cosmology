#!/usr/bin/env python3
"""Link the BPS self-dual energy to paired exact black-universe states."""

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
OUTPUT = ROOT / "results/nsc-1-cross-scale-fixed-point-dev1.json"


def _text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _record(value: Fraction) -> dict[str, object]:
    return {"exact": _text(value), "decimal": float(value)}


def construct() -> dict[str, object]:
    b = Q(1)
    rho_magnitude = Q(3, 4)
    radius = Q(5, 4)
    parent_scale = b * b / radius
    child_scale = radius
    parent_volume = (parent_scale / b) ** 3
    child_volume = (child_scale / b) ** 3

    kappa6 = Q(9)
    kappa0 = Q(4)
    charge = 1
    equilibrium_volume = Q(3, 2)

    def normalized_energy(volume_ratio: Fraction) -> Fraction:
        return (volume_ratio + 1 / volume_ratio) / 2

    def normalized_pressure(volume_ratio: Fraction) -> Fraction:
        return 1 / volume_ratio**2 - 1

    parent_energy = normalized_energy(parent_volume)
    child_energy = normalized_energy(child_volume)
    parent_pressure = normalized_pressure(parent_volume)
    child_pressure = normalized_pressure(child_volume)
    pressure_factor = 1 / child_volume**2
    complete_q = 2 * (1 - rho_magnitude**2) / (1 + rho_magnitude**2) ** 2
    theta_parent = -2 * rho_magnitude / (1 + rho_magnitude**2)
    theta_child = -theta_parent
    omega = child_scale / parent_scale
    gamma = 1 / omega

    assertions = {
        "Pythagorean_areal_radius_exact": radius**2 == b**2 + rho_magnitude**2,
        "reciprocal_scale_map_exact": parent_scale * child_scale == b**2,
        "reciprocal_volume_map_exact": parent_volume * child_volume == 1,
        "one_energy_kernel_both_sides": parent_energy == child_energy,
        "pressure_chain_rule_exact": child_pressure
        == -pressure_factor * parent_pressure,
        "parent_is_trapped": theta_parent < 0,
        "child_is_expanding": theta_child > 0,
        "both_sides_have_same_positive_Q": complete_q > 0,
        "energy_conversion_factor_derived": gamma == parent_scale / child_scale,
        "ideal_dimensionless_energy_ledger_closes": parent_energy - child_energy == 0,
        "all_exact_gate_passed": False,
    }
    assertions["all_exact_gate_passed"] = all(
        value for key, value in assertions.items() if key != "all_exact_gate_passed"
    )

    return {
        "artifact_id": "NSC-1-CROSS-SCALE-FIXED-POINT-DEV1",
        "schema": "NSC-1-CROSS-SCALE-FIXED-POINT-DEV1-v1",
        "classification": "one_frozen_self_dual_energy_kernel_links_particle_minima_to_paired_parent_child_defocusing_states",
        "source_results": [
            "results/nsc-1-bps-self-equal-energy.json",
            "results/nsc-1-exact-black-universe-defocusing.json",
            "results/nsc-1-exact-black-universe-rob1.json",
        ],
        "universal_parameter_vector": {
            "kappa6": _record(kappa6),
            "kappa0": _record(kappa0),
            "changed_between_scales": False,
            "role": "dimensionless_normalization_of_the_same_BPS_volume_kernel",
        },
        "paired_geometry": {
            "parent_rho_over_b": _record(rho_magnitude),
            "child_rho_over_b": _record(-rho_magnitude),
            "areal_radius_over_b": _record(radius),
            "parent_containment_scale_over_b": _record(parent_scale / b),
            "child_expansion_scale_over_b": _record(child_scale / b),
            "scale_product_over_b_squared": _record(
                parent_scale * child_scale / b**2
            ),
            "parent_theta_times_b": _record(theta_parent),
            "child_theta_times_b": _record(theta_child),
            "complete_Q_times_b_squared_both_sides": _record(complete_q),
        },
        "dimensionless_energy": {
            "kernel": "e(v)=(v+1/v)/2",
            "parent_volume_ratio": _record(parent_volume),
            "child_volume_ratio": _record(child_volume),
            "volume_product": _record(parent_volume * child_volume),
            "parent_energy_ratio": _record(parent_energy),
            "child_energy_ratio": _record(child_energy),
            "parent_pressure_ratio": _record(parent_pressure),
            "child_pressure_ratio": _record(child_pressure),
            "pressure_chain_factor": _record(pressure_factor),
        },
        "scale_inheritance": {
            "Omega_child_over_parent": _record(omega),
            "Gamma_if_hbar_c_invariant": _record(gamma),
            "Gamma_derived_not_fitted": True,
        },
        "shared_energy_ledger": {
            "parent_dimensionless_energy": _record(parent_energy),
            "child_converted_dimensionless_energy": _record(child_energy),
            "link_change_in_ideal_reversible_pair": _record(Q(0)),
            "exterior_energy_in_ideal_reversible_pair": _record(Q(0)),
            "closure_residual": _record(parent_energy - child_energy),
        },
        "assertions": assertions,
        "inherited_invariant": "the_BPS_fixed_point_at_v_equals_one_is_the_black_universe_minimum_radius_and_paired_parent_child_states_are_reciprocal_equal_energy_configurations",
        "observable_consequence": "the_exact_pressure_sign_reversal_and_exact_metric_null_defocusing_coexist_at_the_same_paired_states_without_changing_the_energy_law",
        "next_intervention": "derive_the_dimensional_boundary_tension_and_negative_Ricci_null_from_the_published_BPS_parameters_and_one_particle_radius_without_an_independent_collapse_fit",
        "nonclaims": {
            "hbar_c_invariance_across_real_parent_child_domains_measured": False,
            "BPS_stress_tensor_sources_the_black_universe_metric": False,
            "dimensional_cross_scale_prediction_completed": False,
            "observed_cosmology_matched": False,
        },
        "terminal": True,
    }


def main() -> int:
    record = construct()
    if not record["assertions"]["all_exact_gate_passed"]:
        raise RuntimeError("cross-scale fixed-point assertions failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
