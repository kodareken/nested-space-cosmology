#!/usr/bin/env python3
"""Run the first one-action particle/wave reduction gate."""

from __future__ import annotations

from fractions import Fraction
import json
from math import pi, sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.unified_action import (  # noqa: E402
    OPERATOR_BASIS,
    THETA_ORDER,
    OneAction,
    UniversalTheta,
)


OUTPUT = ROOT / "results/nsc-1-s-one-particle-wave-identity.json"


def _compacton_quadrature(order: int) -> tuple[float, float]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    radius = sqrt(2.0)
    radial = radius * (nodes + 1.0) / 2.0
    scaled = radial / radius
    baryon_density = 4.0 / (pi**2 * radius**3) * np.sqrt(1.0 - scaled**2)
    energy_density = 4.0 * (1.0 - scaled**2)
    jacobian = radius / 2.0
    charge = float(np.sum(weights * 4.0 * pi * radial**2 * baryon_density) * jacobian)
    energy = float(np.sum(weights * 4.0 * pi * radial**2 * energy_density) * jacobian)
    return charge, energy


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one particle/wave identity output already exists")
    theta = UniversalTheta(
        length_anchor=1.0,
        coefficients={
            "kappa0": Fraction(1),
            "kappa2": Fraction(0),
            "kappa4": Fraction(0),
            "kappa6": Fraction(1),
            "nonminimal_chi": Fraction(0),
            "degenerate_boundary": Fraction(0),
            "gauge_coupling": Fraction(0),
            "spin_coupling": Fraction(0),
            "mixing_coupling": Fraction(0),
            "boundary_coupling": Fraction(0),
        },
    )
    action = OneAction(theta)
    sector = action.stationary_charge_sector(1)
    wave = action.collective_wave_identity(1)
    quadrature = []
    analytic_energy = 64.0 * sqrt(2.0) * pi / 15.0
    for order in (64, 128, 256):
        charge, energy = _compacton_quadrature(order)
        quadrature.append(
            {
                "order": order,
                "charge": charge,
                "charge_error": abs(charge - 1.0),
                "energy": energy,
                "energy_error": abs(energy - analytic_energy),
                "normalized_energy": 2.0 * energy / analytic_energy,
                "normalized_energy_error": abs(2.0 * energy / analytic_energy - 2.0),
            }
        )
    result = {
        "artifact_id": "NSC-1-S-ONE-PARTICLE-WAVE-IDENTITY",
        "schema": "NSC-1-S-ONE-PARTICLE-WAVE-IDENTITY-v1",
        "classification": "one_charge_localized_configuration_and_collective_wave_amplitude_share_one_field_state",
        "status": "imported_collective_coordinate_reduction_gate_not_project_novelty",
        "one_action": {
            "formula": action.formula,
            "quantum_generator": action.quantum_generator,
            "field_operator_basis": list(OPERATOR_BASIS),
            "theta_order": list(THETA_ORDER),
            "hard_coded_particle_mass_present": False,
            "separate_point_particle_present": False,
            "separate_pilot_wave_present": False,
        },
        "localized_configuration": {
            "field": "BPS_Skyrme_U",
            "topological_charge": 1,
            "profile": "xi(r)=2*acos(r/R) for r<=R; xi=0 outside",
            "dimensionless_radius": sqrt(2.0),
            "normalization": "geometric_compacton_volume_and_analytic_energy_are_scaled_to_stationary_volume_coordinate_1_and_rest_energy_2",
            "stationary_volume_coordinate": sector.stationary_volume,
            "stationary_energy_in_action_units": sector.stationary_energy,
            "positive_volume_hessian": sector.volume_hessian,
            "quadrature_reduction": quadrature,
            "analytic_BPS_energy": analytic_energy,
        },
        "same_configuration_wave": {
            "definition": wave.wave_amplitude,
            "translation_zero_mode": wave.translation_zero_mode,
            "localized_rest_energy": wave.localized_rest_energy,
            "wave_dispersion_rest_energy": wave.wave_dispersion_rest_energy,
            "same_rest_energy": wave.same_rest_energy,
            "no_separate_particle_coordinate": wave.no_separate_particle_coordinate,
            "detector_coupling": wave.detector_coupling,
        },
        "gate": {
            "stationary_localized_energy": sector.volume_hessian > 0.0,
            "charge_converges_to_one": quadrature[-1]["charge_error"] < 1.0e-6,
            "energy_converges_to_analytic_value": quadrature[-1]["normalized_energy_error"] < 1.0e-8,
            "particle_and_wave_use_same_rest_energy": wave.same_rest_energy,
            "particle_and_wave_use_same_field_configuration": wave.no_separate_particle_coordinate,
        },
        "next_result": "activate_the_shared_gauge_spin_and_L2_L4_terms_and_solve_the_electron_proton_neutron_stationary_sectors_with_one_Theta",
        "nonclaims": {
            "electron_solution_obtained": False,
            "spin_half_and_electric_charge_derived": False,
            "Bell_CHSH_reduction_completed": False,
            "complete_Theta_calibrated": False,
            "cross_scale_invariant_closure_completed": False,
        },
        "terminal": True,
    }
    if not all(result["gate"].values()):
        raise RuntimeError("S-one particle/wave identity gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(result))
    return result


def main() -> int:
    result = run()
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
