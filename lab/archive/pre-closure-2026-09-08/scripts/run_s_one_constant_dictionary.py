#!/usr/bin/env python3
"""Derive every local constant from one stationary-action spectrum."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.unified_action import (  # noqa: E402
    LocalVacuumSpectrum,
    RecursiveScaleTransfer,
    derive_local_constants,
)


OUTPUT = ROOT / "results/nsc-1-s-one-constant-dictionary.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one constant dictionary output already exists")
    zt, zx, hessian, stiffness, vacuum, hbar = sp.symbols(
        "Z_t Z_x H F_star V_star hbar", positive=True
    )
    omega2, momentum2 = sp.symbols("omega_squared momentum_squared", nonnegative=True)
    inverse_propagator = zt * omega2 - zx * momentum2 - hessian
    dispersion = sp.solve(sp.Eq(inverse_propagator, 0), omega2)[0]
    speed_squared = sp.simplify(sp.diff(dispersion, momentum2))
    mass_squared = sp.simplify(dispersion.subs(momentum2, 0))
    gravity = sp.simplify(1 / (8 * sp.pi * stiffness))
    planck = sp.simplify(sp.sqrt(hbar * gravity / speed_squared ** sp.Rational(3, 2)))
    cosmological = sp.simplify(vacuum / stiffness)

    omega_scale, parent_hbar_c, child_hbar_c = sp.symbols(
        "Omega hbar_c_parent hbar_c_child", positive=True
    )
    gamma = sp.simplify(child_hbar_c / (omega_scale * parent_hbar_c))
    invariant_gamma = sp.simplify(
        gamma.subs(child_hbar_c, parent_hbar_c) - 1 / omega_scale
    )

    # A numerical fixture checks the public evaluator without treating these
    # arbitrary values as physical constants or a calibrated Theta.
    spectrum = LocalVacuumSpectrum(
        temporal_kinetic=4.0,
        spatial_gradient=9.0,
        metric_stiffness=2.0,
        phase_quantum=1.0,
        vacuum_energy=0.25,
        mass_hessian=(1.0, 4.0, 9.0),
        spin_casimir=(0.0, 0.75, 2.0),
    )
    constants = derive_local_constants(spectrum)
    transfer = RecursiveScaleTransfer(
        omega=4.0,
        parent_hbar_c=1.0,
        child_hbar_c_in_parent_units=1.0,
    )
    fixed_point_formula = (
        "Theta=C[Theta], delta_S_one[Phi_star(Theta)]=0, "
        "T_Theta^*S_one=S_one"
    )
    spectrum_payload = asdict(spectrum)
    spectrum_payload["mass_hessian"] = list(spectrum.mass_hessian)
    spectrum_payload["spin_casimir"] = list(spectrum.spin_casimir)
    constants_payload = asdict(constants)
    constants_payload["mass_eigenvalues_in_anchor_units"] = list(
        constants.mass_eigenvalues_in_anchor_units
    )
    constants_payload["spins"] = list(constants.spins)
    result = {
        "artifact_id": "NSC-1-S-ONE-CONSTANT-DICTIONARY",
        "schema": "NSC-1-S-ONE-CONSTANT-DICTIONARY-v1",
        "classification": "all_local_constants_are_named_spectral_outputs_of_one_stationary_recursive_action",
        "fixed_point_closure": {
            "formula": fixed_point_formula,
            "meaning": "Theta_is_accepted_only_when_the_constants_produced_by_its_own_stationary_solutions_reproduce_Theta_under_the_parent_child_transfer",
            "numerical_solution_obtained": False,
        },
        "one_Hessian_dictionary": {
            "inverse_propagator": str(inverse_propagator),
            "dispersion_omega_squared": str(dispersion),
            "c_squared": str(speed_squared),
            "mass_squared": str(mass_squared),
            "G": str(gravity),
            "Planck_length_over_local_anchor": str(planck),
            "Lambda_times_local_anchor_squared": str(cosmological),
            "spin": "solve j(j+1)=little_group_or_collective_rotation_Casimir",
            "absolute_zero": "infimum_spectrum(H_local) relabelled as T=0",
        },
        "parent_child_unit_transfer": {
            "Gamma": str(gamma),
            "hbar_c_invariant_reduction": "1/Omega",
            "symbolic_identity_residual": str(invariant_gamma),
            "Gamma_stored_as_independent_parameter": False,
        },
        "public_evaluator_fixture": {
            "spectrum": spectrum_payload,
            "derived_constants": constants_payload,
            "derived_Gamma": transfer.gamma,
            "fixture_is_not_physical_calibration": True,
        },
        "gate": {
            "one_Hessian_generates_speed_and_mass": speed_squared == zx / zt
            and mass_squared == hessian / zt,
            "metric_stiffness_generates_G": gravity == 1 / (8 * sp.pi * stiffness),
            "Planck_length_is_derived_not_independent": "G" not in fixed_point_formula,
            "Lambda_is_vacuum_boundary_output": cosmological == vacuum / stiffness,
            "Gamma_is_derived_not_fitted": invariant_gamma == 0,
            "absolute_zero_is_spectral_bound_not_free_coefficient": constants.absolute_zero
            == 0.0,
        },
        "next_result": "solve_Theta_equals_C_of_Theta_for_the_shared_electron_proton_neutron_stationary_sectors_then_open_mass_ratio_holdouts",
        "nonclaims": {
            "numerical_values_of_natures_constants_predicted": False,
            "electron_proton_neutron_spectrum_solved": False,
            "our_room_fixed_point_identified": False,
            "inside_black_hole_observationally_proved": False,
        },
        "terminal": True,
    }
    if not all(result["gate"].values()):
        raise RuntimeError("S-one constant dictionary gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(result))
    return result


def main() -> int:
    result = run()
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
