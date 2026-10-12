#!/usr/bin/env python3
"""Solve the exact fixed point created by an infinite nest of identical rooms."""

from __future__ import annotations

import json
from pathlib import Path
import sys

from scipy.integrate import quad
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.unified_action import RecursiveOutsideKernel  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-recursive-outside-kernel.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("S-one recursive outside-kernel output already exists")

    gamma, local, coupling = sp.symbols("Gamma K b", positive=True)
    equation = sp.Eq(gamma, local - coupling**2 / gamma)
    solutions = sp.solve(equation, gamma)
    physical = (local + sp.sqrt(local**2 - 4 * coupling**2)) / 2
    residual = sp.simplify(physical - local + coupling**2 / physical)
    uncoupled_limit = sp.limit(physical, coupling, 0, dir="+")

    kernel = RecursiveOutsideKernel(boundary_coupling=1.0)
    samples = []
    for local_value in (3.0, 4.0, 8.0, 16.0):
        effective = kernel.effective_inverse(local_value)
        self_energy = kernel.outside_self_energy(local_value)
        fixed_residual = kernel.fixed_point_residual(local_value)
        samples.append(
            {
                "local_inverse": local_value,
                "effective_inverse": effective.real,
                "outside_self_energy": self_energy.real,
                "reconstructed_local_inverse": (effective + self_energy).real,
                "fixed_point_residual": abs(fixed_residual),
            }
        )

    density_integral, density_error = quad(
        kernel.surface_spectral_density,
        -2.0,
        2.0,
        epsabs=1.0e-13,
        epsrel=1.0e-13,
    )
    second_moment, moment_error = quad(
        lambda value: value * value * kernel.surface_spectral_density(value),
        -2.0,
        2.0,
        epsabs=1.0e-13,
        epsrel=1.0e-13,
    )

    record = {
        "artifact_id": "NSC-1-S-ONE-RECURSIVE-OUTSIDE-KERNEL",
        "schema": "NSC-1-S-ONE-RECURSIVE-OUTSIDE-KERNEL-v1",
        "classification": "recursive_self_equality_replaces_an_arbitrary_outside_source_with_one_exact_self_consistent_spectrum",
        "full_functional_fixed_point": "exp(i*Gamma[Phi]/hbar)=integral Dchi exp(i*(S_room[Phi]+S_link[Phi,chi]+Gamma[T_Theta chi])/hbar)",
        "quadratic_derivation": {
            "infinite_nest": "identical room kernels K coupled to the next room by b",
            "one_step_Schur_complement": "Gamma_next=K-b^2/Gamma_tail",
            "self_equal_tail": "Gamma_tail=Gamma_next=Gamma",
            "fixed_point_equation": str(equation),
            "solutions": [str(value) for value in solutions],
            "causal_uncoupled_branch": str(physical),
            "exact_residual": str(residual),
            "uncoupled_limit": str(uncoupled_limit),
        },
        "normalized_unit_link_result": {
            "spectral_support": [-2.0, 2.0],
            "surface_spectral_density": "sqrt(4-E^2)/(2*pi) for abs(E)<2",
            "integrated_spectral_weight": density_integral,
            "integration_error": density_error,
            "second_moment": second_moment,
            "second_moment_error": moment_error,
            "expected_second_moment": 1.0,
            "euclidean_samples": samples,
        },
        "binding_to_known_outside_geometry": {
            "source_artifact": "results/nsc-1-s-one-outside-black-bind.json",
            "RS2_open_term": "the projected Weyl E_mu_nu is undetermined until the outside bulk is solved",
            "recursive_closure": "Pi_outside=-B_dagger*(T_Theta^*Gamma^(2))^-1*B",
            "meaning": "the next room's propagator determines the current room's apparently dark response",
        },
        "observable_routes_from_the_same_kernel": {
            "particles_and_waves": "isolated zeros and pole residues after promoting K and b to the field matrix",
            "dark_response": "the continuous nonlocal spectral contribution over finite momentum",
            "local_constants": "zero-momentum value and principal high-momentum coefficients",
            "black_to_child": "nonlinear saddles of the untruncated functional fixed-point equation",
        },
        "new_result": {
            "arbitrary_outside_function_removed_at_quadratic_scalar_level": True,
            "one_link_coefficient_generates_complete_infinite_tail": True,
            "outside_spectrum_has_unit_total_weight": True,
            "same_equation_repeats_at_every_depth": True,
        },
        "next_result": "promote_K_and_B_to_the_spherical_metric_boundary_matrix_and_test_whether_the_fixed_point_self_energy_reproduces_the_authenticated_black_universe_E_mu_nu_without_an_arbitrary_radial_function",
        "nonclaims": {
            "tensor_bulk_equations_solved": False,
            "particle_masses_predicted": False,
            "dark_matter_lensing_predicted": False,
            "black_universe_stability_proved": False,
            "observed_universe_identified": False,
        },
        "gate": {
            "symbolic_fixed_point_exact": residual == 0,
            "correct_uncoupled_branch": uncoupled_limit == local,
            "numeric_fixed_point_residual_below_1e_14": max(
                item["fixed_point_residual"] for item in samples
            )
            < 1.0e-14,
            "spectral_weight_normalized": abs(density_integral - 1.0) < 1.0e-12,
            "spectral_second_moment_matches_link_squared": abs(second_moment - 1.0)
            < 1.0e-12,
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("recursive outside-kernel gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
