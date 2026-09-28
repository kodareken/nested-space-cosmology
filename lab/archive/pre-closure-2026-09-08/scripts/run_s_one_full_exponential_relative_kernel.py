#!/usr/bin/env python3
"""Evaluate the full exponential flat two-sheet relative-resolution kernel.

The leading doubled-FLRW heat-kernel truncation has no generically stable
positive-de-Sitter parameter overlap.  This calculation goes back to the
untruncated profile fixed by the project, ``Tr exp(-D^2/Lambda^2)``, and asks
the next direct question: is the relative two-sheet conformal mode positive in
the actual exponential spectrum?

An analytic four-dimensional zero-external-momentum result is combined with a
Fourier/Gauss-Legendre finite-momentum calculation.  The determinant part is
exactly invariant under the traceless relative sandwich on the finite spectral
space, so this controlled relative kernel is the anomaly-induced heat response.
The result is Euclidean and flat; curved retarded Lorentzian propagation and a
stationary physical value of Phi/Lambda remain separate next equations.
"""

from __future__ import annotations

import os

# Use the complete CPU deliberately without nesting a threaded BLAS inside
# every process.
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import exp, pi
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-full-exponential-relative-kernel.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-doubled-flrw-closure.json",
    ROOT / "results/nsc-1-s-one-local-two-sheet-anomaly.json",
    ROOT / "results/nsc-1-s-one-spectral-profile-closure.json",
)
MASS_RATIOS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
EXTERNAL_RATIOS = (0.0, 0.5, 1.0, 1.5, 2.0)
RESOLUTIONS = (32, 48, 64)
MODES_PER_CUTOFF = 8
RADIAL_NODES = 20
RADIAL_MAX_OVER_CUTOFF = 5.0
COARSE_EPSILON = 2.0e-3
FINE_EPSILON = 1.0e-3


def _authenticate() -> list[dict[str, object]]:
    authenticated = []
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item["artifact_id"],
            }
        )
    return authenticated


def _gamma_matrices() -> tuple[tuple[np.ndarray, ...], np.ndarray]:
    sigma_1 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    sigma_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
    sigma_3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    identity_2 = np.eye(2, dtype=np.complex128)
    gamma = (
        np.kron(sigma_1, sigma_1),
        np.kron(sigma_1, sigma_2),
        np.kron(sigma_1, sigma_3),
        np.kron(sigma_2, identity_2),
    )
    gamma_5 = np.kron(sigma_3, identity_2)
    identity_4 = np.eye(4, dtype=np.complex128)
    for first in range(4):
        for second in range(4):
            expected = 2.0 * identity_4 if first == second else 0.0 * identity_4
            if not np.allclose(
                gamma[first] @ gamma[second] + gamma[second] @ gamma[first],
                expected,
                rtol=0.0,
                atol=1.0e-13,
            ):
                raise RuntimeError("Euclidean Clifford algebra construction failed")
        if not np.allclose(
            gamma_5 @ gamma[first] + gamma[first] @ gamma_5,
            np.zeros((4, 4), dtype=np.complex128),
            rtol=0.0,
            atol=1.0e-13,
        ):
            raise RuntimeError("gamma_5 anticommutation failed")
    return gamma, gamma_5


def _spectral_momentum(resolution: int) -> np.ndarray:
    momenta = np.fft.fftfreq(resolution, d=1.0 / resolution)
    points = np.arange(resolution)
    fourier = np.exp(
        2.0j * np.pi * np.outer(points, points) / resolution
    ) / np.sqrt(resolution)
    matrix = fourier @ np.diag(momenta) @ fourier.conjugate().T
    return 0.5 * (matrix + matrix.conjugate().T)


def _relative_scalings(
    resolution: int, external_mode: int, epsilon: float
) -> tuple[np.ndarray, np.ndarray]:
    coordinate = 2.0 * np.pi * np.arange(resolution) / resolution
    delta = epsilon * np.cos(external_mode * coordinate)
    identity_4 = np.eye(4, dtype=np.complex128)
    parent = np.kron(identity_4, np.diag(np.exp(-delta / 2.0)))
    child = np.kron(identity_4, np.diag(np.exp(delta / 2.0)))
    return parent, child


def _heat_trace(
    doubled_dirac: np.ndarray,
    parent_scale: np.ndarray,
    child_scale: np.ndarray,
    cutoff: float,
) -> float:
    zeros = np.zeros_like(parent_scale)
    scale = np.block([[parent_scale, zeros], [zeros, child_scale]])
    scaled = scale @ doubled_dirac @ scale
    eigenvalues = np.linalg.eigvalsh(scaled)
    return float(np.sum(np.exp(-(eigenvalues / cutoff) ** 2)).real)


def _one_kernel(task: tuple[int, float, float]) -> dict[str, float | int]:
    resolution, mass_ratio, external_ratio = task
    cutoff = resolution / MODES_PER_CUTOFF
    external_mode = int(round(external_ratio * cutoff))
    if abs(external_mode / cutoff - external_ratio) > 1.0e-15:
        raise RuntimeError("external mode is not common across resolutions")

    gamma, gamma_5 = _gamma_matrices()
    momentum = _spectral_momentum(resolution)
    identity_x = np.eye(resolution, dtype=np.complex128)
    identity_sheet_spin_x = np.eye(4 * resolution, dtype=np.complex128)
    link_block = mass_ratio * cutoff * np.kron(gamma_5, identity_x)

    epsilon_values = (
        0.0,
        COARSE_EPSILON,
        -COARSE_EPSILON,
        FINE_EPSILON,
        -FINE_EPSILON,
    )
    scales = {
        epsilon: _relative_scalings(resolution, external_mode, epsilon)
        for epsilon in epsilon_values
    }

    nodes, weights = leggauss(RADIAL_NODES)
    radial_max = RADIAL_MAX_OVER_CUTOFF * cutoff
    radial_values = radial_max * (nodes + 1.0) / 2.0
    radial_weights = radial_max * weights / 2.0
    integrated = {epsilon: 0.0 for epsilon in epsilon_values}

    for radial, radial_weight in zip(radial_values, radial_weights, strict=True):
        geometric = np.kron(gamma[0], momentum) + radial * np.kron(
            gamma[1], identity_x
        )
        doubled = np.block(
            [[geometric, link_block], [link_block, geometric]]
        )
        if not np.allclose(
            doubled, doubled.conjugate().T, rtol=0.0, atol=1.0e-12
        ):
            raise RuntimeError("doubled Dirac matrix lost self-adjointness")
        measure = radial_weight * radial**2 / (2.0 * np.pi**2)
        for epsilon in epsilon_values:
            parent, child = scales[epsilon]
            integrated[epsilon] += measure * _heat_trace(
                doubled, parent, child, cutoff
            )

    base = integrated[0.0]
    coarse = (
        integrated[COARSE_EPSILON]
        - 2.0 * base
        + integrated[-COARSE_EPSILON]
    ) / COARSE_EPSILON**2
    fine = (
        integrated[FINE_EPSILON]
        - 2.0 * base
        + integrated[-FINE_EPSILON]
    ) / FINE_EPSILON**2
    richardson = (4.0 * fine - coarse) / 3.0
    finite_difference_error = abs(richardson - fine)

    # The Fourier trace is over a circle of coordinate length 2*pi.  A
    # nonconstant cosine has mean square 1/2, so multiply its amplitude
    # Hessian by two to report the Fourier-kernel coefficient.
    fourier_normalization = 1.0 if external_mode == 0 else 2.0
    normalization = fourier_normalization / (2.0 * np.pi * cutoff**4)
    normalized_kernel = richardson * normalization
    normalized_error = finite_difference_error * normalization

    # The determinant is exactly invariant because the relative scaling has
    # unit determinant.  Record a bounded numerical witness on one radial
    # block; the analytic identity is stored in the result owner.
    determinant_residual = 0.0
    if mass_ratio == 0.5 and external_ratio == 1.0:
        radial = cutoff / 3.0
        geometric = np.kron(gamma[0], momentum) + radial * np.kron(
            gamma[1], identity_x
        )
        doubled = np.block(
            [[geometric, link_block], [link_block, geometric]]
        )
        parent, child = _relative_scalings(
            resolution, external_mode, 0.125
        )
        zeros = np.zeros_like(identity_sheet_spin_x)
        scale = np.block([[parent, zeros], [zeros, child]])
        eigenvalues_0 = np.linalg.eigvalsh(doubled)
        eigenvalues_1 = np.linalg.eigvalsh(scale @ doubled @ scale)
        if min(np.min(np.abs(eigenvalues_0)), np.min(np.abs(eigenvalues_1))) <= 0.0:
            raise RuntimeError("determinant witness unexpectedly became singular")
        logdet_0 = float(np.sum(np.log(np.abs(eigenvalues_0))))
        logdet_1 = float(np.sum(np.log(np.abs(eigenvalues_1))))
        determinant_residual = abs(logdet_1 - logdet_0)

    return {
        "resolution": resolution,
        "cutoff": float(cutoff),
        "mass_over_cutoff": float(mass_ratio),
        "external_momentum_over_cutoff": float(external_ratio),
        "external_mode": external_mode,
        "normalized_kernel": float(normalized_kernel),
        "finite_difference_error": float(normalized_error),
        "determinant_residual": float(determinant_residual),
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("full exponential relative-kernel output already exists")
    authenticated = _authenticate()

    tasks = [
        (resolution, mass_ratio, external_ratio)
        for resolution in RESOLUTIONS
        for mass_ratio in MASS_RATIOS
        for external_ratio in EXTERNAL_RATIOS
    ]
    worker_count = min(len(tasks), os.cpu_count() or 1)
    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        rows = list(executor.map(_one_kernel, tasks))
    rows.sort(
        key=lambda item: (
            int(item["resolution"]),
            float(item["mass_over_cutoff"]),
            float(item["external_momentum_over_cutoff"]),
        )
    )

    by_key = {
        (
            int(row["resolution"]),
            float(row["mass_over_cutoff"]),
            float(row["external_momentum_over_cutoff"]),
        ): row
        for row in rows
    }
    finest = [row for row in rows if row["resolution"] == RESOLUTIONS[-1]]
    convergence = []
    for mass_ratio in MASS_RATIOS:
        for external_ratio in EXTERNAL_RATIOS:
            medium = by_key[(RESOLUTIONS[-2], mass_ratio, external_ratio)]
            fine = by_key[(RESOLUTIONS[-1], mass_ratio, external_ratio)]
            scale = max(abs(float(fine["normalized_kernel"])), 1.0e-15)
            convergence.append(
                {
                    "mass_over_cutoff": mass_ratio,
                    "external_momentum_over_cutoff": external_ratio,
                    "medium_to_fine_relative_difference": abs(
                        float(fine["normalized_kernel"])
                        - float(medium["normalized_kernel"])
                    )
                    / scale,
                }
            )

    mu_squared, integration_variable = sp.symbols(
        "mu_squared x", nonnegative=True, real=True
    )
    dimensionless_radial_integrand = (
        integration_variable**2
        * (integration_variable + mu_squared - 1)
        * sp.exp(-integration_variable - mu_squared)
    )
    radial_integral = sp.integrate(
        dimensionless_radial_integrand,
        (integration_variable, 0, sp.oo),
    )
    expected_radial_integral = 2 * (mu_squared + 2) * sp.exp(
        -mu_squared
    )
    analytic_kernel = 4 * (mu_squared + 2) * sp.exp(-mu_squared) / sp.pi**2

    zero_mode_comparison = []
    for mass_ratio in MASS_RATIOS:
        fine = by_key[(RESOLUTIONS[-1], mass_ratio, 0.0)]
        expected = 4.0 * (mass_ratio**2 + 2.0) * exp(
            -(mass_ratio**2)
        ) / pi**2
        zero_mode_comparison.append(
            {
                "mass_over_cutoff": mass_ratio,
                "analytic_kernel": expected,
                "numerical_kernel": fine["normalized_kernel"],
                "relative_error": abs(
                    float(fine["normalized_kernel"]) - expected
                )
                / expected,
            }
        )

    minimum = min(float(row["normalized_kernel"]) for row in finest)
    maximum = max(float(row["normalized_kernel"]) for row in finest)
    max_fd_error = max(float(row["finite_difference_error"]) for row in finest)
    max_convergence = max(
        float(item["medium_to_fine_relative_difference"])
        for item in convergence
    )
    max_zero_error = max(float(item["relative_error"]) for item in zero_mode_comparison)
    max_determinant_residual = max(
        float(row["determinant_residual"]) for row in rows
    )

    record = {
        "artifact_id": "NSC-1-S-ONE-FULL-EXPONENTIAL-RELATIVE-KERNEL",
        "schema": "NSC-1-S-ONE-FULL-EXPONENTIAL-RELATIVE-KERNEL-v1",
        "classification": "the_untruncated_exponential_spectrum_has_a_positive_flat_Euclidean_relative_sheet_kernel_on_the_scanned_mass_momentum_box",
        "authenticated_inputs": authenticated,
        "operator": {
            "geometric_block": "gamma_1*p_1+gamma_2*p_transverse",
            "doubled_Dirac": "[[D_g,gamma_5*Phi],[gamma_5*Phi,D_g]]",
            "relative_rescaling": "V_delta=diag(exp(-delta(x)/2),exp(+delta(x)/2)); D_delta=V_delta*D*V_delta",
            "profile": "Tr exp(-D_delta^2/Lambda^2)",
            "determinant_identity": "det(V_delta*D*V_delta)=det(D) because det(V_delta)=1",
            "fermionic_relative_Hessian_in_controlled_finite_spectrum": "0",
        },
        "analytic_zero_external_momentum": {
            "radial_integrand": str(dimensionless_radial_integrand),
            "radial_integral": str(radial_integral),
            "expected_radial_integral": str(expected_radial_integral),
            "four_spin_component_kernel_over_Lambda4": str(analytic_kernel),
            "strictly_positive_for_all_nonnegative_mass_squared": True,
            "numerical_comparison": zero_mode_comparison,
        },
        "finite_momentum_scan": {
            "mass_over_cutoff": list(MASS_RATIOS),
            "external_momentum_over_cutoff": list(EXTERNAL_RATIOS),
            "resolutions": list(RESOLUTIONS),
            "modes_per_cutoff": MODES_PER_CUTOFF,
            "transverse_Gauss_Legendre_nodes": RADIAL_NODES,
            "transverse_radial_max_over_cutoff": RADIAL_MAX_OVER_CUTOFF,
            "coarse_epsilon": COARSE_EPSILON,
            "fine_epsilon": FINE_EPSILON,
            "parallel_workers": worker_count,
            "finest_minimum_kernel": minimum,
            "finest_maximum_kernel": maximum,
            "maximum_finite_difference_error": max_fd_error,
            "maximum_medium_fine_relative_difference": max_convergence,
            "rows": rows,
            "convergence": convergence,
        },
        "causal_result": {
            "leading_truncation_obstruction_retested_at_full_profile": True,
            "zero_and_sampled_finite_Euclidean_momentum_relative_modes_positive": True,
            "same_Phi_retained": "mass_gap_outside_self_energy_relative_metric_kernel",
            "repair": "return_from_a0_plus_a2_to_the_already_fixed_complete_exponential_profile",
            "new_parameter_added": False,
            "next_equation": "analytically_continue_the_flat_relative_kernel_and_compute_the_curved_two_sheet_retarded_spectrum_then_solve_the_stationary_Phi_over_Lambda_equation",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "analytic_radial_integral_exact": sp.simplify(
                radial_integral - expected_radial_integral
            )
            == 0,
            "analytic_zero_mode_positive": True,
            "all_scanned_finite_momentum_kernels_positive": minimum > 0.0,
            "zero_mode_matches_analytic_below_5e_5": max_zero_error < 5.0e-5,
            "finite_difference_error_below_1e_5": max_fd_error < 1.0e-5,
            "medium_fine_difference_below_5e_3": max_convergence < 5.0e-3,
            "relative_determinant_invariant_below_1e_9": max_determinant_residual
            < 1.0e-9,
            "curved_retarded_kernel_completed": False,
            "stationary_Phi_over_Lambda_solved": False,
            "conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "finite_grid_is_a_continuum_all_momentum_theorem": False,
            "Euclidean_positivity_is_Lorentzian_causality": False,
            "complete_gauge_invariant_graviton_covariance_tested": False,
            "electron_mass_predicted": False,
            "black_child_transition_solved": False,
            "Xi_NSC_predicted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key not in {
            "curved_retarded_kernel_completed",
            "stationary_Phi_over_Lambda_solved",
            "conviction_proof_seed_completed",
        } and value is not True:
            raise RuntimeError(f"full exponential relative-kernel gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
