#!/usr/bin/env python3
"""Test reflection positivity of the selected relative-sheet covariance.

The full exponential heat response is positive pointwise in Euclidean
momentum, but that alone does not imply a positive real-time Hilbert space.
For the child-selected gap, this calculation evaluates every Fourier mode on
a 64-point Euclidean time circle, inverts the kernel, and applies the Gaussian
Osterwalder--Schrader reflection matrix test.

The tested covariance contains the anomaly-induced relative heat response.
If it fails, the result locates the missing induced metric/fermionic observable
block rather than labelling pointwise Euclidean positivity as causality.
"""

from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from scripts.run_s_one_full_exponential_relative_kernel import (  # noqa: E402
    _one_kernel,
)


OUTPUT = ROOT / "results/nsc-1-s-one-relative-reflection-test.json"
KERNEL_RUNNER = ROOT / "scripts/run_s_one_full_exponential_relative_kernel.py"
INPUTS = (
    ROOT / "results/nsc-1-s-one-full-exponential-relative-kernel.json",
    ROOT / "results/nsc-1-s-one-child-orientation.json",
    ROOT / "results/nsc-1-s-one-boundary-retarded-pole.json",
    KERNEL_RUNNER,
)
RESOLUTION = 64
CUTOFF = 8.0
MASS_RATIO = sqrt(1006.0 / 1015.0)


def _authenticate() -> list[dict[str, object]]:
    authenticated = []
    for path in INPUTS:
        raw = path.read_bytes()
        if path.suffix == ".json":
            item = json.loads(raw)
            if item.get("terminal") is not True:
                raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
            artifact_id = item["artifact_id"]
        else:
            artifact_id = None
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": artifact_id,
            }
        )
    return authenticated


def _kernel_for_mode(mode: int) -> dict[str, float | int]:
    return _one_kernel((RESOLUTION, MASS_RATIO, mode / CUTOFF))


def _reflection_spectrum(
    covariance_time: np.ndarray, maximum_positive_time: int
) -> dict[str, object]:
    times = np.arange(1, maximum_positive_time + 1)
    matrix = np.array(
        [
            [covariance_time[(left + right) % RESOLUTION] for right in times]
            for left in times
        ],
        dtype=float,
    )
    matrix = 0.5 * (matrix + matrix.T)
    eigenvalues = np.linalg.eigvalsh(matrix)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    tolerance = 1.0e-10 * scale
    return {
        "maximum_positive_time": maximum_positive_time,
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "negative_eigenvalue_count": int(np.count_nonzero(eigenvalues < -tolerance)),
        "tolerance": tolerance,
        "positive_semidefinite": bool(eigenvalues[0] >= -tolerance),
        "eigenvalues": [float(value) for value in eigenvalues],
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("relative reflection-test output already exists")
    authenticated = _authenticate()

    nonnegative_modes = list(range(RESOLUTION // 2 + 1))
    worker_count = min(len(nonnegative_modes), os.cpu_count() or 1)
    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        rows = list(executor.map(_kernel_for_mode, nonnegative_modes))
    rows.sort(key=lambda item: int(item["external_mode"]))

    kernel_by_mode = {
        int(row["external_mode"]): float(row["normalized_kernel"])
        for row in rows
    }
    fourier_kernel = np.empty(RESOLUTION, dtype=float)
    for index, frequency in enumerate(
        np.fft.fftfreq(RESOLUTION, d=1.0 / RESOLUTION)
    ):
        fourier_kernel[index] = kernel_by_mode[abs(int(frequency))]
    covariance_momentum = 1.0 / fourier_kernel
    covariance_time_complex = np.fft.ifft(covariance_momentum)
    covariance_time = covariance_time_complex.real
    imaginary_residual = float(np.max(np.abs(covariance_time_complex.imag)))

    reflection = [
        _reflection_spectrum(covariance_time, maximum)
        for maximum in (4, 8, 16, 31)
    ]
    full = reflection[-1]
    all_reflection_positive = all(
        bool(item["positive_semidefinite"]) for item in reflection
    )

    # A nearest-neighbour massive lattice scalar has a known positive transfer
    # matrix.  It checks the reflection indexing and finite-circle convention.
    frequencies = np.fft.fftfreq(RESOLUTION, d=1.0 / RESOLUTION)
    control_kernel = 1.0 + 4.0 * np.sin(np.pi * frequencies / RESOLUTION) ** 2
    control_time = np.fft.ifft(1.0 / control_kernel).real
    control_reflection = [
        _reflection_spectrum(control_time, maximum)
        for maximum in (4, 8, 16, 31)
    ]
    control_positive = all(
        bool(item["positive_semidefinite"]) for item in control_reflection
    )

    record = {
        "artifact_id": "NSC-1-S-ONE-RELATIVE-REFLECTION-TEST",
        "schema": "NSC-1-S-ONE-RELATIVE-REFLECTION-TEST-v1",
        "classification": (
            "selected_relative_heat_covariance_passes_the_finite_OS_test"
            if all_reflection_positive
            else "selected_relative_heat_covariance_fails_the_finite_OS_test_and_requires_the_complete_induced_metric_observable"
        ),
        "authenticated_inputs": authenticated,
        "method": {
            "field": "relative_two_sheet_resolution_delta",
            "mass_squared_over_cutoff_squared": "1006/1015",
            "Euclidean_time_resolution": RESOLUTION,
            "cutoff": CUTOFF,
            "Fourier_modes": [0, RESOLUTION // 2],
            "kernel": "complete_exponential_heat_Hessian_K_minus(q)",
            "covariance": "C_minus(q)=1/K_minus(q)",
            "time_covariance": "inverse_discrete_Fourier_transform_of_C_minus(q)",
            "reflection_matrix": "M_ij=C_minus(t_i+t_j) for positive t_i,t_j",
            "parallel_workers": worker_count,
        },
        "kernel": {
            "minimum": float(np.min(fourier_kernel)),
            "maximum": float(np.max(fourier_kernel)),
            "all_positive": bool(np.all(fourier_kernel > 0.0)),
            "rows": rows,
        },
        "covariance": {
            "momentum_minimum": float(np.min(covariance_momentum)),
            "momentum_maximum": float(np.max(covariance_momentum)),
            "time_values": [float(value) for value in covariance_time],
            "maximum_imaginary_DFT_residual": imaginary_residual,
        },
        "reflection_spectra": reflection,
        "free_lattice_scalar_control": {
            "kernel": "1+4*sin(pi*n/N)^2",
            "reflection_spectra": control_reflection,
            "passes": control_positive,
        },
        "causal_result": {
            "finite_OS_reflection_positive": all_reflection_positive,
            "minimum_full_reflection_eigenvalue": full["minimum_eigenvalue"],
            "negative_full_reflection_eigenvalues": full[
                "negative_eigenvalue_count"
            ],
            "owner_if_nonpass": "the_relative_heat_response_is_not_the_complete_gauge_invariant_induced_metric_covariance",
            "next_equation": (
                "continue_the_OS_positive_relative_covariance_to_the_retarded_kernel"
                if all_reflection_positive
                else "add_the_fixed_induced_Einstein_and_anomaly_consistent_fermionic_metric_blocks_then_repeat_the_same_OS_test"
            ),
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "all_Fourier_kernels_positive": bool(np.all(fourier_kernel > 0.0)),
            "real_time_covariance_reconstructed": imaginary_residual < 1.0e-12,
            "finite_OS_reflection_test_executed": True,
            "free_lattice_scalar_control_passes": control_positive,
            "finite_OS_reflection_positive": all_reflection_positive,
            "complete_curved_metric_covariance_tested": False,
            "conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "finite_circle_is_a_continuum_OS_theorem": False,
            "relative_scalar_is_the_complete_metric_tensor": False,
            "curved_black_child_retarded_kernel_proved": False,
            "electron_spectrum_predicted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key in {
            "finite_OS_reflection_positive",
            "complete_curved_metric_covariance_tested",
            "conviction_proof_seed_completed",
        }:
            continue
        if value is not True:
            raise RuntimeError(f"relative reflection test failed to execute: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
