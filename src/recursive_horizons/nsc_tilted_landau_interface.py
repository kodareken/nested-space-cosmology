"""Weighted transmitting maps for the retained tilted Landau Cauchy data.

The spatial boundary owner fixes throat transmission in a common coordinate
spin frame.  The Landau endpoint fixes the positive spin metric
``B_L=exp(-eta*sigma2)``.  These data constrain a channel kernel ``K`` by

    K^dagger G_L K = G_0,

but do not select its canonical unitary ``V``.  This module implements the
entire channel-block family ``K=G_L^-1/2 V G_0^1/2`` and finite witnesses.  It
does not identify ``B_L^1/2`` with the Cauchy map or select a history.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import cosh, isfinite, sinh

import numpy as np


I2 = np.eye(2, dtype=complex)
SIGMA2 = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)


def landau_spin_factors(rapidity: float):
    """Return B_L, its positive square root and its inverse square root."""
    if not isfinite(rapidity):
        raise ValueError("finite Landau rapidity required")
    metric = cosh(rapidity) * I2 - sinh(rapidity) * SIGMA2
    half = cosh(rapidity / 2) * I2 - sinh(rapidity / 2) * SIGMA2
    inverse_half = cosh(rapidity / 2) * I2 + sinh(rapidity / 2) * SIGMA2
    return metric, half, inverse_half


def full_covariance(blocks: np.ndarray) -> np.ndarray:
    """Embed frequency-local 2x2 covariances in one channel matrix."""
    blocks = np.asarray(blocks, dtype=complex)
    if blocks.ndim != 3 or blocks.shape[1:] != (2, 2):
        raise ValueError("frequency-local 2x2 covariance blocks required")
    count = len(blocks)
    result = np.zeros((2 * count, 2 * count), dtype=complex)
    for index, block in enumerate(blocks):
        result[2 * index:2 * index + 2, 2 * index:2 * index + 2] = block
    return result


def frequency_fourier_unitary(count: int) -> np.ndarray:
    """Dense all-frequency witness, acting identically on both spin entries."""
    if isinstance(count, bool) or not isinstance(count, int) or count < 2:
        raise ValueError("at least two frequency nodes required")
    indices = np.arange(count, dtype=float)
    fourier = np.exp(-2j * np.pi * np.outer(indices, indices) / count) / np.sqrt(count)
    return np.kron(fourier, I2)


def frequency_offblock_norm(matrix: np.ndarray) -> float:
    """Frobenius norm outside equal-frequency 2x2 blocks."""
    matrix = np.asarray(matrix)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] % 2:
        raise ValueError("even square channel matrix required")
    masked = np.array(matrix, copy=True)
    for index in range(matrix.shape[0] // 2):
        masked[2 * index:2 * index + 2, 2 * index:2 * index + 2] = 0
    return float(np.linalg.norm(masked))


@dataclass(frozen=True)
class TransmittingTiltedLandauInterface:
    """Admissible weighted-isometry family on one retained channel block."""

    rapidity: float
    tolerance: float = 3e-11

    def __post_init__(self):
        if not isfinite(self.rapidity):
            raise ValueError("finite rapidity required")
        if not isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError("positive finite tolerance required")

    def metrics(self, weights: np.ndarray):
        weights = np.asarray(weights, dtype=float)
        if weights.ndim != 1 or len(weights) < 2 or not np.isfinite(weights).all():
            raise ValueError("finite one-dimensional quadrature weights required")
        if np.min(weights) <= 0:
            raise ValueError("positive quadrature weights required")
        spin_metric, spin_half, spin_inverse_half = landau_spin_factors(self.rapidity)
        root = np.diag(np.sqrt(weights))
        inverse_root = np.diag(1 / np.sqrt(weights))
        g0_half = np.kron(root, I2)
        g0_inverse_half = np.kron(inverse_root, I2)
        gl_half = np.kron(root, spin_half)
        gl_inverse_half = np.kron(inverse_root, spin_inverse_half)
        return {
            "G0": g0_half @ g0_half,
            "GL": gl_half @ gl_half,
            "G0_half": g0_half,
            "G0_inverse_half": g0_inverse_half,
            "GL_half": gl_half,
            "GL_inverse_half": gl_inverse_half,
            "spin_metric": spin_metric,
        }
    def coordinate_kernel(self, weights: np.ndarray, canonical_unitary: np.ndarray):
        """K=G_L^-1/2 V G_0^1/2 for arbitrary channel unitary V."""
        data = self.metrics(weights)
        unitary = np.asarray(canonical_unitary, dtype=complex)
        size = data["G0"].shape[0]
        if unitary.shape != (size, size) or not np.isfinite(unitary).all():
            raise ValueError("finite canonical channel matrix with matching size required")
        return data["GL_inverse_half"] @ unitary @ data["G0_half"]

    def canonical_map(self, weights: np.ndarray, coordinate_kernel: np.ndarray):
        data = self.metrics(weights)
        kernel = np.asarray(coordinate_kernel, dtype=complex)
        return data["GL_half"] @ kernel @ data["G0_inverse_half"]

    def residuals(
        self,
        weights: np.ndarray,
        seed_covariance_blocks: np.ndarray,
        canonical_unitary: np.ndarray,
    ) -> dict[str, float | int]:
        """Evaluate interface variation, rank, CAR and frequency mixing."""
        data = self.metrics(weights)
        unitary = np.asarray(canonical_unitary, dtype=complex)
        kernel = self.coordinate_kernel(weights, unitary)
        reconstructed = self.canonical_map(weights, kernel)
        identity = np.eye(len(unitary), dtype=complex)
        interface = kernel.conj().T @ data["GL"] @ kernel - data["G0"]
        whitened = data["G0_inverse_half"] @ interface @ data["G0_inverse_half"]
        seed = full_covariance(seed_covariance_blocks)
        target = unitary @ seed @ unitary.conj().T
        seed_eigenvalues = np.linalg.eigvalsh(seed)
        target_eigenvalues = np.linalg.eigvalsh(target)
        return {
            "dimension": int(len(unitary)),
            "kernel_rank": int(np.linalg.matrix_rank(kernel)),
            "rank_defect": int(len(unitary) - np.linalg.matrix_rank(kernel)),
            "interface_variation_max_abs": float(np.max(np.abs(interface))),
            "interface_variation_whitened_norm": float(np.linalg.norm(whitened, 2)),
            "canonical_reconstruction_max_abs": float(np.max(np.abs(reconstructed - unitary))),
            "canonical_unitarity_max_abs": float(np.max(np.abs(unitary.conj().T @ unitary - identity))),
            "frequency_offblock_frobenius": frequency_offblock_norm(unitary),
            "target_hermiticity_max_abs": float(np.max(np.abs(target - target.conj().T))),
            "target_minimum_eigenvalue": float(target_eigenvalues.min()),
            "target_maximum_eigenvalue": float(target_eigenvalues.max()),
            "CAR_spectrum_transport_max_abs": float(np.max(np.abs(target_eigenvalues - seed_eigenvalues))),
        }
