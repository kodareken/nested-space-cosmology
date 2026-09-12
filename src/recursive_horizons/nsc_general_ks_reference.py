"""Fourth-order reference covariance on a supplied homogeneous KS history.

This owner generalizes the retained two-level Bloch recursion without
selecting a metric history.  The serialized positive-frequency quadrature is
completed by its axial-parity partner before the shift variation is formed.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi
from typing import Mapping, Sequence

import numpy as np

from .nsc_landau_cauchy_isometry import KSCauchyHistory


FIELDS = ("N", "beta", "q_ADM", "r")


def _first_derivative_matrix(time: np.ndarray, stencil: int) -> np.ndarray:
    """Local polynomial first derivative on an arbitrary ordered grid."""
    time = np.asarray(time, dtype=float)
    count = len(time)
    if stencil < 5 or stencil % 2 == 0 or count < stencil:
        raise ValueError("an odd derivative stencil >=5 within the history is required")
    matrix = np.zeros((count, count), dtype=float)
    half = stencil // 2
    for row in range(count):
        begin = min(max(0, row-half), count-stencil)
        index = np.arange(begin, begin+stencil)
        offset = time[index]-time[row]
        scale = float(np.max(abs(offset)))
        if not isfinite(scale) or scale <= 0:
            raise ValueError("distinct finite history times required")
        x = offset/scale
        vandermonde = np.vstack([x**power for power in range(stencil)])
        target = np.zeros(stencil)
        target[1] = 1.0/scale
        matrix[row, index] = np.linalg.solve(vandermonde, target)
    return matrix


@dataclass(frozen=True)
class ReferenceHistoryEvaluation:
    forces: Mapping[str, np.ndarray]
    residuals: Mapping[str, float]
    bloch_positive_k: np.ndarray
    bloch_negative_k: np.ndarray


@dataclass(frozen=True)
class GeneralKSFourthOrderReferenceHistory:
    """General-KS fourth-order reference for the serialized mode blocks.

    ``forces`` are the reference terms added to ``-Tr(C delta H)`` so that
    the canonical part becomes ``-Tr((C-C_ref) delta H)``.  They are action
    force densities per homogeneous axial coordinate, before any conversion
    into a claimed four-dimensional stress.
    """

    arrays: Mapping[str, np.ndarray]
    channels: Sequence[Mapping[str, object]]
    seed_a_parallel: float
    seed_radius: float = 1.0
    derivative_stencil: int = 9

    def __post_init__(self):
        if not isfinite(self.seed_a_parallel) or self.seed_a_parallel <= 0:
            raise ValueError("positive finite seed axial scale required")
        if not isfinite(self.seed_radius) or self.seed_radius <= 0:
            raise ValueError("positive finite seed radius required")
        required = {
            "hx_seed", "hy_seed", "hz_seed", "quadrature_weight",
            "sample_offsets",
        }
        if not required.issubset(self.arrays):
            raise ValueError("serialized mode-state arrays are incomplete")
        shape = np.asarray(self.arrays["hx_seed"]).shape
        if len(shape) != 1 or any(np.asarray(self.arrays[name]).shape != shape
                                  for name in required if name != "sample_offsets"):
            raise ValueError("mode-state sample arrays must share one flat shape")
        if not all(np.isfinite(np.asarray(self.arrays[name])).all()
                   for name in required):
            raise ValueError("finite serialized mode-state arrays required")

    def _mode_labels(self):
        count = len(np.asarray(self.arrays["hx_seed"]))
        mass = np.empty(count)
        angular = np.empty(count)
        multiplicity = np.empty(count)
        covered = np.zeros(count, dtype=bool)
        for channel in self.channels:
            begin = int(channel["sample_offset"])
            end = begin+int(channel["sample_count"])
            if not 0 <= begin < end <= count:
                raise ValueError("channel slice lies outside the state payload")
            mass[begin:end] = float(channel["compact_mass"])
            angular[begin:end] = float(channel["angular_eigenvalue"])
            multiplicity[begin:end] = (
                int(channel["copy_count"])*int(channel["degeneracy"])/pi
            )
            covered[begin:end] = True
        if not covered.all() or not np.isfinite(mass+angular+multiplicity).all():
            raise ValueError("channel metadata do not cover every sample")
        return mass, angular, multiplicity

    @staticmethod
    def _orders(h: np.ndarray, normal_derivative: np.ndarray):
        energy = np.linalg.norm(h, axis=2)
        if float(np.min(energy)) <= 1e-12:
            raise ValueError("the retained reference requires a nonzero mode gap")
        orders = [-h/energy[:, :, None]]
        recursion_residual = 0.0
        normalization_residual = 0.0
        for order in range(1, 5):
            derivative = np.einsum("ij,jfk->ifk", normal_derivative, orders[-1])
            perpendicular = -np.cross(h, derivative)/(2*energy[:, :, None]**2)
            normalization = np.zeros_like(energy)
            for left in range(1, order):
                normalization += np.einsum(
                    "ifk,ifk->if", orders[left], orders[order-left]
                )
            value = perpendicular-0.5*normalization[:, :, None]*orders[0]
            recurrence = value-perpendicular+0.5*normalization[:, :, None]*orders[0]
            recursion_residual = max(recursion_residual, float(np.max(abs(recurrence))))
            norm_order = (
                2*np.einsum("ifk,ifk->if", orders[0], value)+normalization
            )
            normalization_residual = max(
                normalization_residual, float(np.max(abs(norm_order)))
            )
            orders.append(value)
        return orders, energy, recursion_residual, normalization_residual

    def evaluate(self, history: KSCauchyHistory) -> ReferenceHistoryEvaluation:
        history.validate()
        time = np.asarray(history.time, dtype=float)
        lapse = np.asarray(history.lapse, dtype=float)
        a = np.asarray(history.a_parallel, dtype=float)
        radius = np.asarray(history.radius, dtype=float)
        derivative = _first_derivative_matrix(time, self.derivative_stencil)
        normal_derivative = derivative/lapse[:, None]

        mass, angular, multiplicity = self._mode_labels()
        stored_k = np.asarray(self.arrays["hz_seed"], dtype=float)*self.seed_a_parallel
        k = abs(stored_k)

        def axes(signed_k):
            value = np.empty((len(time), len(k), 3), dtype=float)
            value[:, :, 0] = -mass[None, :]
            value[:, :, 1] = angular[None, :]/radius[:, None]
            value[:, :, 2] = signed_k[None, :]/a[:, None]
            return value

        h_plus, h_minus = axes(k), axes(-k)
        plus, gap_plus, rr_plus, nr_plus = self._orders(h_plus, normal_derivative)
        minus, gap_minus, rr_minus, nr_minus = self._orders(h_minus, normal_derivative)
        bloch_plus = np.sum(plus, axis=0)
        bloch_minus = np.sum(minus, axis=0)

        # The positive-frequency quadrature represents the full +/-k integral
        # as the average of the two parity partners with dk/pi.
        projection_n = 0.5*(
            np.einsum("ifk,ifk->if", bloch_plus, h_plus)
            + np.einsum("ifk,ifk->if", bloch_minus, h_minus)
        )
        projection_beta = np.zeros_like(h_plus[:, :, 0])
        projection_q = 0.5*(
            (-lapse[:, None]*k[None, :]/a[:, None]**2)*bloch_plus[:, :, 2]
            +(lapse[:, None]*k[None, :]/a[:, None]**2)*bloch_minus[:, :, 2]
        )
        projection_r = -0.5*lapse[:, None]*angular[None, :]/radius[:, None]**2*(
            bloch_plus[:, :, 1]+bloch_minus[:, :, 1]
        )
        weighted = multiplicity*np.asarray(self.arrays["quadrature_weight"], dtype=float)
        forces = {
            "N": np.einsum("f,if->i", weighted, projection_n),
            "beta": np.einsum("f,if->i", weighted, projection_beta),
            "q_ADM": np.einsum("f,if->i", weighted, projection_q),
            "r": np.einsum("f,if->i", weighted, projection_r),
        }

        original_axes = np.stack([
            -mass,
            angular/self.seed_radius,
            stored_k/self.seed_a_parallel,
        ], axis=1)
        stored_axes = np.stack([
            np.asarray(self.arrays["hx_seed"]),
            np.asarray(self.arrays["hy_seed"]),
            np.asarray(self.arrays["hz_seed"]),
        ], axis=1)
        seed_axis_residual = float(np.max(abs(original_axes-stored_axes)))
        trace_residual = 0.0  # (I+b.sigma)/2 has trace one identically.
        parity_momentum = float(np.max(abs(projection_beta)))
        residuals = {
            "seed_hamiltonian_axis_residual": seed_axis_residual,
            "maximum_bloch_recursion_residual": max(rr_plus, rr_minus),
            "maximum_order_normalization_residual": max(nr_plus, nr_minus),
            "maximum_projector_trace_residual": trace_residual,
            "maximum_parity_completed_beta_projection": parity_momentum,
            "minimum_reference_gap": float(min(np.min(gap_plus), np.min(gap_minus))),
        }
        return ReferenceHistoryEvaluation(
            forces=forces,
            residuals=residuals,
            bloch_positive_k=bloch_plus,
            bloch_negative_k=bloch_minus,
        )
