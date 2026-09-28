"""Canonical closed-time-path owner for a finite parent--child Dirac system.

The Lorentzian Hamiltonian and its Gaussian covariance carry all causal and
state-dependent response.  State-independent induced forces are supplied by
the same spectral coefficient account and are added exactly once.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Mapping

import numpy as np
from scipy.linalg import expm

from .nsc_boundary_state import GaussianBoundaryState
from .nsc_influence import _covariance, influence


FIELDS = ("N", "beta", "q", "r")


def _history(values, slices, size, name):
    array = np.asarray(values, dtype=complex)
    if array.shape != (slices, size, size) or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a finite square operator history")
    if not all(np.allclose(item, item.conj().T, atol=1e-11, rtol=0) for item in array):
        raise ValueError(f"{name} history must be Hermitian")
    return array


def _vertex_history(values, slices, size, name):
    array = np.asarray(values, dtype=complex)
    if array.shape == (size, size):
        array = np.repeat(array[None, :, :], slices, axis=0)
    if array.shape != (slices, size, size) or not np.isfinite(array).all():
        raise ValueError(f"{name} vertex has the wrong history shape")
    if not all(np.allclose(item, item.conj().T, atol=1e-11, rtol=0) for item in array):
        raise ValueError(f"{name} vertex must be Hermitian")
    return array


@dataclass(frozen=True)
class MetricHistory:
    times: np.ndarray
    parent_hamiltonians: np.ndarray
    child_hamiltonians: np.ndarray
    metric_vertices: Mapping[str, np.ndarray]

    def dimensions(self):
        times = np.asarray(self.times, dtype=float)
        if times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all():
            raise ValueError("at least two finite history times required")
        if not np.all(np.diff(times) > 0):
            raise ValueError("history times must increase strictly")
        parent = np.asarray(self.parent_hamiltonians)
        child = np.asarray(self.child_hamiltonians)
        if parent.ndim != 3 or child.ndim != 3:
            raise ValueError("parent and child histories must have three axes")
        np_, nc = parent.shape[1], child.shape[1]
        _history(parent, len(times), np_, "parent")
        _history(child, len(times), nc, "child")
        if set(self.metric_vertices) != set(FIELDS):
            raise ValueError("metric vertices must contain N, beta, q and r")
        return times, np_, nc


@dataclass(frozen=True)
class GaugeHistory:
    parent_additions: np.ndarray
    child_additions: np.ndarray


@dataclass(frozen=True)
class LinkHistory:
    values: np.ndarray


@dataclass(frozen=True)
class SpectralInducedSource:
    """State-independent source already derived from the common spectrum."""

    forces: Mapping[str, float]
    coefficient_owner: str
    unresolved_terms: tuple[str, ...] = ()

    def __post_init__(self):
        if set(self.forces) != set(FIELDS):
            raise ValueError("one induced force for every metric direction required")
        if not all(np.isfinite(tuple(self.forces.values()))):
            raise ValueError("finite induced forces required")
        if not self.coefficient_owner:
            raise ValueError("spectral coefficient owner required")


@dataclass(frozen=True)
class BoundaryDomain:
    name: str
    parent_dimension: int
    child_dimension: int
    spin_frame: str
    units: str
    stress_factors: Mapping[str, tuple[str, float]] = field(default_factory=dict)
    power_weights: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self):
        if not self.name or not self.spin_frame or not self.units:
            raise ValueError("named domain, spin frame and units required")
        if self.parent_dimension < 1 or self.child_dimension < 1:
            raise ValueError("positive room dimensions required")
        for output, (source, factor) in self.stress_factors.items():
            if source not in FIELDS or not output or not np.isfinite(factor):
                raise ValueError("invalid stress projection")
        if any(source not in FIELDS or not np.isfinite(weight)
               for source, weight in self.power_weights.items()):
            raise ValueError("invalid boundary-power projection")


@dataclass(frozen=True)
class CausalSource:
    Gamma_CTP: complex
    forces: Mapping[str, float]
    matter_forces: Mapping[str, float]
    induced_forces: Mapping[str, float]
    stress: Mapping[str, float]
    boundary_power: float
    retarded_response: Mapping[str, tuple[float, ...]]
    noise_kernel: Mapping[str, tuple[tuple[float, ...], ...]]
    ward_residuals: Mapping[str, float]
    unresolved_terms: tuple[str, ...]
    final_covariance: np.ndarray

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["Gamma_CTP"] = [self.Gamma_CTP.real, self.Gamma_CTP.imag]
        value["retarded_response"] = {
            name: list(rows) for name, rows in self.retarded_response.items()
        }
        value["noise_kernel"] = {
            name: [list(row) for row in rows]
            for name, rows in self.noise_kernel.items()
        }
        value["unresolved_terms"] = list(self.unresolved_terms)
        value["final_covariance"] = [
            [[float(z.real), float(z.imag)] for z in row]
            for row in self.final_covariance
        ]
        return value


class CausalCommonFunctional:
    """Evaluate the normalized physical branch and its causal bilinear kernels."""

    @staticmethod
    def _assemble(metric_history, gauge_history, link_history):
        times, np_, nc = metric_history.dimensions()
        slices = len(times)
        parent = _history(metric_history.parent_hamiltonians, slices, np_, "parent")
        child = _history(metric_history.child_hamiltonians, slices, nc, "child")
        gp = _history(gauge_history.parent_additions, slices, np_, "parent gauge")
        gc = _history(gauge_history.child_additions, slices, nc, "child gauge")
        links = np.asarray(link_history.values, dtype=complex)
        if links.shape != (slices, np_, nc) or not np.isfinite(links).all():
            raise ValueError("finite link history must match both rooms")
        full = []
        for hp, hc, ap, ac, link in zip(parent, child, gp, gc, links):
            full.append(np.block([[hp + ap, link], [link.conj().T, hc + ac]]))
        return times, np_, nc, np.asarray(full)

    @staticmethod
    def _evolve(times, hamiltonians, covariance):
        c = np.array(_covariance(covariance), copy=True)
        identity = np.eye(len(c))
        cumulative = identity.astype(complex)
        histories = [c.copy()]
        unitaries = [cumulative.copy()]
        work = 0.0
        energy_initial = float(np.trace(c @ hamiltonians[0]).real)
        maximum_unitarity_error = 0.0
        for index, dt in enumerate(np.diff(times)):
            step = expm(-1j * dt * hamiltonians[index])
            maximum_unitarity_error = max(
                maximum_unitarity_error,
                float(np.linalg.norm(step.conj().T @ step - identity)),
            )
            cumulative = step @ cumulative
            c = step @ c @ step.conj().T
            work += float(np.trace(c @ (hamiltonians[index + 1] - hamiltonians[index])).real)
            histories.append(c.copy())
            unitaries.append(cumulative.copy())
        energy_final = float(np.trace(c @ hamiltonians[-1]).real)
        return {
            "covariances": histories,
            "unitaries": unitaries,
            "final": c,
            "work": work,
            "energy_initial": energy_initial,
            "energy_final": energy_final,
            "unitarity_error": maximum_unitarity_error,
        }

    def evaluate(
        self,
        metric_history: MetricHistory,
        gauge_history: GaugeHistory,
        link_history: LinkHistory,
        initial_covariance: np.ndarray,
        spectral_coefficients: SpectralInducedSource,
        boundary_domain: BoundaryDomain,
    ) -> CausalSource:
        times, np_, nc, hamiltonians = self._assemble(
            metric_history, gauge_history, link_history
        )
        if (np_, nc) != (boundary_domain.parent_dimension, boundary_domain.child_dimension):
            raise ValueError("boundary domain dimensions differ from the histories")
        partition = np.r_[np.ones(np_, dtype=bool), np.zeros(nc, dtype=bool)]
        boundary_state = GaussianBoundaryState(initial_covariance, partition)
        evolution = self._evolve(times, hamiltonians, boundary_state.covariance)
        covariance = evolution["final"]
        size = np_ + nc
        vertices = {
            name: _vertex_history(values, len(times), size, name)
            for name, values in metric_history.metric_vertices.items()
        }
        matter = {
            name: float(-np.trace(covariance @ values[-1]).real)
            for name, values in vertices.items()
        }
        induced = {name: float(spectral_coefficients.forces[name]) for name in FIELDS}
        forces = {name: matter[name] + induced[name] for name in FIELDS}
        stress = {
            output: float(factor * forces[source])
            for output, (source, factor) in boundary_domain.stress_factors.items()
        }
        boundary_power = float(sum(
            boundary_domain.power_weights.get(name, 0.0) * forces[name] for name in FIELDS
        ))

        retarded, noise = {}, {}
        c0 = _covariance(initial_covariance)
        identity = np.eye(size)
        for name, values in vertices.items():
            heisenberg = [u.conj().T @ vertex @ u
                          for u, vertex in zip(evolution["unitaries"], values)]
            source = heisenberg[0]
            retarded[name] = tuple(
                0.0 if index == 0 else float(
                    (-1j * np.trace(c0 @ (operator @ source - source @ operator))).real
                )
                for index, operator in enumerate(heisenberg)
            )
            matrix = np.array([
                [float(np.trace(c0 @ left @ (identity - c0) @ right).real)
                 for right in heisenberg]
                for left in heisenberg
            ])
            noise[name] = tuple(tuple(float(value) for value in row) for row in matrix)

        equal = influence(c0, evolution["unitaries"][-1], evolution["unitaries"][-1])
        reduced_equal = boundary_state.overlap_schur(
            evolution["unitaries"][-1], evolution["unitaries"][-1]
        )
        gamma = equal["principal_action"]
        if gamma is None:
            raise ArithmeticError("equal-history influence became singular")
        ward = {
            "equal_history_action": float(abs(gamma)),
            "boundary_schur_equal_history_action": float(
                abs(reduced_equal["principal_action"] or 0j)
            ),
            "energy_work_balance": float(
                evolution["energy_final"] - evolution["energy_initial"] - evolution["work"]
            ),
            "covariance_trace": float(np.trace(covariance).real - np.trace(c0).real),
            "unitarity": evolution["unitarity_error"],
            "covariance_hermiticity": float(np.linalg.norm(covariance - covariance.conj().T)),
        }
        return CausalSource(
            Gamma_CTP=complex(gamma),
            forces=forces,
            matter_forces=matter,
            induced_forces=induced,
            stress=stress,
            boundary_power=boundary_power,
            retarded_response=retarded,
            noise_kernel=noise,
            ward_residuals=ward,
            unresolved_terms=tuple(spectral_coefficients.unresolved_terms),
            final_covariance=covariance,
        )
