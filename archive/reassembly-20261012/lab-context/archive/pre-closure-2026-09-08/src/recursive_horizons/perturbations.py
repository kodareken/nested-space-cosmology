"""Deterministic linear-mode transfer through a closed de Sitter core.

This module implements only the selected homogeneous core benchmark.  Its
transfer matrices evolve classical canonical phase-space data ``(y, y')``.
They are *not* a calculation of the Mukhanov--Sasaki curvature spectrum,
Bogoliubov particle production, ``n_s``, ``A_s``, reheating, or a CMB
prediction.  Those require endpoint matching, a quantum state, and a complete
matter/transition model that the repository does not yet provide.

For the closed de Sitter conformal-time patch,
``a(eta) = L sec(eta)``, the two benchmark mode equations are

* spectator/order-parameter scalar (Phi' = 0):
  ``u_n'' + [(n + 1)^2 + ((m L)^2 - 2) sec^2(eta)] u_n = 0``;
* tensor, in the convention stated in the accompanying documentation:
  ``mu_n'' + [n^2 - 2 sec^2(eta)] mu_n = 0``.

The scalar is gauge invariant precisely because the homogeneous background
field is constant.  It must not be relabelled as the curvature perturbation.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, isfinite
from typing import Callable


_HALF_PI = 1.5707963267948966


def _finite(name: str, value: float) -> float:
    value = float(value)
    if not isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _positive_int(name: str, value: int, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer greater than or equal to {minimum}")
    return value


def _core_eta(eta: float) -> float:
    eta = _finite("eta", eta)
    if abs(eta) >= _HALF_PI:
        raise ValueError("eta must satisfy abs(eta) < pi/2 inside the closed de Sitter patch")
    return eta


def closed_desitter_scale_factor(eta: float, curvature_radius: float = 1.0) -> float:
    """Return ``a(eta) = L sec(eta)`` for positive curvature radius ``L``.

    ``eta`` is dimensionless conformal time in the open interval ``(-pi/2, pi/2)``.
    """

    eta = _core_eta(eta)
    curvature_radius = _finite("curvature_radius", curvature_radius)
    if curvature_radius <= 0.0:
        raise ValueError("curvature_radius must be finite and positive")
    return curvature_radius / cos(eta)


def scale_factor_potential(eta: float) -> float:
    """Return ``a''/a = 2 sec^2(eta) - 1`` for ``a = L sec(eta)``."""

    eta = _core_eta(eta)
    sec = 1.0 / cos(eta)
    return 2.0 * sec * sec - 1.0


def spectator_frequency_squared(n: int, eta: float, mass_radius: float = 0.0) -> float:
    """Return the scalar benchmark frequency squared.

    ``mass_radius`` is the dimensionless product ``m L``.  Scalar harmonics on
    ``S^3`` use ``n >= 0``.  This is a spectator/order-parameter fluctuation;
    it is gauge invariant for the stated constant-background-field case.
    """

    n = _positive_int("n", n, 0)
    eta = _core_eta(eta)
    mass_radius = _finite("mass_radius", mass_radius)
    sec = 1.0 / cos(eta)
    return (n + 1) ** 2 + (mass_radius**2 - 2.0) * sec * sec


def tensor_frequency_squared(n: int, eta: float) -> float:
    """Return the tensor benchmark frequency squared for convention ``n >= 3``."""

    n = _positive_int("n", n, 3)
    eta = _core_eta(eta)
    sec = 1.0 / cos(eta)
    return n**2 - 2.0 * sec * sec


@dataclass(frozen=True, slots=True)
class TransferMatrix:
    """A real 2-by-2 canonical transfer matrix for ``(y, y')`` data."""

    m11: float
    m12: float
    m21: float
    m22: float

    @property
    def determinant(self) -> float:
        """Return the Wronskian-preserving determinant (one in exact evolution)."""

        return self.m11 * self.m22 - self.m12 * self.m21

    def apply(self, value: float, derivative: float) -> tuple[float, float]:
        """Apply the matrix to phase-space data ``(y, y')``."""

        return (
            self.m11 * value + self.m12 * derivative,
            self.m21 * value + self.m22 * derivative,
        )

    def max_abs_difference(self, other: "TransferMatrix") -> float:
        """Return the largest componentwise absolute difference."""

        return max(
            abs(self.m11 - other.m11),
            abs(self.m12 - other.m12),
            abs(self.m21 - other.m21),
            abs(self.m22 - other.m22),
        )


@dataclass(frozen=True, slots=True)
class TransferConvergence:
    """One refinement record for a deterministic transfer calculation."""

    steps: int
    matrix: TransferMatrix
    determinant_error: float
    difference_from_previous: float | None


def _rk4_step(
    eta: float,
    step: float,
    y: tuple[float, float, float, float],
    omega_squared: Callable[[float], float],
) -> tuple[float, float, float, float]:
    """Advance two independent solutions of ``y'' + omega^2 y = 0`` once."""

    def rhs(
        time: float, state: tuple[float, float, float, float]
    ) -> tuple[float, float, float, float]:
        w2 = omega_squared(time)
        return (state[1], -w2 * state[0], state[3], -w2 * state[2])

    def add(
        state: tuple[float, float, float, float],
        factor: float,
        increment: tuple[float, float, float, float],
    ) -> tuple[float, float, float, float]:
        return tuple(
            value + factor * delta
            for value, delta in zip(state, increment, strict=True)
        )  # type: ignore[return-value]

    k1 = rhs(eta, y)
    k2 = rhs(eta + step / 2.0, add(y, step / 2.0, k1))
    k3 = rhs(eta + step / 2.0, add(y, step / 2.0, k2))
    k4 = rhs(eta + step, add(y, step, k3))
    return tuple(
        value + step * (a + 2.0 * b + 2.0 * c + d) / 6.0
        for value, a, b, c, d in zip(y, k1, k2, k3, k4, strict=True)
    )  # type: ignore[return-value]


def core_transfer_matrix(
    *,
    sector: str,
    n: int,
    eta_core: float,
    steps: int = 4_096,
    mass_radius: float = 0.0,
) -> TransferMatrix:
    """Integrate the core transfer from ``-eta_core`` to ``+eta_core`` by RK4.

    The output maps input ``(y, y')`` to final ``(y, y')``.  It is a real,
    classical symplectic calculation.  It deliberately makes no particle or
    power-spectrum claim.  ``sector`` is either ``"scalar"`` or ``"tensor"``.
    """

    eta_core = _finite("eta_core", eta_core)
    if eta_core < 0.0 or eta_core >= _HALF_PI:
        raise ValueError("eta_core must satisfy 0 <= eta_core < pi/2")
    steps = _positive_int("steps", steps, 1)
    mass_radius = _finite("mass_radius", mass_radius)

    if sector == "scalar":
        _positive_int("n", n, 0)
        omega_squared = lambda eta: spectator_frequency_squared(n, eta, mass_radius)
    elif sector == "tensor":
        _positive_int("n", n, 3)
        omega_squared = lambda eta: tensor_frequency_squared(n, eta)
    else:
        raise ValueError('sector must be "scalar" or "tensor"')

    if eta_core == 0.0:
        return TransferMatrix(1.0, 0.0, 0.0, 1.0)

    eta = -eta_core
    step = 2.0 * eta_core / steps
    # Column one starts at (1, 0); column two starts at (0, 1).
    state = (1.0, 0.0, 0.0, 1.0)
    for _ in range(steps):
        state = _rk4_step(eta, step, state, omega_squared)
        eta += step
    return TransferMatrix(state[0], state[2], state[1], state[3])


def transfer_convergence(
    *,
    sector: str,
    n: int,
    eta_core: float,
    initial_steps: int = 256,
    refinements: int = 4,
    mass_radius: float = 0.0,
) -> tuple[TransferConvergence, ...]:
    """Return deterministic step-refinement records for a core transfer.

    Each level doubles the number of RK4 steps.  The first record has no prior
    difference; later records report the max component difference from the
    preceding level.  Callers should examine both that difference and the
    Wronskian/determinant error rather than trusting a single grid.
    """

    initial_steps = _positive_int("initial_steps", initial_steps, 1)
    refinements = _positive_int("refinements", refinements, 1)
    previous: TransferMatrix | None = None
    records: list[TransferConvergence] = []
    for level in range(refinements):
        steps = initial_steps * (2**level)
        matrix = core_transfer_matrix(
            sector=sector,
            n=n,
            eta_core=eta_core,
            steps=steps,
            mass_radius=mass_radius,
        )
        records.append(
            TransferConvergence(
                steps=steps,
                matrix=matrix,
                determinant_error=abs(matrix.determinant - 1.0),
                difference_from_previous=(
                    None if previous is None else matrix.max_abs_difference(previous)
                ),
            )
        )
        previous = matrix
    return tuple(records)
