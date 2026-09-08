"""Finite regulated Dirac calculus and normalized recursion.

This is an explicit finite-operator realization, not a completed covariant
quantum action. In particular a spatial Hamiltonian heat trace is not silently
identified with a spacetime determinant. All matrices act in orthonormal bases.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import log, pi, sqrt
import numpy as np
from scipy.linalg import eigh, expm
from scipy.special import erfc, exp1


def hermitian(matrix):
    value = np.array(matrix, dtype=complex, copy=True)
    if value.ndim != 2 or value.shape[0] != value.shape[1] or value.shape[0] == 0:
        raise ValueError("operator must be a nonempty square matrix")
    if not np.isfinite(value).all() or not np.allclose(value, value.conj().T, atol=1e-12, rtol=1e-12):
        raise ValueError("operator must be finite and Hermitian in the declared measure")
    value.setflags(write=False)
    return value


@dataclass(frozen=True)
class OperatorConventions:
    cutoff: float = 1.0
    normalization: float = 1.0
    throat_length: float = 1.0
    child_ratio: float = 1.0
    spatial_dimension: int = 4

    def __post_init__(self):
        for name in ("cutoff", "normalization", "throat_length", "child_ratio"):
            value = getattr(self, name)
            if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if isinstance(self.spatial_dimension, bool) or self.spatial_dimension not in (1, 2, 3, 4):
            raise ValueError("unsupported spatial dimension")

    @property
    def zeta(self):
        return (self.cutoff * self.throat_length) ** 2


@dataclass(frozen=True)
class ClosureSolution:
    """Report mathematical closure separately from unsolved physical equations."""
    scope: str
    residual: float
    unresolved_equations: tuple[str, ...]

    @property
    def physical_stationarity_solved(self):
        # Numerical algebraic closure cannot promote a field-equation claim.
        return self.scope == "full_covariant_stationarity" and not self.unresolved_equations


class RegulatedOperator:
    """Proper-time regulated modulus with an explicit finite normalization.

    Gamma = 1/2 Tr E1(D²/Lambda²) + N[log(M/Lambda)+gamma_E/2].
    At fixed nonzero finite eigenvalues, Lambda -> infinity yields
    -log|det(D/M)|. The last bracket is a declared finite-matrix normalization;
    it is not asserted to be the continuum gravitational counterterm measure.
    """

    def __init__(self, matrix, conventions=None):
        self.matrix = hermitian(matrix)
        self.conventions = conventions or OperatorConventions()
        self.values, self.vectors = eigh(self.matrix)
        if np.min(np.abs(self.values)) < 1e-12 * max(1., np.max(np.abs(self.values))):
            raise ValueError("zero or unresolved near-zero mode: specify reduced determinant and zero-mode measure")

    def action(self):
        c = self.conventions
        return float(.5 * np.sum(exp1((self.values / c.cutoff) ** 2))
                     + len(self.values) * (log(c.normalization / c.cutoff) + np.euler_gamma / 2))

    def negative_logdet_phase(self):
        """Unwrapped eigenvalue-log branch; changes require a zero crossing."""
        return -pi * int(np.sum(self.values < 0))

    def variation(self, first, second=None):
        """First and second total variation for D(t)=D+t D1+t² D2/2.

        The divided-difference Hessian retains noncommuting perturbations and
        coincident eigenvalues; it is not eigenvalue differentiation with fixed
        eigenvectors. A mass, metric, or boundary variation enters through its
        actual D1 and D2 in the common orthonormal measure.
        """
        first = hermitian(first)
        second = np.zeros_like(first) if second is None else hermitian(second)
        if first.shape != self.matrix.shape or second.shape != first.shape:
            raise ValueError("variation dimension differs from operator")
        d = self.values
        cutoff = self.conventions.cutoff
        exponential = np.exp(-(d / cutoff) ** 2)
        fp = -exponential / d
        fpp = exponential * (2 / cutoff**2 + 1 / d**2)
        v = self.vectors.conj().T @ first @ self.vectors
        a = self.vectors.conj().T @ second @ self.vectors
        divided = np.empty((len(d), len(d)))
        for i in range(len(d)):
            for j in range(len(d)):
                if abs(d[i] - d[j]) <= 1e-10 * max(1., abs(d[i]), abs(d[j])):
                    divided[i, j] = (fpp[i] + fpp[j]) / 2
                else:
                    divided[i, j] = (fp[i] - fp[j]) / (d[i] - d[j])
        gradient = float(np.real(np.dot(fp, np.diag(v))))
        hessian = float(np.real(np.sum(divided * np.abs(v)**2) + np.dot(fp, np.diag(a))))
        return gradient, hessian

    def scale_orbit(self, generator, t):
        sigma = hermitian(generator)
        if sigma.shape != self.matrix.shape or not np.isfinite(t):
            raise ValueError("invalid scale orbit")
        v = expm(-t * sigma / 2)
        return RegulatedOperator(v @ self.matrix @ v, self.conventions)

    def scale_defect(self, generator, t):
        """Regulator defect relative to exact finite Grassmann determinant.

        C(t)=Gamma_reg(D_t)-Gamma_reg(D)-t Tr Sigma.
        Removing C restores the raw finite transformation law on THIS orbit.
        This algebra does not declare the relative-sheet field a gauge degree
        of freedom or supply its continuum physical compensator dynamics.
        """
        sigma = hermitian(generator)
        return self.scale_orbit(sigma, t).action() - self.action() - t * float(np.trace(sigma).real)

    def scale_defect_derivative(self, generator, t):
        sigma = hermitian(generator)
        transformed = self.scale_orbit(sigma, t)
        heat = expm(-transformed.matrix @ transformed.matrix / self.conventions.cutoff**2)
        return float(np.trace(sigma @ (heat - np.eye(len(heat)))).real)

    def cutoff_derivative(self):
        """dGamma/dlog Lambda at fixed dimensional D and normalization M."""
        return float(np.sum(np.exp(-(self.values / self.conventions.cutoff)**2)) - len(self.values))

    def ultrastatic_action_per_time(self):
        """Unnormalized proper-time modulus after integrating one frequency.

        Requires a separate ultrastatic Euclidean product operator whose square
        is -partial_tau²+D_spatial². This assumption fails for a generic shifted
        horizon-penetrating Hamiltonian; no continuation is inferred here.
        """
        energies = np.abs(self.values)
        cutoff = self.conventions.cutoff
        return float(np.sum(cutoff * np.exp(-(energies / cutoff)**2) / (2 * sqrt(pi))
                            - energies * erfc(energies / cutoff) / 2))


def recursive_response(local, link, energy: complex, omega: float, depth: int):
    """Dimensionless parent response for a finite chain in one energy frame.

    H_n=Omega^n H_0, B_n=Omega^n b, Lambda_parent=1. Inverse
    response recurses with 1/Omega, not an independently fitted link weight.
    """
    h = hermitian(local)
    b = np.asarray(link, complex)
    if b.shape != h.shape or not np.isfinite(b).all():
        raise ValueError("invalid link")
    if isinstance(depth, bool) or not isinstance(depth, int) or depth < 1:
        raise ValueError("depth must be a positive integer")
    if not np.isfinite(omega) or omega <= 0 or not np.isfinite(energy) or energy.imag <= 0:
        raise ValueError("use positive Omega and upper-half-plane energy")
    gamma = energy / omega**(depth - 1) * np.eye(len(h)) - h
    for n in range(depth - 2, -1, -1):
        gamma = energy / omega**n * np.eye(len(h)) - h - b @ np.linalg.solve(gamma, b.conj().T) / omega
    return gamma


def direct_chain(local, link, energy, omega, depth):
    h = hermitian(local)
    b = np.asarray(link, complex)
    size = len(h)
    chain = np.zeros((depth * size, depth * size), complex)
    for n in range(depth):
        sl = slice(n * size, (n + 1) * size)
        chain[sl, sl] = omega**n * h
        if n + 1 < depth:
            nxt = slice((n + 1) * size, (n + 2) * size)
            chain[sl, nxt] = omega**n * b
            chain[nxt, sl] = omega**n * b.conj().T
    full_inverse = np.linalg.inv(energy * np.eye(len(chain)) - chain)
    return np.linalg.inv(full_inverse[:size, :size])


def pg_principal_witness(rho=0.):
    """Nonelliptic coordinate-time Hamiltonian symbol in the trapped region.

    A=1+3rho+3(1+rho²)(atan(rho)-pi/2), beta=sqrt(1-A).
    Eigenvalues of -beta*k_r I + alpha·k are -beta*k_r +/- |k|.
    Lapse and inverse spatial tetrad cancel their common warp at this order.
    Nonellipticity of this spatial operator does NOT imply ill-posed covariant
    Lorentzian Dirac evolution.
    """
    a = float(1 + 3 * rho + 3 * (1 + rho**2) * (np.arctan(rho) - pi / 2))
    beta = sqrt(1 - a)
    if beta < 1:
        raise ValueError("witness requires horizon or trapped region")
    radial = 1.
    tangent = sqrt(max(0., beta**2 - 1))
    eigenvalues = [-beta - sqrt(radial**2 + tangent**2), -beta + sqrt(radial**2 + tangent**2)]
    return {"rho": rho, "A": a, "beta": beta, "covector": [radial, tangent],
            "principal_eigenvalues": eigenvalues,
            "spatial_Hamiltonian_elliptic": False,
            "covariant_Lorentzian_Dirac_ill_posed_inferred": False}
