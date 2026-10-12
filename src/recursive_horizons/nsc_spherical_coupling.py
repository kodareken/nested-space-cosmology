"""Provisional finite-source coupling of a rank-6 canonical Gaussian to the Weyl chart.

Gamma_one is the induced spherical chart counted once plus the canonical
Gaussian mean field of one finite covariance C0. There is no sea subtraction,
no extra heat-kernel action, and no Gamma_rest. The historical vacuum-matched
quadratic action is not a prerequisite and is not called. The construction is a finite mean-field Cauchy system, not an
absolute continuum vacuum and not a regeneration result.

The discrete generator is the summation-by-parts energy. Geometry uses a
periodic real antisymmetric derivative. Fermions use the existing
antiperiodic momentum. Column forces are the nodal derivatives of
M Tr(C H), M=4 kappa once, and are divided by dx only when a constraint
density or a symplectic momentum equation needs that density.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sympy as sp

from .nsc_conformal_adm_source import conformal_source, sector_multiplicity
from .nsc_covariant_operator import CovariantStaticMetric
from .nsc_nested_qualities import finite_window
from .nsc_spherical_feedback_action import (
    alpha_of,
    feedback_F,
    feedback_V,
    feedback_Z,
    partial_F,
    velocity_rhs,
)

SCHEMA = "NSC-SPHERICAL-COUPLING-CONTROL-v1"
STATUS = "PROVISIONAL"
CAUCHY_OWNED = True
SOURCE_COUPLING_OWNED = True
GLOBAL_REGENERATION = False
ABSOLUTE_VACUUM_CLAIM = False
VACUUM_MATCHED_PREREQUISITE = False
SEA_SUBTRACTED = False
GAMMA_REST_INCLUDED = False
OLD_CIRCULATION_TRANSFERRED = False

# Predeclared before the bounded run. Not retuned to force a pass.
TOL_INITIAL_CONSTRAINT = 1e-8
TOL_FD_RELATIVE = 1e-6
TOL_ENERGY_RELATIVE = 1e-6
TOL_CONSTRAINT_VALIDATION = 1e-3
TOL_RESOLVED_CONSTRAINT = 1e-2
TOL_NEWTON = 1e-10

OCCUPATIONS = np.array([0.75, 0.75, 0.5, 0.5, 0.25, 0.25], dtype=float)
PACKET_COUNT = 3
ELL = 1.0
CARRIER_K = 1.0
OMEGA = 1.5
PERIOD = 8.0
KAPPA = 1
CALIBRATION = {
    "Mweight_executed": 1.533809137604,
    "a0": 0.977957402407,
    "beta0": 0.325985800803,
    "b0": 0.245785078900,
    "phase": -0.842669616901,
    "kappa": KAPPA,
    "Omega": OMEGA,
    "k": CARRIER_K,
    "ell": ELL,
    "length": PERIOD,
}

_LOCKED_NAMES = {
    "A": "A",
    "C_W": "C_Weyl",
    "C_F": "C_gauge",
    "C_E": "C_Euler",
    "C_box": "C_boxR",
    "flux": "magnetic_flux",
    "V_rel": "V_full",
}
_LOCKED_EXPECTED = {
    "A": 0.04501936182826115,
    "C_W": -0.0008441130342798952,
    "C_F": 0.01125484045706527,
    "C_E": 0.0005158468542821582,
    "C_box": -0.0005627420228532635,
    "flux": 4.0,
    "V_rel": 0.0,
}

_LAB_ROOT = Path(__file__).resolve().parents[2]
LOCKED_RECORD = _LAB_ROOT / "results" / "development" / "nsc-subgap-history-response.json"
CONTROL_RECORD = _LAB_ROOT / "results" / "development" / "nsc-spherical-coupling-control-v1.json"


class PositiveChartExit(RuntimeError):
    """Integration left the positive chart. The failed step is not a result."""

    def __init__(self, reason, time_value, last_state):
        super().__init__(reason)
        self.reason = reason
        self.time_value = float(time_value)
        self.last_state = last_state


def locked_coefficients(path=None):
    """Fixed spectral coefficients. V_rel is recorded and is not added to V."""
    record = Path(path) if path is not None else LOCKED_RECORD
    payload = json.loads(record.read_text())
    source = payload["locked_inputs"]
    coefficients = {name: float(source[key]) for name, key in _LOCKED_NAMES.items()}
    for name, expected in _LOCKED_EXPECTED.items():
        if coefficients[name] != expected:
            raise ValueError(f"locked coefficient {name} drifted from the audited literal")
    if coefficients["C_W"] == 0 or coefficients["A"] == 0:
        raise ValueError("positive Weyl chart rejects vanishing A or C_W")
    if coefficients["V_rel"] != 0:
        raise ValueError("this coupling does not carry a nonzero V_rel into V")
    return coefficients


def magnetic_radius_square(coefficients):
    """C_F flux^2 / (4 A). The audited coefficients give 1."""
    return coefficients["C_F"] * coefficients["flux"] ** 2 / (4 * coefficients["A"])


def periodic_derivative(points, length):
    """Real antisymmetric periodic spectral derivative. Nyquist symbol is 0."""
    points = int(points)
    if points < 4 or points % 2:
        raise ValueError("periodic derivative expects an even point count")
    if not np.isfinite(length) or length <= 0:
        raise ValueError("positive period required")
    indices = np.arange(points) - points // 2
    wavenumber = 2 * np.pi * indices / length
    wavenumber[indices == -points // 2] = 0.0
    coordinate = indices.astype(float) * (length / points)
    fourier = np.exp(1j * coordinate[:, None] * (2 * np.pi * indices / length)[None, :])
    fourier /= np.sqrt(points)
    matrix = (fourier * (1j * wavenumber)[None, :]) @ fourier.conj().T
    matrix = np.real_if_close(matrix, tol=1000)
    matrix = np.real(0.5 * (matrix - matrix.T))
    return matrix


def antiperiodic_momentum(points, length):
    """Existing AP momentum on the covariant Fourier grid. Independent of the metric profile."""
    points = int(points)
    if points < 10 or points % 2:
        raise ValueError("AP fermion grid must be even and at least 10 points")
    metric = CovariantStaticMetric(
        float(length),
        np.ones(points),
        np.ones(points),
        np.ones(points),
        eta=0.5,
    )
    return metric.momentum_matrix.copy(), metric


def smooth_step(t):
    """Flat-endpoint transition. S(0+)=0 and S(1-)=1."""
    t = np.asarray(t, dtype=float)
    out = np.empty(t.shape, dtype=float)
    small = t <= 1e-10
    large = t >= 1 - 1e-10
    mid = ~small & ~large
    out[small] = 0.0
    out[large] = 1.0
    if np.any(mid):
        tt = np.clip(t[mid], 1e-10, 1 - 1e-10)
        left = np.exp(-1 / tt)
        right = np.exp(-1 / (1 - tt))
        out[mid] = left / (left + right)
    return out


def profile_coordinate(xi, length=PERIOD):
    """s=x on the packet arc [0, 4], C^infty bridge on [4, length]."""
    xi = np.mod(np.asarray(xi, dtype=float), length)
    arc = xi <= 4.0
    out = np.empty(xi.shape, dtype=float)
    out[arc] = xi[arc]
    bridge = ~arc
    t = (xi[bridge] - 4.0) / (length - 4.0)
    out[bridge] = xi[bridge] - length * smooth_step(t)
    return out


def _bump(u):
    u = np.asarray(u, dtype=float)
    out = np.zeros(u.shape, dtype=float)
    mask = (u > 0.0) & (u < 1.0)
    uu = u[mask]
    out[mask] = np.exp(-1.0 / (uu * (1.0 - uu)))
    return out


def _lobe_norms():
    nodes, weights = np.polynomial.legendre.leggauss(256)
    u = 0.5 * (nodes + 1)
    jacobian = 0.5
    weight = _bump(u)
    even_norm = np.sqrt(np.sum(weights * jacobian * weight ** 2))
    odd_norm = np.sqrt(np.sum(weights * jacobian * ((u - 0.5) * weight) ** 2))
    return float(even_norm), float(odd_norm)


def _sample_lobe(xi, left, kind, dx, norms):
    u = xi - left
    weight = _bump(u)
    if kind == "even":
        raw = weight / norms[0]
    else:
        raw = (u - 0.5) * weight / norms[1]
    # √dx half-density samples of a continuum L2 function.
    return raw * np.sqrt(dx)


def lowdin(columns):
    gram = columns.conj().T @ columns
    evals, evecs = np.linalg.eigh(0.5 * (gram + gram.conj().T))
    if np.min(evals) <= 0:
        raise ValueError("envelope Gram matrix is not positive")
    inverse_sqrt = (evecs * (1 / np.sqrt(evals))) @ evecs.conj().T
    return columns @ inverse_sqrt, gram


@dataclass
class CouplingSystem:
    points: int
    length: float
    dx: float
    xi: np.ndarray
    derivative: np.ndarray
    momentum: np.ndarray
    length_density: np.ndarray
    shift: np.ndarray
    kappa: int
    multiplicity: int
    occupations: np.ndarray
    coefficients: dict
    calibration: dict
    preparation: dict

    @property
    def A(self):
        return self.coefficients["A"]

    @property
    def C_W(self):
        return self.coefficients["C_W"]

    @property
    def C_F(self):
        return self.coefficients["C_F"]

    @property
    def flux(self):
        return self.coefficients["flux"]


@dataclass
class CauchyState:
    Q: np.ndarray
    r: np.ndarray
    chi: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    p_chi: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray

    def copy(self):
        return CauchyState(
            self.Q.copy(), self.r.copy(), self.chi.copy(),
            self.p_Q.copy(), self.p_r.copy(), self.p_chi.copy(),
            self.phi0.copy(), self.phi1.copy(),
        )


@dataclass
class CauchyRate:
    Q: np.ndarray
    r: np.ndarray
    chi: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    p_chi: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray
    force_L: np.ndarray
    force_Q: np.ndarray
    force_beta: np.ndarray
    fieldwork_power: float


def _combine(state, rate, step):
    return CauchyState(
        state.Q + step * rate.Q,
        state.r + step * rate.r,
        state.chi + step * rate.chi,
        state.p_Q + step * rate.p_Q,
        state.p_r + step * rate.p_r,
        state.p_chi + step * rate.p_chi,
        state.phi0 + step * rate.phi0,
        state.phi1 + step * rate.phi1,
    )


def build_system(points=64):
    """Static lapse and shift gauge, AP momentum, periodic geometry derivative."""
    points = int(points)
    coefficients = locked_coefficients()
    if sector_multiplicity(KAPPA) != 4 * KAPPA:
        raise ValueError("multiplicity must remain M=4 kappa once")
    derivative = periodic_derivative(points, PERIOD)
    momentum, _metric = antiperiodic_momentum(points, PERIOD)
    dx = PERIOD / points
    xi = np.arange(points, dtype=float) * dx
    coordinate = profile_coordinate(xi, PERIOD)
    scale = np.exp(np.log(OMEGA) * coordinate)
    length_density = CALIBRATION["b0"] * scale
    shift = CALIBRATION["beta0"] * scale
    q0 = CALIBRATION["b0"] / CALIBRATION["a0"]
    phi0, phi1, preparation = prepare_rank6(xi, dx)
    system = CouplingSystem(
        points=points,
        length=PERIOD,
        dx=dx,
        xi=xi,
        derivative=derivative,
        momentum=momentum,
        length_density=length_density,
        shift=shift,
        kappa=KAPPA,
        multiplicity=int(sector_multiplicity(KAPPA)),
        occupations=OCCUPATIONS.copy(),
        coefficients=coefficients,
        calibration=dict(CALIBRATION),
        preparation=preparation,
    )
    state = CauchyState(
        Q=np.full(points, q0),
        r=np.ones(points),
        chi=np.zeros(points),
        p_Q=np.zeros(points),
        p_r=np.zeros(points),
        p_chi=np.zeros(points),
        phi0=phi0,
        phi1=phi1,
    )
    return system, state


def prepare_rank6(xi, dx):
    """Three real two-lobe packets, then region-local carriers and one minus-column phase.

    The recorded phase multiplies the minus spinor column. It is not a phase of
    the odd spatial lobe. A constant column phase does not change
    C = Φ diag(c) Φ†. Complement occupation 0 is the finite Gaussian support,
    not a physical vacuum identification.
    """
    norms = _lobe_norms()
    phase = CALIBRATION["phase"]
    real_packets = []
    lobe_norms = []
    for region in range(PACKET_COUNT):
        left = region * ELL
        even = _sample_lobe(xi, left, "even", dx, norms)
        odd = _sample_lobe(xi, left + ELL, "odd", dx, norms)
        even_norm = np.linalg.norm(even)
        odd_norm = np.linalg.norm(odd)
        if even_norm == 0 or odd_norm == 0:
            raise ValueError("packet lobe is unresolved on this grid")
        lobe_norms.append({"even_ell2_before": float(even_norm), "odd_ell2_before": float(odd_norm)})
        even = even / even_norm
        odd = odd / odd_norm
        real_packets.append((even + odd) / np.sqrt(2))
    real = np.column_stack(real_packets)
    real_gram = real.conj().T @ real
    envelopes, _ = lowdin(real)
    envelope_imag = float(np.max(np.abs(np.imag(envelopes))))
    if envelope_imag > 1e-10:
        raise ValueError("real packet envelopes did not stay real")
    envelopes = np.real(envelopes)
    phi0 = np.empty((xi.size, 6), dtype=complex)
    phi1 = np.empty((xi.size, 6), dtype=complex)
    minus_phase = np.exp(1j * phase)
    for region in range(PACKET_COUNT):
        local = xi - region * ELL
        carrier = np.exp(1j * CARRIER_K * local) * envelopes[:, region]
        plus0 = carrier / np.sqrt(2)
        plus1 = 1j * carrier / np.sqrt(2)
        phi0[:, 2 * region] = plus0
        phi1[:, 2 * region] = plus1
        phi0[:, 2 * region + 1] = minus_phase * np.conjugate(plus0)
        phi1[:, 2 * region + 1] = minus_phase * np.conjugate(plus1)
    gram = phi0.conj().T @ phi0 + phi1.conj().T @ phi1
    preparation = {
        "occupations": OCCUPATIONS.tolist(),
        "active_rank": 6,
        "complement_dimension": int(2 * xi.size - 6),
        "complement_occupation": 0.0,
        "vacuum_identification": False,
        "sea_subtracted": False,
        "gaussian_form": "C=Phi diag(c) Phi^H",
        "half_density": "sqrt(dx) samples, ell2 orthonormal",
        "real_envelopes_remain_real": True,
        "carrier": "exp(+i k (x-n*ell)) on plus; minus is exp(i phase) times the conjugate",
        "minus_columns": "exp(i*recorded_phase) times conjugate of plus columns",
        "real_envelope_gram_before_phase": {
            "real": np.real(real_gram).tolist(),
            "imag": np.imag(real_gram).tolist(),
        },
        "real_envelope_orthonormal_defect": float(np.max(np.abs(real_gram - np.eye(3)))),
        "real_envelope_imag_after_lowdin": envelope_imag,
        "column_orthonormal_defect": float(np.max(np.abs(gram - np.eye(6)))),
        "lobe_sample_norms": lobe_norms,
        "continuum_lobe_norms": {"even": norms[0], "odd": norms[1]},
        "phase_applied_to_odd_lobe": False,
        "phase_on_minus_spinor_column": True,
    }
    return phi0, phi1, preparation


def apply_dirac(phi0, phi1, length_density, radial_density, shift, kappa, momentum):
    """H Phi = σ2 {a,P}/2 Phi + σ1 κ L Phi - {β,P}/2 Phi, a=L/Q."""
    if np.min(radial_density) <= 0:
        raise PositiveChartExit("Q_left_positive_chart", np.nan, None)
    a = length_density / radial_density

    def anticommutator(values, component):
        return 0.5 * (
            values[:, None] * (momentum @ component)
            + momentum @ (values[:, None] * component)
        )

    kinetic0 = anticommutator(a, phi0)
    kinetic1 = anticommutator(a, phi1)
    shift0 = anticommutator(shift, phi0)
    shift1 = anticommutator(shift, phi1)
    mass = kappa * length_density
    image0 = -1j * kinetic1 + mass[:, None] * phi1 - shift0
    image1 = 1j * kinetic0 + mass[:, None] * phi0 - shift1
    return image0, image1


def column_moments(phi0, phi1, occupations, momentum):
    """Nodal K, S1 and Pmom. K and Pmom are the real parts named by the contract."""
    momentum0 = momentum @ phi0
    momentum1 = momentum @ phi1
    weighted = occupations[None, :]
    k_density = (
        np.conjugate(phi0) * (-1j * momentum1)
        + np.conjugate(phi1) * (1j * momentum0)
    )
    s1_density = np.conjugate(phi0) * phi1 + np.conjugate(phi1) * phi0
    p_density = np.conjugate(phi0) * momentum0 + np.conjugate(phi1) * momentum1
    K = np.sum(k_density * weighted, axis=1).real
    S1 = np.sum(s1_density * weighted, axis=1)
    Pmom = np.sum(p_density * weighted, axis=1).real
    return K, S1, Pmom


def nodal_forces(K, S1, Pmom, length_density, radial_density, kappa, multiplicity):
    """Partial derivatives of M Tr(C H). Not yet densities."""
    s1 = np.real(S1)
    force_l = multiplicity * (K / radial_density + kappa * s1)
    force_q = -multiplicity * length_density * K / radial_density ** 2
    force_beta = -multiplicity * Pmom
    return force_l, force_q, force_beta


def source_from_columns(system, state):
    image0, image1 = apply_dirac(
        state.phi0, state.phi1, system.length_density, state.Q, system.shift,
        system.kappa, system.momentum,
    )
    K, S1, Pmom = column_moments(state.phi0, state.phi1, system.occupations, system.momentum)
    force_l, force_q, force_beta = nodal_forces(
        K, S1, Pmom, system.length_density, state.Q, system.kappa, system.multiplicity,
    )
    return {
        "image0": image0,
        "image1": image1,
        "K": K,
        "S1": S1,
        "Pmom": Pmom,
        "force_L": force_l,
        "force_Q": force_q,
        "force_beta": force_beta,
        "multiplicity_applied_once": True,
    }


def field_energy(system, state, images=None):
    if images is None:
        images = apply_dirac(
            state.phi0, state.phi1, system.length_density, state.Q, system.shift,
            system.kappa, system.momentum,
        )
    amplitude = np.sum(
        np.conjugate(state.phi0) * images[0] + np.conjugate(state.phi1) * images[1],
        axis=0,
    ).real
    return float(system.multiplicity * np.dot(system.occupations, amplitude))


def _pi_and_potentials(system, state):
    force_r, f_chi = partial_F(state.r, system.A, system.C_W)
    f_chi = float(f_chi)
    pi = state.p_r - (force_r / f_chi) * state.p_chi
    potential = feedback_V(state.r, state.chi, system.A, system.C_W, system.C_F, system.flux)
    # V_rel is not added. C_E and C_box stay boundary coefficients on this periodic slice.
    return force_r, f_chi, pi, potential


def gravity_energy(system, state):
    """dx times the SBP sum. DF means D(F), not a silent product rule."""
    D = system.derivative
    radius_derivative = D @ state.r
    chi_derivative = D @ state.chi
    momentum_derivative = D @ state.p_Q
    lapse_derivative = D @ system.length_density
    F = feedback_F(state.r, state.chi, system.A, system.C_W)
    F_derivative = D @ F
    force_r, f_chi, pi, potential = _pi_and_potentials(system, state)
    del force_r
    Z = float(feedback_Z(system.A))
    kinetic = (
        state.p_Q * state.p_chi / (2 * f_chi)
        + pi ** 2 / (4 * Z * state.Q)
        + Z * radius_derivative ** 2 / state.Q
        - state.Q * potential
    )
    density = (
        system.length_density * kinetic
        + 2 * F_derivative * lapse_derivative / state.Q
        + system.shift * (
            state.p_r * radius_derivative
            + state.p_chi * chi_derivative
            - state.Q * momentum_derivative
        )
    )
    return float(system.dx * np.sum(density))


def total_energy(system, state):
    return gravity_energy(system, state) + field_energy(system, state)


def hamilton_constraint(system, state):
    """Geometric lapse constraint density, before the matter density rho."""
    D = system.derivative
    radius_derivative = D @ state.r
    F = feedback_F(state.r, state.chi, system.A, system.C_W)
    F_over_Q = (D @ F) / state.Q
    _force_r, f_chi, pi, potential = _pi_and_potentials(system, state)
    Z = float(feedback_Z(system.A))
    return (
        state.p_Q * state.p_chi / (2 * f_chi)
        + pi ** 2 / (4 * Z * state.Q)
        + Z * radius_derivative ** 2 / state.Q
        - state.Q * potential
        - 2 * (D @ F_over_Q)
    )


def shift_constraint(system, state):
    D = system.derivative
    return (
        state.p_r * (D @ state.r)
        + state.p_chi * (D @ state.chi)
        - state.Q * (D @ state.p_Q)
    )


def constraint_residuals(system, state, forces=None):
    if forces is None:
        forces = source_from_columns(system, state)
    rho = forces["force_L"] / system.dx
    current = forces["force_beta"] / system.dx
    hamilton = hamilton_constraint(system, state) + rho
    momentum = shift_constraint(system, state) + current
    return {
        "rho": rho,
        "current": current,
        "hamilton": hamilton,
        "momentum": momentum,
        "hamilton_max": float(np.max(np.abs(hamilton))),
        "momentum_max": float(np.max(np.abs(momentum))),
        "current_max": float(np.max(np.abs(current))),
    }


def geometric_rates(system, state):
    """Discrete Hamilton equations of the SBP energy. Matter force is not included."""
    D = system.derivative
    radius_derivative = D @ state.r
    chi_derivative = D @ state.chi
    momentum_derivative = D @ state.p_Q
    lapse_derivative = D @ system.length_density
    F = feedback_F(state.r, state.chi, system.A, system.C_W)
    F_derivative = D @ F
    force_r, f_chi, pi, potential = _pi_and_potentials(system, state)
    Z = float(feedback_Z(system.A))
    length_density = system.length_density
    Q = state.Q
    Q_dot = length_density * state.p_chi / (2 * f_chi) + D @ (system.shift * Q)
    r_dot = length_density * pi / (2 * Z * Q) + system.shift * radius_derivative
    chi_dot = (
        length_density * state.p_Q / (2 * f_chi)
        - (force_r / f_chi) * (length_density * pi / (2 * Z * Q))
        + system.shift * chi_derivative
    )
    p_Q_dot = (
        length_density * pi ** 2 / (4 * Z * Q ** 2)
        + length_density * Z * radius_derivative ** 2 / Q ** 2
        + length_density * potential
        + 2 * F_derivative * lapse_derivative / Q ** 2
        + system.shift * momentum_derivative
    )
    pi_radius = (8 * np.pi * system.A / f_chi) * state.p_chi
    weight = 2 * lapse_derivative / Q
    p_r_dot = (
        -length_density * pi * pi_radius / (2 * Z * Q)
        + D @ (2 * length_density * Z * radius_derivative / Q)
        + length_density * Q * (16 * np.pi * system.A * state.r)
        + force_r * (D @ weight)
        + D @ (system.shift * state.p_r)
    )
    alpha = float(alpha_of(system.C_W))
    p_chi_dot = (
        length_density * Q * (-alpha * (4 + 2 * state.chi))
        + f_chi * (D @ weight)
        + D @ (system.shift * state.p_chi)
    )
    return Q_dot, r_dot, chi_dot, p_Q_dot, p_r_dot, p_chi_dot


def chart_failure(system, state):
    arrays = (state.Q, state.r, state.chi, state.p_Q, state.p_r, state.p_chi, state.phi0, state.phi1)
    if not all(np.isfinite(array).all() for array in arrays):
        return "nonfinite_state"
    if np.min(state.Q) <= 0:
        return "Q_left_positive_chart"
    if np.min(state.r) <= 0:
        return "r_left_positive_chart"
    if np.min(system.length_density) <= 0:
        return "L_left_positive_chart"
    return None


def rates(system, state, include_matter_force=True):
    failure = chart_failure(system, state)
    if failure is not None:
        raise PositiveChartExit(failure, np.nan, state.copy())
    source = source_from_columns(system, state)
    Q_dot, r_dot, chi_dot, p_Q_dot, p_r_dot, p_chi_dot = geometric_rates(system, state)
    if include_matter_force:
        p_Q_dot = p_Q_dot - source["force_Q"] / system.dx
    # Representative block evolves by -i H, not by -M i H.
    phi0_dot = -1j * source["image0"]
    phi1_dot = -1j * source["image1"]
    fieldwork = float(np.sum(source["force_Q"] * Q_dot))
    return CauchyRate(
        Q_dot, r_dot, chi_dot, p_Q_dot, p_r_dot, p_chi_dot, phi0_dot, phi1_dot,
        source["force_L"], source["force_Q"], source["force_beta"], fieldwork,
    )


def normal_radius_velocity(system, state, r_dot=None):
    radius_derivative = system.derivative @ state.r
    if r_dot is None:
        _Q, r_dot, _chi, _pQ, _pr, _pc = geometric_rates(system, state)
    normal = (r_dot - system.shift * radius_derivative) / (state.r * system.length_density)
    coordinate = r_dot / (state.r * system.length_density)
    shift_piece = (system.shift * radius_derivative) / (state.r * system.length_density)
    return normal, coordinate, shift_piece


def weyl_proxy(state):
    return state.chi ** 2 / (3 * state.r ** 4)


def _radius_residual(system, r, rho):
    probe = CauchyState(
        Q=np.full(system.points, system.calibration["b0"] / system.calibration["a0"]),
        r=r,
        chi=np.zeros(system.points),
        p_Q=np.zeros(system.points),
        p_r=np.zeros(system.points),
        p_chi=np.zeros(system.points),
        phi0=np.zeros((system.points, 1)),
        phi1=np.zeros((system.points, 1)),
    )
    # Zero columns: hamilton_constraint does not read phi. rho is supplied.
    return hamilton_constraint(system, probe) + rho


def _radius_jacobian(system, r):
    Q = np.full(system.points, system.calibration["b0"] / system.calibration["a0"])
    D = system.derivative
    radius_derivative = D @ r
    force_r = -8 * np.pi * system.A * r
    Z = float(feedback_Z(system.A))
    return (
        np.diag(2 * Z * radius_derivative / Q) @ D
        - np.diag(Q * 16 * np.pi * system.A * r)
        - 2 * D @ np.diag(1 / Q) @ D @ np.diag(force_r)
    )


def solve_initial_radius(system, state):
    """Homotopy Newton of the discrete Hamilton constraint. Momenta and chi stay 0.

    r=1 solves the source-free chart because r_mag^2=1. The actual rho is
    reached by continuation. The residual uses D(F), not the product rule.
    """
    source = source_from_columns(system, state)
    rho = source["force_L"] / system.dx
    rho_variation = float(np.max(np.abs(
        source_from_columns(system, CauchyState(
            state.Q, state.r * 1.7, state.chi, state.p_Q, state.p_r, state.p_chi,
            state.phi0, state.phi1,
        ))["force_L"] - source["force_L"]
    )))
    r = np.ones(system.points)
    steps = []
    for scale in np.linspace(0.0, 1.0, 9)[1:]:
        target = scale * rho
        converged = False
        residual_max = None
        for _iteration in range(12):
            residual = _radius_residual(system, r, target)
            residual_max = float(np.max(np.abs(residual)))
            if residual_max < TOL_NEWTON:
                converged = True
                break
            jacobian = _radius_jacobian(system, r)
            try:
                delta = np.linalg.solve(jacobian, -residual)
            except np.linalg.LinAlgError as error:
                return state, {
                    "converged": False,
                    "blocker": f"singular radius Jacobian: {error}",
                    "residual_max": residual_max,
                    "steps": steps,
                    "rho_independent_of_r_max": rho_variation,
                }
            accepted = False
            step = 1.0
            for _halving in range(16):
                trial = r + step * delta
                if np.min(trial) <= 0 or not np.isfinite(trial).all():
                    step *= 0.5
                    continue
                trial_residual = float(np.max(np.abs(_radius_residual(system, trial, target))))
                if trial_residual < residual_max * (1 - 0.05 * step):
                    r = trial
                    accepted = True
                    break
                step *= 0.5
            if not accepted:
                break
        steps.append({"scale": float(scale), "converged": converged, "residual_max": residual_max})
        if not converged:
            return state, {
                "converged": False,
                "blocker": "Newton left the positive chart or stalled before full rho",
                "residual_max": residual_max,
                "steps": steps,
                "rho_independent_of_r_max": rho_variation,
            }
    solved = state.copy()
    solved.r = r
    final = constraint_residuals(system, solved)
    eight_pi_A = 8 * np.pi * system.A
    bracket = solved.Q ** 2 * magnetic_radius_square(system.coefficients) + solved.Q * final["rho"] / eight_pi_A
    return solved, {
        "converged": True,
        "blocker": None,
        "residual_max": final["hamilton_max"],
        "momentum_residual_max": final["momentum_max"],
        "current_max": final["current_max"],
        "steps": steps,
        "rho_independent_of_r_max": rho_variation,
        "bracket_min": float(np.min(bracket)),
        "bracket_positive": bool(np.min(bracket) > 0),
        "r_min": float(np.min(r)),
        "r_max": float(np.max(r)),
        "r_mean": float(np.mean(r)),
        "imposed_radius": False,
    }


def jacobian_directional_error(system, r, rho, seed=5):
    jacobian = _radius_jacobian(system, r)
    rng = np.random.default_rng(seed)
    direction = rng.normal(size=r.size)
    direction /= np.linalg.norm(direction)
    eps = 1e-6
    plus = _radius_residual(system, r + eps * direction, rho)
    minus = _radius_residual(system, r - eps * direction, rho)
    numeric = (plus - minus) / (2 * eps)
    analytic = jacobian @ direction
    return float(np.max(np.abs(numeric - analytic)))


def covariance_matrix(system, state):
    vectors = np.vstack((state.phi0, state.phi1))
    return (vectors * system.occupations) @ vectors.conj().T


def dense_force_difference(system, state):
    """Column quadratic forms against conformal_source. Small grids only."""
    if system.points > 32:
        raise ValueError("dense covariance cross-check is limited to small grids")
    metric = CovariantStaticMetric(
        system.length,
        system.length_density * state.r,
        state.Q * state.r,
        state.r,
        eta=0.5,
    )
    dense = conformal_source(metric, system.shift, covariance_matrix(system, state), kappa=system.kappa)
    columns = source_from_columns(system, state)
    difference = {
        name: float(np.max(np.abs(dense["nodal"][name] - columns[force])))
        for name, force in (("L", "force_L"), ("Q", "force_Q"), ("beta", "force_beta"))
    }
    difference["energy"] = float(abs(dense["energy"] - field_energy(system, state)))
    difference["multiplicity"] = int(dense["multiplicity"])
    difference["vacuum_subtracted"] = bool(dense["vacuum_subtracted"])
    return difference


def dense_hamiltonian_difference(system, state):
    metric = CovariantStaticMetric(
        system.length,
        system.length_density * np.ones(system.points),
        state.Q,
        np.ones(system.points),
        eta=0.5,
    )
    from .nsc_conformal_adm_source import direct_hamiltonian

    matrix = direct_hamiltonian(metric, system.shift, system.kappa)
    vector = np.vstack((state.phi0, state.phi1))
    dense_image = matrix @ vector
    image0, image1 = apply_dirac(
        state.phi0, state.phi1, system.length_density, state.Q, system.shift,
        system.kappa, system.momentum,
    )
    column_image = np.vstack((image0, image1))
    return float(np.max(np.abs(dense_image - column_image)))


def product_rule_gap(system, state):
    D = system.derivative
    product = D @ (system.shift * state.Q)
    factored = system.shift * (D @ state.Q) + state.Q * (D @ system.shift)
    return float(np.max(np.abs(product - factored)))


def kinematic_module_gap(system, state):
    """Checked velocity map versus the discrete product D(beta Q)."""
    Q_dot, r_dot, chi_dot, _pQ, _pr, _pc = geometric_rates(system, state)
    D = system.derivative
    Q_x = D @ state.Q
    r_x = D @ state.r
    chi_x = D @ state.chi
    beta_x = D @ system.shift
    gaps = []
    for index in (0, system.points // 3, system.points // 2, 2 * system.points // 3):
        checked = velocity_rhs(
            float(system.length_density[index]),
            float(state.Q[index]),
            float(Q_x[index]),
            float(state.r[index]),
            float(r_x[index]),
            float(chi_x[index]),
            float(system.shift[index]),
            float(beta_x[index]),
            float(state.p_Q[index]),
            float(state.p_r[index]),
            float(state.p_chi[index]),
            system.A,
            system.C_W,
        )
        gaps.append((
            abs(Q_dot[index] - float(checked[0])),
            abs(r_dot[index] - float(checked[1])),
            abs(chi_dot[index] - float(checked[2])),
        ))
    gaps = np.array(gaps)
    return {
        "Q_includes_product_rule_gap": float(np.max(gaps[:, 0])),
        "r": float(np.max(gaps[:, 1])),
        "chi": float(np.max(gaps[:, 2])),
        "product_rule_gap": product_rule_gap(system, state),
    }


def _energy_directional(system, state, name, direction, eps):
    def displaced(sign):
        trial = state.copy()
        setattr(trial, name, getattr(trial, name) + sign * eps * direction)
        return total_energy(system, trial)

    return (displaced(1) - displaced(-1)) / (2 * eps)


def _parameter_directional(system, state, name, direction, eps):
    def displaced(sign):
        trial_state = state.copy()
        length_density = system.length_density
        shift = system.shift
        if name == "L":
            length_density = system.length_density + sign * eps * direction
        else:
            shift = system.shift + sign * eps * direction
        cloned = CouplingSystem(
            **{**system.__dict__, "length_density": length_density, "shift": shift}
        )
        return total_energy(cloned, trial_state)

    return (displaced(1) - displaced(-1)) / (2 * eps)


def hamiltonian_fd_errors(system, state, seed=11, eps=1e-6):
    """Directional derivatives of the discrete energy against the analytic RHS."""
    rng = np.random.default_rng(seed)
    rate = rates(system, state, include_matter_force=True)
    errors = {}
    samples = {}
    # E = dx * Σ, so qdot = (1/dx) dE/dp and pdot = -(1/dx) dE/dq.
    state_map = {
        "Q": (rate.p_Q, -system.dx),
        "r": (rate.p_r, -system.dx),
        "chi": (rate.p_chi, -system.dx),
        "p_Q": (rate.Q, system.dx),
        "p_r": (rate.r, system.dx),
        "p_chi": (rate.chi, system.dx),
    }
    for name, (velocity, factor) in state_map.items():
        direction = rng.normal(size=system.points)
        direction /= np.linalg.norm(direction)
        numeric = _energy_directional(system, state, name, direction, eps)
        analytic = factor * float(np.dot(velocity, direction))
        errors[name] = float(abs(numeric - analytic))
        samples[name] = float(abs(analytic))
    for name, velocity in (
        ("L", rate.force_L + system.dx * hamilton_constraint(system, state)),
        ("beta", rate.force_beta + system.dx * shift_constraint(system, state)),
    ):
        direction = rng.normal(size=system.points)
        direction /= np.linalg.norm(direction)
        numeric = _parameter_directional(system, state, name, direction, eps)
        analytic = float(np.dot(velocity, direction))
        errors[name] = float(abs(numeric - analytic))
        samples[name] = float(abs(analytic))
    worst = max(errors.values())
    scale = max(1.0, max(samples.values()))
    return {
        "absolute": errors,
        "relative_worst": float(worst / scale),
        "scale": float(scale),
    }


def mutation_power(system, state, dt=1e-6):
    """Dropping the matter force must leave an uncancelled fieldwork power.

    Centered step so the straight-line O(dt) bias does not hide the identity.
    """
    intact = rates(system, state, include_matter_force=True)
    mutated = rates(system, state, include_matter_force=False)

    def slope(rate):
        forward = total_energy(system, _combine(state, rate, dt))
        backward = total_energy(system, _combine(state, rate, -dt))
        return (forward - backward) / (2 * dt)

    intact_power = float(slope(intact))
    mutated_power = float(slope(mutated))
    fieldwork = intact.fieldwork_power
    return {
        "intact_power": intact_power,
        "mutated_power": mutated_power,
        "fieldwork_power": fieldwork,
        "mutation_minus_intact": mutated_power - intact_power,
        "fieldwork_identity_error": float(abs((mutated_power - intact_power) - fieldwork)),
        "breaks_balance": bool(abs(mutated_power) > 100 * max(abs(intact_power), 1e-10)),
    }


def reduced_dirac(system, state):
    source = source_from_columns(system, state)
    matrix = state.phi0.conj().T @ source["image0"] + state.phi1.conj().T @ source["image1"]
    local = sp.Matrix([[1, sp.I / 5], [-sp.I / 5, 2]])
    link = sp.Matrix([[sp.Rational(1, 4), sp.I / 7], [sp.Rational(1, 9), sp.Rational(1, 6)]])
    window = finite_window(local, link, sp.Rational(3, 2), 0, 3)
    frozen = np.array([[complex(window[i, j]) for j in range(6)] for i in range(6)], dtype=complex)
    difference = matrix - frozen
    onsite = []
    links = []
    for region in range(3):
        block = slice(2 * region, 2 * region + 2)
        onsite.append(float(np.linalg.norm(difference[block, block])))
    for region in range(2):
        row = slice(2 * region, 2 * region + 2)
        col = slice(2 * region + 2, 2 * region + 4)
        links.append(float(np.linalg.norm(difference[row, col])))
    return {
        "matrix_real": matrix.real.tolist(),
        "matrix_imag": matrix.imag.tolist(),
        "frozen_embedded_exactly": False,
        "old_circulation_transferred": False,
        "onsite_frobenius_difference": onsite,
        "link_frobenius_difference": links,
        "max_link_difference": float(max(links)),
        "max_onsite_difference": float(max(onsite)),
        "same_operator_class_not_same_matrix": True,
    }


def _frozen_onsite():
    return np.array([[1, 1j / 5], [-1j / 5, 2]], dtype=complex)


def _old_link():
    """Hand-set window link. Compared, never copied into the state or the Hamiltonian."""
    return np.array(
        [[0.25, 1j / 7], [1 / 9, 1 / 6]],
        dtype=complex,
    )


def continuum_packet_diagnostic(coefficients=None):
    """Real envelopes, region-local carriers, and the minus-column phase.

    Each onsite block is compared with Ω^n H. The old hand-set link is not used.
    """
    del coefficients
    norms = _lobe_norms()
    count = 2 ** 14
    length = PERIOD
    dx = length / count
    xi = np.arange(count) * dx
    phase = CALIBRATION["phase"]
    scale = np.exp(np.log(OMEGA) * xi)
    a = CALIBRATION["a0"] * scale
    beta = CALIBRATION["beta0"] * scale
    length_density = CALIBRATION["b0"] * scale
    symbol = 2 * np.pi * np.fft.fftfreq(count, d=dx)

    def derivative(values):
        return np.fft.ifft(1j * symbol * np.fft.fft(values))

    def momentum(values):
        return -1j * derivative(values)

    def apply_scalar(spatial, spinor):
        antic_a = 0.5 * (a * momentum(spatial) + momentum(a * spatial))
        antic_b = 0.5 * (beta * momentum(spatial) + momentum(beta * spatial))
        sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
        return (
            (sigma2 @ spinor)[:, None] * antic_a
            + (sigma1 @ spinor)[:, None] * (length_density * spatial)
            - spinor[:, None] * antic_b
        )

    plus = np.array([1, 1j]) / np.sqrt(2)
    minus = np.array([1, -1j]) / np.sqrt(2)

    def braket(left_spatial, left_spinor, image):
        density = np.conjugate(left_spatial) * np.einsum("a,an->n", np.conjugate(left_spinor), image)
        return complex(np.sum(density) * dx)

    regions = []
    reference = _frozen_onsite()
    mweight = None
    for region in range(PACKET_COUNT):
        left = region * ELL
        even = np.where((xi > left) & (xi < left + ELL), _bump(xi - left) / norms[0], 0.0)
        odd = np.where(
            (xi > left + ELL) & (xi < left + 2 * ELL),
            ((xi - left - ELL) - 0.5) * _bump(xi - left - ELL) / norms[1],
            0.0,
        )
        envelope = (even + odd) / np.sqrt(2)
        spatial_plus = np.exp(1j * CARRIER_K * (xi - left)) * envelope
        spatial_minus = np.exp(1j * phase) * np.conjugate(spatial_plus)
        if region == 0:
            mweight = float(np.sum(scale * np.abs(spatial_plus) ** 2) * dx)
        image_plus = apply_scalar(spatial_plus, plus)
        image_minus = apply_scalar(spatial_minus, minus)
        block = np.array([
            [braket(spatial_plus, plus, image_plus), braket(spatial_plus, plus, image_minus)],
            [braket(spatial_minus, minus, image_plus), braket(spatial_minus, minus, image_minus)],
        ])
        target = OMEGA ** region * reference
        regions.append({
            "region": region,
            "scale": float(OMEGA ** region),
            "diagonal_abs_error": [
                float(abs(block[0, 0] - target[0, 0])),
                float(abs(block[1, 1] - target[1, 1])),
            ],
            "offdiag_abs_error": float(abs(block[0, 1] - target[0, 1])),
            "onsite_real": block.real.tolist(),
            "onsite_imag": block.imag.tolist(),
        })
    return {
        "Mweight": mweight,
        "Mweight_minus_executed": mweight - CALIBRATION["Mweight_executed"],
        "phase_on_minus_spinor_column": True,
        "phase_applied_to_odd_lobe": False,
        "carrier": "exp(±i k (x-n*ell))",
        "regions": regions,
        "onsite_real": regions[0]["onsite_real"],
        "onsite_imag": regions[0]["onsite_imag"],
        "diagonal_abs_error": regions[0]["diagonal_abs_error"],
        "offdiag_abs_error": regions[0]["offdiag_abs_error"],
        "matches_frozen_offdiag": all(row["offdiag_abs_error"] < 1e-8 for row in regions),
        "omega_scaled_onsite": True,
        "old_B_not_used": True,
    }


def packet_onsite_report(system, state):
    """Onsite blocks versus Ω^n H. Link blocks are reported against the old B and not replaced by it."""
    source = source_from_columns(system, state)
    matrix = state.phi0.conj().T @ source["image0"] + state.phi1.conj().T @ source["image1"]
    reference = _frozen_onsite()
    old = _old_link()
    onsite = []
    for region in range(PACKET_COUNT):
        block = matrix[2 * region: 2 * region + 2, 2 * region: 2 * region + 2]
        target = OMEGA ** region * reference
        onsite.append({
            "region": region,
            "diagonal_abs_error": [
                float(abs(block[0, 0] - target[0, 0])),
                float(abs(block[1, 1] - target[1, 1])),
            ],
            "offdiag_abs_error": float(abs(block[0, 1] - target[0, 1])),
        })
    links = []
    for region in range(PACKET_COUNT - 1):
        block = matrix[2 * region: 2 * region + 2, 2 * region + 2: 2 * region + 4]
        target = OMEGA ** region * old
        links.append({
            "region": region,
            "difference_from_old_B": float(np.linalg.norm(block - target)),
            "copied_from_old_B": False,
        })
    return {"onsite": onsite, "links": links, "old_B_not_used": True}


def observability(system, state, rate=None):
    if rate is None:
        rate = rates(system, state)
    normal, coordinate, shift_piece = normal_radius_velocity(system, state, rate.r)
    residuals = constraint_residuals(system, state, {
        "force_L": rate.force_L,
        "force_beta": rate.force_beta,
    })
    gram = state.phi0.conj().T @ state.phi0 + state.phi1.conj().T @ state.phi1
    number = float(np.sum(system.occupations * np.sum(np.abs(state.phi0) ** 2 + np.abs(state.phi1) ** 2, axis=0)))
    return {
        "time": None,
        "energy": total_energy(system, state),
        "field_energy": field_energy(system, state),
        "gravity_energy": gravity_energy(system, state),
        "hamilton_max": residuals["hamilton_max"],
        "momentum_max": residuals["momentum_max"],
        "current_max": residuals["current_max"],
        "unitarity": float(np.max(np.abs(gram - np.eye(6)))),
        "number": number,
        "fieldwork_power": rate.fieldwork_power,
        "geometric_counterwork_power": -rate.fieldwork_power,
        "normal_radius_velocity_max": float(np.max(np.abs(normal))),
        "coordinate_radius_velocity_max": float(np.max(np.abs(coordinate))),
        "shift_radius_velocity_max": float(np.max(np.abs(shift_piece))),
        "weyl_proxy_max": float(np.max(weyl_proxy(state))),
        "r_min": float(np.min(state.r)),
        "Q_min": float(np.min(state.Q)),
        "chi_max": float(np.max(np.abs(state.chi))),
    }


def rk4_step(system, state, dt):
    k1 = rates(system, state)
    k2 = rates(system, _combine(state, k1, 0.5 * dt))
    k3 = rates(system, _combine(state, k2, 0.5 * dt))
    k4 = rates(system, _combine(state, k3, dt))
    combined = CauchyRate(
        (k1.Q + 2 * k2.Q + 2 * k3.Q + k4.Q) / 6,
        (k1.r + 2 * k2.r + 2 * k3.r + k4.r) / 6,
        (k1.chi + 2 * k2.chi + 2 * k3.chi + k4.chi) / 6,
        (k1.p_Q + 2 * k2.p_Q + 2 * k3.p_Q + k4.p_Q) / 6,
        (k1.p_r + 2 * k2.p_r + 2 * k3.p_r + k4.p_r) / 6,
        (k1.p_chi + 2 * k2.p_chi + 2 * k3.p_chi + k4.p_chi) / 6,
        (k1.phi0 + 2 * k2.phi0 + 2 * k3.phi0 + k4.phi0) / 6,
        (k1.phi1 + 2 * k2.phi1 + 2 * k3.phi1 + k4.phi1) / 6,
        k1.force_L, k1.force_Q, k1.force_beta, k1.fieldwork_power,
    )
    return _combine(state, combined, dt)


def stable_timestep(system, state):
    wavenumber = np.max(np.abs(np.linalg.eigvalsh(system.momentum)))
    a = np.max(system.length_density / state.Q)
    beta = np.max(np.abs(system.shift))
    mass = np.max(np.abs(system.kappa * system.length_density))
    omega = (a + beta) * float(wavenumber) + mass
    # RK4 imaginary stability is near 2√2. Stay at about half of that.
    return float(min(5e-4, 1.4 / max(omega, 1.0))), float(omega)


def evolve(system, state, duration, dt):
    steps = int(np.ceil(duration / dt - 1e-12))
    dt = float(duration / steps)
    current = state.copy()
    samples = [observability(system, current)]
    samples[0]["time"] = 0.0
    fieldwork = 0.0
    counterwork = 0.0
    previous_power = samples[0]["fieldwork_power"]
    try:
        for step in range(steps):
            try:
                current = rk4_step(system, current, dt)
            except PositiveChartExit as exit_chart:
                raise PositiveChartExit(exit_chart.reason, step * dt, current) from exit_chart
            failure = chart_failure(system, current)
            if failure is not None:
                raise PositiveChartExit(failure, (step + 1) * dt, current)
            observation = observability(system, current)
            observation["time"] = (step + 1) * dt
            fieldwork += 0.5 * dt * (previous_power + observation["fieldwork_power"])
            counterwork -= 0.5 * dt * (previous_power + observation["fieldwork_power"])
            previous_power = observation["fieldwork_power"]
            samples.append(observation)
    except PositiveChartExit as exit_chart:
        return {
            "stopped": True,
            "reason": exit_chart.reason,
            "time": exit_chart.time_value,
            "state": current,
            "samples": samples,
            "fieldwork_integral": fieldwork,
            "geometric_counterwork_integral": counterwork,
            "dt": dt,
            "success": False,
        }
    energy0 = samples[0]["energy"]
    energy_drift = float(max(abs(sample["energy"] - energy0) for sample in samples))
    return {
        "stopped": False,
        "reason": None,
        "time": duration,
        "state": current,
        "samples": samples,
        "fieldwork_integral": float(fieldwork),
        "geometric_counterwork_integral": float(counterwork),
        "dt": dt,
        "success": True,
        "energy_drift": energy_drift,
        "hamilton_max": float(max(sample["hamilton_max"] for sample in samples)),
        "momentum_max": float(max(sample["momentum_max"] for sample in samples)),
        "unitarity_max": float(max(sample["unitarity"] for sample in samples)),
        "number_drift": float(max(abs(sample["number"] - samples[0]["number"]) for sample in samples)),
        "final_normal_velocity_max": samples[-1]["normal_radius_velocity_max"],
        "final_weyl_proxy_max": samples[-1]["weyl_proxy_max"],
        "final_coordinate_velocity_max": samples[-1]["coordinate_radius_velocity_max"],
    }


def _compact_run(result):
    samples = result["samples"]
    keep = samples if len(samples) <= 40 else samples[:: max(1, len(samples) // 20)] + samples[-1:]
    compact = {
        "stopped": result["stopped"],
        "reason": result["reason"],
        "time": result["time"],
        "dt": result["dt"],
        "success_interval_completed": result["success"],
        "energy_drift": result.get("energy_drift"),
        "hamilton_max": result.get("hamilton_max"),
        "momentum_max": result.get("momentum_max"),
        "unitarity_max": result.get("unitarity_max"),
        "number_drift": result.get("number_drift"),
        "fieldwork_integral": result["fieldwork_integral"],
        "geometric_counterwork_integral": result["geometric_counterwork_integral"],
        "field_energy_change": samples[-1]["field_energy"] - samples[0]["field_energy"] if samples else None,
        "final_normal_velocity_max": result.get("final_normal_velocity_max"),
        "final_coordinate_velocity_max": result.get("final_coordinate_velocity_max"),
        "final_weyl_proxy_max": result.get("final_weyl_proxy_max"),
        "checkpoints": [
            {key: sample[key] for key in (
                "time", "energy", "hamilton_max", "momentum_max", "unitarity", "number",
                "fieldwork_power", "normal_radius_velocity_max", "weyl_proxy_max",
                "r_min", "Q_min", "chi_max",
            )}
            for sample in keep
        ],
    }
    return compact


def _end_arrays(state):
    return {
        "r": state.r.tolist(),
        "Q": state.Q.tolist(),
        "chi": state.chi.tolist(),
        "p_Q": state.p_Q.tolist(),
        "p_r": state.p_r.tolist(),
        "p_chi": state.p_chi.tolist(),
    }


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imag": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if not np.isfinite(number):
            return None
        return number
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def run_bounded_control(path=None, points=64):
    """N=64, T=0.005 validation. Refine only if that run stays on the chart.

    The checkpointed v1 JSON is not a writable destination. A successor record
    is not created here.
    """
    started = time.perf_counter()
    cpu = time.process_time()
    path = Path(path) if path is not None else CONTROL_RECORD
    if path.resolve() == CONTROL_RECORD.resolve():
        raise ValueError(
            "nsc-spherical-coupling-control-v1.json is checkpointed and cannot be overwritten"
        )
    system, bare = build_system(points)
    solved, solve_info = solve_initial_radius(system, bare)
    initial = None if not solve_info["converged"] else constraint_residuals(system, solved)
    normal = coordinate = shift_piece = None
    if solve_info["converged"]:
        normal, coordinate, shift_piece = normal_radius_velocity(system, solved)
    randomized = None
    fd = None
    mutation = None
    kinematic = None
    jacobian_error = None
    if solve_info["converged"]:
        rng = np.random.default_rng(11)
        randomized = solved.copy()
        randomized.p_Q = 0.05 * rng.normal(size=points)
        randomized.p_r = 0.05 * rng.normal(size=points)
        randomized.p_chi = 0.05 * rng.normal(size=points)
        randomized.chi = 0.02 * np.sin(2 * np.pi * system.xi / system.length)
        fd = hamiltonian_fd_errors(system, randomized)
        mutation = mutation_power(system, solved)
        kinematic = kinematic_module_gap(system, solved)
        jacobian_error = jacobian_directional_error(system, solved.r, initial["rho"])
    small_system, small_state = build_system(32)
    dense = dense_force_difference(small_system, small_state)
    dense["hamiltonian"] = dense_hamiltonian_difference(small_system, small_state)
    dt, omega = stable_timestep(system, solved if solve_info["converged"] else bare)
    runs = {}
    verdict = "FAIL_INITIAL_SOLVE"
    next_blocker = solve_info.get("blocker") or "initial discrete Hamilton constraint"
    if solve_info["converged"] and initial["hamilton_max"] <= TOL_INITIAL_CONSTRAINT and initial["momentum_max"] <= TOL_INITIAL_CONSTRAINT:
        if fd["relative_worst"] > TOL_FD_RELATIVE:
            verdict = "FAIL_DISCRETE_RHS"
            next_blocker = "analytic discrete RHS disagrees with directional energy derivatives"
        else:
            primary = evolve(system, solved, 0.005, dt)
            refined_time = evolve(system, solved, 0.005, 0.5 * dt)
            runs["N64_T0.005"] = _compact_run(primary)
            runs["N64_T0.005_dt_half"] = _compact_run(refined_time)
            runs["N64_T0.005"]["end"] = _end_arrays(primary["state"])
            resolved = (
                not primary["stopped"]
                and not refined_time["stopped"]
                and primary["energy_drift"] <= TOL_ENERGY_RELATIVE * max(1.0, abs(primary["samples"][0]["energy"]))
                and primary["hamilton_max"] <= TOL_RESOLVED_CONSTRAINT
                and primary["momentum_max"] <= TOL_RESOLVED_CONSTRAINT
            )
            passed = (
                resolved
                and primary["hamilton_max"] <= TOL_CONSTRAINT_VALIDATION
                and primary["momentum_max"] <= TOL_CONSTRAINT_VALIDATION
                and refined_time["hamilton_max"] <= TOL_CONSTRAINT_VALIDATION
                and refined_time["momentum_max"] <= TOL_CONSTRAINT_VALIDATION
            )
            if primary["stopped"] or refined_time["stopped"]:
                verdict = "STOP_POSITIVE_CHART"
                next_blocker = primary["reason"] or refined_time["reason"]
            elif not resolved:
                verdict = "PROVISIONAL_CONSTRAINT_DRIFT"
                next_blocker = "validation-interval local constraints exceed the predeclared resolved tolerance"
            else:
                fine_system, fine_bare = build_system(128)
                fine_solved, fine_info = solve_initial_radius(fine_system, fine_bare)
                runs["N128_initial"] = {
                    "converged": fine_info["converged"],
                    "hamilton_max": None if not fine_info["converged"] else constraint_residuals(fine_system, fine_solved)["hamilton_max"],
                    "momentum_max": None if not fine_info["converged"] else constraint_residuals(fine_system, fine_solved)["momentum_max"],
                    "r_min": fine_info.get("r_min"),
                    "r_max": fine_info.get("r_max"),
                }
                if fine_info["converged"]:
                    fine_dt, _fine_omega = stable_timestep(fine_system, fine_solved)
                    fine_run = evolve(fine_system, fine_solved, 0.005, fine_dt)
                    runs["N128_T0.005"] = _compact_run(fine_run)
                    if fine_run["stopped"]:
                        verdict = "STOP_POSITIVE_CHART"
                        next_blocker = fine_run["reason"]
                        passed = False
                    elif fine_run["hamilton_max"] > TOL_CONSTRAINT_VALIDATION or fine_run["momentum_max"] > TOL_CONSTRAINT_VALIDATION:
                        passed = False
                        verdict = "PROVISIONAL_CONSTRAINT_DRIFT"
                        next_blocker = "N=128 local constraints exceed the validation tolerance"
                if passed and verdict != "STOP_POSITIVE_CHART":
                    extended = evolve(system, solved, 0.05, dt)
                    runs["N64_T0.05"] = _compact_run(extended)
                    runs["N64_T0.05"]["end"] = _end_arrays(extended["state"])
                    if extended["stopped"]:
                        verdict = "STOP_POSITIVE_CHART"
                        next_blocker = extended["reason"]
                    elif extended["hamilton_max"] > TOL_CONSTRAINT_VALIDATION or extended["momentum_max"] > TOL_CONSTRAINT_VALIDATION:
                        verdict = "PROVISIONAL_CONSTRAINT_DRIFT"
                        next_blocker = "T=0.05 local constraints exceed the validation tolerance; energy conservation is not treated as a pass"
                    else:
                        verdict = "PROVISIONAL_INITIAL_VALIDATION"
                        next_blocker = "continuum limit of the finite canonical source and independent review"
                elif verdict != "STOP_POSITIVE_CHART":
                    verdict = "PROVISIONAL_CONSTRAINT_DRIFT"
    elif solve_info["converged"]:
        verdict = "FAIL_INITIAL_SOLVE"
        next_blocker = "solved residual above the predeclared initial tolerance"
    record = {
        "schema": SCHEMA,
        "status": STATUS,
        "verdict": verdict,
        "provisional": True,
        "independent_review": False,
        "global_regeneration": False,
        "absolute_vacuum_claim": False,
        "action": "Gamma_one = S_induced[g+] - S_induced[g-] + Gamma_Gaussian(C0), counted once",
        "vacuum_matched_prerequisite": False,
        "sea_subtracted": False,
        "gamma_rest_included": False,
        "V_rel_added_to_potential": False,
        "boundary_coefficients_in_bulk_rhs": False,
        "boundary_coefficients": {"C_E": system.coefficients["C_E"], "C_box": system.coefficients["C_box"]},
        "coefficients": system.coefficients,
        "rmag2": magnetic_radius_square(system.coefficients),
        "multiplicity": system.multiplicity,
        "phi_dot": "-i H Phi",
        "phi_dot_not": "-M i H Phi",
        "static_gauge": ["L", "beta"],
        "calibration": system.calibration,
        "preparation": system.preparation,
        "initial_solve": solve_info,
        "initial_normal_velocity_max": None if normal is None else float(np.max(np.abs(normal))),
        "initial_coordinate_velocity_max": None if coordinate is None else float(np.max(np.abs(coordinate))),
        "initial_shift_velocity_max": None if shift_piece is None else float(np.max(np.abs(shift_piece))),
        "initial_hamilton_max": None if initial is None else initial["hamilton_max"],
        "initial_momentum_max": None if initial is None else initial["momentum_max"],
        "initial_current_max": None if initial is None else initial["current_max"],
        "jacobian_directional_error": jacobian_error,
        "hamiltonian_fd": fd,
        "mutation": mutation,
        "kinematic_gap": kinematic,
        "dense_column_crosscheck_N32": dense,
        "omega_estimate": omega,
        "dt": dt,
        "runs": runs,
        "next_blocker": next_blocker,
        "runtime_seconds": time.perf_counter() - started,
        "cpu_seconds": time.process_time() - cpu,
    }
    if solve_info["converged"]:
        record["reduced_operator"] = reduced_dirac(system, solved)
        record["continuum_packet"] = continuum_packet_diagnostic()
        record["initial_end"] = _end_arrays(solved)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(record), indent=2) + "\n")
    record["path"] = str(path)
    return record
