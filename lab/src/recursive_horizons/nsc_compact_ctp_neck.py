"""Nonlocal CTP covariance of positive compact Dirac levels at the NSC neck.

The fixed free compact reduction gives two four-dimensional Dirac copies at
each m_j=j*pi/(2 L_star).  This module keeps those fields canonical, computes
their massive exterior reflection and transports their Gaussian covariance
through the horizon in the existing Bronnikov geometry.  It returns the
fourth-order superadiabatically subtracted CTP covariance.  The locked
Wilsonian branch restores the corresponding local vacuum response, so the
sum retains the same regulator allocation without an instantaneous-vacuum
subtraction.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from math import atan2, exp, log, pi, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.special import expit
import sympy as sp

from .nsc_unruh_state import ParentDirac, horizon_frame


SIGMA1 = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
SIGMA2 = np.array([[0.0, -1j], [1j, 0.0]], complex)
SIGMA3 = np.diag([1.0, -1.0]).astype(complex)


@dataclass(frozen=True)
class CompactCTPConfig:
    magnetic_flux: int
    omega: float
    cutoff: float
    horizon_rho: float
    surface_gravity: float
    compact_levels: tuple[int, ...] = (1, 2)
    angular_levels: int = 3
    frequency_points: int = 8
    frequency_max: float = 12.0
    evolution_steps: int = 4000
    horizon_offset: float = 1e-10
    scattering_tolerance: float = 3e-9
    outer_floor: float = 60.0

    def __post_init__(self):
        if (isinstance(self.magnetic_flux, bool)
                or not isinstance(self.magnetic_flux, int)
                or self.magnetic_flux == 0):
            raise ValueError("nonzero integer magnetic flux required")
        if not self.compact_levels or any(
            isinstance(level, bool) or not isinstance(level, int) or level < 1
            for level in self.compact_levels
        ):
            raise ValueError("positive compact levels required")
        if len(set(self.compact_levels)) != len(self.compact_levels):
            raise ValueError("compact levels must be unique")
        values = (
            self.omega, self.cutoff, self.horizon_rho, self.surface_gravity,
            self.frequency_max, self.horizon_offset,
            self.scattering_tolerance, self.outer_floor,
        )
        if not all(np.isfinite(values)) or min(
            self.omega, self.cutoff, self.surface_gravity,
            self.frequency_max, self.horizon_offset,
            self.scattering_tolerance, self.outer_floor,
        ) <= 0:
            raise ValueError("finite positive scales required")
        if self.angular_levels < 0 or self.frequency_points < 4:
            raise ValueError("resolved angular/frequency controls required")
        if self.evolution_steps < 100:
            raise ValueError("resolved interior evolution required")


def _su2_step(upper, lower, hx, hy, hz, step):
    """Exact exponential of hx*sigma1+hy*sigma2+hz*sigma3."""
    norm = np.sqrt(hx*hx + hy*hy + hz*hz)
    angle = step * norm
    cosine = np.cos(angle)
    sine_over = np.sinc(angle / pi) * step
    x = sine_over * hx
    y = sine_over * hy
    z = sine_over * hz
    return (
        (cosine - 1j*z) * upper + (-1j*x-y) * lower,
        (-1j*x+y) * upper + (cosine + 1j*z) * lower,
    )


def _frequency_grid(mass, maximum, points):
    edges = sorted(set(float(value) for value in (
        0.0, min(0.25, maximum), min(1.0, maximum),
        min(mass, maximum), min(4.0, maximum),
        min(8.0, maximum), min(12.0, maximum),
        min(16.0, maximum), maximum,
    )))
    edges = [value for index, value in enumerate(edges)
             if index == 0 or value-edges[index-1] > 1e-12]
    nodes, base_weights = leggauss(points)
    frequencies, weights = [], []
    for left, right in zip(edges, edges[1:]):
        frequencies.extend(left + (nodes+1.0)*(right-left)/2.0)
        weights.extend(base_weights*(right-left)/2.0)
    return np.asarray(frequencies), np.asarray(weights), edges


@lru_cache(maxsize=1)
def _adiabatic_bloch_function():
    """Fourth-order superadiabatic negative-energy Bloch vector.

    For compact mass zero, its energy and axial-pressure projections reduce
    exactly to the existing E2+E4 and P2+P4 subtractions.
    """
    q, mass, angular, frequency = sp.symbols(
        "q mass angular frequency", positive=True, real=True
    )
    conformal = (
        3*(sp.pi-q)+sp.Rational(3, 2)*sp.sin(2*q)-sp.sin(q)**2
    )
    h = sp.Matrix([
        -mass/sp.sin(q), angular, -frequency/sp.sqrt(conformal)
    ])
    energy = sp.sqrt(h.dot(h))
    terms = [-h/energy]

    def derivative(vector):
        return vector.applyfunc(
            lambda value: -sp.sqrt(conformal)*sp.diff(value, q)
        )

    for order in range(1, 5):
        previous = derivative(terms[-1])
        perpendicular = -h.cross(previous)/(2*energy**2)
        normalization = sum(
            terms[index].dot(terms[order-index])
            for index in range(1, order)
        )
        terms.append(perpendicular-sp.Rational(1, 2)*normalization*terms[0])
    return sp.lambdify(
        (q, mass, angular, frequency), terms, "numpy", cse=True
    )


def adiabatic_bloch(q, mass, angular, frequency):
    values = _adiabatic_bloch_function()(q, mass, angular, frequency)
    return np.sum([
        np.asarray(value, dtype=float).reshape(3) for value in values
    ], axis=0)


def _massive_reflection(
    background: ParentDirac,
    frequency: float,
    angular: float,
    mass: float,
    *,
    inner_offset: float,
    tolerance: float,
    outer_floor: float,
):
    """Current-normalized scattering for V1*sigma1+V2*sigma2."""
    if frequency > mass:
        momentum = sqrt(frequency*frequency-mass*mass)
        current_kind = "propagating"
    else:
        momentum = 1j*sqrt(max(0.0, mass*mass-frequency*frequency))
        current_kind = "evanescent"
    scale = max(abs(momentum), 0.2)
    end = max(outer_floor, 30.0*max(1.0, angular, mass)/scale)
    metric_a = background.A_from_offset(end-background.horizon_rho)
    radius = sqrt(1.0+end*end)
    v1 = angular*sqrt(metric_a)/radius
    v2 = mass*sqrt(metric_a)
    local_momentum = np.sqrt(complex(frequency*frequency-v1*v1-v2*v2))
    if current_kind == "evanescent" and local_momentum.imag < 0:
        local_momentum = -local_momentum
    if current_kind == "propagating" and local_momentum.real < 0:
        local_momentum = -local_momentum
    denominator = v1-1j*v2
    if abs(denominator) < 1e-14:
        ratio = 0j
    else:
        ratio = (frequency-local_momentum)/denominator
    initial_current = max(0.0, float(1.0-abs(ratio)**2))

    def rhs(value, state):
        offset = exp(value)
        rho = background.horizon_rho+offset
        a = background.A_from_offset(offset)
        r = sqrt(1.0+rho*rho)
        first = angular*sqrt(a)/r
        second = mass*sqrt(a)
        z = state[0]
        dz = offset/a*(
            1j*(first+1j*second)-2j*frequency*z
            +1j*(first-1j*second)*z*z
        )
        dlog = offset/a*(first*z.imag-second*z.real)
        return [dz, dlog]

    run = solve_ivp(
        rhs,
        (log(end-background.horizon_rho), log(inner_offset)),
        [ratio, 0j],
        method="DOP853",
        rtol=tolerance,
        atol=tolerance*0.02,
        max_step=0.18,
    )
    if not run.success:
        raise RuntimeError(run.message)
    raw_ratio = run.y[0, -1]
    log_amplitude = float(run.y[1, -1].real)
    horizon_radius = sqrt(1.0+background.horizon_rho**2)
    distance = sqrt(2.0*inner_offset/background.surface_gravity)/horizon_radius
    effective = sqrt(angular*angular+(mass*horizon_radius)**2)
    phase = atan2(mass*horizon_radius, angular)
    orientation = np.diag([np.exp(-0.5j*phase), np.exp(0.5j*phase)])
    frame = orientation @ horizon_frame(
        frequency, background.surface_gravity, effective, distance,
        background.near_tortoise(inner_offset),
    )
    coefficients = np.linalg.solve(frame, np.array([1.0, raw_ratio]))
    raw_reflection = coefficients[1]/coefficients[0]
    if initial_current > 0:
        transmission = initial_current*exp(-2.0*log_amplitude)/abs(coefficients[0])**2
        transmission = min(1.0, max(0.0, float(transmission)))
        magnitude = sqrt(max(0.0, 1.0-transmission))
        reflection = raw_reflection/abs(raw_reflection)*magnitude
    else:
        transmission = 0.0
        reflection = raw_reflection/abs(raw_reflection)
    return {
        "reflection": complex(reflection),
        "transmission": transmission,
        "current_defect": float(abs(abs(reflection)**2+transmission-1.0)),
        "outer_radius": end,
        "function_evaluations": int(run.nfev),
        "kind": current_kind,
    }


def compact_mode_source(config: CompactCTPConfig, compact_level: int, angular_level: int):
    """State-dependent CTP stress of one compact/angular group at the neck."""
    if compact_level not in config.compact_levels or not 0 <= angular_level <= config.angular_levels:
        raise ValueError("mode is outside the configured sector")
    flux = abs(config.magnetic_flux)
    compact_mass = compact_level*pi/2.0
    angular = sqrt(angular_level*(angular_level+flux))
    degeneracy = flux if angular_level == 0 else 2*(flux+2*angular_level)
    compact_copies = 2
    frequencies, weights, edges = _frequency_grid(
        compact_mass, config.frequency_max, config.frequency_points
    )
    background = ParentDirac(config.horizon_rho, config.surface_gravity)
    scattering = [
        _massive_reflection(
            background, float(frequency), angular, compact_mass,
            inner_offset=config.horizon_offset,
            tolerance=config.scattering_tolerance,
            outer_floor=config.outer_floor,
        )
        for frequency in frequencies
    ]
    reflection = np.asarray([row["reflection"] for row in scattering])
    transmission = np.asarray([row["transmission"] for row in scattering])
    outgoing_occupation = expit(-2*pi*frequencies/config.surface_gravity)
    child_mass = config.omega*compact_mass
    incoming_occupation = expit(
        -2*pi*frequencies/(config.omega*config.surface_gravity)
    )
    incoming_occupation = np.where(frequencies >= child_mass, incoming_occupation, 0.0)

    horizon_radius = sqrt(1.0+background.horizon_rho**2)
    radius2 = horizon_radius*horizon_radius
    signed_offset = (
        -radius2*config.horizon_offset
        + background.horizon_rho*radius2*config.horizon_offset**2
    )
    distance = sqrt(2*config.horizon_offset/config.surface_gravity)/horizon_radius
    tortoise = background.near_tortoise(signed_offset)
    effective = sqrt(angular*angular+(compact_mass*horizon_radius)**2)
    phase = atan2(compact_mass*horizon_radius, angular)
    orientation = np.diag([np.exp(-0.5j*phase), np.exp(0.5j*phase)])
    rotation = np.diag([np.exp(-0.25j*pi), np.exp(0.25j*pi)])
    upper = np.empty(len(frequencies), complex)
    lower = np.empty(len(frequencies), complex)
    incoming_upper = np.empty(len(frequencies), complex)
    incoming_lower = np.empty(len(frequencies), complex)
    covariance_minimum = 1.0
    covariance_maximum = 0.0
    for index, frequency in enumerate(frequencies):
        frame = orientation @ horizon_frame(
            float(frequency), config.surface_gravity, effective,
            distance, tortoise, timelike=True,
        )
        occupied = frame @ np.array([
            sqrt(1.0-outgoing_occupation[index]),
            -1j*sqrt(outgoing_occupation[index])*reflection[index],
        ])
        incoming = frame @ np.array([
            0.0,
            sqrt(incoming_occupation[index]*transmission[index]),
        ])
        covariance = np.outer(occupied, occupied.conj()) + np.outer(incoming, incoming.conj())
        eigenvalues = np.linalg.eigvalsh(covariance)
        covariance_minimum = min(covariance_minimum, float(eigenvalues[0]))
        covariance_maximum = max(covariance_maximum, float(eigenvalues[-1]))
        occupied = rotation @ occupied
        incoming = rotation @ incoming
        upper[index], lower[index] = occupied
        incoming_upper[index], incoming_lower[index] = incoming

    begin = log(config.horizon_offset)
    target = pi/2.0
    end = log(background.horizon_q-target)
    step = (end-begin)/config.evolution_steps
    a = (3.0-2.0*sqrt(3.0))/12.0
    b = (3.0+2.0*sqrt(3.0))/12.0
    c1 = 0.5-sqrt(3.0)/6.0
    c2 = 0.5+sqrt(3.0)/6.0

    def generator(value):
        delta = exp(value)
        q = background.horizon_q-delta
        conformal = background.interior_W(delta)
        radius = 1.0/np.sin(q)
        common = delta/sqrt(conformal)
        return (
            -compact_mass*radius*common,
            angular*common,
            -frequencies*delta/conformal,
        )

    for index in range(config.evolution_steps):
        value = begin+index*step
        x1, y1, z1 = generator(value+c1*step)
        x2, y2, z2 = generator(value+c2*step)
        for hx, hy, hz in (
            (b*x1+a*x2, b*y1+a*y2, b*z1+a*z2),
            (a*x1+b*x2, a*y1+b*y2, a*z1+b*z2),
        ):
            upper, lower = _su2_step(upper, lower, hx, hy, hz, step)
            incoming_upper, incoming_lower = _su2_step(
                incoming_upper, incoming_lower, hx, hy, hz, step
            )

    conformal = sqrt(3*pi/2-1)
    momentum = -frequencies/conformal
    hx = -compact_mass*np.ones_like(momentum)
    hy = angular*np.ones_like(momentum)
    hz = momentum
    energy = np.sqrt(hx*hx+hy*hy+hz*hz)
    differences = []
    for index in range(len(frequencies)):
        incoming_vector = np.array([
            incoming_upper[index], incoming_lower[index]
        ])
        occupied_vector = np.array([upper[index], lower[index]])
        covariance = (
            np.outer(occupied_vector, occupied_vector.conj())
            + np.outer(incoming_vector, incoming_vector.conj())
        )
        bloch = adiabatic_bloch(
            pi/2.0, compact_mass, angular, float(frequencies[index])
        )
        adiabatic = (
            np.eye(2, dtype=complex)
            + bloch[0]*SIGMA1+bloch[1]*SIGMA2+bloch[2]*SIGMA3
        )/2.0
        differences.append(covariance-adiabatic)
    differences = np.asarray(differences)
    trace = lambda matrix: np.einsum("fij,ji->f", differences, matrix).real
    energy_density = np.einsum(
        "fij,fji->f", differences,
        np.asarray([hx[i]*SIGMA1+hy[i]*SIGMA2+hz[i]*SIGMA3 for i in range(len(hx))]),
    ).real
    parallel = momentum*trace(SIGMA3)
    sphere = angular*trace(SIGMA2)/2.0
    factor = compact_copies*degeneracy/(4*pi*pi*conformal)
    rho = float(factor*np.dot(weights, energy_density))
    p_parallel = float(factor*np.dot(weights, parallel))
    p_sphere = float(factor*np.dot(weights, sphere))
    power = float(
        compact_copies*degeneracy/pi
        *np.dot(weights*frequencies, transmission*(outgoing_occupation-incoming_occupation))
    )
    t01 = -power/(4*pi*(3*pi/2-1))
    return {
        "compact_level": compact_level,
        "compact_mass": compact_mass,
        "angular_level": angular_level,
        "angular_mass": angular,
        "angular_degeneracy": int(degeneracy),
        "compact_copies": compact_copies,
        "frequency_edges": edges,
        "frequency_nodes": len(frequencies),
        "rho": rho,
        "T01": t01,
        "p_parallel": p_parallel,
        "p_sphere": p_sphere,
        "radial_null_plus": rho+p_parallel+2*t01,
        "radial_null_minus": rho+p_parallel-2*t01,
        "parent_Killing_power": power,
        "minimum_initial_covariance_eigenvalue": covariance_minimum,
        "maximum_initial_covariance_eigenvalue": covariance_maximum,
        "maximum_scattering_current_defect": max(row["current_defect"] for row in scattering),
        "maximum_outer_radius": max(row["outer_radius"] for row in scattering),
        "total_scattering_function_evaluations": sum(row["function_evaluations"] for row in scattering),
        "incoming_occupation_upper_bound": float(max(incoming_occupation)),
    }


def compact_ctp_source(config: CompactCTPConfig):
    rows = [
        compact_mode_source(config, compact_level, angular_level)
        for compact_level in config.compact_levels
        for angular_level in range(config.angular_levels+1)
    ]
    fields = ("rho", "T01", "p_parallel", "p_sphere",
              "radial_null_plus", "radial_null_minus", "parent_Killing_power")
    totals = {field: float(sum(row[field] for row in rows)) for field in fields}
    totals["both_radial_null_components_negative"] = bool(
        totals["radial_null_plus"] < 0 and totals["radial_null_minus"] < 0
    )
    return {
        "config": {**asdict(config), "compact_levels": list(config.compact_levels)},
        "rows": rows,
        "totals": totals,
        "checks": {
            "minimum_initial_covariance_eigenvalue": min(
                row["minimum_initial_covariance_eigenvalue"] for row in rows
            ),
            "maximum_initial_covariance_eigenvalue": max(
                row["maximum_initial_covariance_eigenvalue"] for row in rows
            ),
            "maximum_scattering_current_defect": max(
                row["maximum_scattering_current_defect"] for row in rows
            ),
        },
    }
