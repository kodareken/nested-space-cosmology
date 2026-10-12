"""Controlled scalar effective-fluid benchmark (DST-1).

DST-1 translates the existing canonical scalar potential
``V = V0 + m^2 varphi^2 / 2`` into elementary homogeneous effective-fluid
identities.  It is not a dark-matter identification, a dark-energy-scale
derivation, reheating microphysics, or a black-hole transition model.

All quantities use ``c = hbar = 1``.  The fixed-H spectator functions are
analytic checks of a damped oscillator; they are not a self-consistent
Friedmann solution.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, exp, isfinite, pi, sin, sqrt
from numbers import Integral, Real
from typing import Final

from .core import C, G


_ZERO_TOLERANCE: Final[float] = 0.0


def _finite_real(name: str, value: Real) -> float:
    """Return a finite real number, rejecting bools and coercive inputs."""

    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _nonnegative(name: str, value: Real) -> float:
    result = _finite_real(name, value)
    if result < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return result


def _positive(name: str, value: Real) -> float:
    result = _finite_real(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _finite_exp(exponent: float, *, name: str) -> float:
    try:
        result = exp(exponent)
    except OverflowError as error:
        raise ValueError(f"{name} is outside the finite evaluation range") from error
    if not isfinite(result) or result == 0.0:
        raise ValueError(f"{name} is outside the finite evaluation range")
    return result


def _equation_of_state(pressure: float, density: float) -> float | None:
    """Return ``p/rho`` or ``None`` for a zero-density component.

    ``None`` is JSON-safe and makes the mathematically undefined zero-over-zero
    case explicit instead of silently emitting NaN or infinity.
    """

    if density == _ZERO_TOLERANCE:
        return None
    return pressure / density


def _finite_mapping(values: dict[str, float | None]) -> dict[str, float | None]:
    """Reject derived overflows so every returned numeric value is JSON-safe."""

    if any(value is not None and not isfinite(value) for value in values.values()):
        raise ValueError("derived quantity is outside the finite evaluation range")
    return values


@dataclass(frozen=True, slots=True)
class DomainScalarSpec:
    """Immutable parameters for ``V=V0+m^2 varphi^2/2``.

    ``vacuum_energy`` has mass dimension four; ``mass`` and ``decay_rate``
    have mass dimension one.  The decay rate is only a phenomenological
    late-time transfer parameter, not an interaction Lagrangian.
    """

    vacuum_energy: float
    mass: float
    decay_rate: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "vacuum_energy", _nonnegative("vacuum_energy", self.vacuum_energy))
        object.__setattr__(self, "mass", _positive("mass", self.mass))
        object.__setattr__(self, "decay_rate", _nonnegative("decay_rate", self.decay_rate))


def instantaneous_split(
    spec: DomainScalarSpec, displacement: Real, velocity: Real
) -> dict[str, float | None]:
    """Return the instantaneous vacuum, oscillatory, and total fluid split.

    A nonzero vacuum offset has ``w_vacuum=-1`` by its stress-tensor form. For
    an identically zero component, the undefined ratio ``p/rho`` is returned
    as ``None`` rather than assigning a continuous-limit convention.
    """

    displacement = _finite_real("displacement", displacement)
    velocity = _finite_real("velocity", velocity)
    kinetic = 0.5 * velocity * velocity
    quadratic = 0.5 * spec.mass * spec.mass * displacement * displacement
    oscillatory_density = kinetic + quadratic
    oscillatory_pressure = kinetic - quadratic
    total_density = spec.vacuum_energy + oscillatory_density
    total_pressure = -spec.vacuum_energy + oscillatory_pressure
    return _finite_mapping({
        "rho_vacuum": spec.vacuum_energy,
        "p_vacuum": -spec.vacuum_energy,
        "w_vacuum": -1.0 if spec.vacuum_energy > 0.0 else None,
        "rho_oscillatory": oscillatory_density,
        "p_oscillatory": oscillatory_pressure,
        "w_oscillatory": _equation_of_state(oscillatory_pressure, oscillatory_density),
        "rho_total": total_density,
        "p_total": total_pressure,
        "w_total": _equation_of_state(total_pressure, total_density),
    })


def quadratic_cycle_average(spec: DomainScalarSpec, amplitude: Real) -> dict[str, float | None]:
    """Return the analytic rapid-quadratic-oscillation average.

    This assumes ``mass >> H`` over an oscillation.  It is an averaging result,
    not a claim that the scalar is observational dark matter.
    """

    amplitude = _finite_real("amplitude", amplitude)
    oscillatory_density = 0.5 * spec.mass * spec.mass * amplitude * amplitude
    total_density = spec.vacuum_energy + oscillatory_density
    total_pressure = -spec.vacuum_energy
    return _finite_mapping({
        "rho_vacuum": spec.vacuum_energy,
        "p_vacuum": -spec.vacuum_energy,
        "w_vacuum": -1.0 if spec.vacuum_energy > 0.0 else None,
        "rho_oscillatory": oscillatory_density,
        "p_oscillatory": 0.0,
        "w_oscillatory": 0.0 if oscillatory_density > 0.0 else None,
        "rho_total": total_density,
        "p_total": total_pressure,
        "w_total": _equation_of_state(total_pressure, total_density),
    })


def monomial_rapid_oscillation(exponent: Real) -> dict[str, float]:
    """Return ``<w>`` and the dilution exponent for ``V proportional |phi|^n``."""

    exponent = _positive("exponent", exponent)
    return _finite_mapping({
        "monomial_exponent": exponent,
        "w_average": (exponent - 2.0) / (exponent + 2.0),
        "density_scale_exponent": 6.0 * exponent / (exponent + 2.0),
    })  # type: ignore[return-value]


def dust_scaled_density(
    initial_density: Real, initial_scale_factor: Real, final_scale_factor: Real
) -> float:
    """Return ``rho_f=rho_i(a_i/a_f)^3`` for an unforced dust component."""

    initial_density = _nonnegative("initial_density", initial_density)
    initial_scale_factor = _positive("initial_scale_factor", initial_scale_factor)
    final_scale_factor = _positive("final_scale_factor", final_scale_factor)
    try:
        result = initial_density * (initial_scale_factor / final_scale_factor) ** 3
    except OverflowError as error:
        raise ValueError("dust-scaled density is outside the finite evaluation range") from error
    if not isfinite(result):
        raise ValueError("dust-scaled density is outside the finite evaluation range")
    return result


def underdamped_spectator_state(
    spec: DomainScalarSpec, amplitude: Real, hubble: Real, time: Real
) -> dict[str, float]:
    """Solve the fixed-H underdamped spectator oscillator exactly.

    It solves ``varphi''+(3H+Gamma)varphi'+m^2 varphi=0`` with
    ``varphi(0)=amplitude`` and ``varphi'(0)=0``, while ``a=exp(H t)``.
    The caller must satisfy ``0 <= 3H+Gamma < 2m``. This endpoint identity is
    not a self-consistent cosmology because the scalar energy does not source
    H. Contracting cases with negative effective damping are outside this
    damped-spectator API.
    """

    amplitude = _finite_real("amplitude", amplitude)
    hubble = _finite_real("hubble", hubble)
    time = _finite_real("time", time)
    damping = 3.0 * hubble + spec.decay_rate
    if damping < 0.0:
        raise ValueError(
            "spectator oscillator requires non-negative effective damping: "
            "3*hubble + decay_rate >= 0"
        )
    alpha = 0.5 * damping
    frequency_squared = spec.mass * spec.mass - alpha * alpha
    if frequency_squared <= 0.0:
        raise ValueError(
            "spectator oscillator must be underdamped: "
            "mass > (3*hubble + decay_rate)/2"
        )
    frequency = sqrt(frequency_squared)
    damping_factor = _finite_exp(-alpha * time, name="damped spectator factor")
    phase = frequency * time
    displacement = amplitude * damping_factor * (
        cos(phase) + alpha * sin(phase) / frequency
    )
    velocity = (
        -amplitude
        * damping_factor
        * (spec.mass * spec.mass / frequency)
        * sin(phase)
    )
    scale_factor = _finite_exp(hubble * time, name="spectator scale factor")
    split = instantaneous_split(spec, displacement, velocity)
    return _finite_mapping({
        "time": time,
        "hubble": hubble,
        "damping": damping,
        "damped_frequency": frequency,
        "scale_factor": scale_factor,
        "displacement": displacement,
        "velocity": velocity,
        "rho_oscillatory": float(split["rho_oscillatory"]),
        "p_oscillatory": float(split["p_oscillatory"]),
    })  # type: ignore[return-value]


def integer_damped_period_benchmark(
    spec: DomainScalarSpec, amplitude: Real, hubble: Real, periods: Integral
) -> dict[str, float]:
    """Check the fixed-H endpoint identity at an integer number of periods.

    At ``t=2*pi*N/omega`` the analytic endpoint obeys
    ``rho_osc a^3 exp(Gamma*t)=rho_osc(0)``.  This is only a fixed-H spectator
    benchmark, not a self-consistent cosmological evolution or reheating proof.
    """

    if isinstance(periods, bool) or not isinstance(periods, Integral):
        raise ValueError("periods must be a non-negative integer")
    periods = int(periods)
    if periods < 0:
        raise ValueError("periods must be non-negative")
    hubble = _finite_real("hubble", hubble)
    amplitude = _finite_real("amplitude", amplitude)
    damping = 3.0 * hubble + spec.decay_rate
    if damping < 0.0:
        raise ValueError(
            "spectator oscillator requires non-negative effective damping: "
            "3*hubble + decay_rate >= 0"
        )
    frequency_squared = spec.mass * spec.mass - 0.25 * damping * damping
    if frequency_squared <= 0.0:
        raise ValueError(
            "spectator oscillator must be underdamped: "
            "mass > (3*hubble + decay_rate)/2"
        )
    frequency = sqrt(frequency_squared)
    time = 2.0 * pi * periods / frequency
    state = underdamped_spectator_state(spec, amplitude, hubble, time)
    initial_density = float(instantaneous_split(spec, amplitude, 0.0)["rho_oscillatory"])
    compensation = _finite_exp(
        damping * time,
        name="volume-and-decay compensation",
    )
    try:
        compensated_density = state["rho_oscillatory"] * compensation
    except OverflowError as error:
        raise ValueError(
            "compensated endpoint density is outside the finite evaluation range"
        ) from error
    if not isfinite(compensated_density):
        raise ValueError("compensated endpoint density is outside the finite evaluation range")
    return _finite_mapping({
        "periods": float(periods),
        "time": time,
        "damped_frequency": frequency,
        "initial_rho_oscillatory": initial_density,
        "compensated_rho_oscillatory": compensated_density,
        "endpoint_residual": compensated_density - initial_density,
    })  # type: ignore[return-value]


def kottler_radial_acceleration(
    mass_kg: Real,
    radius_m: Real,
    lambda_m2: Real,
    gravitational_constant: Real = G,
    causal_speed: Real = C,
) -> float:
    """Return the weak-field Kottler radial acceleration in SI units.

    ``a_r=-GM/r^2 + Lambda*c^2*r/3`` displays competition between compact-mass
    attraction and a positive-Lambda vacuum term. Its zero is an unstable
    radial separator, not a wall, proof that Lambda stops horizon accretion, or
    proof of recursive domains.
    """

    mass_kg = _positive("mass_kg", mass_kg)
    radius_m = _positive("radius_m", radius_m)
    lambda_m2 = _positive("lambda_m2", lambda_m2)
    gravitational_constant = _positive("gravitational_constant", gravitational_constant)
    causal_speed = _positive("causal_speed", causal_speed)
    result = -gravitational_constant * mass_kg / radius_m**2 + lambda_m2 * causal_speed**2 * radius_m / 3.0
    if not isfinite(result):
        raise ValueError("Kottler radial acceleration is outside the finite evaluation range")
    return result


def kottler_balance_radius(
    mass_kg: Real,
    lambda_m2: Real,
    gravitational_constant: Real = G,
    causal_speed: Real = C,
) -> float:
    """Return ``(3GM/(Lambda*c^2))^(1/3)`` in metres for positive inputs.

    This is the unstable weak-field attraction/vacuum-term sign-change radius,
    not a stable physical boundary or a statement about black-hole accretion or
    recursive domains.
    """

    mass_kg = _positive("mass_kg", mass_kg)
    lambda_m2 = _positive("lambda_m2", lambda_m2)
    gravitational_constant = _positive("gravitational_constant", gravitational_constant)
    causal_speed = _positive("causal_speed", causal_speed)
    ratio = 3.0 * gravitational_constant * mass_kg / (lambda_m2 * causal_speed**2)
    if not isfinite(ratio) or ratio <= 0.0:
        raise ValueError("Kottler balance radius is outside the finite evaluation range")
    result = ratio ** (1.0 / 3.0)
    if not isfinite(result) or result <= 0.0:
        raise ValueError("Kottler balance radius is outside the finite evaluation range")
    return result


def exchange_ledger(
    spec: DomainScalarSpec,
    displacement: Real,
    velocity: Real,
    hubble: Real,
    radiation_density: Real,
) -> dict[str, float]:
    """Return the instantaneous ``Q=Gamma*v^2`` exchange and continuity check.

    This is an effective fluid ledger.  A decay rate and transfer law require
    microscopic fields/couplings before they count as reheating microphysics.
    """

    hubble = _finite_real("hubble", hubble)
    radiation_density = _nonnegative("radiation_density", radiation_density)
    split = instantaneous_split(spec, displacement, velocity)
    rho_total = float(split["rho_total"])
    p_total = float(split["p_total"])
    rho_oscillatory = float(split["rho_oscillatory"])
    p_oscillatory = float(split["p_oscillatory"])
    velocity = _finite_real("velocity", velocity)
    transfer = spec.decay_rate * velocity * velocity
    scalar_derivative = -3.0 * hubble * (rho_oscillatory + p_oscillatory) - transfer
    radiation_derivative = -4.0 * hubble * radiation_density + transfer
    total_density = rho_total + radiation_density
    total_pressure = p_total + radiation_density / 3.0
    residual = scalar_derivative + radiation_derivative + 3.0 * hubble * (total_density + total_pressure)
    return _finite_mapping({
        "transfer_q": transfer,
        "scalar_density_derivative": scalar_derivative,
        "radiation_density_derivative": radiation_derivative,
        "total_density": total_density,
        "total_pressure": total_pressure,
        "total_continuity_residual": residual,
    })  # type: ignore[return-value]
