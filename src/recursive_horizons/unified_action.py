"""Minimal owner for the one-shot universal action and its reduced fixed point.

This module does not claim a completed theory.  It owns the finite operator
basis, forbids independent particle masses and per-scale parameter vectors,
and supplies the exact BPS volume sector used by the first particle/wave
identity gate.  Gravity, boundary, gauge-spin, and near-BPS coefficients remain
one shared ``UniversalTheta`` even when a particular reduction does not excite
them.
"""

from __future__ import annotations

from dataclasses import dataclass
import cmath
from fractions import Fraction
from math import isfinite, pi, sqrt
from types import MappingProxyType
from typing import Mapping


Q = Fraction

FIELD_ORDER = ("metric", "connection", "scale", "skyrme", "spinor", "gauge")
OPERATOR_BASIS = (
    "F_chi_times_R",
    "degenerate_beyond_Horndeski_boundary",
    "Skyrme_L0",
    "Skyrme_L2",
    "Skyrme_L4",
    "Skyrme_L6",
    "gauge_spin_kinetic",
    "chi_skyrme_spin_mixing",
    "shared_boundary_kernel",
)
THETA_ORDER = (
    "kappa0",
    "kappa2",
    "kappa4",
    "kappa6",
    "nonminimal_chi",
    "degenerate_boundary",
    "gauge_coupling",
    "spin_coupling",
    "mixing_coupling",
    "boundary_coupling",
)
FORBIDDEN_PARAMETER_FRAGMENTS = (
    "electron_mass",
    "proton_mass",
    "neutron_mass",
    "nuclear_mass",
    "collapse_coupling",
    "child_coupling",
    "gamma_e_fit",
)


def _positive(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer) or answer <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return answer


def _fraction(name: str, value: object, *, nonnegative: bool = False) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an exact integer or Fraction")
    answer = Fraction(value)
    if nonnegative and answer < 0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


@dataclass(frozen=True, slots=True)
class UniversalTheta:
    """One immutable coefficient vector and one dimensional length anchor."""

    length_anchor: float
    coefficients: Mapping[str, Fraction]

    def __post_init__(self) -> None:
        object.__setattr__(self, "length_anchor", _positive("length_anchor", self.length_anchor))
        if not isinstance(self.coefficients, Mapping):
            raise TypeError("coefficients must be a mapping")
        unknown = set(self.coefficients) - set(THETA_ORDER)
        missing = set(THETA_ORDER) - set(self.coefficients)
        if unknown or missing:
            raise ValueError(
                f"UniversalTheta coefficient keys differ: missing={sorted(missing)}, "
                f"unknown={sorted(unknown)}"
            )
        if any(
            fragment in key.lower()
            for key in self.coefficients
            for fragment in FORBIDDEN_PARAMETER_FRAGMENTS
        ):
            raise ValueError("UniversalTheta contains a forbidden hard-coded sector mass")
        normalized = {
            key: _fraction(key, self.coefficients[key], nonnegative=True)
            for key in THETA_ORDER
        }
        if normalized["kappa0"] <= 0 or normalized["kappa6"] <= 0:
            raise ValueError("the stationary BPS core requires positive kappa0 and kappa6")
        object.__setattr__(self, "coefficients", MappingProxyType(normalized))

    def as_tuple(self) -> tuple[Fraction, ...]:
        return tuple(self.coefficients[name] for name in THETA_ORDER)


@dataclass(frozen=True, slots=True)
class StationaryChargeSector:
    charge: int
    stationary_volume: float
    stationary_energy: float
    volume_hessian: float

    def __post_init__(self) -> None:
        if isinstance(self.charge, bool) or not isinstance(self.charge, int) or self.charge <= 0:
            raise ValueError("charge must be a positive integer")
        for name in ("stationary_volume", "stationary_energy", "volume_hessian"):
            object.__setattr__(self, name, _positive(name, getattr(self, name)))


@dataclass(frozen=True, slots=True)
class CollectiveWaveIdentity:
    charge: int
    localized_rest_energy: float
    wave_dispersion_rest_energy: float
    translation_zero_mode: str
    wave_amplitude: str
    detector_coupling: str

    @property
    def same_rest_energy(self) -> bool:
        return self.localized_rest_energy == self.wave_dispersion_rest_energy

    @property
    def no_separate_particle_coordinate(self) -> bool:
        return self.translation_zero_mode == "delta_Phi=-gradient(Phi_star)*delta_X"


@dataclass(frozen=True, slots=True)
class LocalVacuumSpectrum:
    """Hessian data from one local stationary room, before naming constants."""

    temporal_kinetic: float
    spatial_gradient: float
    metric_stiffness: float
    phase_quantum: float
    vacuum_energy: float
    mass_hessian: tuple[float, ...]
    spin_casimir: tuple[float, ...]

    def __post_init__(self) -> None:
        for name in (
            "temporal_kinetic",
            "spatial_gradient",
            "metric_stiffness",
            "phase_quantum",
        ):
            object.__setattr__(self, name, _positive(name, getattr(self, name)))
        if not isfinite(float(self.vacuum_energy)):
            raise ValueError("vacuum_energy must be finite")
        object.__setattr__(self, "vacuum_energy", float(self.vacuum_energy))
        for name in ("mass_hessian", "spin_casimir"):
            values = tuple(float(value) for value in getattr(self, name))
            if not values or any(not isfinite(value) or value < 0.0 for value in values):
                raise ValueError(f"{name} must contain finite nonnegative eigenvalues")
            object.__setattr__(self, name, values)


@dataclass(frozen=True, slots=True)
class LocalConstants:
    """Names assigned to invariant spectral outputs of ``S_one``."""

    characteristic_speed: float
    gravitational_coupling: float
    phase_quantum: float
    planck_length_over_anchor: float
    cosmological_constant_times_anchor_squared: float
    mass_eigenvalues_in_anchor_units: tuple[float, ...]
    spins: tuple[float, ...]
    absolute_zero: float = 0.0


def derive_local_constants(spectrum: LocalVacuumSpectrum) -> LocalConstants:
    if not isinstance(spectrum, LocalVacuumSpectrum):
        raise TypeError("spectrum must be LocalVacuumSpectrum")
    speed = sqrt(spectrum.spatial_gradient / spectrum.temporal_kinetic)
    gravity = 1.0 / (8.0 * pi * spectrum.metric_stiffness)
    planck = sqrt(spectrum.phase_quantum * gravity / speed**3)
    cosmological = spectrum.vacuum_energy / spectrum.metric_stiffness
    masses = tuple(
        sqrt(value / spectrum.temporal_kinetic) for value in spectrum.mass_hessian
    )
    spins = tuple((sqrt(1.0 + 4.0 * value) - 1.0) / 2.0 for value in spectrum.spin_casimir)
    return LocalConstants(
        characteristic_speed=speed,
        gravitational_coupling=gravity,
        phase_quantum=spectrum.phase_quantum,
        planck_length_over_anchor=planck,
        cosmological_constant_times_anchor_squared=cosmological,
        mass_eigenvalues_in_anchor_units=masses,
        spins=spins,
    )


@dataclass(frozen=True, slots=True)
class RecursiveScaleTransfer:
    """Derived parent-to-child unit conversion; no fitted Gamma is stored."""

    omega: float
    parent_hbar_c: float
    child_hbar_c_in_parent_units: float

    def __post_init__(self) -> None:
        for name in ("omega", "parent_hbar_c", "child_hbar_c_in_parent_units"):
            object.__setattr__(self, name, _positive(name, getattr(self, name)))

    @property
    def gamma(self) -> float:
        return (
            self.child_hbar_c_in_parent_units
            / self.parent_hbar_c
            / self.omega
        )


@dataclass(frozen=True, slots=True)
class RecursiveScaleDilation:
    """Unitary scale-and-sheet map for one parent/child spectral step."""

    omega: float
    spatial_dimension: int = 4

    def __post_init__(self) -> None:
        object.__setattr__(self, "omega", _positive("omega", self.omega))
        if (
            isinstance(self.spatial_dimension, bool)
            or not isinstance(self.spatial_dimension, int)
            or self.spatial_dimension <= 0
        ):
            raise ValueError("spatial_dimension must be a positive integer")

    @property
    def wave_amplitude_factor(self) -> float:
        return self.omega ** (self.spatial_dimension / 2.0)

    @property
    def zeta(self) -> float:
        return self.omega**2

    def child_cutoff(self, parent_cutoff: float) -> float:
        return self.omega * _positive("parent_cutoff", parent_cutoff)

    def child_energy(self, parent_energy: float) -> float:
        return self.omega * _positive("parent_energy", parent_energy)

    def child_dimensionless_argument(self, parent_argument: float) -> float:
        value = float(parent_argument)
        if not isfinite(value):
            raise ValueError("parent_argument must be finite")
        return value / self.omega

    def compose(self, other: "RecursiveScaleDilation") -> "RecursiveScaleDilation":
        if not isinstance(other, RecursiveScaleDilation):
            raise TypeError("other must be RecursiveScaleDilation")
        if other.spatial_dimension != self.spatial_dimension:
            raise ValueError("scale dilations must have the same dimension")
        return RecursiveScaleDilation(
            omega=self.omega * other.omega,
            spatial_dimension=self.spatial_dimension,
        )


@dataclass(frozen=True, slots=True)
class ReciprocalGradientEnergy:
    """Derrick-scaled L0/L2/L4/L6 energy with parent/child duality."""

    kappa0: Fraction
    kappa2: Fraction
    kappa4: Fraction
    kappa6: Fraction

    def __post_init__(self) -> None:
        for name in ("kappa0", "kappa2", "kappa4", "kappa6"):
            value = _fraction(name, getattr(self, name), nonnegative=True)
            if value <= 0:
                raise ValueError("reciprocal energy coefficients must be positive")
            object.__setattr__(self, name, value)

    @property
    def self_dual(self) -> bool:
        return self.kappa0 == self.kappa6 and self.kappa2 == self.kappa4

    def energy(self, scale: Fraction) -> Fraction:
        value = _fraction("scale", scale)
        if value <= 0:
            raise ValueError("scale must be positive")
        return (
            self.kappa0 * value**3
            + self.kappa2 * value
            + self.kappa4 / value
            + self.kappa6 / value**3
        )

    def force(self, scale: Fraction) -> Fraction:
        value = _fraction("scale", scale)
        if value <= 0:
            raise ValueError("scale must be positive")
        return -(
            3 * self.kappa0 * value**2
            + self.kappa2
            - self.kappa4 / value**2
            - 3 * self.kappa6 / value**4
        )

    @property
    def fixed_point_hessian(self) -> Fraction:
        return 6 * self.kappa0 + 2 * self.kappa4 + 12 * self.kappa6


@dataclass(frozen=True, slots=True)
class SelfDualScaleCoupling:
    """Local coupling field determined by the gradients it couples."""

    long_gradient_norm_squared: float
    topological_gradient_norm_squared: float
    m0: float
    e0: float

    def __post_init__(self) -> None:
        for name in (
            "long_gradient_norm_squared",
            "topological_gradient_norm_squared",
            "m0",
            "e0",
        ):
            object.__setattr__(self, name, _positive(name, getattr(self, name)))

    @property
    def f(self) -> float:
        return (
            self.topological_gradient_norm_squared
            / (
                self.m0**2
                * self.e0**2
                * self.long_gradient_norm_squared
            )
        ) ** 0.25

    @property
    def long_energy_density(self) -> float:
        return 0.5 * self.m0**2 * self.f**2 * self.long_gradient_norm_squared

    @property
    def topological_energy_density(self) -> float:
        return (
            0.5
            * self.topological_gradient_norm_squared
            / (self.e0**2 * self.f**2)
        )

    @property
    def self_consistent(self) -> bool:
        scale = max(1.0, self.long_energy_density, self.topological_energy_density)
        return (
            abs(self.long_energy_density - self.topological_energy_density)
            <= 64.0 * 2.220446049250313e-16 * scale
        )


@dataclass(frozen=True, slots=True)
class NestedDualPairEnergy:
    """The same dual-balancing law applied to the L06 and L24 pairs."""

    pair_06: float
    pair_24: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "pair_06", _positive("pair_06", self.pair_06))
        object.__setattr__(self, "pair_24", _positive("pair_24", self.pair_24))

    @property
    def scale_state(self) -> float:
        return (self.pair_24 / self.pair_06) ** 0.25

    def energy(self, scale_state: float | None = None) -> float:
        state = self.scale_state if scale_state is None else _positive(
            "scale_state", scale_state
        )
        return 0.5 * (state**2 * self.pair_06 + self.pair_24 / state**2)

    @property
    def pair_contribution(self) -> float:
        return 0.5 * self.scale_state**2 * self.pair_06

    @property
    def self_consistent(self) -> bool:
        other = 0.5 * self.pair_24 / self.scale_state**2
        scale = max(1.0, self.pair_contribution, other)
        return abs(self.pair_contribution - other) <= 64.0 * 2.220446049250313e-16 * scale


@dataclass(frozen=True, slots=True)
class RecursiveOutsideKernel:
    """Exact scalar fixed point obtained by integrating identical nested rooms."""

    boundary_coupling: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "boundary_coupling",
            _positive("boundary_coupling", self.boundary_coupling),
        )

    def effective_inverse(self, local_inverse: float | complex) -> complex:
        """Return the branch continuous with the uncoupled local inverse."""

        value = complex(local_inverse)
        if not isfinite(value.real) or not isfinite(value.imag):
            raise ValueError("local_inverse must be finite")
        discriminant = cmath.sqrt(
            value * value - 4.0 * self.boundary_coupling**2
        )
        first = 0.5 * (value + discriminant)
        second = 0.5 * (value - discriminant)
        return first if abs(first - value) <= abs(second - value) else second

    def fixed_point_residual(self, local_inverse: float | complex) -> complex:
        effective = self.effective_inverse(local_inverse)
        return (
            effective
            - complex(local_inverse)
            + self.boundary_coupling**2 / effective
        )

    def outside_self_energy(self, local_inverse: float | complex) -> complex:
        return self.boundary_coupling**2 / self.effective_inverse(local_inverse)

    def surface_spectral_density(self, spectral_coordinate: float) -> float:
        coordinate = float(spectral_coordinate)
        if not isfinite(coordinate):
            raise ValueError("spectral_coordinate must be finite")
        width = 2.0 * self.boundary_coupling
        if abs(coordinate) >= width:
            return 0.0
        return sqrt(width * width - coordinate * coordinate) / (
            2.0 * pi * self.boundary_coupling**2
        )


@dataclass(frozen=True, slots=True)
class ReciprocalBoundaryLink:
    """Even parent/child link with signed nonlinear transition curvature."""

    baseline: Fraction
    quadratic: Fraction
    quartic: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "baseline",
            _fraction("baseline", self.baseline, nonnegative=True),
        )
        for name in ("quadratic", "quartic"):
            object.__setattr__(
                self,
                name,
                _fraction(name, getattr(self, name)),
            )
        if self.quadratic == 0:
            raise ValueError("the reciprocal boundary link must activate nonlinearly")

    def response(self, gradient_state: Fraction) -> Fraction:
        state = _fraction("gradient_state", gradient_state)
        return (
            self.baseline
            + self.quadratic * state**2
            + self.quartic * state**4
        )

    @property
    def fixed_point_linear_mixing(self) -> Fraction:
        return Fraction(0)

    @property
    def fixed_point_nonlinear_activation(self) -> Fraction:
        return 2 * self.quadratic


@dataclass(frozen=True, slots=True)
class OneAction:
    theta: UniversalTheta

    @property
    def formula(self) -> str:
        return (
            "sum_v int_Mv sqrt(-g)[F_Theta(chi)R/2+L_deg(chi,X,ddchi;Theta)"
            "+L_0246(U;chi,Theta)+L_gauge_spin(Psi,A;chi,U,Theta)"
            "+L_mix(chi,U,Psi;Theta)]"
            "+sum_e int_Sigmae sqrt(|h|)B_one(Phi_minus,Phi_plus,T;Theta)"
        )

    @property
    def quantum_generator(self) -> str:
        return "Z[boundary]=sum_U int Dg DGamma Dchi DPsi DU exp(i*S_one/hbar)"

    def bps_energy(self, charge: int, volume: float) -> float:
        if isinstance(charge, bool) or not isinstance(charge, int) or charge <= 0:
            raise ValueError("charge must be a positive integer")
        size = _positive("volume", volume)
        k0 = float(self.theta.coefficients["kappa0"])
        k6 = float(self.theta.coefficients["kappa6"])
        return k6 * charge * charge / size + k0 * size

    def stationary_charge_sector(self, charge: int) -> StationaryChargeSector:
        if isinstance(charge, bool) or not isinstance(charge, int) or charge <= 0:
            raise ValueError("charge must be a positive integer")
        k0 = float(self.theta.coefficients["kappa0"])
        k6 = float(self.theta.coefficients["kappa6"])
        volume = charge * sqrt(k6 / k0)
        energy = self.bps_energy(charge, volume)
        hessian = 2.0 * k6 * charge * charge / volume**3
        return StationaryChargeSector(charge, volume, energy, hessian)

    def collective_wave_identity(self, charge: int) -> CollectiveWaveIdentity:
        sector = self.stationary_charge_sector(charge)
        return CollectiveWaveIdentity(
            charge=charge,
            localized_rest_energy=sector.stationary_energy,
            wave_dispersion_rest_energy=sector.stationary_energy,
            translation_zero_mode="delta_Phi=-gradient(Phi_star)*delta_X",
            wave_amplitude="psi_q(X,s,phase,t)=<Phi_q_star(X,s,phase)|Psi_universe(t)>",
            detector_coupling="H_int=int_detector J_D(x)*rho_one[Phi_q](x)",
        )

    def reciprocal_pair(self, charge: int, volume: float) -> dict[str, float | bool]:
        sector = self.stationary_charge_sector(charge)
        first = _positive("volume", volume)
        second = sector.stationary_volume**2 / first
        return {
            "first_volume": first,
            "second_volume": second,
            "volume_product": first * second,
            "stationary_volume_squared": sector.stationary_volume**2,
            "first_energy": self.bps_energy(charge, first),
            "second_energy": self.bps_energy(charge, second),
            "energy_invariant": self.bps_energy(charge, first)
            == self.bps_energy(charge, second),
        }


__all__ = [
    "CollectiveWaveIdentity",
    "FIELD_ORDER",
    "FORBIDDEN_PARAMETER_FRAGMENTS",
    "LocalConstants",
    "LocalVacuumSpectrum",
    "NestedDualPairEnergy",
    "OPERATOR_BASIS",
    "OneAction",
    "RecursiveScaleTransfer",
    "RecursiveScaleDilation",
    "ReciprocalGradientEnergy",
    "ReciprocalBoundaryLink",
    "RecursiveOutsideKernel",
    "SelfDualScaleCoupling",
    "StationaryChargeSector",
    "THETA_ORDER",
    "UniversalTheta",
    "derive_local_constants",
]
