"""Small falsification gates for relative cones and literal boundary scalings.

These helpers test only restricted subclasses: a dimensionless tensor/photon
relative-cone shift, separately conserved power-law boundary densities, and
published DESI DR2 Gaussian summaries.  They are not evidence for a variable
locally measured ``c``, an external origin of dark energy, or recursive domains.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import expm1, isfinite, log1p, sqrt
from numbers import Real
from typing import Final


GW170817_DELTA_T_MIN: Final[float] = -3.0e-15
"""Inclusive GW170817 bound on ``delta_T = c_T/c_gamma - 1``."""

GW170817_DELTA_T_MAX: Final[float] = 7.0e-16
"""Inclusive GW170817 bound on ``delta_T = c_T/c_gamma - 1``."""

APPROXIMATE_95_SIGMA: Final[float] = 1.96


def _finite_real(name: str, value: Real) -> float:
    """Return a finite real value while rejecting bools and coercive inputs."""

    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive(name: str, value: Real) -> float:
    result = _finite_real(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _finite_record(values: dict[str, float | bool | str]) -> dict[str, float | bool | str]:
    """Reject non-finite derived numbers before a record becomes JSON output."""

    if any(isinstance(value, float) and not isfinite(value) for value in values.values()):
        raise ValueError("derived quantity is outside the finite evaluation range")
    return values


def tensor_photon_delta(tensor_speed: Real, photon_speed: Real) -> float:
    """Return the observable dimensionless shift ``delta_T=c_T/c_gamma-1``.

    Both speeds must be positive in the same operational units.  This function
    compares cones; it does not posit a variable locally measured ``c``.
    """

    tensor_speed = _positive("tensor_speed", tensor_speed)
    photon_speed = _positive("photon_speed", photon_speed)
    result = tensor_speed / photon_speed - 1.0
    if not isfinite(result) or result <= -1.0:
        raise ValueError("tensor/photon relative-cone shift is outside the finite range")
    return result


def gw170817_cone_gate(delta_tensor_over_photon: Real) -> dict[str, float | bool | str]:
    """Classify a shift against the inclusive GW170817 tensor/photon bound.

    The gate is an observational constraint subject to the source-emission-delay
    assumptions of GW170817.  Passing it is not evidence for a variable local
    causal speed; failing it rejects this low-redshift relative-cone branch.
    """

    delta = _finite_real("delta_tensor_over_photon", delta_tensor_over_photon)
    if delta <= -1.0:
        raise ValueError("delta_tensor_over_photon must exceed -1")
    return _finite_record({
        "delta_tensor_over_photon": delta,
        "lower_inclusive_bound": GW170817_DELTA_T_MIN,
        "upper_inclusive_bound": GW170817_DELTA_T_MAX,
        "within_inclusive_bound": GW170817_DELTA_T_MIN <= delta <= GW170817_DELTA_T_MAX,
        "classification": "relative_cone_constraint_not_variable_local_c_evidence",
    })


def arrival_delay_fraction(delta_tensor_over_photon: Real) -> float:
    """Return ``(t_T-t_gamma)/t_gamma=-delta/(1+delta)`` for a common path."""

    delta = _finite_real("delta_tensor_over_photon", delta_tensor_over_photon)
    if delta <= -1.0:
        raise ValueError("delta_tensor_over_photon must exceed -1")
    result = -delta / (1.0 + delta)
    if not isfinite(result):
        raise ValueError("arrival-delay fraction is outside the finite range")
    return result


def disformal_photon_cone(
    conformal_factor: Real, disformal_factor: Real, kinetic_x: Real
) -> dict[str, float | str]:
    """Evaluate ``g_gamma=C g+D grad(phi)grad(phi)`` in a homogeneous state.

    With ``X=dot(phi)^2/2 >= 0``, returns
    ``c_gamma/c_T=sqrt(1-2DX/C)`` and the numerically stable reciprocal shift
    ``delta_T=c_T/c_gamma-1``.  The requirements ``C>0`` and ``C-2DX>0`` keep
    the photon effective metric Lorentzian.  This is a relative-cone EFT helper,
    not evidence for a variable local ``c``.
    """

    conformal_factor = _positive("conformal_factor", conformal_factor)
    disformal_factor = _finite_real("disformal_factor", disformal_factor)
    kinetic_x = _finite_real("kinetic_x", kinetic_x)
    if kinetic_x < 0.0:
        raise ValueError("kinetic_x must be non-negative for a homogeneous scalar")
    try:
        correction = 2.0 * disformal_factor * kinetic_x
    except OverflowError as error:
        raise ValueError("disformal correction is outside the finite range") from error
    if not isfinite(correction):
        raise ValueError("disformal correction is outside the finite range")
    lorentzian_time_factor = conformal_factor - correction
    if not isfinite(lorentzian_time_factor) or lorentzian_time_factor <= 0.0:
        raise ValueError("conformal_factor - 2*disformal_factor*kinetic_x must be positive")
    ratio_squared = lorentzian_time_factor / conformal_factor
    if not isfinite(ratio_squared) or ratio_squared <= 0.0:
        raise ValueError("photon/tensor cone ratio is outside the finite range")
    photon_over_tensor = sqrt(ratio_squared)
    # c_T/c_gamma - 1 = exp[-log(1-2DX/C)/2] - 1; expm1 avoids
    # cancellation for the extremely small shifts relevant to GW170817.
    delta_tensor_over_photon = expm1(-0.5 * log1p(-correction / conformal_factor))
    if not isfinite(delta_tensor_over_photon):
        raise ValueError("tensor/photon relative-cone shift is outside the finite range")
    return _finite_record({
        "conformal_factor": conformal_factor,
        "disformal_factor": disformal_factor,
        "kinetic_x": kinetic_x,
        "c_gamma_over_c_tensor": photon_over_tensor,
        "delta_tensor_over_photon": delta_tensor_over_photon,
        "classification": "disformal_relative_cone_not_variable_local_c_evidence",
    })  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class BoundaryScalingSpec:
    """A literal-boundary scaling ansatz, not an external-origin inference.

    It uses ``rho/rho_ref=(L/L_ref)^(-q)`` and
    ``L/L_ref=(a/a_ref)^s``.  Only separately conserved components may be
    assigned the effective equation of state derived from this scaling.
    """

    q: float
    s: float
    separately_conserved: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "q", _finite_real("q", self.q))
        object.__setattr__(self, "s", _finite_real("s", self.s))
        if not isinstance(self.separately_conserved, bool):
            raise ValueError("separately_conserved must be a bool")


def boundary_density_ratio(spec: BoundaryScalingSpec, scale_factor_ratio: Real) -> float:
    """Return ``rho/rho_ref=(a/a_ref)^(-q*s)`` for a positive scale ratio."""

    scale_factor_ratio = _positive("scale_factor_ratio", scale_factor_ratio)
    exponent = _boundary_density_exponent(spec)
    try:
        result = scale_factor_ratio ** (-exponent)
    except OverflowError as error:
        raise ValueError("boundary density ratio is outside the finite range") from error
    if not isfinite(result) or result <= 0.0:
        raise ValueError("boundary density ratio is outside the finite range")
    return result


def _boundary_density_exponent(spec: BoundaryScalingSpec) -> float:
    """Return finite ``q*s`` so derived records cannot contain infinity."""

    exponent = spec.q * spec.s
    if not isfinite(exponent):
        raise ValueError("q*s is outside the finite evaluation range")
    return exponent


def _require_separate_conservation(spec: BoundaryScalingSpec) -> None:
    if not spec.separately_conserved:
        raise ValueError(
            "effective equation of state and acceleration gate are not applicable "
            "without separate conservation"
        )


def boundary_equation_of_state(spec: BoundaryScalingSpec) -> float:
    """Return ``w=-1+q*s/3`` for a separately conserved scaling component."""

    _require_separate_conservation(spec)
    result = -1.0 + _boundary_density_exponent(spec) / 3.0
    if not isfinite(result):
        raise ValueError("boundary equation of state is outside the finite range")
    return result


def boundary_acceleration_gate(spec: BoundaryScalingSpec) -> dict[str, float | bool | str]:
    """Classify the literal scaling by the FLRW acceleration condition ``q*s<2``.

    It rejects only this separately conserved literal subclass; it is not
    evidence for or against an external origin of dark energy.
    """

    _require_separate_conservation(spec)
    exponent = _boundary_density_exponent(spec)
    return _finite_record({
        "density_scale_exponent": exponent,
        "w_effective": boundary_equation_of_state(spec),
        "accelerates": exponent < 2.0,
        "classification": "literal_boundary_scaling_only_not_external_origin_evidence",
    })


@dataclass(frozen=True, slots=True)
class PublishedGaussianSummary:
    """Published mean/sigma summary; deliberately not a full likelihood."""

    label: str
    mean: float
    sigma: float

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label:
            raise ValueError("label must be a non-empty string")
        object.__setattr__(self, "mean", _finite_real("mean", self.mean))
        object.__setattr__(self, "sigma", _positive("sigma", self.sigma))


DESI_DR2_FLAT_WCDM_SUMMARIES: Final[tuple[PublishedGaussianSummary, ...]] = (
    PublishedGaussianSummary("DESI+CMB", -1.055, 0.036),
    PublishedGaussianSummary("DESI+CMB+Pantheon+", -0.995, 0.023),
    PublishedGaussianSummary("DESI+CMB+Union3", -0.997, 0.027),
    PublishedGaussianSummary("DESI+CMB+DESY5", -0.971, 0.021),
)


def approximate_gaussian_w_gate(
    proposed_w: Real, summary: PublishedGaussianSummary, sigma_multiplier: Real = APPROXIMATE_95_SIGMA
) -> dict[str, float | bool | str]:
    """Compare a constant ``w`` with one published Gaussian summary.

    This is an approximate 1.96-sigma summary gate, not a DESI likelihood,
    model comparison, or an inference about a dark-energy origin.
    """

    proposed_w = _finite_real("proposed_w", proposed_w)
    sigma_multiplier = _positive("sigma_multiplier", sigma_multiplier)
    z_score = (proposed_w - summary.mean) / summary.sigma
    return _finite_record({
        "proposed_w": proposed_w,
        "published_mean": summary.mean,
        "published_sigma": summary.sigma,
        "sigma_multiplier": sigma_multiplier,
        "z_score": z_score,
        "within_approximate_interval": abs(z_score) <= sigma_multiplier,
        "classification": "approximate_gaussian_summary_not_likelihood_or_origin_inference",
    })
