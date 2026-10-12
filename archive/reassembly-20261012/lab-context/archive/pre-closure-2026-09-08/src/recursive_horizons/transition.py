"""Controlled closed-FLRW transition sector for Recursive Horizons.

This is an exact vacuum toy sector: Einstein gravity plus a canonical scalar
at a positive-potential stationary point. Its closed slicing is global de Sitter space,
``a(t) = L cosh(t/L)``. It is not a black-hole exterior/interior embedding,
does not derive H-SAT, and does not establish a child-universe transition.

Conventions: signature (-,+,+,+), c = hbar = 1, M_pl^2 = (8 pi G)^-1.
"""

from __future__ import annotations

from math import atan2, cosh, isfinite, pi, tan, tanh
from typing import Literal


SurfaceRegion = Literal["trapped", "anti_trapped", "normal", "marginal"]


def _finite(name: str, value: float) -> float:
    value = float(value)
    if not isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _positive(name: str, value: float) -> float:
    value = _finite(name, value)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")
    return value


def _chi(chi: float) -> float:
    chi = _finite("chi", chi)
    if not 0.0 < chi < pi:
        raise ValueError("chi must satisfy 0 < chi < pi")
    return chi


def scale_factor(time: float, length: float) -> float:
    """Return the exact scale factor ``a(t) = L cosh(t/L)``."""

    time = _finite("time", time)
    length = _positive("length", length)
    try:
        value = length * cosh(time / length)
    except OverflowError as error:
        raise ValueError("time / length is outside the finite evaluation range") from error
    if not isfinite(value):
        raise ValueError("time / length is outside the finite evaluation range")
    return value


def hubble_parameter(time: float, length: float) -> float:
    """Return ``H(t) = tanh(t/L)/L``."""

    time = _finite("time", time)
    length = _positive("length", length)
    return tanh(time / length) / length


def hubble_derivative(time: float, length: float) -> float:
    """Return ``dot(H)=sech(t/L)^2/L^2`` without overflow-prone cosh."""

    time = _finite("time", time)
    length = _positive("length", length)
    h_times_l = tanh(time / length)
    return (1.0 - h_times_l * h_times_l) / length**2


def potential_at_stationary_point(
    length: float, reduced_planck_mass: float = 1.0
) -> float:
    """Return the constant stationary-point potential ``V0 = 3 M_pl^2/L^2``."""

    length = _positive("length", length)
    reduced_planck_mass = _positive("reduced_planck_mass", reduced_planck_mass)
    return 3.0 * reduced_planck_mass**2 / length**2


def friedmann_residual(
    time: float, length: float, reduced_planck_mass: float = 1.0
) -> float:
    """Return ``H^2 + 1/a^2 - V0/(3 M_pl^2)`` for the exact solution."""

    reduced_planck_mass = _positive("reduced_planck_mass", reduced_planck_mass)
    a = scale_factor(time, length)
    h = hubble_parameter(time, length)
    v0 = potential_at_stationary_point(length, reduced_planck_mass)
    return h * h + 1.0 / (a * a) - v0 / (3.0 * reduced_planck_mass**2)


def curvature_invariants(length: float) -> dict[str, float]:
    """Return the constant 4D de Sitter invariants for curvature radius ``L``."""

    length = _positive("length", length)
    return {
        "ricci_scalar": 12.0 / length**2,
        "ricci_tensor_squared": 36.0 / length**4,
        "riemann_tensor_squared": 24.0 / length**4,
    }


def future_null_expansions(time: float, chi: float, length: float) -> tuple[float, float]:
    """Return ``(theta_plus, theta_minus)`` for future ``k_±=d_t±a^-1d_chi``.

    With areal radius ``R_A=a sin(chi)``, the result is
    ``theta_± = 2[H ± cot(chi)/a]`` and ``k_+ . k_- = -2``.
    """

    chi = _chi(chi)
    a = scale_factor(time, length)
    h = hubble_parameter(time, length)
    cotangent = 1.0 / tan(chi)
    return (2.0 * (h + cotangent / a), 2.0 * (h - cotangent / a))


def classify_surface(
    time: float, chi: float, length: float, *, tolerance: float = 1e-12
) -> SurfaceRegion:
    """Classify a sphere from signs of the two future null expansions."""

    tolerance = _finite("tolerance", tolerance)
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    theta_plus, theta_minus = future_null_expansions(time, chi, length)
    if abs(theta_plus) <= tolerance or abs(theta_minus) <= tolerance:
        return "marginal"
    if theta_plus > 0.0 and theta_minus > 0.0:
        return "anti_trapped"
    if theta_plus < 0.0 and theta_minus < 0.0:
        return "trapped"
    return "normal"


def apparent_horizon_chi(time: float, length: float) -> dict[str, float]:
    """Return locations where respectively ``theta_plus`` and ``theta_minus`` vanish."""

    a_h = scale_factor(time, length) * hubble_parameter(time, length)
    return {
        "theta_plus_zero": atan2(1.0, -a_h),
        "theta_minus_zero": atan2(1.0, a_h),
    }


def closed_curvature_discriminator(
    scale_factor_value: float, hubble_value: float, causal_speed: float = 1.0
) -> float:
    """Return ``Omega_K=-(c/(aH))^2`` for the strict closed branch.

    It is undefined at the exact bounce, so this helper requires nonzero H.
    """

    scale_factor_value = _positive("scale_factor_value", scale_factor_value)
    hubble_value = _finite("hubble_value", hubble_value)
    if hubble_value == 0.0:
        raise ValueError("hubble_value must be nonzero")
    causal_speed = _positive("causal_speed", causal_speed)
    return -(causal_speed / (scale_factor_value * hubble_value)) ** 2
