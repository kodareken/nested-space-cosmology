"""NGS-0 finite-scale, recurrence, and eternal-opportunity controls.

These functions make three pieces of the Nested Gradient Spectrum discussion
executable: the preferred scale of the ``-alpha*k^2 + beta*k^4`` control, the
positive eigenvalue of a two-level linear recurrence, and the probability of
at least one occurrence in a declared finite independent-trial ensemble.

They are mathematical controls.  They do not construct nested spacetimes,
derive a covariant restoring action, identify either dark sector, establish a
cosmological golden ratio, guarantee observers, or turn a metaphysical reading
into physical evidence.
"""

from __future__ import annotations

from math import exp, expm1, hypot, isfinite, log1p, sqrt
from numbers import Integral, Real


def _finite_real(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _positive(name: str, value: Real) -> float:
    result = _finite_real(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def finite_band_selection(alpha: Real, beta: Real) -> dict[str, float | bool | str]:
    """Return the nonzero scale selected by ``-alpha*k^2 + beta*k^4``.

    The stationary nonzero wavenumber is ``k_star^2=alpha/(2*beta)`` and
    ``ell_star=1/k_star``.  This is a spatial-energy control only; it makes no
    claim about a relativistic principal symbol or higher-time derivatives.
    """

    alpha_value = _positive("alpha", alpha)
    beta_value = _positive("beta", beta)
    k_squared = 0.5 * (alpha_value / beta_value)
    if not isfinite(k_squared) or k_squared <= 0.0:
        raise ValueError("selected wavenumber is outside the finite positive range")
    k_star = sqrt(k_squared)
    ell_star = 1.0 / k_star
    stationary_residual = -2.0 * alpha_value * k_star + 4.0 * beta_value * k_star**3
    energy_coefficient = -alpha_value * k_squared + beta_value * k_squared**2
    if not all(
        isfinite(value)
        for value in (k_star, ell_star, stationary_residual, energy_coefficient)
    ):
        raise ValueError("finite-band control is outside the finite range")
    return {
        "alpha": alpha_value,
        "beta": beta_value,
        "k_star_squared": k_squared,
        "k_star": k_star,
        "ell_star": ell_star,
        "stationary_derivative_residual": stationary_residual,
        "energy_coefficient_at_k_star": energy_coefficient,
        "nonzero_finite_scale_selected": True,
        "classification": "spatial_scale_selection_control_not_covariant_gravity",
    }


def two_level_recurrence_ratio(a: Real, b: Real) -> dict[str, float | bool | str]:
    """Return the positive fixed ratio for ``L[n+1]=a*L[n]+b*L[n-1]``.

    The positive root satisfies ``q^2=a*q+b``.  ``a=b=1`` yields the golden
    ratio, but no coefficients are privileged by this calculation.
    """

    a_value = _finite_real("a", a)
    b_value = _positive("b", b)
    if a_value < 0.0:
        raise ValueError("a must be non-negative for this positive-scale control")
    sqrt_b = sqrt(b_value)
    discriminant_root = hypot(a_value, 2.0 * sqrt_b)
    q = 0.5 * a_value + 0.5 * discriminant_root
    scale = max(a_value, sqrt_b)
    scaled_q = q / scale
    scaled_a = a_value / scale
    scaled_sqrt_b = sqrt_b / scale
    scaled_residual = (
        scaled_q * scaled_q
        - scaled_a * scaled_q
        - scaled_sqrt_b * scaled_sqrt_b
    )
    golden = a_value == 1.0 and b_value == 1.0
    if not isfinite(q) or q <= 0.0 or not isfinite(scaled_residual):
        raise ValueError("recurrence ratio is outside the finite positive range")
    return {
        "a": a_value,
        "b": b_value,
        "positive_ratio": q,
        "scaled_characteristic_residual": scaled_residual,
        "golden_ratio_fixture": golden,
        "coefficients_derived_from_physics": False,
        "classification": "conditional_recurrence_eigenvalue_not_cosmological_prediction",
    }


def eventual_occurrence_probability(
    probability_per_trial: Real,
    trial_count: Integral,
) -> dict[str, float | int | bool | str]:
    """Return ``1-(1-p)^N`` with stable arithmetic and explicit assumptions.

    The finite calculation assumes independent trials with one fixed ``p``.
    Its ``N -> infinity`` limit is one only for ``p>0``.  This is not a proof
    of independence, ergodicity, a domain measure, or logical necessity.
    """

    probability = _finite_real("probability_per_trial", probability_per_trial)
    if probability < 0.0 or probability > 1.0:
        raise ValueError("probability_per_trial must lie in the closed interval [0, 1]")
    if isinstance(trial_count, bool) or not isinstance(trial_count, Integral):
        raise ValueError("trial_count must be a non-negative integer")
    trials = int(trial_count)
    if trials < 0:
        raise ValueError("trial_count must be a non-negative integer")
    if trials > 10**15:
        raise ValueError("trial_count exceeds the deterministic control range")

    if trials == 0 or probability == 0.0:
        at_least_one = 0.0
        none = 1.0
    elif probability == 1.0:
        at_least_one = 1.0
        none = 0.0
    else:
        log_none = trials * log1p(-probability)
        at_least_one = -expm1(log_none)
        none = exp(log_none)
    normalization_residual = at_least_one + none - 1.0
    if not all(isfinite(value) for value in (at_least_one, none, normalization_residual)):
        raise ValueError("occurrence probability is outside the finite range")
    return {
        "probability_per_trial": probability,
        "trial_count": trials,
        "probability_at_least_one": at_least_one,
        "probability_none": none,
        "normalization_residual": normalization_residual,
        "infinite_independent_trial_limit": 1.0 if probability > 0.0 else 0.0,
        "fixed_probability_assumed": True,
        "independence_assumed": True,
        "ergodicity_derived": False,
        "logical_necessity_proven": False,
        "classification": "conditional_independent_trial_control_not_inevitable_life",
    }


def required_nonclaims() -> dict[str, bool]:
    """Return the scientific conclusions NGS-0 explicitly does not supply."""

    return {
        "covariant_fgc_action_derived": False,
        "nested_spacetime_constructed": False,
        "finite_cross_domain_spectrum_derived": False,
        "dark_matter_identified": False,
        "dark_energy_identified": False,
        "universal_hertz_increment_derived": False,
        "golden_ratio_cosmology_derived": False,
        "observers_guaranteed": False,
        "metaphysical_plan_proven": False,
    }
