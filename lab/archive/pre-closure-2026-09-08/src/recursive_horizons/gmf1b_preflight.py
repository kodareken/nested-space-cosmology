"""GMF-1B-PF1 focusing and thermal-spin-fluid characteristic preflight.

The artifact rejects only two specified routes: a canonical scalar plus
Lambda cannot defocus a twist-free future null congruence under the stated
Raychaudhuri assumptions, and the naive thermal averaged EC perfect-fluid
closure is not a healthy causal material law through its EC-1 bounce.  It
does not reject full Einstein--Cartan--Dirac theory, the EC-1 background
algebra, noncanonical actions, modified gravity, or topology-changing
physics.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt
from numbers import Real


def _real(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _positive(name: str, value: Real) -> float:
    result = _real(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _nonnegative(name: str, value: Real) -> float:
    result = _real(name, value)
    if result < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return result


def canonical_scalar_null_energy(k_dot_phi: Real) -> float:
    """Return ``T_kk=(k dot dPhi)^2`` for a canonical scalar plus Lambda."""

    derivative = _real("k_dot_phi", k_dot_phi)
    result = derivative * derivative
    if not isfinite(result):
        raise ValueError("canonical null energy is outside the finite range")
    return result


def canonical_scalar_raychaudhuri_derivative(
    gravitational_constant: Real,
    expansion: Real,
    shear_squared: Real,
    k_dot_phi: Real,
) -> dict[str, float | bool]:
    """Evaluate twist-free metric-null focusing for canonical scalar + Lambda.

    For affine future null generators and vanishing twist,
    ``dtheta/dlambda=-theta^2/2-sigma^2-8*pi*G*Tkk``.  A cosmological
    constant cancels from the null contraction.  This local derivative gate
    is not a global no-go outside its assumptions.
    """

    G = _positive("gravitational_constant", gravitational_constant)
    theta = _real("expansion", expansion)
    sigma_squared = _nonnegative("shear_squared", shear_squared)
    Tkk = canonical_scalar_null_energy(k_dot_phi)
    derivative = -0.5 * theta * theta - sigma_squared - 8.0 * 3.141592653589793 * G * Tkk
    if not isfinite(derivative):
        raise ValueError("Raychaudhuri derivative is outside the finite range")
    return {
        "T_kk": Tkk,
        "raychaudhuri_derivative": derivative,
        "twist_free": True,
        "canonical_null_convergence": True,
        "negative_expansion_can_rise_to_zero_before_caustic_or_endpoint": False,
    }


def canonical_scalar_focusing_gate(
    gravitational_constant: Real,
    initial_expansion: Real,
    shear_squared: Real,
    k_dot_phi: Real,
) -> dict[str, float | bool | str]:
    """Classify the canonical-scalar future-null focusing implication.

    The stated conclusion requires an affine, future-directed, twist-free
    metric-null generator with a canonical scalar, Einstein gravity, and no
    caustic or endpoint before the claimed zero crossing.
    """

    theta = _real("initial_expansion", initial_expansion)
    derivative = canonical_scalar_raychaudhuri_derivative(
        gravitational_constant, theta, shear_squared, k_dot_phi
    )
    return {
        **derivative,
        "initial_expansion": theta,
        "negative_initial_expansion": theta < 0.0,
        "same_generator_defocusing_to_zero_excluded_under_stated_assumptions": theta < 0.0,
        "classification": "conditional_canonical_scalar_null_focusing_gate",
    }


def thermal_ec_state(z: Real) -> dict[str, float | None | str | bool]:
    """Classify the naive thermal averaged EC perfect-fluid closure.

    ``z=alpha*h_n^2*T^2/h_star``.  The returned density and pressure are
    normalized by ``h_star*T^4``.  ``cs2`` is ``None`` at the singular
    enthalpy/derivative point, avoiding JSON NaN/infinity.
    """

    z = _nonnegative("z", z)
    rho_normalized = 1.0 - z
    pressure_normalized = 1.0 / 3.0 - z
    # Threshold-centred forms avoid cancellation to an accidental zero for
    # the representable floats immediately beside z=2/9 and z=2/3.
    denominator = 6.0 * (2.0 / 3.0 - z)
    numerator = 6.0 * (2.0 / 9.0 - z)
    enthalpy_normalized = 2.0 * (2.0 / 3.0 - z)
    if z < 2.0 / 9.0:
        classification = "causal_barotrope"
        cs2: float | None = numerator / denominator
    elif z == 2.0 / 9.0:
        classification = "degenerate_zero_sound_speed"
        cs2 = 0.0
    elif z < 2.0 / 3.0:
        classification = "gradient_unstable"
        cs2 = numerator / denominator
    elif z == 2.0 / 3.0:
        classification = "singular_enthalpy_and_density_derivative"
        cs2 = None
    else:
        classification = "superluminal_relative_to_metric_cone"
        cs2 = numerator / denominator
    derived = (
        rho_normalized,
        pressure_normalized,
        enthalpy_normalized,
        denominator,
        numerator,
    )
    if not all(isfinite(value) for value in derived) or (cs2 is not None and not isfinite(cs2)):
        raise ValueError("thermal EC state is outside the finite range")
    return {
        "z": z,
        "rho_over_hT4": rho_normalized,
        "p_over_hT4": pressure_normalized,
        "rho_plus_p_over_hT4": enthalpy_normalized,
        "density_temperature_derivative_over_hT3": denominator,
        "pressure_temperature_derivative_over_hT3": numerator,
        "cs2": cs2,
        "classification": classification,
        "healthy_causal_barotrope": classification == "causal_barotrope",
    }


@dataclass(frozen=True, slots=True)
class EC1BounceSpec:
    """Positive EC-1 background parameters with two distinct turning roots."""

    A: float
    B: float

    def __post_init__(self) -> None:
        A = _positive("A", self.A)
        B = _positive("B", self.B)
        discriminant = A * A - 4.0 * B
        if not isfinite(discriminant) or discriminant <= 0.0:
            raise ValueError("A**2 must be finite and strictly greater than 4*B")
        object.__setattr__(self, "A", A)
        object.__setattr__(self, "B", B)


def ec1_bounce_thermal_gate(spec: EC1BounceSpec) -> dict[str, float | bool | str | None]:
    """Map an EC-1 bounce to the naive thermal closure's characteristic gate.

    At the lower EC-1 turning point, ``z_bounce=B/(A*x_min)=x_max/A``.
    Since ``A**2>4B``, it lies strictly in ``(1/2,1)``.  It therefore cannot
    lie in the closure's healthy causal interval ``0<cs2<=1``.
    """

    discriminant = spec.A * spec.A - 4.0 * spec.B
    root = sqrt(discriminant)
    x_max = 0.5 * (spec.A + root)
    x_min = spec.B / x_max
    # ``x_max/A`` is mathematically below one, but can round to 1.0 when
    # B/A^2 is below binary64 resolution. Retain its positive complement.
    one_minus_z_bounce = spec.B / (spec.A * x_max)
    if (
        not isfinite(x_min)
        or x_min <= 0.0
        or not isfinite(one_minus_z_bounce)
        or one_minus_z_bounce <= 0.0
    ):
        raise ValueError("EC-1 bounce roots are outside the finite nonzero range")
    z_bounce = 1.0 - one_minus_z_bounce
    thermal = thermal_ec_state(z_bounce)
    return {
        "A": spec.A,
        "B": spec.B,
        "discriminant": discriminant,
        "x_min": x_min,
        "x_max": x_max,
        "z_bounce": z_bounce,
        "one_minus_z_bounce": one_minus_z_bounce,
        "z_bounce_strictly_between_one_half_and_one": True,
        "thermal_classification_at_bounce": thermal["classification"],
        "cs2_at_bounce": thermal["cs2"],
        "healthy_causal_barotrope_at_bounce": thermal["healthy_causal_barotrope"],
        "continuous_low_z_to_bounce_branch_healthy_causal": False,
        "naive_thermal_averaged_perfect_fluid_rejected_as_gmf1b_material_law": True,
        "full_einstein_cartan_dirac_rejected": False,
        "ec1_background_algebra_rejected": False,
    }
