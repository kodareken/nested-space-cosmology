"""Convention-locked algebra for the FGC-1 candidate tournament.

The module implements only quantities that follow directly from the declared
four-dimensional action

``L = (M_pl^2 + beta*phi^2) R / 2
     - (nabla phi)^2 / 2 - V(phi) - (nabla chi)^2 / 2 + f(phi) G``

with ``V = mu^2 phi^2/2 + g4 phi^4/4``.  It deliberately does not implement
a spherical reduction, principal-symbol calculation, or collapse evolution.
Those omissions are machine-readable in every certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from math import isfinite
from numbers import Real
from typing import Any


class ModelID(str, Enum):
    """Frozen entries in the FGC-1 candidate tournament."""

    GR_0 = "GR-0"
    SGB_L = "SGB-L"
    FGC_QR = "FGC-QR"


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


def _finite_output(name: str, value: float) -> float:
    if not isfinite(value):
        raise ValueError(f"{name} is outside the finite range")
    return value


@dataclass(frozen=True, slots=True)
class ActionParameters:
    """Parameters for one exact action branch in natural units.

    Dimensions are ``[M_pl]=[mu]=[phi]=mass``, ``[alpha]=mass^-1``,
    ``[eta]=mass^-2``, ``[L0]=mass^-1``, and ``[g4]=[beta]=1``.
    The tournament uses ``g4>0`` for bounded quartic saturation, ``eta>0``,
    and a nonzero signed ``beta`` in its main FGC-QR entry.  Positivity of
    ``F=M_pl^2+beta*phi^2`` is a solution-level hard boundary rather than a
    parameter-sign shortcut.  These choices are not a nonlinear stability
    result.
    """

    model_id: ModelID
    planck_mass: float
    scalar_mass: float
    quartic_coupling: float
    pulse_width: float
    scalar_field: float = 0.0
    ricci_coupling: float = 0.0
    linear_gb_coupling: float = 0.0
    quadratic_gb_coupling: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.model_id, ModelID):
            raise ValueError("model_id must be a ModelID")
        _positive("planck_mass", self.planck_mass)
        _positive("scalar_mass", self.scalar_mass)
        _positive("quartic_coupling", self.quartic_coupling)
        _positive("pulse_width", self.pulse_width)
        _finite_real("scalar_field", self.scalar_field)
        beta = _finite_real("ricci_coupling", self.ricci_coupling)
        alpha = _finite_real("linear_gb_coupling", self.linear_gb_coupling)
        eta = _finite_real("quadratic_gb_coupling", self.quadratic_gb_coupling)

        if self.model_id is ModelID.GR_0:
            if any(value != 0.0 for value in (beta, alpha, eta)):
                raise ValueError("GR-0 requires beta=alpha=eta=0")
        elif self.model_id is ModelID.SGB_L:
            if beta != 0.0 or eta != 0.0 or alpha == 0.0:
                raise ValueError("SGB-L requires alpha!=0 and beta=eta=0")
        else:
            if alpha != 0.0 or beta == 0.0 or eta <= 0.0:
                raise ValueError("FGC-QR requires beta!=0, eta>0, and alpha=0")


def coupling_functions(
    parameters: ActionParameters, phi: Real | None = None
) -> dict[str, float]:
    """Evaluate ``F``, ``V``, and ``f`` and their first two derivatives."""

    field = _finite_real("phi", parameters.scalar_field if phi is None else phi)
    planck = _positive("planck_mass", parameters.planck_mass)
    mass = _positive("scalar_mass", parameters.scalar_mass)
    quartic = _positive("quartic_coupling", parameters.quartic_coupling)
    beta = float(parameters.ricci_coupling)

    try:
        effective_planck = planck * planck + beta * field * field
        effective_planck_prime = 2.0 * beta * field
        effective_planck_second = 2.0 * beta
        potential = 0.5 * mass * mass * field * field + 0.25 * quartic * field**4
        potential_prime = mass * mass * field + quartic * field**3
        potential_second = mass * mass + 3.0 * quartic * field * field

        if parameters.model_id is ModelID.SGB_L:
            alpha = float(parameters.linear_gb_coupling)
            gb_coupling = alpha * field
            gb_coupling_prime = alpha
            gb_coupling_second = 0.0
        elif parameters.model_id is ModelID.FGC_QR:
            eta = float(parameters.quadratic_gb_coupling)
            gb_coupling = eta * field * field / 8.0
            gb_coupling_prime = eta * field / 4.0
            gb_coupling_second = eta / 4.0
        else:
            gb_coupling = 0.0
            gb_coupling_prime = 0.0
            gb_coupling_second = 0.0
    except OverflowError as exc:
        raise ValueError("coupling functions are outside the finite range") from exc

    values = {
        "F": effective_planck,
        "F_prime": effective_planck_prime,
        "F_second": effective_planck_second,
        "V": potential,
        "V_prime": potential_prime,
        "V_second": potential_second,
        "f_gb": gb_coupling,
        "f_gb_prime": gb_coupling_prime,
        "f_gb_second": gb_coupling_second,
    }
    if not all(isfinite(value) for value in values.values()):
        raise ValueError("coupling functions are outside the finite range")
    return values


def effective_planck_mass_squared(
    parameters: ActionParameters, phi: Real | None = None
) -> float:
    """Return the Jordan-frame coefficient ``F(phi)=M_pl^2+beta*phi^2``."""

    return coupling_functions(parameters, phi)["F"]


def scalar_equation_algebraic_residual(
    parameters: ActionParameters,
    ricci_scalar: Real = 0.0,
    gauss_bonnet: Real = 0.0,
    phi: Real | None = None,
) -> float:
    """Return the scalar equation residual after setting ``Box(phi)=0``.

    The convention-locked scalar equation is
    ``Box(phi)-V'(phi)+F'(phi)R/2+f'(phi)G=0``.
    """

    ricci = _finite_real("ricci_scalar", ricci_scalar)
    invariant = _finite_real("gauss_bonnet", gauss_bonnet)
    functions = coupling_functions(parameters, phi)
    return _finite_output(
        "scalar equation residual",
        -functions["V_prime"]
        + 0.5 * functions["F_prime"] * ricci
        + functions["f_gb_prime"] * invariant,
    )


def linear_effective_mass_squared(
    parameters: ActionParameters,
    ricci_scalar: Real = 0.0,
    gauss_bonnet: Real = 0.0,
    phi: Real | None = None,
) -> float:
    """Return the local linear scalar coefficient fixed by the action.

    It is ``V''-F''R/2-f''G``.  At ``phi=0`` in FGC-QR this becomes
    ``mu^2-beta*R-eta*G/4``.  This coefficient is not a full perturbative
    stability or principal-symbol calculation.
    """

    ricci = _finite_real("ricci_scalar", ricci_scalar)
    invariant = _finite_real("gauss_bonnet", gauss_bonnet)
    functions = coupling_functions(parameters, phi)
    return _finite_output(
        "linear effective mass squared",
        functions["V_second"]
        - 0.5 * functions["F_second"] * ricci
        - functions["f_gb_second"] * invariant,
    )


def schwarzschild_gauss_bonnet(
    schwarzschild_radius: Real, areal_radius: Real
) -> float:
    """Return ``G=12*r_s^2/r^6`` for the Ricci-flat Schwarzschild control."""

    radius_s = _positive("schwarzschild_radius", schwarzschild_radius)
    radius = _positive("areal_radius", areal_radius)
    try:
        result = 12.0 * radius_s * radius_s / radius**6
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError(
            "Schwarzschild Gauss-Bonnet invariant is outside the finite range"
        ) from exc
    result = _finite_output("Schwarzschild Gauss-Bonnet invariant", result)
    if result <= 0.0:
        raise ValueError("Schwarzschild Gauss-Bonnet invariant must be positive")
    return result


def schwarzschild_activation_radius(
    parameters: ActionParameters, schwarzschild_radius: Real
) -> float:
    """Return the FGC-QR zero-mass radius in the Schwarzschild control.

    Solving ``mu^2-eta*G/4=0`` gives
    ``r_*=(3*eta*r_s^2/mu^2)^(1/6)``.  It is a linear-instability scale, not a
    demonstrated collapse transition, EFT cutoff, or bounce radius.
    """

    if parameters.model_id is not ModelID.FGC_QR:
        raise ValueError("Schwarzschild activation radius is defined only for FGC-QR")
    radius_s = _positive("schwarzschild_radius", schwarzschild_radius)
    eta = _positive("quadratic_gb_coupling", parameters.quadratic_gb_coupling)
    mass = _positive("scalar_mass", parameters.scalar_mass)
    try:
        result = (3.0 * eta * radius_s * radius_s / (mass * mass)) ** (1.0 / 6.0)
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("activation radius is outside the finite range") from exc
    result = _finite_output("activation radius", result)
    if result <= 0.0:
        raise ValueError("activation radius must be positive")
    return result


def dimensionless_scan_coordinates(
    parameters: ActionParameters,
    *,
    phi: Real | None = None,
    ricci_scalar: Real = 0.0,
    gauss_bonnet: Real = 0.0,
    areal_radius: Real | None = None,
    schwarzschild_radius: Real | None = None,
) -> dict[str, float]:
    """Return the dimensionless coordinates used by a candidate scan."""

    field = _finite_real("phi", parameters.scalar_field if phi is None else phi)
    ricci = _finite_real("ricci_scalar", ricci_scalar)
    invariant = _finite_real("gauss_bonnet", gauss_bonnet)
    planck = _positive("planck_mass", parameters.planck_mass)
    width = _positive("pulse_width", parameters.pulse_width)
    try:
        output = {
            "phi_over_planck_mass": field / planck,
            "mu_times_L0": float(parameters.scalar_mass) * width,
            "alpha_over_L0": float(parameters.linear_gb_coupling) / width,
            "eta_over_L0_squared": float(parameters.quadratic_gb_coupling)
            / (width * width),
            "g4": float(parameters.quartic_coupling),
            "beta": float(parameters.ricci_coupling),
            "R_times_L0_squared": ricci * width * width,
            "G_times_L0_fourth": invariant * width**4,
            "effective_planck_ratio": effective_planck_mass_squared(parameters, field)
            / (planck * planck),
        }
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("scan coordinates are outside the finite range") from exc
    if areal_radius is not None or schwarzschild_radius is not None:
        if areal_radius is None or schwarzschild_radius is None:
            raise ValueError(
                "areal_radius and schwarzschild_radius must be provided together"
            )
        radius = _positive("areal_radius", areal_radius)
        radius_s = _positive("schwarzschild_radius", schwarzschild_radius)
        output["areal_radius_over_schwarzschild_radius"] = radius / radius_s
        if parameters.model_id is ModelID.FGC_QR:
            output["areal_radius_over_activation_radius"] = (
                radius / schwarzschild_activation_radius(parameters, radius_s)
            )
    if not all(isfinite(value) for value in output.values()):
        raise ValueError("scan coordinates are outside the finite range")
    return output


def required_nonclaims() -> dict[str, bool]:
    """Return every conclusion that this algebraic certificate withholds."""

    return {
        "full_covariant_metric_variation_expanded_and_cross_checked": False,
        "coupled_metric_chi_background_verified": False,
        "spherical_field_equations_derived": False,
        "spherical_principal_symbol_derived": False,
        "nonlinear_strong_hyperbolicity_proven": False,
        "nonlinear_ghost_freedom_proven": False,
        "constraint_propagation_proven": False,
        "metric_null_defocusing_derived": False,
        "finite_collapse_solution_constructed": False,
        "singularity_resolution_proven": False,
        "child_spacetime_constructed": False,
        "shear_robustness_proven": False,
        "global_flux_entropy_closure_derived": False,
        "dark_matter_derived": False,
        "dark_energy_derived": False,
        "variable_speed_of_light_derived": False,
        "particle_spectrum_derived": False,
        "theory_of_everything_derived": False,
        "observational_confirmation_obtained": False,
    }


def action_certificate(
    parameters: ActionParameters,
    *,
    ricci_scalar: Real = 0.0,
    gauss_bonnet: Real = 0.0,
    schwarzschild_radius: Real | None = None,
    areal_radius: Real | None = None,
) -> dict[str, Any]:
    """Build a JSON-safe action/background certificate with explicit gates."""

    ricci = _finite_real("ricci_scalar", ricci_scalar)
    invariant = _finite_real("gauss_bonnet", gauss_bonnet)
    if schwarzschild_radius is not None or areal_radius is not None:
        if schwarzschild_radius is None or areal_radius is None:
            raise ValueError(
                "schwarzschild_radius and areal_radius must be provided together"
            )
        invariant = schwarzschild_gauss_bonnet(schwarzschild_radius, areal_radius)

    functions = coupling_functions(parameters)
    planck_squared = functions["F"]
    if planck_squared <= 0.0:
        raise ValueError("effective Planck mass squared must be positive")
    mass_squared = linear_effective_mass_squared(parameters, ricci, invariant)
    residual = scalar_equation_algebraic_residual(parameters, ricci, invariant)
    scan = dimensionless_scan_coordinates(
        parameters,
        ricci_scalar=ricci,
        gauss_bonnet=invariant,
        areal_radius=areal_radius,
        schwarzschild_radius=schwarzschild_radius,
    )

    certificate: dict[str, Any] = {
        "schema_version": 1,
        "model_id": parameters.model_id.value,
        "classification": "covariant_action_and_linear_background_certificate_not_collapse_solution",
        "convention": {
            "metric_signature": "-+++",
            "units": "natural_units_c_hbar_1",
            "gauss_bonnet": "G=R^2-4*R_ab*R^ab+R_abcd*R^abcd",
            "action_density": "(M_Pl^2+beta*phi^2)*R/2-(nabla_phi)^2/2-V(phi)-(nabla_chi)^2/2+f(phi)*G",
            "potential": "V(phi)=mu^2*phi^2/2+g4*phi^4/4",
            "scalar_equation": "Box(phi)-V_prime+beta*phi*R+f_prime*G=0",
            "linear_effective_mass_squared": "V_second-beta*R-f_second*G",
        },
        "parameters": {
            "planck_mass": float(parameters.planck_mass),
            "scalar_mass_mu": float(parameters.scalar_mass),
            "quartic_coupling_g4": float(parameters.quartic_coupling),
            "pulse_width_L0": float(parameters.pulse_width),
            "reference_scalar_field": float(parameters.scalar_field),
            "ricci_coupling_beta": float(parameters.ricci_coupling),
            "linear_gb_coupling_alpha": float(parameters.linear_gb_coupling),
            "quadratic_gb_coupling_eta": float(parameters.quadratic_gb_coupling),
        },
        "couplings_at_reference_field": functions,
        "background": {
            "ricci_scalar": ricci,
            "gauss_bonnet": invariant,
            "scalar_equation_algebraic_residual": residual,
            "linear_effective_mass_squared": mass_squared,
            "linear_tachyonic_at_reference_point": mass_squared < 0.0,
        },
        "dimensionless_scan": scan,
        "action_gate": {
            "candidate_covariant_action_declared": True,
            "scalar_equations_derived": True,
            "metric_equation_variationally_defined": True,
            "separate_matter_and_regulator_fields": True,
        },
        "low_gradient_reference_gate": {
            "effective_planck_mass_squared": planck_squared,
            "effective_planck_ratio": planck_squared
            / (float(parameters.planck_mass) ** 2),
            "effective_planck_mass_positive": planck_squared > 0.0,
            "zero_regulator_reference": float(parameters.scalar_field) == 0.0,
            "zero_regulator_solves_declared_scalar_equation": residual == 0.0,
            "coupled_metric_chi_background_verified": False,
            "exact_gr0_action": parameters.model_id is ModelID.GR_0,
            "solution_level_gr_recovery_demonstrated": False,
        },
        "weak_field_kinetic_reference": {
            "canonical_phi_kinetic_coefficient": 1.0,
            "canonical_chi_kinetic_coefficient": 1.0,
            "positive_einstein_coefficient_at_reference": planck_squared > 0.0,
            "full_coupled_quadratic_action_diagonalized": False,
            "nonlinear_ghost_freedom_proven": False,
        },
        "open_gates": [
            "expand and independently cross-check the full metric variation",
            "derive the spherical equations and constraint propagation system",
            "derive the unredefined nonlinear principal symbol and hyperbolicity domain",
            "compute the physical affine Raychaudhuri sign and magnitude",
            "construct convergent independent-matter collapse inside EFT control",
            "close flux, entropy, shear, and any global continuation separately",
        ],
        "nonclaims": required_nonclaims(),
    }
    if parameters.model_id is ModelID.FGC_QR and schwarzschild_radius is not None:
        activation = schwarzschild_activation_radius(parameters, schwarzschild_radius)
        activation_invariant = schwarzschild_gauss_bonnet(
            schwarzschild_radius, activation
        )
        activation_mass_squared = linear_effective_mass_squared(
            parameters, 0.0, activation_invariant, phi=0.0
        )
        certificate["schwarzschild_linear_activation_control"] = {
            "activation_radius": activation,
            "gauss_bonnet_at_activation": activation_invariant,
            "linear_effective_mass_squared_at_activation": activation_mass_squared,
            "sample_areal_radius": float(areal_radius),
            "sample_is_inside_activation_radius": float(areal_radius) < activation,
            "sample_is_on_activation_radius": float(areal_radius) == activation,
            "sample_is_outside_activation_radius": float(areal_radius) > activation,
            "scaling": "r_star proportional to r_s^(1/3)*eta^(1/6)*mu^(-1/3)",
            "interpretation": "linear_instability_scale_not_transition_or_eft_cutoff",
        }
    return certificate
