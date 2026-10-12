"""Small, dependency-free calculations with explicit epistemic labels.

Nothing in this module solves the proposed transition.  The functions separate
standard cosmological definitions from the conditional horizon-saturation
postulate so an algebraic closure cannot be mistaken for an independent
prediction.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asinh, isfinite, pi, sqrt
from typing import Any

# SI constants.  c is exact; G and hbar use the 2018 CODATA recommended values.
C = 299_792_458.0
G = 6.67430e-11
HBAR = 1.054_571_817e-34
MPC = 3.085_677_581_491_367_3e22
SOLAR_MASS = 1.988_47e30
YEAR = 31_557_600.0
PLANCK_LENGTH = sqrt(HBAR * G / C**3)


def _positive(name: str, value: float) -> float:
    value = float(value)
    if not isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


@dataclass(frozen=True, slots=True)
class CosmologyInputs:
    """Rounded parameters for a transparent flat matter-plus-Lambda example.

    These are not a likelihood object.  The defaults reproduce the rounded
    Planck 2018 base-flat-LambdaCDM values used in the manuscript.
    """

    h0_km_s_mpc: float = 67.4
    omega_m: float = 0.315
    omega_lambda: float = 0.685

    def validate(self) -> None:
        _positive("h0_km_s_mpc", self.h0_km_s_mpc)
        _positive("omega_m", self.omega_m)
        _positive("omega_lambda", self.omega_lambda)
        if abs((self.omega_m + self.omega_lambda) - 1.0) > 1e-10:
            raise ValueError(
                "flat_lcdm_age requires omega_m + omega_lambda = 1 "
                "in the matter-plus-Lambda approximation"
            )


def hubble_parameter_si(h0_km_s_mpc: float) -> float:
    """Convert H0 from km s^-1 Mpc^-1 to s^-1."""

    return _positive("h0_km_s_mpc", h0_km_s_mpc) * 1_000.0 / MPC


def critical_density(hubble_si: float) -> float:
    """Return the FLRW critical mass density, 3 H^2/(8 pi G), in kg m^-3."""

    hubble_si = _positive("hubble_si", hubble_si)
    return 3.0 * hubble_si**2 / (8.0 * pi * G)


def hubble_radius(hubble_si: float) -> float:
    """Return c/H in metres; this is not generally an event-horizon radius."""

    return C / _positive("hubble_si", hubble_si)


def hubble_mass(hubble_si: float) -> float:
    """Critical-density mass inside a sphere of radius c/H, in kilograms."""

    hubble_si = _positive("hubble_si", hubble_si)
    radius = hubble_radius(hubble_si)
    return (4.0 * pi / 3.0) * critical_density(hubble_si) * radius**3


def compactness(mass_kg: float, radius_m: float) -> float:
    """Return the dimensionless spherical ratio 2GM/(c^2 R)."""

    mass_kg = _positive("mass_kg", mass_kg)
    radius_m = _positive("radius_m", radius_m)
    return 2.0 * G * mass_kg / (C**2 * radius_m)


def lambda_from_h0_omega_lambda(
    h0_km_s_mpc: float, omega_lambda: float
) -> float:
    """Return Lambda = 3 H0^2 Omega_Lambda / c^2 in m^-2.

    This is a standard flat-LambdaCDM parameter conversion, not a prediction of
    Lambda by the Recursive Horizons hypothesis.
    """

    omega_lambda = _positive("omega_lambda", omega_lambda)
    hubble_si = hubble_parameter_si(h0_km_s_mpc)
    return 3.0 * hubble_si**2 * omega_lambda / C**2


def de_sitter_radius(lambda_m2: float) -> float:
    """Return sqrt(3/Lambda) in metres for positive Lambda."""

    return sqrt(3.0 / _positive("lambda_m2", lambda_m2))


def horizon_entropy_over_kb(radius_m: float) -> float:
    """Return A/(4 l_P^2) = pi R^2/l_P^2 for a spherical horizon."""

    radius_m = _positive("radius_m", radius_m)
    return pi * radius_m**2 / PLANCK_LENGTH**2


def parent_inference_from_lambda(lambda_m2: float) -> dict[str, float | str]:
    """Infer a spherical parent horizon conditional on the H-SAT postulate.

    The result is explicitly an inverse inference from an input Lambda.  It must
    not be reported as an independent prediction of that same Lambda.
    """

    lambda_m2 = _positive("lambda_m2", lambda_m2)
    radius_m = de_sitter_radius(lambda_m2)
    mass_kg = C**2 * radius_m / (2.0 * G)
    lambda_reconstructed = 3.0 / radius_m**2
    return {
        "classification": "conditional_inference_not_independent_prediction",
        "assumption": "S_child_dS_equals_S_parent_BH_for_spherical_horizons",
        "parent_radius_m": radius_m,
        "parent_mass_kg": mass_kg,
        "parent_mass_solar": mass_kg / SOLAR_MASS,
        "lambda_reconstructed_m^-2": lambda_reconstructed,
        "relative_closure_error": abs(lambda_reconstructed - lambda_m2)
        / lambda_m2,
    }


def flat_lcdm_age(inputs: CosmologyInputs) -> float:
    """Age in seconds for a flat matter-plus-Lambda FLRW approximation.

    Radiation, massive-neutrino details, and likelihood uncertainty are omitted.
    The analytic expression is therefore a transparent benchmark, not a full
    cosmological parameter inference.
    """

    inputs.validate()
    hubble_si = hubble_parameter_si(inputs.h0_km_s_mpc)
    ratio = sqrt(inputs.omega_lambda / inputs.omega_m)
    return 2.0 * asinh(ratio) / (3.0 * hubble_si * sqrt(inputs.omega_lambda))


def reproduce_core(inputs: CosmologyInputs | None = None) -> dict[str, Any]:
    """Return the JSON-serializable identity/conditional-inference record."""

    inputs = inputs or CosmologyInputs()
    inputs.validate()

    hubble_si = hubble_parameter_si(inputs.h0_km_s_mpc)
    rho_critical = critical_density(hubble_si)
    radius_hubble = hubble_radius(hubble_si)
    mass_hubble = hubble_mass(hubble_si)
    hubble_compactness = compactness(mass_hubble, radius_hubble)
    lambda_m2 = lambda_from_h0_omega_lambda(
        inputs.h0_km_s_mpc, inputs.omega_lambda
    )
    radius_ds = de_sitter_radius(lambda_m2)
    age_seconds = flat_lcdm_age(inputs)
    parent = parent_inference_from_lambda(lambda_m2)

    return {
        "schema_version": 1,
        "project_version": "0.11.0",
        "scope": (
            "Standard definitions, exact identities, and a conditional H-SAT "
            "inverse inference; no transition dynamics or empirical detection."
        ),
        "inputs": {
            "classification": "rounded_Planck_2018_base_flat_LCDM_example",
            "h0_km_s_mpc": inputs.h0_km_s_mpc,
            "omega_m": inputs.omega_m,
            "omega_lambda": inputs.omega_lambda,
        },
        "constants_si": {
            "c_m_s": C,
            "G_m3_kg_s2": G,
            "hbar_J_s": HBAR,
            "Mpc_m": MPC,
            "solar_mass_kg": SOLAR_MASS,
            "planck_length_m": PLANCK_LENGTH,
        },
        "standard_cosmology": {
            "classification": "definitions_and_flat_LCDM_approximation",
            "H0_s^-1": hubble_si,
            "critical_density_kg_m^-3": rho_critical,
            "hubble_radius_m": radius_hubble,
            "hubble_mass_kg": mass_hubble,
            "hubble_mass_solar": mass_hubble / SOLAR_MASS,
            "hubble_sphere_compactness": hubble_compactness,
            "lambda_m^-2": lambda_m2,
            "de_sitter_radius_m": radius_ds,
            "de_sitter_entropy_over_kB": horizon_entropy_over_kb(radius_ds),
            "matter_plus_lambda_age_s": age_seconds,
            "matter_plus_lambda_age_Gyr": age_seconds / (1e9 * YEAR),
        },
        "conditional_horizon_saturation": parent,
        "logical_status": {
            "hubble_compactness": "exact_identity_from_shared_definitions",
            "lambda_conversion": "standard_parameter_conversion",
            "parent_mass": "inferred_from_input_lambda_under_H_SAT",
            "lambda_closure": "algebraic_resubstitution_not_prediction",
            "age": "flat_matter_plus_lambda_analytic_benchmark",
        },
    }
