"""Pre-step gate for coupled charged-CTP/metric evolution.

The source-selected neck changes the Dirac Hamiltonian.  A stored integrated
stress tensor is not a Gaussian covariance and cannot be evolved by the CTP
equation.  This owner quantifies that initial metric vertex and refuses to
manufacture a fluid closure from four stress moments.
"""
from __future__ import annotations

from math import isfinite, log


def source_revaluation_preflight(
    *,
    seed_radius: float,
    selected_radius: float,
    rho_landau: float,
    p_sphere_landau: float,
) -> dict[str, float]:
    values = (seed_radius, selected_radius, rho_landau, p_sphere_landau)
    if not all(isfinite(value) for value in values):
        raise ValueError("finite source and geometry data required")
    if seed_radius <= 0 or selected_radius <= 0:
        raise ValueError("positive radii required")
    log_radius_change = log(selected_radius/seed_radius)
    density_vertex = -2.0*(rho_landau+p_sphere_landau)
    return {
        "seed_radius": seed_radius,
        "selected_radius": selected_radius,
        "log_radius_change": log_radius_change,
        "angular_Hamiltonian_coefficient_ratio": seed_radius/selected_radius,
        "angular_Hamiltonian_fractional_change": seed_radius/selected_radius-1.0,
        "fixed_covariance_density_log_radius_vertex": density_vertex,
        "linearized_density_change_diagnostic": density_vertex*log_radius_change,
    }
