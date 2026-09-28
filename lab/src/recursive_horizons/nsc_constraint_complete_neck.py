"""Source-selected constraint-complete NSC neck initial data.

This uses the existing homogeneous Kantowski--Sachs ADM constraint surface.
The completed tensor selects its Landau normal and the areal radius.  No
coupling, state occupation, or scale parameter is varied.
"""
from __future__ import annotations

from math import atanh, isfinite, sqrt

from .nsc_background_projection import ChildFrameTensor


def landau_frame(tensor: ChildFrameTensor) -> dict[str, object]:
    """Diagonalize the radial 1+1 stress by its unique subluminal boost."""
    rho, flux, pressure = tensor.rho, tensor.T01, tensor.p_parallel
    enthalpy = rho+pressure
    discriminant = enthalpy*enthalpy-4.0*flux*flux
    if discriminant <= 0 or flux == 0 or enthalpy >= 0:
        raise ValueError("the recorded branch requires type-I rho+p<0 and nonzero flux")
    # Stable form of the small root of
    # flux*(1+v^2)+v*(rho+pressure)=0.
    velocity = 2.0*flux/(-enthalpy+sqrt(discriminant))
    if not 0 < abs(velocity) < 1:
        raise ArithmeticError("no unique subluminal Landau boost")
    gamma = 1.0/sqrt(1.0-velocity*velocity)
    landau_rho = gamma*gamma*(
        rho+2.0*velocity*flux+velocity*velocity*pressure
    )
    landau_flux = gamma*gamma*(
        flux*(1.0+velocity*velocity)+velocity*enthalpy
    )
    landau_pressure = gamma*gamma*(
        pressure+2.0*velocity*flux+velocity*velocity*rho
    )
    old_null_plus = rho+pressure+2.0*flux
    old_null_minus = rho+pressure-2.0*flux
    new_null_plus = landau_rho+landau_pressure+2.0*landau_flux
    new_null_minus = landau_rho+landau_pressure-2.0*landau_flux
    return {
        "velocity_old_child_frame": velocity,
        "rapidity": atanh(velocity),
        "gamma": gamma,
        "stress": {
            "rho": landau_rho,
            "T01": landau_flux,
            "p_parallel": landau_pressure,
            "p_sphere": tensor.p_sphere,
        },
        "type_I_discriminant": discriminant,
        "residuals": {
            "Landau_flux": landau_flux,
            "trace_invariance": (
                landau_rho-landau_pressure-2.0*tensor.p_sphere
                -(rho-pressure-2.0*tensor.p_sphere)
            ),
            "null_product_invariance": (
                new_null_plus*new_null_minus
                -old_null_plus*old_null_minus
            ),
        },
        "radial_nulls": {
            "old_plus": old_null_plus,
            "old_minus": old_null_minus,
            "Landau_plus": new_null_plus,
            "Landau_minus": new_null_minus,
        },
    }


def source_selected_neck(
    tensor: ChildFrameTensor,
    *,
    einstein_coefficient: float,
    retained_a_parallel: float,
    retained_H_parallel: float,
) -> dict[str, object]:
    """Select a homogeneous minimal-radius Cauchy slice from the tensor."""
    values = (einstein_coefficient, retained_a_parallel, retained_H_parallel)
    if not all(isfinite(value) for value in values):
        raise ValueError("finite source/geometry values required")
    if einstein_coefficient <= 0 or retained_a_parallel <= 0:
        raise ValueError("positive Einstein coefficient and axial scale required")
    frame = landau_frame(tensor)
    stress = frame["stress"]
    radius = sqrt(2.0*einstein_coefficient/stress["rho"])
    H_sphere = 0.0
    hamiltonian = (
        H_sphere*H_sphere+2.0*retained_H_parallel*H_sphere
        +1.0/(radius*radius)-stress["rho"]/(2.0*einstein_coefficient)
    )
    momentum = -stress["T01"]/(2.0*einstein_coefficient)
    return {
        "frame": frame,
        "geometry": {
            "metric_class": "Kantowski-Sachs: dT_L^2-a_parallel^2 dz_L^2-r^2 dOmega_2^2",
            "normal": "future Landau eigenvector of the completed radial stress block",
            "a_parallel": retained_a_parallel,
            "r": radius,
            "H_parallel": retained_H_parallel,
            "H_sphere": H_sphere,
            "is_temporal_minimum_sphere": True,
        },
        "constraints": {
            "Hamiltonian": hamiltonian,
            "momentum": momentum,
        },
    }
