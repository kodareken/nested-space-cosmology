"""Initial-data gate for the locked NSC child backreaction problem.

The Einstein/ADM constraints and Bronnikov neck tensor are imported.  This
module only compares the completed charged CTP tensor with that locked Cauchy
surface.  Time evolution is deliberately not started when a constraint fails.
"""
from __future__ import annotations

from math import isfinite, pi, sqrt

from .nsc_background_projection import ChildFrameTensor


def locked_neck_constraint_gate(
    tensor: ChildFrameTensor,
    *,
    einstein_coefficient: float,
    a_parallel: float,
    sphere_radius: float,
    H_parallel: float,
    H_sphere: float,
) -> dict[str, object]:
    """Evaluate the existing homogeneous ADM constraints at the neck.

    In the stored unit-radius Bronnikov frame the imported two-derivative
    source requirement is rho=2*A and j=0.  The normalized Hamiltonian form is
    H_sphere^2+2 H_parallel H_sphere+1/r^2-rho/(2 A)=0.
    """
    values = (
        einstein_coefficient, a_parallel, sphere_radius,
        H_parallel, H_sphere,
    )
    if not all(isfinite(value) for value in values):
        raise ValueError("finite locked neck data required")
    if min(einstein_coefficient, a_parallel, sphere_radius) <= 0:
        raise ValueError("positive coefficient and metric scales required")

    required_rho = 2.0*einstein_coefficient/(sphere_radius*sphere_radius)
    required_pressure = (
        2.0*einstein_coefficient*(1.0-3.0*pi)
        /(sphere_radius*sphere_radius)
    )
    geometric_hamiltonian = (
        H_sphere*H_sphere+2.0*H_parallel*H_sphere
        +1.0/(sphere_radius*sphere_radius)
    )
    source_hamiltonian = tensor.rho/(2.0*einstein_coefficient)
    hamiltonian_residual = geometric_hamiltonian-source_hamiltonian
    momentum_residual = -tensor.T01/(2.0*einstein_coefficient)

    discriminant = (
        H_parallel*H_parallel
        -1.0/(sphere_radius*sphere_radius)
        +source_hamiltonian
    )
    nearest_constrained_H_sphere = (
        -H_parallel+sqrt(discriminant) if discriminant >= 0 else None
    )
    return {
        "required_two_derivative_neck_tensor": {
            "rho": required_rho,
            "T01": 0.0,
            "p_parallel": required_pressure,
            "p_sphere": required_pressure,
        },
        "component_residual_required_minus_CTP": {
            "rho": required_rho-tensor.rho,
            "T01": -tensor.T01,
            "p_parallel": required_pressure-tensor.p_parallel,
            "p_sphere": required_pressure-tensor.p_sphere,
        },
        "normalized_constraints": {
            "hamiltonian_geometry_minus_source": hamiltonian_residual,
            "momentum_geometry_minus_source": momentum_residual,
            "geometric_hamiltonian": geometric_hamiltonian,
            "source_hamiltonian": source_hamiltonian,
        },
        "diagnostic_repairs_not_applied": {
            "Einstein_coefficient_required_by_density": tensor.rho/2.0,
            "Einstein_coefficient_fractional_change": (
                tensor.rho/(2.0*einstein_coefficient)-1.0
            ),
            "nearest_H_sphere_satisfying_only_Hamiltonian": (
                nearest_constrained_H_sphere
            ),
            "corresponding_r_dot_at_unit_radius": (
                nearest_constrained_H_sphere*sphere_radius
                if nearest_constrained_H_sphere is not None else None
            ),
            "keeps_a_minimum_sphere": bool(
                nearest_constrained_H_sphere is not None
                and abs(nearest_constrained_H_sphere) < 1e-14
            ),
            "momentum_constraint_repaired": False,
        },
        "radial_nulls": {
            "plus": tensor.rho+tensor.p_parallel+2.0*tensor.T01,
            "minus": tensor.rho+tensor.p_parallel-2.0*tensor.T01,
        },
    }
