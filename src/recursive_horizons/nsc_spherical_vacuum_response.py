"""Constrained linear vacuum response of the owned spherical auxiliary action.

Physical areal gauge: ds^2=N^2 dt^2-a^2 dR^2-R^2 dOmega^2, N=a=1
in the reference. Maxwell and Dirac states vanish. The zero-mass perturbation
has c=delta(chi)=R psi; lapse and radial metric traces are constrained by
psi. This is one admissible response column, not three free Dirichlet inputs.
All face coefficients use the minimal unreferenced boundary representative.
"""
from __future__ import annotations

import cmath
from math import isfinite, pi, sqrt

import sympy as sp

from .nsc_spherical_feedback_action import alpha_of, require_weyl_chart


def vacuum_coefficients(A, C_W):
    """Coefficients derived by eliminating the linear lapse/shift constraints.

    S2_bulk=K int [psi_t^2-psi_R^2-(6/R^2+m2) psi^2] dt dR.
    The sign of K refers to this owned action, before any EFT order reduction.
    """
    require_weyl_chart(C_W, A)
    symbolic = isinstance(A, sp.Expr) or isinstance(C_W, sp.Expr)
    if symbolic:
        A, C_W = sp.sympify(A), sp.sympify(C_W)
        if A.is_positive is not True:
            raise ValueError("vacuum reference requires explicitly positive A")
    elif not (isfinite(float(A)) and float(A) > 0 and isfinite(float(C_W))):
        raise ValueError("vacuum reference requires finite A>0 and nonzero real C_W")
    alpha = alpha_of(C_W)
    circle_pi = sp.pi if symbolic else pi
    k = 4 * circle_pi * A / alpha
    return {"alpha": alpha, "constraint_k": k,
            "mass_squared": -k / 6, "quadratic_kinetic_coefficient": -6 * alpha / k}


def _laplace_domain(laplace_s, radius, A, C_W):
    coefficients = vacuum_coefficients(float(A), float(C_W))
    s = complex(laplace_s)
    R = float(radius)
    if not (isfinite(R) and R > 0 and isfinite(s.real) and isfinite(s.imag)):
        raise ValueError("finite positive radius and finite Laplace parameter required")
    growth = sqrt(max(0.0, -coefficients["mass_squared"]))
    if s.real <= growth:
        raise ValueError("retarded transform requires Re(s) greater than the vacuum growth bound")
    q = cmath.sqrt(s*s + coefficients["mass_squared"])
    if q.real <= 0:
        raise ValueError("no decaying radial branch in the admitted Laplace half-plane")
    return coefficients, s, R, q, growth


def vacuum_curvature_profile(laplace_s, radius, boundary_radius, A, C_W):
    """c(R)/c(R0), zero initial/mass data, decaying at infinity for Re(s)>growth."""
    _coefficients, _s, R0, q, _growth = _laplace_domain(laplace_s, boundary_radius, A, C_W)
    R = float(radius)
    if not (isfinite(R) and R >= R0):
        raise ValueError("profile radius must lie in the exterior R>=R0")
    polynomial = lambda z: z*z + 3*z + 3
    return cmath.exp(-q*(R-R0)) * polynomial(q*R)/polynomial(q*R0) * R0/R


def vacuum_exterior_response(laplace_s, boundary_radius, A, C_W, *, orientation=-1):
    """Constrained metric and face responses per unit boundary c in Laplace space.

    The exterior is R>=R0, so its physical outward normal at R0 has s=-1.
    R derivatives point into the exterior; the orientation is applied only
    to outward normal derivatives and face responses. N approaches 1 at
    infinity and the ADM mass perturbation is zero. The time transform is
    causal with compatible zero initial data. No intrinsic face term is added.
    """
    if orientation not in (-1, 1):
        raise ValueError("orientation must be +1 or -1")
    coefficients, s, R, q, growth = _laplace_domain(laplace_s, boundary_radius, A, C_W)
    alpha, k = coefficients["alpha"], coefficients["constraint_k"]
    z = q*R
    polynomial = z*z + 3*z + 3
    dpsi = -q - 3*(z+2)/(R*polynomial)
    dc = 1/R + dpsi
    b = -(dpsi+2/R)/(k*R)
    n = b + (z+3)/(6*polynomial)
    nr = -(s*s+2*dpsi/R+4/R**2)/(k*R)
    return {
        **coefficients, "laplace_s": s, "radius": R, "decay_q": q,
        "growth_bound": growth, "orientation": orientation,
        "psi_R_over_psi": dpsi, "chi_R_per_chi": dc,
        "outward_chi_normal_per_chi": orientation*dc,
        "radial_metric_b_per_chi": b, "lapse_n_per_chi": n,
        "lapse_n_R_per_chi": nr,
        "unreferenced_energy_per_chi": -4*orientation*alpha/R,
        "unreferenced_area_pressure_per_chi": -orientation*alpha*(s*s+1/R**2)/(2*pi*R),
        "auxiliary_proper_response_per_chi": (
            4*orientation*alpha*(s*s+3*dpsi/R+6/R**2)/(k*R)
        ),
    }
