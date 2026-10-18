"""Orthogonal timelike Dirichlet faces of the owned auxiliary bulk action.

The representative is sqrt|h| [F R_h + Z (grad r)^2 + V], h=(+L^2,-Q^2),
beta=0. Variations vanish near the time caps; no moving-face/corner law or
matter boundary state is supplied. The minimal completion adds no intrinsic
face potential. Its charges are unreferenced gravitational responses, not a
parent mass, a bulk density, or a cosmological constant.
"""
from __future__ import annotations

from math import isfinite, pi as math_pi

import sympy as sp

from .nsc_spherical_feedback_action import (
    feedback_F, feedback_Z, partial_F, require_weyl_chart,
)


def _positive(name, value):
    if isinstance(value, sp.Expr):
        if value.is_positive is not True:
            raise ValueError(name + " requires an explicitly positive symbol/value")
    else:
        try:
            valid = isfinite(float(value)) and value > 0
        except (TypeError, ValueError, OverflowError):
            valid = False
        if not valid:
            raise ValueError(name + " must be finite and positive")


def _finite_real(name, value):
    if isinstance(value, sp.Expr):
        valid = value.is_real is not False and value.is_finite is not False and value is not sp.nan
    else:
        try:
            valid = isfinite(float(value))
        except (TypeError, ValueError, OverflowError):
            valid = False
    if not valid:
        raise ValueError(name + " must be finite and real")


def _domain(L, Q, r, A, C_W, orientation):
    for name, value in (("L", L), ("Q", Q), ("r", r)):
        _positive(name, value)
    require_weyl_chart(C_W, A)
    for name, value in (("A", A), ("C_W", C_W)):
        _finite_real(name, value)
    if orientation not in (-1, 1):
        raise ValueError("outward orientation must be +1 or -1")


def timelike_completion(L, Q, r, chi, L_x, A, C_W, *, orientation=1):
    """Face density -2 s F L_x/Q, integrated against coordinate time.

    This cancels the spatial curvature divergence +2(F L_x/Q)_x.
    It is the epsilon=-1 term 2 epsilon sqrt|gamma_h| F K_h, with
    K_h=s L_x/(L Q). Other bulk total divergences need their own completion.
    """
    _domain(L, Q, r, A, C_W, orientation)
    _finite_real("chi", chi)
    _finite_real("L_x", L_x)
    return -2 * orientation * feedback_F(r, chi, A, C_W) * L_x / Q


def timelike_gravitational_response(
    L, Q, r, chi, L_x, r_x, chi_x, A, C_W, *, orientation=1,
):
    """Coefficients of the completed face variation at fixed x and beta=0.

    delta S_face = dt [Pi_N delta N + Pi_r delta r + Pi_chi delta chi],
    N=rL. Define E=-Pi_N and P=Pi_r/(8 pi r N), so the same form is
    dt [-E delta N + N P delta(4 pi r^2) + Pi_chi delta chi].
    This definition fixes no reference subtraction or exterior physical state.
    L is varied independently of Q before any conformal gauge substitution.
    """
    _domain(L, Q, r, A, C_W, orientation)
    for name, value in (("chi", chi), ("L_x", L_x), ("r_x", r_x), ("chi_x", chi_x)):
        _finite_real(name, value)
    F_r, F_chi = partial_F(r, A, C_W)
    F_x = F_r * r_x + F_chi * chi_x
    Z = feedback_Z(A)
    sign = orientation
    coordinate_L = -2 * sign * F_x / Q
    coordinate_r = -2 * sign * (F_r * L_x + Z * L * r_x) / Q
    coordinate_chi = -2 * sign * F_chi * L_x / Q
    physical_N = coordinate_L / r
    physical_r = coordinate_r - L * coordinate_L / r
    N = r * L
    pi = sp.pi if any(isinstance(v, sp.Expr) for v in (L, Q, r, A, C_W)) else math_pi
    return {
        "orientation": sign,
        "N": N,
        "coordinate": {"L": coordinate_L, "r": coordinate_r, "chi": coordinate_chi, "Q": 0},
        "physical": {"N": physical_N, "r": physical_r, "chi": coordinate_chi},
        "unreferenced_energy": -physical_N,
        "unreferenced_area_pressure": physical_r / (8 * pi * r * N),
        "auxiliary_response_per_proper_time": coordinate_chi / N,
    }
