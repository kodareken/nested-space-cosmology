"""Flat Euclidean TT response of a fixed matched Dirac spectral profile.

For Gamma=(1/2) Tr H(D4^2), remove its already selected volume term.
The quadratic density is (A/4) h_TT P_H(-partial^2) h_TT. Compact scales,
warp and matching parameters are fixed during this external metric variation.
This is not a curved retarded response or a Lorentzian state prescription.
"""
from __future__ import annotations

from math import isfinite, pi, sqrt
from typing import Callable

import numpy as np
from scipy.special import roots_jacobi


def profile_local_coefficients(value_at_zero: float, first_moment: float) -> dict[str, float]:
    """The existing compact-matching normalization, Q1=int_0^infty H(u) du."""
    q0, q1 = float(value_at_zero), float(first_moment)
    if not (isfinite(q0) and isfinite(q1) and q1 > 0):
        raise ValueError("finite H(0) and positive finite first spectral moment required")
    norm = (4*pi)**2
    return {"A": q1/(6*norm), "C_W": -q0/(40*norm),
            "local_quartic_per_quadratic": -3*q0/(10*q1)}


def flat_spin2_profile_response(
    squared_momentum: float,
    profile: Callable[[float], float],
    first_moment: float,
    *,
    quadrature_order: int = 32,
) -> dict[str, float | int | None]:
    """Evaluate P_H(x) without inverse Laplace reconstruction or moment refits.

    P_H=x [1-x/(4 Q1) int_0^1 (1+2v) sqrt(1-v) H(xv/4) dv].
    This follows from the linear heat response whenever its spectral
    functional calculus and integrations exist. H must have the supplied
    first moment and a finite continuous limit at zero. Positivity or a
    positive Laplace measure is not assumed. Sample extrema are not bounds.
    """
    x, q1 = float(squared_momentum), float(first_moment)
    if not (isfinite(x) and x >= 0 and isfinite(q1) and q1 > 0):
        raise ValueError("finite nonnegative squared momentum and positive first moment required")
    if (isinstance(quadrature_order, bool) or not isinstance(quadrature_order, int)
            or quadrature_order < 2):
        raise ValueError("Jacobi quadrature requires an integer order of at least two")
    if x == 0:
        return {"squared_momentum": x, "response": 0.0, "response_over_Einstein": 1.0,
                "correction_integral": 0.0, "quadrature_order": quadrature_order,
                "sampled_profile_min": None, "sampled_profile_max": None}
    nodes, weights = roots_jacobi(quadrature_order, 0.5, 0.0)
    v = (nodes+1)/2
    values = np.asarray([float(profile(float(x*y/4))) for y in v])
    if not np.isfinite(values).all():
        raise ValueError("profile returned a nonfinite value inside the integration domain")
    # The Jacobi weight (1-node)^(1/2) includes the square-root endpoint.
    integral = float(x/4*np.dot(weights/(2*sqrt(2)), (1+2*v)*values))
    ratio = 1-integral/q1
    return {"squared_momentum": x, "response": x*ratio,
            "response_over_Einstein": ratio, "correction_integral": integral,
            "quadrature_order": quadrature_order,
            "sampled_profile_min": float(np.min(values)),
            "sampled_profile_max": float(np.max(values))}
