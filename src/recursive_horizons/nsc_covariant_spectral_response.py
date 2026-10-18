"""Curvature-squared profile kernels and conserved static weak-field response.

Gamma=(1/2) Tr_4 H(D4^2), with fixed compact data and the existing volume
subtraction. The curvature action is [C C_H(-nabla^2) C + R R_H(-nabla^2) R]
/[2(4*pi)^2], modulo Euler/boundary terms in a perturbative fixed topology
and O(curvature^3). This supplies no real-time state or strong-curvature wall.
"""
from __future__ import annotations

from math import isfinite, sqrt
from typing import Callable

import numpy as np
from scipy.special import roots_jacobi


def curvature_profile_kernels(
    squared_momentum: float,
    profile: Callable[[float], float],
    value_at_zero: float,
    *,
    quadrature_order: int = 32,
) -> dict[str, float | int]:
    """The covariant Weyl/scalar functions; a single set of profile samples.

    The scalar integral's weight has exactly zero mean. Subtracting H(0)
    makes its small-momentum cancellation explicit without changing it.
    The same regular spectral-functional calculus as the TT owner is assumed.
    """
    x, h0 = float(squared_momentum), float(value_at_zero)
    if not (isfinite(x) and x >= 0 and isfinite(h0)):
        raise ValueError("finite nonnegative squared momentum and finite H(0) required")
    if (isinstance(quadrature_order,bool) or not isinstance(quadrature_order,int)
            or quadrature_order < 2):
        raise ValueError("Jacobi quadrature requires an integer order of at least two")
    if x == 0:
        return {'squared_momentum':x,'weyl_kernel':-h0/20,'scalar_kernel':0.0,
                'quadrature_order':quadrature_order}
    nodes, weights = roots_jacobi(quadrature_order,0.5,0.0)
    v, weights = (nodes+1)/2, weights/(2*sqrt(2))
    values = np.asarray([float(profile(float(x*y/4))) for y in v])
    if not np.isfinite(values).all():
        raise ValueError("profile returned a nonfinite curvature-kernel sample")
    return {'squared_momentum':x,
            'weyl_kernel':-float(np.dot(weights,(1+2*v)*values))/24,
            'scalar_kernel':float(np.dot(weights,(1-2.5*v)*(values-h0)))/72,
            'quadrature_order':quadrature_order}


def static_gravitational_response(
    squared_momentum: float,
    weyl_kernel: float,
    scalar_kernel: float,
    first_moment: float,
) -> dict[str,float | None]:
    """Potential ratios for a conserved, static T00=rho probe at linear order.

    ds^2=(1+2 Phi)dt^2-(1-2 Psi)dx^2; Phi_Einstein=-rho/(4 A p^2).
    At x=0 only the continuous transfer-ratio limit is returned, not the
    inverse Poisson response to a homogeneous density. No causal kernel is
    selected by this static calculation.
    """
    x, fc, fr, q1 = map(float,(squared_momentum,weyl_kernel,scalar_kernel,first_moment))
    if not (all(isfinite(v) for v in (x,fc,fr,q1)) and x >= 0 and q1 > 0):
        raise ValueError("finite curvature kernels, nonnegative momentum and positive moment required")
    a2, a0 = 1+6*x*fc/q1, 1-18*x*fr/q1
    if a2 == 0 or a0 == 0:
        raise ValueError("the static transfer is singular at this momentum")
    phi, psi = 4/(3*a2)-1/(3*a0), 2/(3*a2)+1/(3*a0)
    return {'spin2_factor':a2,'scalar_trace_factor':a0,
            'Newton_potential_ratio':phi,'spatial_potential_ratio':psi,
            'lensing_potential_ratio':(phi+psi)/2,
            'slip_Psi_over_Phi':None if phi == 0 else psi/phi}
