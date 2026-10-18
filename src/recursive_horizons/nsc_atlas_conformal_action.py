"""Original Jordan action in unit-sphere conformal variables, per 4*pi.

g_J=r^2*(gbar2+dOmega2), chi=r*phi_J. Boundary terms are retained.
These local action/junction/Legendre interfaces are not a new evolved
solution, a correction force or a replacement for the P001 bulk tests.
"""
from __future__ import annotations

import math

import numpy as np


def _configuration(radius, field, kappa):
    r, chi, k = float(radius), float(field), float(kappa)
    if not all(math.isfinite(v) for v in (r, chi, k)) or r <= 0 or k <= 0:
        raise ValueError('require finite r>0, kappa>0 and invariant field')
    F = 1-k*chi*chi/(6*r*r)
    if F <= 0:
        raise ValueError('the original positive-coupling branch F>0 is required')
    return r, chi, k, F


def reduced_bulk_density(*, radius, field, kappa, vacuum_density,
                         radial_curvature, grad_radius_squared, grad_field_squared):
    """sqrt(-gbar2) is external; grad squares use that Lorentzian two-metric."""
    r, chi, k, _ = _configuration(radius, field, kappa)
    V0, R2, dr2, dc2 = (float(v) for v in (vacuum_density, radial_curvature,
                                         grad_radius_squared, grad_field_squared))
    if not all(math.isfinite(v) for v in (V0, R2, dr2, dc2)):
        raise ValueError('finite one-sided vacuum, curvature and gradient scalars are required')
    X = r*r/(2*k)-chi*chi/12
    return X*(R2+2)+3*dr2/k-dc2/2-V0*r**4


def junction_matrix(*, radius, field, kappa):
    """Rows act on ([nbar r],[nbar chi],[Kbar2]), with one common normal."""
    r, chi, k, _ = _configuration(radius, field, kappa)
    return np.array([[r, -k*chi/6, 0.], [3., 0., r], [0., 1., chi/3]])


def tension_junction(*, radius, field, kappa, wall_k):
    r, chi, k, F = _configuration(radius, field, kappa)
    b = float(wall_k)
    if not math.isfinite(b) or b < 0:
        raise ValueError('require finite nonnegative original tension k_wall=kappa*sigma/2')
    matrix = junction_matrix(radius=r, field=chi, kappa=k)
    rhs = np.array([-b*r**3, -3*b*r*r, 0.])
    jumps = np.linalg.solve(matrix, rhs)
    return dict(normal_radius_jump=float(jumps[0]), normal_field_jump=float(jumps[1]),
        extrinsic_curvature_jump=float(jumps[2]), F=F,
        determinant=-r*r*F, residual_max=float(np.max(np.abs(matrix @ jumps-rhs))))


def time_kinetic_matrix(*, radius, field, kappa, radial_length, lapse):
    """p=M*v per 4*pi, q=(a,r,chi), v=(a_t-(a*b)_y,r_t-b*r_y,chi_t-b*chi_y).

    a,N belong to gbar2. This gravitational matrix need not be positive:
    lapse/shift constraints still have to be imposed. It is nonsingular
    on the admitted F>0 branch at r>0, with det=-4*a*r^2*F/(kappa^2*N^3).
    """
    r, chi, k, _ = _configuration(radius, field, kappa)
    a, N = float(radial_length), float(lapse)
    if not all(math.isfinite(v) and v > 0 for v in (a, N)):
        raise ValueError('positive finite conformal radial length and lapse are required')
    return np.array([[0., -2*r/k, chi/3], [-2*r/k, -6*a/k, 0.], [chi/3, 0., a]])/N


def velocities_from_momenta(momentum, **configuration):
    """Recover actual advected field/metric velocities; no momentum-as-velocity alias."""
    p = np.asarray(momentum, dtype=float)
    if p.shape != (3,) or not np.all(np.isfinite(p)):
        raise ValueError('require three finite local canonical momenta per 4*pi')
    return np.linalg.solve(time_kinetic_matrix(**configuration), p)
