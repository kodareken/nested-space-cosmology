"""Common supernova likelihood for inherited parent and FLRW controls.

Only the existing protected-future source is compared. An observer epoch
can be calibrated; neither the scalar normalization nor the curvature is fit.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

from .nsc_atlas_inherited_matter import exact_initial_geometry, uniform_bound


@lru_cache(maxsize=4)
def _gauss(order):
    if order not in (32, 64, 96):
        raise ValueError('quadrature order must be 32, 64 or 96')
    return np.polynomial.legendre.leggauss(order)


def flrw_distance(z_hd, z_hel, *, matter, curvature=0., order=64):
    """D_L/(c/H_o) with release heliocentric prefactor; radiation omitted.

    Omega_k>0 is open. The closed branch is restricted before its antipode.
    """
    z, hel = np.broadcast_arrays(np.asarray(z_hd, float), np.asarray(z_hel, float))
    if (z.ndim != 1 or not np.isfinite(z).all() or not np.isfinite(hel).all()
            or np.any(z <= 0) or np.any(hel <= -1)
            or not math.isfinite(matter) or not math.isfinite(curvature)
            or matter < 0 or 1-matter-curvature < 0):
        raise ValueError('finite positive redshifts and nonnegative matter/vacuum required')
    nodes, weights = _gauss(order)
    zz = z[:, None]*(nodes+1)/2
    E2 = matter*(1+zz)**3+curvature*(1+zz)**2+1-matter-curvature
    if np.any(E2 <= 0):
        raise ValueError('light cone leaves expanding branch')
    integral = z/2*((1/np.sqrt(E2)) @ weights)
    if curvature < 0:
        angle = math.sqrt(-curvature)*integral
        if np.any(angle >= math.pi/2):
            raise ValueError('closed ray leaves the monotone sine branch')
        radial = integral*np.sinc(angle/np.pi)
    elif curvature > 0:
        angle = math.sqrt(curvature)*integral
        radial = np.sinh(angle)/math.sqrt(curvature)
    else:
        radial = integral
    return (1+hel)*radial


def parent_fractions(inverse_observer_scale):
    """s=1/A_o; s=0 is the infinite-future vacuum control, not a finite age."""
    s = float(inverse_observer_scale)
    if not math.isfinite(s) or not 0 <= s <= 1:
        raise ValueError('protected-future inverse scale must be in [0,1]')
    d = exact_initial_geometry()
    vacuum, matter, C = float(d['vacuum']/3), float(d['excitation']/3), float(d['curvature'])
    H2 = vacuum+matter*s**3-C*s*s
    return dict(matter=matter*s**3/H2, vacuum=vacuum/H2,
                curvature=-C*s*s/H2, H_squared=H2)


def parent_distance(z_hd, z_hel, inverse_observer_scale, *, order=64, minimum_scale=1.):
    z = np.asarray(z_hd, float)
    s = float(inverse_observer_scale)
    if minimum_scale not in (1.,1/3):
        raise ValueError('only the proved future or A=1/3 continuation is available')
    if z.size == 0 or s*(1+float(np.max(z))) > 1/minimum_scale+2e-14:
        raise ValueError('emitter precedes the proved initial-scale boundary')
    fractions = parent_fractions(s)
    return flrw_distance(z_hd, z_hel, matter=fractions['matter'],
                         curvature=fractions['curvature'], order=order)


def magnitude_error_bound(omega=10**9, *, minimum_scale=1.):
    """Bound on unnormalised physical D_L magnitudes, up to a common offset.

    The corrected-energy theorem holds on A>=1. Sine is increasing and
    concave on the verified short-ray branch. The normalising actual H_o
    differs only by a common magnitude offset, removed by the likelihood.
    """
    if minimum_scale==1.:
        bound=uniform_bound(omega)
    elif minimum_scale==1/3:
        from .nsc_atlas_parent_past import past_bound
        bound=past_bound(omega)
    else:
        raise ValueError('unproved minimum scale')
    r = float(bound['hubble_squared_relative_error'])
    return 5/(2*math.log(10))*(-math.log1p(-r))


def parent_ray_angle_upper(z_max, inverse_observer_scale, omega=10**9, *, minimum_scale=1.):
    """Conservative absolute angle using H>=H_- and A_e>=1."""
    s = float(inverse_observer_scale)
    if z_max < 0 or s < 0 or s*(1+z_max) > 1/minimum_scale+2e-14:
        raise ValueError('future endpoints required')
    # uniform_bound verifies H_->0.266 for the same initial geometry.
    uniform_bound(omega)
    C = float(exact_initial_geometry()['curvature'])
    return math.sqrt(C)*s*z_max/.266


class SupernovaShapeLikelihood:
    """Fixed covariance, analytic common-offset projection, no diagonal jitter."""

    def __init__(self, magnitude, covariance):
        self.magnitude = np.asarray(magnitude, float)
        C = np.asarray(covariance, float)
        n = self.magnitude.size
        if (self.magnitude.ndim != 1 or n < 3 or C.shape != (n,n)
                or not np.isfinite(C).all() or not np.isfinite(self.magnitude).all()
                or np.max(abs(C-C.T)) > 1e-12*max(1., float(np.max(abs(C))))):
            raise ValueError('finite aligned magnitudes and symmetric covariance required')
        self.factor = cho_factor((C+C.T)/2, lower=True)
        self.ones = np.ones(n)
        self.precision_ones = cho_solve(self.factor, self.ones)
        self.normalization = float(self.ones @ self.precision_ones)

    def evaluate(self, distance):
        distance = np.asarray(distance, float)
        if (distance.shape != self.magnitude.shape or np.any(distance <= 0)
                or not np.isfinite(distance).all()):
            raise ValueError('aligned finite positive distances required')
        residual = self.magnitude-5*np.log10(distance)
        offset = float(self.precision_ones @ residual/self.normalization)
        centered = residual-offset
        return dict(chi_squared=float(centered @ cho_solve(self.factor, centered)),
                    magnitude_intercept=offset, offset_score=float(self.precision_ones @ centered))


def fit_shape(likelihood, distance, upper):
    """Bounded one-parameter calibration plus explicit endpoint controls."""
    evaluations = 0
    def objective(value):
        nonlocal evaluations
        evaluations += 1
        if evaluations > 200:
            raise RuntimeError('shape fit exceeded 200 objective evaluations')
        return likelihood.evaluate(distance(value))['chi_squared']
    fit = minimize_scalar(objective, bounds=(0., upper), method='bounded',
                          options=dict(xatol=2e-11, maxiter=180))
    if not fit.success:
        raise RuntimeError('bounded shape calibration failed')
    choices = [(float(fit.fun), float(fit.x)), (objective(0.), 0.), (objective(upper), upper)]
    _, selected = min(choices)
    return dict(parameter=selected, at_boundary=selected in (0., upper),
                objective_evaluations=evaluations, **likelihood.evaluate(distance(selected)))
