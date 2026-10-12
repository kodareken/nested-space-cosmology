"""Derivative-coherent local (w,U) histories for the approved incoming gate.

The finite basis parameterizes geometry. It does not fit source data, choose
a cosmological duration, or turn a collocation residual into an existence
certificate. Every candidate still requires the owned evolved source.
"""
from dataclasses import dataclass
from math import comb, fsum

import numpy as np
from numpy.polynomial import Chebyshev
from scipy.special import expit

from .nsc_compatible_history_geometry import CompatibleIncomingMetric, CompatibleRadiusDirection
from .nsc_ks_spacetime_variation import chart_coordinates


def plateau_derivatives(z, center, inner=.03, outer=.06):
    """One C-infinity plateau and its actual derivatives through order four."""
    z, center, inner, outer = map(float, (z, center, inner, outer))
    if not np.isfinite([z, center, inner, outer]).all() or not 0 < inner < outer:
        raise ValueError('finite axial coordinates and 0 < inner < outer required')
    distance = abs(z-center)
    if distance <= inner:
        return np.array([1., 0., 0., 0., 0.])
    if distance >= outer:
        return np.zeros(5)
    width = outer-inner
    u = (distance-inner)/width
    # q = 1/(1-u)-1/u. Evaluate both logistic factors independently:
    # 1-eta rounds to zero prematurely on the inner edge.
    q = 1/(1-u)-1/u
    eta, complement = expit(-q), expit(q)
    product = eta*complement
    q1 = (1-u)**-2+u**-2
    q2 = 2*((1-u)**-3-u**-3)
    q3 = 6*((1-u)**-4+u**-4)
    q4 = 24*((1-u)**-5-u**-5)
    f1 = -product
    f2 = product*(1-2*eta)
    f3 = -product*(1-6*eta+6*eta**2)
    f4 = product*(1-14*eta+36*eta**2-24*eta**3)
    values = np.array([eta, f1*q1, f2*q1*q1+f1*q2,
                       f3*q1**3+3*f2*q1*q2+f1*q3,
                       f4*q1**4+6*f3*q1*q1*q2+3*f2*q2*q2+4*f2*q1*q3+f1*q4])
    values *= (np.sign(z-center)/width)**np.arange(5)
    if not np.isfinite(values).all():
        raise ArithmeticError('axial plateau derivative overflow')
    return values


@dataclass(frozen=True)
class LocalAxialFunction:
    coefficients: tuple
    center: float
    inner: float = .03
    outer: float = .06

    def __post_init__(self):
        raw = np.asarray(self.coefficients)
        if (raw.ndim != 1 or len(raw) == 0 or np.iscomplexobj(raw)
                or not np.isfinite(raw).all()):
            raise ValueError('finite real Chebyshev coefficients required')
        plateau_derivatives(self.center, self.center, self.inner, self.outer)
        coeff = tuple(map(float, raw))
        poly = Chebyshev(coeff, domain=[self.center-self.outer, self.center+self.outer])
        object.__setattr__(self, 'coefficients', coeff)
        object.__setattr__(self, '_derivatives', tuple(poly.deriv(j) for j in range(5)))

    def __call__(self, z, derivative_order):
        if isinstance(derivative_order, bool) or derivative_order not in range(5):
            raise ValueError('owned coherent derivative orders zero through four required')
        z = float(z)
        window = plateau_derivatives(z, self.center, self.inner, self.outer)
        if not np.any(window):
            return 0.
        n = derivative_order
        return float(fsum(comb(n, k)*window[k]*self._derivatives[n-k](z)
                          for k in range(n+1)))


@dataclass(frozen=True)
class LocalIncomingFamily:
    """Two coefficient vectors, one actual history, and a declared local I.

    Amplitude order is all w coefficients, then all U coefficients. The
    existing geometry owner computes every incoming mixed derivative from
    these same functions. Source preparation is not part of this object.
    """
    coefficients: object
    center: float | None = None
    axial_inner: float = .03
    axial_outer: float = .06
    normal_inner: float = .007
    normal_outer: float = .03

    def __post_init__(self):
        raw = np.asarray(self.coefficients)
        if (raw.ndim != 2 or raw.shape[0] != 2 or raw.shape[1] not in (8, 16, 32, 64)
                or np.iscomplexobj(raw) or not np.isfinite(raw).all()):
            raise ValueError('two finite real vectors with 8, 16, 32 or 64 coefficients required')
        center = chart_coordinates(1.)[1]+.15 if self.center is None else float(self.center)
        plateau_derivatives(center, center, self.axial_inner, self.axial_outer)
        if not np.isfinite([self.normal_inner,self.normal_outer]).all() or not 0 < self.normal_inner < self.normal_outer:
            raise ValueError('positive ordered normal support radii required')
        object.__setattr__(self, 'center', center)
        object.__setattr__(self, 'coefficients', tuple(tuple(map(float,row)) for row in raw))

    @property
    def interval(self):
        return (self.center-self.axial_inner, self.center+self.axial_inner)

    @property
    def support(self):
        return (self.center-self.axial_outer, self.center+self.axial_outer)

    @property
    def functions(self):
        return tuple(LocalAxialFunction(row,self.center,self.axial_inner,self.axial_outer)
                     for row in self.coefficients)

    def metric(self):
        n = len(self.coefficients[0])
        basis = [LocalAxialFunction(tuple(row),self.center,self.axial_inner,self.axial_outer)
                 for row in np.eye(n)]
        zero = LocalAxialFunction((0.,),self.center,self.axial_inner,self.axial_outer)
        directions = tuple(CompatibleRadiusDirection(w,zero,self.normal_inner,self.normal_outer) for w in basis)
        directions += tuple(CompatibleRadiusDirection(zero,U,self.normal_inner,self.normal_outer) for U in basis)
        return CompatibleIncomingMetric(np.asarray(self.coefficients).ravel(), directions)

    def radius_lower_bound(self):
        """Conservative global radius bound from |T_j|<=1 on the cutoff.

        The reference radius is >=1. Directed upward cushions cover each
        nonnegative binary64 sum/product; this bounds this finite geometry,
        not the state or the residual between collocation nodes.
        """
        up = lambda v: float(np.nextafter(float(v), np.inf))
        w, U = (up(fsum(abs(v) for v in row)) for row in self.coefficients)
        s = float(self.normal_outer)
        s3 = up(up(s*s)*s)
        perturbation = up(up(s*w)+up(up(s3/6)*U))
        return float(np.nextafter(1.-perturbation, -np.inf))

    def collocation_nodes(self, count=None):
        count = 2*len(self.coefficients[0])+1 if count is None else count
        if isinstance(count,bool) or not isinstance(count,(int,np.integer)) or count < 3:
            raise ValueError('at least three local collocation nodes required')
        return self.center+self.axial_inner*np.cos(np.linspace(np.pi,0,int(count)))

    def description(self):
        return {'class': 'local source-fixed pure-radius compatible histories',
                'interval': list(self.interval), 'axial_support': list(self.support),
                'normal_window': [self.normal_inner,self.normal_outer],
                'coefficients': [list(row) for row in self.coefficients],
                'coefficient_order': 'w then U; increasing Chebyshev degree',
                'radius_lower_bound': self.radius_lower_bound(),
                'duration_selected': False, 'state_supplied_by_family': False,
                'physical_gate': 'OPEN; actual evolved full source and error enclosure required'}
