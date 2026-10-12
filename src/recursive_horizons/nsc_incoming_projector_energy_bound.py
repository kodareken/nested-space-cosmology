"""Incoming energy-specific error of the same pure vacuum projector.

The rank-one relation retains the quadratic longitudinal error instead of
bounding every projector direction by the large instantaneous energy.
Thermal and numerical quadrature errors remain separate.
"""
import math
import mpmath as mp
import sympy as sp

from .nsc_incoming_vacuum_tail_bound import (
    _precision, _lo, _hi, _range, _up_float, _geometry_jets, _mul, _derivative,
    IncomingVacuumTailBoundGeometry,
)
from .nsc_incoming_middle_bound import _parameters


def rank_one_energy_identity():
    """Exact geometry of two rank-one 2x2 projectors, no state interpolation."""
    d, longitudinal, transverse = sp.symbols('d longitudinal transverse', real=True)
    # |Delta n|=2d and n.Delta n=-2d^2 for |n|=|n+Delta n|=1.
    squared_transverse = sp.expand(4*d*d-(-2*d*d)**2)
    residual = sp.expand(squared_transverse-4*d*d*(1-d*d))
    if residual != 0: raise ArithmeticError('rank-one projector identity failed')
    return {'transverse_squared_residual': str(residual),
            'delta_n_longitudinal': '-2*d^2',
            'energy_error': '2*norm(cross(h,n))*d + 2*abs(dot(h,n))*d^2',
            'domain': 'two normalized rank-one vacuum projectors; d=operator_norm(P_exact-P_approx)'}


def incoming_riccati_interval_coefficients(channel, *, order=16):
    """The owned recurrence at rho=1, with directed coefficient arithmetic."""
    A, sphere = _geometry_jets(mp.iv.mpf(1), order)
    _, _, mass, angular, _ = _parameters(channel)
    upper = [v*angular for v in sphere]; upper[0] += mp.iv.j*mass
    lower = [mp.iv.mpc(v.real, -v.imag) for v in upper]
    terms = [upper]; derivative_A = _derivative(A)
    for n in range(1, order):
        first = _mul(A, _derivative(terms[-1])); second = _mul(derivative_A, terms[-1])
        product = [mp.iv.mpc(0) for _ in A]
        for j in range(1, n):
            value = _mul(terms[j-1], terms[n-j-1])
            product = [a+b for a, b in zip(product, value)]
        nonlinear = _mul(_mul(A, lower), product)
        terms.append([mp.iv.j*(a+b/2)+c for a, b, c in zip(first, second, nonlinear)])
    return [row[0] for row in terms]


def local_cross_product_coefficients(channel, *, order=16, precision=40):
    """Bound |h cross n_order| by sum_j B_j/(2E)^j, j=1,...,2*order.

    Exact c1=L+i*m cancellation is applied before arithmetic. The positive
    normalization denominator 1+a^2|S|^2 is bounded below by1.
    """
    if order not in (16, 24):
        raise ValueError('explicit existing16 or improved24 numerical order required')
    with _precision(precision):
        a, r, mass, ell, _ = _parameters(channel); L = ell/r
        c = [mp.iv.mpc(0)]+incoming_riccati_interval_coefficients(channel, order=order)
        u = [mp.iv.mpf(0) for _ in range(2*order+1)]
        for j in range(2, 2*order+1):
            u[j] = a*a*sum((c[i].real*c[j-i].real+c[i].imag*c[j-i].imag
                for i in range(1, order+1) if 1 <= j-i <= order), mp.iv.mpf(0))
        coefficients = []
        for j in range(1, 2*order+1):
            x, y, z = -L*u[j], -mass*u[j], mp.iv.mpf(0)
            if j+1 <= order:
                x -= c[j+1].real; y -= c[j+1].imag
            if 2 <= j <= order:
                z = 2*a*(mass*c[j].real-L*c[j].imag)
            coefficients.append(_up_float(mp.iv.sqrt(x*x+y*y+z*z)))
        return {'first_inverse_power': 1, 'last_inverse_power': 2*order,
                'cross_numerator_coefficient_upper': coefficients,
                'physical_Riccati_order': order, 'normalization_denominator_lower': 1.,
                'constant_coefficient_exact': 0., 'c1': 'L+i*m',
                'angular_sign_norm_symmetry': 'c_n(-lambda)=(-1)^n conjugate(c_n(lambda))'}


def _primitive(power, moment, left, right):
    p = power-moment-1
    if p <= 0 or not 0 < left < right:
        raise ValueError('integrable positive finite-energy powers required')
    return (mp.iv.mpf(left)**(-p)-mp.iv.mpf(right)**(-p))/p


def projected_energy_bound(channel, radial, cross, left, right, *, precision=40):
    """Reuse radial d<=sum I_n/(2E)^n in the rank-one energy bound."""
    if radial['group'] != channel['index'] or radial['lower'] != right:
        raise ValueError('owned group and finite-band endpoint required')
    order = cross['physical_Riccati_order']
    if order not in (16, 24) or radial.get('physical_Riccati_order', 16) != order or cross['constant_coefficient_exact'] != 0.:
        raise ValueError('matching explicit numerical order and exact leading cancellation required')
    B = cross['cross_numerator_coefficient_upper']
    if len(B) != 2*order or any(not math.isfinite(v) or v < 0 for v in B):
        raise ValueError('all finite nonnegative local cross coefficients required')
    signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
    if tuple(row['sign'] for row in radial['per_sign']) != signs:
        raise ValueError('actual angular signs required exactly once')
    with _precision(precision):
        a, r, m, ell, factor = _parameters(channel); factor /= len(signs)
        M = mp.iv.sqrt(m*m+(ell/r)**2)
        linear, quadratic = mp.iv.mpf(0), mp.iv.mpf(0)
        for row in radial['per_sign']:
            values = row['I16_to_I32_upper'] if order == 16 else row['coefficient_integral_upper']
            if len(values) != order+1 or any(not math.isfinite(v) or v < 0 for v in values):
                raise ValueError('owned radial defect integral bounds required')
            for p, integral in enumerate(values, order):
                for j, local in enumerate(B, 1):
                    linear += 2*factor*mp.iv.mpf(local)*mp.iv.mpf(integral)/mp.iv.mpf(2)**(p+j)*_primitive(p+j, 0, left, right)
                for q, other in enumerate(values, order):
                    quadratic += 2*factor*mp.iv.mpf(integral)*mp.iv.mpf(other)/mp.iv.mpf(2)**(p+q)*(
                        _primitive(p+q, 1, left, right)/a+M*_primitive(p+q, 0, left, right))
        total = linear+quadratic
        return {'density_error_upper': _up_float(total),
                'linear_transverse_density_upper': _up_float(linear),
                'quadratic_longitudinal_density_upper': _up_float(quadratic),
                'lapse_action_error_upper': _up_float(4*mp.iv.pi*a*r*r*total),
                'vacuum_current_error': 0., 'physical_Riccati_order': order,
                'thermal_quadrature_low_band_errors_included': False}


class RefinedDefectGeometry(IncomingVacuumTailBoundGeometry):
    """Explicit bounded refinement of the same directed collar partition."""
    def __init__(self, config, *, intervals=64, precision=40):
        if intervals not in (32, 64, 128):
            raise ValueError('new bounded partition32,64,128 required')
        super().__init__(config, intervals=16, precision=precision)
        self.intervals = intervals
        with _precision(precision):
            self.cells = []
            for i in range(intervals):
                t = _range(mp.mpf(i)/intervals, mp.mpf(i+1)/intervals)
                rho = 1+(self.horizon-1)*(1-t*t)
                self.cells.append((rho, *_geometry_jets(rho, 16)))
