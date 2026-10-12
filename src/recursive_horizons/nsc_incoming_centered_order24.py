"""Trimmed Taylor enclosure of the unchanged group14 order24 projector."""
import mpmath as mp

from .nsc_incoming_vacuum_tail_bound import (
    _precision, _lo, _hi, _range, _geometry_jets, _up_float,
    IncomingVacuumTailBoundGeometry,
)
from .nsc_incoming_defect_taylor_bound import centered_coefficient_enclosure
from .nsc_incoming_middle_bound import _parameters


def _product(a, b, degree):
    """Required low convolution only; no replacement of missing derivatives."""
    if min(len(a), len(b)) <= degree:
        raise ValueError('complete low product jets required')
    return [sum((a[j]*b[n-j] for j in range(n+1)), mp.iv.mpc(0))
            for n in range(degree+1)]


def _differentiate(a, degree):
    if len(a) <= degree+1:
        raise ValueError('next derivative coefficient required')
    return [(j+1)*a[j+1] for j in range(degree+1)]


def trimmed_defect_jets(A, sphere, mass, angular, *, order=24, depth=4):
    """f_order..f_2order normalized coordinate jets through degree depth.

    c_j requires degree order+depth-j+1. The last c_order derivative
    therefore has degree depth, without any invented terminal derivative.
    Order16 is available solely for comparison to the existing owner.
    """
    if order not in (16, 24) or depth != 4:
        raise ValueError('order16 ownership control or order24, depth4 required')
    if min(len(A), len(sphere)) <= order+depth:
        raise ValueError('geometry through degree order+depth required')
    degree = order+depth
    upper = [angular*v for v in sphere[:degree+1]]
    upper[0] += mp.iv.j*mass
    lower = [mp.iv.mpc(v.real, -v.imag) for v in upper]
    terms = [upper]
    for n in range(1, order):
        degree = order+depth-n
        first = _product(A, _differentiate(terms[-1], degree), degree)
        second = _product(_differentiate(A, degree), terms[-1], degree)
        product = [mp.iv.mpc(0) for _ in range(degree+1)]
        for j in range(1, n):
            value = _product(terms[j-1], terms[n-j-1], degree)
            product = [a+b for a, b in zip(product, value)]
        nonlinear = _product(_product(A, lower, degree), product, degree)
        terms.append([mp.iv.j*(a+b/2)+c for a, b, c in zip(first, second, nonlinear)])
    answer = []
    for power in range(order, 2*order+1):
        square = [mp.iv.mpc(0) for _ in range(depth+1)]
        for j in range(1, order+1):
            if 1 <= power-j <= order:
                value = _product(terms[j-1], terms[power-j-1], depth)
                square = [a+b for a, b in zip(square, value)]
        coefficient = [-v for v in _product(_product(A, lower, depth), square, depth)]
        if power == order:
            first = _product(A, _differentiate(terms[-1], depth), depth)
            second = _product(_differentiate(A, depth), terms[-1], depth)
            coefficient = [c-mp.iv.j*(a+b/2) for c, a, b in zip(coefficient, first, second)]
        answer.append(coefficient)
    return answer


class CenteredOrder24Geometry:
    """One fixed128-cell, degree4 certificate; no adaptive refinement loop."""
    def __init__(self, config, *, precision=40):
        self.precision = precision
        inherited = IncomingVacuumTailBoundGeometry(config, intervals=16, precision=precision)
        self.horizon, self.weight = inherited.horizon, inherited.weight
        self.cells = []
        with _precision(precision):
            for i in range(128):
                t = _range(mp.mpf(i)/128, mp.mpf(i+1)/128)
                rho = 1+(self.horizon-1)*(1-t*t)
                center = mp.iv.mpf((_lo(rho)+_hi(rho))/2)
                if not _lo(rho) <= _lo(center) <= _hi(center) <= _hi(rho):
                    raise ArithmeticError('Taylor center outside cell')
                self.cells.append((rho, center))

    def bound_channel(self, channel, *, progress=None):
        if channel.get('index') != 14:
            raise ValueError('this centered order24 certificate owns group14 only')
        with _precision(self.precision):
            _, _, mass, angular, _ = _parameters(channel)
            totals = {s: [mp.iv.mpf(0) for _ in range(25)] for s in (1, -1)}
            natural = {s: [mp.iv.mpf(0) for _ in range(25)] for s in (1, -1)}
            for index, (rho, center) in enumerate(self.cells):
                geometry = _geometry_jets(rho, 28)
                centered = _geometry_jets(center, 28)
                for sign in (1, -1):
                    whole = trimmed_defect_jets(*geometry, mass, sign*angular)
                    point = trimmed_defect_jets(*centered, mass, sign*angular)
                    for j, (w, p) in enumerate(zip(whole, point)):
                        enclosure = centered_coefficient_enclosure(w, p, rho-center, 4)
                        totals[sign][j] += self.weight*mp.iv.mpf(_hi(abs(enclosure)))/128
                        natural[sign][j] += self.weight*mp.iv.mpf(_hi(abs(w[0])))/128
                if progress is not None and (index+1)%16 == 0:
                    progress(index+1)
            return {
                'group': 14, 'lower': 160., 'intervals': 128,
                'physical_Riccati_order': 24, 'centered_remainder_depth': 4,
                'first_inverse_power': 24, 'last_inverse_power': 48,
                'directed_precision': self.precision,
                'per_sign': [{'sign': s,
                              'coefficient_integral_upper': [_up_float(v) for v in totals[s]],
                              'natural_coefficient_integral_upper': [_up_float(v) for v in natural[s]]}
                             for s in (1, -1)],
                'endpoint_weight_upper': _up_float(self.weight),
                'horizon_enclosure': [mp.nstr(_lo(self.horizon), 55), mp.nstr(_hi(self.horizon), 55)],
                'initial_projector_difference': 0.,
                'initial_scope': 'same declared affine-horizon vacuum; unchanged explicit order24 projector',
                'finite_offset_archive_accuracy_certified': False,
                'physical_state_or_approximation_changed': False,
                'bound_method': 'trimmed normalized jets, centered degree4 interval remainder intersected with natural rectangles',
            }
