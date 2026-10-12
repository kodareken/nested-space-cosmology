"""Centered interval enclosure of the unchanged order16 Riccati defect.

Taylor depth controls the enclosure only. The physical projector still has
exactly sixteen inverse-energy coefficients. Existing radial weights and
positive finite-energy primitives retain ownership of the final bound.
"""
import mpmath as mp

from .nsc_incoming_vacuum_tail_bound import (
    IncomingVacuumTailBoundGeometry, _precision, _lo, _hi, _range,
    _geometry_jets, _mul, _derivative, _up_float,
    recurrence_cancellation_identity,
)


def defect_coefficient_jets(A, sphere, mass, angular, *, depth):
    """f16,...,f32 jets with coefficients f^(j)/j!, j=0,...,depth.

    This is the existing Riccati recurrence with its derivative bookkeeping
    extended. No c17 term is inserted into the approximate mode.
    """
    if depth < 0 or len(A) != len(sphere) or len(A) < 17+depth:
        raise ValueError('geometry jets through degree16+depth required')
    upper = [angular*v for v in sphere]; upper[0] += mp.iv.j*mass
    lower = [mp.iv.mpc(v.real, -v.imag) for v in upper]
    terms = [upper]; derivative_A = _derivative(A)
    zero = lambda: [mp.iv.mpc(0) for _ in A]
    for n in range(1, 16):
        first = _mul(A, _derivative(terms[-1]))
        second = _mul(derivative_A, terms[-1])
        product = zero()
        for j in range(1, n):
            value = _mul(terms[j-1], terms[n-j-1])
            product = [a+b for a, b in zip(product, value)]
        nonlinear = _mul(_mul(A, lower), product)
        terms.append([mp.iv.j*(a+b/2)+c for a, b, c in zip(first, second, nonlinear)])
    answer = []
    for power in range(16, 33):
        square = zero()
        for j in range(1, 17):
            if 1 <= power-j <= 16:
                value = _mul(terms[j-1], terms[power-j-1])
                square = [a+b for a, b in zip(square, value)]
        coefficient = [-v for v in _mul(_mul(A, lower), square)]
        if power == 16:
            first = _mul(A, _derivative(terms[-1]))
            second = _mul(derivative_A, terms[-1])
            coefficient = [c-mp.iv.j*(a+b/2) for c, a, b in zip(coefficient, first, second)]
        answer.append(coefficient[:depth+1])
    return answer


def intersect_complex(left, right):
    """Intersect valid rectangular enclosures, rejecting an empty result."""
    def intersect(a, b):
        lo, hi = max(_lo(a), _lo(b)), min(_hi(a), _hi(b))
        if lo > hi:
            raise ArithmeticError('empty intersection of defect enclosures')
        return _range(lo, hi)
    return mp.iv.mpc(intersect(left.real, right.real), intersect(left.imag, right.imag))


def centered_coefficient_enclosure(interval_jets, center_jets, displacement, depth):
    """Taylor theorem: center coefficients below p, whole-cell remainder p."""
    if not 1 <= depth < len(interval_jets) or len(center_jets) <= depth:
        raise ValueError('center and interval jets through positive remainder degree required')
    polynomial = sum((center_jets[j]*displacement**j for j in range(depth)), mp.iv.mpc(0))
    centered = polynomial+interval_jets[depth]*displacement**depth
    return intersect_complex(centered, interval_jets[0])


class CenteredDefectGeometry:
    """One derivative-enclosure refinement on the existing radial cell map."""
    def __init__(self, config, *, intervals=16, depth=4, precision=40):
        if depth not in (2, 4, 6, 8):
            raise ValueError('bounded centered Taylor depth2,4,6,8 required')
        self.depth, self.intervals, self.precision = depth, intervals, precision
        original = IncomingVacuumTailBoundGeometry(config, intervals=intervals, precision=precision)
        self.weight, self.horizon = original.weight, original.horizon
        self.horizon_derivative, self.kappa = original.horizon_derivative, original.kappa
        self.cells = []
        with _precision(precision):
            for rho, _, _ in original.cells:
                # The midpoint is a selected exact floating number. The
                # displacement interval encloses its difference from ALL rho.
                center = mp.iv.mpf((_lo(rho)+_hi(rho))/2)
                if _lo(center) < _lo(rho) or _hi(center) > _hi(rho):
                    raise ArithmeticError('center is outside its radial enclosure')
                cell_geometry = _geometry_jets(rho, 16+depth)
                center_geometry = _geometry_jets(center, 16+depth)
                self.cells.append((rho, center, cell_geometry, center_geometry))

    def bound_channel(self, channel, lower):
        group = channel.get('index')
        if group not in range(1, 33) or lower != (320. if group in (10, 11, 12, 31, 32) else 160.):
            raise ValueError('existing retained group and archived endpoint required')
        with _precision(self.precision):
            definition = channel['compact_level']*mp.iv.pi/2
            mass = _range(min(_lo(definition), mp.mpf(float(channel['compact_mass']))),
                          max(_hi(definition), mp.mpf(float(channel['compact_mass']))))
            definition = mp.iv.sqrt(channel['angular_level']*(channel['angular_level']+4))
            angular = _range(min(_lo(definition), mp.mpf(float(channel['angular_eigenvalue']))),
                             max(_hi(definition), mp.mpf(float(channel['angular_eigenvalue']))))
            signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
            details = []
            for sign in signs:
                integrals = [mp.iv.mpf(0) for _ in range(17)]
                natural = [mp.iv.mpf(0) for _ in range(17)]
                for rho, center, cell_geometry, center_geometry in self.cells:
                    interval_jets = defect_coefficient_jets(*cell_geometry, mass, sign*angular, depth=self.depth)
                    center_jets = defect_coefficient_jets(*center_geometry, mass, sign*angular, depth=self.depth)
                    for i, (whole, point) in enumerate(zip(interval_jets, center_jets)):
                        enclosure = centered_coefficient_enclosure(whole, point, rho-center, self.depth)
                        integrals[i] += self.weight*mp.iv.mpf(_hi(abs(enclosure)))/self.intervals
                        natural[i] += self.weight*mp.iv.mpf(_hi(abs(whole[0])))/self.intervals
                details.append({'sign': sign, 'I16_to_I32_upper': [_up_float(v) for v in integrals],
                                'natural_I16_to_I32_upper': [_up_float(v) for v in natural]})
            return {'group': group, 'lower': float(lower), 'intervals': self.intervals,
                    'physical_Riccati_order': 16, 'centered_remainder_depth': self.depth,
                    'directed_precision': self.precision, 'per_sign': details,
                    'lower_order_cancellation': recurrence_cancellation_identity(),
                    'horizon_enclosure': [mp.nstr(_lo(self.horizon), 55), mp.nstr(_hi(self.horizon), 55)],
                    'initial_projector_difference': 0.,
                    'initial_scope': 'declared affine-horizon limit; finite-delta input artifacts unchanged',
                    'finite_offset_initial_projector_term_bound': None,
                    'finite_offset_archive_accuracy_certified': False,
                    'physical_state_or_approximation_changed': False,
                    'bound_method': 'centered Taylor theorem with directed whole-cell remainder; intersection with natural interval'}
