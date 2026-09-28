"""Directed-interval defect bound for the declared affine-horizon vacuum.

The bound compares its occupied projector with the existing order16 Riccati
projector on the SAME stationary collar. It uses the imported unitary defect
estimate, an analytic endpoint weight, positive interval upper sums, and exact
energy-power integration. No field ODE, scattering solve or finite-start
state replacement occurs. Frozen finite-delta data remain numerical controls.
"""
from contextlib import contextmanager
from functools import lru_cache
import math

import mpmath as mp
import numpy as np
import sympy as sp


@contextmanager
def _precision(dps):
    old = mp.iv.dps
    mp.iv.dps = dps
    try:
        with mp.workdps(max(90, dps+30)):
            yield
    finally:
        mp.iv.dps = old


def _lo(value): return mp.make_mpf(value._mpi_[0])
def _hi(value): return mp.make_mpf(value._mpi_[1])
def _range(lower, upper): return mp.iv.mpf([lower, upper])
def _up_float(value): return float(np.nextafter(float(_hi(value)), np.inf))
def _atan(value): return mp.iv.atan2(value, mp.iv.mpf(1))


def _A(rho): return 1+3*rho+3*(1+rho*rho)*(_atan(rho)-mp.iv.pi/2)
def _Aprime(rho): return 6+6*rho*(_atan(rho)-mp.iv.pi/2)
def _Ahalfsecond(rho): return 3*(_atan(rho)-mp.iv.pi/2+rho/(1+rho*rho))


def _horizon_enclosure(hint, kappa_hint, dps):
    """Enclose the existing profile root; do not change any stored parameter."""
    h = mp.mpf(float(hint)); left, right = h-mp.mpf('1e-11'), h+mp.mpf('1e-11')
    if _hi(_A(mp.iv.mpf(left))) >= 0 or _lo(_A(mp.iv.mpf(right))) <= 0:
        raise ValueError('declared horizon is outside the bounded numerical root bracket')
    for _ in range(80):
        middle = (left+right)/2; value = _A(mp.iv.mpf(middle))
        if _hi(value) < 0: left = middle
        elif _lo(value) > 0: right = middle
        else: break
    horizon = _range(left, right)
    derivative = _Aprime(horizon)
    if _lo(derivative) <= 0:
        raise ArithmeticError('positive horizon derivative not enclosed')
    kappa = derivative/2
    # State/geometry receipts store binary floats. Report a rounding enclosure
    # including both, and reject any discrepancy beyond numerical rounding.
    supplied = mp.mpf(float(kappa_hint))
    distance = max(mp.mpf(0), _lo(kappa)-supplied, supplied-_hi(kappa))
    if distance > mp.mpf('3e-14'):
        raise ValueError('declared surface gravity does not match the owned profile')
    rounded = _range(min(_lo(kappa), supplied), max(_hi(kappa), supplied))
    return horizon, derivative, rounded


def _mul(a, b):
    return [sum((a[j]*b[n-j] for j in range(n+1)), mp.iv.mpc(0)) for n in range(len(a))]


def _derivative(a):
    return [(n+1)*a[n+1] for n in range(len(a)-1)]+[mp.iv.mpc(0)]


def _geometry_jets(rho, order):
    """Interval Taylor coefficients of the same rational/arctangent profile."""
    size = order+2
    zero = lambda: [mp.iv.mpf(0) for _ in range(size)]
    r2 = zero(); r2[:3] = [1+rho*rho, 2*rho, mp.iv.mpf(1)]
    inverse = zero(); inverse[0] = 1/r2[0]
    for n in range(1, size):
        inverse[n] = -sum((r2[j]*inverse[n-j] for j in range(1, n+1)), mp.iv.mpf(0))/r2[0]
    sphere = zero(); sphere[0] = mp.iv.sqrt(inverse[0])
    for n in range(1, size):
        sphere[n] = (inverse[n]-sum((sphere[j]*sphere[n-j] for j in range(1, n)), mp.iv.mpf(0)))/(2*sphere[0])
    # Monotonic scalar ranges reduce dependency overestimation; higher jets
    # use the exact identity A'''=12/(1+rho^2)^2 of the owned profile.
    left, right = mp.iv.mpf(_lo(rho)), mp.iv.mpf(_hi(rho))
    A = zero()
    A[0] = _range(_lo(_A(left)), _hi(_A(right)))
    A[1] = _range(_lo(_Aprime(right)), _hi(_Aprime(left)))
    A[2] = _range(_lo(_Ahalfsecond(left)), _hi(_Ahalfsecond(right)))
    inverse_squared = _mul(inverse, inverse)
    for n in range(3, size):
        A[n] = 12*inverse_squared[n-3]/(n*(n-1)*(n-2))
    return A, sphere


def _defect_coefficients(A, sphere, mass, angular, order=16):
    """Coefficients f16..f32; lower coefficients cancel by the recurrence."""
    upper = [angular*v for v in sphere]; upper[0] += mp.iv.j*mass
    # Explicit rectangular conjugation avoids ctx_iv's scalar-conjugate
    # implementation, which does not accept its own interval endpoint pair.
    lower = [mp.iv.mpc(v.real, -v.imag) for v in upper]
    terms = [upper]; derivative_A = _derivative(A)
    zero = lambda: [mp.iv.mpc(0) for _ in A]
    for n in range(1, order):
        first, second = _mul(A, _derivative(terms[-1])), _mul(derivative_A, terms[-1])
        product = zero()
        for j in range(1, n):
            value = _mul(terms[j-1], terms[n-j-1])
            product = [a+b for a, b in zip(product, value)]
        nonlinear = _mul(_mul(A, lower), product)
        terms.append([mp.iv.j*(a+b/2)+c for a, b, c in zip(first, second, nonlinear)])
    answer = []
    for power in range(order, 2*order+1):
        square = sum((terms[j-1][0]*terms[power-j-1][0]
                      for j in range(1, order+1) if 1 <= power-j <= order), mp.iv.mpc(0))
        coefficient = -A[0]*lower[0]*square
        if power == order:
            coefficient -= mp.iv.j*(A[0]*terms[-1][1]+A[1]*terms[-1][0]/2)
        answer.append(coefficient)
    return answer


@lru_cache(maxsize=1)
def recurrence_cancellation_identity():
    """Exact scalar check of all fifteen lower coefficient substitutions."""
    A, Ap, U = sp.symbols('A Ap U')
    c = sp.symbols('c1:18'); dc = sp.symbols('dc1:17')
    residuals = []
    for p in range(1, 16):
        rhs = sp.I*(A*dc[p-1]+Ap*c[p-1]/2)+A*U*sum(c[j-1]*c[p-j-1] for j in range(1, p))
        residuals.append(str(sp.expand((c[p]-rhs).subs(c[p], rhs))))
    if residuals != ['0']*15:
        raise ArithmeticError('lower inverse-energy terms did not cancel')
    return residuals


def vacuum_tail_from_coefficients(channel, lower, per_sign, *, precision=40):
    """Replay the bound from outward-rounded radial coefficient integrals.

    This does not construct geometry cells or repeat the interval recurrence.
    Its inputs are upper bounds, not samples or fitted coefficients.
    """
    group = channel.get('index')
    if group not in range(1, 33) or lower != (320. if group in (10, 11, 12, 31, 32) else 160.):
        raise ValueError('existing retained group and archived energy split required')
    expected = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
    if tuple(row['sign'] for row in per_sign) != expected:
        raise ValueError('each actual angular sign must be included exactly once')
    with _precision(precision):
        mass_definition = channel['compact_level']*mp.iv.pi/2
        angular_definition = mp.iv.sqrt(channel['angular_level']*(channel['angular_level']+4))
        mass = _range(min(_lo(mass_definition), mp.mpf(float(channel['compact_mass']))),
                      max(_hi(mass_definition), mp.mpf(float(channel['compact_mass']))))
        angular = _range(min(_lo(angular_definition), mp.mpf(float(channel['angular_eigenvalue']))),
                         max(_hi(angular_definition), mp.mpf(float(channel['angular_eigenvalue']))))
        factor = (channel['copy_count']*channel['degeneracy']
                  /(4*mp.iv.pi**2*2*mp.iv.sqrt(3*mp.iv.pi/2-4)*len(expected)))
        a1 = mp.iv.sqrt(3*mp.iv.pi/2-4); r1 = mp.iv.sqrt(2)
        M1 = mp.iv.sqrt(mass*mass+(angular/r1)**2)
        total = [mp.iv.mpf(0) for _ in range(4)]
        for row in per_sign:
            values = row['I16_to_I32_upper']
            if len(values) != 17 or any(not math.isfinite(v) or v < 0 for v in values):
                raise ValueError('seventeen finite nonnegative directed integral bounds required')
            bounds = [mp.iv.mpf(0) for _ in range(4)]
            for power, integral in zip(range(16, 33), values):
                primitive0 = mp.iv.mpf(lower)**(1-power)/(power-1)
                primitive1 = mp.iv.mpf(lower)**(2-power)/(power-2)
                common = mp.iv.mpf(integral)/(mp.iv.mpf(2)**power)
                bounds[0] += 2*common*(primitive1/a1+M1*primitive0)
                bounds[1] += 2*common*primitive1/a1
                bounds[3] += common*angular/r1*primitive0
            for i in range(4): total[i] += factor*bounds[i]
        return [_up_float(v) if i != 2 else 0. for i, v in enumerate(total)]


def aggregate_vacuum_tail_bounds(groups, *, precision=40):
    """Outward-rounded aggregate; the budget belongs to all32 groups."""
    if [row['group'] for row in groups] != list(range(1, 33)):
        raise ValueError('all32 retained non-LLL groups in owned order required')
    if any(row['initial_scope'] != 'declared affine-horizon limit; finite-delta input artifacts unchanged'
           or row['initial_projector_difference'] != 0.
           or row['finite_offset_initial_projector_term_bound'] is not None
           or row['finite_offset_archive_accuracy_certified'] for row in groups):
        raise ValueError('affine-horizon bound cannot silently certify a finite-offset archive')
    with _precision(precision):
        total = [sum((mp.iv.mpf(row['vacuum_tail_error_upper'][i]) for row in groups), mp.iv.mpf(0))
                 for i in range(4)]
        upper = [_up_float(v) if i != 2 else 0. for i, v in enumerate(total)]
    budget = 32*3e-12
    return {'vacuum_tail_error_upper': upper, 'component_budget': budget,
            'group_component_budget': 3e-12, 'group_count': 32,
            'within_component_budget': max(upper) <= budget,
            'finite_offset_initial_projector_term_bound': None,
            'finite_offset_archive_accuracy_certified': False}


class IncomingVacuumTailBoundGeometry:
    """Cache interval geometry once for actual retained channels/signs."""
    def __init__(self, config, *, intervals=8, precision=40):
        if intervals not in (8, 16) or precision < 30:
            raise ValueError('bounded partition8 or16 and directed precision>=30 required')
        if config.get('magnetic_flux') != 4:
            raise ValueError('the locked magnetic q=4 branch is required')
        self.intervals, self.precision = intervals, precision
        with _precision(precision):
            self.horizon, self.horizon_derivative, self.kappa = _horizon_enclosure(
                config['horizon_rho'], config['surface_gravity'], precision)
            # A''<0 for rho>0, so A'(rho)>=A'(rho_h)>0. For
            # rho=1+(rho_h-1)(1-t^2), |d rho|/a <=weight dt globally.
            alpha = mp.iv.mpf(_lo(self.horizon_derivative))
            self.weight = 2*mp.iv.sqrt((self.horizon-1)/alpha)
            self.cells = []
            for i in range(intervals):
                # Form the cell map with directed operations as well. Scalar
                # endpoint arithmetic, even at higher precision, is not an
                # enclosure and can leave tiny gaps between adjacent cells.
                t = _range(mp.mpf(i)/intervals, mp.mpf(i+1)/intervals)
                rho = 1+(self.horizon-1)*(1-t*t)
                self.cells.append((rho, *_geometry_jets(rho, 16)))

    def bound_channel(self, channel, lower):
        group = channel.get('index')
        if group not in range(1, 33) or lower != (320. if group in (10, 11, 12, 31, 32) else 160.):
            raise ValueError('existing retained group and archived energy split required')
        with _precision(self.precision):
            signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
            # Enclose the stored label's binary rounding and its already-owned
            # compact/magnetic definition; do not choose a new parameter.
            mass_definition = channel['compact_level']*mp.iv.pi/2
            angular_definition = mp.iv.sqrt(channel['angular_level']*(channel['angular_level']+4))
            mass = _range(min(_lo(mass_definition), mp.mpf(float(channel['compact_mass']))),
                          max(_hi(mass_definition), mp.mpf(float(channel['compact_mass']))))
            angular = _range(min(_lo(angular_definition), mp.mpf(float(channel['angular_eigenvalue']))),
                             max(_hi(angular_definition), mp.mpf(float(channel['angular_eigenvalue']))))
            details = []
            for sign in signs:
                integrals = [mp.iv.mpf(0) for _ in range(17)]
                for _, A, sphere in self.cells:
                    coefficients = _defect_coefficients(A, sphere, mass, sign*angular)
                    for i, value in enumerate(coefficients):
                        # A positive directed upper sum, not a quadrature
                        # sample of the transport defect.
                        integrals[i] += self.weight*mp.iv.mpf(_hi(abs(value)))/self.intervals
                details.append({'sign': sign, 'I16_to_I32_upper': [_up_float(v) for v in integrals]})
            identity = recurrence_cancellation_identity()
            # Store the actual outward-rounded coefficient inputs used in
            # this final contraction, so verification replays exactly.
            upper = vacuum_tail_from_coefficients(channel, lower, details, precision=self.precision)
            return {'group': group, 'lower': float(lower), 'intervals': self.intervals,
                'directed_precision': self.precision, 'vacuum_tail_error_upper': upper,
                'per_sign': details, 'lower_order_cancellation': identity,
                'horizon_enclosure': [mp.nstr(_lo(self.horizon), 55), mp.nstr(_hi(self.horizon), 55)],
                'kappa_rounding_enclosure': [mp.nstr(_lo(self.kappa), 55), mp.nstr(_hi(self.kappa), 55)],
                'endpoint_weight_upper': _up_float(self.weight),
                'bound_method': 'directed coefficient enclosures and positive interval upper sums, then analytic energy-power integration',
                'initial_projector_difference': 0.,
                'initial_scope': 'declared affine-horizon limit; finite-delta input artifacts unchanged',
                'finite_offset_initial_projector_term_bound': None,
                'finite_offset_archive_accuracy_certified': False,
                'physical_scope': 'vacuum projector on the fixed incoming collar, same operator and horizon state',
                'state_or_metric_evolution': False,
                'status': 'PASS' if max(upper) < 3e-12 else 'OPEN: valid enclosure is too broad for the component budget',
                'component_budget': 3e-12}
