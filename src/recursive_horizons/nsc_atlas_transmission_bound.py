"""Exact comparison algebra for the connected continuum scalar problem.

The caller supplies proved uniform coefficient/residual bounds on one common
geometry. This does not estimate them from samples or certify metric feedback.
See the connected scalar characteristic inverse in the finite-amplitude owner.
"""
from collections.abc import Mapping
from fractions import Fraction as Q


COMPONENTS = ('cap_outgoing', 'wall_value', 'exterior_transverse',
              'exterior_value', 'exterior_longitudinal', 'cap_incoming',
              'cap_value')
COEFFICIENTS = ('width', 'exterior_potential', 'cap_reaction', 'cap_drift',
                'cap_speed', 'wall_speed', 'incoming_gain', 'outgoing_gain',
                'incoming_value_gain', 'outgoing_value_gain')
N = len(COMPONENTS)


def _rational(value, name, *, positive=False):
    # Float conversion does not turn a measured maximum into an upper bound.
    if isinstance(value, bool) or not isinstance(value, (int, str, Q)):
        raise ValueError(name+' requires an exact integer, rational or rational string')
    try:
        result = Q(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(name+' requires a finite rational') from error
    if result < 0 or (positive and result == 0):
        raise ValueError(name+' is outside its nonnegative/positive domain')
    return result


def _identity():
    return [[Q(i == j) for j in range(N)] for i in range(N)]


def canonical_source_error_bound(*, alpha, radius_floor, u_bound, c_x_bound,
                                r_x_bound, errors, budget):
    """Convert certified invariant source bounds to canonical value/first jet.

    The common segment joining exact and approximate inputs must have r above
    radius_floor, |c/r|<=u_bound and the two supplied first-jet upper bounds.
    No second source derivative or floating elementary function is evaluated.
    """
    if not callable(getattr(budget, 'check', None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(512)
    al = _rational(alpha, 'alpha')
    r = _rational(radius_floor, 'radius floor', positive=True)
    u = _rational(u_bound, 'u bound')
    cx = _rational(c_x_bound, 'c_x bound')
    rx = _rational(r_x_bound, 'r_x bound')
    if not isinstance(errors, Mapping) or set(errors) != {'r','c','r_x','c_x'}:
        raise ValueError('require four invariant source error bounds')
    e = {k:_rational(v, k+' error') for k,v in errors.items()}
    j = al*u**3/18
    derivative = al*u*(cx+u*rx)/(6*r)
    budget.check()
    return dict(value=e['c']+j*e['r'],
                first=e['c_x']+j*e['r_x']+derivative*(e['c']+u*e['r']),
                frame_gap_bound=j, combined_derivative_bound=derivative,
                scope='canonical source representation error on a supplied common tube')


def _product(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(N)), Q(0))
             for j in range(N)] for i in range(N)]


def scalar_inverse_bound(coefficients, data, *, duration, budget, weights=None):
    """Bound seven running suprema from u <= d+B*u+integral(A*u).

    All entries and returned bounds are exact Fractions. The ten coefficients
    are uniform upper bounds, not field values. data contains nonnegative
    upper bounds for initial/source/boundary errors and integrated residuals.
    weights optionally scales the seven error components. A failed short-time
    criterion raises; it is not evidence of an unstable physical solution.
    """
    if not callable(getattr(budget, 'check', None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(5000)
    if not isinstance(coefficients, Mapping) or set(coefficients) != set(COEFFICIENTS):
        raise ValueError('require exactly the ten coefficient bounds')
    if not isinstance(data, Mapping) or set(data) != set(COMPONENTS):
        raise ValueError('require exactly the seven forcing bounds')
    c = {k: _rational(v, k) for k, v in coefficients.items()}
    d = [_rational(data[k], k) for k in COMPONENTS]
    t = _rational(duration, 'duration')
    if weights is None:
        scales = [Q(1)]*N
    else:
        if not isinstance(weights, Mapping) or set(weights) != set(COMPONENTS):
            raise ValueError('require exactly seven positive component scales')
        scales = [_rational(weights[k], k+' scale', positive=True) for k in COMPONENTS]
    B = [[Q(0) for _ in range(N)] for _ in range(N)]
    A = [[Q(0) for _ in range(N)] for _ in range(N)]
    B[2][0], B[2][1] = c['outgoing_gain'], c['outgoing_value_gain']
    B[3][2] = c['width']
    B[4][3] = c['width']*c['exterior_potential']
    B[5][4], B[5][1], B[6][1] = c['incoming_gain'], c['incoming_value_gain'], Q(1)
    A[0][0], A[0][6] = c['cap_drift'], c['cap_reaction']
    A[1][2], A[1][4] = c['wall_speed'], Q(1)
    A[2][3] = c['exterior_potential']
    A[5][5], A[5][6] = c['cap_drift'], c['cap_reaction']
    A[6][0] = A[6][5] = c['cap_speed']/2
    G = _identity()
    power = _identity()
    for _ in range(4):
        power = _product(power, B)
        G = [[G[i][j]+power[i][j] for j in range(N)] for i in range(N)]
    if any(any(row) for row in _product(power, B)):
        raise AssertionError('instantaneous transmission graph is no longer acyclic')
    K = _product(G, A)
    forcing = [sum((G[i][j]*d[j] for j in range(N)), Q(0)) for i in range(N)]
    rate = max(sum((K[i][j]*scales[j] for j in range(N)), Q(0))/scales[i]
               for i in range(N))
    theta = t*rate
    if theta >= 1:
        raise ValueError('short-time inverse margin is not positive; no bound established')
    size = max(forcing[i]/scales[i] for i in range(N))/(1-theta)
    budget.check()
    return dict(instantaneous=B, evolution=A, inverse=G, propagated=K,
                forcing=dict(zip(COMPONENTS, forcing)), rate=rate, theta=theta,
                bounds=dict(zip(COMPONENTS, (size*w for w in scales))),
                scope='conditional continuum scalar inverse; metric feedback and coefficient bounds external')
