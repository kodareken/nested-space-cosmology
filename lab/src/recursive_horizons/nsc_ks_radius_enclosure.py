"""Exact bounds for the finite-history reciprocal-radius potential remainder."""
from fractions import Fraction as Q

from .nsc_ks_residual_error import nonnegative
from .nsc_local_incoming_family import LocalAxialFunction


def axial_profile_bounds(profile):
    """Sup bounds through two derivatives of the represented analytic profile.

    Binary Chebyshev-map coefficients are treated exactly. This bounds the
    analytic function; floating evaluation errors still require accounting.
    Power-coefficient absolute sums are deliberately conservative.
    """
    if not isinstance(profile, LocalAxialFunction):
        raise TypeError('owned LocalAxialFunction required')
    off, scale = map(Q, profile._derivatives[0].mapparms())
    center, outer, inner = map(Q, (profile.center, profile.outer, profile.inner))
    M = max(abs(off + scale*(center-outer)), abs(off + scale*(center+outer)))
    basis = [[Q(1)], [Q(0), Q(1)]]
    for n in range(2, len(profile.coefficients)):
        row = [Q(0)] + [2*x for x in basis[-1]]
        for j, value in enumerate(basis[-2]):
            row[j] -= value
        basis.append(row)
    polynomial = [Q(0)] * len(profile.coefficients)
    for coefficient, row in zip(profile.coefficients, basis):
        for j, value in enumerate(row):
            polynomial[j] += Q(coefficient) * value
    norms = []
    for derivative in range(3):
        norms.append(abs(scale)**derivative * sum(abs(c)*M**j for j, c in enumerate(polynomial)))
        polynomial = [j*c for j, c in enumerate(polynomial)][1:]
    p0, p1, p2 = norms
    width = outer - inner
    # For eta=expit(1/u-1/(1-u)), |eta'|<=8, |eta''|<=105.
    return p0, p1 + 8*p0/width, p2 + 16*p1/width + 105*p0/width**2


def reciprocal_radius_tail(w_norms, U_norms, *, normal_support,
                           reference_radius_lower, axial_lower, absolute_angular, order):
    """Sup bounds for V-V_order and its first two axial derivatives.

    V=(i ell/a)(1/(r_ref+delta r)-1/r_ref) sigma2. V_order uses
    powers 1..order of -delta r/r_ref. No approximation is made to the
    physical generator; this function bounds the validation remainder.
    """
    if isinstance(order, bool) or not isinstance(order, int) or order < 1:
        raise ValueError('positive integer series order required')
    if len(w_norms) != 3 or len(U_norms) != 3:
        raise ValueError('three derivative bounds required')
    W, U = (tuple(nonnegative(v, name) for v in values)
            for name, values in (('w norms', w_norms), ('U norms', U_norms)))
    s, r, a, ell = [nonnegative(v, name) for v, name in (
        (normal_support, 'normal support'), (reference_radius_lower, 'reference radius'),
        (axial_lower, 'axial scale'), (absolute_angular, 'angular label'))]
    if min(s, r, a) <= 0:
        raise ValueError('positive support and geometry lower bounds required')
    eta, eta1, eta2 = ((s*w + s**3*u/6)/r for w, u in zip(W, U))
    if eta >= 1:
        raise ValueError('reciprocal series requires a strict radius margin eta<1')
    den, n = 1-eta, order+1
    f0 = eta**n/den
    f1 = n*eta**(n-1)/den + eta**n/den**2
    f2 = n*(n-1)*eta**(n-2)/den + 2*n*eta**(n-1)/den**2 + 2*eta**n/den**3
    factor = ell/(a*r)
    return {'radius_ratio_upper': eta, 'radius_lower': r*(1-eta),
            'potential_derivative_tail_bounds': (factor*f0, factor*f1*eta1,
                                                 factor*(f2*eta1**2+f1*eta2)),
            'order': order,
            'scope': 'exact conditional geometry remainder; no field or full-source certificate'}
