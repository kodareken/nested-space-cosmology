"""Interval Chebyshev profile jets by midpoint evaluation and endpoint remainder.

This owner generalizes the Chebyshev endpoint recurrence of
`nsc_ks_current_history_bounds` to every derivative of the represented
polynomial. Production evaluation, the hash-bound module, and the v2
field adapter are unchanged.

Mathematical basis
------------------
Chebyshev polynomials of the first kind obey the three-term recurrence
(DLMF 18.9.1, Table 18.9.1)

    T_0(x) = 1,  T_1(x) = x,
    T_{n+1}(x) = 2 x T_n(x) - T_{n-1}(x).

Differentiating j times and applying Leibniz to the 2x T_n term gives

    T_{n+1}^{(j)}(x) = 2 j T_n^{(j-1)}(x) + 2 x T_n^{(j)}(x) - T_{n-1}^{(j)}(x).

For M >= 1 and |x| <= M, |T_n^{(j)}(x)| <= T_n^{(j)}(M). On [-1, 1] this
is the Chebyshev endpoint bound: DLMF 18.9.21 has T_n' = n U_{n-1}, and
U_{n-1}(1) = n, so T_n'(1) = n^2. Higher endpoint values follow the same
recurrence; Markov's inequality identifies T_n as the extremal degree-n
polynomial. For x >= 1, T_n(x) = cosh(n arcosh x) (DLMF 18.5.1 continued
off the cut), so T_n^{(j)} is positive and increasing on [1, ∞). Parity
covers [-M, -1]. Ultraspherical derivatives (DLMF 18.9.19)
d/dx C_n^{(λ)} = 2 λ C_{n-1}^{(λ+1)} with U_n = C_n^{(1)} (DLMF 18.7.4)
recover the same first-derivative identity. For j>=1, iterating these
identities gives a positive multiple of C_{n-j}^{(j)}; the endpoint
inequality in DLMF 18.14.4 supplies the bound throughout [-1,1].

The axial map is the exact binary Chebyshev domain map
x = offset + scale z from numpy `mapparms()`, retained as dyadic
rationals. Then P^{(j)}(z) = scale^j sum_n a_n T_n^{(j)}(x(z)), and

    |P^{(j)}(z)| <= |scale|^j sum_n |a_n| T_n^{(j)}(M)

on the mapped interval. Expanding a high-degree profile into monomials
and summing absolute power coefficients is a different, looser bound.

On a real interval with midpoint c and radius r, two valid enclosures of
P^{(j)} are formed and intersected:

1. Mean-value remainder (parent form): P^{(j)}(c) + [-r, r] B_{j+1},
   with B_{j+1} the global endpoint bound of |P^{(j+1)}|.
2. Exact finite Taylor of a degree-N polynomial,
   P^{(j)}(c+h) = sum_{m=0}^{N-j} P^{(j+m)}(c) h^m / m!,
   plus the Lagrange remainder r^k / k! B_{j+k} for any truncation
   length k. The full sum is an identity, not an extra hypothesis.
   Absolute sums of the Taylor terms avoid Horner wrapping in h.

The endpoint/Lagrange intersections are used only when the mapped input
ball lies inside their declared [-M,M]. Outside that domain, including
intervals inflated slightly beyond a cutoff endpoint by ball arithmetic,
only the domain-independent complete Taylor enclosure is used.

Interval Clenshaw of T_n on a broad x-ball is a valid natural interval
extension, but dependency explosion makes it unusable for Chebyshev-32
eighth derivatives. Midpoint evaluation has no such wrapping.

The owned plateau cutoff jet (`plateau_series`) multiplies the
polynomial jet through the Leibniz / `arb_series` product, preserving
the Taylor convention f^{(n)}/n! and the declared inner/outer domains.
On the plateau interior the cutoff is identically 1.

The represented (w, U) coefficient vectors are the same analytic
geometry as `one_direction_metric(family)`. Summing 64 separate
broad-interval basis enclosures is not used.

This is a continuous proof method for compact profile derivatives. It
does not close the local incoming gate and does not bound source or
whole-field error.
"""
from fractions import Fraction as Q
from math import factorial

from flint import acb, arb, arb_series, ctx

from . import nsc_ks_ball_geometry as G
from .nsc_ks_profile_fourier_bound import ProfileFourierBound, alias_and_tail_bounds
from .nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily


MATHEMATICAL_JUSTIFICATION = __doc__
DEFAULT_MAX_ENDPOINT_ORDER = 9
DEFAULT_TRANSITION_PANELS = 16
DEFAULT_INTERIOR_PANELS = 4
DEFAULT_DERIVATIVE_ORDER = 8


class _PolyWork:
    """Raise flint series cap for a polynomial jet; restore on exit."""

    def __init__(self, order):
        if isinstance(order, bool) or not isinstance(order, int) or order < 0:
            raise ValueError('nonnegative polynomial jet order required')
        self.order = order

    def __enter__(self):
        self._cap, self._prec = ctx.cap, ctx.prec
        ctx.cap = self.order + 1
        ctx.prec = max(self._prec, 80)
        return self.order

    def __exit__(self, *exc):
        ctx.cap, ctx.prec = self._cap, self._prec


def binary_domain_map(profile):
    """Exact dyadic Chebyshev map x = offset + scale z from mapparms()."""
    if not isinstance(profile, LocalAxialFunction):
        raise TypeError('owned LocalAxialFunction required')
    offset, scale = map(Q, profile._derivatives[0].mapparms())
    center, inner, outer = map(Q, (profile.center, profile.inner, profile.outer))
    if outer <= inner:
        raise ValueError('positive axial transition width required')
    return offset, scale, center, inner, outer


def _endpoint_M(offset, scale, center, outer):
    return max(Q(1), abs(offset + scale * (center - outer)),
               abs(offset + scale * (center + outer)))


def _chebyshev_endpoint_table(degree, M, max_order):
    """T_n^{(j)}(M) for n <= degree, j <= max_order, as exact rationals.

    Recurrence: T_{n+1}^{(j)} = 2 j T_n^{(j-1)} + 2 M T_n^{(j)} - T_{n-1}^{(j)},
    from differentiating DLMF 18.9.1. All entries are nonnegative for M >= 1.
    """
    if M < 1:
        raise ValueError('Chebyshev endpoint bound requires M >= 1')
    zero = tuple(Q(0) for _ in range(max_order + 1))
    rows = [tuple(Q(1) if j == 0 else Q(0) for j in range(max_order + 1))]
    if degree >= 1:
        first = [Q(0)] * (max_order + 1)
        first[0] = M
        if max_order >= 1:
            first[1] = Q(1)
        rows.append(tuple(first))
    for _ in range(2, degree + 1):
        current, previous = rows[-1], rows[-2]
        nxt = []
        for j in range(max_order + 1):
            value = 2 * M * current[j] - previous[j]
            if j:
                value += (2 * j) * current[j - 1]
            nxt.append(value)
        rows.append(tuple(nxt))
    if any(value < 0 for row in rows for value in row):
        raise ArithmeticError('Chebyshev endpoint derivative became negative')
    return tuple(rows)


def chebyshev_endpoint_derivative_bounds(profile, max_order=DEFAULT_MAX_ENDPOINT_ORDER):
    """Global sup bounds of the polynomial (no cutoff) through max_order.

    polynomial_bounds[j] = |scale|^j sum |a_n| T_n^{(j)}(M). Orders 0, 1, 2
    match `chebyshev_profile_bounds` polynomial_bounds on the same profile.
    """
    if isinstance(max_order, bool) or not isinstance(max_order, int) or max_order < 0:
        raise ValueError('nonnegative endpoint derivative order required')
    offset, scale, center, inner, outer = binary_domain_map(profile)
    degree = len(profile.coefficients) - 1
    order = max(max_order, degree)
    M = _endpoint_M(offset, scale, center, outer)
    table = _chebyshev_endpoint_table(degree, M, order)
    bounds = []
    abs_scale = abs(scale)
    for j in range(order + 1):
        total = sum((abs(Q(coefficient)) * row[j]
                     for coefficient, row in zip(profile.coefficients, table)), Q())
        bounds.append((abs_scale ** j) * total)
    return {
        'map_endpoint_upper': M,
        'map_offset': offset,
        'map_scale': scale,
        'polynomial_bounds': tuple(bounds),
        'max_order': order,
    }


def _arb_rational(value):
    return arb(value.numerator) / arb(value.denominator)


def _symmetric(center, radius):
    radius = arb(radius)
    if not radius >= 0:
        raise ValueError('nonnegative remainder radius required')
    return center + radius * arb(0, 1)


def _cover_radius(z_ball, center):
    return abs(z_ball - center).upper()


def _polynomial_point_jet(profile, z_mid, order):
    """Taylor coefficients P^{(k)}(z_mid)/k! by Clenshaw AD at a point."""
    offset, scale = profile._derivatives[0].mapparms()
    with _PolyWork(order):
        z = G._identity(G._binary_arb(z_mid, 'z midpoint'), order)
        x = G._binary_arb(float(offset), 'mapparms off') + G._binary_arb(
            float(scale), 'mapparms scale') * z
        return G._finite_series(G._chebyshev_series(profile.coefficients, x, order),
                                'Chebyshev midpoint', order)


def _derivative_from_point_jet(coeffs, derivative_order, radius, global_bounds,
                               *, endpoint_valid=True):
    """P^{(j)}(c+h) for |h|<=r from a midpoint jet and endpoint remainders.

    Truncation k=1 is the parent mean-value form P^{(j)}(c)+[-r,r] B_{j+1}.
    Truncation k=N-j+1 is the exact finite Taylor identity. Every k is a
    valid enclosure; the returned ball is their intersection.
    """
    degree = len(coeffs) - 1
    if derivative_order > degree:
        return arb(0)
    point = factorial(derivative_order) * coeffs[derivative_order]
    remaining = degree - derivative_order
    if not radius > 0:
        return point

    def enclosure_for(lagrange_order):
        extra = arb(0)
        power = radius
        for m in range(1, lagrange_order):
            extra += abs(factorial(derivative_order + m) * coeffs[derivative_order + m]
                         / factorial(m)) * power
            power *= radius
        global_index = derivative_order + lagrange_order
        if global_index <= degree:
            extra += _arb_rational(global_bounds[global_index]) * power / factorial(lagrange_order)
        return _symmetric(point, extra)

    # The endpoint majorants cover only their declared mapped interval.
    # The complete finite Taylor sum remains valid on any finite interval.
    # Intersecting that sum with a majorant outside its domain can remove
    # the true value, even for T_2 at z=2.
    if not endpoint_valid:
        return enclosure_for(remaining + 1)
    result = enclosure_for(1)
    for lagrange_order in range(2, remaining + 2):
        candidate = enclosure_for(lagrange_order)
        if not result.overlaps(candidate):
            raise ArithmeticError('polynomial midpoint and endpoint remainders disagree')
        result = result.intersection(candidate)
    if not result.is_finite():
        raise ArithmeticError('nonfinite polynomial derivative enclosure')
    return result


def _endpoint_domain_covers(profile, z_ball, endpoint_bounds):
    offset, scale, *_ = binary_domain_map(profile)
    mapped = _arb_rational(offset) + _arb_rational(scale)*z_ball
    return bool(abs(mapped).upper() <= _arb_rational(endpoint_bounds['map_endpoint_upper']))


def polynomial_derivative_on_interval(profile, z_ball, derivative_order, endpoint_bounds=None):
    """Enclosure of P^{(j)} on a real z-interval."""
    if not isinstance(profile, LocalAxialFunction):
        raise TypeError('owned LocalAxialFunction required')
    if (isinstance(derivative_order, bool) or not isinstance(derivative_order, int)
            or derivative_order < 0):
        raise ValueError('nonnegative polynomial derivative order required')
    z_ball = G._binary_arb(z_ball, 'z')
    degree = len(profile.coefficients) - 1
    if endpoint_bounds is None:
        endpoint_bounds = chebyshev_endpoint_derivative_bounds(profile, max(derivative_order + 1, degree))
    center = z_ball.mid()
    radius = _cover_radius(z_ball, center)
    jet = _polynomial_point_jet(profile, center, degree)
    return _derivative_from_point_jet(
        G._coeffs(jet, degree), derivative_order, radius, endpoint_bounds['polynomial_bounds'],
        endpoint_valid=_endpoint_domain_covers(profile, z_ball, endpoint_bounds))


def polynomial_range_series(profile, z_ball, order, endpoint_bounds=None):
    """arb_series whose k-th coefficient encloses P^{(k)}(I)/k! on the ball."""
    z_ball = G._binary_arb(z_ball, 'z')
    degree = len(profile.coefficients) - 1
    if endpoint_bounds is None:
        endpoint_bounds = chebyshev_endpoint_derivative_bounds(profile, max(order + 1, degree))
    jet = _polynomial_point_jet(profile, z_ball.mid(), degree)
    coeffs = G._coeffs(jet, degree)
    radius = _cover_radius(z_ball, z_ball.mid())
    covered = _endpoint_domain_covers(profile, z_ball, endpoint_bounds)
    terms = [_derivative_from_point_jet(coeffs, j, radius, endpoint_bounds['polynomial_bounds'],
                                      endpoint_valid=covered) / factorial(j)
             for j in range(order + 1)]
    return arb_series(terms, prec=order + 1)


def local_axial_interval_series(profile, z_ball, order, endpoint_bounds=None):
    """Jet of χ P on a real z-ball, polynomial part by midpoint Taylor.

    Cutoff jets and domain classification are the owned `plateau_series`
    path. Derivative factorials stay in the arb_series convention.
    """
    if not isinstance(profile, LocalAxialFunction):
        raise TypeError('owned LocalAxialFunction required')
    if endpoint_bounds is None:
        endpoint_bounds = chebyshev_endpoint_derivative_bounds(profile, max(order + 1, 9))
    with G._JetWork(order) as n:
        z0 = G._binary_arb(z_ball, 'z')
        center = G._binary_arb(profile.center, 'center')
        inner = G._binary_arb(profile.inner, 'inner')
        if abs(z0 - center).upper() <= inner:
            window = G._constant_series(arb(1), n)
        elif z0 >= center:
            window = G.plateau_series(G._identity(z0, n) - center,
                                      profile.inner, profile.outer, n)
        elif z0 <= center:
            window = G.plateau_series(center - G._identity(z0, n),
                                      profile.inner, profile.outer, n)
        else:
            raise G.SubdivisionNeeded(
                'z ball straddles the plateau center; |z-center| is not analytic')
        poly = polynomial_range_series(profile, z0, n, endpoint_bounds)
        return G._finite_series(window * poly, 'interval LocalAxialFunction', n)


def actual_profile_functions(family):
    """The represented (w, U) vectors; same analytic geometry as one_direction_metric."""
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned local incoming family required')
    return family.functions


def _common_window(family):
    w, U = actual_profile_functions(family)
    windows = {(float(p.center), float(p.inner), float(p.outer)) for p in (w, U)}
    if len(windows) != 1:
        raise ValueError('the approved common axial cutoff is required')
    return next(iter(windows))


def monomial_interval_series(family, z_ball, keys, order, endpoint_cache=None):
    """Jets of w^p U^q on a z-ball from the actual 32-coefficient profiles."""
    w_profile, u_profile = actual_profile_functions(family)
    if endpoint_cache is None:
        endpoint_cache = (
            chebyshev_endpoint_derivative_bounds(w_profile, max(order + 1, 9)),
            chebyshev_endpoint_derivative_bounds(u_profile, max(order + 1, 9)),
        )
    with G._JetWork(order):
        w = local_axial_interval_series(w_profile, z_ball, order, endpoint_cache[0])
        U = local_axial_interval_series(u_profile, z_ball, order, endpoint_cache[1])
        results = {}
        for key in keys:
            series = G._constant_series(arb(1), order)
            if key[0]:
                series *= w ** key[0]
            if key[1]:
                series *= U ** key[1]
            results[key] = G._finite_series(series, 'monomial', order)
        return results


def monomial_derivative_l1(family, keys, order=DEFAULT_DERIVATIVE_ORDER,
                           transition_panels=DEFAULT_TRANSITION_PANELS,
                           interior_panels=DEFAULT_INTERIOR_PANELS, bits=120):
    """L1 bounds of (w^p U^q)^{(order)} over the compact cutoff support.

    Panel layout matches the v2 Fourier owner: `transition_panels` on each
    cutoff ramp and `interior_panels` on the plateau. The derivative on each
    cell is the midpoint/endpoint jet, not interval Clenshaw of the basis.
    """
    if (not keys or len(set(keys)) != len(keys)
            or any(len(k) != 2 or min(k) < 0 or sum(k) < 1 for k in keys)):
        raise ValueError('distinct positive-total-degree profile powers required')
    if (isinstance(transition_panels, bool) or not isinstance(transition_panels, int)
            or transition_panels < 1):
        raise ValueError('positive transition panel count required')
    if (isinstance(interior_panels, bool) or not isinstance(interior_panels, int)
            or interior_panels < 1):
        raise ValueError('positive interior panel count required')
    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError('nonnegative derivative order required')
    with ctx.workprec(bits):
        center, inner, outer = map(arb, _common_window(family))
        w_profile, u_profile = actual_profile_functions(family)
        endpoint_cache = (
            chebyshev_endpoint_derivative_bounds(w_profile, max(order + 1, 9)),
            chebyshev_endpoint_derivative_bounds(u_profile, max(order + 1, 9)),
        )
        sums = {key: arb(0) for key in keys}
        intervals = []
        for left, right, count in ((-outer, -inner, transition_panels),
                                   (-inner, inner, interior_panels),
                                   (inner, outer, transition_panels)):
            for i in range(count):
                intervals.append((center + left + (right - left) * i / count,
                                  center + left + (right - left) * (i + 1) / count, 0))
        cells = 0
        while intervals:
            left, right, depth = intervals.pop()
            try:
                series = monomial_interval_series(
                    family, left.union(right), keys, order, endpoint_cache)
                values = {key: G._coeffs(item, order)[order] * factorial(order)
                          for key, item in series.items()}
            except G.SubdivisionNeeded:
                if depth >= 12:
                    raise
                middle = (left + right) / 2
                intervals.extend(((left, middle, depth + 1), (middle, right, depth + 1)))
                continue
            width = right - left
            for key, value in values.items():
                if not value.is_finite():
                    raise ArithmeticError('nonfinite monomial derivative enclosure')
                sums[key] += width * abs(value).upper()
            cells += 1
        return {key: value.upper() for key, value in sums.items()}, cells


def enclose_profile_fourier_jet(family, keys, period_origin, period_length, *,
                                derivative_order=DEFAULT_DERIVATIVE_ORDER,
                                transition_panels=DEFAULT_TRANSITION_PANELS,
                                interior_panels=DEFAULT_INTERIOR_PANELS,
                                quadrature_points=1024, retained_index=64,
                                max_derivative=2, bits=120):
    """Retained DFT balls with jet L1 alias/tail bounds for w^p U^q.

    Point samples reuse owned order-0 axial series. The L1 bound that feeds
    the classical alias identity is the new midpoint/endpoint jet, not the
    v2 interval-Clenshaw `_derivative_l1`.
    """
    from .nsc_ks_ball_operator import AnalyticRadiusFamily
    from .nsc_ks_value_evaluator import one_direction_metric

    if (not keys or len(set(keys)) != len(keys)
            or any(len(k) != 2 or min(k) < 0 or sum(k) < 1 for k in keys)):
        raise ValueError('distinct positive-total-degree profile powers required')
    with ctx.workprec(bits):
        origin, length = arb(period_origin), arb(period_length)
        center, inner, outer = map(arb, _common_window(family))
        if not origin < center - outer or not origin + length > center + outer:
            raise ValueError('compact support must lie strictly inside the numerical period')
        l1, cells = monomial_derivative_l1(
            family, keys, order=derivative_order, transition_panels=transition_panels,
            interior_panels=interior_panels, bits=bits)
        errors = {key: alias_and_tail_bounds(
            l1[key], length, derivative_order, quadrature_points, retained_index,
            max_derivative) for key in keys}
        model = AnalyticRadiusFamily(one_direction_metric(family))
        samples = {key: [] for key in keys}
        for node in range(quadrature_points):
            z = origin + length * node / quadrature_points
            if z <= center - outer or z >= center + outer:
                values = {key: arb(0) for key in keys}
            else:
                w, U = model.profile_series(z, 0)
                values = {key: w[0] ** key[0] * U[0] ** key[1] for key in keys}
            for key in keys:
                samples[key].append(acb(values[key]))
        results = {}
        for key in keys:
            transformed = acb.dft(samples[key])
            alias, tails = errors[key]
            error = arb((0, 0), alias.man_exp())
            coefficients = tuple(
                transformed[k % quadrature_points] / quadrature_points + acb(error, error)
                for k in range(-retained_index, retained_index + 1))
            results[key] = ProfileFourierBound(
                key, coefficients, l1[key], alias, tails, origin, length,
                derivative_order, quadrature_points, retained_index, cells)
        return results
