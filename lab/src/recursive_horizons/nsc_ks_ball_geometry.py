"""Ball-arithmetic Taylor jets for the owned KS chart, shift, and axial plateau.

Coefficients enclose derivatives/n!. They are not a field or gate certificate.
"""
from dataclasses import dataclass
from fractions import Fraction

from flint import acb, arb, arb_series, ctx

from .nsc_local_incoming_family import LocalAxialFunction

MAX_JET_ORDER = 16
FLAT_EPS = Fraction(1, 256)
RHO_MIN = Fraction(1, 1)
RHO_MAX = Fraction(33, 32)
AXIAL_LOWER = Fraction(4, 5)
SHIFT_DERIVATIVE_BOUND = Fraction(5, 4)


class SubdivisionNeeded(ValueError):
    """The ball is too broad to classify without crossing a non-analytic wall."""


class _JetWork:
    __slots__ = ('order', '_cap', '_prec')

    def __init__(self, order):
        self.order = _order(order)

    def __enter__(self):
        self._cap, self._prec = ctx.cap, ctx.prec
        ctx.cap = self.order + 1
        ctx.prec = max(self._prec, 80)
        return self.order

    def __exit__(self, *exc):
        ctx.cap, ctx.prec = self._cap, self._prec


def _order(order):
    if isinstance(order, bool) or not isinstance(order, int) or not 0 <= order <= MAX_JET_ORDER:
        raise ValueError('jet order must be an integer in [0, 16]')
    return order


def _binary_arb(value, name):
    if isinstance(value, arb):
        ball = value
    elif isinstance(value, Fraction):
        ball = arb(value.numerator) / arb(value.denominator)
    elif isinstance(value, bool) or isinstance(value, str):
        raise ValueError(f'binary real {name} required; decimal strings are not used')
    elif isinstance(value, int):
        ball = arb(value)
    else:
        try:
            ball = arb(float(value))
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(f'binary real {name} required') from error
    if not ball.is_finite():
        raise ValueError(f'finite {name} required')
    return ball


def _rho_ball(value):
    rho = _binary_arb(value, 'rho')
    if not (rho >= arb(1) and rho <= arb(33) / 32):
        raise ValueError('rho must lie in [1, 33/32]')
    return rho


def _identity(ball, order):
    if order == 0:
        return arb_series([ball], prec=1)
    return arb_series([ball, arb(1)], prec=order + 1)


def _coeffs(series, order):
    values = list(series.coeffs())
    if len(values) < order + 1:
        values.extend(arb(0) for _ in range(order + 1 - len(values)))
    return values[:order + 1]


def _finite_series(series, name, order):
    values = _coeffs(series, order)
    if any(not c.is_finite() for c in values):
        raise ArithmeticError(f'nonfinite {name} jet')
    return arb_series(values, prec=order + 1)


def _zero_const(series, order):
    values = _coeffs(series, order)
    values[0] = arb(0)
    return arb_series(values, prec=order + 1)


def _constant_series(value, order):
    return arb_series([value], prec=order + 1)


def restricted_axial_lower_bound():
    """a>=4/5 on [1, 33/32] from pi>3.14 and atan(1/65)<1/65.

    a^2 is decreasing on the slab because
    d/drho[(1+rho^2)(pi/2-atan rho)-rho]=2(rho*angle-1) and
    rho*angle<=(33/32)(pi/4)<(33/32)(22/28)=363/448<1. The minimum is
    therefore at rho=33/32, where angle=pi/4-atan(1/65).
    """
    angle_lo = Fraction(314, 100) / 4 - Fraction(1, 65)
    rho = RHO_MAX
    a2_lo = 3 * ((1 + rho * rho) * angle_lo - rho) - 1
    if a2_lo < AXIAL_LOWER * AXIAL_LOWER:
        raise ArithmeticError('rational axial lower bound failed')
    return AXIAL_LOWER


def inv_a_integrand(z, analytic):
    """1/a with the python-flint analytic flag forwarded to acb.sqrt."""
    # atan is analytic throughout the open right half-plane. Reject a
    # requested complex analyticity region that crosses its branch domain.
    if analytic and not z.real > 0:
        return acb('nan')
    angle = acb.pi() / 2 - z.atan()
    inside = 3 * ((1 + z * z) * angle - z) - 1
    return 1 / inside.sqrt(analytic=analytic)


@dataclass(frozen=True)
class BackgroundSeries:
    r: object
    a: object
    inv_a: object
    inv_a2: object


def background_series(rho_ball, order):
    """Jets of r=sqrt(1+rho^2), a, 1/a and 1/a^2 on the restricted slab."""
    with _JetWork(order) as n:
        rho = _identity(_rho_ball(rho_ball), n)
        radius = (1 + rho * rho).sqrt()
        angle = arb.pi() / 2 - rho.atan()
        inside = 3 * ((1 + rho * rho) * angle - rho) - 1
        axial = inside.sqrt()
        if any(not c.is_finite() for c in _coeffs(axial, n)) or not (axial[0] > 0):
            raise SubdivisionNeeded(
                'axial scale is not proved positive on this rho ball; subdivide')
        inv_a = axial.inv()
        return BackgroundSeries(
            r=_finite_series(radius, 'r', n),
            a=_finite_series(axial, 'a', n),
            inv_a=_finite_series(inv_a, '1/a', n),
            inv_a2=_finite_series(inv_a * inv_a, '1/a^2', n),
        )


def shift_series(rho_ball, order):
    """s=T(rho)-T(1)=-int_1^rho d rho/a, with interval value from |T'|<=5/4."""
    restricted_axial_lower_bound()
    with _JetWork(order) as n:
        rho = _rho_ball(rho_ball)
        inv_a = background_series(rho, n).inv_a
        higher = _finite_series((-inv_a).integral(), 'shift tail', n)
        center = rho.mid()
        integral = acb.integral(inv_a_integrand, acb(1), acb(center))
        if not integral.is_finite() or not integral.imag.contains(0):
            raise ArithmeticError('acb.integral of 1/a did not return a real enclosure')
        value = -integral.real + (arb(5) / 4) * rho.rad() * arb(0, 1)
        coeffs = _coeffs(higher, n)
        coeffs[0] = value
        return _finite_series(arb_series(coeffs, prec=n + 1), 'shift', n)


def _analytic_eta(u):
    q = (arb(1) - u).inv() - u.inv()
    if u[0] >= arb(1)/2:
        small = (-q).exp()
        return small/(arb(1)+small)
    return (arb(1) + q.exp()).inv()


def _flat_complement_bounds(order):
    eps = arb(1) / 256
    if order and not (eps <= arb(2) / (3 * order)):
        raise ValueError('flat-strip Cauchy bound needs eps<=2/(3n)')
    t = (1 / (1 - 3 * eps / 2) - 2 / (3 * eps)).exp()
    if not (t < 1):
        raise ArithmeticError('flat-strip |f| bound is not strictly less than 1')
    scale = t / (1 - t)
    bound, growth = [scale], arb(1)
    two_over_eps = arb(2) / eps
    for _ in range(order):
        growth *= two_over_eps
        bound.append(growth * scale)
    return bound


def _compose_strip(u, order, left):
    # Cauchy balls bound eta^{(k)}/k! on the strip; compose with nonconstant u.
    mag = _flat_complement_bounds(order)
    const = arb(1).union(arb(1) - mag[0]) if left else arb(0).union(mag[0])
    acc = _constant_series(const, order)
    delta, power = _zero_const(u, order), _constant_series(arb(1), order)
    for k in range(1, order + 1):
        power *= delta
        acc += mag[k] * arb(0, 1) * power
    return _finite_series(acc, 'plateau strip', order)


def plateau_series(distance_series, inner, outer, order):
    """C-infinity eta(u)=expit(1/u-1/(1-u)), u=(distance-inner)/(outer-inner)."""
    if not isinstance(distance_series, arb_series):
        raise TypeError('distance_series must be an arb_series')
    with _JetWork(order) as n:
        inner, outer = _binary_arb(inner, 'inner'), _binary_arb(outer, 'outer')
        width = outer - inner
        if not (inner > 0 and width > 0):
            raise ValueError('0 < inner < outer required')
        distance = _finite_series(distance_series, 'distance', n)
        d0 = distance[0]
        if d0 <= inner:
            return _constant_series(arb(1), n)
        if d0 >= outer:
            return _constant_series(arb(0), n)
        u = _finite_series((distance - inner) / width, 'u', n)
        u0, eps = u[0], arb(1) / 256
        if u0 <= eps:
            return _compose_strip(u, n, left=True)
        if u0 >= arb(1) - eps:
            return _compose_strip(u, n, left=False)
        if u0 > 0 and u0 < 1:
            try:
                return _finite_series(_analytic_eta(u), 'plateau', n)
            except ArithmeticError as error:
                raise SubdivisionNeeded('plateau interval loses a denominator margin; subdivide') from error
        raise SubdivisionNeeded(
            'distance ball is too broad to classify against the flat endpoints')


def _chebyshev_series(coefficients, x, order):
    terms = [_binary_arb(c, 'Chebyshev coefficient') for c in coefficients]
    prec = order + 1
    if len(terms) == 1:
        return arb_series([terms[0]], prec=prec)
    b1 = b2 = arb_series([0], prec=prec)
    two_x = 2 * x
    for ak in reversed(terms[1:]):
        b1, b2 = arb_series([ak], prec=prec) + two_x * b1 - b2, b1
    return arb_series([terms[0]], prec=prec) + x * b1 - b2


def local_axial_series(profile, z_ball, order):
    """Jet of an owned LocalAxialFunction at a real z ball."""
    if not isinstance(profile, LocalAxialFunction):
        raise TypeError('owned LocalAxialFunction required')
    with _JetWork(order) as n:
        z0 = _binary_arb(z_ball, 'z')
        center = _binary_arb(profile.center, 'center')
        z = _identity(z0, n)
        if abs(z0-center).upper() <= _binary_arb(profile.inner, 'inner'):
            window = _constant_series(arb(1), n)
        elif z0 >= center:
            distance = z - center
            window = plateau_series(distance, profile.inner, profile.outer, n)
        elif z0 <= center:
            distance = center - z
            window = plateau_series(distance, profile.inner, profile.outer, n)
        else:
            raise SubdivisionNeeded(
                'z ball straddles the plateau center; |z-center| is not analytic')
        off, scale = profile._derivatives[0].mapparms()
        x = _binary_arb(float(off), 'mapparms off') + _binary_arb(
            float(scale), 'mapparms scale') * z
        poly = _finite_series(_chebyshev_series(profile.coefficients, x, n),
                              'Chebyshev', n)
        return _finite_series(window * poly, 'LocalAxialFunction', n)
