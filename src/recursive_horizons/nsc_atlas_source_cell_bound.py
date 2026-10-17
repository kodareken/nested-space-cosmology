"""Whole-cell source representation error, with no sampling or ODE solve.

The cubic's binary coefficients and the source parameters define exact
mathematical inputs. Arb encloses those objects, not floating evaluation
of SourceSnapshot or a coupled spacetime evolution.
"""
from fractions import Fraction as Q
import math

from flint import arb, ctx, fmpq

from .nsc_atlas_background_cap import _budget


def _exact(value):
    if isinstance(value, bool):
        raise ValueError('finite real input required')
    value = Q(value)
    if not math.isfinite(float(value)):
        raise ValueError('finite real input required')
    return value


def _ball(value):
    value = _exact(value)
    return arb(fmpq(value.numerator, value.denominator))


def _upper(value):
    return abs(value).upper()


def _record(value):
    value = _upper(value)
    if not value.is_finite():
        raise ValueError('nonfinite source bound')
    mantissa, exponent = value.man_exp()
    display = 0. if value == 0 else math.nextafter(float(value), math.inf)
    if not math.isfinite(display):
        raise ValueError('unrepresentable source bound')
    return dict(upper_float=display, upper_dyadic=dict(mantissa=str(mantissa), exponent=int(exponent)))


def _coordinate(x, geometry):
    scale = _ball(.8)
    u, va, vb = (_ball(geometry[k]) for k in ('U0', 'Va', 'Vb'))
    e = (x/scale).exp()
    z = 4*((1-x/scale)*e/u-(va+vb)/2)/(vb-va)
    zx = -4*x*e/(scale**2*u*(vb-va))
    zxx = -4*e*(1+x/scale)/(scale**2*u*(vb-va))
    return z, zx, zxx


def _bump_jets(z, zx, zxx):
    a = 1-z*z
    bz = -8*z*a**3
    return a**4, bz*zx, 8*a*a*(7*z*z-1)*zx*zx+bz*zxx


def _pulse_jets(left, right, geometry):
    """Monotone coordinate endpoint hull, intersected with compact support."""
    lo, hi = _coordinate(left, geometry), _coordinate(right, geometry)
    if hi[0].upper() <= -1 or lo[0].lower() >= 1:
        return (arb(0),)*3
    zl, zh = max(arb(-1), lo[0].lower()), min(arb(1), hi[0].upper())
    z = zl.union(zh)
    zx = lo[1].lower().union(hi[1].upper())
    zxx = lo[2].lower().union(hi[2].upper())
    values = _bump_jets(z, zx, zxx)
    if lo[0].lower() <= -1 or hi[0].upper() >= 1:
        values = tuple(v.union(arb(0)) for v in values)
    return values


def _polynomial_jets(power, t, width):
    a, b, c, d = power
    return a+t*(b+t*(c+t*d)), (b+t*(2*c+3*d*t))/width, (2*c+6*d*t)/width**2


def _Q(x, h, w, g, gx, gxx, alpha):
    s, v = x+alpha*h, 1+alpha*w
    return (g*gxx-2*gx*gx)/(6*s)+g*gx*v/(3*s*s)


def radial_error_bounds(length, residual, Lh, Lw):
    """Integral-remainder estimate for a shared zero left value/first jet."""
    if not length > 0 or any(v.lower() < 0 for v in (residual, Lh, Lw)):
        raise ValueError('positive length and nonnegative error inputs required')
    theta = (Lh*length**2/2+Lw*length).upper()
    if not theta < 1:
        raise ValueError('source residual inverse is not certified: theta >= 1')
    forcing = residual/(1-theta)
    return theta, (length**2*forcing/2).upper(), (length*forcing).upper()


def frame_jacobian(s, v, g, gx, alpha):
    """Exact derivatives in the independent variables (s,v,g,gx)."""
    def positive_power(value, n):
        # Repeated ball multiplication can include zero for a wide positive
        # input. Monotonic endpoint powers retain its proved positive floor.
        if not value > 0:
            raise ValueError('frame derivative denominator is not positive')
        return (value.lower()**n).union(value.upper()**n)
    s2, s3, s4 = (positive_power(s, n) for n in (2, 3, 4))
    F = 1-alpha*g*g/(6*s2)
    r = (s2-alpha*g*g/6).sqrt()
    r3, F2 = positive_power(r, 3), positive_power(F, 2)
    n = s*v-alpha*g*gx/6
    z = g*gx/s2-g*g*v/s3
    f = (alpha*g*g/(3*s3), arb(0), -alpha*g/(3*s2), arb(0))
    dz = (-2*g*gx/s3+3*g*g*v/s4, -g*g/s3,
          gx/s2-2*g*v/s3, g/s2)
    return dict(r=(s/r, arb(0), -alpha*g/(6*r), arb(0)),
        log_C=tuple(a/F for a in f), c=(arb(0), arb(0), arb(1), arb(0)),
        r_x=(v/r-n*s/r3, s/r, -alpha*gx/(6*r)+n*alpha*g/(6*r3), -alpha*g/(6*r)),
        log_C_x=tuple(-alpha*(a/F-z*b/F2)/3 for a, b in zip(dz, f)),
        c_x=(arb(0), arb(0), arb(0), arb(1)))


def source_cell_error_bound(geometry, *, alpha, left, right, power, budget):
    """Enclose the first active source cell using its actual cubic powers.

    power is in increasing powers of t=(x-left)/(right-left). The caller
    extracts it exactly from the stored binary PPoly coefficients. The
    exact target uses the analytic compact source; the snapshot's additional
    mask at the declared front is bounded separately. No derivatives are
    replaced by their ODE targets and no adaptive subdivision is performed.
    """
    budget = _budget(budget); before = budget.evaluations
    left, right, alpha = map(_exact, (left, right, alpha))
    power = tuple(map(_exact, power))
    front = _exact(geometry['support_x'][0])
    if (len(power) != 4 or power[0] or power[1] or alpha < 0
            or not 0 < left < front < right < _exact(.8)
            or not _exact(geometry['U0']) < 0
            or not _exact(geometry['Va']) < _exact(geometry['Vb']) < 0):
        raise ValueError('require a zero-left-jet first source cell and original IV signs')
    budget.check(2048)
    with ctx.workprec(160):
        a, b, al, cut = map(_ball, (left, right, alpha, front))
        width = b-a; coeff = tuple(map(_ball, power))
        if not _coordinate(a, geometry)[0] < -1 or not _coordinate(b, geometry)[0] < 0:
            raise ValueError('cell must start before analytic support and remain in its first half')
        panels = []
        for lo, hi in ((a, cut), (cut, b)):
            budget.check(2048)
            t = ((lo-a)/width).union((hi-a)/width)
            h, w, hxx = _polynomial_jets(coeff, t, width)
            g, gx, gxx = _pulse_jets(lo, hi, geometry)
            x = lo.union(hi)
            panels.append(dict(lo=lo, hi=hi, x=x, h=h, w=w, hxx=hxx, g=g, gx=gx, gxx=gxx))
        budget.check(4096)
        G0, G1, G2 = (max(_upper(p[k]) for p in panels) for k in ('g', 'gx', 'gxx'))
        smin, vmax, Fmin = _ball(Q(1, 5)), _ball(Q(3, 2)), _ball(Q(1, 2))
        qmax = ((G0*G2+2*G1*G1)/(6*smin)+G0*G1*vmax/(3*smin*smin)).upper()
        exact_slo = a-al*width**2*qmax/2
        exact_shi = b+al*width**2*qmax/2
        exact_vabs = 1+al*width*qmax
        if not exact_slo > smin or not exact_vabs < vmax or not 1-al*G0*G0/(6*smin*smin) > Fmin:
            raise ValueError('exact source continuation tube is not certified')
        approximate_s = [(p['x']+al*p['h']) for p in panels]
        if any(not s > smin for s in approximate_s) or any(not abs(1+al*p['w']) < vmax for p in panels):
            raise ValueError('source polynomial exits the common comparison tube')
        Lh = (al*(G0*G2+2*G1*G1)/(6*smin*smin)+2*al*G0*G1*vmax/(3*smin**3)).upper()
        Lw = (al*G0*G1/(3*smin*smin)).upper()
        defects = [_upper(p['hxx']-_Q(p['x'], p['h'], p['w'], p['g'], p['gx'], p['gxx'], al)) for p in panels]
        D = max(defects)
        theta, H, W = radial_error_bounds(width, D, Lh, Lw)
        # The declared front can lie slightly after the exact analytic root.
        # Hull that omitted sliver; never declare its analytic field zero.
        zfront = _coordinate(cut, geometry)[0]
        mask0 = mask1 = arb(0)
        if zfront.upper() > -1:
            zsliver = arb(-1).union(zfront.upper())
            _, zx, zxx = _coordinate(b, geometry)
            mask = _bump_jets(zsliver, arb(0).union(zx.upper()), arb(0).union(zxx.upper()))
            mask0, mask1 = _upper(mask[0]), _upper(mask[1])
        smax = max(exact_shi.upper(), *(s.upper() for s in approximate_s))
        derivatives = frame_jacobian(smin.union(smax), (-vmax).union(vmax),
                                     (-G0).union(G0), (-G1).union(G1), al)
        differences = (al*H, al*W, mask0, mask1)
        frames = {name: _record(sum(_upper(v)*e for v, e in zip(row, differences)))
                  for name, row in derivatives.items()}
        result = dict(precision_bits=ctx.prec, panel_count=2,
            radial=dict(residual=_record(D), Lh=_record(Lh), Lw=_record(Lw), theta=_record(theta),
                        h_error=_record(H), h_x_error=_record(W)),
            analytic_source_bounds={k: _record(v) for k, v in zip(('g', 'g_x', 'g_xx'), (G0, G1, G2))},
            front_mask_difference=dict(c=_record(mask0), c_x=_record(mask1)),
            frame_error=frames,
            panels=[dict(left=str(left), right=str(front), residual=_record(defects[0])),
                    dict(left=str(front), right=str(right), residual=_record(defects[1]))],
            tube=dict(rJ_lower='1/5', absolute_rJ_x_upper='3/2', F_lower='1/2',
                      exact_continuation_proved=True, polynomial_inside=True),
            charged_values=budget.evaluations-before,
            scope='exact binary source cubic and analytic pulse; floating evaluation and coupled exterior errors excluded')
    budget.check()
    return result
