"""Certified exact-reference gauge bounds, separate from floating root errors.

Two fixed endpoint integrals enclose the analytic inverse. No root search,
source solve or trajectory is performed. Constants are the rational reference
a=4/5, R0=9/20, x0=2/5; the numerical gauge's input/roundoff errors are external.
"""
from fractions import Fraction as Q
import math

from flint import acb, arb, ctx, fmpq

from .nsc_atlas_background_cap import _budget


def _ball(value):
    if isinstance(value, (bool, float)):
        raise ValueError('supply an explicit exact rational input')
    value = Q(value)
    return arb(fmpq(value.numerator, value.denominator))


def _dyadic(value):
    mantissa, exponent = value.man_exp()
    return dict(mantissa=str(mantissa), exponent=int(exponent))


def _record(value):
    if not value.is_finite():
        raise ValueError('nonfinite reference enclosure')
    lo, hi = value.lower(), value.upper()
    return dict(lower=_dyadic(lo), upper=_dyadic(hi),
                lower_float=math.nextafter(float(lo), -math.inf),
                upper_float=math.nextafter(float(hi), math.inf))


def _state(R):
    bm, bp = 1/(5*R**2)+3*R/4, 1/(5*R**2)-5*R/4
    q = 25*R**2/16-1+3/(10*R)+1/(25*R**4)
    return bm, bp, q


def _advanced(R, analytic):
    # Forward the branch flag even though the real integration path is regular.
    _, bp, q = _state(R)
    return -5*R/(4-5*R)*(1+bp/q.sqrt(analytic=analytic))


def _integral(endpoint, budget):
    def integrand(R, analytic):
        budget.check(64)
        return _advanced(R, analytic)
    result = acb.integral(integrand, acb(_ball(Q(9, 20))), acb(_ball(endpoint)),
                          rel_tol=arb(2)**-80, abs_tol=arb(2)**-80,
                          deg_limit=64, eval_limit=512, depth_limit=20)
    budget.check()
    if not result.is_finite() or not result.imag.contains(0):
        raise ValueError('endpoint integral did not enclose a finite real value')
    return result.real


def reference_gauge_bound(left, right, radius_lower, radius_upper, *, budget):
    """Prove an all-x radius bracket and enclose its analytic wall coefficients."""
    budget = _budget(budget)
    before = budget.evaluations
    # Refuse implicit binary/decimal reinterpretation at the public boundary.
    for value in (left, right, radius_lower, radius_upper):
        _ball(value)
    left, right, radius_lower, radius_upper = map(Q, (left, right, radius_lower, radius_upper))
    if not Q(2, 5) < left < right < Q(4, 5):
        raise ValueError('incoming cell must be strictly inside the reference horizon')
    if not Q(9, 20) < radius_lower < radius_upper < Q(11, 20):
        raise ValueError('fixed radius bracket must lie inside the reference wall segment')
    budget.check(2048)
    with ctx.workprec(160):
        a, x0 = _ball(Q(4, 5)), _ball(Q(2, 5))
        xl, xr = _ball(left), _ball(right)
        rl, rh = _ball(radius_lower), _ball(radius_upper)
        path = _ball(Q(9, 20)).union(rh)
        _, _, path_q = _state(path)
        if not path_q > 0 or not a-path > 0:
            raise ValueError('real integral path not certified regular')
        # q-bp^2=(a-R)/R>0 proves I'<0 along the entire real path.
        def F(r):
            return (1-r/a)*(r/a).exp()
        targets = [2*a*(F(x)/F(x0)).log() for x in (xl, xr)]
        integrals = [_integral(endpoint, budget) for endpoint in (radius_lower, radius_upper)]
        margins = [integrals[0]-targets[0], targets[1]-integrals[1]]
        if not all(value > 0 for value in margins):
            raise ValueError('predeclared inverse bracket was not proved; no retry')
        x, R = xl.union(xr), rl.union(rh)
        bm, bp, q = _state(R)
        v = q.sqrt()
        difference = (a-R)/(R*(v+bp))  # v-bp, using the exact junction identity
        Rx = 2*x*v*difference/(a-x)
        sx = 2*x*difference/(a-x)
        C = x/R*((x-R)/a).exp()
        W = 4*a*(1-F(R)/F(x))
        Wx = sx**2/C
        vR = (25*R/8-3/(10*R**2)-4/(25*R**5))/(2*v)
        bpR = -2/(5*R**3)-_ball(Q(5, 4))
        Lx_total = 1/x+1/a-(1/R+1/a)*Rx
        log_sx_x = 1/x+1/(a-x)+(vR-bpR)*Rx/difference
        eta_x = -log_sx_x+Lx_total/2
        A0 = 1-R**2
        fields = dict(R=R,A0=A0,v=v,beta_minus=bm,beta_plus=bp,R_x=Rx,
                      s_x=sx,C_wall=C,W=W,W_x=Wx,eta_x=eta_x,
                      Udot=sx/C,Vdot=1/sx,T_x=bm*sx/A0,Z_x=Rx/A0,
                      T_x_minus_Z_x=sx/(bm+v),
                      T_x_plus_Z_x=sx*(bm+v)/A0)
        if not all(fields[key] > 0 for key in fields if key != 'eta_x'):
            raise ValueError('positive reference wall coefficient not enclosed')
        # On the vacuum wedge, r is between x and R(x), C<=1 and r_x<=1.
        Cmin = xl/rh*((xl-rh)/a).exp()
        rxmin = xl*(a-rh)/(rh*(a-xl))
        wedge = dict(r=xl.union(rh),C=Cmin.union(arb(1)),
                     r_x=rxmin.union(arb(1)),
                     r_U=(Cmin*(a-xr)/(4*xr)).union((a-xl)/(4*xl)),
                     K=arb(0).union(1/(5*xl**3)))
        budget.check()
        return dict(precision_bits=160,normalization='exact rational corner x0=2/5,R0=9/20,a=4/5',
                    inputs=dict(left=str(left),right=str(right),radius_lower=str(radius_lower),
                                radius_upper=str(radius_upper)),
                    integral_path_q=_record(path_q),
                    endpoint_integrals=[_record(value) for value in integrals],
                    endpoint_targets=[_record(value) for value in targets],
                    bracket_margins=[_record(value) for value in margins],
                    wall={key:_record(value) for key,value in fields.items()},
                    vacuum_wedge={key:_record(value) for key,value in wedge.items()},
                    charged_values=budget.evaluations-before,
                    scope='exact reference only; numerical gauge errors and sourced coefficient tube remain external')
