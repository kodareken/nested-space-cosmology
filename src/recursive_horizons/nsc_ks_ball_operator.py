"""Taylor enclosures of the existing KS time coefficients and spatial profiles."""
from math import comb, factorial

from flint import arb, arb_series, ctx

from . import nsc_ks_ball_geometry as G
from .nsc_ks_source_envelope import _as_metric
from .nsc_local_incoming_family import LocalAxialFunction


def slab_ball(lower, upper):
    """Enclose a slab, using an exact dyadic boundary ball when lower=1."""
    lower, upper = arb(lower), arb(upper)
    if not lower >= 1 or not upper <= arb(33)/32 or not upper >= lower:
        raise ValueError('ordered slab in [1,33/32] required')
    result = lower.union(upper)
    if not result >= 1 and lower == 1:
        radius = arb(1)/64
        while 1+radius >= upper:
            radius /= 2
        # Tuple input preserves this power-of-two radius exactly. Passing an
        # arb as rad may inflate it by a mag ulp and cross rho=1 again.
        result = arb((1+radius).man_exp(), radius.man_exp())
    if not result >= 1 or not result <= arb(33)/32:
        raise G.SubdivisionNeeded('outward slab ball crosses the declared chart domain')
    return result


class AnalyticRadiusFamily:
    """Same amplitude/basis data as the prepared owner; no fitted coefficients."""
    def __init__(self, family):
        _, self.metric = _as_metric(family)
        if not self.metric.directions:
            raise ValueError('an explicit compatible radius basis is required')
        windows = {(float(d.inner_radius), float(d.outer_radius)) for d in self.metric.directions}
        if len(windows) != 1:
            raise ValueError('this enclosure owns the common normal-window family')
        self.inner, self.outer = next(iter(windows))
        if not all(isinstance(p, LocalAxialFunction) for d in self.metric.directions for p in (d.w, d.U)):
            raise TypeError('analytic LocalAxialFunction profiles required; sampled callbacks are insufficient')
        self.w_zero = all(a == 0 or all(c == 0 for c in d.w.coefficients)
                          for a, d in zip(self.metric.amplitudes, self.metric.directions))
        self.U_zero = all(a == 0 or all(c == 0 for c in d.U.coefficients)
                          for a, d in zip(self.metric.amplitudes, self.metric.directions))

    def profile_series(self, z, order):
        with G._JetWork(order):
            w, U = arb_series([0], prec=order+1), arb_series([0], prec=order+1)
            for amplitude, direction in zip(self.metric.amplitudes, self.metric.directions):
                if amplitude:
                    if any(direction.w.coefficients):
                        w += arb(float(amplitude))*G.local_axial_series(direction.w, z, order)
                    if any(direction.U.coefficients):
                        U += arb(float(amplitude))*G.local_axial_series(direction.U, z, order)
            return w, U

    def time_series(self, rho, angular, order, reciprocal_order):
        with G._JetWork(order):
            bg = G.background_series(rho, order)
            shift = G.shift_series(rho, order)
            window = G.plateau_series(-shift, self.inner, self.outer, order)
            first, third = shift*window/bg.r, shift**3*window/(6*bg.r)
            inv_ar = bg.inv_a/bg.r
            result = {'inv_a2': bg.inv_a2, 'inv_a': bg.inv_a, 'inv_ar': inv_ar}
            for n in range(1, reciprocal_order+1):
                for q in range(n+1):
                    p = n-q
                    if (p and self.w_zero) or (q and self.U_zero):
                        continue
                    result[(p, q)] = arb(float(angular))*(-1)**n*comb(n, q)*inv_ar*first**p*third**q
            return result


def time_operator_enclosure(model, rho_start, rho_end, angular, *, degree=8, reciprocal_order=4, bits=120):
    """Power polynomials in step fraction plus uniform continuous remainders."""
    if degree+1 > G.MAX_JET_ORDER or degree < 0:
        raise ValueError('unsupported time Taylor degree')
    with ctx.workprec(bits):
        start, end = arb(rho_start), arb(rho_end)
        if not start > end:
            raise ValueError('decreasing rho step required')
        h, center = end-start, (start+end)/2
        interval = slab_ball(end, start)
        at_center = model.time_series(center, angular, degree, reciprocal_order)
        on_interval = model.time_series(interval, angular, degree+1, reciprocal_order)
        polynomials, errors = {}, {}
        for name, series in at_center.items():
            c = G._coeffs(series, degree)
            polynomials[name] = [sum((c[j]*h**j*comb(j,k)*(-arb(1)/2)**(j-k)
                                     for j in range(k,degree+1)), arb(0)) for k in range(degree+1)]
            derivative = G._coeffs(on_interval[name], degree+1)[degree+1]
            errors[name] = (abs(derivative)*(abs(h)/2)**(degree+1)).upper()
        return polynomials, errors


def spatial_potential_jets(model, z, keys, order, *, bits=120, absolute=False):
    """Ordinary z derivatives of w^p U^q at a point or over an interval."""
    with ctx.workprec(bits), G._JetWork(order):
        w, U = model.profile_series(z, order)
        result = {}
        for p, q in keys:
            series = w**p*U**q
            values = [c*factorial(j) for j, c in enumerate(G._coeffs(series,order))]
            result[(p,q)] = [abs(v).upper() for v in values] if absolute else values
        return result
