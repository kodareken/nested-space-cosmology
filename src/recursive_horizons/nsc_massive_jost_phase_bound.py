"""Directed subgap Jost phase error on the entire exterior r >= R.

The original binary ratio polynomial is an approximation, not an exact
asymptotic solution. Its complete phase defect is enclosed, including rounded
coefficients. See docs/nsc-massive-jost-phase-bound.md for the comparison proof.
This does not enclose subsequent ODE integration, horizon sewing or rho=1 fields.
"""
from contextlib import contextmanager
import math

import numpy as np
from flint import arb, arb_series, ctx

from .nsc_ks_ball_trajectory import exact_upper


@contextmanager
def _series_context(bits, degree):
    previous = ctx.cap
    try:
        with ctx.workprec(bits):
            ctx.cap = degree + 1
            yield
    finally:
        ctx.cap = previous


def _variable(center, degree):
    return arb_series([center, 1], prec=degree + 1)


def _polynomial(coefficients, x, degree):
    value = arb_series([0], prec=degree + 1)
    for coefficient in reversed(coefficients):
        value = value*x + coefficient
    return value


def _metric(x, upper_x, terms, degree, *, enclose_tail):
    # A(r)=1-3*((1+r^2)*atan(1/r)-r), x=1/r. The omitted
    # powers start at 2*terms+1; their derivatives are dominated by
    # the positive coefficient majorant 6*x^(2*terms+1)/(1-x).
    value = arb_series([1], prec=degree + 1)
    for j in range(terms):
        value -= arb(6)*(-1)**j*x**(2*j+1)/((2*j+1)*(2*j+3))
    if enclose_tail:
        endpoint = _variable(upper_x, degree)
        majorant = 6*endpoint**(2*terms+1)/(1-endpoint)
        value = arb_series([
            value[k] + arb(0, majorant[k].abs_upper())
            for k in range(degree + 1)
        ], prec=degree + 1)
    return value


def _phase_defect(x, real, imag, energy, mass, angular, upper_x, terms, degree,
                  *, enclose_tail):
    a = _metric(x, upper_x, terms, degree, enclose_tail=enclose_tail)
    p = _polynomial(real, x, degree)
    q = _polynomial(imag, x, degree)
    dp = _polynomial([k*real[k] for k in range(1, len(real))], x, degree)
    dq = _polynomial([k*imag[k] for k in range(1, len(imag))], x, degree)
    norm2 = p*p + q*q
    norm = norm2.sqrt()
    phi_r = -x*x*(p*dq-q*dp)/norm2
    rhs = 2/a*(a.sqrt()*(angular*x/(1+x*x).sqrt()*p+mass*q)/norm-energy)
    return phi_r-rhs


def _endpoint(value, lower=False):
    point = value.lower() if lower else value.upper()
    if not point.is_finite():
        raise ArithmeticError('finite directed endpoint required')
    mantissa, exponent = point.man_exp()
    return {'mantissa': str(mantissa), 'exponent': int(exponent)}


def subgap_phase_bound(energy, mass, angular, radius, coefficients,
                       initializer_phase, *, bits=192, degree=12,
                       metric_terms=12, tube='0.00001'):
    """Bound |theta_decaying(R)-initializer_phase| for the fixed NSC metric.

    Scalars may be exact binary floats or Arb parameter intervals. Coefficients
    are the supplied numerical ratio polynomial, stored as exact binary64 values.
    The proof requires a positive real part and a positive phase derivative on
    the whole exterior tube. Failure of those checks is not a physical obstruction.
    The E=ell=0 case is admitted for an exact manufactured control; production
    source rows still use the original positive-energy preparation contract.
    """
    if (isinstance(degree, bool) or not isinstance(degree, int) or degree < 2
            or isinstance(metric_terms, bool) or not isinstance(metric_terms, int)
            or metric_terms < 1 or 2*metric_terms+1 < degree
            or isinstance(bits, bool) or not isinstance(bits, int) or bits < 80):
        raise ValueError('valid Taylor degree, metric tail order and precision required')
    coefficients = np.asarray(coefficients, dtype=complex)
    if (coefficients.ndim != 1 or not len(coefficients)
            or not np.isfinite(coefficients).all()
            or not math.isfinite(float(initializer_phase))):
        raise ValueError('finite polynomial and numerical initializer required')
    with _series_context(bits, degree):
        energy, mass, angular, radius = map(arb, (energy, mass, angular, radius))
        eta = arb(tube)
        if (not all(v.is_finite() for v in (energy, mass, angular, radius, eta))
                or not energy >= 0 or not mass > energy or not radius >= 8
                or not eta > 0 or not eta < arb(1)/4):
            raise ValueError('subgap positive mass, R>=8 and small positive phase tube required')
        real = [arb(float(z.real)) for z in coefficients]
        imag = [arb(float(z.imag)) for z in coefficients]
        upper_x = (1/radius).upper()
        domain = arb(0).union(upper_x)
        x = _variable(domain, degree)
        p = _polynomial(real, x, degree)[0]
        q = _polynomial(imag, x, degree)[0]
        a = _metric(x, upper_x, metric_terms, degree, enclose_tail=True)[0]
        norm = (p*p+q*q).sqrt()
        if not a > 0 or not p > 0 or not norm > 0:
            raise ArithmeticError('whole-exterior metric/phase chart is not enclosed')
        # |cos(phi+e)-cos(phi)| and |sin(phi+e)-sin(phi)| <= |e|.
        bracket = mass.lower()*((p/norm).lower()-eta.upper())
        bracket -= angular.abs_upper()*upper_x*((q/norm).abs_upper()+eta.upper())
        lam = (2/a.upper().sqrt()*bracket).lower()
        if not lam > 0:
            raise ArithmeticError('positive exterior phase derivative not established')
        # Derivatives at zero below degree are exact for the truncated metric:
        # the omitted metric series starts at power >= degree.
        center = _phase_defect(_variable(arb(0), degree), real, imag, energy,
                               mass, angular, upper_x, metric_terms, degree,
                               enclose_tail=False)
        interval = _phase_defect(x, real, imag, energy, mass, angular,
                                 upper_x, metric_terms, degree, enclose_tail=True)
        residual = sum((center[k].abs_upper()*upper_x**k
                        for k in range(degree)), arb(0))
        residual += interval[degree].abs_upper()*upper_x**degree
        residual = residual.upper()
        error = (residual/lam).upper()
        limit_gap = abs((imag[0]/real[0]).atan()-(energy/mass).asin())
        if (not residual.is_finite() or not error < eta.lower()
                or not limit_gap < eta.lower()):
            raise ArithmeticError('phase tube or outgoing asymptotic branch is unresolved')
        at_r = _variable(1/radius, degree)
        pr = _polynomial(real, at_r, degree)[0]
        qr = _polynomial(imag, at_r, degree)[0]
        rounding = abs((qr/pr).atan()-arb(float(initializer_phase))).upper()
        return {
            'scope': 'decaying subgap phase at R; entire exterior half-line controlled',
            'metric': 'A=1-3*((1+r^2)*atan(1/r)-r)',
            'parameter_intervals': {
                name: {'lower': _endpoint(v, True), 'upper': _endpoint(v)}
                for name, v in [('energy', energy), ('mass', mass), ('angular', angular),
                                ('radius', radius)]},
            'bits': bits, 'taylor_degree': degree, 'metric_terms': metric_terms,
            'phase_tube_upper': exact_upper(eta),
            'metric_lower': _endpoint(a, True),
            'phase_derivative_lower': _endpoint(lam, True),
            'phase_defect_upper': exact_upper(residual),
            'asymptotic_branch_gap_upper': exact_upper(limit_gap),
            'asymptotic_phase_error_upper': exact_upper(error),
            'initializer_rounding_upper': exact_upper(rounding),
            'combined_phase_initial_error_upper': exact_upper(error+rounding),
            'all_exterior_radii': True,
            'radial_transport_error': None,
            'horizon_sewing_error': None,
            'physical_rho1_source_error': None,
        }
