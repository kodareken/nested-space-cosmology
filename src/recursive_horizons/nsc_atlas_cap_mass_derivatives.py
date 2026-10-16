"""Derivatives of the declared cap mass integrands on positive radii.

This diagnostic companion leaves the retained metric producer unchanged.
It differentiates the actual small-radius and potential series formulas.
There is no field preparation, evolution, extrapolation or finite differencing.
"""
import math

import numpy as np

from .nsc_atlas_characteristic_background import (
    CharacteristicProfile, MAX_POINTS, SMALL_RADIUS, _array, _budget, _moment,
)
from .nsc_atlas_polar_cap import _potential, _scalar


def third_moments(profile, T, distance, *, budget):
    """Integrate z**k F'''(T+distance*z), k=0,1,2, on the actual pieces.

    C2 matching cancels moving-knot boundary terms in first derivatives of
    integrals of F' and F''. A discontinuous F''' is integrated piecewise.
    """
    budget = _budget(budget)
    if not isinstance(profile, CharacteristicProfile):
        raise ValueError("one CharacteristicProfile is required")
    T, distance = _scalar(T, "T"), _scalar(distance, "distance", positive=True)
    profile._domain(np.array([T-distance, T+distance]))
    labels = profile.breakpoints
    internal = labels[(labels > T-distance) & (labels < T+distance)]
    cuts = np.r_[-1., (internal-T)/distance, 1.]
    result = np.zeros(3)
    budget.check(9*(len(cuts)-1))
    for low, high in zip(cuts[:-1], cuts[1:]):
        midpoint = T+distance*(low+high)/2
        index = int(np.searchsorted(labels, midpoint, side="right")-1)
        if index < 0:
            continue
        index = min(index, len(labels)-2)
        coefficient = profile._spline.c[:, index]
        at = T-labels[index]
        degree=len(coefficient)-1
        budget.check(3*max(0,degree-5))
        moments = [_moment(low, high, power) for power in range(degree)]
        for order in range(degree-2):
            factor = (float(np.polyval(np.polyder(coefficient, order+3), at))
                      * distance**order/math.factorial(order))
            result += factor*np.asarray(moments[order:order+3])
    if not np.all(np.isfinite(result)):
        raise ValueError("nonfinite background moment derivative")
    budget.check()
    return result


def _radial_coefficients(r):
    """Actual truncated Phi coefficient and stable derivatives of d/r,d/(r*A).

    The two derivative evaluation series include terms through r**11.
    Their omitted terms are O(r**13) on 0<r<=.01, below binary64 rounding
    of these coefficients. They are not substituted for the Phi polynomial.
    """
    z = r*r
    c = r*(2/3+z*(4/5+z*(6/7+z*(8/9+z*10/11))))
    c_r = 2/3+z*(12/5+z*(30/7+z*(56/9+z*90/11)))
    rho_r = c+12*r**11/13
    gamma_r = r*(8/3+z*(4/5+z*(16/35+z*(20/63+z*(8/33+z*28/143)))))/(1-z)**2
    return c, c_r, rho_r, gamma_r


def background_partials(profile, T, r, *, budget):
    """Return f,Pi,Phi and their T/r partials at up to 513 positive radii.

    Algebraic endpoint partials are used only along admissible directions by
    the caller. The profile still rejects every query outside its domain.
    No assertion about Pi_r at a merely C2 axis knot is made here.
    """
    budget = _budget(budget)
    if not isinstance(profile, CharacteristicProfile):
        raise ValueError("one CharacteristicProfile is required")
    times, radii = np.broadcast_arrays(_array(T, "T"), _array(r, "r"))
    if times.size > MAX_POINTS or np.any(radii <= 0) or np.any(radii >= 1):
        raise ValueError("mass derivatives require at most 513 points with 0<r<1")
    shape = times.shape
    times, radii = times.ravel(), radii.ravel()
    base = profile.background(times, radii, budget=budget)
    out = {name: np.empty_like(radii) for name in
           ("f_T", "f_r", "Pi_T", "Pi_r", "Phi_T", "Phi_r")}
    ordinary = radii > SMALL_RADIUS
    if np.any(ordinary):
        rr, tt = radii[ordinary], times[ordinary]
        A, d = 1-rr*rr, np.arctanh(rr)
        jets = profile._evaluate(np.stack((tt-d, tt+d)), budget)
        pm, pp = jets["F_prime"]
        sm, sp = jets["F_second"]
        Pi, Phi, q = (base[name][ordinary] for name in ("Pi", "Phi", "q"))
        out["f_T"][ordinary] = (pp-pm)/rr
        out["f_r"][ordinary] = Phi
        out["Pi_T"][ordinary] = (sp-sm)/(rr*A)
        out["Pi_r"][ordinary] = (sp+sm)/(rr*A*A)-Pi*(1-3*rr*rr)/(rr*A)
        out["Phi_T"][ordinary] = (sp+sm)/(A*rr)-(pp-pm)/(rr*rr)
        out["Phi_r"][ordinary] = (sp-sm)/(rr*A*A)+2*q/A-2*Phi/rr
    for i in np.flatnonzero(~ordinary):
        rr, tt = float(radii[i]), float(times[i])
        A, d = 1-rr*rr, float(np.arctanh(rr))
        first, second, weighted = profile._split_integrals(tt, d, budget)
        j0, j1, j2 = third_moments(profile, tt, d, budget=budget)
        rho, gamma = d/rr, d/(rr*A)
        c, c_r, rho_r, gamma_r = _radial_coefficients(rr)
        out["f_T"][i] = rho*second
        out["f_r"][i] = rho_r*first+gamma*weighted
        out["Pi_T"][i] = gamma*j0
        out["Pi_r"][i] = gamma_r*second+gamma*j1/A
        out["Phi_T"][i] = c*second+gamma*j1
        out["Phi_r"][i] = c_r*first+c*weighted/A+gamma_r*weighted+gamma*j2/A
    budget.check(9*len(radii))
    if any(not np.all(np.isfinite(values)) for values in out.values()):
        raise ValueError("nonfinite background partial derivative")
    return {**{name: base[name].reshape(shape) for name in ("f", "Pi", "Phi")},
            **{name: value.reshape(shape) for name, value in out.items()}}


def potential_partial(f, alpha, *, budget):
    """Pdiff and d(Pdiff)/df, including its actual small-q polynomial."""
    budget = _budget(budget)
    f = _array(f, "f")
    alpha = _scalar(alpha, "alpha", nonnegative=True)
    _, difference, derivative = _potential(f, alpha)
    q = np.sqrt(alpha/6)*f
    small = abs(q) < 1e-4
    s = 1+q*q/6+q**4/120
    # s+q*s_q is this polynomial, not cosh(q).
    value = f*s*((np.cosh(q)**2+1)*(1+q*q/2+q**4/24)
                 +q*s*np.cosh(q)*np.sinh(q))
    derivative = np.where(small, value, derivative)
    budget.check(2*f.size)
    if not np.all(np.isfinite(difference)) or not np.all(np.isfinite(derivative)):
        raise ValueError("nonfinite potential partial derivative")
    return {"Pdiff": difference, "Pdiff_f": derivative}
