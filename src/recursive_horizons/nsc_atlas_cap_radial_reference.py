"""Smooth-potential reference arithmetic and one moving affine radial map.

These diagnostic components do not change the cap producer or evolve a wall.
A panel caller supplies one smooth owner/formula piece and its fixed-r jets.
Solver tolerances are requested local tolerances, not an error enclosure.
"""
import numpy as np
from scipy.integrate import solve_ivp

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_characteristic_background import _array
from .nsc_atlas_polar_cap import _scalar


def reference_potential(f, alpha, *, budget):
    """Analytic potential and its positive gap from the native small-q series.

    The gap is evaluated as a sinh-series remainder, not by subtracting
    rounded potentials. Terms through q**10 retain its relative accuracy
    on |q|<1e-4. Returned series bounds cover the omitted real-arithmetic
    terms only; floating-point/elementary-function error is separate.
    Each returned point charges 12 coefficient evaluations to the caller.
    """
    budget = _budget(budget)
    f, alpha = _array(f, "field"), _scalar(alpha, "alpha", nonnegative=True)
    budget.check(12*f.size)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            scale = np.sqrt(alpha/6)
            if alpha > 0 and scale == 0:
                raise ValueError("reference alpha scale underflows")
            q = scale*f
            ch, sh = np.cosh(q), np.sinh(q)
            sinc = np.ones_like(q)
            np.divide(sh, q, out=sinc, where=q != 0)
            factor = (ch*ch+1)/2
            P, U = f*f*sinc*sinc*factor, 2*f*sinc*ch**3
            # Use safe small arguments for the remainder on all array slots.
            small = abs(q) < 1e-4
            z = np.where(small, q, 0.)
            z2 = z*z
            series = 1+z2/6+z2*z2/120
            tail = z**6*(1/5040+z2*(1/362880+z2/39916800))
            qtail_prime = z**6*(6/5040+z2*(8/362880+10*z2/39916800))
            qseries_prime = z2/3+z2*z2/30
            small_factor = (np.cosh(z)**2+1)/2
            qfactor_prime = z*np.cosh(z)*np.sinh(z)
            d = small_factor*tail*(2*series+tail)
            qd_prime = (qfactor_prime*tail*(2*series+tail)
                +small_factor*(2*(series+tail)*qtail_prime+2*tail*qseries_prime))
            gap, gap_f = f*f*d, f*(2*d+qd_prime)
            # First omitted term z**12/13!, followed by positive terms with
            # ratios <=z**2/210; for q*d(tail)/dq the ratio is <=z**2/180.
            rest = z**12/6227020800/(1-z2/210)
            qrest_prime = 12*z**12/6227020800/(1-z2/180)
            h = rest*(2*(series+tail)+rest)
            d_bound = small_factor*h
            qd_bound = abs(qfactor_prime)*h+small_factor*(
                2*rest*qtail_prime+2*(series+tail+rest)*qrest_prime+2*rest*qseries_prime)
            p_bound, u_bound = f*f*d_bound, abs(f)*(2*d_bound+qd_bound)
    except FloatingPointError as error:
        raise ValueError("nonfinite reference potential") from error
    result = dict(Pdiff=P, Pdiff_f=U, implementation_gap=gap,
                  implementation_gap_f=gap_f, gap_series_bound=p_bound,
                  gap_f_series_bound=u_bound)
    if not all(np.all(np.isfinite(v)) for v in result.values()):
        raise ValueError("nonfinite reference potential")
    return result


def affine_radial_map(coefficients, *, left, right, left_rate, right_rate,
                      budget, rtol=2e-11, atol=2e-13):
    """Integrate (S,T,S_x,T_x) on one declared smooth radial piece.

    coefficients(r) returns scalar a,b,a_x,b_x,a_r,b_r for the reference
    equations; it owns any field-sampling charges. Here x derivatives are
    fixed-radius partials. Endpoint rates map them to fixed panel coordinates.
    Calls at the closed endpoints need explicit owner/formula limits. A
    caller must split discontinuities first; this function discovers none.
    Each RHS charges 16 coefficient/state components, in addition to callback
    work. A zero-width backward-born panel is handled analytically.
    """
    budget = _budget(budget)
    left, right = (_scalar(v, name, nonnegative=True) for v,name in
                   ((left,"left"),(right,"right")))
    left_rate, right_rate = (_scalar(v, name) for v,name in
                             ((left_rate,"left rate"),(right_rate,"right rate")))
    rtol, atol = (_scalar(v, name, positive=True) for v,name in
                 ((rtol,"rtol"),(atol,"atol")))
    if right < left or right >= 1:
        raise ValueError("reference panel requires 0 <= left <= right < 1")
    if rtol < 1e-13 or rtol > 1e-5 or atol < 1e-15 or atol > 1e-7:
        raise ValueError("reference solver tolerances outside the declared range")
    if not callable(coefficients):
        raise ValueError("reference coefficient callback required")
    width, width_rate = right-left, right_rate-left_rate
    if width == 0 and width_rate > 0:
        raise ValueError("zero panel has negative width on the backward ray")
    before, calls = budget.evaluations, 0

    def fields(r):
        nonlocal calls
        budget.check(16)
        calls += 1
        row = coefficients(float(r))
        if not isinstance(row, dict) or set(row) != {"a","b","a_x","b_x","a_r","b_r"}:
            raise ValueError("reference panel needs exactly six fixed-radius coefficients")
        return {key:_scalar(value,key) for key,value in row.items()}

    if width == 0:
        row = fields(left)
        values = np.array([1.,0.,-width_rate*row["a"],width_rate*row["b"]])
    else:
        def rhs(u, state):
            # Keep exact endpoints instead of a last-bit excursion caused
            # by l+(h-l). The intermediate panel coordinate stays in [0,1].
            r = left if u == 0 else right if u == 1 else left+u*width
            row = fields(r)
            velocity = left_rate+u*width_rate
            a, b = width*row["a"], width*row["b"]
            ax = width_rate*row["a"]+width*(row["a_x"]+velocity*row["a_r"])
            bx = width_rate*row["b"]+width*(row["b_x"]+velocity*row["b_r"])
            S,T,Sx,Tx = state
            value = np.array([-a*S,b-a*T,-a*Sx-ax*S,bx-ax*T-a*Tx])
            if not np.all(np.isfinite(value)):
                raise ValueError("nonfinite reference panel variation")
            return value
        solution = solve_ivp(rhs,(0.,1.),[1.,0.,0.,0.],method="DOP853",
                             rtol=rtol,atol=atol,max_step=.25)
        if not solution.success or solution.t[-1] != 1.:
            raise RuntimeError("reference panel integration failed: "+solution.message)
        values = solution.y[:,-1]
    if not np.all(np.isfinite(values)) or values[0] <= 0:
        raise ValueError("nonfinite or nonpositive reference affine map")
    budget.check()
    return dict(zip(("S","T","S_x","T_x"),map(float,values)),
        coefficient_calls=calls,evaluations=budget.evaluations-before,
        born=width == 0,requested_rtol=rtol,requested_atol=atol,
        uniform_error_enclosure=False)
