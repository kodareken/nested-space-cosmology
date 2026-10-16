"""Fixed-cap mass gradients for regular, primitive-conforming P1 fields.

Forward radial sensitivities and independently weighted gradient quadrature
check one constrained mass functional. The weak RHS uses its common gradient;
neither its projected wall momentum nor its normal derivative is silently
identified with a reconstructed trace. No time evolution is performed here.
"""
from fractions import Fraction

import numpy as np
from flint import arb, ctx, fmpq
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_geometric_reference import _positive_gap
from .nsc_atlas_cap_radial_reference import reference_potential
from .nsc_atlas_polar_cap import _scalar, _vector

GRADIENT_ABSOLUTE_TOLERANCE = 2e-8


def _exact_ball(value):
    ratio = Fraction(float(value))
    return arb(fmpq(ratio.numerator, ratio.denominator))


def _metric_floor_from_slopes(Z, lc, lk, al, *, scope):
    """Apply the radial inequality to certified global derivative bounds.

    Arb inputs and the working precision belong to the field-space caller.
    """
    R, A0 = Z.tanh(), 1/Z.cosh()**2
    if not (R > 0 and R < 1):
        raise ValueError("radial interval is not resolved inside 0<R<1")
    F, pmax = lc*Z/R, lk*Z/A0
    hmax = lc/A0+F
    argument = (al/6).sqrt()*F
    if not argument <= 100:
        raise ValueError("potential argument is outside the finite reference domain")
    if al == 0 or lc == 0:
        potential = F*F
    else:
        sinc = argument.sinh()/argument
        potential = F*F*sinc*sinc*(argument.cosh()**2+1)/2
    bmax = R*R*potential/2+(pmax*pmax+hmax*hmax)/4
    bound = A0-2*al*bmax
    if not bound.is_finite() or not bound > 0:
        raise ValueError("sufficient all-radius metric bound is not positive")
    lower = bound.lower()
    mantissa, exponent = lower.man_exp()
    down = float(np.nextafter(float(lower), -np.inf))
    if not np.isfinite(down) or down <= 0:
        raise ValueError("positive metric margin is not representable")
    return dict(A_lower=dict(mantissa=str(mantissa), exponent=int(exponent)),
                A_lower_float_down=down, precision_bits=ctx.prec,
                scope=scope)


def _metric_domain(z, chi, k, alpha, budget):
    """Sufficient all-radius bound for the exact declared P1 solution."""
    budget.check(64*len(z)+64)
    with ctx.workprec(160):
        zb, cb, kb = [[_exact_ball(v) for v in values] for values in (z, chi, k)]
        Z, al = zb[-1], _exact_ball(alpha)
        lc = max(abs((b-a)/(v-u)).upper() for a,b,u,v in zip(cb[:-1],cb[1:],zb[:-1],zb[1:]))
        lk = max(abs((b-a)/(v-u)).upper() for a,b,u,v in zip(kb[:-1],kb[1:],zb[:-1],zb[1:]))
        return _metric_floor_from_slopes(
            Z,lc,lk,al,scope="exact radial solution for the declared P1 fields; numerical ODE error excluded")


def _radial_variations(breaks, coefficients, *, centre_lapse, derivative_count,
                       nonzero_field, budget, rtol, atol):
    """Shared radial mass/sensitivity engine on a caller's fixed coordinate.

    coefficients returns (r, dr/dcoordinate, A, a, rho, source_partials).
    A shape partial includes its Jacobian derivative in source_partials.
    Forward and adjoint routes share the metric, not sensitivity values.
    """
    state = np.zeros(2+derivative_count)
    state[1] = centre_lapse
    pieces = []
    samples = 0
    sampled_minimum_A = 1.
    def evaluate(cell, point, mass):
        nonlocal samples, sampled_minimum_A
        row = coefficients(cell,point,mass)
        samples += 1
        _,jacobian,A,a,rho,partials = row
        if (not np.all(np.isfinite([jacobian,A,a,rho])) or jacobian <= 0 or A <= 0
                or len(partials) != derivative_count or not np.all(np.isfinite(partials))):
            raise ValueError("invalid radial coefficient or derivative source")
        sampled_minimum_A = min(sampled_minimum_A,A)
        return row
    for cell in range(len(breaks)-1):
        def rhs(point, values):
            _,jacobian,_,a,rho,partials = evaluate(cell,point,values[0])
            budget.check(len(values))
            return jacobian*np.r_[rho,a,partials-a*values[2:]]
        solution = solve_ivp(rhs,(breaks[cell],breaks[cell+1]),state,method='DOP853',
                             rtol=rtol,atol=atol,dense_output=True)
        if not solution.success or not np.all(np.isfinite(solution.y)):
            raise ValueError("radial mass-gradient integration failed")
        state = solution.y[:,-1].copy()
        pieces.append(solution)
    mass, lapse = float(state[0]), float(state[1])
    if mass < 0 or (mass == 0 and nonzero_field):
        raise ValueError("positive reference mass is unresolved")
    nodes, weights = leggauss(32)
    adjoint = np.zeros(derivative_count)
    for cell, solution in enumerate(pieces):
        lo, hi = breaks[cell], breaks[cell+1]
        points = (lo+hi)/2+(hi-lo)*nodes/2
        values = solution.sol(points)
        for j, point in enumerate(points):
            _,jacobian,_,_,_,partials = evaluate(cell,point,values[0,j])
            budget.check(derivative_count)
            adjoint += (hi-lo)*weights[j]*jacobian*np.exp(values[1,j]-lapse)*partials/2
    gradient = state[2:]
    discrepancy = float(np.max(abs(gradient-adjoint)))
    if discrepancy > GRADIENT_ABSOLUTE_TOLERANCE:
        raise ValueError("forward and adjoint mass gradients exceed the declared tolerance")
    return dict(mass=mass,lapse=lapse,gradient=gradient.copy(),adjoint_gradient=adjoint,
                gradient_discrepancy=discrepancy,field_sites=samples,
                sampled_minimum_A=sampled_minimum_A,adjoint_order=32)


def variational_cap_reference(z_nodes, chi, canonical_density, *, alpha, T_x,
                              incoming_w_plus, budget, rtol=2e-11, atol=2e-13):
    """One fixed-radius weak RHS, with forward and adjoint mass gradients.

    k=A0*p is the canonical density used in the z coordinate; it equals
    chi_T only at alpha=0. The returned coefficient arrays hold centre-zero
    continuous P1 fields. All parameter derivatives keep the mesh and alpha
    fixed. Solver tolerances and gradient discrepancies are not error bounds.
    """
    budget = _budget(budget)
    z = _vector(z_nodes, "z_nodes").copy()
    chi = _vector(chi, "chi").copy()
    k = _vector(canonical_density, "canonical_density").copy()
    alpha = _scalar(alpha, "alpha", nonnegative=True)
    tx = _scalar(T_x, "T_x", positive=True)
    incoming = _scalar(incoming_w_plus, "incoming_w_plus")
    rtol, atol = _scalar(rtol, "rtol", positive=True), _scalar(atol, "atol", positive=True)
    if (not 2 <= len(z) <= 129 or z[0] != 0 or np.any(np.diff(z) <= 0)
            or chi.shape != z.shape or k.shape != z.shape):
        raise ValueError("require 2..129 ordered nodes from zero and matching field arrays")
    if chi[0] != 0 or k[0] != 0:
        raise ValueError("regular centre requires chi(0)=k(0)=0")
    if rtol < 64*np.finfo(float).eps:
        raise ValueError("requested radial tolerance is below the reference floating-point domain")
    before = budget.evaluations
    certificate = _metric_domain(z, chi, k, alpha, budget)
    m = len(z)-1
    widths = np.diff(z)
    if not 0 < float(np.tanh(z[-1])) < 1:
        raise ValueError("wall radius is not representable inside 0<R<1")
    if alpha > 0 and np.sqrt(alpha/6) == 0:
        raise ValueError("reference alpha scale underflows")

    def coefficients(cell, point, mass):
        budget.check(32+16*m)
        if not z[cell] <= point <= z[cell+1]:
            raise ValueError("radial evaluation escaped its declared element")
        r, A0 = float(np.tanh(point)), float(1/np.cosh(point)**2)
        # Basis values and derivatives for free nodes; centre is eliminated.
        left, right, width = z[cell], z[cell+1], widths[cell]
        basis = np.zeros(m)
        slopes = np.zeros(m)
        intercepts = np.zeros(m)
        if cell > 0:
            basis[cell-1] = (right-point)/width
            slopes[cell-1], intercepts[cell-1] = -1/width, right/width
        basis[cell] = (point-left)/width
        slopes[cell], intercepts[cell] = 1/width, -left/width
        if point == 0:
            # f_theta has a finite nonzero centre limit, but all radial
            # mass/lapse and sensitivity RHS values vanish at this endpoint.
            return r,A0,1.,0.,0.,np.zeros(2*m)
        f_basis = basis/r
        p_basis = basis/A0
        h_basis = slopes*(point*_positive_gap(point)/A0)-intercepts/r
        f, p, h = float(f_basis@chi[1:]), float(p_basis@k[1:]), float(h_basis@chi[1:])
        potential = reference_potential(np.array([f]), alpha, budget=budget)
        P, U = float(potential['Pdiff'][0]), float(potential['Pdiff_f'][0])
        A = A0-2*alpha*mass/r
        if not np.isfinite(A) or A <= 0:
            raise ValueError("numerical radial solution left the positive metric domain")
        square = p*p+h*h
        a = alpha*square/(2*r)
        rho = r*r*P/2+A*square/4
        partials = np.r_[r*r*U*f_basis/2+A*h*h_basis/2, A*p*p_basis/2]
        return r,A0,A,a,rho,partials

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            centre_f = float(chi[1]/z[1])
            centre_q = np.sqrt(alpha/6)*centre_f
            centre_lapse = -float(np.log1p(2*np.sinh(centre_q/2)**2))
            reference = _radial_variations(
                z,coefficients,centre_lapse=centre_lapse,derivative_count=2*m,
                nonzero_field=bool(np.any(chi != 0) or np.any(k != 0)),
                budget=budget,rtol=rtol,atol=atol)
            mass,lapse = reference['mass'],reference['lapse']
            gradient,adjoint = reference['gradient'],reference['adjoint_gradient']
            gradient_discrepancy = reference['gradient_discrepancy']
            samples = reference['field_sites']
            sampled_minimum_A = reference['sampled_minimum_A']
            gc, gk = gradient[:m], gradient[m:]
            budget.check(32*m+32)
            diagonal = np.r_[(widths[:-1]+widths[1:])/3,widths[-1]/3]
            off = widths[1:]/6
            band = np.zeros((3,m))
            band[1], band[0,1:], band[2,:-1] = diagonal, off, off
            scale = float(np.exp(lapse))
            velocity = 2*scale*solve_banded((1,1),band,gk)
            R, A0R = float(np.tanh(z[-1])), float(1/np.cosh(z[-1])**2)
            A = A0R-2*alpha*mass/R
            B = A*scale
            forcing = -2*scale*gc.copy()
            forcing[-1] += B*R*incoming-velocity[-1]
            k_rate_T = solve_banded((1,1),band,forcing)
            mass_rate = float(tx*(gc@velocity+gk@k_rate_T))
            wall_velocity = float(velocity[-1])
            flux = wall_velocity*(B*R*incoming-wall_velocity)/(2*scale)
            p_trace = k[-1]/A0R
            slope_wall = (chi[-1]-chi[-2])/widths[-1]
            q_trace = slope_wall/A0R
            f_wall = chi[-1]/R
            h_trace = (slope_wall-chi[-1]/z[-1]+_positive_gap(z[-1])*chi[-1])/A0R
            p_weak = wall_velocity/B
            q_weak = f_wall+R*incoming-p_weak
            result = dict(
                z_nodes=z,chi=chi,canonical_density=k,alpha=alpha,radius=R,radius_rate=0.,T_x=tx,
                mass=mass,lapse=lapse,gradient=dict(chi=gc.copy(),canonical_density=gk.copy()),
                adjoint_gradient=dict(chi=adjoint[:m],canonical_density=adjoint[m:]),
                gradient_discrepancy=gradient_discrepancy,
                gradient_absolute_tolerance=GRADIENT_ABSOLUTE_TOLERANCE,
                rates=dict(chi=np.r_[0.,tx*velocity],canonical_density=np.r_[0.,tx*k_rate_T]),
                mass_rate=mass_rate,physical_weak_flux=float(flux),balance_defect=float(mass_rate-tx*flux),
                incoming_flux_bound=float(tx*(B*R*incoming)**2/(8*scale)),
                wall=dict(A=A,B=B,f=float(f_wall),velocity=wall_velocity,
                          p_weak=float(p_weak),q_weak=float(q_weak),p_trace=float(p_trace),q_trace=float(q_trace),
                          momentum_discrepancy=float(p_trace-p_weak),
                          normal_discrepancy=float(h_trace-R*incoming+p_weak),
                          incoming_discrepancy=float((p_trace+h_trace-R*incoming)/R)),
                centre=dict(f=centre_f,Pi=float(k[1]/z[1]),Phi=0.,delta=centre_lapse),
                metric_certificate=certificate,sampled_minimum_A=float(sampled_minimum_A),
                requested_tolerances=dict(rtol=rtol,atol=atol),adjoint_order=32,
                coefficient_evaluations=budget.evaluations-before,field_sites=samples,
                numerical_error_enclosure=False,
                domain="fixed-radius primitive-conforming weak reference; no time evolution",
            )
    except (FloatingPointError,np.linalg.LinAlgError) as error:
        raise ValueError("unresolved variational cap arithmetic") from error
    def finite(value):
        if isinstance(value,dict):
            return all(finite(v) for v in value.values())
        return isinstance(value,str) or bool(np.all(np.isfinite(value)))
    if not finite(result):
        raise ValueError("nonfinite variational cap reference")
    budget.check()
    return result
