"""Exact-model wall momentum projection bound for the fixed P1 cap.

This encloses the spatial Galerkin projection only. It excludes numerical
radial/gradient errors, the normal boundary discrepancy and mesh motion.
"""
import numpy as np
from flint import arb, ctx

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_variational_reference import _exact_ball, _metric_domain
from .nsc_atlas_polar_cap import _scalar, _vector


def _endpoint_record(value, *, upper):
    endpoint = value.upper() if upper else value.lower()
    if not endpoint.is_finite():
        raise ValueError('projection bound is not finite')
    mantissa, exponent = endpoint.man_exp()
    rounded = float(endpoint)
    if endpoint != 0:
        rounded = float(np.nextafter(rounded,np.inf if upper else -np.inf))
    if not np.isfinite(rounded):
        raise ValueError('projection bound display is not representable')
    if not upper and rounded <= 0:
        raise ValueError('positive wall lower bound display is not representable')
    return dict(mantissa=str(mantissa),exponent=int(exponent),
                **{'float_up' if upper else 'float_down':rounded})


def cap_momentum_projection_bound(z_nodes, chi, canonical_density, *, alpha, budget):
    """Enclose |V_W-B_W p_W| and |V_W/B_W-p_W| for exact mass gradients.

    Inputs are exact binary nodal values defining continuous P1 chi,k in z.
    The bound uses real-arithmetic field norms and the existing radial
    constraints; it executes no numerical radial or time integration.
    """
    budget = _budget(budget)
    z = _vector(z_nodes,'z_nodes').copy()
    chi = _vector(chi,'chi').copy()
    k = _vector(canonical_density,'canonical_density').copy()
    alpha = _scalar(alpha,'alpha',nonnegative=True)
    if (not 2 <= len(z) <= 129 or z[0] != 0 or np.any(np.diff(z) <= 0)
            or chi.shape != z.shape or k.shape != z.shape or chi[0] != 0 or k[0] != 0):
        raise ValueError('require ordered nodes and regular centre-zero P1 fields')
    before = budget.evaluations
    certificate = _metric_domain(z,chi,k,alpha,budget)
    budget.check(192*len(z)+128)
    with ctx.workprec(160):
        zb, cb, kb = [[_exact_ball(v) for v in values] for values in (z,chi,k)]
        Z, al = zb[-1],_exact_ball(alpha)
        R,A0 = Z.tanh(),1/Z.cosh()**2
        widths = [b-a for a,b in zip(zb[:-1],zb[1:])]
        lc = max(abs((b-a)/h).upper() for a,b,h in zip(cb[:-1],cb[1:],widths))
        lk = max(abs((b-a)/h).upper() for a,b,h in zip(kb[:-1],kb[1:],widths))
        centre_f = cb[1]/zb[1]
        delta0 = -((al/6).sqrt()*centre_f).cosh().log()
        floor = certificate['A_lower']
        eta = arb(int(floor['mantissa']))*arb(2)**floor['exponent']
        lower_B = eta*delta0.exp()
        if not lower_B > 0:
            raise ValueError('wall lapse/metric lower bound is unresolved')
        if alpha == 0 or lk == 0:
            second_derivative = arb(0)
        else:
            F, hmax = lc*Z/R,lc/A0+lc*Z/R
            argument = (al/6).sqrt()*F
            if lc == 0:
                potential, Umax = arb(0),arb(0)
            else:
                sinc = argument.sinh()/argument
                potential = F*F*sinc*sinc*(argument.cosh()**2+1)/2
                Umax = 2*F*sinc*argument.cosh()**3
            pi_max = lk*Z/(R*A0)
            # The first-element regular Phi kernel has positive r-series
            # coefficients. Other elements are bounded using r>=r_1.
            phi_first = lc*abs((R/A0-Z)/(R*R))
            global_phi = max(phi_first.upper(),(hmax/zb[1].tanh()).upper())
            # On each element Phi=a_e*K(r)-b_e/r^2, where the positive
            # kernel K increases with r. The first intercept is exactly
            # zero; every other element has a positive lower radius.
            local_phi = []
            for j,width in enumerate(widths):
                slope = (cb[j+1]-cb[j])/width
                intercept = cb[j]-slope*zb[j]
                hi = zb[j+1].tanh()
                ahi = 1/zb[j+1].cosh()**2
                bound = abs(slope)*abs((hi/ahi-zb[j+1])/(hi*hi))
                if j > 0:
                    bound += abs(intercept)/zb[j].tanh()**2
                local_phi.append(bound.upper())
            phi_max = min(global_phi,max(local_phi))
            square = pi_max*pi_max+phi_max*phi_max
            cbound = potential/2+square/4
            nu = cbound/3
            Cmax = R*(2*nu+potential+4*nu*R*R/A0)
            Crmax = (2*cbound+4*nu+potential+R*Umax*phi_max
                     +4*R*R*cbound/A0+8*nu*R**4/A0**2)
            amax = al*R*square/2
            delta_upper = delta0+al*R*R*square/4
            beta_first = al*delta_upper.exp()*Cmax
            beta_second = al*delta_upper.exp()*(amax*Cmax+Crmax)
            second_derivative = beta_second*lk*Z+2*beta_first*lk
        largest_width = max(h.upper() for h in widths)
        interpolation = largest_width**2*second_derivative/8
        velocity = 3*interpolation
        momentum = velocity/lower_B
        if not (second_derivative >= 0 and velocity >= 0 and momentum >= 0):
            raise ValueError('nonnegative projection bounds are unresolved')
        result = dict(
            z_nodes=z,chi=chi,canonical_density=k,alpha=alpha,precision_bits=160,
            metric_certificate=certificate,
            B_wall_lower=_endpoint_record(lower_B,upper=False),
            g_second_derivative_upper=_endpoint_record(second_derivative,upper=True),
            interpolation_upper=_endpoint_record(interpolation,upper=True),
            wall_velocity_projection_upper=_endpoint_record(velocity,upper=True),
            wall_momentum_projection_upper=_endpoint_record(momentum,upper=True),
            coefficient_evaluations=budget.evaluations-before,
            scope='exact fixed-cap Galerkin spatial projection only',
            includes_radial_solver_error=False,includes_gradient_error=False,
            includes_normal_trace_error=False,includes_mesh_motion=False,
        )
    budget.check()
    return result
