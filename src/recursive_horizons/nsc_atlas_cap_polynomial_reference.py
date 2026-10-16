"""Regular polynomial cap mass gradients and one moving-domain weak RHS.

The radial solver is shared with the P1 reference. The radius derivative is
integrated independently of the shape identity; boundary trace errors remain
explicit. This module performs no physical time evolution.
"""
import numpy as np
from flint import ctx

from . import nsc_atlas_cap_variational_reference as radial
from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_geometric_reference import _positive_gap
from .nsc_atlas_cap_radial_reference import reference_potential
from .nsc_atlas_polar_cap import _scalar, _vector
from .nsc_atlas_cap_polynomial_basis import odd_legendre_operators, odd_legendre_values
from .nsc_atlas_cap_background_jets import RetainedCapBackground, SMALL_Z
from numpy.polynomial.legendre import leggauss


SHAPE_ABSOLUTE_TOLERANCE = 2e-8


def _background_projection(background, T, Z, count, breaks, budget):
    """Coefficients of projected background time jets on existing pieces."""
    nodes,weights = leggauss(19)
    budget.check(19*19)
    result = np.zeros((2,count))
    for left,right in zip(breaks[:-1],breaks[1:]):
        half = (right-left)/2
        for x,weight in zip(left+half*(nodes+1),weights):
            jet = background.reduced_jets(T,Z*float(x),budget=budget)
            phi = odd_legendre_values(count,float(x),budget=budget)['value']
            result += half*weight*np.array([jet['k'],jet['k_T']])[:,None]*phi
            budget.check(8*count)
    return result*(2*np.arange(1,2*count,2)+1)


def polynomial_cap_reference(Z, chi_coefficients, density_coefficients, *,
                             alpha, T_x, radius_rate, incoming_w_plus, budget,
                             rtol=2e-11, atol=2e-13, basis='monomial',
                             background=None, background_time=None):
    """One RHS for regular odd polynomials in xi=z/Z; no time integration.

    radius_rate is dR/dx; T_x converts cap time to that same wall label.
    Use 1..3 monomials or 1..16 Legendre modes in the same field space.
    With a RetainedCapBackground, Legendre coefficients describe departures;
    the explicit-time mass remainder is reported, not silently removed.
    All numerical discrepancy thresholds concern these dimensionless inputs;
    they are acceptance checks, not uniform numerical error enclosures.
    """
    budget = _budget(budget)
    Z = _scalar(Z,'Z',positive=True)
    chi = _vector(chi_coefficients,'chi_coefficients').copy()
    k = _vector(density_coefficients,'density_coefficients').copy()
    alpha = _scalar(alpha,'alpha',nonnegative=True)
    tx = _scalar(T_x,'T_x',positive=True)
    rx = _scalar(radius_rate,'radius_rate')
    incoming = _scalar(incoming_w_plus,'incoming_w_plus')
    rtol,atol = _scalar(rtol,'rtol',positive=True),_scalar(atol,'atol',positive=True)
    if basis not in ('monomial','legendre'):
        raise ValueError('basis must be monomial or legendre')
    if background is not None:
        if not isinstance(background,RetainedCapBackground) or basis != 'legendre':
            raise ValueError('retained background requires Legendre departure coordinates')
        time = _scalar(background_time,'background_time')
    elif background_time is not None:
        raise ValueError('background_time requires a retained background')
    maximum = 3 if basis == 'monomial' else 16
    if not 1 <= len(chi) <= maximum or k.shape != chi.shape:
        raise ValueError('reference requires matching coefficients within its declared mode count')
    if rtol < 64*np.finfo(float).eps:
        raise ValueError('requested radial tolerance is below the reference floating-point domain')
    if alpha > 0 and np.sqrt(alpha/6) == 0:
        raise ValueError('reference alpha scale underflows')
    before = budget.evaluations
    count = len(chi)
    powers = np.arange(1,2*count,2,dtype=float)
    if basis == 'legendre':
        operators = odd_legendre_operators(count,Z,budget=budget)
        slope_bounds = operators['derivative_bound_xi']
        dilation = operators['dilation']
        centre_weights = operators['centre_derivative_xi']
        gram = operators['gram']
    else:
        slope_bounds = powers
        dilation = np.diag(powers)
        centre_weights = np.eye(1,count).ravel()
        gram = Z/(powers[:,None]+powers[None,:]+1)
    budget.check(256*count+256)
    breaks = np.array([0.,1.])
    if background is None:
        with ctx.workprec(160):
            zb,al = radial._exact_ball(Z),radial._exact_ball(alpha)
            lc = sum(int(p)*abs(radial._exact_ball(v)) for p,v in zip(slope_bounds,chi))/zb
            lk = sum(int(p)*abs(radial._exact_ball(v)) for p,v in zip(slope_bounds,k))/zb
            certificate = radial._metric_floor_from_slopes(
                zb,lc,lk,al,scope='exact radial solution for declared odd polynomials; numerical ODE error excluded')
    else:
        from .nsc_atlas_cap_background_bound import background_metric_bound
        certificate = background_metric_bound(background,time,Z,chi,k,alpha=alpha,budget=budget)
        knots = background.profile.breakpoints
        mapped = np.r_[(knots-time)/Z,(time-knots)/Z,SMALL_Z/Z]
        breaks = np.unique(np.r_[0.,mapped[(mapped>0)&(mapped<1)],1.])
        a_h,b_h = _background_projection(background,time,Z,count,breaks,budget)
    R,A0R = float(np.tanh(Z)),float(1/np.cosh(Z)**2)
    if not 0 < R < 1 or not A0R > 0:
        raise ValueError('wall radius is not representable inside 0<R<1')
    condition = float(np.linalg.cond(gram))
    if not np.isfinite(condition) or condition*np.finfo(float).eps > 1e-9:
        raise ValueError('polynomial Gram conditioning is outside the reference domain')

    def fields(xi):
        z = Z*xi
        r,A0 = float(np.tanh(z)),float(1/np.cosh(z)**2)
        if xi == 0:
            return r,A0,np.zeros(count),np.zeros(count),np.zeros(count)
        if basis == 'legendre':
            values = odd_legendre_values(count,xi,budget=budget)
            Q,Qt = values['quotient'],values['quotient_derivative_t']
            fb = Q*z/(Z*r)
            pb = values['value']/A0
            hb = (2*xi*xi*Qt+z*_positive_gap(z)*Q)/(Z*A0)
        else:
            base = xi**(powers-1)/Z
            fb = base*z/r
            pb = xi**powers/A0
            hb = base*((powers-1)+z*_positive_gap(z))/A0
        return r,A0,fb,pb,hb

    nonzero_background_seen = False
    def coefficients(cell,xi,mass):
        nonlocal nonzero_background_seen
        budget.check(40+20*count)
        if not 0 <= cell < len(breaks)-1 or not breaks[cell] <= xi <= breaks[cell+1]:
            raise ValueError('polynomial radial evaluation escaped the unit interval')
        r,A0,fb,pb,hb = fields(xi)
        J = Z*A0
        if xi == 0:
            return r,J,1.,0.,0.,np.zeros(2*count+1+(background is not None))
        fd,pd,hd = float(fb@chi),float(pb@k),float(hb@chi)
        f,p,h = fd,pd,hd
        if background is not None:
            bg = background.reduced_jets(time,Z*xi,budget=budget)
            nonzero_background_seen |= any(bg[name] != 0 for name in ('f','p','h'))
            f,p,h = f+bg['f'],p+bg['p'],h+bg['h']
        potential = reference_potential(np.array([f]),alpha,budget=budget)
        P,U = float(potential['Pdiff'][0]),float(potential['Pdiff_f'][0])
        A = A0-2*alpha*mass/r
        if not np.isfinite(A) or A <= 0:
            raise ValueError('numerical radial solution left the positive metric domain')
        square = p*p+h*h
        a = alpha*square/(2*r)
        rho = r*r*P/2+A*square/4
        partials = np.r_[r*r*U*fb/2+A*h*hb/2,A*p*pb/2]
        # Shape partials hold xi and coefficients fixed, while the running
        # mass is held fixed here. Its derivative is the shared -a*G term.
        rZ,A0Z = xi*A0,-2*r*xi*A0
        fZ,pZ = -fd*rZ/r,2*r*xi*pd
        # q_Z-f_Z cancels two nonzero centre limits. Use their exact
        # gap identity so the O(z^2) normal shape derivative stays resolved.
        hZ = (-1/Z+2*r*xi)*hd+xi*(2*r-_positive_gap(Z*xi))*fd
        if background is not None:
            fZ,pZ,hZ = fZ+xi*bg['f_z'],pZ+xi*bg['p_z'],hZ+xi*bg['h_z']
        AZ = A0Z+2*alpha*mass*rZ/(r*r)
        shape_source = (r*rZ*P+r*r*U*fZ/2+AZ*square/4
                        +A*(p*pZ+h*hZ)/2+(1/Z-2*r*xi)*rho)
        sources = np.r_[partials,shape_source]
        if background is not None:
            time_source = r*r*U*bg['f_T']/2+A*(p*bg['p_T']+h*bg['h_T'])/2
            sources = np.r_[sources,time_source]
        return r,J,A,a,rho,sources

    try:
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            centre_f = float(centre_weights@chi/Z)
            centre_pi = float(centre_weights@k/Z)
            if background is not None:
                axis = background.reduced_jets(time,0.,budget=budget)
                centre_f += axis['f']
                centre_pi += axis['f_T']
            cq = np.sqrt(alpha/6)*centre_f
            delta0 = -float(np.log1p(2*np.sinh(cq/2)**2))
            reference = radial._radial_variations(
                breaks,coefficients,centre_lapse=delta0,
                derivative_count=2*count+1+(background is not None),
                nonzero_field=bool(np.any(chi) or np.any(k)),
                budget=budget,rtol=rtol,atol=atol)
            mass,lapse = reference['mass'],reference['lapse']
            if mass == 0 and nonzero_background_seen:
                raise ValueError('positive retained-background mass is unresolved')
            gradient = reference['gradient']
            gc,gk,shape = gradient[:count],gradient[count:2*count],float(gradient[2*count])
            explicit_time = float(gradient[-1]) if background is not None else 0.
            _,_,A,_,rho,_ = coefficients(len(breaks)-2,1.,mass)
            _,_,fb,pb,hb = fields(1.)
            f,p,h = float(fb@chi),float(pb@k),float(hb@chi)
            if background is not None:
                wall_bg = background.reduced_jets(time,Z,budget=budget)
                f,p,h = f+wall_bg['f'],p+wall_bg['p'],h+wall_bg['h']
            scale = float(np.exp(lapse))
            B = A*scale
            if not abs(rx) < tx*B:
                raise ValueError('prescribed wall speed is not timelike')
            dilation_chi,dilation_k = dilation@chi,dilation@k
            shape_from_dilation = A0R*rho-(gc@dilation_chi+gk@dilation_k)/Z
            shape_discrepancy = shape-shape_from_dilation
            if abs(shape_discrepancy) > SHAPE_ABSOLUTE_TOLERANCE:
                raise ValueError('independent radial shape and dilation identity disagree')
            budget.check(64*count+64)
            normalized_residuals = []
            def gram_solve(rhs):
                value = np.linalg.solve(gram,rhs)
                denominator = np.linalg.norm(gram,np.inf)*np.linalg.norm(value,np.inf)+np.linalg.norm(rhs,np.inf)
                residual = float(np.linalg.norm(gram@value-rhs,np.inf)/denominator) if denominator else 0.
                if not np.all(np.isfinite(value)) or residual > 64*np.finfo(float).eps:
                    raise ValueError('polynomial Gram solve residual is unresolved')
                normalized_residuals.append(residual)
                return value
            V = 2*scale*gram_solve(gk)
            departure_velocity = V if background is None else V-a_h
            wall_velocity = float(np.sum(departure_velocity))
            if background is not None:
                wall_velocity += wall_bg['k']
            drive = B*R*incoming
            fixed_k_rate = gram_solve(-2*scale*gc+drive-wall_velocity)
            if background is not None:
                fixed_k_rate -= b_h
            zx = rx/A0R
            chi_rate = tx*departure_velocity+zx*dilation_chi/Z
            k_rate = tx*fixed_k_rate+zx*dilation_k/Z
            mass_rate = float(tx*explicit_time+gc@chi_rate+gk@k_rate+shape*zx)
            weak_flux = wall_velocity*(drive-wall_velocity)/(2*scale)
            moving_density = rx*rho
            remainder = 0.
            if background is not None:
                remainder = float(explicit_time-gc@a_h-gk@b_h
                    -(wall_bg['k']-np.sum(a_h))*(drive-wall_velocity)/(2*scale))
            p_weak = wall_velocity/B
            h_weak = R*incoming-p_weak
            normal_discrepancy = h-h_weak
            scalar_tangent = (tx*wall_velocity+rx*h)/R
            weak_tangent = (tx*wall_velocity+rx*h_weak)/R
            # The existing hybrid interface combines its imposed incoming
            # characteristic with the reconstructed outgoing characteristic.
            # That pair is not in general (p_weak,h_weak).
            outgoing = (p-h)/R
            p_interface = R*(incoming+outgoing)/2
            h_interface = R*(incoming-outgoing)/2
            transmitted_tangent = (tx*B*p_interface+rx*h_interface)/R
            interface_flux = A*B*p_interface*h_interface/2
            density_gap = A*(p*p+h*h-p_interface*p_interface-h_interface*h_interface)/4
            flux_gap = weak_flux-interface_flux
            result = dict(
                Z=Z,chi_coefficients=chi,density_coefficients=k,powers=powers,basis=basis,
                alpha=alpha,T_x=tx,radius=R,radius_rate=rx,Z_rate=zx,
                mass=mass,lapse=lapse,
                gradient=dict(chi=gc.copy(),canonical_density=gk.copy(),Z=shape),
                adjoint_gradient=dict(chi=reference['adjoint_gradient'][:count].copy(),
                                      canonical_density=reference['adjoint_gradient'][count:2*count].copy(),
                                      Z=float(reference['adjoint_gradient'][2*count])),
                gradient_discrepancy=reference['gradient_discrepancy'],
                shape_from_dilation=float(shape_from_dilation),shape_discrepancy=float(shape_discrepancy),
                shape_absolute_tolerance=SHAPE_ABSOLUTE_TOLERANCE,
                rates=dict(chi=chi_rate,canonical_density=k_rate),
                mass_rate=mass_rate,physical_weak_flux=float(weak_flux),
                moving_wall_density=float(moving_density),
                balance_defect=float(mass_rate-tx*weak_flux-moving_density),
                interface_flux_trace_gap=float(tx*flux_gap),
                interface_moving_density_gap=float(rx*density_gap),
                interface_balance_defect=float(mass_rate-tx*interface_flux-rx*(rho-density_gap)),
                wall=dict(A=A,B=B,f=f,p_trace=p,q_trace=f+h,p_weak=p_weak,q_weak=f+h_weak,
                          velocity=wall_velocity,rho=rho,momentum_discrepancy=p-p_weak,
                          normal_discrepancy=normal_discrepancy,
                          incoming_discrepancy=(p+h-R*incoming)/R,
                          outgoing_w_minus=outgoing,p_interface=p_interface,q_interface=f+h_interface,
                          weak_characteristic_scalar_tangent=weak_tangent,
                          weak_scalar_tangent_discrepancy=scalar_tangent-weak_tangent,
                          scalar_tangent=scalar_tangent,transmitted_scalar_tangent=transmitted_tangent,
                          scalar_tangent_discrepancy=scalar_tangent-transmitted_tangent),
                centre=dict(f=centre_f,Pi=centre_pi,Phi=0.,delta=delta0),
                metric_certificate=certificate,sampled_minimum_A=reference['sampled_minimum_A'],
                gram_condition=condition,gram_normalized_residual=max(normalized_residuals),
                requested_tolerances=dict(rtol=rtol,atol=atol),adjoint_order=reference['adjoint_order'],
                coefficient_evaluations=budget.evaluations-before,field_sites=reference['field_sites']+1,
                numerical_error_enclosure=False,
                domain='one regular polynomial moving-cap weak RHS; no physical time evolution',
            )
            if background is not None:
                result['domain'] = 'one retained-background and polynomial-departure moving-cap weak RHS; no physical time evolution'
                result['gradient']['T'] = explicit_time
                result['adjoint_gradient']['T'] = float(reference['adjoint_gradient'][-1])
                result['background'] = dict(time=time,radial_intervals=len(breaks)-1,
                    projected_time_chi=a_h.copy(),projected_time_density=b_h.copy(),
                    wall_time_projection_gap=float(wall_bg['k']-np.sum(a_h)))
                result['explicit_time_remainder'] = remainder
                result['remainder_identity_defect'] = float(result['balance_defect']-tx*remainder)
    except (FloatingPointError,np.linalg.LinAlgError) as error:
        raise ValueError('unresolved polynomial reference arithmetic') from error
    def finite(value):
        if isinstance(value,dict):
            return all(finite(v) for v in value.values())
        return isinstance(value,str) or bool(np.all(np.isfinite(value)))
    if not finite(result):
        raise ValueError('nonfinite polynomial cap reference')
    budget.check()
    return result
