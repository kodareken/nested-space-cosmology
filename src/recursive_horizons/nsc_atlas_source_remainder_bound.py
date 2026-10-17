"""Quadratic source-geometry remainders from the original uniform source tube.

The invariant trace g is independent of alpha. Its Jordan radius obeys
rJ=x+alpha*h with quiet h,h_x anchors. Only existing ceilings are consumed;
no source integration or amplitude-derivative solution is evaluated.
"""
from fractions import Fraction as Q

from .nsc_atlas_source_forcing_bound import _up
from .nsc_atlas_transmission_bound import _rational


def source_remainder_comparison(*, g, g_x, g_xx, h, h_x, duration, right,
                                r_min, F_min, budget):
    """Exact rational sufficient constants for |slot-alpha*slot1|<=alpha^2*C."""
    budget.check(512)
    values=dict(g=g,g_x=g_x,g_xx=g_xx,h=h,h_x=h_x,duration=duration,
                right=right,r_min=r_min,F_min=F_min)
    v={k:_rational(x,k) for k,x in values.items()}
    g,gx,gxx,H,Hx,dt,xhi,r,F=(v[k] for k in
        ('g','g_x','g_xx','h','h_x','duration','right','r_min','F_min'))
    if not 0<r<=xhi<Q(4,5) or dt<=0 or not Q(1,2)<F<=1:
        raise ValueError('require original positive Jordan/frame tube and source cell')
    A=(g*gxx+2*gx*gx)/6
    B=g*gx/3
    H2=A*H/r**2+B*(Hx/r**2+2*H/r**3)
    Ch,Chx=dt*dt*H2/2,dt*H2
    # sqrt(F)>=F_min: positive rational bounds on
    # f(u)=1/(1+sqrt(1-u)) and f'(u) along 0<=u<=1-F_min.
    fmax=1/(1+F)
    fpmax=1/(2*F*(1+F)**2)
    Cr=Ch+g*g/(6*r*r)*(H*fmax+g*g*fpmax/(6*r))
    CP=Ch+xhi*Chx+H*Hx
    CL=g*g*H/(3*r**3)+g**4/(72*r**4*F**2)
    KQ=2*gx*H/r**3+g*Hx/r**3+3*g*H/r**4
    CLx=g/(3*F)*(KQ+(gx/r**2+g/r**3)*g*g/(6*r*r))
    budget.check()
    return dict(coefficients=dict(r=Cr,log_C=CL,log_C_x=CLx,r_times_r_x_minus_x=CP),
                comparison=dict(h_xx_difference=H2,h_difference=Ch,h_x_difference=Chx,
                                frame_quotient=fmax,frame_quotient_derivative=fpmax,
                                lapse_slope_numerator_difference=KQ))


def source_remainder_bound(parameters, source, *, budget):
    if not callable(getattr(budget,'check',None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(2048)
    left,right,alpha=(Q(parameters[k]) for k in ('left','right','alpha_max'))
    dt=right-left
    if (dt!=Q(source['duration']) or alpha!=Q(source['alpha_max']) or alpha<=0
            or not 0<left<right<Q(4,5)):
        raise ValueError('require original source domain and coupling')
    rmin=Q(source['tube']['rJ_lower'])
    if not 0<rmin<=left:
        raise ValueError('Jordan floor must also bound the reference radius x')
    values={k:Q(parameters[k]) for k in ('g','g_x','g_xx')}
    values.update(h=Q(source['radial']['h']),h_x=Q(source['radial']['h_x']),
        duration=dt,right=right,r_min=rmin,F_min=Q(source['tube']['F_lower']))
    result=source_remainder_comparison(**values,budget=budget)
    coefficients={k:_up(v) for k,v in result['coefficients'].items()}
    comparison={k:_up(v) for k,v in result['comparison'].items()}
    budget.check(1024)
    return dict(alpha_max=alpha,left=left,right=right,duration=dt,inputs=values,
        coefficients=coefficients,comparison=comparison,dyadic_bits=128,
        maximum_amplitude_remainders={k:alpha*alpha*v for k,v in coefficients.items()},
        first_coefficients=dict(r='h0-g^2/(12*x)',log_C='-g^2/(6*x^2)',
            log_C_x='-g/3*(g_x/x^2-g/x^3)',r_times_r_x_minus_x='h0+x*h0_x-g*g_x/6'),
        scope='uniform quadratic remainders of four source geometry slots in original affine '
              'label; h0 is the quiet first Jordan response; no source integration, '
              'pointwise h0 evaluation, coupled bulk/wall remainder or observer error')
