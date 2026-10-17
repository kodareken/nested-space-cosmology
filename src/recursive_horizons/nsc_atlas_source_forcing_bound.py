"""Actual prescribed source relative to vacuum, from certified pulse bounds.

This is source forcing, not interpolation error. The zero-left-jet source
constraint is enclosed for the whole coupling interval 0<=alpha<=alpha_max.
Every returned number is an exact rational bound; fixed dyadic rounding
keeps the certificate small without using floating-point arithmetic.
"""
from fractions import Fraction as Q

from .nsc_atlas_transmission_bound import _rational

SCALE=2**128


def _up(value):
    value=Q(value)*SCALE
    return Q(-((-value.numerator)//value.denominator),SCALE)


def _down(value):
    value=Q(value)*SCALE
    return Q(value.numerator//value.denominator,SCALE)


def source_forcing_bound(*,left,right,alpha_max,g,g_x,g_xx,budget):
    """Uniform absolute forcing for the original quiet-initial source ODE."""
    if not callable(getattr(budget,'check',None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(2048)
    left=_rational(left,'left',positive=True);right=_rational(right,'right',positive=True)
    al=_rational(alpha_max,'maximum coupling')
    g,gx,gxx=(_rational(v,k) for v,k in ((g,'g'),(g_x,'g_x'),(g_xx,'g_xx')))
    if not left<right<Q(4,5):raise ValueError('source cell must lie inside the horizon')
    width=right-left;smin=Q(1,5);vmax=Q(3,2)
    # h''=(g*g''-2*g'^2)/(6*sJ)+g*g'*vJ/(3*sJ^2), with h=h'=0 at left.
    acceleration=_up((g*gxx+2*gx*gx)/(6*smin)+g*gx*vmax/(3*smin*smin))
    h=_up(width*width*acceleration/2);hx=_up(width*acceleration)
    d0,d1=_up(al*h),_up(al*hx)
    slo,shi=_down(left-d0),_up(right+d0);vhi=_up(1+d1)
    if not slo>smin or not vhi<vmax:
        raise ValueError('source bootstrap leaves its common Jordan tube')
    loss=_up(al*g*g/(6*slo*slo));Flo=_down(1-loss)
    if not Q(1,2)<Flo<=1:raise ValueError('source frame loses the proved positive margin')
    rlo=_down(Flo*slo)  # sqrt(F)>=F on 0<F<=1
    shrink=_up(al*g*g/(6*(slo+rlo)))
    er=_up(d0+shrink)
    erx=_up((shrink+shi*d1+al*g*gx/6)/rlo)
    eL=_up(loss/Flo)  # -log(1-loss)<=loss/(1-loss)
    eLx=_up(al*(g*gx/(slo*slo)+g*g*vhi/(slo**3))/(3*Flo))
    product=_up(d0+shi*d1+al*g*gx/6)
    u=_up(g/rlo);gap=_up(al*u**3/18)
    Qbound=_up(gx+gap*(1+erx))
    budget.check()
    return dict(
        alpha_max=al,duration=width,dyadic_bits=128,
        radial=dict(h_xx=acceleration,h=h,h_x=hx,rJ_minus_x=d0,rJ_x_minus_one=d1),
        tube=dict(rJ_lower=slo,rJ_upper=shi,absolute_rJ_x_upper=vhi,F_lower=Flo,
                  r_lower=rlo,r_upper=shi,absolute_r_x_upper=_up(1+erx)),
        forcing=dict(r=er,r_x=erx,log_C=eL,log_C_x=eLx,c=_up(g),c_x=_up(gx),
                     X=_up(g),Q=Qbound,r_times_r_x_minus_x=product),
        frame=dict(one_minus_F=loss,rJ_minus_r=shrink,canonical_gap=gap),
        scope='actual source relative to vacuum on 0<=alpha<=alpha_max; no coupled evolution or interpolation-error claim')
