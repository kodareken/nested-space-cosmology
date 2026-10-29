"""HM lifetime ordering and regular compact O(4) CDL exclusion tests.

The completed-square test excludes four candidates analytically. A bounded
two-pole search in the reference reports its finite scope, not a proof of
absence or an unconditional selection theorem.
"""
import math,time
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_volume_family import member
from .nsc_atlas_volume_routes import FRONTIER


def exclusion_margin(model):
    """Global minimum of e^(d*s)*(V_phiphi+1.3V+.65V_phi)."""
    d=model.d;c=model.coefficients
    weights=[d/math.sqrt(model.K),(d+2)/math.sqrt(model.K),2*d/math.sqrt(model.K)]
    a,b,t=[float(c[k])*(p*p-13*p/20+13/10) for k,p in zip(('A','B','T'),weights)]
    x=(2*b/(d*t))**(2/(d-2))
    return dict(reduced_global_minimum=a-(d-2)*b*x/d,x_at_minimum=x,
                criterion='G=V_phiphi+13*V/10+13*V_phi/20')


def exclusion_lower_bound(model):
    """Exact rational lower bound on the global certificate, not sampling."""
    d=model.d;c=model.coefficients;K=c['kinetic'];scale=10**12
    n=math.isqrt(K.numerator*scale*scale//K.denominator)
    rootlo,roothi=F(n,scale),F(n+1,scale)
    # p^2-.65p+1.3 is increasing for all three exponents here.
    poly=lambda p:p*p-F(13,20)*p+F(13,10)
    al=c['A']*poly(F(d)/roothi)
    bu=c['B']*poly(F(d+2)/rootlo)
    tl=c['T']*poly(F(2*d)/roothi)
    power=(2*bu/(d*tl))**2
    lo,hi=F(0),F(2)
    assert hi**(d-2)>power
    for _ in range(48):
        mid=(lo+hi)/2
        if mid**(d-2)<power:lo=mid
        else:hi=mid
    lower=al-F(d-2,d)*bu*hi
    decimal_floor=lower.numerator*10**10//lower.denominator
    return dict(lower_bound=f'{decimal_floor//10**10}.{decimal_floor%10**10:010d}',strictly_positive=lower>0)


def pole_shot(model,s0,*,rtol=2e-10):
    """Regular pole to first equator, with local-Hubble rescaling."""
    H=math.sqrt(model.potential(s0)/3);eps=1e-5
    force=model.potential(s0,1)/(model.K*H*H)
    calls=0
    def rhs(t,y):
        nonlocal calls
        calls+=1
        if calls>10000:raise RuntimeError('single pole allowance exhausted')
        s,u,a,b=y
        return u,model.potential(s,1)/(model.K*H*H)-3*b*u/a,b,-a*(model.K*u*u+model.potential(s)/(H*H))/3
    def equator(t,y):return y[3]
    equator.terminal=True;equator.direction=-1
    initial=[s0+force*eps**2/8,force*eps/4,eps-eps**3/6,1-eps**2/2]
    sol=solve_ivp(rhs,(eps,30),initial,events=equator,method='DOP853',rtol=rtol,atol=rtol/100)
    if not sol.success or not len(sol.t_events[0]):raise RuntimeError('equator not reached')
    s,u,a,b=sol.y[:,-1]
    physical_u,physical_a=u*H,a/H
    residual=b*b-1-physical_a**2*(model.K*physical_u**2/2-model.potential(s))/3
    return dict(s=s,derivative=physical_u,radius=physical_a,constraint_residual=residual,calls=calls)


def reference_search():
    model=VolumeParent();points=member(12,2)['stationary_points']
    top,false=points[0]['s'],points[1]['s'];width=top-false
    H=math.sqrt(points[0]['vacuum']/3);calls=0;maximum_constraint=0.;start=time.process_time()
    def residual(z):
        nonlocal calls,maximum_constraint
        calls+=1
        if calls>1600 or time.process_time()-start>120:raise RuntimeError('bounded CDL search allowance exhausted')
        a=pole_shot(model,top+width*z[0]);b=pole_shot(model,false+width*z[1])
        maximum_constraint=max(maximum_constraint,abs(a['constraint_residual']),abs(b['constraint_residual']))
        return [(a['s']-b['s'])/width,(a['derivative']+b['derivative'])/(H*width)]
    rows=[]
    for seed in ((.01,.95),(.05,.7),(.1,.6),(.3,.5),(.5,.3),(1.,.1),(2.,.05),(5.,.01)):
        sol=least_squares(residual,seed,bounds=([1e-5,1e-6],[20,1-1e-6]),max_nfev=100,
                          xtol=2e-9,ftol=2e-9,gtol=2e-9)
        error=residual(sol.x)
        rows.append(dict(seed=list(seed),pole_coordinates=sol.x.tolist(),matching_residual=error,
                         admitted_nonconstant_bounce=bool(np.linalg.norm(error)<1e-8 and sum((sol.x[0],1-sol.x[1]))>1e-3)))
    return dict(rows=rows,residual_calls=calls,maximum_equator_constraint_residual=maximum_constraint,
        nonconstant_bounces_found=sum(r['admitted_nonconstant_bounce'] for r in rows),
        status='Eight bounded seeds approach HM; this is not an exhaustive root count. Reference CDL absence remains unresolved.',
        domain='s_left=top+width*l, 1e-5<=l<=20; s_right=false+width*r, 1e-6<=r<=1-1e-6')


def calculation(*,search=True):
    rows=[]
    for d,j in FRONTIER:
        model=VolumeParent(d,j);m=member(d,j);top,false=m['stationary_points']
        hm=24*math.pi**2*(1/false['vacuum']-1/top['vacuum'])
        margin=exclusion_margin(model)
        certificate=exclusion_lower_bound(model)
        rows.append(dict(d=d,j=j,Lambda_Lref0_squared=m['Lambda_Lref0_squared'],
            HM_exponent_times_g=hm,top_negative_curvature_ratio=-3*top['stiffness']/(model.K*top['vacuum']),
            CDL_exclusion=margin,exact_exclusion_certificate=certificate,
            regular_compact_CDL_excluded=certificate['strictly_positive'],
            CDL_exponent_times_g=None,
            CDL_status='excluded by global completed-square identity' if margin['reduced_global_minimum']>0 else 'absence unresolved; bounded reference search below'))
    return dict(schema='nsc-atlas-volume-lifetime-v1',rows=rows,
        reference_search=reference_search() if search else None,
        HM_order=[dict(d=r['d'],j=r['j']) for r in sorted(rows,key=lambda r:r['HM_exponent_times_g'],reverse=True)],
        conditional_survival_theorem='For HM-only rates A_i(g) exp(-b_i/g), subexponential positive prefactors, common sufficiently small g, positive initial populations and no replenishment, conditioning on survival at late times selects (12,2).',
        unconditional_selection_derived=False,
        remaining_lifetime_question='A possible lower-action CDL or omitted channel in (12,2), prefactors and the population/measure assumptions prevent upgrading the HM ordering to unconditional selection.',
        scope='Regular compact O(4) bounces of the declared Einstein/canonical-volume truncation. Singular instantons and omitted internal modes are not classified.')
