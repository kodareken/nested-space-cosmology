"""Smooth finite-coupling cap characteristics and a common-value wall port.

These components reuse the reduced polar PDE and hybrid wall law. They do
not reconstruct a history, integrate radial constraints, advance time, or
establish a numerical mass balance. The caller supplies one common metric.
"""
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_characteristic_background import _array, MAX_POINTS
from .nsc_atlas_cap_radial_reference import reference_potential
from .nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from .nsc_atlas_polar_cap import _scalar
from .nsc_atlas_reduced_cap import _potential_reaction_quotient


def _frame_gap(f,alpha):
    """u-f*du/df without subtracting its equal leading terms.

    For |q|<=.1 the first omitted real-arithmetic term is bounded by
    |f|*q^12/518918400/(1-q^2/180). Floating function error is separate.
    """
    q=np.sqrt(alpha/6)*f
    if abs(q)<=.1:
        z=q*q
        return -f*z*(1/3+z*(1/30+z*(1/840+z*(1/45360+z/3991680))))
    return f*(np.sinh(q)/q-np.cosh(q))


def _metric(r,f,mass,lapse,alpha,budget):
    """Pointwise constrained coefficients with analytic vacuum cancellation."""
    A0=1-r*r
    A=A0-2*alpha*mass/r
    mu=np.exp(lapse)
    if np.any(A<=0) or np.any(mu<=0):
        raise ValueError('characteristics require positive A and lapse')
    P=reference_potential(f,alpha,budget=budget)['Pdiff']
    reaction=_potential_reaction_quotient(f,alpha)
    B=A*mu;beta=B/A0
    beta_z=alpha*mu*(2*mass/(r*r)-r*P-4*mass/A0)
    source=alpha*mu*(r*reaction+(2*mass/(r*r)-r*P)*f)
    return dict(A0=A0,A=A,mu=mu,B=B,beta=beta,beta_z=beta_z,
                source=source,Pdiff=P)


def characteristic_point_rhs(radius,chi,E_plus,E_minus,chi_z,E_plus_z,E_minus_z,
                             *,mass,lapse,alpha,clock_rate,coordinate_rate=0.,budget):
    """Evaluate up to 513 smooth points strictly inside 0<r<1.

    chi=r*f, E_plus=A0*(p+q), E_minus=A0*(p-q), p=r*Pi, q=chi_r.
    Spatial derivatives use z=atanh(r). Returned rates are T_x*partial_T
    plus coordinate_rate*partial_z, where coordinate_rate is the motion
    dz/dx of the query point. Its metric mass is vacuum-subtracted N.
    chi_z is independently supplied, so primitive inconsistency is reported.
    Axis limits, discontinuities and moving interfaces belong to the caller.
    """
    budget=_budget(budget);before=budget.evaluations
    al=_scalar(alpha,'alpha',nonnegative=True)
    tx=_scalar(clock_rate,'clock rate',positive=True)
    names=('r','chi','Ep','Em','chiz','Epz','Emz','mass','lapse','zx')
    values=(radius,chi,E_plus,E_minus,chi_z,E_plus_z,E_minus_z,mass,lapse,coordinate_rate)
    arrays=[_array(v,k) for k,v in zip(names,values)]
    try: data=dict(zip(names,np.broadcast_arrays(*arrays)))
    except ValueError as error: raise ValueError('characteristic point shapes differ') from error
    r=data['r'];count=r.size
    if count>MAX_POINTS or np.any(r<=0) or np.any(r>=1):
        raise ValueError('smooth characteristic points require 0<r<1')
    budget.check(80*count)
    try:
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            f=data['chi']/r
            metric=_metric(r,f,data['mass'],data['lapse'],al,budget)
            A0,A,B,beta,bz,S=(metric[k] for k in ('A0','A','B','beta','beta_z','source'))
            ep,em,epz,emz=(data[k] for k in ('Ep','Em','Epz','Emz'))
            plus_T=beta*epz+bz*ep-A0*S
            minus_T=-beta*emz-bz*em-A0*S
            chi_T=beta*(ep+em)/2
            p=(ep+em)/(2*A0);q=(ep-em)/(2*A0);h=q-f
            rho=r*r*metric['Pdiff']/2+A*(p*p+h*h)/4
            target=A*B*p*h/2
            result=dict(chi_rate=tx*chi_T+data['zx']*data['chiz'],
                        E_plus_rate=tx*plus_T+data['zx']*epz,
                        E_minus_rate=tx*minus_T+data['zx']*emz,
                        primitive_residual=data['chiz']-(ep-em)/2,
                        f=f,p=p,q=q,h=h,mass_density=rho,
                        temporal_mass_target=target,**metric)
    except FloatingPointError as error:
        raise ValueError('nonfinite smooth characteristic arithmetic') from error
    if not all(np.all(np.isfinite(v)) for v in result.values()):
        raise ValueError('nonfinite smooth characteristic result')
    budget.check()
    return dict(result,evaluations=budget.evaluations-before,
                scope='smooth point PDE and conditional mass targets; no radial solve or evolution',
                numerical_error_enclosure=False)


def characteristic_wall_port(radius,chi,E_minus,*,mass,lapse,alpha,exterior_wall,
                             C_wall,eta_plus,eta_plus_x,budget):
    """Derive the incoming trace and common wall rates from the existing hybrid.

    chi and E_minus are the cap's actual endpoint traces. The radial mass
    and lapse must come from that same cap. exterior_wall receives immutable
    {r,c} and returns r_V,c_V,log_C_V on the common invariant field value.
    The existing boost gauge determines the evolving exterior lapse rate.
    No independent incident waveform or wall-value evolution is supplied.
    The scalar tangent residual is a local conversion check, not an
    independently evolved regional matching residual.
    """
    budget=_budget(budget);before=budget.evaluations
    R=_scalar(radius,'wall radius',positive=True)
    if R>=1: raise ValueError('wall radius must be below one')
    chi=_scalar(chi,'wall primitive');em=_scalar(E_minus,'outgoing characteristic')
    n=_scalar(mass,'wall mass');delta=_scalar(lapse,'wall lapse')
    al=_scalar(alpha,'alpha',nonnegative=True)
    C=_scalar(C_wall,'wall C',positive=True)
    boost=_scalar(eta_plus,'wall boost');boost_x=_scalar(eta_plus_x,'wall boost derivative')
    if not callable(exterior_wall): raise ValueError('exterior wall must be callable')
    budget.check(256)
    try:
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            f=chi/R
            metric=_metric(np.asarray(R),np.asarray(f),np.asarray(n),np.asarray(delta),al,budget)
            A0,A,B=(float(metric[k]) for k in ('A0','A','B'))
            af=np.sqrt(al/6)*f
            u=f*(np.sinh(af)/af if af else 1.)
            rootD=np.cosh(af);c=R*u
            if not np.isfinite(c): raise ValueError('nonfinite common wall scalar')
            common=MappingProxyType(dict(r=R,c=float(c)))
            raw=exterior_wall(common)
            names={'r_V','c_V','log_C_V'}
            if not isinstance(raw,Mapping) or set(raw)!=names:
                raise ValueError('exterior callback requires three outgoing first jets')
            outgoing={k:_scalar(raw[k],'exterior '+k) for k in names}
            state=dict(common,eta_plus=boost)
            cap=dict(A=A,delta=delta,w_minus=(em/A0+f)/R)
            first=hybrid_wall_boundary(state,outgoing,cap,dict(C_wall=C,log_C_s=0.),alpha=al)
            xs=float(first['Vdot_plus'])
            Lx=2*boost_x-2*float(first['k_plus'])/xs+2*outgoing['log_C_V']
            hybrid=hybrid_wall_boundary(state,outgoing,cap,dict(C_wall=C,log_C_s=Lx*xs),alpha=al)
            tx=float(hybrid['T_s'])/xs;rx=float(hybrid['r_s'])/xs
            ep=A0*(R*float(hybrid['w_plus'])+f)
            p=(ep+em)/(2*A0);q=(ep-em)/(2*A0);h=q-f
            chi_x=tx*B*p+rx*q
            gap=_frame_gap(f,al)
            c_x=rootD*chi_x+gap*rx
            target=float(hybrid['c_s'])/xs
            rho=R*R*float(metric['Pdiff'])/2+A*(p*p+h*h)/4
            temporal=A*B*p*h/2
            rates=dict(r=rx,L=Lx,T=tx,s=1/xs,c=target,chi=chi_x)
            scalars=dict(E_plus=ep,E_minus=em,p=p,q=q,h=h,cap_c_x=c_x,frame_gap=gap,
                scalar_tangent_residual=c_x-target,mass_density=rho,
                temporal_mass_target=temporal,moving_mass_target=tx*temporal+rx*rho)
    except FloatingPointError as error:
        raise ValueError('nonfinite characteristic wall arithmetic') from error
    if not all(np.isfinite(v) for v in (*rates.values(),*scalars.values())):
        raise ValueError('nonfinite characteristic wall result')
    budget.check()
    return dict(common=dict(common),outgoing=outgoing,hybrid=hybrid,wall_rates=rates,
                **scalars,evaluations=budget.evaluations-before,
                scope='local common-value port; conditional mass targets, no regional evolution')
