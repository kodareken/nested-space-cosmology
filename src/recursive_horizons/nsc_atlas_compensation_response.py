"""Local spatial response of the retained mass-compensated initial geometry.

The transform uses the areal geometric mass measure, not proper energy.
Its separate envelope measure retains the omitted positive gradient mass.
No initial field, reservoir radius or inherited parameter is retuned.
"""
from functools import cached_property
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre, spherical_jn

from .nsc_atlas_scale_hierarchy import ScaleHierarchy
from .nsc_atlas_core_reservoir import compact_bump


def sinc_minus_one(x):
    """Cancellation-safe j0(x)-1, including its value at zero."""
    x=np.asarray(x,dtype=float);y=x*x
    series=y*(-1/6+y*(1/120+y*(-1/5040+y/362880)))
    return np.where(abs(x)<.05,series,np.sinc(x/np.pi)-1)


def nodes(order,left,right):
    x,w=roots_legendre(order)
    return left+(x+1)*(right-left)/2,w*(right-left)/2


class CompensationResponse:
    def __init__(self,model=None,*,order=160,rtol=2e-11):
        self.model=model or ScaleHierarchy(1e9,.001,.5)
        if order<24 or order>1024:raise ValueError('quadrature order outside admitted range')
        self.order=order;self.rtol=rtol
        self.preparation=self.model.solve(rtol=rtol)
        self.core,self.core_solution=self.model.core_data(rtol=rtol)
        self.mass=self.core['mass'];self.outer=self.preparation['reservoir_outer_radius']
        self.inner=self.model.radius
        _,self.reservoir_solution=self.model.reservoir(self.outer,self.mass,rtol=rtol)

    def cumulative(self,s):
        s=np.asarray(s,dtype=float);r=s*self.outer
        w=self.model.frequency;z=self.model.weak_amplitude
        x=np.clip(w*z*r,0,20)
        core=np.array([1.,1/w,1/w**2])@self.core_solution.sol(x)
        core=z*core/(w*self.mass)
        a=np.clip((r-self.inner)/(self.outer-self.inner),0,1)
        reservoir=self.reservoir_solution.sol(a)[0]
        return np.where(r<=self.inner,core,np.where(r<=self.outer,reservoir,0.))

    def densities(self,s):
        """Return (dM_GR/ds, dM_env/ds, dM_grad/ds)/Mc."""
        s=np.asarray(s,dtype=float);r=s*self.outer;m=self.model
        w=m.frequency;d=m.ambient_coefficient;k=m.reference.scalar.kappa
        f,fp=m.core.field(np.minimum(r,self.inner))
        in_core=r<=self.inner;f=np.where(in_core,f,0.);fp=np.where(in_core,fp,0.)
        bump,bp=compact_bump(r,self.inner,self.outer)
        delta=f-m.depletion*d/w*bump
        derivative=fp-m.depletion*d/w*bp
        B=1-m.curvature*r*r-np.divide(2*self.mass*self.cumulative(s),r,
                                    out=np.zeros_like(r),where=r>0)
        # U-U0 at the coherent turning point equals the leading envelope
        # mass density. Its cancellation-safe expression is exact here.
        du=w*d*delta+w*w*delta*delta/2
        env=k*r*r*du*self.outer/(2*self.mass)
        grad=k*r*r*B*derivative*derivative*self.outer/(4*self.mass)
        inside=(s>=0)&(s<=1)
        return tuple(np.where(inside,v,0.) for v in (env+grad,env,grad))

    @cached_property
    def quadrature(self):
        cut=self.inner/self.outer
        parts=[nodes(self.order,a,b) for a,b in ((0,.75*cut),(.75*cut,cut),(cut,1))]
        s=np.concatenate([p[0] for p in parts]);weights=np.concatenate([p[1] for p in parts])
        gr,env,grad=self.densities(s)
        return s,weights,gr,env,grad

    def transform(self,K,*,channel='geometry',core_only=False):
        """K=k Rd. Never silently subtract the envelope's total residual."""
        K=np.asarray(K,dtype=float)
        if np.any(~np.isfinite(K)) or np.any(K<0):raise ValueError('K must be finite and nonnegative')
        s,q,gr,env,grad=self.quadrature
        values={'geometry':gr,'envelope':env,'gradient':grad}[channel]
        weights=q*values
        if core_only:weights=weights*(s<self.inner/self.outer)
        # Separately retain zeroth moment; the exact constraint makes it
        # zero for GR, but not for the unmodified leading envelope.
        return np.sum(sinc_minus_one(K[...,None]*s)*weights,axis=-1)+np.sum(weights)

    def by_parts(self,K):
        # Integrating j0(Ks) dM gives j0(K) M(1) + K int M j1(Ks) ds.
        s,q,*_=self.quadrature;K=np.asarray(K,dtype=float)
        endpoint=self.reservoir_solution.y[0,-1]
        return np.sinc(K/np.pi)*endpoint+K*np.sum(
            q*self.cumulative(s)*spherical_jn(1,K[...,None]*s),axis=-1)

    def half_power(self,*,channel='geometry'):
        # First crossing in the declared finite bracket. The producer
        # checks monotonicity on this interval, not on all wavenumbers.
        def f(K):return float(self.transform(K,channel=channel)/self.transform(K,channel=channel,core_only=True)-1/math.sqrt(2))
        return brentq(f,.2,8.,xtol=2e-12)

    def projected(self,b,*,channel='geometry'):
        """Projected areal mass, surface contrast and mean-minus-local.

        b is projected radius/Rd. Surface units are Mc/Rd^2. This is a
        spatial projection, not a full relativistic null observable.
        Line-of-sight integration removes the Abel endpoint singularity.
        """
        if not math.isfinite(b) or b<=0:raise ValueError('projected radius must be positive')
        if b>=1:
            total=0. if channel=='geometry' else float(self.transform(0.,channel=channel))
            return dict(cylindrical_mass=total,surface=0.,mean_minus_local=total/(np.pi*b*b))
        index={'geometry':0,'envelope':1,'gradient':2}[channel]
        cuts=[b]+[r for r in (.75*self.inner/self.outer,self.inner/self.outer) if r>b]+[1.]
        surface=0.;outside=0.
        for left,right in zip(cuts[:-1],cuts[1:]):
            z,q=nodes(self.order,math.sqrt(max(0,left*left-b*b)),math.sqrt(right*right-b*b))
            s=np.sqrt(b*b+z*z);density=self.densities(s)[index]
            surface+=float(np.sum(q*density/(s*s)))/(2*np.pi)
            outside+=float(np.sum(q*density*z*z/(s*s)))
        # Direct GR zero is exact in the mathematical boundary condition;
        # retain the measured zeroth moment for all numerical channels.
        total=float(self.transform(0.,channel=channel))
        cylinder=total-outside
        return dict(cylindrical_mass=cylinder,surface=surface,
                    mean_minus_local=cylinder/(np.pi*b*b)-surface)

    def record(self):
        s,q,gr,env,grad=self.quadrature
        moments=[float(np.sum(q*gr*s**p)) for p in (0,2,4)]
        a2=-moments[1]/6;remainder=float(np.sum(q*abs(gr)*s**4)/120)
        K=self.half_power();Ke=self.half_power(channel='envelope')
        sample=np.geomspace(.0001,30,151)
        direct=self.transform(sample);parts=self.by_parts(sample)
        gap=float(np.sum(q*grad))
        ratio=np.array([self.transform(k)/self.transform(k,core_only=True) for k in np.linspace(.2,K,81)])
        return dict(frequency=self.model.frequency,weak_amplitude=self.model.weak_amplitude,
            depletion=self.model.depletion,quadrature_order=self.order,constraint_rtol=self.rtol,
            core_mass=self.mass,core_radius=self.inner,reservoir_radius=self.outer,
            mass_moments_normalized=moments,low_K_quadratic_coefficient=a2,
            low_K_fourth_remainder_bound=remainder,
            gradient_mass_over_core=gap,envelope_total_mass_over_core=float(np.sum(q*env)),
            half_power_K=K,envelope_half_power_K=Ke,
            half_power_wavenumber_atlas=K/self.outer,
            Lambda_over_core_mass_k_half_cubed=self.model.vacuum*self.outer**3/(self.mass*K**3),
            root_bracket_monotone=bool(np.min(np.diff(ratio))>0),
            integration_by_parts_max_gap=float(np.max(abs(direct-parts))),
            all_K_envelope_error_bound=gap,
            sampled_envelope_difference=float(np.max(abs(direct-self.transform(sample,channel='envelope')))),
            low_K_envelope_residual_crossover=math.sqrt(gap/a2),
            initial_parent_curvature_Rd_squared=self.model.curvature*self.outer**2,
            minimum_core_B=self.core['minimum_core_beta_squared'])


def point_core_limit(depletion=.5,order=160):
    """Derived omega->infinity shape limit at fixed zeta and depletion."""
    s,q=nodes(order,0,1);b,_=compact_bump(s,0.,1.)
    weight=q*s*s*(2*depletion*b-depletion**2*b*b)
    J=float(np.sum(weight));weight/=J
    def F(K):return -float(np.sum(weight*sinc_minus_one(K*s)))
    K=brentq(lambda x:F(x)-1/math.sqrt(2),.2,8.,xtol=2e-13)
    Q=1.5/J
    return dict(depletion=depletion,J=J,half_power_K=K,
        low_K_quadratic_coefficient=float(np.sum(weight*s*s)/6),
        Lambda_Rd_cubed_over_mass=Q,Lambda_over_core_mass_k_half_cubed=Q/K**3)


def retained_evolution_response(radius,initial_u,final_u,scale_factor,ambient,scaled_initial_Rd):
    """Postprocess an authenticated retained spherical evolution, no steps.

    K is physical k times the ORIGINAL Rd; the profile need not remain
    compact at that radius. The finite domain's monopole is kept explicitly.
    """
    from .nsc_atlas_hierarchy_evolution import enclosed
    r=radius;dr=float(r[1]-r[0]);A=scale_factor
    first=r*r*(abs(initial_u/r)**2+2*ambient*np.real(initial_u/r))
    initial_core=enclosed(r,first,dr,20.)
    radial=r*r*(abs(final_u/r)**2+2*ambient*np.real(final_u/r))
    weights=dr*radial/initial_core;s=A*r/scaled_initial_Rd
    def full(K):return float(np.sum(weights*sinc_minus_one(K*s))+np.sum(weights))
    def core(K):return enclosed(r,radial*np.sinc(K*s/np.pi),dr,20/A)/initial_core
    K=brentq(lambda k:full(k)/core(k)-1/math.sqrt(2),.2,8.,xtol=1e-11)
    core_ratio=core(0.)
    moments=[float(np.sum(weights*s**p)) for p in (0,2,4)]
    return dict(final_scale_factor=A,core_mass_over_initial=core_ratio,
        total_excess_over_initial_core=moments[0],moments_in_original_Rd_units=moments,
        low_K_quadratic_coefficient=-moments[1]/6,half_power_K_in_original_Rd_units=K,
        envelope_mass_inside_original_Rd_over_initial_core=enclosed(r,radial,dr,scaled_initial_Rd/A)/initial_core,
        finite_outer_radius_over_initial_Rd=A*(r[-1]+dr)/scaled_initial_Rd,
        sampling_scope='Physical wavenumber, finite retained domain and final event. Monopole retained; no exact k^4 limit or compact support at original Rd is imposed.')
