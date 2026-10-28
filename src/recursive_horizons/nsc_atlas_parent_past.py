"""Finite backward Cauchy continuation of the same inherited scalar parent.

The bootstrap extends source control to A=1/3; a two-way causal bound keeps
observer rays in the domain determined by the original homogeneous cap.
"""
from fractions import Fraction as F
import math

from .nsc_atlas_inherited_matter import exact_initial_geometry, uniform_bound


def past_bound(omega=10**9):
    w=F(str(omega))
    d=exact_initial_geometry()
    Amin,hmin,hmax,Kmax=F(1,3),F(1,3),F(3),F(57,2500)
    C,rv,E0=d['curvature'],d['vacuum'],d['excitation']
    # Bootstrap 1/2 < A^3 E/E0 < 2. The minimum of H_lower^2
    # lies at A=1 since 3*(E0/6)>2*C throughout A<=1.
    assert 3*(E0/6)>2*C
    assert hmin*hmin<rv/3+E0/6-C
    assert hmax*hmax>rv/3+2*E0/(3*Amin**3)
    delta=3*hmax/(2*w)
    variation=Kmax/(2*hmin)*(Amin**-2-1)+2*E0/(3*hmin)*(Amin**-3-1)
    if w<=0 or delta>=1:
        raise ValueError('positive corrected energy required')
    gamma=3*variation/(2*w*(1-delta))
    if gamma>=1:
        raise ValueError('rational backward comparison unavailable')
    lo=(1-gamma)/(1+delta)
    hi=1/((1-gamma)*(1-delta))
    if not F(1,2)<lo<1<hi<2:
        raise ValueError('bootstrap does not close for this frequency')
    future=uniform_bound(omega)
    lower=min(lo,future['lower']); upper=max(hi,future['upper'])
    error=max(1-lower,upper-1)
    # For A>=1/3 this denominator is increasing.
    assert rv*Amin*Amin-C>0
    ratio=(E0/3)/(rv*Amin**3/3+E0/3-C*Amin)
    r=ratio*error
    # Bound pi-asin(R_outer sqrt(C)) from below and the future
    # conformal horizon from above without evaluating transcendental endpoints.
    assert F(19,10)**2*C<F(3,10)**2
    assert F(31,100)-F(31,100)**3/6>F(3,10)
    initial_angle_lower=F(314,100)-F(31,100)
    assert 4*C<rv/3   # arcsin(sqrt(C/(rv/3))) < pi/6 < 11/21
    future_angle_upper=F(11,21)
    assert C<F(19,125)**2
    past_angle_upper=F(19,125)*(1/Amin-1)/hmin
    two_way_angle_upper=future_angle_upper+2*past_angle_upper
    assert two_way_angle_upper<initial_angle_lower
    return dict(omega=w,minimum_scale=Amin,lower=lower,upper=upper,
                density_relative_error=error,hubble_squared_relative_error=r,
                delta=delta,Gamma=gamma,H_variation_upper=variation,
                H_past_lower=hmin,H_past_upper=hmax,
                initial_angle_lower=initial_angle_lower,
                future_angle_upper=future_angle_upper,past_angle_upper=past_angle_upper,
                two_way_angle_upper=two_way_angle_upper,
                causal_margin_lower=initial_angle_lower-two_way_angle_upper)


def solve_past(*,clock='proper',omega=200):
    """Independent finite numerical control, not the high-frequency proof."""
    import time
    import numpy as np
    from scipy.integrate import solve_ivp
    from .nsc_atlas_flux_conversion import ScalarFlux
    from .nsc_atlas_shared_scalar_slice import SharedScalarSlice
    started=time.process_time()
    atlas=SharedScalarSlice(scalar=ScalarFlux(omega=omega))
    d=atlas.parent_causal_domain()
    H0=d['initial_H']; a0=d['initial_scale_factor']; rv=d['ground_density']
    chi0=atlas.field-atlas.scalar.minimum(atlas.parent_flux)
    E0=omega**2*chi0**2/2
    target=-math.log(3); calls=0
    def count():
        nonlocal calls
        calls+=1
        if calls>80000 or time.process_time()-started>60:
            raise RuntimeError('backward scalar control exhausted finite allowance')
    if clock=='proper':
        def rhs(t,y):
            count()
            cc,pp,N,H=y
            return atlas.parent_rates(cc,pp,N,H)[:4]
        def event(t,y):return y[2]-target
        event.terminal=True
        sol=solve_ivp(rhs,(0.,-4.),(chi0,0.,0.,H0),events=event,
                      dense_output=True,method='DOP853',rtol=2e-12,atol=2e-14)
        if not sol.success or not len(sol.t_events[0]):
            raise RuntimeError('backward proper-time endpoint not reached')
        final=float(sol.t_events[0][0]);cc,pp,N,H=sol.sol(np.linspace(0,final,1001))
    elif clock=='efolds':
        def rhs(N,y):
            count()
            cc,pp=y
            H=math.sqrt((rv+(pp*pp+omega**2*cc*cc)/2)/3-math.exp(-2*N)/a0**2)
            return [pp/H,(-3*H*pp-omega**2*cc)/H]
        sol=solve_ivp(rhs,(0.,target),(chi0,0.),dense_output=True,
                      method='DOP853',rtol=2e-12,atol=2e-14)
        if not sol.success:raise RuntimeError(sol.message)
        N=np.linspace(0,target,1001);cc,pp=sol.sol(N)
        H=np.sqrt((rv+(pp*pp+omega**2*cc*cc)/2)/3-np.exp(-2*N)/a0**2)
    else:
        raise ValueError('proper or efolds clock required')
    R=np.exp(3*N)*(pp*pp+omega**2*cc*cc)/(2*E0)
    constraint=H*H-((rv+(pp*pp+omega**2*cc*cc)/2)/3-np.exp(-2*N)/a0**2)
    return dict(clock=clock,omega=omega,minimum_scale=1/3,
                sampled_density_ratio_range=[float(R.min()),float(R.max())],
                endpoint_density_ratio=float(R[-1]),maximum_Friedmann_residual=float(np.max(abs(constraint))),
                RHS_evaluations=calls,process_CPU_seconds=time.process_time()-started)
