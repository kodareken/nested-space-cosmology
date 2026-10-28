"""Conditional Planck distance diagnostic; no invented NSC acoustic physics.

A supplied decoupling redshift and profiled sound ratio produce an optimistic
geometric diagnostic. The unchanged parent has no predicted acoustic horizon.
"""
import math

import numpy as np
from scipy.integrate import quad

from .nsc_atlas_inherited_matter import exact_initial_geometry,uniform_bound
from .nsc_atlas_shared_ruler import fractions,C_OVER_100


# Chen--Huang--Wang 1808.05724v1, Table 1. Use the same curvature-relaxed
# compression in BOTH columns. This reconstructs its rounded covariance.
PRIOR_MEAN=np.array([1.7429,301.409])
PRIOR_SIGMA=np.array([.0051,.091])
PRIOR_COVARIANCE=np.diag(PRIOR_SIGMA)@np.array([[1.,.54],[.54,1.]])@np.diag(PRIOR_SIGMA)


def early_parent_bound(observer_scale,z_star,omega=10**9):
    """Bootstrap on A_min<=A<=1 using H>=.24 A^(-3/2).

    No evolution grid and no assertion that the scalar produces recombination.
    """
    d=exact_initial_geometry();C=float(d['curvature']);rv=float(d['vacuum']);E0=float(d['excitation'])
    Ao=float(observer_scale);Amin=Ao/(1+z_star);h=6/25
    if not 1<=Ao<=20 or not 0<Amin<1:raise ValueError('declared observer/redshift range required')
    # Minimum of E0/6-C*A+rv*A^3/3 is E0/6-2*C*sqrt(C/rv)/3.
    assert E0/6-2*C*math.sqrt(C/rv)/3>h*h
    Hmax=math.ceil(math.sqrt(rv/3+2*E0/(3*Amin**3)))+1
    delta=3*Hmax/(2*omega)
    variation=2*C/h*(Amin**-.5-1)+4*E0/(3*h)*(Amin**-1.5-1)
    gamma=3*variation/(2*omega*(1-delta))
    lower=(1-gamma)/(1+delta);upper=1/((1-gamma)*(1-delta))
    if not .5<lower<1<upper<2:raise ValueError('early scalar bootstrap failed')
    error=max(1-lower,upper-1)
    r=(E0/3)/(E0/3-C)*error
    # A point-specific two-way initial-cap footprint.
    past_angle=2*math.sqrt(C)/h*(1-math.sqrt(Amin))
    observer_angle=math.sqrt(C)/.266*(1-1/Ao)
    margin=2.83-observer_angle-2*past_angle
    if margin<=0:raise ValueError('CMB probe leaves the proved initial-cap domain')
    return dict(minimum_scale=Amin,source_density_relative_error=error,
                hubble_squared_relative_error=r,past_angle_upper=past_angle,
                observer_angle_upper=observer_angle,causal_margin_lower=margin,
                H_variation_upper=variation,H_upper=Hmax)


def dust_transverse_distance(model,theta,z_star):
    """H0 D_M/c on the dust/vacuum/curvature continuation.

    x=1/sqrt(1+z) removes the high-z integration boundary layer.
    """
    f=fractions(model,theta);m,k,v=f['Omega_m'],f['Omega_k'],f['Omega_Lambda']
    if m<=0 or z_star<=0:raise ValueError('positive matter and source redshift required')
    lo=1/math.sqrt(1+z_star)
    integrand=lambda x:2/math.sqrt(m+k*x*x+v*x**6)
    I,error=quad(integrand,lo,1,epsabs=2e-12,epsrel=2e-12)
    if k<0:
        angle=math.sqrt(-k)*I
        if angle>=math.pi/2:raise ValueError('CMB probe outside monotone sine branch')
        DM=math.sin(angle)/math.sqrt(-k)
    else:
        angle=0.;DM=I
    x,w=np.polynomial.legendre.leggauss(96)
    xx=(1+lo)/2+(1-lo)*x/2
    control=(1-lo)*np.sum(w/np.sqrt(m+k*xx**2+v*xx**6))
    if abs(I-control)>2e-10:raise ValueError('independent high-z quadrature failed')
    return DM,dict(conformal_angle=angle,quadrature_difference=abs(I-control),quad_error_indicator=error)


def profiled_distance_prior(R):
    """Profile unknown l_A, rather than silently fixing baryon/sound physics."""
    shift=R-PRIOR_MEAN[0]
    ell=PRIOR_MEAN[1]+PRIOR_COVARIANCE[1,0]/PRIOR_COVARIANCE[0,0]*shift
    residual=np.array([R,ell])-PRIOR_MEAN
    score=float(residual@np.linalg.solve(PRIOR_COVARIANCE,residual))
    assert abs(score-(shift/PRIOR_SIGMA[0])**2)<1e-8*(1+score)
    return dict(R=R,profiled_lA=ell,conditional_chi_squared_lower_bound=score,
                profiled_acoustic_parameters=1,remaining_geometric_contrasts=1)


def diagnostic(model,theta,rd_h,z_star):
    f=fractions(model,theta)
    DM,numerical=dust_transverse_distance(model,theta,z_star)
    R=math.sqrt(f['Omega_m'])*DM
    result=dict(**profiled_distance_prior(R),H0_DM_over_c=DM,
                implied_rs_star_over_rd=math.pi*(C_OVER_100/rd_h)*DM/profiled_distance_prior(R)['profiled_lA'],
                numerical=numerical)
    if model=='parent':
        b=early_parent_bound(math.exp(theta),z_star)
        future=uniform_bound(10**9);r=b['hubble_squared_relative_error']
        Rlo=R*math.sqrt(float(future['lower'])/(1+r))
        Rhi=R*math.sqrt(float(future['upper'])/(1-r))
        closest=min(max(PRIOR_MEAN[0],Rlo),Rhi)
        result.update(scalar_past_control=b,unaveraged_R_enclosure=[Rlo,Rhi],
                      phase_optimistic_chi_squared_lower_bound=((closest-PRIOR_MEAN[0])/PRIOR_SIGMA[0])**2)
    return result
