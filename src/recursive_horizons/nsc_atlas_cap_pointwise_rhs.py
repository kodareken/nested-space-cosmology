"""Smooth cap scalar equations on supplied local field and constrained metric jets.

One point away from the centre, in the original H*ell=1 cap. No state,
metric, material junction, incoming history or evolution is solved here.
The defect of a weak RHS is not a source to be added to the physical model.
"""
from __future__ import annotations

from math import exp,isfinite


def scalar_rhs(*,radius,chi,k,chi_z,chi_zz,k_z,mass,lapse,alpha,
               potential_difference,potential_derivative):
    """Return chi_T, k_T and partial_z(chi_T) at fixed r=tanh(z).

    P= (Vhat-3)/alpha and U=dP/df use the existing continuous alpha=0
    definitions. p=r*Pi=k/A0; h=r*Phi=chi_z/A0-chi/r; B=A*exp(delta).
    Mass is vacuum-subtracted: m_full=r^3/2+alpha*mass.
    """
    vals=tuple(map(float,(radius,chi,k,chi_z,chi_zz,k_z,mass,lapse,alpha,
                          potential_difference,potential_derivative)))
    if not all(isfinite(v) for v in vals) or not 0<radius<1 or alpha<0:
        raise ValueError('finite cap data with 0<r<1 and alpha>=0 required')
    r=radius;A0=1-r*r;A=A0-2*alpha*mass/r
    if A<=0: raise ValueError('pointwise cap metric must remain positive')
    mu=exp(lapse);B=A*mu
    p=k/A0;h=chi_z/A0-chi/r
    h_z=chi_zz/A0+2*r*chi_z/A0-chi_z/r+chi*A0/(r*r)
    p_z=(k_z+2*r*k)/A0
    rho=r*r*potential_difference/2+A*(p*p+h*h)/4
    A_r=-2*r-2*alpha*(rho/r-mass/(r*r))
    delta_r=alpha*(p*p+h*h)/(2*r)
    B_z=A0*mu*(A_r+A*delta_r)
    result=dict(chi_T=B*p,k_T=A0*B*h/r+B_z*h+B*h_z-A0*r*mu*potential_derivative,
        chi_T_z=B_z*p+B*p_z,p=p,h=h,h_z=h_z,p_z=p_z,A0=A0,A=A,B=B,
        normalized_mass_density=rho,A_r=A_r,delta_r=delta_r,B_z=B_z)
    if not all(isfinite(v) for v in result.values()): raise ValueError('nonfinite pointwise cap RHS')
    return result


def incoming_trace_defect(*,radius,continuum,represented_chi_T,
                          represented_k_T,represented_chi_T_z):
    """Difference of cap contributions to d_T[(p+h)/r] at fixed radius.

    This excludes the exterior incoming derivative and wall-position terms.
    It cannot be identified with the complete coupled transmission rate.
    """
    r=float(radius)
    if not 0<r<1: raise ValueError('incoming trace requires 0<r<1')
    A0=1-r*r
    v=float(represented_chi_T)-continuum['chi_T']
    d=float(represented_k_T)-continuum['k_T']
    vz=float(represented_chi_T_z)-continuum['chi_T_z']
    density=d/(r*A0);slope=vz/(r*A0);value=-v/(r*r)
    result=dict(chi_equation_defect=v,k_equation_defect=d,chi_equation_defect_z=vz,
        incoming_density_channel=density,incoming_velocity_slope_channel=slope,
        incoming_velocity_value_channel=value,incoming_cap_rate_defect=density+slope+value)
    if not all(isfinite(v) for v in result.values()): raise ValueError('nonfinite represented cap RHS')
    return result
