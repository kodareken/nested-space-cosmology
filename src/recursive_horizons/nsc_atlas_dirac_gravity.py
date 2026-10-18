"""Spherical Einstein/clock/wall variations for the actual quiet Dirac probe.

ds2=-C dU dV+r² dOmega²; L=log C, U=T-x,V=T+x.
The supplied nonvacuum stress excludes the regional constant Lambda term.
gravity_coupling is 8*pi*G in the caller's consistent units. For a Dirac
stress expressed in hbar=c=ell=1, restore 8*pi*(l_P/ell)^2 explicitly.

This is an equation interface, not an evolution, a selected state or a slip.
Linear equations about a vacuum region require a conserved complete source.
A metric variation of a nonzero quantum mean has additional Ward/tadpole
terms; a matter-only commutator cannot silently replace that source.
Controls: independent 4D Einstein contraction, constraint identity, actual
quiet linearization and moving-wall/gauge checks, 60 process-CPU seconds.
No physical field evolution, new P001 grid or reference JSON is admitted.
"""
from __future__ import annotations

import numpy as np

from .nsc_atlas_cap_wall_response import (
    H_ELL, K_ELL, M_OVER_ELL, reference_coefficients)


def einstein_from_null_jets(*, r, C, ru, rv, ruu, ruv, rvv, Lu, Lv, Luv):
    """Covariant UU, UV, VV and orthonormal angular Einstein components.

    Pure algebra also accepts symbolic jets for an independent calculation.
    The metric domain has r>0,C>0. No vacuum term is included in G.
    """
    return dict(uu=-2*(ruu-Lu*ru)/r,
        uv=2*ruv/r+2*ru*rv/r**2+C/(2*r**2),
        vv=-2*(rvv-Lv*rv)/r,
        pt=-(2*Luv+4*ruv/r)/C)


def null_stress_from_orthonormal(*, C, rho, radial_pressure, outward_energy_flux,
                                angular_pressure):
    """Transform the existing physical Dirac stress; no new degeneracy factor."""
    return dict(uu=C*(rho+radial_pressure+2*outward_energy_flux)/4,
        uv=C*(rho-radial_pressure)/4,
        vv=C*(rho+radial_pressure-2*outward_energy_flux)/4,
        pt=angular_pressure)


def quiet_null_jets(fields):
    """Exact vacuum jets from the existing actual radar field interface."""
    r,C = fields['radius'],fields['conformal_factor']
    ru,rv = (fields['radius_t']-fields['radius_x'])/2,(
        fields['radius_t']+fields['radius_x'])/2
    Lu,Lv = fields['log_C_u'],fields['log_C_v']
    return dict(r=r,C=C,ru=ru,rv=rv,Lu=Lu,Lv=Lv,
        ruu=Lu*ru,rvv=Lv*rv,ruv=fields['radius_uv'],Luv=fields['log_C_uv'],
        Lambda=np.where(fields['is_cap'],3*H_ELL**2,0.))


def linearized_null_rhs(background, *, s, su, sv, q, qu, qv, stress,
                        gravity_coupling=1.):
    """Linear mixed evolution and the two null constraints.

    s=delta r, q=delta log C. The outputs are delta r_UU, delta r_UV,
    delta r_VV and delta L_UV. Residual null-coordinate freedom and
    constrained initial data are the caller's responsibility.
    """
    r,C,ru,rv,Lu,Lv,lam = (background[k] for k in
        ('r','C','ru','rv','Lu','Lv','Lambda'))
    k = gravity_coupling
    mixed = rv*su+ru*sv
    return dict(
        ruu=Lu*su+ru*qu-k*r*stress['uu']/2,
        rvv=Lv*sv+rv*qv-k*r*stress['vv']/2,
        ruv=-mixed/r+(ru*rv/r**2+C/(4*r**2)+lam*C/4)*s
            +C*(-1/(4*r)+lam*r/4)*q+k*r*stress['uv']/2,
        Luv=2*mixed/r**2-(4*ru*rv+C)*s/r**3+C*q/(2*r**2)
            -k*stress['uv']-k*C*stress['pt']/2)


def constraint_transport(*, r, ru, rv, constraint_u, constraint_v,
                         conservation_u, conservation_v, gravity_coupling=1.):
    """Exact compatibility identity on mixed Einstein evolution.

    E_U=r_UU-L_U*r_U+k*r*t_UU/2, and
    Q_U=t_UU,V+t_UV,U+2*r_V*t_UU/r
        +(2*r_U/r-L_U)*t_UV+C*r_U*p_t/r.
    Q_V swaps U,V. Q=0 is stress conservation. Thus r*E_U is constant
    along V, and r*E_V along U, on the constraint surface.
    """
    return dict(constraint_u_v=-rv*constraint_u/r
                +gravity_coupling*r*conservation_u/2,
        constraint_v_u=-ru*constraint_v/r
                +gravity_coupling*r*conservation_v/2)


def quiet_wall_jets(radius):
    """Common-normal wall jets, ordered cap then exterior; ell=1."""
    R = float(radius)
    g = reference_coefficients(R)
    beta = np.array([g['beta_minus'],g['beta_plus']],dtype=float)
    beta_R = np.array([g['beta_minus_x'],g['beta_plus_x']],dtype=float)
    fp = np.array([-2*H_ELL**2*R,2*M_OVER_ELL/R**2])
    fpp = np.array([-2*H_ELL**2,-4*M_OVER_ELL/R**3])
    return dict(radius=R,v=float(g['v']),acceleration=float(g['radial_acceleration']),
        beta=beta,beta_R=beta_R,beta_RR=float(g['beta_second']),
        Lx=2*beta_R,Lxx=fpp,rxx=float(g['radial_acceleration'])+fp,
        Ktime=beta_R,Kangle=beta/R,k_wall=K_ELL)


def wall_metric_variation(radius, *, s, sx, q, qx, displacement,
                          velocity, acceleration, normal_stress_jump=0.,
                          gravity_coupling=1.):
    """Linear geometry at x=w(T), in the same conformal-coordinate gauge.

    Physical continuity requires the two surface radii Z and log-lapses Q
    below to agree. Their jump and both Israel jumps must vanish.
    The force residual is a Gauss/Israel compatibility condition, not an
    independent extra junction datum. Baseline vacuum normal pressure has
    constant jump; normal_stress_jump is only its supplied perturbation.
    """
    wall = quiet_wall_jets(radius)
    arrays = [np.asarray(x,dtype=float) for x in (s,sx,q,qx)]
    if any(x.shape != (2,) or not np.isfinite(x).all() for x in arrays):
        raise ValueError('two finite real cap/exterior traces required')
    if not all(np.isfinite(x) for x in (displacement,velocity,acceleration,
                                       normal_stress_jump,gravity_coupling)):
        raise ValueError('finite displacement, stress and coupling required')
    s,sx,q,qx = arrays
    R,beta,Lx = wall['radius'],wall['beta'],wall['Lx']
    w,wd,wdd = displacement,velocity,acceleration
    Z,Q = s+beta*w,q+Lx*w
    time = wdd+qx/2+w*wall['Lxx']/2-Q*Lx/4
    angle = (sx+w*wall['rxx']+wd*wall['v']-Q*beta/2)/R-beta*Z/R**2
    return dict(surface_radius=Z,surface_log_lapse=Q,
        proper_clock_rate=Q/2,Ktime=time,Kangle=angle,
        radius_continuity=Z[1]-Z[0],metric_continuity=Q[1]-Q[0],
        Israel_time=time[1]-time[0],Israel_angle=angle[1]-angle[0],
        normal_force=np.mean(time+2*angle)+gravity_coupling*normal_stress_jump/(2*K_ELL))


def radar_observer_variation(background, *, s, q, wall_u, wall_v):
    """Read a supplied metric perturbation at fixed proper-wall radar labels.

    wall_u/v supply the wall time_shift, displacement, clock_rate and velocity
    at the two null intersections. If the input gauge has a moving wall,
    Uhat=U+A(U), Vhat=V+B(V), A=delta_tau+w, B=delta_tau-w.
    These shifts and their Jacobians must be included before reading a clock.
    The radial light protocol emits phase -Uhat from the unit proper wall
    clock and receives it on a fixed radar-position observer.
    """
    A = wall_u['time_shift']+wall_u['displacement']
    B = wall_v['time_shift']-wall_v['displacement']
    Au = wall_u['clock_rate']+wall_u['velocity']
    Bv = wall_v['clock_rate']-wall_v['velocity']
    radar_q = q-A*background['Lu']-B*background['Lv']-Au-Bv
    radar_s = s-A*background['ru']-B*background['rv']
    fractional_frequency = -radar_q/2
    return dict(radius=radar_s,log_conformal_factor=radar_q,
        received_over_wall_frequency_fraction=fractional_frequency,
        received_over_wall_frequency_change=fractional_frequency/np.sqrt(background['C']))
