"""Exact initial shear/internal-coupling relation for the inherited kicks.

The shear is the full contraction norm sqrt(sigma_ij sigma^ij), without
the alternative factor1/2 convention. H_star is the reference vacuum
curvature scale; the initial trace expansion itself is zero.
"""
import math

import numpy as np

from .nsc_atlas_finite_parent_data import source_profile
from .nsc_atlas_kick_york import solve_york_kick
from .nsc_atlas_parent_circulation import Parent


def _source(value):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or value not in (1,2):
        raise ValueError('source_harmonic must be integer1 or2')
    return int(value)


def geometry_factor_squared(z,source_harmonic):
    """Return the closed York geometry factor squared for a matched source.

    With t=q^2=(z+8)/(6(z+3)), the two rational expressions are
      F1^2=(3t^2-45t+211)/(450(1+t)^2),
      F2^2=(75t^4-1215t^3+7813t^2-564t+12)/(22050(1-t^2)^2).
    They follow from the ambient polynomial York solution at
    X=(q,sqrt((1-q^2)/2),0,sqrt((1-q^2)/2)), not an observational fit.
    """
    source=_source(source_harmonic)
    if isinstance(z,(bool,np.bool_)):raise ValueError('positive finite z=mu^2/H_star^2 required')
    try:z=float(z)
    except (TypeError,ValueError) as error:raise ValueError('positive finite z=mu^2/H_star^2 required') from error
    if not math.isfinite(z) or z<=0:raise ValueError('positive finite z=mu^2/H_star^2 required')
    # Equivalent ratio, without overflowing six times a large finite z.
    t=(1+5/(z+3))/6
    if source==1:return (3*t*t-45*t+211)/(450*(1+t)**2)
    return (75*t**4-1215*t**3+7813*t*t-564*t+12)/(22050*(1-t*t)**2)


def observer_relation(*,epsilon=.25,source_harmonic=1,parent=None):
    """Compare the closed ratio with an independent reference York solve.

    For any nonzero kick amplitude lambda and positive conformal factor,
      ||sigma||=|lambda| psi^-6 ||Abar||,
      partial_tau^3 ln(g_KK)=-kappa mu^2 chi lambda f psi^-6.
    Thus the displayed dimensionless ratio is F/z. The current protocol
    fixes the freely specified conformal TT datum to zero. Changing it can
    change the result. No SI scale or quantum state is selected.
    """
    source=_source(source_harmonic)
    if isinstance(epsilon,(bool,np.bool_)):
        raise ValueError('the ratio requires finite0<epsilon<=1')
    try:epsilon=float(epsilon)
    except (TypeError,ValueError) as error:raise ValueError('the ratio requires finite0<epsilon<=1') from error
    if not math.isfinite(epsilon) or not 0<epsilon<=1:
        raise ValueError('the ratio requires finite0<epsilon<=1; epsilon0 gives0/0')
    p=parent or Parent.retained()
    if not math.isfinite(p.kappa) or p.kappa<=0 or not math.isfinite(p.mass_squared) or p.mass_squared<=0:
        raise ValueError('positive finite gravitational coupling and mass required')
    H=p.H;z=p.mass_squared/H**2
    mode=p.initial();a=math.sqrt(2)/math.pi*mode[2];b=math.sqrt(12)/math.pi*mode[4]
    q=a/b;t=q*q;matched=(1+5/(z+3))/6
    if not 0<q<1 or not math.isclose(t,matched,rel_tol=2e-13,abs_tol=2e-15):
        raise ValueError('the source is outside the matched two-harmonic family')
    v=math.sqrt((1-t)/2);point=np.array([[q,v,0.,v]])
    kick=solve_york_kick(epsilon=epsilon,source_harmonic=source,parent=p)
    tensor=kick.values(point)[0];norm=float(np.linalg.norm(tensor))
    chi=epsilon*float(source_profile(point,p)[0][0]);f=float(kick.momentum_profile(point)[0])
    if chi*f==0:raise ValueError('the ratio is undefined at a node of chi*f')
    factor_york=norm/(p.kappa*abs(chi*f))
    factor_squared=geometry_factor_squared(z,source);factor=math.sqrt(factor_squared)
    jerk_coefficient=-p.kappa*p.mass_squared*chi*f
    ratio_york=H*H*norm/abs(jerk_coefficient);ratio=factor/z
    gap=abs(ratio_york-ratio)
    return dict(schema='nsc-atlas-kick-relation-v1',source_harmonic=source,epsilon=epsilon,
        z=z,amplitude_ratio_a_over_b=q,amplitude_ratio_squared=t,
        equal_energy_ratio_squared=matched,observer_quaternion=point[0].tolist(),
        reference_H_star=H,reference_mu_squared=p.mass_squared,kappa=p.kappa,
        closedfactor=factor,geometry_factor_squared=factor_squared,
        geometry_factor_from_york=factor_york,
        shear_to_logcoupling_jerk_ratio=ratio,
        independent_york_shear_to_logcoupling_jerk_ratio=ratio_york,
        ratio_absolute_check_error=gap,ratio_relative_check_error=gap/abs(ratio),
        chi_at_observer=chi,momentum_profile_at_observer=f,
        shear_norm_over_abs_lambda_psi_minus6=norm,
        logcoupling_jerk_over_lambda_psi_minus6=jerk_coefficient,
        shear_norm_definition='sqrt(sigma_ij sigma^ij), without a factor1/2',
        jerk_identity='partial_tau^3 ln(g_KK)=-kappa mu^2 chi Pi_chi; K_s=d(d+2)/(2kappa), mu^2=omega0^2 exp(-d s), ln g_KK=constant-(d+2)s/2',
        hypotheses=['Initially constant s and normal(s)=0','Initial trK=0',
            'Observer initially follows the unit normal','Matched two-harmonic chi profile',
            'Nonzero lambda and chi*f','Zero freely specified conformal TT datum'],
        cancellation='The factors psi^-6, |lambda|, epsilon and overall kick-profile normalization cancel. The source shape and matched amplitude ratio remain specified.',
        TT_scope='An added conformal transverse-traceless tensor can change the shear and this ratio; no radiation-state selection is claimed.',
        reference_scope='H_star is the vacuum reference curvature scale, not initial expansion (which is zero).',
        SI_scale_selected=False,quantum_state_selected=False,
        scope='Exact initial scalar relation conditional on the Hamiltonian constraint and stated conformal-TT protocol. No late-time universal ratio or physical-state selection is asserted.')
