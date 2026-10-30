"""Coordinate-independent observer comparisons for the atlas map.

Exact Kasner-de Sitter is a same-vacuum anisotropic background with chosen
initial shear, not a selected cosmological history. Circular geodesics and
wave redshift are defined controls, not fitted Galactic trajectories.
"""
import math
import numpy as np
from scipy.special import logsumexp
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_golden_vortex import PHI,kasner


def log_scale(u,exponents):
    if u<=0:raise ValueError('positive u=3Ht required')
    log_sinh=u-math.log(2)+math.log(-math.expm1(-2*u))
    log_tanh=math.log(-math.expm1(-u))-math.log1p(math.exp(-u))
    return log_sinh/3+(np.asarray(exponents)-1/3)*log_tanh


def redshift(emit,receive,direction,exponents):
    n=np.asarray(direction,dtype=float);n=n/np.linalg.norm(n)
    keep=n!=0;growth=log_scale(receive,exponents)-log_scale(emit,exponents)
    log_ratio=.5*logsumexp(np.log(n[keep]**2)+2*growth[keep])
    return math.expm1(float(log_ratio))


def circular_orbit(radius_over_mass=10.):
    r=radius_over_mass
    if r<=3:raise ValueError('timelike Schwarzschild circular orbit requires r/m>3')
    L2=r*r/(r-3);E2=(r-2)**2/(r*(r-3))
    return dict(radius_over_mass=r,L_squared_over_mass_squared=L2,E_squared=E2,
        radial_velocity=0.,angular_frequency_times_mass=r**-1.5,
        radial_potential_second_derivative=2*(r-6)/(r**3*(r-3)),
        stable=r>6,spacetime_Kretschmann_times_mass4=48/r**6,
        external_balancing_masses=0,
        scope='Exact test-particle geodesic in the isolated Schwarzschild metric; excludes self-force, radiation and other bodies, so it is not a perpetual-motion claim for the real Galaxy.')


def calculation():
    minimum=VolumeParent().vacuum_roots()[1];p=kasner(PHI);H=math.sqrt(minimum['vacuum']/3)
    rows=[]
    for initial,final in ((.1,.2),(1.,2.),(4.,5.),(8.,9.)):
        axis=[redshift(initial,final,n,p) for n in np.eye(3)]
        rows.append(dict(u_emit=initial,u_receive=final,axis_redshifts=axis,
            directional_ratio_max_over_min=max(1+z for z in axis)/min(1+z for z in axis),
            late_shear_fraction=1/math.cosh(final)**2))
    return dict(schema='nsc-atlas-observer-geometry-v1',
        invariant='1+z=(-k_mu u_emit^mu)/(-k_mu u_receive^mu)',
        anisotropic=dict(Lambda=minimum['vacuum'],H=H,Kasner_exponents=p.tolist(),
            metric='ds^2=-dt^2+sum_i a_i(t)^2 dx_i^2',
            solution='a_i=A_i sinh(u)^(1/3) tanh(u/2)^(p_i-1/3), u=3Ht',
            rates='H_i/H=coth(u)+3(p_i-1/3)csch(u)',
            shear='sigma^2=3H^2 csch(u)^2; sigma^2/(3 Hmean^2)=sech(u)^2',
            optical_law='1+z(n)=sqrt(sum_i n_i^2 [a_i(receive)/a_i(emit)]^2)',
            rows=rows,
            shape_conclusion='Directional expansion isotropizes; this does not select a spherical boundary or finite topology.',
            state='The golden Kasner exponents are a declared comparison inherited from the earlier BKL motif, not selected by the late-vacuum action.',
            past='The generic Kasner-de Sitter past still has a singularity; no birth/reset or origin-of-time proof.',
            observer_epoch_selected=False),
        orbital=circular_orbit(),
        wave=dict(model='Massless geometric-optics probe on isotropic FLRW; local orthonormal speed c.',
            relation='f proportional to 1/a, wavelength proportional to a, f*wavelength=c; dr/dt=c/a for a radial comoving coordinate.',
            twofold_expansion=dict(frequency_ratio=.5,wavelength_ratio=2.,local_speed_ratio=1.,radiation_energy_density_ratio=1/16),
            scope='Radiation-density ratio assumes conserved photon number in comoving volume and no absorption/emission. A change in the photon operator or medium can change propagation and must be specified separately.'),
        galactic_centre=dict(source='https://arxiv.org/abs/1807.09409',
            observation='S2 near 2018 pericentre: combined gravitational and transverse-Doppler shift about 200 km/s divided by c, in addition to line-of-sight motion.',
            model_fit_performed=False),
        physical_cosmology_fit=False)
