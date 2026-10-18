"""Local spherical curvature readout in a specified radial observer frame.

These are Riemann sectional contrasts, not a replacement definition for
the nonlocal Newtonian potentials Phi and Psi. The caller supplies the
observer's orthonormal stress and mass excess above Lambda*r^3/6.
"""
from __future__ import annotations


def tidal_contrasts(*,radius,mass_excess,energy_density,radial_pressure,
                    angular_pressure,gravity_coupling=1.):
    """A=R_0r0r-R_0theta0theta, B=R_rtheta rtheta-R_thetaphi thetaphi.

    Convention: vacuum Schwarzschild R_0r0r=-2M/r^3-Lambda/3.
    Lambda cancels from both contrasts. Weyl w is invariant under radial
    boosts; A and B generally are not when matter has anisotropic stress.
    """
    k=gravity_coupling
    w=mass_excess/radius**3-k*(energy_density-radial_pressure+angular_pressure)/6
    anisotropy=k*(radial_pressure-angular_pressure)
    return dict(weyl=w,temporal=-3*w-anisotropy/2,spatial=-3*w+anisotropy/2)


def mean_tidal_response(*,radius,background_mass,mass_change,radius_change,
                        stress,gravity_coupling=1.):
    """Finite prepared stress, first order at fixed physical radar labels.

    Background is vacuum plus fixed Lambda and a constant central mass.
    The areal-radius variation in M/r^3 is retained. A change of the radial
    observer frame acting on that vacuum adds no anisotropic Ricci term.
    Return the leading change of B/A only when background A is nonzero.
    """
    k=gravity_coupling
    rho,pr,pt=(stress[x] for x in ('energy_density','radial_pressure','angular_pressure'))
    w=(mass_change/radius**3-3*background_mass*radius_change/radius**4
       -k*(rho-pr+pt)/6)
    anisotropy=k*(pr-pt)
    a0=-3*background_mass/radius**3
    return dict(weyl_change=w,temporal_change=-3*w-anisotropy/2,
        spatial_change=-3*w+anisotropy/2,background_contrast=a0,
        ratio_change=None if background_mass==0 else anisotropy/a0)
