"""Relative-state bounds for the fixed free Dirac tower in the NSC child.

The homogeneous mode covariance is compared on the SAME geometry/domain.
Finite weighted trace norms are required; no absolute vacuum is selected.
"""
from math import atan, isfinite, pi, sqrt


PARALLEL_RATIO_MINIMUM = sqrt(3*pi/2-1)


def interior_factors(rho):
    """Existing interior, with s=a_parallel/r and angularly integrated volume.

    Natural units and the stored unit throat profile. The volume factor
    multiplies one coordinate length dz, not a compact spatial volume.
    """
    if not isfinite(rho) or rho > 0:
        raise ValueError("finite expanding-child coordinate rho<=0 required")
    radius=sqrt(1+rho*rho)
    ratio2=3*(pi/2-atan(rho))-(1+3*rho)/(1+rho*rho)
    ratio=sqrt(ratio2)
    return {"rho":rho,"radius":radius,"parallel_over_radius":ratio,
            "a_parallel":ratio*radius,"proper_volume_per_dz":4*pi*ratio*radius**3}


def relative_energy_bound(radius, mass_moment, momentum_moment):
    """|Delta rho| <= [M+K/r]/[4*pi*s_min*r^3].

    M=sum/integral |m_j| ||Delta C_j||_1 and
    K=sum/integral (|k_j|/s_min+|kappa_j|) ||Delta C_j||_1.
    The mode measure includes physical degeneracies exactly once.
    This is a homogeneous angular-average density or per-dz density bound,
    not a pointwise bound on arbitrary inhomogeneous Hadamard states.
    """
    if not all(isfinite(x) for x in (radius,mass_moment,momentum_moment)):
        raise ValueError("finite radius and weighted state moments required")
    if radius < 1 or mass_moment < 0 or momentum_moment < 0:
        raise ValueError("expanding branch r>=1 and nonnegative moments required")
    volume=4*pi*PARALLEL_RATIO_MINIMUM*radius**3
    return {"massive_part":mass_moment/volume,
            "redshifting_part":momentum_moment/(volume*radius),
            "total":(mass_moment+momentum_moment/radius)/volume}


def channel_energy_norm(radius, parallel_scale, mass, momentum, angular):
    """Norm of m*beta+(k/a_parallel)*alpha_z+(kappa/r)*alpha_sphere."""
    if radius<=0 or parallel_scale<=0:
        raise ValueError("positive metric scales required")
    return sqrt(mass**2+(momentum/parallel_scale)**2+(angular/radius)**2)


def relative_null_state_source(rho, delta_tu, delta_tv):
    """Existing free LLL state constants mapped into the child proper frame.

    Flux is the covariant orthonormal T_hatT_hatz component. The common
    geometric anomaly cancels in the difference; no absolute vacuum chosen.
    """
    g=interior_factors(rho)
    factor=1/(4*pi*g["a_parallel"]**2*g["radius"]**2)
    return {"rho_coordinate":rho,"radius":g["radius"],
            "delta_density":(delta_tu+delta_tv)*factor,
            "delta_p_parallel":(delta_tu+delta_tv)*factor,
            "delta_p_sphere":0.,"delta_T_hatT_hatz":(delta_tv-delta_tu)*factor,
            "absolute_density_bound":(abs(delta_tu)+abs(delta_tv))/(4*pi*PARALLEL_RATIO_MINIMUM**2*g["radius"]**4)}
