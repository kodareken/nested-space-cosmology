"""Apply the existing compact weight to the published four-sphere spectrum.

The Euclidean sphere is an invariant endpoint control, not the global
Lorentzian black-universe state. No Einstein or vacuum term is added to
the full Dirac functional a second time.
"""
from math import ceil, exp, pi, sqrt

import numpy as np
from scipy.special import exp1

from .nsc_lorentzian import geometry


def geometry_scales(rho):
    """Existing curvature contractions in the quasi-global metric chart.

    The sectional expressions extend through A=0 as scalar-contraction
    inputs; they are not a static orthonormal observer inside the horizon.
    All curvature values scale as L^-2 under a physical profile dilation.
    """
    x = float(rho)
    beta, beta_prime, lapse = map(float, geometry(x))
    radius = sqrt(1+x*x)
    prime = -2*beta*beta_prime
    second = 6*(np.arctan(x)-pi/2)+6*x/(1+x*x)
    radial_prime, radial_second = x/radius, 1/radius**3
    a = second/2
    b = prime*radial_prime/(2*radius)
    c = lapse*radial_second/radius+b
    d = (1-lapse*radial_prime**2)/radius**2
    scalar = 2*(a+2*b+2*c-d)
    riemann2 = 4*(a*a+2*b*b+2*c*c+d*d)
    return {"rho": x, "A": lapse, "radius": radius,
            "sections": [a,b,c,d], "R_L": scalar,
            "Riemann_squared": riemann2,
            "curvature_scale_squared": sqrt(riemann2/24),
            "largest_section_magnitude": max(map(abs,(a,b,c,d)))}


def sphere_spectrum(radius, last_level):
    """D has signs +/-m/a; each sign has 2(m^3-m)/3 states, m>=2."""
    if not np.isfinite(radius) or radius <= 0 or last_level < 2:
        raise ValueError("positive sphere radius and at least level two required")
    levels = np.arange(2, last_level+1, dtype=float)
    return levels**2/radius**2, 4*(levels**3-levels)/3


def endpoint(weight, radius, last_level=None):
    """Gamma=1/2 Tr4 h(D4^2) and its log-radius derivative at fixed inputs.

    Angular tail bounds use the established nonpositive-warp min-max
    comparison. Retained compact Galerkin errors are assessed separately.
    """
    cutoff_radius = weight.cutoff*radius
    maximum = max(2,ceil(8*cutoff_radius)+2) if last_level is None else last_level
    y, multiplicity = sphere_spectrum(radius, maximum)
    responses = [weight.response(z) for z in y]
    action = .5*sum(d*r["weight"] for d,r in zip(multiplicity,responses))
    slope = -sum(d*z*r["d_squared"] for d,z,r in zip(multiplicity,y,responses))
    # E1(z)<=exp(-z)/z, compact heat trace <=1+ell*Lambda/sqrt(pi).
    # |y*h'(y)|<=exp(s*a)*Tr_compact exp(-epsilon/Lambda^2), by
    # Cauchy-Schwarz in the same generalized quadratic form.
    if maximum <= sqrt(1.5)*cutoff_radius:
        raise ValueError("sphere tail must start past both decreasing envelopes")
    heat = 1+weight.interval*weight.cutoff/sqrt(pi)
    tail_exp = exp(-(maximum/cutoff_radius)**2)
    action_bound = heat*cutoff_radius**4*tail_exp/3
    slope_bound = (2/3)*exp(weight.path*weight.amplitude)*heat*cutoff_radius**2*(maximum**2+cutoff_radius**2)*tail_exp
    volume = 8*pi*pi*radius**4/3
    return {"radius": radius, "last_level": maximum, "action": float(action),
            "d_log_radius": float(slope), "four_volume": volume,
            "isotropic_action_response": float(slope/(4*volume)),
            "sphere_action_tail_bound": action_bound,
            "sphere_log_derivative_tail_bound": slope_bound}


def light_endpoint(radius, matching_cutoff, last_level):
    y, multiplicity = sphere_spectrum(radius,last_level)
    return {"action": float(.5*np.dot(multiplicity,exp1(y/matching_cutoff**2))),
            "d_log_radius": float(np.dot(multiplicity,np.exp(-y/matching_cutoff**2)))}


def local_endpoint(radius, coefficients, last_level=80):
    """Exact light field plus its already-matched local complement through a4.

    On S4, C^2=box R=0, integral E4=64*pi^2. The Dirac R^2 coefficient
    is zero; no claim about the full theory's finite R^2 term is made.
    """
    light = light_endpoint(radius,coefficients["matching_cutoff"],last_level)
    volume = 8*pi*pi*radius**4/3
    V,A,E = (coefficients[k] for k in ("V_Dirac","A_Dirac","C_Euler_Dirac"))
    return {"action": light["action"]+volume*(V-12*A/radius**2)+64*pi*pi*E,
            "d_log_radius": light["d_log_radius"]+volume*(4*V-24*A/radius**2),
            "light_action": light["action"], "light_d_log_radius": light["d_log_radius"]}
