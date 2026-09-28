"""A massless reference specification and its finite-cutoff conversion.

The reference is canonical Hadamard, not a new independent gravitational
source. The full compact determinant is reused from its authenticated record.
"""
from math import atan, ceil, exp, isfinite, pi

from .nsc_spectral_endpoint import light_endpoint


def conformal_ratio(u):
    """Analytic continuation of the existing future conformal spatial factor."""
    return 3*(pi-atan(u))+(3*u-u*u)/(1+u*u)


def reference_extension(u, width=.1):
    """Actual metric for u>=-width/2, ultrastatic spatial metric for u<=-width.

    Lapse is converted to conformal proper time. In the constant region
    that time extends freely, irrespective of the auxiliary u chart range.
    Width is a reference choice; it is not selected by the physical action.
    """
    if not isfinite(u) or not isfinite(width) or not 0<width<=.25:
        raise ValueError("finite u and reference width in (0,.25] required")
    if u>=-width/2:
        return conformal_ratio(u)
    if u<=-width:
        return 3*pi
    t=2*(u+width)/width
    first,second=exp(-1/t),exp(-1/(1-t))
    switch=first/(first+second)
    return 3*pi+switch*(conformal_ratio(u)-3*pi)


def canonical_reference_density(radius):
    """One massless four-component Dirac, in the declared +--- anomaly scheme.

    Exact on the invariant de Sitter reference, and the asymptotic value
    of the conformally regular untwisted reference class described in the note.
    It is not the finite-radius KS stress or the full NSC vacuum source.
    """
    if not isfinite(radius) or radius<=0:
        raise ValueError("positive curvature radius required")
    return 11/(960*pi*pi*radius**4)


def light_conversion(radius, matching_cutoff):
    """Metric projection of C_nu,mu = Gamma_light,nu - Gamma_light,ren.

    The Euler logarithm is independent of the S4 radius under variation at
    fixed nu,mu. The finite remainder is retained, including at nu*a<1.
    """
    if not isfinite(matching_cutoff) or matching_cutoff<=0:
        raise ValueError("positive matching cutoff required")
    x=radius*matching_cutoff
    maximum=max(2,ceil(10*x)+2)
    heat=light_endpoint(radius,matching_cutoff,maximum)["d_log_radius"]
    volume=8*pi*pi*radius**4/3
    power=(2/3)*(x**4-x*x)
    anomaly=11/90
    remainder=heat-power-anomaly
    # A decreasing m^3 exp(-m²/x²) envelope bounds omitted sphere levels.
    tail=(2/3)*x*x*(maximum*maximum+x*x)*exp(-(maximum/x)**2)
    return {"matching_cutoff":matching_cutoff,"nu_times_radius":x,"last_level":maximum,
            "finite_light_log_derivative":heat,"power_log_derivative":power,
            "renormalized_log_derivative":anomaly,"finite_remainder_log_derivative":remainder,
            "conversion_log_derivative":heat-anomaly,
            "finite_light_density":heat/(4*volume),"power_density":power/(4*volume),
            "canonical_reference_density":canonical_reference_density(radius),
            "finite_remainder_density":remainder/(4*volume),
            "conversion_density":(heat-anomaly)/(4*volume),
            "sphere_log_derivative_tail_bound":tail}


def matched_reference_budget(full_density, radius, matching_cutoff):
    """Move the same contribution between reference and complement exactly once."""
    light=light_conversion(radius,matching_cutoff)
    old_complement=full_density-light["finite_light_density"]
    converted=old_complement+light["conversion_density"]
    total=light["canonical_reference_density"]+converted
    return {"light":light,"old_complement_density":old_complement,
            "converted_complement_density":converted,"unchanged_full_density":full_density,
            "matched_total_density":total,"matching_residual":total-full_density,
            "unmatched_replacement_density":light["canonical_reference_density"]+old_complement}
