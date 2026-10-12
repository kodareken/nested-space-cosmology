"""Exact analytic radius bounds for the current two-function history.

Chebyshev derivative extrema retain the represented basis. Expanding a high
degree profile into powers and then summing absolute coefficients can destroy
the radius margin even when the represented function is small.
"""
from fractions import Fraction as Q
import math

from .nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily


def chebyshev_profile_bounds(profile):
    """Global sup bounds of the analytic profile and its first two derivatives.

    For M>=1, |T_n^(j)(x)| <= T_n^(j)(M) on [-M,M], j=0,1,2.
    On [-1,1] these are the Chebyshev endpoint derivative bounds; outside
    that interval the derivative polynomials increase in absolute value.
    The actual binary domain-map coefficients and profile coefficients are
    treated as exact rationals. Evaluation roundoff is a separate error.
    """
    if not isinstance(profile, LocalAxialFunction):
        raise TypeError("owned LocalAxialFunction required")
    offset, scale = map(Q, profile._derivatives[0].mapparms())
    center, inner, outer = map(Q, (profile.center, profile.inner, profile.outer))
    width = outer - inner
    if width <= 0:
        raise ValueError("positive axial transition width required")
    M = max(Q(1), abs(offset+scale*(center-outer)), abs(offset+scale*(center+outer)))
    values = [(Q(1), Q(0), Q(0))]
    if len(profile.coefficients) > 1:
        values.append((M, Q(1), Q(0)))
    for _ in range(2, len(profile.coefficients)):
        f, d, dd = values[-1]
        old_f, old_d, old_dd = values[-2]
        values.append((2*M*f-old_f, 2*f+2*M*d-old_d, 4*d+2*M*dd-old_dd))
    polynomial = tuple(abs(scale)**j * sum(
        (abs(Q(coefficient))*row[j] for coefficient, row in zip(profile.coefficients, values)),
        Q()) for j in range(3))
    if any(value < 0 for value in polynomial):
        raise ArithmeticError("Chebyshev derivative upper became negative")
    p0, p1, p2 = polynomial
    # The same owned smooth cutoff bounds as axial_profile_bounds.
    bounds = (p0, p1 + 8*p0/width, p2 + 16*p1/width + 105*p0/width**2)
    return {"map_endpoint_upper": M, "polynomial_bounds": polynomial,
            "profile_bounds": bounds}


def radius_bounds(family, rho_up):
    """Reuse the proved background slab and bound this history's radius."""
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("owned local incoming family required")
    rho_up = Q(rho_up)
    rho_upper, axial_lower, reference_lower = Q(33, 32), Q(4, 5), Q(7, 5)
    if not Q(1) < rho_up <= rho_upper:
        raise ValueError("upstream lies outside the existing background slab")
    a2_lower = 3*((1+rho_upper**2)*(Q(157, 200)-Q(1, 65))-rho_upper)-1
    if a2_lower <= axial_lower**2 or reference_lower**2 >= 2:
        raise ArithmeticError("existing background lower bound failed")
    sigma = Q(family.normal_outer)
    if not Q(0) < sigma <= Q(float(0.03)):
        raise ValueError("normal support lies outside the declared history class")
    w, u = (chebyshev_profile_bounds(profile) for profile in family.functions)
    W, U = w["profile_bounds"], u["profile_bounds"]
    delta = tuple(sigma*W[j] + sigma**3*U[j]/6 for j in range(3))
    lower = reference_lower-delta[0]
    if lower <= 0:
        raise ValueError("analytic bounds do not establish a positive radius")
    return {
        "rho_upper": rho_upper, "rho_up": rho_up, "axial_lower": axial_lower,
        "reference_radius_lower": reference_lower, "normal_support": sigma,
        "w": w, "U": u, "delta_radius_bounds": delta, "radius_lower": lower,
        "profile_evaluation_roundoff_included": False,
        "physical_local_gate": "OPEN",
    }


def value_integral_bounds(bounds, mass, angular):
    """K0,D,K1 for existing fixed-history value energy majorants.

    No history-tangent majorant is established here. Both w_z and U_z enter
    the spatial derivative of the actual-radius potential.
    """
    duration = bounds["rho_up"]-1
    a, r = bounds["axial_lower"], bounds["radius_lower"]
    ell, mass = abs(Q(angular)), abs(Q(mass))
    delta_z = bounds["delta_radius_bounds"][1]
    return {
        "K0": duration/a*(mass+ell/r),
        "D": duration/a**2,
        "K1": duration*ell*delta_z/(a*r**2),
        "history_tangent_bound": None,
    }


def rational_record(value):
    """Serialize exact values and outward binary upper bounds for inspection."""
    if isinstance(value, Q):
        rounded = float(value)
        if Q(rounded) < value:
            rounded = math.nextafter(rounded, math.inf)
        return {"exact_rational": str(value), "binary64_upper": rounded}
    if isinstance(value, dict):
        return {key: rational_record(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [rational_record(item) for item in value]
    return value
