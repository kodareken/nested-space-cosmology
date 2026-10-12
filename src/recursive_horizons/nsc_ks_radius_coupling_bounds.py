"""Continuous A-to-D coupling bounds for the owned pure-radius history.

These bounds close M and M_z inputs only. They neither bound the numerical
reference residual nor the physical upstream preparation error.
"""
from copy import deepcopy
from fractions import Fraction as Q

from .nsc_ks_current_history_bounds import radius_bounds, rational_record
from .nsc_ks_residual_error import nonnegative


def radius_coupling_integrals(w_bounds, u_bounds, normal_support,
                              reference_radius_lower, radius_lower, angular):
    """Integrate induced row norms continuously over the upstream half slab.

    M = -(i ell/a) delta_r/[r0 (r0+delta_r)] sigma_2,
    M_z = -(i ell/a) delta_r,z/(r0+delta_r)^2 sigma_2.
    The reference radius is homogeneous in z and sigma_2 has row-sum norm 1.
    |d rho| = a |ds| cancels 1/a. The support lies in -sigma <= s <= 0;
    integrating |s| W_j + |s|^3 U_j/6 yields sigma^2 W_j/2+sigma^4 U_j/24.
    Integrating the full support is valid also when preparation begins later.
    """
    if len(w_bounds) < 2 or len(u_bounds) < 2:
        raise ValueError("profile value and axial derivative bounds required")
    w0, w1 = (nonnegative(v, "w bound") for v in w_bounds[:2])
    u0, u1 = (nonnegative(v, "U bound") for v in u_bounds[:2])
    s, r0, r = (nonnegative(v, name) for v, name in (
        (normal_support, "normal support"),
        (reference_radius_lower, "reference radius lower"),
        (radius_lower, "radius lower")))
    if isinstance(angular, bool):
        raise ValueError("finite angular label required")
    ell = nonnegative(abs(Q(angular)), "absolute angular")
    if s == 0 or r0 == 0 or r == 0:
        raise ValueError("positive support and radius lower bounds required")
    integral0 = s**2*w0/2 + s**4*u0/24
    integral1 = s**2*w1/2 + s**4*u1/24
    return {"M_integral": ell*integral0/(r0*r),
            "Mz_integral": ell*integral1/r**2}


def history_coupling_bounds(family, rho_up, angular):
    bounds = radius_bounds(family, rho_up)
    return radius_coupling_integrals(
        bounds["w"]["profile_bounds"], bounds["U"]["profile_bounds"],
        bounds["normal_support"], bounds["reference_radius_lower"],
        bounds["radius_lower"], angular)


def with_radius_coupling_bounds(inputs, family, rho_up, angular, *, bits=160):
    """Add exactly two owned slots to the earlier cone input interface.

    Recompute the relevant geometry and require its binding to agree. Do not
    rewrite archived cone records or fill other unknown inputs with zero.
    """
    from .nsc_ks_ball_trajectory import exact_upper
    from .nsc_ks_current_field_cone import PROPAGATION_INPUT_NAMES
    from flint import arb, ctx

    bounds = radius_bounds(family, rho_up)
    expected = {
        "rho_up": bounds["rho_up"], "radius_lower": bounds["radius_lower"],
        "normal_support": bounds["normal_support"],
        "w_profile_bounds": bounds["w"]["profile_bounds"],
        "U_profile_bounds": bounds["U"]["profile_bounds"],
    }
    for name, value in expected.items():
        if inputs.get("geometry", {}).get(name) != rational_record(value):
            raise ValueError("coupling geometry binding mismatch: " + name)
    result = deepcopy(inputs)
    # The angular channel is bound independently of the radius profile.
    from .nsc_ks_ball_trajectory import restored_upper
    channel = result["matter"]["absolute_angular"]["value"]
    with ctx.workprec(bits):
        ell = abs(Q(angular))
        exact_ell = arb(ell.numerator)/arb(ell.denominator)
        if not restored_upper(channel).contains(exact_ell):
            raise ValueError("coupling angular channel mismatch")
        for name, value in history_coupling_bounds(family, rho_up, angular).items():
            result["propagation"][name] = {
                "value": exact_upper(arb(value.numerator)/arb(value.denominator)),
                "owner": "nsc_ks_radius_coupling_bounds.history_coupling_bounds",
                "reason": None,
            }
    result["complete_for_propagate_difference_error"] = all(
        result["propagation"][name]["value"] is not None
        for name in PROPAGATION_INPUT_NAMES)
    return result
