"""Normalized mixed Einstein/scalar equations at strictly positive radius.

The fixed cap has h=1 and the exterior h=0. The invariant c=chi/epsilon
survives alpha=0. This pointwise operator solves no PDE, chooses no gauge,
and is not a regular-centre discretization or a constraint certificate.
"""
from collections.abc import Mapping

import numpy as np


def _fields(value, names, label):
    if not isinstance(value, Mapping) or set(value) != set(names):
        raise ValueError(label + " must contain exactly " + ", ".join(names))
    return value


def _arrays(**values):
    arrays = []
    for name, value in values.items():
        try:
            if isinstance(value, (list, tuple)) and any(
                isinstance(item, (bool, np.bool_))
                for item in np.asarray(value, dtype=object).flat
            ):
                raise ValueError
            raw = np.asarray(value)
            if raw.dtype.kind not in "iuf" or raw.size == 0:
                raise ValueError
            array = np.asarray(raw, dtype=float)
            if not np.all(np.isfinite(array)):
                raise ValueError
        except (ValueError, TypeError, OverflowError) as error:
            raise ValueError(name + " must contain finite real numbers, not booleans") from error
        arrays.append(array)
    try:
        return dict(zip(values, np.broadcast_arrays(*arrays)))
    except ValueError as error:
        raise ValueError("input shapes are not broadcast compatible") from error


def bulk_operator(state, first_jets, *, alpha, region):
    """Return mixed evolution and explicitly qualified constraint data.

    ``state`` contains r, C, c; ``first_jets`` contains r_U, r_V, c_U, c_V.
    Numeric inputs broadcast. r and C must be positive, alpha nonnegative.
    ``region`` is "cap" (h_squared=1) or "exterior" (h_squared=0), with no
    adjustable vacuum parameter. Coordinates and Einstein radius are
    dimensionless; the two-metric is -C dU dV and L=log(C).

    psi_*_normalized are Psi_*/sqrt(alpha), evaluated by their continuous
    formulas at alpha=0. potential_increment is P-3h_squared, computed
    without subtracting nearly equal terms. m includes the cap vacuum mass.

    constraint_a_rhs means r_aa-L_a*r_a on the constraint surface. Pure
    second derivatives and L_a are not supplied, so no constraint residual
    is evaluated. For E_a=r_aa-L_a*r_a+r*Psi_a**2/2, mixed evolution implies
    d_V E_U=constraint_U_transport_coefficient*E_U and symmetrically in V.

    The actual differentiated mass satisfies
    m_U=mass_U_on_constraints+mass_U_residual_coefficient*E_U, and similarly
    in V. The named constrained rates alone are not derivatives of arbitrary
    input jets. The caller owns compatible characteristic data, the wall,
    regular-centre treatment and all physical domain bounds.
    """
    state = _fields(state, ("r", "C", "c"), "state")
    first_jets = _fields(first_jets, ("r_U", "r_V", "c_U", "c_V"), "first_jets")
    if not isinstance(region, str) or region not in ("cap", "exterior"):
        raise ValueError("region must be cap or exterior")
    data = _arrays(**state, **first_jets, alpha=alpha)
    if np.any(data["r"] <= 0) or np.any(data["C"] <= 0) or np.any(data["alpha"] < 0):
        raise ValueError("require r>0, C>0 and alpha>=0")
    try:
        with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
            r, C, c, alpha = (data[key] for key in ("r", "C", "c", "alpha"))
            rU, rV, cU, cV = (data[key] for key in ("r_U", "r_V", "c_U", "c_V"))
            u = c/r
            uU, uV = (cU-u*rU)/r, (cV-u*rV)/r
            q = alpha*u*u/6
            D = 1+q
            root_D = np.sqrt(D)
            h_squared = 1 if region == "cap" else 0
            increment = 3*q*(2+q) if h_squared else np.zeros_like(D)
            potential = 3*h_squared+increment
            sigmaUU, sigmaUV, sigmaVV = (alpha*uU*uU/D, alpha*uU*uV/D, alpha*uV*uV/D)
            rr = r*r
            radial_product = rU*rV
            rUV = -(radial_product+C/4)/r+3*h_squared*C*r/4+C*r*increment/4
            LUV = 2*(radial_product+C/4)/rr-sigmaUV
            # Separate the vacuum kernel from the scalar potential increment.
            # This preserves the alpha=0 cap cancellation without constructing
            # and subtracting its two vacuum mass/potential contributions.
            K = (radial_product+C*(1-h_squared*rr)/4)/rr-C*increment/12-sigmaUV/6
            result = {
                "r_UV": rUV, "log_C_UV": LUV, "c_UV": -K*c,
                "K": K, "m": r*(1+4*radial_product/C)/2,
                "D": D, "F": 1/D, "potential": potential,
                "potential_increment": increment,
                "psi_U_normalized": uU/root_D, "psi_V_normalized": uV/root_D,
                "constraint_U_rhs": -r*sigmaUU/2, "constraint_V_rhs": -r*sigmaVV/2,
                "constraint_U_transport_coefficient": -rV/r,
                "constraint_V_transport_coefficient": -rU/r,
                "mass_U_on_constraints": rr*potential*rU/2-rr*rV*sigmaUU/C,
                "mass_V_on_constraints": rr*potential*rV/2-rr*rU*sigmaVV/C,
                "mass_U_residual_coefficient": 2*r*rV/C,
                "mass_V_residual_coefficient": 2*r*rU/C,
            }
    except FloatingPointError as error:
        raise ValueError("bulk arithmetic overflow or invalid derived quantity") from error
    if not all(np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError("derived bulk quantities must be finite")
    return result
