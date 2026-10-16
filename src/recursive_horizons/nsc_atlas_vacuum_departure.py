"""Exact Schwarzschild geometry and nonlinear exterior departure prefixes.

The background is expressed in the original affine source labels (Ubar,x).
Subtracting its known integrals changes the quadrature, not the field equations.
No source, moving wall, cap, or time evolution is constructed here.
"""
from collections.abc import Mapping

import numpy as np
from scipy.integrate import cumulative_simpson, cumulative_trapezoid
from scipy.special import lambertw

from .nsc_atlas_exterior_characteristic_slice import (
    MAX_NODES, STATE_KEYS, JET_KEYS, _array, _mapping, _check_tube,
)


SOURCE_KEYS = ("r", "log_C", "c", "r_x", "log_C_x", "c_x")
HORIZON_RADIUS = .8
MASS = .4


def _budget(value):
    if not callable(getattr(value, "check", None)):
        raise ValueError("budget must expose a live cumulative check(count)")
    value.check()
    return value


def _coordinates(Ubar, x):
    U, x = _array(Ubar, "Ubar"), _array(x, "x")
    try:
        U, x = np.broadcast_arrays(U, x)
    except ValueError as error:
        raise ValueError("vacuum coordinates must broadcast") from error
    if U.ndim > 1 or U.size > MAX_NODES:
        raise ValueError("vacuum queries require at most 513 scalar or one-dimensional points")
    if np.any((x <= 0) | (x >= HORIZON_RADIUS)) or np.any((U < 0) | (U >= 4*HORIZON_RADIUS)):
        raise ValueError("vacuum coordinates leave 0<x<.8, 0<=Ubar<3.2")
    return U, x


def vacuum_slice(Ubar, x, *, budget):
    """Return exact trapped Schwarzschild values and analytic first/second jets.

    The original affine source is Ubar=0, r=x, C=1. The principal Lambert-W
    branch gives 0<r<.8 and mass .4. The algebraic coordinate guard is wider
    than the admitted sourced wedge; its narrower domain belongs to the caller.
    Scalar or broadcast one-dimensional queries have at most 513 points.
    One charge per point covers this algebraic evaluation, with no quadrature.
    """
    U, x = _coordinates(Ubar, x)
    budget = _budget(budget)
    budget.check(int(U.size))
    try:
        with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
            a = HORIZON_RADIUS
            product = (1-U/(4*a))*(1-x/a)*np.exp(x/a)
            r = a*(1+lambertw(-product/np.e, 0).real)
            # Preserve the exact affine anchor, including its first jets.
            r = np.where(U == 0, x, r)
            if np.any((r <= 0) | (r >= a)):
                raise ValueError("vacuum radius leaves the representable trapped branch")
            L = np.log(x/r)+(x-r)/a
            C = np.exp(L)
            p = C*(a-x)/(4*x)
            R = (1-U/(4*a))*C
            B = 1/r+1/a
            LU = -B*p
            Lx = 1/x+1/a-B*R
            rUU, rxx = LU*p, Lx*R
            values = {
                "Ubar": U.copy(), "x": x.copy(), "r": r, "L": L, "C": C,
                "r_Ubar": p, "r_x": R, "log_C_Ubar": LU, "log_C_x": Lx,
                "r_Ubar_Ubar": rUU, "r_xx": rxx,
                "log_C_Ubar_Ubar": p*p/(r*r)-B*rUU,
                "log_C_xx": -1/(x*x)+R*R/(r*r)-B*rxx,
                "r_Ubar_x": -MASS*C/(2*r*r),
                "log_C_Ubar_x": MASS*C/(r*r*r),
            }
    except (FloatingPointError, OverflowError) as error:
        raise ValueError("nonfinite vacuum arithmetic") from error
    if not all(np.all(np.isfinite(value)) for value in values.values()):
        raise ValueError("nonfinite vacuum output")
    budget.check()
    return values


def departure_slice(Ubar, states, departure_first_jets, source, *, x, alpha,
                    tube, budget, quadrature_order=2):
    """Reconstruct full nonlinear longitudinal jets using vacuum differences.

    states contains physical r,L,c on a uniform 3..513-node Ubar grid. The
    supplied transverse jets differentiate physical-minus-vacuum values;
    the source offset is spatially constant and therefore has no Ubar jet.
    The vacuum scalar is zero, including when alpha=0 keeps a normalized
    test scalar present. Physical source values must exactly match node zero.
    A full SourceSnapshot row is accepted; its six required values/first jets
    are preserved. This operator neither replaces its curvature nor enforces
    a source constraint. x is one fixed incoming label.

    Subtracted radial and lapse integrals retain all nonlinear terms. The
    canonical scalar integral is the full equation. Seven N-point passes are
    charged: vacuum, normalization, three evolution prefixes, and two vacuum
    quadratures used only to report the removed discretization error. Order
    2/4 uses cumulative trapezoid/Simpson as in the original prefix operator.
    The quadrature correction vanishes under refinement; no error enclosure
    or pure-null constraint result follows from these mixed-equation jets.
    """
    U = _array(Ubar, "Ubar")
    x = float(_array(x, "x", ()))
    if (U.ndim != 1 or not 3 <= U.size <= MAX_NODES or U[0] != 0
            or np.any(np.diff(U) <= 0)
            or not np.allclose(np.diff(U), U[-1]/(U.size-1), rtol=5e-13, atol=0)):
        raise ValueError("Ubar must be 3..513 uniform increasing nodes starting at zero")
    _coordinates(U, x)
    if type(quadrature_order) is not int or quadrature_order not in (2, 4):
        raise ValueError("quadrature_order must be 2 or 4")
    alpha = float(_array(alpha, "alpha", ()))
    if alpha < 0:
        raise ValueError("alpha must be nonnegative")
    states = _mapping(states, STATE_KEYS, "states", U.shape)
    jets = _mapping(departure_first_jets, JET_KEYS, "departure_first_jets", U.shape)
    if not isinstance(source, Mapping) or not set(SOURCE_KEYS) <= set(source):
        raise ValueError("source must include values and first x jets")
    src = {key: float(_array(source[key], "source."+key, ())) for key in SOURCE_KEYS}
    if any(states[key][0] != src["log_C" if key == "L" else key] for key in STATE_KEYS):
        raise ValueError("physical node zero must exactly preserve the supplied source")
    _check_tube(states, tube)
    budget = _budget(budget)
    N = int(U.size)
    vacuum = vacuum_slice(U, x, budget=budget)

    def prefix(values):
        budget.check(N)
        if quadrature_order == 2:
            return cumulative_trapezoid(values, x=U, initial=0.)
        return cumulative_simpson(values, x=U, initial=0.)

    try:
        with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
            budget.check(N)
            r, L, c = (states[key] for key in STATE_KEYS)
            rv, Cv, pv, Rv = (vacuum[key] for key in ("r", "C", "r_Ubar", "r_x"))
            dr, dL = r-rv, L-vacuum["L"]
            dC = Cv*np.expm1(dL)
            C = np.exp(L)
            p = pv+jets["r_Ubar"]
            LU = vacuum["log_C_Ubar"]+jets["log_C_Ubar"]
            cU = jets["c_Ubar"]
            u = c/r
            D = 1+alpha*u*u/6
            root_D = np.sqrt(D)
            fU = (cU-u*p)/(r*root_D)
            source_fx = (src["c_x"]-u[0]*src["r_x"])/(r[0]*root_D[0])

            integral_dC = prefix(dC/4)
            source_radial_offset = (src["r"]-x)+src["r"]*(src["r_x"]-1)
            dR = (source_radial_offset-integral_dC-Rv*dr)/r
            R = Rv+dR
            R[0] = src["r_x"]
            dR[0] = src["r_x"]-1
            integral_scalar = prefix(R*fU)
            fx = (r[0]*source_fx-integral_scalar)/r
            cx = u*R+r*root_D*fx
            cx[0] = src["c_x"]

            dp = jets["r_Ubar"]
            Av = pv*Rv+Cv/4
            dA = pv*dR+Rv*dp+dp*dR+dC/4
            lapse_difference = (2*dA/(r*r)
                -2*Av*dr*(r+rv)/(r*r*rv*rv)-alpha*fU*fx)
            integral_lapse_difference = prefix(lapse_difference)
            dLx = src["log_C_x"]+integral_lapse_difference
            Lx = vacuum["log_C_x"]+dLx
            Lx[0] = src["log_C_x"]

            # These two extra quadratures expose the change to the old scheme;
            # they never enter the physical prefixes above.
            vacuum_C_quadrature = prefix(Cv/4)
            vacuum_lapse_quadrature = prefix(vacuum["log_C_Ubar_x"])
            vacuum_C_exact = x-rv*Rv
            removed_radial = (vacuum_C_quadrature-vacuum_C_exact)/r
            removed_lapse = vacuum["log_C_x"]-vacuum_lapse_quadrature
            full_lapse_rhs = 2*(p*R+C/4)/(r*r)-alpha*fU*fx
            radial_rhs = -(p*R+C/4)/r
            kernel = (p*R+C/4)/(r*r)-alpha*fU*fx/6
            arrays = {
                "Ubar": U, "r_Ubar": p, "log_C_Ubar": LU, "c_Ubar": cU,
                "r_x": R, "log_C_x": Lx, "c_x": cx,
                "delta_r_x": dR, "delta_log_C_x": dLx, "delta_c_x": cx.copy(),
                "delta_r": dr, "delta_log_C": dL, "delta_C": dC,
                "D": D, "f_Ubar": fU, "f_x": fx,
                "integral_delta_C_over_four": integral_dC,
                "integral_r_x_f_Ubar": integral_scalar,
                "integral_lapse_difference": integral_lapse_difference,
                "lapse_difference_integrand": lapse_difference,
                "removed_radial_quadrature_error": removed_radial,
                "removed_lapse_quadrature_error": removed_lapse,
                "r_Ubar_x": radial_rhs, "log_C_Ubar_x": full_lapse_rhs,
                "c_Ubar_x": -kernel*c,
            }
    except (FloatingPointError, OverflowError) as error:
        raise ValueError("nonfinite vacuum departure arithmetic") from error
    if not all(np.all(np.isfinite(value)) for value in arrays.values()):
        raise ValueError("nonfinite vacuum departure output")
    budget.check()
    return {**arrays, "vacuum": vacuum, "source": src,
            "settings": {"nodes": N, "quadrature_order": quadrature_order,
                         "charged_evaluations": 7*N},
            "scope": {"source_integrations": 0, "time_steps": 0,
                      "pure_constraints_evaluated": False,
                      "quadrature_error_is_enclosed": False,
                      "vacuum_quadrature_drives_evolution": False}}
