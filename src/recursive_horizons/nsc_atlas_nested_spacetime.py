"""Regular quiet Kottler bulk charts and a specified central test flash.

The chart is the algebraic form of Lake gr-qc/0507031v2 (5)--(10).
It does not select global horizon identifications or a matter preparation.
The flash is confined to the de Sitter cap and has negligible backreaction.
"""
from __future__ import annotations

import math

import numpy as np


def regular_kottler_chart(u, w, *, mass, H, horizon):
    """Return the radial metric and jets, with the root error kept explicit.

    u,w have length units; ds2 = F du2 + 2 du dw + R2 dOmega2.
    A simple positive horizon is required. At finite coordinates the chart
    remains regular for R>0. Another chart may be needed at coordinate infinity.
    """
    mass, H, C = map(float, (mass, H, horizon))
    if (not all(math.isfinite(x) for x in (mass, H, C))
            or mass < 0 or H < 0 or C <= 0):
        raise ValueError("finite nonnegative mass/H and positive horizon required")
    root_residual = 1-2*mass/C-H*H*C*C
    if abs(root_residual) > 128*np.finfo(float).eps*max(1., 2*mass/C, H*H*C*C):
        raise ValueError("horizon is not a root of the supplied mass and H")
    derivative_factor = 1-3*H*H*C*C
    if abs(derivative_factor) <= 128*np.finfo(float).eps*max(1.,3*H*H*C*C):
        raise ValueError("a numerically resolved simple horizon is required")
    kappa = derivative_factor/(2*C)
    u, w = np.broadcast_arrays(np.asarray(u, dtype=float), np.asarray(w, dtype=float))
    if not np.all(np.isfinite(u)) or not np.all(np.isfinite(w)):
        raise ValueError("chart coordinates must be finite")
    R = C+kappa*u*w
    if not np.all(np.isfinite(R)) or np.any(R <= 0):
        raise ValueError("regular chart is restricted to finite R>0")
    B = 1/C-H*H*C
    F = w*w*(H*H+B/R)
    F_u = -kappa*B*w**3/R**2
    F_w = 2*w*(H*H+B/R)-kappa*B*u*w*w/R**2
    F_ww = 2*H*H+2*B*C*C/R**3
    return {
        "R": R, "F": F, "F_u": F_u, "F_w": F_w, "F_ww": F_ww,
        "R_u": kappa*w, "R_w": kappa*u, "R_uw": kappa,
        "radial_metric_determinant": -1., "kappa": kappa,
        "supplied_mass": mass, "chart_mass": C*(1-H*H*C*C)/2,
        "horizon_root_residual": root_residual,
        "mass_roundoff_difference": C*(1-H*H*C*C)/2-mass,
    }


def static_slice_frame(radius, V, *, mass, H, orientation):
    """Oriented s,N in a single static Killing coordinate (T,R).

    The left branch has future N(T)<0. Its locally future coordinate -T
    is a different coordinate and reverses all T derivatives together.
    """
    R, V, mass, H = map(float, (radius, V, mass, H))
    if (not all(math.isfinite(x) for x in (R, V, mass, H))
            or R <= 0 or mass < 0 or H < 0 or orientation not in (-1, 1)):
        raise ValueError("finite static data and orientation +/-1 required")
    f = 1-2*mass/R-H*H*R*R
    if f <= 0:
        raise ValueError("static frame requires f>0; use a regular chart at a horizon")
    Rp = orientation*math.sqrt(f+V*V)
    return {"f": f, "Rprime": Rp,
            "s_T": V/f, "s_R": Rp, "N_T": Rp/f, "N_R": V,
            "T_l": V/f, "T_R": V/(f*Rp),
            "outgoing_l": -1/(V+Rp), "advanced_N": (V+Rp)/f}


def central_inward_flash(*, radius, H, beta_minus, Rdot):
    """Exact cap-side arrival and frequency ratio for a wall-emitted flash.

    t_emit=0 in the flat cap chart, central receiver, inward geometric optics,
    prescribed wall proper frequency. Returned time is a length (c=1).
    Input/model uncertainties are not supplied by floating-point precision.
    """
    R, H, beta, speed = map(float, (radius, H, beta_minus, Rdot))
    if (not all(math.isfinite(x) for x in (R, H, beta, speed))
            or R <= 0 or H < 0 or not H*R < 1):
        raise ValueError("finite data with R>0,H>=0 and H*R<1 required")
    if not math.isclose(beta*beta-speed*speed, 1-(H*R)**2,
                        rel_tol=2e-13, abs_tol=2e-13):
        raise ValueError("wall factors do not match the declared de Sitter cap")
    a = (beta+speed)/(1+H*R)
    if a <= 0:
        raise ValueError("wall emitter must be future timelike in the cap frame")
    gamma, gamma_v = (a+1/a)/2, (a-1/a)/2
    time_to_centre = R if H == 0 else -math.log1p(-H*R)/H
    frequency_ratio = a/(1-H*R)
    return {"central_proper_time": time_to_centre,
            "one_plus_redshift": frequency_ratio, "redshift": frequency_ratio-1,
            "emission_doppler_factor": a, "emitter_gamma": gamma,
            "emitter_local_speed": gamma_v/gamma,
            "cap_clock_origin": "flat-slicing emission event t=0",
            "probe": "inward test flash from material wall to central geodesic"}
