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


def expanding_cap_wall_clock(*,radius,H,beta_minus,Rdot):
    """Future wall clock in a flat expanding de Sitter cap, through HR=1.

    ds2=-(1-H2*R2)dt2-2*HR*dt*dR+dR2, c=1. The supplied signed
    wall factors must belong to this cap. No wall trajectory is integrated.
    The clock formula uses the regular local boost, with no division by
    1-H2*R2; t is central proper time. Its origin is left unspecified.
    """
    R,H,beta,speed=map(float,(radius,H,beta_minus,Rdot))
    if (not all(math.isfinite(x) for x in (R,H,beta,speed)) or R<=0 or H<0
            or not math.isclose(beta*beta-speed*speed,1-(H*R)**2,
                rel_tol=2e-13,abs_tol=2e-13*max(1.,beta*beta,speed*speed,(H*R)**2))):
        raise ValueError('finite wall factors matching the expanding cap required')
    doppler=(beta+speed)/(1+H*R)
    if doppler<=0:
        raise ValueError('wall must lie in the future expanding flat patch')
    gamma=(doppler+1/doppler)/2
    return dict(d_cap_time_d_wall_time=gamma,wall_speed_in_cap_time=speed/gamma,
        local_radial_velocity=(doppler-1/doppler)/(2*gamma),
        inward_null_speed_in_cap_time=H*R-1,
        outside_central_horizon=H*R>1)


def central_radar_clock(*, radius, H, beta_minus, Rdot, Rddot):
    """Clock observables for an ideal transponder inside the static cap.

    All times are central proper lengths (c=1); the wall clock is proper
    tau. The test transponder has zero delay and negligible stress. It is
    not a reflective boundary condition for the model's transparent field.
    Echo stretch is dt_receive/dt_send; redshift is dt_receive/dtau.
    """
    R, H, beta, speed, acceleration = map(float,
        (radius, H, beta_minus, Rdot, Rddot))
    if (not all(math.isfinite(x) for x in (R, H, beta, speed, acceleration))
            or R <= 0 or H < 0 or H*R >= 1 or beta <= 0):
        raise ValueError('finite future wall inside the static cap, H*R<1, required')
    f = 1-(H*R)**2
    if not math.isclose(beta*beta-speed*speed, f,
                        rel_tol=2e-13, abs_tol=2e-13*max(1., beta*beta)):
        raise ValueError('wall tangent and signed cap normal disagree')
    distance = R if H == 0 else math.atanh(H*R)/H
    # Use the two null clock derivatives directly. Their product is 1/f.
    receive_rate = 1/(beta-speed)
    send_rate = 1/(beta+speed)
    beta_tau = speed*(acceleration-H*H*R)/beta
    return dict(radar_distance=distance, one_plus_redshift=receive_rate,
        echo_stretch=receive_rate/send_rate,
        log_redshift_rate=acceleration-beta_tau,
        send_clock_rate=send_rate, receive_clock_rate=receive_rate,
        cap_lapse_squared=f,
        scope='Calibrated test clocks on a quiet cap; no uncalibrated source inverse.')


def _radar_shape(lapse_squared):
    x = math.sqrt(1-lapse_squared)
    if x == 0:
        return 0., 1.
    optical = math.atanh(x)
    return optical, x/optical


def radar_clock_metric(*, one_plus_redshift, echo_stretch):
    """Measure the radial two-metric ON a timed transponder's worldline.

    For unique radial send/return null paths, radar coordinates are
    T=(t_plus+t_minus)/2, D=(t_plus-t_minus)/2. ANY radial metric is
    ds2=C(T,D)*(-dT2+dD2), including a time-dependent one. Proper wall
    normalization gives C=K/Z^2. No field equation, angular-radius law,
    static metric, or wall material law is imposed by this readout.
    """
    Z, K = map(float, (one_plus_redshift, echo_stretch))
    if not all(math.isfinite(x) and x > 0 for x in (Z, K)):
        raise ValueError('positive finite calibrated clock redshift and echo stretch required')
    C = (K/Z)/Z
    if not math.isfinite(C) or C <= 0:
        raise ValueError('radar metric factor must be finite and positive')
    return dict(radial_metric_factor=C, centre_clock_rate=(Z+Z/K)/2,
        radar_distance_rate=(Z-Z/K)/2, radar_velocity=(K-1)/(K+1),
        scope='Measured radial metric on the beacon, not a full spacetime or angular-radius reconstruction.')


def invert_central_radar_clock(*, radar_distance, one_plus_redshift,
                              echo_stretch, log_redshift_rate):
    """Recover a declared quiet cap and its local junction combinations.

    Radar and calibrated proper-clock data first recover F=K/Z^2 without
    knowing the wall motion. Identifying this F with a regular vacuum cap
    then gives H and R. General radial metrics need a separate areal-radius
    law; the clocks alone do not establish a de Sitter interior.
    """
    D, Z, K, drift = map(float, (radar_distance, one_plus_redshift,
                               echo_stretch, log_redshift_rate))
    if (not all(math.isfinite(x) for x in (D, Z, K, drift))
            or min(D, Z, K) <= 0):
        raise ValueError('positive finite distance, redshift and echo stretch required')
    f = radar_clock_metric(one_plus_redshift=Z,echo_stretch=K)['radial_metric_factor']
    if not 0 < f <= 1:
        raise ValueError('clock data incompatible with a nonnegative-curvature static cap')
    optical, radius_factor = _radar_shape(f)
    H, R = optical/D, D*radius_factor
    beta = (Z*f+1/Z)/2
    speed = (Z*f-1/Z)/2
    beta_radius = Z*(drift-H*H*R)
    a = R*R*(beta-R*beta_radius)/3
    c = (beta_radius+2*beta/R)/3
    return dict(H=H, radius=R, Lambda_cap=3*H*H,
        beta_minus=beta, Rdot=speed, beta_radius=beta_radius,
        mass_over_tension=a, junction_radius_coefficient=c,
        cap_lapse_squared=f,
        scope='Conditional quiet-cap inverse; Lambda_cap is reconstructed, not generated.')


def radar_junction_branches(inferred, *, host_H):
    """Keep every positive constant-tension branch for a specified host H.

    The cap mass is zero. Host curvature is an independent model input;
    without it the inverse determines M/k and (k^2+H_host^2-H^2)/(2k),
    not k, M and H_host separately.
    """
    H, a, c, host = map(float, (inferred['H'], inferred['mass_over_tension'],
                               inferred['junction_radius_coefficient'], host_H))
    if not all(math.isfinite(x) for x in (H, a, c, host)) or min(H, host, a) < 0:
        raise ValueError('nonnegative finite cap/host H and inferred M/k required')
    difference = host*host-H*H
    discriminant = c*c-difference
    if discriminant < 0:
        return []
    if discriminant == 0:
        roots = [c]
    else:
        q = c+math.copysign(math.sqrt(discriminant), c)
        roots = [q, difference/q] if q != 0 else []
    return [dict(tension=k, mass=k*a, host_H=host)
            for k in sorted(set(roots)) if k > 0]


def radar_cap_input_bounds(*, radar_distance, one_plus_redshift, echo_stretch):
    """Analytic box bounds on H^2 and R under the declared cap model.

    Input intervals are assumed instrument bounds, not observed errors.
    Correlations may narrow the box. Endpoints use float64; this is not an
    interval-arithmetic certification of transcendental rounding or a bound
    on an incorrect cap, a drifting standard, or a delayed transponder.
    """
    intervals = [tuple(map(float, pair)) for pair in
                 (radar_distance, one_plus_redshift, echo_stretch)]
    if any(len(pair) != 2 or not all(math.isfinite(x) for x in pair)
           or not 0 < pair[0] <= pair[1] for pair in intervals):
        raise ValueError('three positive finite ordered input intervals required')
    (D0, D1), (Z0, Z1), (K0, K1) = intervals
    f0, f1 = (K0/Z1)/Z1, (K1/Z0)/Z0
    if not 0 < f0 <= f1 <= 1:
        raise ValueError('entire input box must remain inside the admitted cap domain')
    y0, r0 = _radar_shape(f0)
    y1, r1 = _radar_shape(f1)
    curvature = [(y1/D1)**2, (y0/D0)**2]
    return dict(H_squared=curvature, Lambda_cap=[3*x for x in curvature],
        radius=[D0*r0, D1*r1], cap_lapse_squared=[f0, f1],
        scope='Conditional analytic input-box bounds evaluated in float64; no observational calibration.')
