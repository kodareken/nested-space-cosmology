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


class TwoWallHistoryBudget:
    """Count scalar geometry/RHS work across one finite application."""

    def __init__(self, limit=20000, cpu_seconds=60.):
        import time
        self.limit = int(limit)
        self.cpu_seconds = float(cpu_seconds)
        self.count = 0
        self.started = time.process_time()

    def add(self, count=1):
        import time
        self.count += int(count)
        if self.count > self.limit:
            raise RuntimeError('two-wall scalar evaluation allowance exhausted')
        if time.process_time()-self.started > self.cpu_seconds:
            raise RuntimeError('two-wall process-CPU allowance exhausted')


class TwoWallHostHistory:
    """The existing quiet two-wall state in one oriented host atlas.

    The bridge has T=0 and r_star(.5*ell)=0; Lake's scale is L0=ell.
    Clock origins come from the original slice, not independent wall resets.
    Only the admitted tau/ell in [0,.04] is continued. Equal proper times on
    different walls are not asserted simultaneous. No new field is evolved.
    """

    def __init__(self, *, order=64, length_scale=1., budget=None):
        from recursive_horizons.nsc_atlas_nested_slice import kottler_horizons
        if order not in (32, 64):
            raise ValueError('two-wall quadrature order must be 32 or 64')
        self.ell = float(length_scale)
        if not math.isfinite(self.ell) or self.ell <= 0:
            raise ValueError('positive finite length scale required')
        self.order = order
        self.budget = budget if budget is not None else TwoWallHistoryBudget()
        self.nodes, self.weights = np.polynomial.legendre.leggauss(order)
        self.mass, self.H, self.tension = .01*self.ell, 1/self.ell, 2/self.ell
        roots = kottler_horizons(self.mass, self.H, self.budget)
        self.rb, self.rc = roots['rb'], roots['rc']
        self.ro = -self.rb-self.rc
        self.kappa = (1-3*(self.H*self.rb)**2)/(2*self.rb)
        self.reference = .5*self.ell
        self.initial_walls = {name: self.wall(name, radius*self.ell)
                              for name, radius in (('inner', .8), ('outer', 1.5))}
        self._initial_events = {}
        self._outer_label = None

    def _integral(self, function, lower, upper):
        if lower == upper:
            return 0.
        points = (upper+lower)/2+(upper-lower)/2*self.nodes
        self.budget.add(points.size)
        return float((upper-lower)/2*self.weights @ np.asarray(function(points)))

    def lapse(self, radius):
        radius = np.asarray(radius)
        return 1-2*self.mass/radius-(self.H*radius)**2

    def wall(self, name, radius):
        from recursive_horizons.nsc_atlas_nested_slice import tension_wall
        self.budget.add()
        if name == 'inner':
            result = tension_wall(radius, 0., self.H, self.mass, self.H, self.tension)
            result['host_beta'] = result['beta_plus']
        elif name == 'outer':
            result = tension_wall(radius, self.mass, self.H, .4*self.ell, 0., self.tension)
            result['host_beta'] = result['beta_minus']
        else:
            raise ValueError('wall must be inner or outer')
        return result

    def _wall_arrays(self, name, radius):
        """Vectorized signed motion used only inside counted quadrature."""
        radius = np.asarray(radius)
        if name == 'inner':
            minus_mass, minus_H, plus_mass, plus_H = 0., self.H, self.mass, self.H
        elif name == 'outer':
            minus_mass, minus_H, plus_mass, plus_H = self.mass, self.H, .4*self.ell, 0.
        else:
            raise ValueError('wall must be inner or outer')
        f_minus = 1-2*minus_mass/radius-(minus_H*radius)**2
        f_plus = 1-2*plus_mass/radius-(plus_H*radius)**2
        beta_minus = (f_minus-f_plus+(self.tension*radius)**2)/(2*self.tension*radius)
        speed = np.sqrt(beta_minus**2-f_minus)
        beta_host = beta_minus-self.tension*radius if name == 'inner' else beta_minus
        return speed, beta_host

    def initial_frame(self, branch, radius):
        """Original (Rprime,V), never a newly chosen taper."""
        from recursive_horizons.nsc_atlas_nested_slice import (
            positive_sheet_velocity, taper_velocity, transition_velocity,
        )
        radius = np.asarray(radius)
        if branch == 'inner':
            if np.any(radius < .1*self.ell) or np.any(radius > .8*self.ell):
                raise ValueError('inner source outside original left taper')
            V, _ = taper_velocity(radius, .1*self.ell, .8*self.ell,
                                  self.initial_walls['inner']['V_plus'])
            Rp = -np.sqrt(self.lapse(radius)+V*V)
        elif branch == 'outer':
            if np.any(radius < .5*self.ell) or np.any(radius > 1.5*self.ell):
                raise ValueError('outer source outside original right host')
            transition, _ = transition_velocity(radius, .5*self.ell, .9*self.ell,
                                                 self.mass, self.H)
            sheet, _ = positive_sheet_velocity(radius, self.mass, self.H)
            V = np.where(radius < .9*self.ell, transition, sheet)
            Rp = np.sqrt(self.lapse(radius)+V*V)
        else:
            raise ValueError('branch must be inner or outer')
        return Rp, V

    def tortoise_regular_part(self, radius):
        """Smooth remainder at rb from the exact partial fractions of 1/f."""
        radius = float(radius)
        value = 0.
        for root in (self.rc, self.ro):
            derivative = 2*self.mass/root**2-2*self.H**2*root
            value += math.log(abs((radius-root)/(self.reference-root)))/derivative
        return value

    def tortoise(self, radius):
        radius = float(radius)
        if radius == self.rb:
            return -math.inf
        if radius == self.rc:
            return math.inf
        return (math.log(abs((radius-self.rb)/(self.reference-self.rb)))/(2*self.kappa)
                + self.tortoise_regular_part(radius))

    def _initial_null(self, branch, radius):
        """p on the left and q on the right, with the shared bridge origin."""
        radius = float(radius)
        if branch == 'inner':
            def integrand(points):
                Rp, V = self.initial_frame('inner', points)
                return V/(self.lapse(points)*Rp)
            T = self._integral(integrand, .1*self.ell, radius)
            return T+self.tortoise(radius)
        if branch == 'outer':
            def integrand(points):
                Rp, V = self.initial_frame('outer', points)
                return -1/(Rp*(V+Rp))
            split = min(radius, .9*self.ell)
            return (self._integral(integrand, .5*self.ell, split)
                    + self._integral(integrand, split, radius))
        raise ValueError('branch must be inner or outer')

    def lake_event(self, radius, q, sign, *, speed=None, beta=None):
        """Common Lake coordinates, optionally with a wall proper tangent."""
        self.budget.add()
        u = sign*self.ell*math.exp(-self.kappa*q)
        w = (radius-self.rb)/(self.kappa*u)
        chart = regular_kottler_chart(u, w, mass=self.mass, H=self.H, horizon=self.rb)
        result = dict(u=u, w=w, F=float(chart['F']),
                      radius_residual=float(chart['R'])-radius)
        if speed is not None:
            q_tau = 1/(beta+speed)
            u_tau = -self.kappa*u*q_tau
            w_tau = speed/(self.kappa*u)-w*u_tau/u
            result.update(u_tau=u_tau, w_tau=w_tau,
                normalization_residual=result['F']*u_tau*u_tau+2*u_tau*w_tau+1,
                radial_rate_residual=self.kappa*(w*u_tau+u*w_tau)-speed)
        return result

    def initial_event(self, branch, radius):
        self.budget.add()
        radius = float(radius)
        Rp, V = map(float, self.initial_frame(branch, radius))
        null = self._initial_null(branch, radius)
        star = self.tortoise(radius)
        q = null-2*star if branch == 'inner' else null
        T = null-star if branch == 'inner' else null+star
        lake = self.lake_event(radius, q, -1 if branch == 'inner' else 1)
        q_l = -1/(Rp+V)
        u_l = -self.kappa*lake['u']*q_l
        w_l = Rp/(self.kappa*lake['u'])-lake['w']*u_l/lake['u']
        N_u, N_w = -u_l, lake['F']*u_l+w_l
        lake.update(u_l=u_l, w_l=w_l, N_u=N_u, N_w=N_w,
            spatial_norm_residual=lake['F']*u_l*u_l+2*u_l*w_l-1,
            normal_norm_residual=lake['F']*N_u*N_u+2*N_u*N_w+1,
            orthogonality_residual=lake['F']*u_l*N_u+u_l*N_w+w_l*N_u,
            normal_radius_residual=self.kappa*(lake['w']*N_u+lake['u']*N_w)-V)
        return dict(radius=radius, T=T if math.isfinite(T) else None,
                    null_coordinate=null, q=q, Rprime=Rp, V=V, lake=lake)

    def bridge_event(self, y):
        """Smooth signed static bridge, including its bifurcation sphere."""
        from recursive_horizons.nsc_atlas_nested_slice import bridge_chart
        self.budget.add()
        y = float(y)
        radius = self.rb+y*y
        remainder = self.tortoise_regular_part(radius)
        factor = self.ell/math.sqrt(self.reference-self.rb)*math.exp(self.kappa*remainder)
        u, w = y*factor, y/(self.kappa*factor)
        remainder_R = sum(1/((2*self.mass/root**2-2*self.H**2*root)*(radius-root))
                          for root in (self.rc, self.ro))
        u_y = factor*(1+2*self.kappa*y*y*remainder_R)
        w_y = (1-2*self.kappa*y*y*remainder_R)/(self.kappa*factor)
        slice_row = bridge_chart(y, self.rb, self.H, self.mass)
        dl_dy = slice_row['dl_dy']
        chart = regular_kottler_chart(u, w, mass=self.mass, H=self.H, horizon=self.rb)
        F = float(chart['F'])
        u_l, w_l = u_y/dl_dy, w_y/dl_dy
        N_u, N_w = -u_l, F*u_l+w_l
        return dict(radius=radius, u=u, w=w, T=0., F=F, dl_dy=dl_dy,
            induced_dy_squared=F*u_y*u_y+2*u_y*w_y,
            spatial_norm_residual=F*u_l*u_l+2*u_l*w_l-1,
            normal_radius_residual=self.kappa*(w*N_u+u*N_w),
            future_N_u=N_u, future_N_w=N_w)

    def initial_wall(self, name):
        if name not in self._initial_events:
            radius = self.initial_walls[name]['radius']
            self._initial_events[name] = self.initial_event(name, radius)
        return self._initial_events[name]

    def wall_event(self, name, proper_time):
        from scipy.optimize import brentq
        tau = float(proper_time)
        if not math.isfinite(tau) or not 0 <= tau <= .04*self.ell:
            raise ValueError('wall proper time outside admitted [0,.04 ell]')
        initial = self.initial_wall(name)
        R0 = initial['radius']
        upper = (.85 if name == 'inner' else 1.6)*self.ell
        def elapsed(radius):
            return self._integral(lambda x: 1/self._wall_arrays(name, x)[0], R0, radius)
        radius = R0 if tau == 0 else brentq(lambda r: elapsed(r)-tau, R0, upper,
                                           xtol=5e-15*self.ell, rtol=8e-15)
        def null_rate(points):
            speed, beta = self._wall_arrays(name, points)
            denominator = beta-speed if name == 'inner' else beta+speed
            return 1/(speed*denominator)
        null = initial['null_coordinate']+self._integral(null_rate, R0, radius)
        star = self.tortoise(radius)
        q = null-2*star if name == 'inner' else null
        T = null-star if name == 'inner' else null+star
        wall = self.wall(name, radius)
        lake = self.lake_event(radius, q, -1 if name == 'inner' else 1,
                               speed=wall['Rdot'], beta=wall['host_beta'])
        return dict(proper_time=tau, radius=radius, T=T, q=q,
            null_coordinate=null, null_kind='p=T+rstar' if name == 'inner' else 'q=T-rstar',
            Rdot=wall['Rdot'], host_beta=wall['host_beta'], lake=lake,
            junction_residuals={key:wall[key] for key in
                ('mass_map_residual', 'normal_force_residual',
                 'minus_rate_residual', 'plus_rate_residual')},
            angular_jump_residual=wall['angular_jump']+self.tension,
            temporal_jump_residual=wall['temporal_jump']+self.tension)

    def outer_initial_label(self):
        from recursive_horizons.nsc_atlas_nested_slice import bridge_chart
        if self._outer_label is None:
            left = self._integral(lambda x: -1/self.initial_frame('inner', x)[0],
                                  .1*self.ell, .8*self.ell)
            def bridge(points):
                return bridge_chart(points, self.rb, self.H, self.mass)['dl_dy']
            neck = (self._integral(bridge, -math.sqrt(.1*self.ell-self.rb), 0.)
                    + self._integral(bridge, 0., math.sqrt(.5*self.ell-self.rb)))
            right = (self._integral(lambda x: 1/self.initial_frame('outer', x)[0],
                                    .5*self.ell, .9*self.ell)+.6*self.ell)
            self._outer_label = .8*self.ell+left+neck+right
        return self._outer_label

    def characteristic_origin(self, name, event):
        """Unique original-slice origin of the incoming host characteristic."""
        from scipy.optimize import brentq
        R0 = self.initial_walls[name]['radius']
        lower = (.7 if name == 'inner' else 1.3)*self.ell
        target = event['null_coordinate']
        if event['proper_time'] == 0:
            radius = R0
        else:
            radius = brentq(lambda r: self._initial_null(name, r)-target, lower, R0,
                             xtol=5e-15*self.ell, rtol=8e-15)
        initial = self.initial_event(name, radius)
        Rp, V = initial['Rprime'], initial['V']
        if name == 'inner':
            label = .8*self.ell+self._integral(
                lambda x: -1/self.initial_frame(name, x)[0], radius, R0)
            derivative = (event['host_beta']-event['Rdot'])/(Rp-V)
        else:
            label = self.outer_initial_label()-self._integral(
                lambda x: 1/self.initial_frame(name, x)[0], radius, R0)
            derivative = -(event['host_beta']+event['Rdot'])/(Rp+V)
        return dict(radius=radius, proper_label=label, arrival_d_tau_d_label=derivative,
            null_matching_residual=initial['null_coordinate']-target,
            initial_event=initial)
