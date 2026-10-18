"""Wall-pinned time chart for the quiet cap, material wall and Schwarzschild exterior.

The construction pulls the exact de Sitter cap and the outgoing
Eddington--Finkelstein exterior back to one wall-pinned coordinate. That
coordinate equals the initial proper distance minus the initial wall
distance; it is not future exterior proper distance. This is a coordinate
chart on the fixed quiet sources, not a new stress, a state, a Hadamard
condition or a Dirac evolution.

Model units are fixed at ell=1, H=1, k_wall=2 and M=.4. The wall radius stays
on the original strip .45 <= R <= .55. Quadrature is only for the retarded
coordinate and the wall clocks; every producer derivative is analytic.
"""
from __future__ import annotations

import math
import time

import numpy as np

from .nsc_atlas_cap_wall_response import reference_coefficients
from .nsc_atlas_connection_slice import smooth_step
from .nsc_atlas_dirac_evolution import radial_time_coefficients
from .nsc_atlas_dirac_slice import source_anchor_slice


class ChartDomainError(ValueError):
    """Fixed chart sample is not spacelike or not future-directed.

    The details locate the sampled failure. They are not a gauge search and
    not a certified boundary of the domain.
    """

    def __init__(self, message, details):
        super().__init__(message)
        self.details = details


CPU_SECONDS = 120
MAX_QUADRATURE_POINTS = 200000
ELL = 1.0
H = 1.0
K_WALL = 2.0
M = 0.4
R_ANCHOR = 0.45
R_STRIP = (0.45, 0.55)
TAPER = (1.0, 1.2)
DISPLAY_RADIUS = 4.0
CUTOFF = 1.0
UNDERFLOW = 1e-10
ALPHA = np.array([[0.0, -1j], [1j, 0.0]])
_EYE = np.eye(2)
_GAUSS = {}


def _order(order):
    if isinstance(order, bool) or not isinstance(order, int) or order not in (32, 64):
        raise ValueError('quadrature order must be the fixed comparison value 32 or 64')
    return order


def _finite(value, name):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.size == 0:
        raise ValueError(name + ' must contain finite real numbers')
    array = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(name + ' must contain finite real numbers')
    return array


def _scalar(value, name='value'):
    """One Python float. Size-one arrays are accepted without a ndim>0 cast."""
    array = _finite(value, name)
    if array.size != 1:
        raise ValueError(name + ' must be one finite real number')
    return array.reshape(()).item()


def _one_wall_radius(radius, name):
    value = _scalar(radius, name)
    if not R_STRIP[0] <= value <= R_STRIP[1]:
        raise ValueError(name + ' requires one wall radius in [.45, .55]')
    return value


def _gauss(order):
    if order not in _GAUSS:
        _GAUSS[order] = np.polynomial.legendre.leggauss(order)
    return _GAUSS[order]


def _integrate(func, lower, upper, order, counter):
    """Gauss integral of an analytic scalar integrand on each [lower, upper]."""
    lower, upper = np.broadcast_arrays(_finite(lower, 'lower'), _finite(upper, 'upper'))
    shape = lower.shape
    lower, upper = lower.reshape(-1), upper.reshape(-1)
    result = np.zeros(lower.size, dtype=float)
    span = upper - lower
    active = span != 0.0
    if np.any(active):
        nodes, weights = _gauss(order)
        half = 0.5 * span[active]
        mid = 0.5 * (lower[active] + upper[active])
        sample = mid[:, None] + half[:, None] * nodes
        if counter[0]+int(sample.size) > MAX_QUADRATURE_POINTS:
            raise RuntimeError('fixed coordinate-chart quadrature budget exhausted')
        values = np.asarray(func(sample), dtype=float)
        if values.shape != sample.shape or not np.all(np.isfinite(values)):
            raise ValueError('quadrature integrand must return finite samples of the same shape')
        result[active] = half * (values * weights).sum(axis=1)
        counter[0] += int(sample.size)
    return result.reshape(shape)


def _step_second(x):
    """Second derivative of the logistic step on the open interval away from its ends."""
    L = -1.0 / x + 1.0 / (1.0 - x)
    Lp = 1.0 / x**2 + 1.0 / (1.0 - x)**2
    Lpp = -2.0 / x**3 + 2.0 / (1.0 - x)**3
    weight = np.exp(-np.abs(L))
    product = weight / (1.0 + weight)**2
    one_minus_two = np.where(L >= 0.0, -(1.0 - weight) / (1.0 + weight), (1.0 - weight) / (1.0 + weight))
    return product * (one_minus_two * Lp**2 + Lpp)


def _cutoff(y):
    """S=1-smooth_step, supported on [0, 1], with analytic first and second derivatives.

    The step is flat at both ends. Within 1e-10 of either end the positive-order
    derivatives underflow to zero before the rational prefactors overflow.
    """
    y = _finite(y, 'y')
    flat = np.atleast_1d(y)
    step, step_y = smooth_step(flat)
    step_y = np.array(step_y, dtype=float, copy=True)
    step_yy = np.zeros(flat.shape, dtype=float)
    band = (flat <= UNDERFLOW) | (flat >= 1.0 - UNDERFLOW)
    step_y[band] = 0.0
    inner = (flat > UNDERFLOW) & (flat < 1.0 - UNDERFLOW)
    if np.any(inner):
        step_yy[inner] = _step_second(flat[inner])
    return ((1.0 - step).reshape(y.shape), (-step_y).reshape(y.shape), (-step_yy).reshape(y.shape))


def _wall_geometry(radius):
    """Exact quiet-wall jets. beta_plus=0 is an ordinary point: no division by it or by f."""
    coeff = reference_coefficients(radius)
    R, v, bm, bp = coeff['x'], coeff['v'], coeff['beta_minus'], coeff['beta_plus']
    bm_R = coeff['beta_minus_x']
    A = 1.0 - R * R
    root = np.sqrt(A)
    gamma, gamma_v = bm / root, v / root
    slope = R / root
    slope_R = 1.0 / A**1.5
    v_R = (bm * bm_R + R) / v
    E = K_WALL * slope * v
    J = root - K_WALL * slope * bm
    E_R = K_WALL * (v_R * slope + v * slope_R)
    J_R = -R / root - K_WALL * (bm_R * slope + bm * slope_R)
    speed = E + J
    outgoing = v * (bp + v)
    if np.any(speed <= 0.0) or np.any(outgoing <= 0.0):
        raise ValueError('wall tangent or outgoing retarded derivative is not positive on this radius')
    return dict(R=R, v=v, v_R=v_R, beta_minus=bm, beta_plus=bp, beta_minus_R=bm_R,
        A=A, root=root, gamma=gamma, gamma_v=gamma_v, E=E, J=J, E_R=E_R, J_R=J_R,
        speed=speed, U_R=1.0 / outgoing, s=gamma_v / gamma)


def wall_time_data(radius, chart, order=32):
    """Clocks and outgoing-EF wall embedding at one or more radii on the original strip.

    t=q0*(R-.45) normalizes the initial wall lapse to one. U, cap time and wall
    proper time are integrals from .45; their t and R derivatives are exact.
    """
    order = _order(order)
    geometry = _wall_geometry(radius)
    q0 = chart['q0']
    R = geometry['R']
    t = q0 * (R - R_ANCHOR)
    anchor = float(R_ANCHOR)
    U = _integrate(lambda sample: _wall_geometry(np.clip(sample, R_ANCHOR, R_STRIP[1]))['U_R'],
        anchor, R, order, chart['quadrature_points'])
    Tcap = _integrate(lambda sample: _cap_time_density(np.clip(sample, R_ANCHOR, R_STRIP[1])),
        anchor, R, order, chart['quadrature_points'])
    tau = _integrate(lambda sample: 1.0 / _wall_geometry(np.clip(sample, R_ANCHOR, R_STRIP[1]))['v'],
        anchor, R, order, chart['quadrature_points'])
    U_t = geometry['U_R'] / q0
    R_t = np.full(R.shape, 1.0 / q0)
    E_t, J_t = geometry['E_R'] / q0, geometry['J_R'] / q0
    return dict(radius=R, t=t, U=U, R=R, E=geometry['E'], J=geometry['J'], E_t=E_t, J_t=J_t,
        E_R=geometry['E_R'], J_R=geometry['J_R'], U_t=U_t, R_t=R_t, U_R=geometry['U_R'],
        tau=tau, tau_R=1.0 / geometry['v'], tau_t=1.0 / (q0 * geometry['v']),
        Tcap=Tcap, Tcap_R=geometry['beta_minus'] / (geometry['A'] * geometry['v']),
        Tcap_t=geometry['beta_minus'] / (q0 * geometry['A'] * geometry['v']),
        lambda_w=np.arcsin(R), lambda_w_t=1.0 / (q0 * geometry['root']),
        N_wall=geometry['gamma'] / (q0 * geometry['v']), shift_wall=1.0 / (q0 * geometry['root']),
        gamma=geometry['gamma'], gamma_v=geometry['gamma_v'], s=geometry['s'], v=geometry['v'],
        beta_minus=geometry['beta_minus'], beta_plus=geometry['beta_plus'], A=geometry['A'],
        ew_U=-1.0 / geometry['speed'], ew_R=geometry['J'],
        ew_U_t=(geometry['E_R'] + geometry['J_R']) / (geometry['speed']**2 * q0),
        ew_R_t=geometry['J_R'] / q0)


def _cap_time_density(radius):
    geometry = _wall_geometry(radius)
    return geometry['beta_minus'] / (geometry['A'] * geometry['v'])


def _tail_U(radius):
    return -(radius + 2.0 * M * np.log(radius - 2.0 * M))


def _tail_proper(radius):
    return np.sqrt(radius * (radius - 2.0 * M)) + 2.0 * M * np.arccosh(np.sqrt(radius / (2.0 * M)))


def _coordinate_velocity(N0, N0_y, N0_yy, b0, b0_y, b0_yy, V, V_y, V_yy, rp, rp_y, rp_yy):
    """Initial Xdot=(u_dot, r_dot) and its proper-distance jets. No division by f."""
    c0 = V + rp
    if np.any(c0 <= 0.0):
        raise ValueError('initial outgoing combination V+rprime is not positive')
    c0_y = V_y + rp_y
    c0_yy = V_yy + rp_yy
    num, num_y, num_yy = N0 - b0, N0_y - b0_y, N0_yy - b0_yy
    u = num / c0
    u_y = (num_y * c0 - num * c0_y) / c0**2
    u_yy = (num_yy * c0 - num * c0_yy) / c0**2 - 2.0 * u_y * c0_y / c0
    r = N0 * V + b0 * rp
    r_y = N0_y * V + N0 * V_y + b0_y * rp + b0 * rp_y
    r_yy = N0_yy * V + 2.0 * N0_y * V_y + N0 * V_yy + b0_yy * rp + 2.0 * b0_y * rp_y + b0 * rp_yy
    return u, u_y, u_yy, r, r_y, r_yy


def _lapse_shift(y, chart):
    S, S_y, S_yy = _cutoff(y)
    f1 = y * S
    f1_y = S + y * S_y
    f1_yy = 2.0 * S_y + y * S_yy
    a_plus, s0 = chart['a_plus0'], chart['s0']
    N0 = np.exp(a_plus * f1)
    N0_y = N0 * a_plus * f1_y
    N0_yy = N0 * ((a_plus * f1_y)**2 + a_plus * f1_yy)
    return dict(S=S, S_y=S_y, S_yy=S_yy, f1=f1, f1_y=f1_y, f1_yy=f1_yy,
        N0=N0, N0_y=N0_y, N0_yy=N0_yy, b0=s0 * S, b0_y=s0 * S_y, b0_yy=s0 * S_yy)


def _radius_jets(R, V, V_R, V_RR, rp=None):
    # The bridge descends through a negative slope; the outgoing square root is the positive branch.
    if rp is None:
        rp = np.sqrt(1.0 - 2.0 * M / R + V * V)
    rp_y = M / R**2 + V * V_R
    rp_yy = (-2.0 * M / R**3 + V_R**2 + V * V_RR) * rp
    V_y = V_R * rp
    V_yy = V_RR * rp**2 + V_R * rp_y
    return rp, rp_y, rp_yy, V_y, V_yy


def _base_from_geometry(y, U0, R, V, V_R, V_RR, chart, rp=None, **extra):
    rp, rp_y, rp_yy, V_y, V_yy = _radius_jets(R, V, V_R, V_RR, rp)
    profile = _lapse_shift(y, chart)
    u, u_y, u_yy, r, r_y, r_yy = _coordinate_velocity(
        profile['N0'], profile['N0_y'], profile['N0_yy'], profile['b0'], profile['b0_y'],
        profile['b0_yy'], V, V_y, V_yy, rp, rp_y, rp_yy)
    speed = V + rp
    base = dict(y=y, U0=U0, R0=R, U0_y=-1.0 / speed, R0_y=rp,
        U0_yy=(V_y + rp_y) / speed**2, R0_yy=rp_y, V0=V, V0_y=V_y, V0_yy=V_yy,
        V_R=V_R, V_RR=V_RR, Xdot0_U=u, Xdot0_U_y=u_y, Xdot0_U_yy=u_yy,
        Xdot0_R=r, Xdot0_R_y=r_y, Xdot0_R_yy=r_yy)
    base.update(profile)
    base.update(extra)
    return base


def _bridge_U(q, chart, order):
    q_wall = chart['q_wall']
    rmin, slope, E = chart['rmin'], chart['slope'], chart['E0']

    def integrand(sample):
        sample = np.clip(sample, q_wall, chart['q_taper'])
        lam = 2.0 * rmin / slope * np.cosh(sample)**2
        rp = slope * np.tanh(sample)
        return -lam / (E + rp)

    return _integrate(integrand, q_wall, q, order, chart['quadrature_points'])


def initial_bridge_points(q, chart, order=32):
    """Initial slice on q_wall <= q <= q(R=1), including the minimum and horizon."""
    order = _order(order)
    q = _finite(q, 'q')
    if np.any(q < chart['q_wall'] - 1e-12) or np.any(q > chart['q_taper'] + 1e-12):
        raise ValueError('bridge parameter is outside the wall-to-taper q range')
    q = np.clip(q, chart['q_wall'], chart['q_taper'])
    rmin, slope, E = chart['rmin'], chart['slope'], chart['E0']
    R = rmin * np.cosh(q)**2
    V = np.full(q.shape, E)
    zeros = np.zeros(q.shape, dtype=float)
    y = chart['l_minimum'] + rmin / slope * (q + np.sinh(q) * np.cosh(q)) - chart['l_wall']
    U0 = _bridge_U(q, chart, order)
    return _base_from_geometry(y, U0, R, V, zeros, zeros, chart, rp=slope * np.tanh(q), q=q)


def _outgoing_profile(radius, E):
    width = TAPER[1] - TAPER[0]
    flat = np.atleast_1d(radius)
    x = (flat - TAPER[0]) / width
    step, step_y = smooth_step(x)
    step_y = np.array(step_y, dtype=float, copy=True)
    step_yy = np.zeros(flat.shape, dtype=float)
    band = (x <= UNDERFLOW) | (x >= 1.0 - UNDERFLOW)
    step_y[band] = 0.0
    inner = (x > UNDERFLOW) & (x < 1.0 - UNDERFLOW)
    if np.any(inner):
        step_yy[inner] = _step_second(x[inner])
    V = E * (1.0 - step)
    V_R = -E * step_y / width
    V_RR = -E * step_yy / width**2
    tail = flat >= TAPER[1]
    values = (np.where(tail, 0.0, V), np.where(tail, 0.0, V_R), np.where(tail, 0.0, V_RR))
    return tuple(np.reshape(value, np.shape(radius)) for value in values)


def _outgoing_increments(radius, chart, order):
    start, stop = TAPER
    E = chart['E0']
    shape = np.shape(radius)
    flat = np.atleast_1d(radius)

    def densities(sample):
        sample = np.clip(sample, start, stop)
        V, _, _, rp = _outgoing_rp(sample, E)
        return -1.0 / (rp * (V + rp)), 1.0 / rp

    upper = np.minimum(flat, stop)
    dU = np.atleast_1d(_integrate(lambda sample: densities(sample)[0], start, upper, order, chart['quadrature_points']))
    dy = np.atleast_1d(_integrate(lambda sample: densities(sample)[1], start, upper, order, chart['quadrature_points']))
    tail = flat > stop
    if np.any(tail):
        dU = np.array(dU, copy=True)
        dy = np.array(dy, copy=True)
        dU[tail] = dU[tail] + _tail_U(flat[tail]) - _tail_U(stop)
        dy[tail] = dy[tail] + _tail_proper(flat[tail]) - _tail_proper(stop)
    return dU.reshape(shape), dy.reshape(shape)


def _outgoing_rp(radius, E):
    V, V_R, V_RR = _outgoing_profile(radius, E)
    rp = np.sqrt(1.0 - 2.0 * M / radius + V * V)
    return V, V_R, V_RR, rp


def initial_outgoing_points(radius, chart, order=32):
    """Initial outgoing vacuum from R=1 through the taper and the exact AE tail to R=4."""
    order = _order(order)
    R = _finite(radius, 'radius')
    if np.any(R < TAPER[0] - 1e-12) or np.any(R > DISPLAY_RADIUS + 1e-12):
        raise ValueError('outgoing radius must lie in the original taper-to-display interval [1, 4]')
    R = np.clip(R, TAPER[0], DISPLAY_RADIUS)
    E = chart['E0']
    V, V_R, V_RR = _outgoing_profile(R, E)
    dU, dy = _outgoing_increments(R, chart, order)
    return _base_from_geometry(chart['y_taper'] + dy, chart['U0_taper'][order] + dU, R, V, V_R, V_RR, chart, radius=R)


def make_time_chart():
    """Fixed quiet chart context, including the prepared initial wall anchor."""
    anchor = source_anchor_slice(ELL)
    if not math.isclose(anchor['ell'], ELL) or not math.isclose(anchor['R'], R_ANCHOR):
        raise RuntimeError('time chart requires the ell=1 source anchor at R=.45')
    root = math.sqrt(1.0 - R_ANCHOR**2)
    q0 = anchor['beta_minus'] / (anchor['Rdot'] * root)
    chart = dict(ell=ELL, H=H, k_wall=K_WALL, M=M, R_anchor=R_ANCHOR, R_strip=R_STRIP,
        q0=q0, q_wall=anchor['q_wall'], q_minimum=0.0, q_range=(anchor['q_wall'], None),
        taper=TAPER, taper_start=TAPER[0], taper_stop=TAPER[1], display_radius=DISPLAY_RADIUS,
        cutoff_width=CUTOFF, a_plus0=-R_ANCHOR / root - K_WALL / anchor['gamma'],
        s0=anchor['gamma_v'] / anchor['gamma'], E0=anchor['V_plus'], J0=anchor['Rprime_plus'],
        gamma0=anchor['gamma'], gamma_v0=anchor['gamma_v'], v0=anchor['Rdot'],
        beta_minus0=anchor['beta_minus'], beta_plus0=anchor['beta_plus'], A0=1.0 - R_ANCHOR**2,
        rmin=anchor['minimum_radius'], l_wall=anchor['wall_distance'],
        l_minimum=anchor['minimum_distance'], slope=math.sqrt(1.0 + anchor['V_plus']**2),
        angular_k=1, quadrature_points=[0],
        source='quiet de Sitter cap, tension wall and Schwarzschild exterior; scalar amplitude 0')
    chart['q_taper'] = math.acosh(math.sqrt(TAPER[0] / chart['rmin']))
    chart['q_range'] = (chart['q_wall'], chart['q_taper'])
    chart['q_horizon'] = math.acosh(math.sqrt(2.0 * M / chart['rmin']))
    if not chart['q_wall'] < chart['q_horizon'] < chart['q_taper']:
        raise RuntimeError('event horizon is not on the constant-rate bridge; the fixed chart stops')
    chart['y_taper'] = chart['l_minimum'] + chart['rmin'] / chart['slope'] * (
        chart['q_taper'] + math.sinh(chart['q_taper']) * math.cosh(chart['q_taper'])) - chart['l_wall']
    chart['U0_taper'] = {}
    for order in (32, 64):
        chart['U0_taper'][order] = _scalar(_bridge_U(chart['q_taper'], chart, order), 'U0_taper')
    anchor_point = initial_bridge_points(chart['q_wall'], chart, order=32)
    chart['anchor'] = dict(y=0.0, U0=0.0, R0=R_ANCHOR,
        X0_U=0.0, X0_R=R_ANCHOR,
        X0_y_U=_scalar(anchor_point['U0_y'], 'U0_y'), X0_y_R=_scalar(anchor_point['R0_y'], 'R0_y'),
        Xdot0_U=_scalar(anchor_point['Xdot0_U'], 'Xdot0_U'), Xdot0_R=_scalar(anchor_point['Xdot0_R'], 'Xdot0_R'),
        Xdot0_y_U=_scalar(anchor_point['Xdot0_U_y'], 'Xdot0_U_y'),
        Xdot0_y_R=_scalar(anchor_point['Xdot0_R_y'], 'Xdot0_R_y'),
        Xdot0_yy_U=_scalar(anchor_point['Xdot0_U_yy'], 'Xdot0_U_yy'),
        Xdot0_yy_R=_scalar(anchor_point['Xdot0_R_yy'], 'Xdot0_R_yy'),
        N0=_scalar(anchor_point['N0'], 'N0'), b0=_scalar(anchor_point['b0'], 'b0'),
        N0_y=_scalar(anchor_point['N0_y'], 'N0_y'))
    wall = wall_time_data(R_ANCHOR, chart, order=32)
    de_t_U = _scalar(wall['ew_U_t'], 'ew_U_t') - chart['anchor']['Xdot0_y_U']
    de_t_R = _scalar(wall['ew_R_t'], 'ew_R_t') - chart['anchor']['Xdot0_y_R']
    velocity_U = _scalar(wall['U_t'], 'U_t') - chart['anchor']['Xdot0_U']
    velocity_R = _scalar(wall['R_t'], 'R_t') - chart['anchor']['Xdot0_R']
    chart['de_t0'] = (de_t_U, de_t_R)
    if max(abs(de_t_U), abs(de_t_R), abs(velocity_U), abs(velocity_R), abs(_scalar(wall['N_wall'], 'N_wall') - 1.0)) > 1e-9:
        raise RuntimeError('initial wall velocity or lapse normalization failed: de_t0='
            + repr(chart['de_t0']) + ' velocity_gap=' + repr((velocity_U, velocity_R)))
    chart['prepared_initial_anchor'] = True
    return chart


def _embed(base, wall, chart):
    anchor = chart['anchor']
    t = _scalar(wall['t'], 't')
    dXw_U = _scalar(wall['U'], 'U') - t * anchor['Xdot0_U']
    dXw_R = _scalar(wall['R'], 'R') - anchor['X0_R'] - t * anchor['Xdot0_R']
    de_U = _scalar(wall['ew_U'], 'ew_U') - anchor['X0_y_U'] - t * anchor['Xdot0_y_U']
    de_R = _scalar(wall['ew_R'], 'ew_R') - anchor['X0_y_R'] - t * anchor['Xdot0_y_R']
    dXw_t_U = _scalar(wall['U_t'], 'U_t') - anchor['Xdot0_U']
    dXw_t_R = _scalar(wall['R_t'], 'R_t') - anchor['Xdot0_R']
    de_t_U = _scalar(wall['ew_U_t'], 'ew_U_t') - anchor['Xdot0_y_U']
    de_t_R = _scalar(wall['ew_R_t'], 'ew_R_t') - anchor['Xdot0_y_R']
    S, S_y, f1, f1_y = base['S'], base['S_y'], base['f1'], base['f1_y']
    S_yy, f1_yy = base['S_yy'], base['f1_yy']
    U = base['U0'] + t * base['Xdot0_U'] + dXw_U * S + de_U * f1
    R = base['R0'] + t * base['Xdot0_R'] + dXw_R * S + de_R * f1
    U_y = base['U0_y'] + t * base['Xdot0_U_y'] + dXw_U * S_y + de_U * f1_y
    R_y = base['R0_y'] + t * base['Xdot0_R_y'] + dXw_R * S_y + de_R * f1_y
    U_yy = base['U0_yy'] + t * base['Xdot0_U_yy'] + dXw_U * S_yy + de_U * f1_yy
    R_yy = base['R0_yy'] + t * base['Xdot0_R_yy'] + dXw_R * S_yy + de_R * f1_yy
    U_t = base['Xdot0_U'] + dXw_t_U * S + de_t_U * f1
    R_t = base['Xdot0_R'] + dXw_t_R * S + de_t_R * f1
    U_ty = base['Xdot0_U_y'] + dXw_t_U * S_y + de_t_U * f1_y
    R_ty = base['Xdot0_R_y'] + dXw_t_R * S_y + de_t_R * f1_y
    return dict(U=U, R=R, U_y=U_y, R_y=R_y, U_yy=U_yy, R_yy=R_yy, U_t=U_t, R_t=R_t, U_ty=U_ty, R_ty=R_ty)


def chart_fields(radius, base, chart, order=32):
    """Exact pullback metric, constraint jets and shifted radial Dirac coefficients.

    Signed lapse and g_yy are rejected when nonpositive. Absolute value is not
    used to hide a caustic or an orientation reversal.
    """
    order = _order(order)
    Rwall = _one_wall_radius(radius, 'chart_fields radius')
    for key in ('y', 'U0', 'R0', 'U0_y', 'R0_y', 'U0_yy', 'R0_yy', 'S', 'S_y', 'S_yy',
                'f1', 'f1_y', 'f1_yy', 'Xdot0_U', 'Xdot0_R', 'Xdot0_U_y', 'Xdot0_R_y',
                'Xdot0_U_yy', 'Xdot0_R_yy'):
        if key not in base:
            raise ValueError('base is missing prepared initial jet ' + key)
    wall = wall_time_data(Rwall, chart, order=order)
    embedded = _embed(base, wall, chart)
    R = embedded['R']
    f_safe = np.isfinite(R) & (R > 0.0)
    f = np.where(f_safe, 1.0 - 2.0 * M / np.where(f_safe, R, 1.0), np.nan)
    f_R = np.where(f_safe, 2.0 * M / np.where(f_safe, R, 1.0)**2, np.nan)
    U_y, R_y = embedded['U_y'], embedded['R_y']
    U_t, R_t = embedded['U_t'], embedded['R_t']
    U_yy, R_yy = embedded['U_yy'], embedded['R_yy']
    U_ty, R_ty = embedded['U_ty'], embedded['R_ty']
    h = -f * U_y**2 - 2.0 * U_y * R_y
    D = U_t * R_y - R_t * U_y
    if (not np.all(f_safe) or not np.all(np.isfinite(h)) or not np.all(np.isfinite(D))
            or np.any(h <= 0.0) or np.any(D <= 0.0)):
        y = np.asarray(base['y'], dtype=float)
        finite_h = np.where(np.isfinite(h), h, np.inf)
        finite_D = np.where(np.isfinite(D), D, np.inf)
        iH, iD = int(np.argmin(finite_h)), int(np.argmin(finite_D))
        details = dict(radius=float(Rwall), min_spatial_norm=float(np.min(finite_h)),
            min_jacobian=float(np.min(finite_D)), y_at_min_spatial_norm=float(y.reshape(-1)[iH]),
            y_at_min_jacobian=float(y.reshape(-1)[iD]),
            areal_radius_at_min_jacobian=float(np.asarray(R, dtype=float).reshape(-1)[iD]))
        raise ChartDomainError('quiet time chart leaves the spacelike future domain at Rwall='
            + repr(Rwall) + ': min spatial_norm=' + repr(details['min_spatial_norm'])
            + ' at y=' + repr(details['y_at_min_spatial_norm'])
            + ', min jacobian=' + repr(details['min_jacobian'])
            + ' at y=' + repr(details['y_at_min_jacobian']), details)
    h_y = -f_R * R_y * U_y**2 - 2.0 * f * U_y * U_yy - 2.0 * (U_yy * R_y + U_y * R_yy)
    h_t = -f_R * R_t * U_y**2 - 2.0 * f * U_y * U_ty - 2.0 * (U_ty * R_y + U_y * R_ty)
    a = np.sqrt(h)
    a_y, a_t = h_y / (2.0 * a), h_t / (2.0 * a)
    B = -f * U_t * U_y - U_t * R_y - R_t * U_y
    B_y = (-f_R * R_y * U_t * U_y - f * (U_ty * U_y + U_t * U_yy)
        - (U_ty * R_y + U_t * R_yy) - (R_ty * U_y + R_t * U_yy))
    D_y = U_ty * R_y + U_t * R_yy - R_ty * U_y - R_t * U_yy
    shift = B / h
    shift_y = (B_y * h - B * h_y) / h**2
    c = D / h
    c_y = (D_y * h - D * h_y) / h**2
    N = D / a
    N_y = (D_y * a - D * a_y) / a**2
    if np.any(N <= 0.0):
        raise ValueError('signed lapse is nonpositive at Rwall=' + repr(Rwall))
    Kll = -(a_t / a - shift_y - shift * a_y / a) / N
    normal_V = (R_t - shift * R_y) / N
    Ktheta = -normal_V / R
    rprime = R_y / a
    rsecond = R_yy / a**2 - R_y * a_y / a**3
    scalar = -4.0 * rsecond / R + 2.0 * (1.0 - rprime**2) / R**2
    normal_V_y = (R_ty - shift_y * R_y - shift * R_yy) / N - normal_V * N_y / N
    Ktheta_y = -normal_V_y / R + normal_V * R_y / R**2
    momentum = 2.0 * rprime / R * (Kll - Ktheta) - 2.0 * Ktheta_y / a
    gtt_y = -f_R * R_y * U_t**2 - 2.0 * f * U_t * U_ty - 2.0 * (U_ty * R_t + U_t * R_ty)
    dirac = radial_time_coefficients(lapse=N, lapse_l=N_y, radial_length=a, radial_length_l=a_y,
        radius=R, angular_k=chart['angular_k'])
    dirac['shift'] = shift
    dirac['shift_companion'] = 0.5 * shift_y
    return dict(y=base['y'], U=embedded['U'], R=R, N=N, a=a, shift=shift, c=c, c_y=c_y, shift_y=shift_y,
        N_y=N_y, a_y=a_y, Kll=Kll, Ktheta=Ktheta, mass=R / 2.0 * (1.0 + normal_V**2 - rprime**2),
        hamiltonian_residual=scalar + 4.0 * Kll * Ktheta + 2.0 * Ktheta**2, momentum_residual=momentum,
        jacobian=D, spatial_norm=h, normal_V=normal_V, rprime=rprime, gtt_y=gtt_y,
        U_t=U_t, R_t=R_t, U_y=U_y, R_y=R_y, U_ty=U_ty, R_ty=R_ty, U_yy=U_yy, R_yy=R_yy, dirac=dirac)


def cap_chart_fields(radius, y, chart):
    """Cap metric on -lambda_w <= y <= 0. The angular coefficient is omitted at r=0.

    y>0 would be a smooth continuation across the material wall and is rejected.
    """
    Rwall = _one_wall_radius(radius, 'cap radius')
    y = _finite(y, 'y')
    geometry = _wall_geometry(Rwall)
    q0 = chart['q0']
    lambda_w = math.asin(Rwall)
    if np.any(y < -lambda_w - 1e-12) or np.any(y > 1e-12):
        raise ValueError('cap chart is only defined for -lambda_w <= y <= 0 inside the material wall')
    y = np.clip(y, -lambda_w, 0.0)
    T_t = _scalar(geometry['beta_minus'] / (q0 * geometry['A'] * geometry['v']), 'Tcap_t')
    shift = _scalar(1.0 / (q0 * geometry['root']), 'shift')
    lam = y + lambda_w
    r = np.sin(lam)
    r_y = np.cos(lam)
    N = np.cos(lam) * T_t
    N_y = -np.sin(lam) * T_t
    zeros = np.zeros(y.shape, dtype=float)
    ones = np.ones(y.shape, dtype=float)
    fields = dict(y=y, radius=np.full(y.shape, Rwall), lambda_coordinate=lam, R=r, N=N, a=ones,
        shift=np.full(y.shape, shift), c=N, c_y=N_y, shift_y=zeros, N_y=N_y, a_y=zeros,
        R_y=r_y, R_yy=-r, Kll=zeros, Ktheta=zeros, normal_V=zeros, gtt=-N**2 + shift**2,
        gtt_y=-2.0 * N * N_y, spatial_norm=ones, jacobian=N, mass=r**3 / 2.0,
        angular_coefficient_evaluated=bool(np.all(r > 0.0)))
    if fields['angular_coefficient_evaluated']:
        fields['dirac'] = radial_time_coefficients(lapse=N, lapse_l=N_y, radial_length=ones,
            radial_length_l=zeros, radius=r, angular_k=chart['angular_k'])
        fields['dirac']['shift'] = fields['shift']
        fields['dirac']['shift_companion'] = zeros
    return fields


def fixed_wall_transmission(radius, chart, trace, order=32):
    """Fixed-wall metric match, temporal Israel jump and parallel-frame spinor transmission.

    chi=sqrt(a)*R*psi. No surface spinor force is added. The result is a domain
    identity of this chart, not an eta invariant or a selected state.
    """
    psi = np.asarray(trace, dtype=complex)
    if psi.shape != (2,) or not np.all(np.isfinite(psi)):
        raise ValueError('trace must be a finite complex radial two-column')
    order = _order(order)
    base = initial_bridge_points(chart['q_wall'], chart, order=order)
    exterior = chart_fields(radius, base, chart, order=order)
    cap = cap_chart_fields(radius, 0.0, chart)
    wall = wall_time_data(radius, chart, order=order)
    N = _scalar(exterior['N'], 'N')
    a = _scalar(exterior['a'], 'a')
    shift = _scalar(exterior['shift'], 'shift')
    R = _scalar(exterior['R'], 'R')
    gamma = _scalar(wall['gamma'], 'gamma')
    gamma_v = _scalar(wall['gamma_v'], 'gamma_v')
    c_y = _scalar(exterior['c_y'], 'c_y') - _scalar(cap['c_y'], 'cap c_y')
    shift_y = _scalar(exterior['shift_y'], 'shift_y') - _scalar(cap['shift_y'], 'cap shift_y')
    a_y = _scalar(exterior['a_y'], 'a_y') - _scalar(cap['a_y'], 'cap a_y')
    R_y = _scalar(exterior['R_y'], 'R_y') - _scalar(cap['R_y'], 'cap R_y')
    N_y = _scalar(exterior['N_y'], 'N_y') - _scalar(cap['N_y'], 'cap N_y')
    Kll = _scalar(exterior['Kll'], 'Kll') - _scalar(cap['Kll'], 'cap Kll')
    B = _scalar(exterior['c'], 'c') * ALPHA - shift * _EYE
    B_y = c_y * ALPHA - shift_y * _EYE
    chi = math.sqrt(a) * R * psi
    try:
        chi_y = -0.5 * np.linalg.solve(B, B_y @ chi)
    except np.linalg.LinAlgError as error:
        raise ValueError('fixed-wall principal matrix is singular') from error
    psi_y = chi_y / R - (0.5 * a_y + R_y / R) * psi
    omega = gamma_v * (N_y / N) - gamma * Kll
    parallel = psi_y / gamma + 0.5 * omega * (ALPHA @ psi)
    gtt_jump = _scalar(exterior['gtt_y'], 'gtt_y') - _scalar(cap['gtt_y'], 'cap gtt_y')
    Ktau = -gamma / (2.0 * N**2) * gtt_jump
    beta_plus = gamma_v * _scalar(exterior['normal_V'], 'normal_V') + gamma * _scalar(exterior['rprime'], 'rprime')
    scale = max(1.0, float(np.max(np.abs(psi))))
    n_wall = _scalar(wall['N_wall'], 'N_wall')
    shift_wall = _scalar(wall['shift_wall'], 'shift_wall')
    radius_wall = _scalar(wall['radius'], 'radius')
    return dict(radius=radius_wall,
        radial_length_error=abs(a - 1.0),
        lapse_error=abs(N - n_wall),
        shift_error=abs(shift - shift_wall),
        cap_lapse_error=abs(_scalar(cap['N'], 'cap N') - n_wall),
        cap_shift_error=abs(_scalar(cap['shift'], 'cap shift') - shift_wall),
        cap_radius_error=abs(_scalar(cap['R'], 'cap R') - radius_wall),
        tangent_U_error=abs(_scalar(exterior['U_y'], 'U_y') - _scalar(wall['ew_U'], 'ew_U')),
        tangent_R_error=abs(_scalar(exterior['R_y'], 'R_y') - _scalar(wall['ew_R'], 'ew_R')),
        position_U_error=abs(_scalar(exterior['U'], 'U') - _scalar(wall['U'], 'U')),
        position_R_error=abs(R - radius_wall),
        israel_Ktau_error=abs(Ktau + K_WALL),
        beta_plus_error=abs(beta_plus - _scalar(wall['beta_plus'], 'beta_plus')),
        spinor_parallel_error=float(np.max(np.abs(parallel - 1.5 * K_WALL * psi))) / scale,
        exterior_Kll=_scalar(exterior['Kll'], 'Kll'),
        exterior_a_y=_scalar(exterior['a_y'], 'a_y'),
        cap_Kll=0.0, cap_a_y=0.0, surface_spinor_force=0.0, nested_Psi_over_Phi=None)


def _max_abs(value):
    return float(np.max(np.abs(value)))


def _metric_pullback_error(fields):
    f = 1.0 - 2.0 * M / fields['R']
    gtt = -f * fields['U_t']**2 - 2.0 * fields['U_t'] * fields['R_t']
    gty = -f * fields['U_t'] * fields['U_y'] - fields['U_t'] * fields['R_y'] - fields['R_t'] * fields['U_y']
    gyy = -f * fields['U_y']**2 - 2.0 * fields['U_y'] * fields['R_y']
    adm_tt = -fields['N']**2 + fields['a']**2 * fields['shift']**2
    return max(_max_abs(gtt - adm_tt), _max_abs(gty - fields['a']**2 * fields['shift']),
        _max_abs(gyy - fields['spatial_norm']))


def reference_chart_report():
    """One compact sampled check. Positivity on these nodes is not a global enclosure."""
    started = time.process_time()
    chart = make_time_chart()
    spent_setup = chart['quadrature_points'][0]
    radii = np.linspace(R_STRIP[0], R_STRIP[1], 9)
    bridge_q = np.concatenate((np.linspace(chart['q_wall'], chart['q_taper'], 65),
        [chart['q_minimum'], chart['q_horizon']]))
    outgoing_R = np.concatenate((np.linspace(TAPER[0], TAPER[1], 65), [1.3, 1.5, 2.0, 3.0, 4.0]))
    bridge = initial_bridge_points(bridge_q, chart, order=32)
    outgoing = initial_outgoing_points(outgoing_R, chart, order=32)
    bridge_U_32 = _scalar(initial_bridge_points(chart['q_taper'], chart, order=32)['U0'], 'bridge U')
    bridge_U_64 = _scalar(initial_bridge_points(chart['q_taper'], chart, order=64)['U0'], 'bridge U')
    out32 = initial_outgoing_points(DISPLAY_RADIUS, chart, order=32)
    out64 = initial_outgoing_points(DISPLAY_RADIUS, chart, order=64)
    wall32 = wall_time_data(R_STRIP[1], chart, order=32)
    wall64 = wall_time_data(R_STRIP[1], chart, order=64)
    trace = np.array([0.2 + 0.1j, -0.4 + 0.3j])
    checks = dict(de_t0_abs=max(abs(chart['de_t0'][0]), abs(chart['de_t0'][1])),
        bulk_mass_max_error=0.0, bulk_hamiltonian_max_abs=0.0, bulk_momentum_max_abs=0.0,
        metric_pullback_max_error=0.0, initial_lapse_max_error=0.0, initial_shift_max_error=0.0,
        initial_radial_length_max_error=0.0, initial_Kll_max_error=0.0, initial_velocity_max_error=0.0,
        wall_match_max_error=0.0, israel_Ktau_max_error=0.0, beta_plus_max_error=0.0,
        spinor_parallel_max_error=0.0)
    minima = dict(lapse=math.inf, spatial_norm=math.inf, jacobian=math.inf)
    failure = None
    successful = []
    for radius in radii:
        completed = True
        for name, base in (('bridge', bridge), ('outgoing', outgoing)):
            try:
                fields = chart_fields(float(radius), base, chart, order=32)
            except ChartDomainError as error:
                failure = dict(error.details)
                failure['segment'] = name
                failure['gauge_was_changed'] = False
                completed = False
                break
            checks['bulk_mass_max_error'] = max(checks['bulk_mass_max_error'], _max_abs(fields['mass'] - M))
            checks['bulk_hamiltonian_max_abs'] = max(checks['bulk_hamiltonian_max_abs'], _max_abs(fields['hamiltonian_residual']))
            checks['bulk_momentum_max_abs'] = max(checks['bulk_momentum_max_abs'], _max_abs(fields['momentum_residual']))
            checks['metric_pullback_max_error'] = max(checks['metric_pullback_max_error'], _metric_pullback_error(fields))
            minima['lapse'] = min(minima['lapse'], float(np.min(fields['N'])))
            minima['spatial_norm'] = min(minima['spatial_norm'], float(np.min(fields['spatial_norm'])))
            minima['jacobian'] = min(minima['jacobian'], float(np.min(fields['jacobian'])))
            if math.isclose(float(radius), R_ANCHOR):
                checks['initial_lapse_max_error'] = max(checks['initial_lapse_max_error'], _max_abs(fields['N'] - base['N0']))
                checks['initial_shift_max_error'] = max(checks['initial_shift_max_error'], _max_abs(fields['shift'] - base['b0']))
                checks['initial_radial_length_max_error'] = max(checks['initial_radial_length_max_error'], _max_abs(fields['a'] - 1.0))
                checks['initial_Kll_max_error'] = max(checks['initial_Kll_max_error'], _max_abs(fields['Kll'] + base['V_R']))
                if name == 'bridge':
                    checks['initial_velocity_max_error'] = max(checks['initial_velocity_max_error'],
                        _max_abs(fields['U_t'] - base['Xdot0_U']), _max_abs(fields['R_t'] - base['Xdot0_R']))
        if completed:
            successful.append(float(radius))
        junction = fixed_wall_transmission(float(radius), chart, trace, order=32)
        checks['wall_match_max_error'] = max(checks['wall_match_max_error'],
            junction['radial_length_error'], junction['lapse_error'], junction['shift_error'],
            junction['cap_lapse_error'], junction['cap_shift_error'], junction['cap_radius_error'],
            junction['tangent_U_error'], junction['tangent_R_error'], junction['position_U_error'],
            junction['position_R_error'])
        checks['israel_Ktau_max_error'] = max(checks['israel_Ktau_max_error'], junction['israel_Ktau_error'])
        checks['beta_plus_max_error'] = max(checks['beta_plus_max_error'], junction['beta_plus_error'])
        checks['spinor_parallel_max_error'] = max(checks['spinor_parallel_max_error'], junction['spinor_parallel_error'])
        if failure is not None:
            break
    quadrature = dict(
        bridge_U_32_64_abs_difference=abs(bridge_U_32 - bridge_U_64),
        outgoing_U_32_64_abs_difference=abs(_scalar(out32['U0'], 'U0') - _scalar(out64['U0'], 'U0')),
        outgoing_proper_distance_32_64_abs_difference=abs(_scalar(out32['y'], 'y') - _scalar(out64['y'], 'y')),
        wall_U_32_64_abs_difference=abs(_scalar(wall32['U'], 'U') - _scalar(wall64['U'], 'U')),
        wall_Tcap_32_64_abs_difference=abs(_scalar(wall32['Tcap'], 'Tcap') - _scalar(wall64['Tcap'], 'Tcap')),
        wall_tau_32_64_abs_difference=abs(_scalar(wall32['tau'], 'tau') - _scalar(wall64['tau'], 'tau')),
        is_rigorous_error_enclosure=False, compares_pde_refinement=False)
    failed = [name for name, value in checks.items() if value > 1e-8]
    if failed:
        raise RuntimeError('time-chart identity failed: ' + ', '.join(name + '=' + repr(checks[name]) for name in failed))
    cpu = time.process_time() - started
    if cpu > CPU_SECONDS:
        raise RuntimeError('time-chart report exceeded the 120 process-CPU-second reference budget')
    points = chart['quadrature_points'][0]
    return dict(inputs=dict(ell=ELL, H_ell=H, k_wall_ell=K_WALL, M_over_ell=M, scalar_amplitude=0.0,
            Rwall_over_ell=[float(value) for value in radii], taper_R_over_ell=[TAPER[0], TAPER[1]],
            right_display_R_over_ell=DISPLAY_RADIUS, q0=chart['q0'], angular_k=1),
        sources=dict(cap='de Sitter H*ell=1', wall='tension k_wall*ell=2', exterior='Schwarzschild M/ell=.4',
            scalar_amplitude=0.0, beta_plus_zero_uses_singular_division=False, division_by_f=False),
        gauge=dict(time='t=q0*(Rwall-.45), Nwall(t=0)=1',
            coordinate='y is the wall-pinned coordinate equal to initial proper distance minus initial wall distance; it is not future exterior proper distance',
            cutoff='S=1-smooth_step(y/ell) on [0,1]', initial_lapse='exp(a_plus0*y*S)', initial_shift='s0*S',
            matching='S and y*S fix the wall value and spatial tangent', a_plus0=chart['a_plus0'], s0=chart['s0'],
            q_range=[chart['q_wall'], chart['q_taper']], q_minimum=0.0, q_horizon=chart['q_horizon']),
        checks=checks,
        sampled_minima=dict(lapse=minima['lapse'], spatial_norm=minima['spatial_norm'], jacobian=minima['jacobian'],
            on_successful_metric_calls_only=True, is_uniform_interval_enclosure=False),
        numerical_comparison=quadrature,
        domain=dict(status='exited' if failure else 'sampled_positive',
            successful_wall_radii=successful, failure=failure,
            note='Sampled sign change is a stop, not a cue to change the cutoff, lapse or shift.'),
        counts=dict(planned_wall_radius_cuts=9,
            attempted_wall_radius_cuts=len(successful)+int(failure is not None),
            successful_wall_evaluations=len(successful),
            bridge_nodes=int(bridge_q.size), outgoing_nodes=int(outgoing_R.size),
            quadrature_points=int(points), setup_quadrature_points=int(spent_setup),
            P001_evolution_steps=0, Dirac_evolution_steps=0, state_selections=0, loop_polarizations=0),
        quadrature_budget=dict(reference_process_cpu_seconds=CPU_SECONDS, process_cpu_seconds=cpu,
            maximum_reference_quadrature_points=MAX_QUADRATURE_POINTS, points_used=int(points), blas_threads=1),
        nested_Psi_over_Phi=None,
        validity=('The fixed chart exits at the recorded sample; remaining bulk cuts were not evaluated.' if failure
                  else 'Sampled h>0 and signed lapse>0 only; not a certified global enclosure.'),
        scope='Quiet joined coordinate chart and its Dirac coefficient jets. No quantum state, Hadamard condition, polarization, P001 step or evolved spinor.')
