"""Shared invariant scalar action on an exact quiet wall-normal collar.

Units ell=1; prescribed H=1, k_wall=2, M=.4. This is a finite-action
control, not a backreacted P001 evolution. Timelike GHY and the scalar
time-end canonical transformation are retained. No source is fitted.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .nsc_atlas_cap_wall_response import reference_coefficients

ETA = .02
MASS = .4
RMIN, RMAX = .45, .55
AMPLITUDE = .01


@dataclass
class Budget:
    started: float
    field_rhs: int = 0
    background_rhs: int = 0

    def check(self):
        if (time.process_time()-self.started > 180 or self.field_rhs > 6000
                or self.background_rhs > 400000):
            raise RuntimeError('shared-collar reference compute budget exhausted')


def positivity_bounds():
    """Conservative bounds from the geodesic/Jacobi equations, before sampling."""
    lo, hi = reference_coefficients(RMIN), reference_coefficients(RMAX)
    cap_r = RMIN*math.cos(ETA)-float(lo['beta_minus'])*math.sin(ETA)
    ext_r = RMIN+min(float(hi['beta_plus']), 0.)*ETA
    ext_N = 1+float(lo['beta_plus_x'])*ETA-MASS*ETA**2/ext_r**3
    return dict(radius_lower=min(cap_r, ext_r), lapse_lower=min(math.cos(ETA), ext_N),
                cap_radius_lower=cap_r, exterior_radius_lower=ext_r)


class QuietCollar:
    """Exact normal-geodesic solution; exterior integration is only spatial.

    r_R=N, N_R=w, r_tau=v*N. The signed r_n is retained through beta_plus=0.
    Cap and exterior share wall proper time, but the radius has its real kink.
    """
    def __init__(self, radius, budget=None):
        self.radius = float(radius)
        self.wall = {k: float(v) for k, v in reference_coefficients(radius).items()}
        self.v = self.wall['v']
        self.v_R = self.wall['radial_acceleration']/self.v

        def rhs(n, y):
            if budget is not None:
                budget.background_rhs += 1
                budget.check()
            r, rn, N, Nn, w, wn = y
            if r <= 0:
                raise ValueError('normal collar reached the areal-radius singularity')
            return rn, MASS/r**2, Nn, -2*MASS*N/r**3, wn, (
                -2*MASS*w/r**3+6*MASS*N*N/r**4)

        solution = solve_ivp(rhs, (0., ETA), [self.radius, self.wall['beta_plus'],
            1., self.wall['beta_plus_x'], 0., self.wall['beta_second']],
            method='DOP853', rtol=2e-12, atol=2e-14, dense_output=True)
        if not solution.success:
            raise RuntimeError(solution.message)
        self.exterior = solution.sol

    def at(self, normal):
        n = np.atleast_1d(np.asarray(normal, dtype=float))
        if n.ndim != 1 or not np.all(np.isfinite(n)) or np.any(np.abs(n) > ETA):
            raise ValueError('finite normal coordinates in the fixed |n|<=.02 collar required')
        y = np.empty((6, n.size))
        inside = n <= 0
        c, s = np.cos(n[inside]), np.sin(n[inside])
        R, b, bR, bRR = (self.radius, self.wall['beta_minus'],
                        self.wall['beta_minus_x'], self.wall['beta_second'])
        y[:, inside] = np.array([R*c+b*s, -R*s+b*c, c+bR*s,
                                -s+bR*c, bRR*s, bRR*c])
        if np.any(~inside):
            y[:, ~inside] = self.exterior(n[~inside])
        if np.min(y[0]) <= 0 or np.min(y[2]) <= 0:
            raise ValueError('normal collar lost positive radius or future lapse')
        return y

    def wall_sides(self):
        return np.array([[self.radius, self.wall['beta_minus'], 1.,
            self.wall['beta_minus_x'], 0., self.wall['beta_second']],
            [self.radius, self.wall['beta_plus'], 1.,
            self.wall['beta_plus_x'], 0., self.wall['beta_second']]]).T


def hermite_basis(normal):
    """q=(chi_left,eta*chi_n_left,chi_wall,eta*chi_n_wall,chi_right,eta*chi_n_right)."""
    n = np.atleast_1d(np.asarray(normal, dtype=float))
    if n.ndim != 1 or np.any(np.abs(n) > ETA) or not np.all(np.isfinite(n)):
        raise ValueError('basis coordinates must lie in the fixed collar')
    B, D = np.zeros((n.size, 6)), np.zeros((n.size, 6))
    for left, indices in ((True, (0, 1, 2, 3)), (False, (2, 3, 4, 5))):
        mask = n <= 0 if left else n > 0
        t = (n[mask]+ETA)/ETA if left else n[mask]/ETA
        B[np.ix_(mask, indices)] = np.column_stack((2*t**3-3*t*t+1,
            t**3-2*t*t+t, -2*t**3+3*t*t, t**3-t*t))
        D[np.ix_(mask, indices)] = np.column_stack((6*t*t-6*t,
            3*t*t-4*t+1, -6*t*t+6*t, 3*t*t-2*t))/ETA
    return B, D


def boundary_coefficients(profile):
    """Scalar GHY stiffness and its Rwall derivative; outward signs (-,+)."""
    r, rn, N, Nn, w, wn = profile
    eps = np.array([-1., 1.])
    value = eps*(Nn-N*rn/r)/3
    derivative = eps*(wn-(w*rn+N*Nn)/r+N*N*rn/r**2)/3
    return value, derivative


@dataclass
class Assembly:
    geometry: QuietCollar
    n: np.ndarray
    weights: np.ndarray
    B: np.ndarray
    D: np.ndarray
    profile: np.ndarray
    endpoint_profile: np.ndarray
    boundary: np.ndarray
    M: np.ndarray
    K: np.ndarray
    S: np.ndarray
    MR: np.ndarray
    KR: np.ndarray
    SR: np.ndarray

    def velocity(self, q, momentum):
        return np.linalg.solve(self.M, momentum-self.S @ q)

    def energy(self, q, velocity):
        return float((velocity @ self.M @ velocity+q @ self.K @ q)/2)


@lru_cache(maxsize=4)
def _fixed_quadrature(order):
    xi, wi = leggauss(order)
    n = np.concatenate((ETA*(xi-1)/2, ETA*(xi+1)/2))
    weights = np.tile(wi*ETA/2, 2)
    B, D = hermite_basis(n)
    for array in (n, weights, B, D):
        array.setflags(write=False)
    return n, weights, B, D


def assemble(radius, order=32, budget=None):
    if not isinstance(order, int) or not 4 <= order <= 64:
        raise ValueError('Gauss order must be an integer from 4 through 64')
    geometry = QuietCollar(radius, budget)
    n, weights, B, D = _fixed_quadrature(order)
    profile = geometry.at(n)
    r, rn, N, Nn, w, wn = profile

    def mass_matrix(coefficient, basis=B):
        return basis.T @ ((weights*coefficient)[:, None]*basis)

    exterior = n > 0
    V = np.where(exterior, 2*MASS*N/r**3, 0.)
    VR = np.where(exterior, 2*MASS*(w/r**3-3*N*N/r**4), 0.)
    M = mass_matrix(1/N)
    K = mass_matrix(N, D)+mass_matrix(V)
    MR = mass_matrix(-w/N**2)
    KR = mass_matrix(w, D)+mass_matrix(VR)
    S = -geometry.v*mass_matrix(1/r)/3
    SR = -geometry.v_R*mass_matrix(1/r)/3+geometry.v*mass_matrix(N/r**2)/3
    endpoint = geometry.at([-ETA, ETA])
    boundary, boundary_R = boundary_coefficients(endpoint)
    K[(0, 4), (0, 4)] += boundary
    KR[(0, 4), (0, 4)] += boundary_R
    return Assembly(geometry, n, weights, B, D, profile, endpoint, boundary,
                    M, K, S, MR, KR, SR)


def manufactured(radius, normal):
    """Independent continuum field/derivatives, not an assembled-matrix source."""
    n = np.asarray(normal, dtype=float)
    angle = math.pi*(radius-RMIN)/(RMAX-RMIN)
    A = math.sin(angle)**2
    AR = math.pi*math.sin(2*angle)/(RMAX-RMIN)
    ARR = 2*(math.pi/(RMAX-RMIN))**2*math.cos(2*angle)
    g = reference_coefficients(radius)
    v, acceleration = float(g['v']), float(g['radial_acceleration'])
    P, Pn, Pnn = 1+.3*n+.2*n*n+.1*n**3, .3+.4*n+.3*n*n, .4+.6*n
    return dict(chi=AMPLITUDE*A*P, normal=AMPLITUDE*A*Pn,
        normal_second=AMPLITUDE*A*Pnn, time=AMPLITUDE*v*AR*P,
        time_second=AMPLITUDE*(acceleration*AR+v*v*ARR)*P,
        normal_time=AMPLITUDE*v*AR*Pn,
        normal_time_second=AMPLITUDE*(acceleration*AR+v*v*ARR)*Pn)


def manufactured_coefficients(radius):
    field = manufactured(radius, np.array([-ETA, 0., ETA]))
    q, u, acceleration = np.empty(6), np.empty(6), np.empty(6)
    q[::2], q[1::2] = field['chi'], ETA*field['normal']
    u[::2], u[1::2] = field['time'], ETA*field['normal_time']
    acceleration[::2] = field['time_second']
    acceleration[1::2] = ETA*field['normal_time_second']
    return q, u, acceleration


def continuum_force(action):
    g = action.geometry
    f = manufactured(g.radius, action.n)
    r, rn, N, Nn, w, wn = action.profile
    V = np.where(action.n > 0, 2*MASS*N/r**3, 0.)
    J = (f['time_second']/N-g.v*w*f['time']/N**2-Nn*f['normal']
         -N*f['normal_second']+V*f['chi'])
    force = action.B.T @ (action.weights*J)
    end = manufactured(g.radius, np.array([-ETA, ETA]))
    flux = np.array([-1., 1.])*action.endpoint_profile[2]*end['normal']+action.boundary*end['chi']
    force[[0, 4]] += flux
    return force


def reference_report():
    """One fixed-space manufactured evolution; no P001 production field forcing."""
    budget = Budget(time.process_time())
    bounds = positivity_bounds()
    if min(bounds['radius_lower'], bounds['lapse_lower']) <= 0:
        raise ValueError('fixed collar has no positive pre-sampling bound')
    a16, a32 = assemble(.5, 16, budget), assemble(.5, 32, budget)
    quadrature = {name: float(np.max(np.abs(getattr(a16, name)-getattr(a32, name))))
                  for name in ('M', 'K', 'S', 'MR', 'KR', 'SR')}

    def rhs(radius, state):
        budget.field_rhs += 1
        budget.check()
        a = assemble(radius, budget=budget)
        q, p = state[:6], state[6:12]
        u = a.velocity(q, p)
        force = continuum_force(a)
        qR = u/a.geometry.v
        pR = a.SR @ q+a.S @ qR+(force-a.K @ q)/a.geometry.v
        metric_work = float((-u @ a.MR @ u+q @ a.KR @ q)/2)
        source_work = float(u @ force/a.geometry.v)
        return np.concatenate((qR, pR, [metric_work, source_work]))

    sol = solve_ivp(rhs, (RMIN, RMAX), np.zeros(14), method='DOP853',
        rtol=1e-10, atol=1e-12, max_step=(RMAX-RMIN)/80, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    errors = dict(field=0., proper_time_velocity=0., energy_closure=0.,
                  jordan_transmission=0., internal_field_power=0., internal_GHY=0.)
    scales = dict(field=0., proper_time_velocity=0., energy=0.)
    max_canonical_shift = 0.
    max_radius_kink = 0.
    endpoint = None
    for R in np.linspace(RMIN, RMAX, 33):
        budget.check()
        a = assemble(float(R), budget=budget)
        y = sol.sol(R)
        q, p = y[:6], y[6:12]
        u = a.velocity(q, p)
        qe, ue, _ = manufactured_coefficients(float(R))
        for key, actual, exact in (('field', q, qe), ('proper_time_velocity', u, ue)):
            errors[key] = max(errors[key], float(np.max(np.abs(actual-exact))))
            scales[key] = max(scales[key], float(np.max(np.abs(exact))))
        energy = a.energy(q, u)
        errors['energy_closure'] = max(errors['energy_closure'], abs(energy-y[12]-y[13]))
        scales['energy'] = max(scales['energy'], abs(energy), abs(y[12])+abs(y[13]))
        max_canonical_shift = max(max_canonical_shift, float(np.max(np.abs(a.S @ q))))
        r, rn, N, Nn, w, wn = a.geometry.wall_sides()
        phi = q[2]/r
        normal_phi = q[3]/(ETA*r)-q[2]*rn/r**2
        errors['jordan_transmission'] = max(errors['jordan_transmission'],
                                           abs(float(np.diff(normal_phi)[0]-2*phi[0])))
        powers = N*u[2]*q[3]/ETA
        errors['internal_field_power'] = max(errors['internal_field_power'], abs(float(np.diff(powers)[0])))
        robin = (Nn-N*rn/r)/3
        errors['internal_GHY'] = max(errors['internal_GHY'], abs(float(np.diff(robin)[0])))
        max_radius_kink = max(max_radius_kink, abs(float(np.diff(rn)[0])))
        endpoint = dict(radius=float(R), field_max=float(np.max(np.abs(q))),
                        velocity_max=float(np.max(np.abs(u))), energy=energy,
                        metric_work=float(y[12]), source_and_boundary_work=float(y[13]))
    normalized = dict(field=errors['field']/scales['field'],
        proper_time_velocity=errors['proper_time_velocity']/scales['proper_time_velocity'],
        energy_closure=errors['energy_closure']/scales['energy'])
    status = 'shared_action_manufactured_control_passed'
    if max(normalized.values()) > 1e-6:
        status = 'shared_action_manufactured_control_failed'
    return dict(status=status, domain=dict(wall_radius=[RMIN, RMAX], normal=[-ETA, ETA],
            nodes=3, degrees_of_freedom=6, samples=33, quadrature_order=32),
        source=dict(amplitude=AMPLITUDE, temporal='sin(pi*(R-.45)/.1)^2',
                    spatial='1+.3*n+.2*n^2+.1*n^3', interpretation='manufactured numerical control'),
        pre_sampling_bounds=bounds, quadrature_16_32_absolute_differences=quadrature,
        absolute_errors=errors, normalization_scales=scales, normalized_errors=normalized,
        max_ADM_momentum_shift=max_canonical_shift, max_preserved_radius_kink=max_radius_kink,
        endpoint=endpoint, counts=dict(field_rhs=budget.field_rhs,
            background_rhs=budget.background_rhs, accepted_steps=sol.t.size-1),
        limits='prescribed quiet geometry; not backreacted P001, Einstein mass or a sky prediction')
