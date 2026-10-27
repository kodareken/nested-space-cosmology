"""Geometry-derived axial current vertices in a finite conditional quantum sector.

The probes have zero background and do not change the sourced classical body.
One regular Dirichlet axial mode and two positive harmonics on each wall are
retained. The computed four-point object is the propagating exchange/pole
contribution of this projection, not the full constrained Einstein functional.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import time

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_gravitational_response import NestedBody, nested_body, stiff_inside


def body_with_child(radius=.25, scale=1.):
    base = nested_body(scale)
    inner = stiff_inside(base.parent, inside_mass=0., radius=radius*scale)
    if inner.radial_potential(inner.radius)[2] <= 0:
        raise ValueError("inner control failed radial stability")
    return NestedBody(inner, base.outer, base.observer)


def _cap_pair(core, radius, lam, ell=2):
    """Regular cap field and x derivative, without an arbitrary clock reset."""
    r = mp.mpf(radius)
    H = mp.sqrt(core.H2)
    w = mp.sqrt(lam)/(mp.mpf(core.clock)*H)
    z = H*H*r*r
    a, b, c = (ell+2-1j*w)/2, (ell+1-1j*w)/2, mp.mpf(ell)+mp.mpf("1.5")
    F = mp.hyp2f1(a, b, c, z)
    value = r**(ell+1)*(1-z)**(-1j*w/2)*F
    deriv = value*((ell+1)/r+1j*w*H*H*r/(1-z)
                    +2*H*H*r*a*b/c*mp.hyp2f1(a+1, b+1, c+1, z)/F)
    return mp.re(value), mp.re(mp.mpf(core.clock)*(1-z)*deriv)


def _cap_phase(body, lam, ell, k):
    with mp.workdps(40):
        def angle(lmb):
            value, p = _cap_pair(body.core, body.inner.radius, lmb, ell)
            return mp.atan2(value, p/mp.mpf(k))
        theta = angle(mp.mpf(lam))
        sensitivity = mp.diff(angle, mp.mpf(lam))
    return np.array([float(theta), 0., float(sensitivity)])


def _jump(state, q, k):
    theta, logamp, sensitivity = state
    s, c = math.sin(theta), math.cos(theta)
    ratio = q/k
    D = s*s+(c+ratio*s)**2
    return np.array([theta+math.atan2(-ratio*s*s, 1+ratio*s*c),
                     logamp+.5*math.log(D), sensitivity/D])


@dataclass
class Mode:
    body: NestedBody
    omega: float
    ell: int
    k: float
    cap_start: np.ndarray
    walls: list
    end: np.ndarray
    parts: list
    rayleigh_upper: float
    rhs_calls: int
    cpu_seconds: float

    def normalized_pair(self, phase_state):
        theta, amplitude, _ = phase_state
        factor = math.exp(amplitude-self.end[1])/math.sqrt(self.k*self.end[2])
        return np.array([factor*math.sin(theta), self.k*factor*math.cos(theta)])

    def wall_pair(self, index, side="plus"):
        return self.normalized_pair(self.walls[index][side])

    def norm_quadrature(self, order=96):
        """Independent field-square integral, including the analytic cap."""
        x, weight = np.polynomial.legendre.leggauss(order)
        R = self.body.inner.radius
        cap_boundary = self.normalized_pair(self.cap_start)
        with mp.workdps(30):
            reference = float(_cap_pair(self.body.core, R, self.omega**2, self.ell)[0])
            r = R*(x+1)/2
            values = np.array([float(_cap_pair(self.body.core, ri, self.omega**2, self.ell)[0])
                               for ri in r])*cap_boundary[0]/reference
        total = float(np.sum(weight*(R/2)*values**2/self.body.core.speed(r)))
        for region, lo, hi, sol in self.parts:
            r = (lo+hi)/2+(hi-lo)*x/2
            state = sol.sol(r)
            factors = np.exp(state[1]-self.end[1])/math.sqrt(self.k*self.end[2])
            values = factors*np.sin(state[0])
            total += float(np.sum(weight*((hi-lo)/2)*values**2/region.speed(r)))
        return total


def first_axial_mode(body, *, ell=2, rtol=1e-10):
    if ell < 2:
        raise ValueError("physical axial multipole ell>=2 required")
    started = time.process_time()
    M = body.exterior.mass
    k = 1/M
    calls = 0

    def shoot(lam, keep=False):
        nonlocal calls
        initial = _cap_phase(body, lam, ell, k)
        state = initial.copy()
        walls, parts = [], []
        for index, (shell, region, stop) in enumerate([
                (body.inner, body.parent, body.outer.radius),
                (body.outer, body.exterior, body.observer)]):
            before = state.copy()
            state = _jump(state, shell.jump_strength, k)
            walls.append(dict(minus=before, plus=state.copy()))
            def rhs(r, y):
                nonlocal calls
                calls += 1
                if calls > 400_000 or time.process_time()-started > 180:
                    raise RuntimeError("mode allowance exhausted")
                theta, _, sensitivity = y
                s, c = math.sin(theta), math.cos(theta)
                U = region.potential(r, ell)
                speed = region.speed(r)
                return [(k*c*c+(lam-U)*s*s/k)/speed,
                        (k*k+U-lam)*s*c/(k*speed),
                        (2*(lam-U-k*k)*s*c*sensitivity/k+s*s/k)/speed]
            sol = solve_ivp(rhs, (shell.radius, stop), state, method="DOP853",
                            rtol=rtol, atol=rtol*.02, dense_output=keep)
            if not sol.success:
                raise RuntimeError(sol.message)
            state = sol.y[:, -1]
            if keep:
                parts.append((region, shell.radius, stop, sol))
        return state, initial, walls, parts

    # Rayleigh trial obeys both endpoint conditions; shell terms are positive.
    R = body.observer
    trial = lambda r: r**(ell+1)*(1-r/R)
    derivative = lambda r: (ell+1)*r**ell-(ell+2)*r**(ell+1)/R
    intervals = [(body.core, 0., body.inner.radius),
                 (body.parent, body.inner.radius, body.outer.radius),
                 (body.exterior, body.outer.radius, R)]
    norm = sum(quad(lambda r: trial(r)**2/g.speed(r), lo, hi,
                    epsabs=1e-10, epsrel=1e-11)[0] for g, lo, hi in intervals)
    form = sum(quad(lambda r: g.speed(r)*derivative(r)**2
                     +g.potential(r, ell)*trial(r)**2/g.speed(r), lo, hi,
                     epsabs=1e-10, epsrel=1e-11)[0] for g, lo, hi in intervals)
    form += sum(s.jump_strength*trial(s.radius)**2 for s in (body.inner, body.outer))
    upper = form/norm
    lo, hi = shoot(0.)[0][0]-math.pi, shoot(upper)[0][0]-math.pi
    if not lo < 0 < hi:
        raise RuntimeError("first-mode Rayleigh/phase bracket failed")
    lam = brentq(lambda x: shoot(x)[0][0]-math.pi, 0., upper,
                 xtol=rtol*.01/M**2, rtol=max(rtol*.01, 1e-14))
    end, initial, walls, parts = shoot(lam, True)
    if end[2] <= 0:
        raise RuntimeError("eigenvalue sensitivity/norm is not positive")
    return Mode(body, math.sqrt(lam), ell, k, initial, walls, end, parts,
                math.sqrt(upper), calls, time.process_time()-started)


def wall_vertex(mode, index, G=1.):
    """Cubic h-current coupling from canonical spectator and CPM master actions."""
    if mode.ell != 2 or index not in (0, 1) or not math.isfinite(G) or G <= 0:
        raise ValueError("this selected angular vertex requires ell=2, wall0/1 and G>0")
    shell = (mode.body.inner, mode.body.outer)[index]
    R = shell.radius
    lapse = shell.outside.clock*math.sqrt(shell.outside.f(R))
    frequency1, frequency2 = lapse*math.sqrt(2)/R, lapse*math.sqrt(6)/R
    field, p = mode.wall_pair(index)
    angular = 3/(2*math.sqrt(math.pi))
    multipole_factor = math.factorial(mode.ell+2)/math.factorial(mode.ell-2)
    cpm_per_canonical = math.sqrt(64*math.pi*G/multipole_factor)
    h0_per_oscillator = cpm_per_canonical*(shell.outside.speed(R)*field+R*p)/2
    current_coupling = angular*h0_per_oscillator/R**2
    pair_numerator = current_coupling*(frequency1-frequency2)
    return dict(radius=R, lapse=lapse, frequencies=[frequency1, frequency2],
                normalized_master=field, normalized_master_x=p,
                h0_per_canonical_oscillator=h0_per_oscillator,
                angular_integral=angular, current_coupling=current_coupling,
                vacuum_pair_numerator=pair_numerator)


def interaction(mode, G=1.):
    vertices = [wall_vertex(mode, i, G) for i in (0, 1)]
    left, right = [sum(v["frequencies"]) for v in vertices]
    total = left+right
    K4 = 1/(total*(left+mode.omega)*(right+mode.omega))
    # H3=q_g sum c_s(P_1 Q_2-P_2 Q_1). Its first wavefunction term is imaginary;
    # the mixed connected term in log(Psi) is minus the two real pair numerators.
    coefficient = -vertices[0]["vacuum_pair_numerator"]*vertices[1]["vacuum_pair_numerator"]*K4
    product = math.prod(w for v in vertices for w in v["frequencies"])
    normalized_correlation = coefficient/(2*math.sqrt(product))
    M = mode.body.exterior.mass
    return dict(gravity_frequency=mode.omega, vertices=vertices, Et=total,
                EL=left+mode.omega, ER=right+mode.omega, kernel=K4,
                log_wavefunction_exchange_coefficient=coefficient,
                coefficient_shape=coefficient*M**4/G,
                selected_exchange_correlation=normalized_correlation,
                correlation_shape=normalized_correlation*M*M/G,
                sector="one retained axial-mode exchange; constraint/contact and other modes excluded")
