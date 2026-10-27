"""Finite Einstein/four-form/stiff-scalar shells and axial gravitational response.

G=c=M_outer=1 unless a scale is explicitly supplied. This is a separate
static candidate, not a rescaling of the formed scalar source used by P004.
The axial junction is derived from the general metric conditions of Pani
et al. 0909.0287, with a common global clock retained throughout.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import time

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq


@dataclass(frozen=True)
class Region:
    mass: float
    H2: float
    clock: float = 1.0

    def __post_init__(self):
        if (not all(math.isfinite(x) for x in (self.mass, self.H2, self.clock))
                or self.mass < 0 or self.H2 < 0 or self.clock <= 0):
            raise ValueError("finite nonnegative mass/vacuum and positive clock required")

    def f(self, r):
        return 1-2*self.mass/r-self.H2*r*r

    def fp(self, r):
        return 2*self.mass/r**2-2*self.H2*r

    def fpp(self, r):
        return -4*self.mass/r**3-2*self.H2

    def speed(self, r):
        return self.clock*self.f(r)

    def potential(self, r, ell):
        return self.clock**2*self.f(r)*(ell*(ell+1)/r**2-6*self.mass/r**3)


@dataclass(frozen=True)
class Interface:
    radius: float
    inside: Region
    outside: Region

    def __post_init__(self):
        R = self.radius
        if not math.isfinite(R) or R <= 0 or min(self.inside.f(R), self.outside.f(R)) <= 0:
            raise ValueError("static interface requires positive radius and two positive lapses")
        a = self.inside.clock**2*self.inside.f(R)
        b = self.outside.clock**2*self.outside.f(R)
        if abs(a-b) > 1e-12*max(a, b):
            raise ValueError("interface clocks do not match")

    @property
    def sigma(self):
        R = self.radius
        return (math.sqrt(self.inside.f(R))-math.sqrt(self.outside.f(R)))/(4*math.pi*R)

    @property
    def pressure(self):
        R = self.radius
        def curvature_sum(g):
            return g.fp(R)/(2*math.sqrt(g.f(R)))+math.sqrt(g.f(R))/R
        return (curvature_sum(self.outside)-curvature_sum(self.inside))/(8*math.pi)

    @property
    def jump_strength(self):
        """Psi continuous, [dPsi/dx]=q Psi; dx/dr=1/(A f)."""
        R = self.radius
        return (self.inside.speed(R)-self.outside.speed(R))/R

    def radial_potential(self, r):
        # No exchange with a cosmological vacuum: sigma*r^4 is conserved.
        k = 4*math.pi*self.sigma*(self.radius/r)**4
        w = k*r
        wm, wmm = -3*w/r, 12*w/r**2
        D = self.inside.f(r)-self.outside.f(r)
        Dp = self.inside.fp(r)-self.outside.fp(r)
        Dpp = self.inside.fpp(r)-self.outside.fpp(r)
        beta = D/(2*w)+w/2
        bp = Dp/(2*w)-D*wm/(2*w*w)+wm/2
        bpp = (Dpp/(2*w)-Dp*wm/w**2-D*wmm/(2*w*w)
               +D*wm*wm/w**3+wmm/2)
        return np.array([self.inside.f(r)-beta*beta,
                         self.inside.fp(r)-2*beta*bp,
                         self.inside.fpp(r)-2*(bp*bp+beta*bpp)])

    def record(self):
        R = self.radius
        fm, fp = self.inside.f(R), self.outside.f(R)
        lapse = self.outside.clock*math.sqrt(fp)
        Kminus = np.diag([self.inside.fp(R)/(2*math.sqrt(fm)),
                          math.sqrt(fm)/R, math.sqrt(fm)/R])
        Kplus = np.diag([self.outside.fp(R)/(2*math.sqrt(fp)),
                         math.sqrt(fp)/R, math.sqrt(fp)/R])
        jump = Kplus-Kminus
        stress = -(jump-np.trace(jump)*np.eye(3))/(8*math.pi)
        expected = np.diag([-self.sigma, self.pressure, self.pressure])
        # Contract S^a_b with the mean mixed extrinsic curvature.
        normal_contraction = np.trace(expected @ ((Kplus+Kminus)/2))
        bulk_normal_jump = -3*(self.outside.H2-self.inside.H2)/(8*math.pi)
        V, Vp, Vpp = self.radial_potential(R)
        return dict(radius=R, sigma=self.sigma, pressure=self.pressure,
                    scalar_proper_rate=math.sqrt(2*self.sigma),
                    clock_lapse=lapse, f_minus=fm, f_plus=fp,
                    V=float(V), Vprime=float(Vp), Vsecond=float(Vpp),
                    radial_proper_frequency=math.sqrt(Vpp/2) if Vpp > 0 else None,
                    q=self.jump_strength,
                    induced_clock_gap=self.inside.clock**2*fm-self.outside.clock**2*fp,
                    israel_gap=float(np.max(np.abs(stress-expected))),
                    stiff_gap=self.pressure-self.sigma,
                    force_gap=float(normal_contraction-bulk_normal_jump))


def stiff_inside(outside, *, inside_mass, radius):
    """Solve both static junctions for H_inside with p=sigma, ordinary normals."""
    R, m = float(radius), float(inside_mass)
    if not (R > 0 and 0 <= m < outside.mass and outside.f(R) > 0):
        raise ValueError("ordered positive masses and an untrapped interface required")
    b = math.sqrt(outside.f(R))
    C = (1-outside.mass/R-2*outside.H2*R*R)/b+2*b
    d = 1-3*m/R
    discriminant = C*C+16*d
    if discriminant <= 0:
        raise ValueError("stiff-interface equation has no admitted positive root")
    a = (C+math.sqrt(discriminant))/8
    H2 = (1-2*m/R-a*a)/(R*R)
    if not (H2 >= 0 and a > b):
        raise ValueError("positive vacuum and positive stiff surface density required")
    inside = Region(m, H2, outside.clock*b/a)
    return Interface(R, inside, outside)


@dataclass(frozen=True)
class NestedBody:
    inner: Interface
    outer: Interface
    observer: float = 10.0

    @property
    def core(self):
        return self.inner.inside

    @property
    def parent(self):
        return self.inner.outside

    @property
    def exterior(self):
        return self.outer.outside

    def causal_time(self, *, center=False):
        """Proper receiving-clock two-way travel time, not an echo-peak estimate."""
        a = quad(lambda r: 1/self.exterior.speed(r), self.outer.radius, self.observer,
                 epsabs=2e-12, epsrel=2e-12)[0]
        b = quad(lambda r: 1/self.parent.speed(r), self.inner.radius, self.outer.radius,
                 epsabs=2e-12, epsrel=2e-12)[0]
        c = 0.0
        if center:
            H = math.sqrt(self.core.H2)
            c = math.atanh(H*self.inner.radius)/(self.core.clock*H)
        return 2*math.sqrt(self.exterior.f(self.observer))*(a+b+c)


def nested_body(scale=1.0):
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("positive finite exterior mass scale required")
    exterior = Region(scale, 0.0)
    outer = stiff_inside(exterior, inside_mass=.1*scale, radius=2.5*scale)
    inner = stiff_inside(outer.inside, inside_mass=0.0, radius=.25*scale)
    if min(outer.radial_potential(outer.radius)[2],
           inner.radial_potential(inner.radius)[2]) <= 0:
        raise ValueError("candidate failed its radial stability preflight")
    return NestedBody(inner, outer, 10*scale)


@dataclass
class Budget:
    limit: int = 500_000
    cpu_seconds: float = 180.0
    calls: int = 0
    started: float = 0.0

    def __post_init__(self):
        self.started = time.process_time()

    def charge(self):
        self.calls += 1
        if self.calls > self.limit or time.process_time()-self.started > self.cpu_seconds:
            raise RuntimeError("gravitational response compute allowance exhausted")


def _advance(region, start, stop, state, omega, ell, budget, *, rtol, method="DOP853"):
    if ell < 2 or omega <= 0 or region.f(start) <= 0 or region.f(stop) <= 0:
        raise ValueError("ell>=2, positive real frequency and positive regional f required")
    def rhs(r, y):
        budget.charge()
        speed = region.speed(r)
        return np.array([y[1]/speed,
                         (region.potential(r, ell)-omega*omega)*y[0]/speed])
    sol = solve_ivp(rhs, (start, stop), np.asarray(state, dtype=complex),
                    rtol=rtol, atol=rtol*.02, method=method)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1]


def transfer(region, start, stop, omega, ell, budget, *, rtol=1e-11):
    columns = [_advance(region, start, stop, np.eye(2)[:, j], omega, ell,
                        budget, rtol=rtol) for j in range(2)]
    return np.column_stack(columns)


def interface_transfer(interface):
    return np.array([[1., 0.], [interface.jump_strength, 1.]])


def dirichlet_kernel(matrix):
    """Two-endpoint canonical gravitational action, away from a Dirichlet pole.

    Its boundary tractions are (-p_left,p_right). A unit-determinant transfer
    gives a symmetric kernel, so eliminating a shared endpoint is a Schur step.
    """
    if abs(matrix[0, 1]) < 1e-14:
        raise ValueError("Dirichlet pole: retain the transfer representation")
    return np.array([[matrix[0, 0], -1.], [-1., matrix[1, 1]]])/matrix[0, 1]


def axial_metric_data(region, r, state, omega, ell=2):
    """Metric amplitudes for exp(-i omega t), plus the temporal traction data.

    The canonical field is psi=A*f*h1/r and p=dpsi/dx, dx/dr=1/(A*f).
    This also checks the shell's remaining axial Israel equation directly.
    """
    psi, p = state
    c = region.speed(r)
    h0 = 1j/omega*(c*psi+r*p)
    h1 = r*psi/c
    h0p = 1j/omega*(region.clock*region.fp(r)*psi+2*p
                    +r*(region.potential(r, ell)-omega**2)*psi/c)
    return dict(h0=h0, h1=h1, h0p=h0p)


def analytic_cap_state(core, radius, omega, ell):
    """Independent hypergeometric regular de Sitter solution; no ODE history."""
    if core.mass != 0 or core.H2 <= 0 or core.f(radius) <= 0:
        raise ValueError("regular positive-vacuum static cap required")
    with mp.workdps(50):
        r, H = mp.mpf(radius), mp.sqrt(core.H2)
        z = H*H*r*r
        w = mp.mpf(omega)/(mp.mpf(core.clock)*H)
        aa, bb, cc = (ell+2-1j*w)/2, (ell+1-1j*w)/2, mp.mpf(ell)+mp.mpf("1.5")
        hyper = mp.hyp2f1(aa, bb, cc, z)
        value = r**(ell+1)*(1-z)**(-1j*w/2)*hyper
        derivative = value*((ell+1)/r+1j*w*H*H*r/(1-z)
                            +2*H*H*r*aa*bb/cc
                            *mp.hyp2f1(aa+1, bb+1, cc+1, z)/hyper)
        result = np.array([complex(value), complex(core.speed(radius)*derivative)])
    return result/max(np.abs(result))


def cap_state(body, omega, ell, budget, *, rtol, fraction, method="DOP853"):
    core, R = body.core, body.inner.radius
    r0 = fraction*R
    s = ell+1
    c = (core.H2*s*(s+1)-(omega/core.clock)**2)/(4*s+2)
    slope = s/r0+2*c*r0/(1+c*r0*r0)
    result = _advance(core, r0, R, [1., core.speed(r0)*slope], omega, ell,
                      budget, rtol=rtol, method=method)
    return result/max(np.abs(result))


def horizon_state(body, omega, ell, budget, *, rtol, fraction, method="DOP853"):
    parent = body.parent
    # The retained collar is below the cosmological maximum of f.
    rh = brentq(parent.f, parent.mass, 2.2*parent.mass, xtol=1e-15*parent.mass)
    delta = fraction*rh
    r0 = rh+delta
    f1, f2 = parent.fp(rh), parent.fpp(rh)
    q = omega/(parent.clock*f1)
    s = -1j*q
    W = ell*(ell+1)/rh**2-6*parent.mass/rh**3
    c = ((q*q-s/2)*f2+W)/(f1*(1+2*s))
    state0 = np.array([1., parent.speed(r0)*(s/delta+c/(1+c*delta))])
    initial_flux = float(np.imag(np.conj(state0[0])*state0[1]))
    state = _advance(parent, r0, body.inner.radius, state0, omega, ell, budget,
                     rtol=rtol, method=method)
    return state, initial_flux, rh


def _reflection(state, omega):
    incident = (state[0]-state[1]/(1j*omega))/2
    outgoing = (state[0]+state[1]/(1j*omega))/2
    return outgoing/incident, incident


def response(body, omega, *, ell=2, termination="cap", rtol=1e-11,
             origin_fraction=.001, horizon_fraction=4e-6, analytic=False,
             method="DOP853"):
    budget = Budget()
    if termination == "cap":
        seed = (analytic_cap_state(body.core, body.inner.radius, omega, ell) if analytic
                else cap_state(body, omega, ell, budget, rtol=rtol,
                               fraction=origin_fraction, method=method))
        seed = interface_transfer(body.inner) @ seed
        seed_flux = float(np.imag(np.conj(seed[0])*seed[1]))
        rh = None
    elif termination == "horizon":
        seed, seed_flux, rh = horizon_state(body, omega, ell, budget, rtol=rtol,
                                           fraction=horizon_fraction, method=method)
    else:
        raise ValueError("termination must be cap or horizon")
    state = _advance(body.parent, body.inner.radius, body.outer.radius, seed,
                     omega, ell, budget, rtol=rtol, method=method)
    state = interface_transfer(body.outer) @ state
    state = _advance(body.exterior, body.outer.radius, body.observer, state,
                     omega, ell, budget, rtol=rtol, method=method)
    reflection, incident = _reflection(state, omega)
    flux = float(np.imag(np.conj(state[0])*state[1]))
    absorption = -seed_flux/(omega*abs(incident)**2)
    lapse = math.sqrt(body.exterior.f(body.observer))
    return dict(omega=omega, observer_proper_frequency=omega/lapse, ell=ell,
                termination=termination, reflection=reflection,
                reflected_power=abs(reflection)**2, absorbed_power=absorption,
                energy_balance_gap=1-abs(reflection)**2-absorption,
                flux_relative_gap=abs(flux-seed_flux)/max(abs(seed_flux), 1e-100)
                    if termination == "horizon" else abs(flux),
                impedance=state[1]/state[0], horizon_radius=rh,
                rhs_calls=budget.calls, cpu_seconds=time.process_time()-budget.started)
