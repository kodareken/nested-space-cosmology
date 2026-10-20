"""Exact cosmological attachment and explicitly separate bulk-gravity examples.

G=c=1 in the four-dimensional attachment; lengths are in the original ell.
The five-dimensional examples keep their own L5 and kappa5_squared. They
are imported candidate families, not a change to the NSC shell action.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


@dataclass(frozen=True)
class DustParent:
    """Flat expanding Einstein--Straus parent of a Schwarzschild collar."""

    mass: float = .4
    initial_boundary: float = 4.

    def __post_init__(self):
        if not (math.isfinite(self.mass) and math.isfinite(self.initial_boundary)
                and self.mass > 0 and self.initial_boundary > 2*self.mass):
            raise ValueError("positive mass and an initially untrapped boundary required")

    @property
    def H0(self):
        return math.sqrt(2*self.mass/self.initial_boundary**3)

    def state(self, tau):
        x = 1+1.5*self.H0*tau
        if not math.isfinite(x) or x <= 0:
            raise ValueError("cosmic time must be later than the dust big bang")
        a = x**(2/3)
        R = self.initial_boundary*a
        H = self.H0/x
        v = H*R
        f = 1-2*self.mass/R
        if f <= 0:
            raise ValueError("this static-chart evaluator requires R_boundary>2M")
        rho = 3*H*H/(8*math.pi)
        return dict(tau=tau, scale_factor=a, radius=R, H=H, density=rho,
                    radial_velocity=v, radial_acceleration=-self.mass/R**2,
                    f=f, T_rate=1/f, mass_removed=4*math.pi*rho*R**3/3,
                    K_angular=1/R, K_time=0., deceleration=.5,
                    parent_Lambda=0., Ricci_trace_over_four=2*math.pi*rho)

    def static_time(self, tau):
        return quad(lambda t: self.state(t)["T_rate"], 0., tau,
                    epsabs=2e-12, epsrel=2e-12)[0]

    def tortoise(self, radius):
        if radius <= 2*self.mass:
            raise ValueError("tortoise evaluation requires R>2M")
        return radius+2*self.mass*math.log(radius/(2*self.mass)-1)

    def arrival(self, emitter_radius=2.):
        """Static source at T=T_boundary(0); reception by the moving join.

        This experiment is entirely in the parent-side exterior. It is
        not a flash transmitted from the causally separated nested cap.
        """
        if not 2*self.mass < emitter_radius < self.initial_boundary:
            raise ValueError("emitter must lie between the horizon and initial join")
        distance = self.tortoise(self.initial_boundary)-self.tortoise(emitter_radius)

        def outgoing_phase(t):
            return quad(lambda s: 1/(1+self.state(s)["radial_velocity"]),
                        0., t, epsabs=2e-12, epsrel=2e-12)[0]-distance

        # The expanding boundary slows down, so this upper bound is sufficient.
        upper = distance*(1+self.state(0.)["radial_velocity"])
        tau = brentq(outgoing_phase, 0., upper, xtol=2e-13)
        state = self.state(tau)
        fe = 1-2*self.mass/emitter_radius
        ratio = (1+state["radial_velocity"])/math.sqrt(fe)
        return dict(emitter_radius=emitter_radius, reception_time=tau,
                    reception_radius=state["radius"], one_plus_redshift=ratio,
                    received_over_emitted_frequency=1/ratio,
                    null_equation_residual=self.static_time(tau)
                    -self.tortoise(state["radius"])+self.tortoise(emitter_radius))

    def initial_slice(self, radius, inner=1., outer=3.5):
        """A C2 spacelike splice from the retained static collar to cosmic time.

        Its slope is w(R)*sqrt(2M/R)/f. The quintic w is constant on
        collars; C2 is sufficient for the classical initial constraints.
        No physical surface source is inserted at either splice.
        """
        if not 2*self.mass < inner < outer < self.initial_boundary:
            raise ValueError("slice splice must fit outside the horizon and inside join")
        r = np.asarray(radius, dtype=float)
        if np.any(~np.isfinite(r)) or np.any(r <= 2*self.mass):
            raise ValueError("slice is evaluated on the exterior branch R>2M")
        s = np.clip((r-inner)/(outer-inner), 0., 1.)
        w = s**3*(10+s*(-15+6*s))
        f = 1-2*self.mass/r
        v = np.sqrt(2*self.mass/r)
        slope = w*v/f
        spatial_metric = (1-w*w*v*v)/f
        return dict(weight=w, T_radius=slope, radial_metric=spatial_metric,
                    normal_radius_rate=w*v/np.sqrt(1-w*w*v*v)*np.sqrt(f))


def quadratic_brane_stress(stress):
    """SMS pi_ab in an orthonormal (-+++) frame, for symmetric covariant T."""
    T = np.asarray(stress, dtype=float)
    if T.shape != (4, 4) or not np.all(np.isfinite(T)) or not np.allclose(T, T.T):
        raise ValueError("finite symmetric covariant 4x4 stress required")
    eta = np.diag([-1., 1., 1., 1.])
    tr = np.trace(eta @ T)
    square_trace = np.trace(eta @ T @ eta @ T)
    return (-T @ eta @ T/4+tr*T/12
            +eta*square_trace/8-eta*tr*tr/24)


def projected_brane_source(stress, *, tension, kappa5_squared, L5, weyl=None):
    """Gauss/Israel projection with geometric Lambda5=-6/L5^2.

    The returned Weyl tensor is supplied bulk data; this projection does
    not solve its bulk evolution or establish a bulk completion for T.
    """
    if min(tension, kappa5_squared, L5) <= 0:
        raise ValueError("positive tension, coupling and bulk length required")
    T = np.asarray(stress, dtype=float)
    pi = quadratic_brane_stress(T)
    eta = np.diag([-1., 1., 1., 1.])
    E = np.zeros((4, 4)) if weyl is None else np.asarray(weyl, dtype=float)
    if (E.shape != (4, 4) or not np.all(np.isfinite(E))
            or not np.allclose(E, E.T) or abs(np.trace(eta @ E)) > 1e-12):
        raise ValueError("projected bulk Weyl tensor must be symmetric and traceless")
    k4 = kappa5_squared**2*tension/6
    lam4 = -3/L5**2+kappa5_squared**2*tension**2/12
    K = -kappa5_squared/2*(T+(tension-np.trace(eta @ T))*eta/3)
    Ktr = np.trace(eta @ K)
    extrinsic = Ktr*K-K @ eta @ K-eta/2*(Ktr*Ktr-np.trace(eta @ K @ eta @ K))
    direct = 3/L5**2*eta+extrinsic-E
    effective = -lam4*eta+k4*T+kappa5_squared**2*pi-E
    return dict(kappa4_squared=k4, Lambda4=lam4, extrinsic_curvature=K,
                quadratic_stress=pi, Einstein_direct=direct,
                Einstein_effective=effective,
                projection_residual=float(np.max(np.abs(direct-effective))))


def de_sitter_bulk_slice(y, *, detuning=1.1, L5=1.):
    """Exact pure-AdS5 warp with a positive-tension dS4 slice at y=0.

    detuning = kappa5_squared*tension*L5/6 > 1 is supplied, not selected.
    """
    if not (math.isfinite(detuning) and math.isfinite(L5)
            and detuning > 1 and L5 > 0):
        raise ValueError("this dS branch requires finite detuning>1 and L5>0")
    z = np.asarray(y, dtype=float)/L5
    horizon = math.atanh(1/detuning)
    if np.any(~np.isfinite(z)) or np.any(z < 0) or np.any(z >= horizon):
        raise ValueError("evaluate the positive warp side before its coordinate horizon")
    A = np.cosh(z)-detuning*np.sinh(z)
    Ap = (np.sinh(z)-detuning*np.cosh(z))/L5
    App = A/L5**2
    curvature4 = (detuning**2-1)/L5**2
    return dict(warp=A, derivative=Ap, second_derivative=App,
                H4_squared=curvature4, Lambda4_L5_squared=3*(detuning**2-1),
                coordinate_horizon=horizon*L5,
                tangential_bulk_curvature=(curvature4-Ap*Ap)/(A*A),
                normal_bulk_curvature=-App/A)


def shadow_lensing(shadow_fraction, *, separation_over_L5=5.):
    """Two-brane massless-radion, zero-mode, static weak-field prediction.

    f is the fraction of the *dynamical potential*, with common radial
    profiles for the two components. It is not a cosmological mass fraction.
    The observer is on the positive-tension brane.
    """
    if not (math.isfinite(shadow_fraction) and 0 <= shadow_fraction <= 1
            and math.isfinite(separation_over_L5) and separation_over_L5 > 0):
        raise ValueError("0<=shadow fraction<=1 and finite positive separation required")
    warp_squared = math.exp(-2*separation_over_L5)
    gamma_visible = (3-warp_squared)/(3+warp_squared)
    gamma = (1-shadow_fraction)*gamma_visible+.5*shadow_fraction
    return dict(shadow_fraction=shadow_fraction, gamma_visible=gamma_visible,
                Psi_over_Phi=gamma, lensing_over_GR_same_dynamical_mass=(1+gamma)/2)


def population_balance(a, *, exponent, initial_fraction=.1):
    """Open pressureless BH sector and a conserving pressureless reservoir.

    A transparent counterexample to inferring physical vacuum pressure
    solely from M proportional to a**exponent. Total initial density is 1.
    The reservoir must remain nonnegative on the declared interval.
    """
    if (not all(math.isfinite(x) for x in (a, exponent, initial_fraction))
            or a <= 0 or not 0 < initial_fraction < 1):
        raise ValueError("finite positive scale and 0<initial fraction<1 required")
    bh = initial_fraction*a**(exponent-3)
    total = a**-3
    reservoir = total-bh
    if reservoir < -1e-14:
        raise ValueError("the specified finite reservoir has been exhausted")
    return dict(scale_factor=a, black_hole_density=bh, reservoir_density=reservoir,
                total_density=total, total_pressure=0.,
                transfer_Q_over_H=exponent*bh,
                closed_effective_w_if_separately_conserved=-exponent/3,
                total_deceleration=.5)


def tension_potential(x, *, mu, u):
    """Original constant-tension wall: x=kR, mu=kM, u=H^2/k^2.

    dR/dtau squared + V = 0; derivatives here are with respect to x.
    Signed beta factors still need checking when reconstructing a geometry.
    """
    x = np.asarray(x, dtype=float)
    if (np.any(~np.isfinite(x)) or np.any(x <= 0)
            or not math.isfinite(mu) or mu <= 0
            or not math.isfinite(u) or u < 0):
        raise ValueError("finite x>0, mu>0 and u>=0 required")
    C = (1+u)**2/4
    V = 1-mu**2/x**4-mu*(1-u)/x-C*x*x
    Vx = 4*mu**2/x**5+mu*(1-u)/x**2-2*C*x
    Vxx = -20*mu**2/x**6-2*mu*(1-u)/x**3-2*C
    return dict(V=V, derivative=Vx, second_derivative=Vxx,
                beta_minus=mu/x**2+(1-u)*x/2,
                beta_plus=mu/x**2-(1+u)*x/2)


def critical_tension_wall(u):
    """The unstable static family of the original wall law, without a scan."""
    if not math.isfinite(u) or u < 0:
        raise ValueError("finite u>=0 required")
    q = (u-1+math.sqrt(9*u*u+14*u+9))/8
    x = 1/math.sqrt(3*((1+u)**2/4-q*q))
    mu = q*x**3
    state = tension_potential(x, mu=mu, u=u)
    return dict(u=u, x=x, mu=mu, beta_minus=float(state["beta_minus"]),
                beta_plus=float(state["beta_plus"]), V=float(state["V"]),
                derivative=float(state["derivative"]),
                second_derivative=float(state["second_derivative"]),
                Lambda_Lref_squared=12*u,
                static_is_radially_stable=False)
