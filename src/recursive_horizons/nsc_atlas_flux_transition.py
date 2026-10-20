"""Einstein/four-form membranes: a finite inherited vacuum and shared parent.

This is an explicit alternative action, not a Schur term in the original
conformal-field model. c=hbar=1; kappa=8*pi*G; M is geometric mass.
Normal points from inside to outside; K_ab=+h_a^c nabla_c n_b.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad, solve_ivp


@dataclass(frozen=True)
class FluxLaw:
    kappa: float = 1.
    bare_lambda: float = 0.
    charge: float = .5
    tension: float = .2

    def __post_init__(self):
        if (not all(math.isfinite(x) for x in (
                self.kappa, self.bare_lambda, self.charge, self.tension))
                or min(self.kappa, self.tension) <= 0):
            raise ValueError("finite data, positive gravity and membrane tension required")

    @property
    def k(self):
        return self.kappa*self.tension/2

    def vacuum(self, flux):
        if not math.isfinite(flux):
            raise ValueError("finite flux required")
        return self.bare_lambda+self.kappa*flux*flux/2

    def inherit(self, parent_flux, generations):
        if type(generations) is not int or generations < 0:
            raise ValueError("nonnegative integer transition count required")
        self.vacuum(parent_flux)
        flux = parent_flux-generations*self.charge
        return dict(generation=generations, flux=flux, Lambda=self.vacuum(flux))


@dataclass(frozen=True)
class KottlerRegion:
    mass: float
    Lambda: float

    def __post_init__(self):
        if not all(math.isfinite(x) for x in (self.mass, self.Lambda)):
            raise ValueError("finite geometric mass and vacuum curvature required")

    def f(self, radius):
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError("finite positive radius required")
        return 1-2*self.mass/radius-self.Lambda*radius*radius/3

    def f_prime(self, radius):
        self.f(radius)
        return 2*self.mass/radius**2-2*self.Lambda*radius/3

    def tidal(self, radius):
        """Initial relative accelerations of freely falling, locally resting probes.

        Units are acceleration per separation; a laboratory clock in SI
        measures c^2 times these values. Static orthonormal frame requires f>0.
        """
        if self.f(radius) <= 0:
            raise ValueError("tidal laboratory must lie in the static region")
        radial = 2*self.mass/radius**3+self.Lambda/3
        transverse = -self.mass/radius**3+self.Lambda/3
        return dict(radial=radial, transverse=transverse,
                    trace=radial+2*transverse, Ricci_scalar=4*self.Lambda)


@dataclass(frozen=True)
class ChargedWall:
    law: FluxLaw
    outside_flux: float
    inside: KottlerRegion
    outside: KottlerRegion

    def __post_init__(self):
        for actual, expected in (
            (self.inside.Lambda, self.law.vacuum(self.outside_flux-self.law.charge)),
            (self.outside.Lambda, self.law.vacuum(self.outside_flux)),
        ):
            if not math.isclose(actual, expected, rel_tol=2e-14, abs_tol=2e-14):
                raise ValueError("regional vacuum does not follow the declared flux law")

    @classmethod
    def from_inside(cls, law, outside_flux, inside_mass, radius, velocity=0.):
        """Solve the UNSQUARED ordinary-orientation angular junction for outer M."""
        inner = KottlerRegion(inside_mass, law.vacuum(outside_flux-law.charge))
        if not math.isfinite(velocity) or inner.f(radius) <= 0:
            raise ValueError("reference data require a timelike static collar")
        beta = math.sqrt(velocity*velocity+inner.f(radius))
        if beta-law.k*radius <= 0:
            raise ValueError("ordinary outward orientation is not admissible")
        outer_lambda = law.vacuum(outside_flux)
        mass = (inside_mass+(inner.Lambda-outer_lambda)*radius**3/6
                +law.k*radius**2*beta-law.k**2*radius**3/2)
        outer = KottlerRegion(mass, outer_lambda)
        if outer.f(radius) <= 0:
            raise ValueError("outer reference collar must be static")
        return cls(law, outside_flux, inner, outer)

    def coefficients(self):
        A = (self.outside.mass-self.inside.mass)/self.law.k
        B = (self.outside.Lambda-self.inside.Lambda)/(6*self.law.k)+self.law.k/2
        return A, self.inside.mass+A*B, self.inside.Lambda/3+B*B

    def potential(self, radius):
        self.inside.f(radius)
        A, D, C = self.coefficients()
        return 1-A*A/radius**4-2*D/radius-C*radius*radius

    def acceleration(self, radius):
        self.inside.f(radius)
        A, D, C = self.coefficients()
        return -2*A*A/radius**5-D/radius**2+C*radius

    def geometry(self, radius, velocity):
        fm, fp = self.inside.f(radius), self.outside.f(radius)
        if min(fm, fp) <= 0 or not math.isfinite(velocity):
            raise ValueError("evaluation left the declared static collars")
        bm, bp = math.sqrt(velocity**2+fm), math.sqrt(velocity**2+fp)
        acc = self.acceleration(radius)
        ktm = (acc+self.inside.f_prime(radius)/2)/bm
        ktp = (acc+self.outside.f_prime(radius)/2)/bp
        km, kp = ktm+2*bm/radius, ktp+2*bp/radius
        pressure_drop = self.law.charge*(self.outside_flux-self.law.charge/2)
        return dict(radius=radius, velocity=velocity, acceleration=acc,
                    f_inside=fm, f_outside=fp, beta_inside=bm, beta_outside=bp,
                    K_inside=km, K_outside=kp, vacuum_pressure_drop=pressure_drop,
                    angular_residual=bm-bp-self.law.k*radius,
                    temporal_residual=ktp-ktm+self.law.k,
                    force_residual=self.law.tension*(km+kp)/2-pressure_drop,
                    energy_residual=velocity**2+self.potential(radius))


@dataclass(frozen=True)
class DustFluxParent:
    """Exact flat dust+vacuum FLRW outside a mass-compensated Kottler vacuole."""
    region: KottlerRegion
    initial_boundary: float = 1.2

    def __post_init__(self):
        if (self.region.mass <= 0 or self.region.Lambda <= 0
                or self.region.f(self.initial_boundary) <= 0):
            raise ValueError("positive mass/vacuum and initially static join required")

    def state(self, tau):
        if not math.isfinite(tau) or tau < 0:
            raise ValueError("finite nonnegative parent cosmic time required")
        h = math.sqrt(self.region.Lambda/3)
        dust = 2*self.region.mass/self.initial_boundary**3
        H0 = math.sqrt(h*h+dust)
        x = 1.5*h*tau
        u = math.cosh(x)+H0/h*math.sinh(x)
        a = u**(2/3)
        H = h*(math.sinh(x)+H0/h*math.cosh(x))/u
        R = self.initial_boundary*a
        velocity = H*R
        f = self.region.f(R)
        density = 3*dust/a**3  # kappa*rho_m
        acceleration = (h*h-density/6)*R
        return dict(tau=tau, scale_factor=a, H=H, radius=R, velocity=velocity,
                    f=f, kappa_rho_m=density, acceleration=acceleration,
                    deceleration=(density/6-h*h)/(H*H),
                    dust_tidal_trace=3*acceleration/R,
                    vacuum_from_expansion=H*H+2*acceleration/R,
                    Lambda=self.region.Lambda,
                    compensation_residual=density*R**3/6-self.region.mass,
                    angular_residual=velocity**2+f-1,
                    temporal_residual=acceleration+self.region.f_prime(R)/2,
                    Friedmann_residual=3*H*H-density-self.region.Lambda)

    def static_time(self, tau):
        def rate(t):
            f = self.state(t)["f"]
            if f <= 0:
                raise ValueError("static join chart ended at its horizon")
            return 1/f
        return quad(rate, 0., tau, epsabs=2e-12, epsrel=2e-12)[0]


def reference_geometry():
    law = FluxLaw()
    inner = ChargedWall.from_inside(law, 1., 0., .4)
    outer = ChargedWall.from_inside(law, 1.5, inner.outside.mass, .8)
    assert outer.inside == inner.outside
    parent = DustFluxParent(outer.outside)
    return law, outer, inner, parent


def evolve_reference(*, end=.25, rtol=2e-11, atol=2e-13):
    """Both walls on ONE middle Kottler clock, with the exterior cosmic join.

    y = R_outer,v_outer,tau_outer,T_parent,R_inner,v_inner,tau_inner,T_core,tau_FLRW.
    Velocities v are individual proper-radius derivatives, not common-clock rates.
    No field grids or parameter search. A second solve controls integration error.
    """
    if not 0 < end <= .25 or not 0 < atol <= rtol < 1e-5:
        raise ValueError("admitted finite common-clock interval is 0<T<=.25")
    law, outer, inner, parent = reference_geometry()
    evaluations = 0

    def rhs(t, y):
        nonlocal evaluations
        evaluations += 1
        if evaluations > 12000:
            raise RuntimeError("finite ODE budget exceeded")
        go, gi = outer.geometry(y[0], y[1]), inner.geometry(y[4], y[5])
        if y[4] >= y[0]:
            raise ValueError("nested walls collided")
        qo = go["f_inside"]/go["beta_inside"]
        qi = gi["f_outside"]/gi["beta_outside"]
        To = go["beta_outside"]/go["f_outside"]*qo
        Ti = gi["beta_inside"]/gi["f_inside"]*qi
        p = parent.state(y[8])
        if p["f"] <= 0 or y[0] >= p["radius"]:
            raise ValueError("parent join left its admitted collar")
        return (y[1]*qo, go["acceleration"]*qo, qo, To,
                y[5]*qi, gi["acceleration"]*qi, qi, Ti, p["f"]*To)

    initial = (.8, 0., 0., 0., .4, 0., 0., 0., 0.)
    solution = solve_ivp(rhs, (0., end), initial, method="DOP853",
                         rtol=rtol, atol=atol, dense_output=True, max_step=.025)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution, dict(evaluations=evaluations, common_clock_end=end)
