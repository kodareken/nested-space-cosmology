"""Radial test-body transport through the prepared scalar/gravity geometry.

The cap profiles include their solved gravity and wall. The receiving open
child uses its vacuum comparison metric, with a separately recorded scalar
energy bound. This is a kinematic response, not scalar production or a
solution of the perturbed Einstein/membrane system.
"""
from dataclasses import dataclass
from functools import cached_property
import math
from typing import Callable

from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq


@dataclass(frozen=True)
class RadialCapTransport:
    child_length: float
    parent_length: float
    child_a: Callable[[float], float]
    parent_a: Callable[[float], float]
    parent_ap: Callable[[float], float]
    child_h: float
    pole_h_squared: float
    quadrature_tolerance: float = 2e-11

    @classmethod
    def from_preparation(cls, model, solution):
        Li, Lo = map(float, solution.p)
        E0 = model.pole_transfer(solution)["scalar_energy_at_pole"]
        return cls(Li, Lo,
                   lambda r: float(solution.sol(r/Li)[0]),
                   lambda r: float(solution.sol(r/Lo)[4]),
                   lambda r: float(solution.sol(r/Lo)[5]),
                   math.sqrt(model.scalar.kappa*model.child_vacuum/3),
                   model.scalar.kappa*(model.child_vacuum+E0)/3)

    def integrate(self, f, left, right):
        return quad(f, left, right, epsabs=self.quadrature_tolerance,
                    epsrel=self.quadrature_tolerance, limit=200)[0]

    @cached_property
    def wall_x(self):
        # x-ln(xi) is regular. The omitted pole remainder is O(eps^4).
        eps = min(1e-4, self.child_length*1e-4)
        correction = self.pole_h_squared*eps*eps/12
        return (math.log(self.child_length)+correction
                +self.integrate(lambda r: (r-self.child_a(r))/(r*self.child_a(r)),
                                eps, self.child_length))

    def source(self, depth, eta=0.):
        if not (math.isfinite(depth) and 0 <= depth < self.parent_length
                and math.isfinite(eta) and 0 <= eta <= 10):
            raise ValueError("source requires 0<=depth<parent length and 0<=eta<=10")
        r = self.parent_length-depth
        a = self.parent_a(r)
        x = self.wall_x+self.integrate(lambda t: 1/self.parent_a(t),
                                      r, self.parent_length)
        # The positive-density caps are concave; the ordinary child cap is
        # increasing and the parent has at most one maximum on this route.
        peak = max(self.child_a(self.child_length), a)
        if self.parent_ap(r) > 0 > self.parent_ap(self.parent_length):
            equator = brentq(self.parent_ap, r, self.parent_length)
            peak = max(peak, self.parent_a(equator))
        return dict(depth=float(depth), eta=float(eta), parent_pole_coordinate=r,
                    a_source=a, x_source=x, maximum_a_on_route=peak)

    def null(self, depth, eta=0.):
        src = self.source(depth, eta)
        log_label = eta+src["x_source"]
        limit = math.log(2/self.child_h)
        row = dict(source=src, log_label=log_label,
                   central_arrival=log_label < limit,
                   latest_emission_eta=limit-src["x_source"])
        if row["central_arrival"]:
            tau = 2/self.child_h*math.atanh(math.exp(log_label-limit))
            b = math.sinh(self.child_h*tau)/self.child_h
            row.update(arrival_time=tau, one_plus_redshift=b/src["a_source"])
        return row

    def massive(self, depth, gamma_source, eta=0.):
        if not math.isfinite(gamma_source) or gamma_source <= 1:
            raise ValueError("strictly inward source requires finite gamma>1")
        src = self.source(depth, eta)
        K = src["a_source"]*gamma_source
        row = dict(source=src, gamma_source=gamma_source, radial_charge=K,
                   central_arrival=False)
        if K <= src["maximum_a_on_route"]:
            row["status"] = "turning_or_separatrix_before_child_cone"
            return row

        def delay(a):
            root = math.sqrt(1-(a/K)**2)
            return a/(K*K*root*(1+root))

        D = (self.integrate(lambda r: delay(self.child_a(r)), 0, self.child_length)
             +self.integrate(lambda r: delay(self.parent_a(r)),
                             src["parent_pole_coordinate"], self.parent_length))
        log_label = eta+src["x_source"]+D
        limit = 2*K/(1+self.child_h*K)
        row.update(delay=D, log_label=log_label, central_label_limit=limit,
                   central_arrival=log_label < math.log(limit))
        if not row["central_arrival"]:
            row["status"] = "enters_child_but_misses_central_observer"
            return row
        tau = de_sitter_arrival(self.child_h, K, log_label)
        b = math.sinh(self.child_h*tau)/self.child_h
        gamma_r = math.hypot(1., K/b)
        row.update(status="central_arrival", arrival_time=tau,
                   child_scale_factor=b, gamma_received=gamma_r,
                   energy_ratio=gamma_r/gamma_source,
                   speed_received=K/math.hypot(b, K))
        return row

    def central_threshold(self, depth, eta=0.):
        if not self.null(depth, eta)["central_arrival"]:
            raise ValueError("source is outside the central null horizon")
        src = self.source(depth, eta)
        lo = max(1., src["maximum_a_on_route"]/src["a_source"])*(1+1e-7)

        def gap(g):
            row = self.massive(depth, g, eta)
            return row["log_label"]-math.log(row["central_label_limit"])

        hi = max(2., 2*lo)
        for _ in range(16):
            if gap(hi) < 0:
                return brentq(gap, lo, hi, xtol=2e-12)
            hi *= 2
        raise RuntimeError("finite threshold bracket budget exceeded")

    def null_collar(self, eta=0.):
        """Parent portion of the eta slice in the central receiver's causal past."""
        target = math.log(2/self.child_h)-eta
        if target <= self.wall_x:
            return dict(depth=0., parent_slice_volume_fraction=0.)
        depth = brentq(lambda d: self.source(d, eta)["x_source"]-target,
                       0., self.parent_length*(1-1e-7), xtol=2e-11)
        # cosh(eta)^2 and 4pi cancel from this proper-volume fraction.
        volume = self.integrate(lambda r: self.parent_a(r)**2, 0, self.parent_length)
        collar = self.integrate(lambda r: self.parent_a(r)**2,
                                self.parent_length-depth, self.parent_length)
        return dict(depth=depth, parent_slice_volume_fraction=collar/volume)


def de_sitter_arrival(h, K, log_label):
    """Invert the regularized radial timelike integral, within its central cone."""
    if min(h, K) <= 0 or not all(map(math.isfinite, (h, K, log_label))):
        raise ValueError("positive finite h,K and finite logarithmic label required")
    if log_label >= math.log(2*K/(1+h*K)):
        raise ValueError("no finite central arrival")
    lam = 2*math.exp(-log_label)/h
    c = 1/(h*K)**2-1
    z = (lam*lam-c)/(2*lam)
    return math.atanh(1/z)/h


def independent_child_arrival(h, K, log_label, *,max_time=12.):
    """Integrate the first-order trajectory, without the closed-form primitive.

    y=chi+ln(tau) removes the regular-cone logarithm. The event chi=0
    determines the reception time independently of analytic inversion.
    """
    start = 1e-5
    y0 = log_label+(h*h/12+1/(4*K*K))*start*start

    def rhs(t, _):
        b = math.sinh(h*t)/h
        return [1/t-K/(b*math.hypot(b, K))]

    def arrival(t, y):
        return y[0]-math.log(t)

    arrival.terminal = True
    arrival.direction = -1
    sol = solve_ivp(rhs, (start, max_time), [y0], events=arrival,
                    rtol=3e-12, atol=3e-13, max_step=.04)
    if not sol.success or not len(sol.t_events[0]):
        raise RuntimeError("no central arrival within independent integration budget")
    return float(sol.t_events[0][0])


def de_sitter_embedding_control(h, source_angle, gamma):
    """Independent geodesic X(s)=X0 cosh(hs)+U0 sinh(hs)/h in R^(1,4)."""
    if not 0 < source_angle < math.pi/2 or gamma <= 1/math.cos(source_angle):
        raise ValueError("embedding source must be in the central timelike cone")
    hs = math.atanh(math.tan(source_angle)/math.sqrt(gamma*gamma-1))
    tau = math.asinh(gamma*math.sinh(hs))/h
    return dict(arrival_time=tau, gamma_received=gamma*math.cos(source_angle),
                threshold=1/math.cos(source_angle))
