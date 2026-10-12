"""Quadratic null-constraint data forced by the existing conformal pulse.

This integrates one initial characteristic, not the field or Einstein evolution.
The affine parameter is the reference radius x=R0/ell, anchored at the corner;
holding its lapse fixed does not also hold the perturbed radius fixed.
"""
from dataclasses import dataclass, field
import math
import time

import numpy as np

from .nsc_atlas_causal_reference import reference_wall_events
from .nsc_atlas_field_preparation import (
    radius_from_null, radius_product, validate_preparation_domain,
    validate_preparation_inputs,
)
from .nsc_atlas_field_transfer import compact_bump

SCHEMA = "nsc-atlas-characteristic-initialization-v1"
MAX_SAMPLES = 129
MAX_EVALUATIONS = 20000
CPU_SECONDS = 60
TOLERANCE_ABSOLUTE = 1e-8
TOLERANCE_RELATIVE = 1e-8


@dataclass
class Budget:
    """One creation owns the cumulative quadrature and CPU budget."""
    maximum_evaluations: int = MAX_EVALUATIONS
    cpu_seconds: float = CPU_SECONDS
    evaluations: int = 0
    started: float = field(default_factory=time.process_time)

    def check(self, count=0):
        self.evaluations += int(count)
        if (self.evaluations > self.maximum_evaluations
                or time.process_time()-self.started > self.cpu_seconds):
            raise RuntimeError("characteristic initialization budget exhausted")


def validate_model(*, xi=1/6, M_over_ell=.4, H_ell=1., k_ell=2.,
                   scalar_kinetic_coefficient=1., scalar_potential=0.,
                   surface_scalar_coupling=0., corner_radius_correction=0.,
                   corner_derivative_correction=0., bottom_scalar=0.,
                   lapse_fixed=True, radius_fixed=False):
    values = (xi, M_over_ell, H_ell, k_ell, scalar_kinetic_coefficient,
              scalar_potential, surface_scalar_coupling, corner_radius_correction,
              corner_derivative_correction, bottom_scalar)
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("model data must be finite")
    validate_preparation_inputs(xi=xi, surface_scalar_coupling=surface_scalar_coupling,
                                source_corner=bottom_scalar, bottom_corner=bottom_scalar)
    if values[:6] != (1/6, .4, 1., 2., 1., 0.):
        raise ValueError("changed theory requires a separately admitted candidate")
    if corner_radius_correction or corner_derivative_correction:
        raise ValueError("both corner corrections are fixed to zero")
    if lapse_fixed is not True or radius_fixed is not False:
        raise ValueError("hold the affine lapse fixed, not both lapse and radius")


def source_geometry(order=64, budget=None):
    """Derive source endpoints from the original geometric adapter each time."""
    if order not in (32, 64):
        raise ValueError("only admitted geometric orders 32 and 64 are available")
    # In the unchanged borrowed adapter: one anchor interval, then 4+3
    # intervals for each of v and tau at the two fixed wall radii.
    if budget is not None:
        budget.check(15*order)
    wall = reference_wall_events(np.array([.45, .55]), order)
    Ua, Ub = (float(v) for v in wall["U"])
    Va, Vb = (float(v) for v in wall["V"])
    U0 = float(radius_product(.4))/Va
    validate_preparation_domain(U0, Ua, Ub, Va, Vb)
    support_V = [Va+(Vb-Va)/4, Va+3*(Vb-Va)/4]
    return {"order": order, "U0": U0, "Ua": Ua, "Ub": Ub, "Va": Va, "Vb": Vb,
            "x_corner": .4, "x_last": float(radius_from_null(U0, Vb)),
            "support_x": radius_from_null(U0, support_V).tolist(),
            "support_V": support_V}


def source_jet(x, geometry, amplitude=1.):
    """Analytic p, dp/dx, d2p/dx2 with phi=(epsilon/ell)p; no saved differences."""
    x = np.asarray(x, dtype=float)
    amplitude = float(amplitude)
    if (not math.isfinite(amplitude) or not np.all(np.isfinite(x))
            or np.any(x < geometry["x_corner"]-2e-15)
            or np.any(x > geometry["x_last"]+2e-15)):
        raise ValueError("field jet exits the admitted initial characteristic")
    U0, Va, Vb = (geometry[key] for key in ("U0", "Va", "Vb"))
    V = radius_product(x)/U0
    z = 4*(V-(Va+Vb)/2)/(Vb-Va)
    b, bz = compact_bump(z)
    active = abs(z) < 1
    safe = np.where(active, z, 0.)
    bzz = np.where(active, 8*(1-safe**2)**2*(7*safe**2-1), 0.)
    first = -x*np.exp(x/.8)/.8**2
    second = -np.exp(x/.8)*(1+x/.8)/.8**2
    zx, zxx = 4*first/(U0*(Vb-Va)), 4*second/(U0*(Vb-Va))
    bx, bxx = bz*zx, bzz*zx*zx+bz*zxx
    return amplitude*b/x, amplitude*(bx/x-b/x**2), amplitude*(bxx/x-2*bx/x**2+2*b/x**3)


def validate_geometry(geometry):
    validate_preparation_domain(*(geometry[key] for key in ("U0", "Ua", "Ub", "Va", "Vb")))
    Va, Vb, U0 = (geometry[key] for key in ("Va", "Vb", "U0"))
    support = np.array([Va+(Vb-Va)/4, Va+3*(Vb-Va)/4])
    if (geometry["order"] not in (32, 64) or geometry["x_corner"] != .4
            or not math.isclose(geometry["x_last"], float(radius_from_null(U0, Vb)), rel_tol=0, abs_tol=1e-14)
            or not .4 < geometry["x_last"] < .55
            or not np.allclose(geometry["support_V"], support, rtol=0, atol=1e-14)
            or not np.allclose(geometry["support_x"], radius_from_null(U0, support), rtol=0, atol=1e-14)):
        raise ValueError("changed source domain or support requires a new admission")


def constraint_samples(geometry, *, samples=MAX_SAMPLES, order=64, budget=None,
                       amplitude=1.):
    """Direct stress integrals and integration-by-parts moments on one source leg."""
    if not isinstance(samples, int) or samples < 5 or samples > MAX_SAMPLES or (samples-1) % 4:
        raise ValueError("use at most 129 samples including both exact support edges")
    if order not in (32, 64):
        raise ValueError("only admitted constraint orders 32 and 64 are available")
    validate_geometry(geometry)
    budget = Budget() if budget is None else budget
    budget.check()
    V = np.linspace(geometry["Va"], geometry["Vb"], samples)
    V[(samples-1)//4], V[3*(samples-1)//4] = geometry["support_V"]
    x = radius_from_null(geometry["U0"], V)
    x[0], x[-1] = geometry["x_corner"], geometry["x_last"]
    # Endpoint roundoff should not disguise the exact C3 zeros.
    p, dp, ddp = source_jet(x, geometry, amplitude)
    for i in (0, (samples-1)//4, 3*(samples-1)//4, samples-1):
        p[i] = dp[i] = ddp[i] = 0.
    nodes, weights = np.polynomial.legendre.leggauss(order)
    data = {name: np.zeros(samples) for name in (
        "h", "h_prime", "h_second", "h_direct", "h_prime_direct", "h_einstein",
        "integral_p_squared", "integral_x_dp_squared", "weighted_integral_x_dp_squared")}
    lower, support_last = geometry["support_x"]
    for i, endpoint in enumerate(x):
        upper = min(endpoint, support_last)
        if upper <= lower:
            continue
        s = (lower+upper)/2+(upper-lower)*nodes/2
        w = (upper-lower)*weights/2
        budget.check(len(s))
        q, dq, ddq = source_jet(s, geometry, amplitude)
        constraint = -s*dq*dq/3+s*q*ddq/6
        I0 = float(w @ (q*q))
        I1 = float(w @ (s*dq*dq))
        I2 = float(w @ ((endpoint-s)*s*dq*dq))
        data["integral_p_squared"][i] = I0
        data["integral_x_dp_squared"][i] = I1
        data["weighted_integral_x_dp_squared"][i] = I2
        data["h"][i] = -.5*I2+endpoint*p[i]**2/12-I0/6
        data["h_prime"][i] = -.5*I1+endpoint*p[i]*dp[i]/6-p[i]**2/12
        data["h_second"][i] = -endpoint*dp[i]**2/3+endpoint*p[i]*ddp[i]/6
        data["h_direct"][i] = float(w @ ((endpoint-s)*constraint))
        data["h_prime_direct"][i] = float(w @ constraint)
        data["h_einstein"][i] = -.5*I2
    budget.check()
    return {"V": V, "x": x, "p": p, "p_prime": dp, "p_second": ddp, **data}


def amplitude_parameter(epsilon, kappa_4, length_scale):
    """alpha is dimensionless; a length-unit change must transform all inputs."""
    epsilon, kappa_4, ell = map(float, (epsilon, kappa_4, length_scale))
    if not all(math.isfinite(v) for v in (epsilon, kappa_4, ell)) or kappa_4 <= 0 or ell <= 0:
        raise ValueError("finite epsilon and positive kappa_4 and length scale required")
    return kappa_4*epsilon**2/ell**2


def exact_algebra_checks():
    """Check the derivative/frame identities, independently of quadrature values."""
    import sympy as sp
    if sp.__version__ != "1.14.0":
        raise RuntimeError("exact checks require pinned SymPy 1.14.0")
    x = sp.symbols("x", positive=True)
    p = sp.Function("p")(x)
    S, I, W = (sp.Function(name)(x) for name in ("S", "I", "W"))
    source = -x*sp.diff(p, x)**2/3+x*p*sp.diff(p, x, 2)/6
    hp = -S/2+x*p*sp.diff(p, x)/6-p**2/12
    h = -W/2+x*p**2/12-I/6
    substitutions = {sp.diff(S, x): x*sp.diff(p, x)**2,
                     sp.diff(I, x): p**2, sp.diff(W, x): S}
    residuals = {
        "improved_null_stress_constraint": source+x*(sp.Rational(2, 3)*sp.diff(p, x)**2
                                                        -p*sp.diff(p, x, 2)/3)/2,
        "first_moment_derivative": sp.diff(hp, x).subs(substitutions)-source,
        "second_moment_derivative": sp.diff(h, x).subs(substitutions)-hp,
        "einstein_affine_focusing": source-sp.diff(x*p**2, x, 2)/12
                                    +sp.diff(p**2, x)/6+x*sp.diff(p, x)**2/2,
        "einstein_frame_moment": h-x*p**2/12+I/6+W/2,
    }
    result = {name: str(sp.simplify(value)) for name, value in residuals.items()}
    if any(value != "0" for value in result.values()):
        raise ValueError("exact characteristic identity failed")
    return result


def initialization_calculation(budget=None):
    """One bounded 32/64 calculation; no field re-solve or parameter selection."""
    budget = Budget() if budget is None else budget
    validate_model()
    geometry32, geometry64 = source_geometry(32, budget), source_geometry(64, budget)
    geometric_evaluations = 15*(32+64)
    fine = constraint_samples(geometry64, order=64, budget=budget)
    coarse_quadrature = constraint_samples(geometry64, order=32, budget=budget)
    coarse_geometry = constraint_samples(geometry32, order=64, budget=budget)
    refinements = {}
    for comparison, coarse in (("constraint_quadrature", coarse_quadrature),
                               ("geometry_quadrature", coarse_geometry)):
        refinements[comparison] = {}
        for name in fine:
            difference = abs(fine[name]-coarse[name])
            refinements[comparison][name] = float(np.max(difference))
            if np.any(difference > TOLERANCE_ABSOLUTE+TOLERANCE_RELATIVE*abs(fine[name])):
                raise RuntimeError("32/64 "+comparison+" acceptance failed: "+name)
    direct = {name: float(np.max(abs(fine[name]-fine[name+"_direct"])))
              for name in ("h", "h_prime")}
    for name in direct:
        if np.any(abs(fine[name]-fine[name+"_direct"])
                  > TOLERANCE_ABSOLUTE+TOLERANCE_RELATIVE*abs(fine[name])):
            raise RuntimeError("direct/parts characteristic constraint disagreement")
    frame_error = float(np.max(abs(fine["h"]-fine["x"]*fine["p"]**2/12
                                  +fine["integral_p_squared"]/6-fine["h_einstein"])))
    if frame_error > TOLERANCE_ABSOLUTE:
        raise RuntimeError("Einstein-frame affine comparison failed")
    if not fine["h"][-1] < 0 or not fine["h_prime"][-1] < 0:
        raise RuntimeError("nonzero compact pulse violated analytic post-support sign")
    budget.check()
    return {
        "geometry": geometry64,
        "samples": {key: value.tolist() for key, value in fine.items()},
        "last": {key: float(value[-1]) for key, value in fine.items()},
        "checks": {"quadrature_32_64_max_differences": refinements,
                   "geometry_comparison": "same fractional V positions and independent geometric endpoints; constraint order 64 fixed",
                   "constraint_comparison": "geometric order 64 fixed; constraint orders 32 versus 64",
                   "direct_vs_parts_max_differences": direct,
                   "einstein_affine_moment_max_difference": frame_error,
                   "tolerance_absolute": TOLERANCE_ABSOLUTE,
                   "tolerance_relative": TOLERANCE_RELATIVE,
                   "refinement_is_rigorous_enclosure": False},
        "execution": {"displayed_source_samples": MAX_SAMPLES,
                      "scalar_integrand_evaluations": budget.evaluations-geometric_evaluations,
                      "geometric_integrand_evaluations": geometric_evaluations,
                      "total_integrand_evaluations": budget.evaluations,
                      "maximum_scalar_integrand_evaluations": MAX_EVALUATIONS,
                      "budget_guard_counts_geometry_too": True,
                      "cpu_seconds_max": CPU_SECONDS, "quadrature_orders": [32, 64],
                      "field_solves": 0, "bulk_einstein_evolutions": 0},
    }
