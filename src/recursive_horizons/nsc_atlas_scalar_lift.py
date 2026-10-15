"""Continuum boundary jets and an exactly integrated finite Hermite lift.

Derivative identities use first scalar jets and mass conservation.  They are
not exact derivatives of the numerical field/clock interpolants; independent
withheld boundary queries must measure that consistency before a PDE uses the
lift.  No exterior scalar correction is evolved here.
"""
import math
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from .nsc_atlas_cap_response import radial_metric
from .nsc_atlas_cap_wall_response import M_OVER_ELL, reference_coefficients
from .nsc_atlas_scalar_boundary import (
    MAX_QUERY, _number, _order, _quadrature, _ray_domain, _shared_budget, _time,
    cap_potential, wall_boundary_coefficients,
)


SCHEMA = "nsc-atlas-scalar-lift-v1"
SCOPE = "finite Hermite boundary and exact polynomial lift; withheld sensitivity is not a uniform error enclosure"
KNOT_COUNTS = (5, 9, 17, 33)
MAX_SAMPLES = 129


def _array(value, name):
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.dtype.kind == "b":
        raise ValueError(name + " must contain finite real values")
    values = np.asarray(raw, dtype=float)
    if values.ndim > 1 or not 1 <= values.size <= MAX_QUERY or not np.all(np.isfinite(values)):
        raise ValueError(name + " needs 1 to 513 finite scalar or one-dimensional values")
    return values


def cap_potential_u(r, nu_E, p, p_T, p_r, f, f_prime):
    """At fixed q_+, compute K1,u from both mass flux and its quiet reduction."""
    cap_potential(r, nu_E, p, p_T, p_r, f, f_prime)  # Shared finite-radius validation.
    r, nu, p, pT, pr, f, fp = np.broadcast_arrays(
        *[np.asarray(value, dtype=float) for value in (r, nu_E, p, p_T, p_r, f, f_prime)])
    A = 1-r*r
    nu_T = r*r*A*pT*pr/2
    nu_r = r*r*(pT*pT/A+A*pr*pr)/4+r*r*p*p/2
    nu_u = (nu_T-A*nu_r)/2
    quiet_nu_u = -A*(1+r*r)*f*f/(8*r*r)
    # K=a(r)*nu+b(r)*f*f'+c(r)*f^2, and dr/du=-A/2.
    a, ar = A/(2*r**3), (r*r-3)/(2*r**4)
    br, cr = (3-r*r)/(12*r**4), (2*r*r-1)/(6*r**5)
    general = a*nu_u-A*(ar*nu+br*f*fp+cr*f*f)/2
    quiet = (A*(3-r*r)*nu/(4*r**4)-A*(3-r*r)*f*fp/(24*r**4)
             +A*(1-8*r*r+3*r**4)*f*f/(48*r**5))
    return {"K1_u": general, "K1_u_quiet": quiet, "nu_u": nu_u,
            "mass_null_identity_defect": nu_u-quiet_nu_u,
            "potential_u_identity_defect": general-quiet}


def outgoing_cap_jet(cap_field, q_minus, q_plus, *, budget,
                     radial_order=16, ray_order=16):
    """Compute J_minus and its partial u derivative in the same q quadrature."""
    budget = _shared_budget(cap_field, budget)
    radial_order, ray_order = _order(radial_order, "radial order"), _order(ray_order, "ray order")
    u, upper, margins = _ray_domain(cap_field, q_minus, q_plus)
    edges = np.asarray(cap_field.edges, dtype=float)
    if edges.ndim != 1 or not 1 <= len(edges) <= 16 or not np.all(np.isfinite(edges)):
        raise ValueError("boundary jets need finite bounded profile edges")
    labels, radii, weights = _quadrature(u, float(cap_field.domain["qa"]), upper,
                                        edges, ray_order, "q")
    positive = negative = integral_u = 0.
    maxima = dict.fromkeys(("potential_identity_defect", "potential_u_identity_defect",
                           "mass_null_identity_defect"), 0.)
    for first in range(0, len(labels), MAX_QUERY):
        qs, rs, ws = [array[first:first+MAX_QUERY] for array in (labels, radii, weights)]
        metric = radial_metric(cap_field, (u+qs)/2, rs, radial_order=radial_order, budget=budget)
        profile = cap_field.profile(qs)
        args = (rs, metric["nu_E"], metric["p"], metric["p_T"], metric["p_r"],
                profile["f"], profile["f_prime"])
        potential, derivative = cap_potential(*args), cap_potential_u(*args)
        source = -potential["K1"]*profile["f"]
        positive += float(ws @ np.maximum(source, 0.))
        negative += float(ws @ np.minimum(source, 0.))
        integral_u += float(ws @ (derivative["K1_u"]*profile["f"]))
        for key in maxima:
            values = potential[key] if key in potential else derivative[key]
            maxima[key] = max(maxima[key], float(np.max(abs(values))))
        budget.check()
    radius = math.tanh((upper-u)/2)
    endpoint = radial_metric(cap_field, (u+upper)/2, radius,
                             radial_order=radial_order, budget=budget)
    profile = cap_field.profile(np.array([upper]))
    f, fp = float(profile["f"][0]), float(profile["f_prime"][0])
    endpoint_potential = cap_potential(radius, endpoint["nu_E"], endpoint["p"],
                                       endpoint["p_T"], endpoint["p_r"], f, fp)
    maxima["potential_identity_defect"] = max(maxima["potential_identity_defect"],
                                               abs(float(endpoint_potential["potential_identity_defect"])))
    result = {"J_minus": positive+negative, "J_minus_u": -integral_u,
              "integral_Ku_f": integral_u, "K1_wall": float(endpoint_potential["K1"]),
              "wall_profile": f, "positive_integral": positive, "negative_integral": negative,
              "ray_nodes": len(labels), **maxima, **margins}
    if not all(math.isfinite(float(value)) for value in result.values()):
        raise ValueError("boundary jet produced nonfinite values")
    return result


def wall_boundary_derivative(wall, J_minus, J_minus_s):
    """Material derivatives of all wall pullbacks, including the changing Vdot."""
    base = wall_boundary_coefficients(wall, J_minus)
    keys = ("x", "v", "beta_plus", "beta_minus", "p", "w", "j", "A", "rho_s",
            "V", "V1", "V1_s", "Vdot")
    x, v, bp, bm, p, w, j, p_s, rho_s, V, V1, V1_s, Vs = [_number(wall[key], key) for key in keys]
    if V >= 0:
        raise ValueError("boundary derivative requires the retained negative V branch")
    J_minus_s = _number(J_minus_s, "outgoing material derivative")
    reference = reference_coefficients(x)
    v_s = float(reference["radial_acceleration"])
    bp_s, bm_s = float(reference["beta_plus_x"])*v, float(reference["beta_minus_x"])*v
    d, A0, qm_s = base["d"], 1-x*x, base["q_minus_s"]
    d_s, A0_s = v_s-bp_s, -2*x*v
    qm_ss = (bm_s-v_s)/A0-(bm-v)*A0_s/A0**2
    qp_s = (bm+v)/A0
    multiplier, multiplier_s = 2*qm_s/d, 2*qm_ss/d-2*qm_s*d_s/d**2
    Vs_s = Vs*(Vs/V-d_s/d)
    e_s = rho_s-(v*p*p+2*x*p*p_s)/12
    E_s = e_s-v_s*V1/Vs-v*V1_s/Vs+v*V1*Vs_s/Vs**2
    C_s = 2*M_OVER_ELL*(j-w*(3*v/x+2*d_s/d))/(x**3*d*d)
    derivative = (multiplier_s*J_minus+multiplier*J_minus_s
                  -E_s*base["c0_rr"]-base["displacement_at_fixed_V"]*C_s)
    rU = -4*M_OVER_ELL**2*V*math.exp(-x/(2*M_OVER_ELL))/x
    rU_s = rU*(Vs/V-v*(1/(2*M_OVER_ELL)+1/x))
    b = rU*base["c1_r"]
    b_s = rU_s*base["c1_r"]+rU*derivative
    return {**base, "c1_r_s": derivative, "b": b, "b_s": b_s, "b_V": b_s/Vs,
            "r_U": rU, "r_U_s": rU_s, "d_s": d_s, "q_minus_ss": qm_ss,
            "q_plus_s": qp_s, "V0_ss": Vs_s, "e_s": e_s,
            "displacement_at_fixed_V_s": E_s, "c0_rr_s": C_s}


class BoundaryJet:
    """Shared cache of continuum boundary values/slopes for nested lift fits."""

    def __init__(self, boundary):
        self.boundary = boundary
        self.cap_field, self.wall_solution = boundary.cap_field, boundary.wall_solution
        self.budget = _shared_budget(self.cap_field, boundary.budget)
        if self.wall_solution.budget is not self.budget:
            raise ValueError("boundary jets require the same cumulative wall budget")
        if self.cap_field.prepared_field is not self.wall_solution.wedge.prepared_field:
            raise ValueError("boundary jets require the same prepared field")
        self.domain, self.settings = dict(boundary.domain), dict(boundary.settings)
        self.geometry = self.cap_field.geometry
        self.cache, self.ray_nodes = {}, 0

    def sample(self, t):
        t = _time(t)
        self.budget.check()
        if t in self.cache:
            return dict(self.cache[t])
        if len(self.cache) >= MAX_SAMPLES:
            raise ValueError("retain at most 129 boundary jet samples")
        wall = self.wall_solution.sample(t)
        V = float(wall["V"])
        q = float(self.geometry["arrival"](V))
        u = q-2*math.atanh(wall["x"])
        jet = outgoing_cap_jet(self.cap_field, u, q, budget=self.budget, **self.settings)
        x, v, bm = [float(wall[key]) for key in ("x", "v", "beta_minus")]
        qm_s, qp_s = (bm-v)/(1-x*x), (bm+v)/(1-x*x)
        Js = -qp_s*jet["K1_wall"]*jet["wall_profile"]-qm_s*jet["integral_Ku_f"]
        derivative = wall_boundary_derivative(wall, jet["J_minus"], Js)
        U, U_V = float(self.geometry["wall"](V)), float(self.geometry["wall"](V, 1))
        if not math.isfinite(U_V) or U_V <= 0:
            raise ValueError("boundary Hermite map needs a positive wall derivative")
        stored = self.boundary.cache.get(t)
        row = {"t": t, "U": U, "V": V, "V0_s": float(wall["Vdot"]),
               "x": x, "q_minus": u, "q_plus": q,
               **jet, **derivative, "J_minus_s": Js, "U_V": U_V,
               "b_U": derivative["b_V"]/U_V,
               "arrival_rate_defect": float(self.geometry["arrival"](V, 1))*wall["Vdot"]-qp_s,
               "wall_coordinate_rate_defect": U_V*wall["Vdot"]-derivative["d"]/(2*derivative["r_U"]),
               "wall_profile_defect": jet["wall_profile"]-wall["w"],
               "boundary_cache_checked": stored is not None,
               "boundary_cache_value_defect": 0. if stored is None else derivative["c1_r"]-stored["c1_r"]}
        self.ray_nodes += jet["ray_nodes"]
        self.cache[t] = row
        return dict(row)


def boundary_derivative(boundary, t):
    """One uncached query; use BoundaryJet for multiple or nested fit queries."""
    return BoundaryJet(boundary).sample(t)


class BoundaryLift:
    """Hermite boundary in U_W and its exactly integrated polynomial lift.

    Pass the same BoundaryJet to coarse/fine fits to share physical queries.
    Validation points are never inserted into a fit.  Reported withheld
    residuals are observed numerical differences, not certified sup bounds.
    """

    def __init__(self, boundary, *, knots=9, validation_times=()):
        if isinstance(knots, (bool, np.bool_)) or not isinstance(knots, (int, np.integer)) or knots not in KNOT_COUNTS:
            raise ValueError("lift knot count must be one of " + str(KNOT_COUNTS))
        self.jets = boundary if isinstance(boundary, BoundaryJet) else BoundaryJet(boundary)
        self.budget, self.domain = self.jets.budget, dict(self.jets.domain)
        self.geometry = self.jets.geometry
        self.knot_count = int(knots)
        self.knot_times = np.linspace(0., 1., knots)
        self.rows = [self.jets.sample(float(t)) for t in self.knot_times]
        self.U_knots = np.array([row["U"] for row in self.rows])
        self.b_knots = np.array([row["b"] for row in self.rows])
        self.b_U_knots = np.array([row["b_U"] for row in self.rows])
        if (not all(np.all(np.isfinite(a)) for a in (self.U_knots, self.b_knots, self.b_U_knots))
                or np.any(np.diff(self.U_knots) <= 0)):
            raise ValueError("lift knots must be finite and strictly increasing in U")
        tolerance = 8*np.finfo(float).eps*max(1., abs(self.domain["Ua"]), abs(self.domain["Ub"]))
        if (abs(self.U_knots[0]-self.domain["Ua"]) > tolerance
                or abs(self.U_knots[-1]-self.domain["Ub"]) > tolerance):
            raise ValueError("lift knots must cover the full finite wall domain")
        self.U_knots[[0, -1]] = [self.domain["Ua"], self.domain["Ub"]]
        if abs(self.b_knots[0]) > 1e-24 or abs(self.b_U_knots[0]) > 1e-24:
            raise ValueError("lift requires the quiet value and derivative anchor")
        self.b_knots[0] = self.b_U_knots[0] = 0.
        self.polynomial = CubicHermiteSpline(self.U_knots, self.b_knots, self.b_U_knots, extrapolate=False)
        self.antiderivative = self.polynomial.antiderivative()
        self.anchor = float(self.antiderivative(self.domain["Ua"]))
        self.validation_rows = []
        self.validate(validation_times)
        self.budget.check()

    def _U(self, U):
        values = _array(U, "lift U query")
        if np.any(values < self.domain["U0"]) or np.any(values > self.domain["Ub"]):
            raise ValueError("lift query exits the retained U domain")
        self.budget.check()
        return values

    def _V(self, V):
        values = _array(V, "boundary V query")
        if np.any(values < self.domain["Va"]) or np.any(values > self.domain["Vb"]):
            raise ValueError("boundary query exits the finite V domain")
        return values

    def _wall_U(self, V):
        values = self._V(V)
        U = np.asarray(self.geometry["wall"](values), dtype=float)
        tolerance = 8*np.finfo(float).eps*max(1., abs(self.domain["Ua"]), abs(self.domain["Ub"]))
        if (U.shape != values.shape or not np.all(np.isfinite(U))
                or np.any(U < self.domain["Ua"]-tolerance) or np.any(U > self.domain["Ub"]+tolerance)):
            raise ValueError("wall map exits the finite lift domain")
        return np.clip(U, self.domain["Ua"], self.domain["Ub"])

    def b_at_U(self, U):
        """Boundary function parametrized by the wall's U coordinate."""
        values = self._U(U)
        active = values > self.domain["Ua"]
        return np.where(active, self.polynomial(np.maximum(values, self.domain["Ua"])), 0.)

    def b_U(self, U):
        """Derivative of the interpolated boundary with respect to wall U."""
        values = self._U(U)
        active = values > self.domain["Ua"]
        return np.where(active, self.polynomial(np.maximum(values, self.domain["Ua"]), 1), 0.)

    def B(self, U):
        values = self._U(U)
        active = values > self.domain["Ua"]
        integral = self.antiderivative(np.maximum(values, self.domain["Ua"]))-self.anchor
        return np.where(active, integral, 0.)

    def b(self, V):
        return self.b_at_U(self._wall_U(V))

    def B_at_wall(self, V):
        return self.B(self._wall_U(V))

    def validate(self, times):
        times = tuple(times)
        if len(times) > MAX_SAMPLES or len(self.validation_rows)+len(times) > MAX_SAMPLES:
            raise ValueError("retain at most 129 withheld queries")
        for raw in times:
            t = _time(raw)
            if np.any(abs(self.knot_times-t) <= 4*np.finfo(float).eps):
                raise ValueError("withheld query must not be a fitted knot")
            if any(row["t"] == t for row in self.validation_rows):
                continue
            row = self.jets.sample(t)
            U, V = row["U"], row["V"]
            index = min(int(np.searchsorted(self.U_knots, U, side="right")-1), len(self.U_knots)-2)
            width = self.U_knots[index+1]-self.U_knots[index]
            fitted = float(self.polynomial(U))
            slope = float(self.polynomial(U, 1))
            self.validation_rows.append({"t": t, "U": U, "V": V, "b_reference": row["b"],
                "b_interpolated": fitted, "b_defect": fitted-row["b"],
                "b_U_reference": row["b_U"], "b_U_interpolated": slope,
                "scaled_slope_defect": width*(slope-row["b_U"]),
                "B_lift": float(self.B(U)), "B_wall": float(self.B_at_wall(V)),
                "lift_derivative_defect": float(self.antiderivative(U, 1))-fitted,
                "wall_transport_defect": float(self.antiderivative(U, 1))*row["U_V"]-fitted*row["U_V"]})
        self.budget.check()
        return [dict(row) for row in self.validation_rows]

    def summary(self):
        defects = ("b_defect", "scaled_slope_defect", "lift_derivative_defect", "wall_transport_defect")
        return {"schema": SCHEMA, "scope": SCOPE, "domain": dict(self.domain),
                **self.jets.settings, "interpolation_variable": "U_W",
                "knot_count": self.knot_count, "knot_samples": [dict(row) for row in self.rows],
                "validation_samples": [dict(row) for row in self.validation_rows],
                "withheld_evaluations": len(self.validation_rows),
                "maximum_absolute_defects": {key: max((abs(row[key]) for row in self.validation_rows), default=0.)
                                               for key in defects},
                "terminal_lift": float(self.B(self.domain["Ub"])),
                "quiet_anchor": {"b": float(self.b_at_U(self.domain["Ua"])),
                                 "b_U": float(self.b_U(self.domain["Ua"])),
                                 "B": float(self.B(self.domain["Ua"]))},
                "boundary_jet_evaluations": len(self.jets.cache), "ray_nodes": self.jets.ray_nodes,
                "metric_values": self.budget.metric_values, "scalar_kernel_values": self.jets.cap_field.budget.kernel_values,
                "process_cpu_seconds": time.process_time()-self.budget.started}
