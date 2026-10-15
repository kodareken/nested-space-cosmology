"""Known cubic outgoing cap forcing and displaced exterior wall derivative.

This adapter evaluates a boundary input on the established no-return window.
It does not evolve the scalar correction or construct a boundary interpolant.
"""
import math
import time

import numpy as np

from .nsc_atlas_bulk_response import BulkBudget
from .nsc_atlas_cap_response import radial_metric
from .nsc_atlas_cap_wall_response import M_OVER_ELL


SCHEMA = "nsc-atlas-scalar-boundary-v1"
SCOPE = "formal cubic outgoing boundary input; no propagated scalar correction or finite-amplitude remainder"
QUADRATURE_ORDERS = (4, 8, 16, 32)
MAX_QUERY = 513
DOMAIN_KEYS = ("U0", "Ua", "Ub", "Va", "Vb", "qa", "qb", "Ta", "Tb", "qmin", "xa", "xb")
DEFECT_KEYS = ("coordinate_integral_defect", "potential_identity_defect",
               "wall_flux_defect", "wall_profile_defect")
MARGIN_KEYS = ("no_return_margin", "quiet_centre_time_margin", "initial_time_margin",
               "last_profile_margin", "first_profile_margin", "active_radius_lower")


def _number(value, name):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(name + " must be a finite real scalar")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(name + " must be a finite real scalar")
    return value


def _order(value, name):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
            or value not in QUADRATURE_ORDERS):
        raise ValueError(name + " must be one of " + str(QUADRATURE_ORDERS))
    return int(value)


def _time(value):
    value = _number(value, "normalized wall time")
    if not 0 <= value <= 1:
        raise ValueError("wall time exits the preparation domain")
    return value


def _shared_budget(cap_field, budget):
    if not isinstance(budget, BulkBudget) or budget.external is not cap_field.budget:
        raise ValueError("boundary queries require the shared prepared-field budget")
    if budget.started != cap_field.budget.started:
        raise ValueError("boundary queries require the cumulative CPU clock")
    budget.check()
    return budget


def cap_potential(r, nu_E, p, p_T, p_r, f, f_prime):
    """General first-jet potential and its independent quiet-profile form.

    The second expression assumes T<q_a and p=f(T+atanh(r))/r.  Comparing
    them checks the scalar coordinate conversion, not the radial quadrature.
    """
    arrays = np.broadcast_arrays(*[np.asarray(value, dtype=float)
                                   for value in (r, nu_E, p, p_T, p_r, f, f_prime)])
    r, nu_E, p, p_T, p_r, f, f_prime = arrays
    if (not all(np.all(np.isfinite(value)) for value in arrays)
            or np.any(r <= 0) or np.any(r >= 1)):
        raise ValueError("outgoing potential needs finite data at 0<r<1")
    A0 = 1-r*r
    general = A0*(2*nu_E/r**3-p*p/3-p_T*p_T/(6*A0)+A0*p_r*p_r/6)/4
    quiet = A0*nu_E/(2*r**3)-A0*f*f_prime/(12*r**3)+A0*(1-3*r*r)*f*f/(24*r**4)
    return {"K1": general, "K1_quiet": quiet, "potential_identity_defect": general-quiet}


def wall_boundary_coefficients(wall, J_minus):
    """Known backscatter plus material-to-areal wall pullback, without Hessians."""
    keys = ("x", "v", "beta_plus", "beta_minus", "p", "w", "rho", "V1", "Vdot")
    x, v, bp, bm, p, w, rho, V1, Vdot = [_number(wall[key], key) for key in keys]
    J_minus = _number(J_minus, "outgoing correction derivative")
    d, A0 = v-bp, 1-x*x
    if not .45 <= x <= .55 or d <= 0 or Vdot <= 0 or bm-v <= 0:
        raise ValueError("wall coefficients leave the admitted signed reference branch")
    qm_s = (bm-v)/A0
    e = rho-x*p*p/12
    displacement_at_V = e-v*V1/Vdot
    c0_rr = 2*M_OVER_ELL*w/(x**3*d*d)
    backscatter = 2*qm_s*J_minus/d
    displacement = -displacement_at_V*c0_rr
    c1_r = backscatter+displacement
    flux_defect = d*(c1_r+displacement_at_V*c0_rr)/2-qm_s*J_minus
    return {"c1_r": c1_r, "backscatter_term": backscatter,
            "displacement_term": displacement, "e": e,
            "displacement_at_fixed_V": displacement_at_V, "c0_rr": c0_rr,
            "q_minus_s": qm_s, "d": d, "wall_flux_defect": flux_defect}


def _ray_domain(cap_field, u, q):
    u, q = _number(u, "outgoing label"), _number(q, "incoming wall label")
    domain = cap_field.domain
    qa, qb, Ta, Tb, qmin = [float(domain[key]) for key in ("qa", "qb", "Ta", "Tb", "qmin")]
    tolerance = 8*np.finfo(float).eps*max(1., abs(Ta), abs(qmin), abs(qb))
    if not (Ta < Tb < qa < qb) or not u < qa:
        raise ValueError("outgoing boundary requires the certified quiet-centre no-return domain")
    if u < qmin-tolerance or q < qa-tolerance or q > qb+tolerance:
        raise ValueError("outgoing ray exits the finite profile domain")
    # Clamp only floating endpoint roundoff; retain actual labels in the result.
    upper = min(qb, max(qa, q))
    times = np.array([(u+qa)/2, (u+upper)/2])
    radii = np.tanh((np.array([qa, upper])-u)/2)
    if np.any(times < Ta-tolerance) or np.any(times > Tb+tolerance):
        raise ValueError("outgoing ray exits the known pre-arrival cap times")
    limit = np.asarray(cap_field.radial_limit(times), dtype=float)
    if limit.shape != times.shape or np.any(radii > limit+tolerance):
        raise ValueError("outgoing ray exits the retained cap radius")
    margins = dict(no_return_margin=qa-u, quiet_centre_time_margin=qa-times[-1],
                   initial_time_margin=times[0]-Ta, last_profile_margin=qb-q,
                   first_profile_margin=u-qmin, active_radius_lower=float(radii[0]))
    return u, upper, margins


def _quadrature(u, lower, upper, edges, order, coordinate):
    cuts = np.unique([lower, *[float(edge) for edge in edges if lower < edge < upper], upper])
    if coordinate == "r":
        cuts = np.tanh((cuts-u)/2)
    elif coordinate != "q":
        raise ValueError("outgoing integral coordinate must be q or r")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    points, weighted = [], []
    for first, last in zip(cuts[:-1], cuts[1:]):
        half = (last-first)/2
        points.extend(first+half*(1+nodes))
        weighted.extend(half*weights)
    points, weighted = np.asarray(points), np.asarray(weighted)
    if coordinate == "q":
        labels, radii = points, np.tanh((points-u)/2)
    else:
        radii, labels = points, u+2*np.arctanh(points)
        weighted = weighted*2/(1-radii*radii)
    return labels, radii, weighted


def outgoing_cap_flux(cap_field, q_minus, q_plus, *, budget,
                      radial_order=16, ray_order=16, coordinate="q"):
    """Evaluate J_minus=-integral K1*c0 dq on one no-return outgoing ray.

    The r-coordinate rule uses independently placed Gauss nodes and dq/dr.
    Its discrepancy from the q rule is a numerical diagnostic, not an error
    enclosure.  No centre derivative or new physical evolution is required.
    """
    budget = _shared_budget(cap_field, budget)
    radial_order, ray_order = _order(radial_order, "radial order"), _order(ray_order, "ray order")
    if coordinate not in ("q", "r"):
        raise ValueError("outgoing integral coordinate must be q or r")
    u, upper, margins = _ray_domain(cap_field, q_minus, q_plus)
    edges = np.asarray(cap_field.edges, dtype=float)
    if edges.ndim != 1 or not 1 <= len(edges) <= 16 or not np.all(np.isfinite(edges)):
        raise ValueError("outgoing integration needs finite bounded profile edges")
    labels, radii, weights = _quadrature(u, float(cap_field.domain["qa"]), upper, edges, ray_order, coordinate)
    positive = negative = defect = 0.
    min_K = max_K = min_source = max_source = 0.
    if len(labels):
        min_K = min_source = math.inf
        max_K = max_source = -math.inf
    for first in range(0, len(labels), MAX_QUERY):
        qs, rs, ws = (array[first:first+MAX_QUERY] for array in (labels, radii, weights))
        metric = radial_metric(cap_field, (u+qs)/2, rs, radial_order=radial_order, budget=budget)
        profile = cap_field.profile(qs)
        values = cap_potential(rs, metric["nu_E"], metric["p"], metric["p_T"], metric["p_r"],
                               profile["f"], profile["f_prime"])
        source = -values["K1"]*profile["f"]
        if source.shape != qs.shape or not np.all(np.isfinite(source)):
            raise ValueError("invalid outgoing cap forcing")
        positive += float(ws @ np.maximum(source, 0.))
        negative += float(ws @ np.minimum(source, 0.))
        defect = max(defect, float(np.max(abs(values["potential_identity_defect"]))))
        min_K, max_K = min(min_K, float(np.min(values["K1"]))), max(max_K, float(np.max(values["K1"])))
        min_source, max_source = min(min_source, float(np.min(source))), max(max_source, float(np.max(source)))
        budget.check()
    return {"J_minus": positive+negative, "positive_integral": positive,
            "negative_integral": negative, "absolute_integral": positive-negative,
            "potential_identity_defect": defect, "K1_sample_min": min_K,
            "K1_sample_max": max_K, "source_sample_min": min_source,
            "source_sample_max": max_source, "ray_nodes": len(labels), **margins}


class ScalarBoundary:
    """Cached boundary coefficients from one cap field and cumulative wall solve."""

    def __init__(self, cap_field, wall_solution, *, radial_order=16, ray_order=16, budget=None):
        self.settings = {"radial_order": _order(radial_order, "radial order"),
                         "ray_order": _order(ray_order, "ray order")}
        budget = wall_solution.budget if budget is None else budget
        self.budget = _shared_budget(cap_field, budget)
        if budget is not wall_solution.budget:
            raise ValueError("boundary and wall require the same cumulative metric budget")
        if (not hasattr(wall_solution, "wedge")
                or cap_field.prepared_field is not wall_solution.wedge.prepared_field):
            raise ValueError("boundary and wall require the same prepared field solution")
        if any(not math.isclose(float(cap_field.domain[key]), float(wall_solution.geometry[key]),
                               rel_tol=0, abs_tol=1e-13) for key in DOMAIN_KEYS[:5]):
            raise ValueError("cap and wall preparation domains differ")
        self.cap_field, self.wall_solution = cap_field, wall_solution
        self.domain = {key: float(cap_field.domain[key]) for key in DOMAIN_KEYS}
        certificate = cap_field.quiet_certificate
        if (certificate.get("quiet_cap") is not True
                or not self.domain["Ta"] < self.domain["Tb"] < self.domain["qa"] < self.domain["qb"]):
            raise ValueError("boundary adapter requires the certified pre-arrival quiet domain")
        self.cache = {}
        self.maxima = dict.fromkeys(DEFECT_KEYS, 0.)
        self.minimum_margins = dict.fromkeys(MARGIN_KEYS, math.inf)
        self.ray_nodes_q = self.ray_nodes_r = 0

    def sample(self, t):
        t = _time(t)
        self.budget.check()
        if t in self.cache:
            return dict(self.cache[t])
        if len(self.cache) >= 129:
            raise ValueError("retain at most 129 boundary samples")
        wall = self.wall_solution.sample(t)
        q = float(self.cap_field.geometry["arrival"](wall["V"]))
        u = q-2*math.atanh(wall["x"])
        q_integral = outgoing_cap_flux(self.cap_field, u, q, budget=self.budget, **self.settings)
        r_integral = outgoing_cap_flux(self.cap_field, u, q, budget=self.budget, coordinate="r", **self.settings)
        boundary = wall_boundary_coefficients(wall, q_integral["J_minus"])
        profile = self.cap_field.profile(np.array([q]))
        defects = {"coordinate_integral_defect": q_integral["J_minus"]-r_integral["J_minus"],
                   "potential_identity_defect": max(q_integral["potential_identity_defect"],
                                                       r_integral["potential_identity_defect"]),
                   "wall_flux_defect": boundary["wall_flux_defect"],
                   "wall_profile_defect": float(profile["f"][0])-wall["w"]}
        for key, value in defects.items():
            self.maxima[key] = max(self.maxima[key], abs(value))
        for key in MARGIN_KEYS:
            self.minimum_margins[key] = min(self.minimum_margins[key], q_integral[key])
        self.ray_nodes_q += q_integral["ray_nodes"]
        self.ray_nodes_r += r_integral["ray_nodes"]
        result = {"t": t, "V": float(wall["V"]), "x": float(wall["x"]),
                  "v": float(wall["v"]), "V0_s": float(wall["Vdot"]),
                  "beta_plus": float(wall["beta_plus"]), "beta_minus": float(wall["beta_minus"]),
                  "q_minus": u, "q_plus": q,
                  "T_wall": (u+q)/2, "rho": float(wall["rho"]), "V1": float(wall["V1"]),
                  "w": float(wall["w"]), "p_wall": float(wall["p"]),
                  **q_integral, **boundary, **defects, "J_minus_r": r_integral["J_minus"],
                  "ray_nodes_r": r_integral["ray_nodes"]}
        self.cache[t] = result
        return dict(result)

    def summary(self, times=tuple(np.linspace(0., 1., 9))):
        times = tuple(times)
        if not 1 <= len(times) <= 129:
            raise ValueError("display 1 to 129 boundary samples")
        samples = [self.sample(value) for value in times]
        return {"schema": SCHEMA, "scope": SCOPE, "domain": dict(self.domain), **self.settings,
                "boundary_samples": samples, "maximum_absolute_defects": dict(self.maxima),
                "minimum_domain_margins": dict(self.minimum_margins),
                "boundary_evaluations": len(self.cache), "ray_nodes_q": self.ray_nodes_q,
                "ray_nodes_r": self.ray_nodes_r, "metric_values": self.budget.metric_values,
                "rhs_calls": self.budget.rhs_calls, "source_nodes": self.budget.source_nodes,
                "scalar_kernel_values": self.cap_field.budget.kernel_values,
                "profile_counts": dict(self.cap_field.evaluation_counts),
                "process_cpu_seconds": time.process_time()-self.budget.started}
