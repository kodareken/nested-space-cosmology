"""Quadratic cap metric and central null-arrival coefficients.

The regular scalar is supplied by the same prepared field. Radial Einstein
constraints retain the transformed cap potential and normalize static time
to the Jordan central clock. The scalar amplitude correction is not solved.
"""
import math
import time

import numpy as np

from .nsc_atlas_bulk_response import BulkBudget
from .nsc_atlas_wall_evolution import _gauss_antiderivatives


SCHEMA = "nsc-atlas-cap-response-v1"
SCOPE = "formal cap metric and central null-arrival coefficients; no corrected scalar amplitude or finite-amplitude remainder"
QUADRATURE_ORDERS = (4, 8, 16, 32, 64)
TIME_PANELS = (4, 8, 16, 32)
MAX_QUERY = 513
DOMAIN_KEYS = ("U0", "Ua", "Ub", "Va", "Vb", "Ta", "Tb", "qa", "qb", "qmin", "xa", "xb")
DEFECT_KEYS = ("wall_einstein_mass_defect", "wall_jordan_mass_defect",
              "wall_frame_mass_defect", "wall_scalar_value_defect",
              "wall_scalar_tangent_defect", "wall_scalar_normal_defect",
              "reference_arrival_rate_defect",
              "wall_normalization_defect", "central_clock_normalization_defect",
              "fixed_V_conversion_defect")


def _number(value, label):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(label+" must be a finite real number")
    try:
        value = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(label+" must be a finite real number") from error
    if not math.isfinite(value):
        raise ValueError(label+" must be a finite real number")
    return value


def _choice(value, choices, label):
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer)) or value not in choices):
        raise ValueError(label+" must be one of "+str(choices))
    return int(value)


def _time(t):
    t = _number(t, "normalized time")
    if not 0 <= t <= 1:
        raise ValueError("cap wall query exits the finite preparation window")
    return t


def _budget(cap_field, budget):
    if budget is None:
        shared = cap_field.budget
        budget = BulkBudget(started=shared.started, cpu_seconds=shared.cpu_seconds,
                            external=shared)
    if not isinstance(budget, BulkBudget) or budget.external is not cap_field.budget:
        raise ValueError("cap response requires the prepared field's shared compute budget")
    budget.check()
    return budget


def _jet(cap_field, T, r, budget):
    T, r = np.broadcast_arrays(np.asarray(T, dtype=float), np.asarray(r, dtype=float))
    shape = T.shape
    if (T.size == 0 or not np.all(np.isfinite(T)) or not np.all(np.isfinite(r))
            or np.any(r < 0) or np.any(r >= 1)):
        raise ValueError("cap metric queries require finite times and 0<=r<1")
    values = {key: [] for key in ("p", "p_T", "p_r", "p_center")}
    for first in range(0, T.size, MAX_QUERY):
        ts, rs = T.ravel()[first:first+MAX_QUERY], r.ravel()[first:first+MAX_QUERY]
        row = cap_field.jet(ts, rs)
        for key in values:
            array = np.asarray(row[key], dtype=float)
            if array.shape != ts.shape or not np.all(np.isfinite(array)):
                raise ValueError("invalid regular cap field jet: "+key)
            values[key].append(array)
        budget.check(metric=len(ts))
    return {key: np.concatenate(parts).reshape(shape) for key, parts in values.items()}


def _radial_cuts(T, upper, edges):
    """Both incoming and reflected profile edges on a constant-T radius."""
    cuts = [0., float(upper)]
    for edge in edges:
        radius = math.tanh(abs(float(edge)-T))
        if 0 < radius < upper:
            cuts.append(radius)
    return np.unique(cuts)


def radial_metric(cap_field, T, r, *, radial_order=8, budget=None):
    """Integrate the regular, centre-anchored radial metric at paired points.

    Every quadrature endpoint is checked by the leading-field domain guard.
    Queries are batched; source/arrival edges split the radial panels. The
    potential contribution is kept separate for independent omission checks.
    """
    radial_order = _choice(radial_order, QUADRATURE_ORDERS, "cap radial order")
    budget = _budget(cap_field, budget)
    T, r = np.broadcast_arrays(np.asarray(T, dtype=float), np.asarray(r, dtype=float))
    if T.size > MAX_QUERY:
        raise ValueError("request at most 513 paired cap metric endpoints")
    endpoint = _jet(cap_field, T, r, budget)
    shape = T.shape
    times, radii = T.ravel(), r.ravel()
    nodes, weights = np.polynomial.legendre.leggauss(radial_order)
    points, point_times, point_weights, owners = [], [], [], []
    for index, (value, upper) in enumerate(zip(times, radii)):
        cuts = _radial_cuts(float(value), float(upper), cap_field.edges)
        for lower, last in zip(cuts[:-1], cuts[1:]):
            half = (last-lower)/2
            if half == 0:
                continue
            points.extend(lower+half*(1+nodes))
            point_times.extend(np.full(radial_order, value))
            point_weights.extend(half*weights)
            owners.extend(np.full(radial_order, index))
    kinetic, potential, lapse = (np.zeros(T.size) for _ in range(3))
    if points:
        points, point_times, point_weights = (np.asarray(array) for array in (points, point_times, point_weights))
        owners = np.asarray(owners, dtype=int)
        jets = _jet(cap_field, point_times, points, budget)
        p, pT, pr = (jets[key] for key in ("p", "p_T", "p_r"))
        A0 = 1-points*points
        density = points*points*(pT*pT/A0+A0*pr*pr)/4
        kinetic = np.bincount(owners, point_weights*density, minlength=T.size)
        potential = np.bincount(owners, point_weights*points*points*p*p/2, minlength=T.size)
        lapse = np.bincount(owners, point_weights*points*(pr*pr+pT*pT/A0**2)/2, minlength=T.size)
    p, pr, pc = (endpoint[key].ravel() for key in ("p", "p_r", "p_center"))
    nuE = kinetic+potential
    delta0 = -pc*pc/12
    deltaE = delta0+lapse
    frame = -radii*radii*(1-radii*radii)*p*pr/6-radii**3*p*p/12
    result = {"nu_E": nuE, "nu_J": nuE+frame, "delta_E": deltaE,
              "delta_E_center": delta0, "frame_mass": frame,
              "mass_kinetic": kinetic, "mass_potential": potential,
              "central_clock_normalization_defect": 2*delta0+pc*pc/6}
    return {"T": T, "r": r, **endpoint,
            **{key: value.reshape(shape) for key, value in result.items()}}


def wall_clock_coefficients(wall, nu_E, delta_E):
    """Cap coordinate-time coefficient and independent induced-normalization residual."""
    x, v, beta = (float(wall[key]) for key in ("x", "v", "beta_minus"))
    p, A, B = (float(wall[key]) for key in ("p", "A", "B_minus"))
    nu_E, delta_E = _number(nu_E, "cap Einstein mass"), _number(delta_E, "cap Einstein lapse")
    A0 = 1-x*x
    if not .45 <= x <= .55:
        raise ValueError("cap wall exits the admitted radius interval")
    e = wall["rho"]-x*p*p/12
    es = wall["rho_s"]-(v*p*p+2*x*p*A)/12
    A1 = -2*x*e-2*nu_E/x
    T1s = (wall["delta_beta_minus"]-x*p*B/6-beta*p*p/12)/A0-beta*delta_E/A0-beta*A1/A0**2
    # This is independent of the radial/timelike equality of the cap mass.
    # It therefore measures that mismatch instead of cancelling identically.
    normalization = (2*beta*T1s+2*delta_E*beta*beta/A0
        +A1*(beta*beta+v*v)/A0**2-2*v*es/A0+p*p/6)
    return {"e": float(e), "e_s": float(es), "A1": float(A1),
            "T1_s": float(T1s), "q0_s": (beta+v)/A0,
            "wall_normalization_defect": float(normalization)}


def null_integrand(r, nu_E, delta_E):
    r, mass, lapse = np.broadcast_arrays(np.asarray(r, dtype=float),
        np.asarray(nu_E, dtype=float), np.asarray(delta_E, dtype=float))
    if (not all(np.all(np.isfinite(value)) for value in (r, mass, lapse))
            or np.any(r < 0) or np.any(r >= 1)):
        raise ValueError("null integrand requires regular finite cap data")
    if np.any((r == 0) & (mass != 0)):
        raise ValueError("regular centre requires zero enclosed Einstein excess mass")
    mass_term = np.divide(2*mass, r, out=np.zeros_like(mass), where=r > 0)
    return mass_term/(1-r*r)**2-lapse/(1-r*r)


def ray_integral(cap_field, q0, x, *, radial_order=8, ray_order=8, budget=None):
    """Metric correction along the unperturbed incoming ray T=q0-atanh(r)."""
    q0, x = _number(q0, "arrival label"), _number(x, "emission radius")
    radial_order = _choice(radial_order, QUADRATURE_ORDERS, "cap radial order")
    ray_order = _choice(ray_order, QUADRATURE_ORDERS, "cap ray order")
    budget = _budget(cap_field, budget)
    if not .45 <= x <= .55:
        raise ValueError("arrival ray must start on the admitted wall strip")
    # Validate the centre and wall before creating any quadrature points.
    _jet(cap_field, np.array([q0, q0-np.arctanh(x)]), np.array([0., x]), budget)
    cuts = [0., x]
    for edge in cap_field.edges:
        if q0 > edge:
            # Reflected q-=q0-2*r* and centre-gauge T=q0-r* transitions.
            for factor in (1., .5):
                point = math.tanh(factor*(q0-float(edge)))
                if 0 < point < x:
                    cuts.append(point)
    cuts = np.unique(cuts)
    nodes, weights = np.polynomial.legendre.leggauss(ray_order)
    points, point_weights = [], []
    for first, last in zip(cuts[:-1], cuts[1:]):
        half = (last-first)/2
        points.extend(first+half*(1+nodes))
        point_weights.extend(half*weights)
    points, point_weights = np.asarray(points), np.asarray(point_weights)
    values = np.empty_like(points)
    for first in range(0, len(points), MAX_QUERY):
        radii = points[first:first+MAX_QUERY]
        metric = radial_metric(cap_field, q0-np.arctanh(radii), radii,
                               radial_order=radial_order, budget=budget)
        values[first:first+MAX_QUERY] = null_integrand(radii, metric["nu_E"], metric["delta_E"])
    positive = float(point_weights @ np.maximum(values, 0.))
    negative = float(point_weights @ np.minimum(values, 0.))
    return {"integral": positive+negative, "positive_integral": positive,
            "negative_integral": negative,
            "integrand_sample_min": float(np.min(values)),
            "integrand_sample_max": float(np.max(values)), "ray_nodes": len(points)}


class CapSolution:
    def __init__(self, cap_field, wall_solution, settings, budget):
        self.cap_field, self.wall_solution = cap_field, wall_solution
        self.settings, self.budget = settings, budget
        self.domain = {key: float(cap_field.domain[key]) for key in DOMAIN_KEYS}
        self.wall_cache, self.ray_cache = {}, {}
        self.maxima = dict.fromkeys(DEFECT_KEYS, 0.)
        self.ray_nodes = 0

    def _wall(self, t):
        t = _time(t)
        if t in self.wall_cache:
            return self.wall_cache[t]
        wall = self.wall_solution._wall(t)
        q0 = float(self.cap_field.geometry["arrival"](wall["V"]))
        T0 = q0-float(np.arctanh(wall["x"]))
        metric = radial_metric(self.cap_field, T0, wall["x"],
                               radial_order=self.settings["radial_order"], budget=self.budget)
        metric = {key: float(value) for key, value in metric.items()}
        clock = wall_clock_coefficients(wall, metric["nu_E"], metric["delta_E"])
        A0 = 1-wall["x"]**2
        tangent = wall["beta_minus"]*metric["p_T"]/A0+wall["v"]*metric["p_r"]
        normal = wall["v"]*metric["p_T"]/A0+wall["beta_minus"]*metric["p_r"]
        defects = {
            "wall_einstein_mass_defect": metric["nu_E"]-wall["integral_L_minus"],
            "wall_jordan_mass_defect": metric["nu_J"]-wall["nu_minus"],
            "wall_frame_mass_defect": metric["frame_mass"]-wall["Q_minus"],
            "wall_scalar_value_defect": metric["p"]-wall["p"],
            "wall_scalar_tangent_defect": tangent-wall["A"],
            "wall_scalar_normal_defect": normal-wall["B_minus"],
            "reference_arrival_rate_defect":
                float(self.cap_field.geometry["arrival"](wall["V"], 1))*wall["Vdot"]-clock["q0_s"],
            "wall_normalization_defect": clock["wall_normalization_defect"],
            "central_clock_normalization_defect": metric["central_clock_normalization_defect"],
        }
        for key, value in defects.items():
            self.maxima[key] = max(self.maxima[key], abs(value))
        row = {"t": t, "V": wall["V"], "x": wall["x"], "T0_wall": T0, "q0": q0,
            "nu_E_wall": metric["nu_E"], "nu_J_wall": metric["nu_J"],
            "delta_E_wall": metric["delta_E"], "nu_E_potential_wall": metric["mass_potential"],
            "integral_L_minus": wall["integral_L_minus"],
            "ds_dt": wall["ds_dt"], "V0_s": wall["Vdot"], **clock, **defects}
        self.wall_cache[t] = row
        return row

    def sample(self, t):
        t = _time(t)
        self.budget.check()
        row = dict(self._wall(t))
        panels = self.settings["time_panels"]
        panel = min(int(t*panels), panels-1)
        local = 2*(t-self.time_cuts[panel])*panels-1
        weights = np.array([polynomial(local) for polynomial in self.time_basis])/(2*panels)
        T1 = self.time_prefix[panel]+float(weights @ self.time_values[panel])
        if t not in self.ray_cache:
            self.ray_cache[t] = ray_integral(self.cap_field, row["q0"], row["x"],
                radial_order=self.settings["radial_order"], ray_order=self.settings["ray_order"], budget=self.budget)
            self.ray_nodes += self.ray_cache[t]["ray_nodes"]
        integral = self.ray_cache[t]["integral"]
        ray_diagnostics = {"ray_"+key: self.ray_cache[t][key] for key in
            ("positive_integral", "negative_integral", "integrand_sample_min", "integrand_sample_max")}
        delta_s = T1+row["e"]/(1-row["x"]**2)+integral
        exterior = self.wall_solution.sample(t)
        shift = row["q0_s"]*exterior["V1"]/row["V0_s"]
        delta_V = delta_s-shift
        conversion = delta_V-delta_s+shift
        self.maxima["fixed_V_conversion_defect"] = max(self.maxima["fixed_V_conversion_defect"], abs(conversion))
        return {**row, **ray_diagnostics, "T1": float(T1), "ray_integral": integral,
            "Delta_C_fixed_s": float(delta_s), "Delta_C_fixed_V": float(delta_V),
            "exterior_label_shift": float(shift), "V1": float(exterior["V1"]),
            "fixed_V_conversion_defect": float(conversion)}

    def summary(self, times=tuple(np.linspace(0, 1, 9))):
        times = tuple(times)
        if not 1 <= len(times) <= 129:
            raise ValueError("display at most 129 cap response time cuts")
        keys = ("t", "V", "x", "T0_wall", "q0", "nu_E_wall", "nu_J_wall",
            "delta_E_wall", "nu_E_potential_wall", "integral_L_minus", "T1", "T1_s",
            "e", "e_s", "A1", "q0_s", "V0_s", "ray_integral", "Delta_C_fixed_s",
            "ray_positive_integral", "ray_negative_integral",
            "ray_integrand_sample_min", "ray_integrand_sample_max",
            "Delta_C_fixed_V", "exterior_label_shift", "V1",
            "wall_einstein_mass_defect", "wall_jordan_mass_defect", "wall_normalization_defect")
        rows = [{key: row[key] for key in keys} for row in (self.sample(t) for t in times)]
        return {"schema": SCHEMA, "scope": SCOPE, "domain": dict(self.domain),
            **self.settings, "clock_time_order": 4, "maximum_absolute_defects": dict(self.maxima),
            "cap_samples": rows, "wall_metric_evaluations": len(self.wall_cache),
            "ray_evaluations": len(self.ray_cache), "ray_nodes": self.ray_nodes,
            "rhs_calls": self.budget.rhs_calls, "metric_values": self.budget.metric_values,
            "source_nodes": self.budget.source_nodes,
            "process_cpu_seconds": time.process_time()-self.budget.started}


def solve_cap_response(cap_field, wall_solution, *, radial_order=8, ray_order=8,
                       time_panels=8, budget=None):
    """Integrate the cap wall clock; central rays are evaluated on requested cuts."""
    radial_order = _choice(radial_order, QUADRATURE_ORDERS, "cap radial order")
    ray_order = _choice(ray_order, QUADRATURE_ORDERS, "cap ray order")
    time_panels = _choice(time_panels, TIME_PANELS, "cap time panels")
    budget = _budget(cap_field, wall_solution.budget if budget is None else budget)
    if budget is not wall_solution.budget:
        raise ValueError("cap and wall must share the cumulative metric budget")
    if wall_solution.budget.external is not cap_field.budget:
        raise ValueError("cap and wall must use the same prepared field budget")
    if (not hasattr(cap_field, "prepared_field") or not hasattr(wall_solution, "wedge")
            or not hasattr(wall_solution.wedge, "prepared_field")
            or cap_field.prepared_field is not wall_solution.wedge.prepared_field):
        raise ValueError("cap and wall must use the same prepared field solution")
    if any(not math.isclose(float(cap_field.domain[key]), float(wall_solution.geometry[key]),
                           rel_tol=0, abs_tol=1e-13) for key in ("U0", "Ua", "Ub", "Va", "Vb")):
        raise ValueError("cap and wall preparation domains differ")
    settings = dict(radial_order=radial_order, ray_order=ray_order, time_panels=time_panels)
    result = CapSolution(cap_field, wall_solution, settings, budget)
    first = result._wall(0.)
    if any(abs(first[key]) > 1e-24 for key in ("nu_E_wall", "delta_E_wall", "T1_s", "e")):
        raise ValueError("cap time evolution requires the quiet zero anchor")
    nodes, result.time_basis = _gauss_antiderivatives()
    result.time_cuts = np.linspace(0., 1., time_panels+1)
    result.time_values = np.zeros((time_panels, 4))
    for panel, lower in enumerate(result.time_cuts[:-1]):
        for index, point in enumerate(lower+(1+nodes)/(2*time_panels)):
            row = result._wall(float(point))
            result.time_values[panel, index] = row["T1_s"]*row["ds_dt"]
    weights = np.array([polynomial(1.) for polynomial in result.time_basis])/(2*time_panels)
    result.time_prefix = np.r_[0., np.cumsum(result.time_values @ weights)]
    budget.check()
    return result
