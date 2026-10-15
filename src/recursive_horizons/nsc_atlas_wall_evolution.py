"""Quadratic wall masses, displacement and exterior null-clock coefficient.

The scalar is the unchanged prepared field. Radial Einstein-frame constraints
are integrated only in the retained wedge. The independent timelike mass laws
drive the signed material junction. All outputs are coefficients of alpha,
not a finite-amplitude solution or the cap's central-observer clock.
"""
import math
import time

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import solve_ivp

from .nsc_atlas_bulk_response import A_RADIUS, M, BulkBudget, background, source_traces
from .nsc_atlas_cap_wall_response import (
    K_ELL, reference_coefficients, field_mass_ingredients, wall_response_coefficients,
)
from .nsc_atlas_characteristic_initialization import validate_geometry


SCHEMA = "nsc-atlas-wall-evolution-v1"
TIME_PANELS = (4, 8, 16, 32)
RADIAL_ORDERS = (2, 4)
MAX_QUERY = 513


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
        raise ValueError("wall query exits the finite preparation window")
    return t


def _budget(wedge, budget):
    if budget is None:
        shared = wedge.budget
        budget = BulkBudget(started=shared.started,
                            cpu_seconds=shared.cpu_seconds, external=shared)
    if not isinstance(budget, BulkBudget) or budget.external is not wedge.budget:
        raise ValueError("wall calculations require the field's shared compute budget")
    budget.check()
    return budget


def _domain(wedge, geometry):
    validate_geometry(geometry)
    if any(not math.isclose(float(wedge.domain[key]), geometry[key], rel_tol=0, abs_tol=1e-13)
           for key in ("U0", "Ua", "Ub", "Va", "Vb")):
        raise ValueError("field and wall-response domains do not agree")
    grid = np.asarray(wedge.grid, dtype=float)
    if (grid.ndim != 1 or not 2 <= len(grid) <= MAX_QUERY
            or not np.all(np.isfinite(grid)) or np.any(np.diff(grid) <= 0)
            or grid[0] != geometry["U0"] or grid[-1] != geometry["Ub"]):
        raise ValueError("radial integration requires the original finite field grid")


def _field(wedge, points, t, order, budget):
    rows = []
    keys = ("p", "p_U", "p_V", "chi", "chi_U", "chi_V",
            "chi_U_lower_limit_bound")
    for first in range(0, len(points), MAX_QUERY):
        current = points[first:first+MAX_QUERY]
        row = wedge.bulk(current, t, order=order)
        for key in keys:
            value = np.asarray(row[key], dtype=float)
            if value.shape != current.shape or not np.all(np.isfinite(value)):
                raise ValueError("invalid wedge field jet: "+key)
            row[key] = value
        rows.append(row)
        budget.check(metric=len(current))
    return {key: np.concatenate([row[key] for row in rows]) for key in keys}


def radial_response(wedge_field, geometry, t, *, upper_U=None, radial_order=2,
                    field_order=32, source_order=64, source_scale=1., budget=None):
    """Integrate radial mass/lapse, optionally ending on the old rectangle.

    The original scalar interpolation cells are quadrature panels. Jordan
    reconstruction uses the generic endpoint jets, so an interior rectangle
    endpoint is never silently given the wall condition chi_U=0.
    """
    _domain(wedge_field, geometry)
    t = _time(t)
    radial_order = _choice(radial_order, RADIAL_ORDERS, "radial order")
    field_order = _choice(field_order, (32, 64), "field derivative order")
    source_order = _choice(source_order, (32, 64), "source order")
    budget = _budget(wedge_field, budget)
    V = min(geometry["Vb"], max(geometry["Va"],
            geometry["Va"]+(geometry["Vb"]-geometry["Va"])*t))
    wall = float(wedge_field.geometry["wall"](V))
    upper = wall if upper_U is None else _number(upper_U, "radial endpoint")
    lower = geometry["U0"]
    if not lower <= upper <= wall:
        raise ValueError("radial integration exits the retained exterior wedge")
    grid = np.asarray(wedge_field.grid)
    cuts = np.r_[lower, grid[(grid > lower) & (grid < upper)], upper]
    node, weight = np.polynomial.legendre.leggauss(radial_order)
    half = np.diff(cuts)/2
    points = ((cuts[:-1]+half)[:, None]+half[:, None]*node).ravel()
    weights = (half[:, None]*weight).ravel()
    all_points = np.r_[lower, points, upper]
    sf = _field(wedge_field, all_points, t, field_order, budget)
    source = source_traces(V, geometry, order=source_order,
                           source_scale=source_scale, budget=budget)
    if not np.allclose([sf["p"][0], sf["p_V"][0]],
                       [source["p"], source["p_V"]], rtol=2e-9, atol=2e-10):
        raise ValueError("prepared field is incompatible with the fixed source gauge/data")
    g = background(all_points, V)
    r, rU, rV, C = (g[key] for key in ("r", "r_U", "r_V", "C"))
    p, pU, pV = (sf[key] for key in ("p", "p_U", "p_V"))
    muE = source["mu_E"]-float(weights @ (r[1:-1]**2*rV[1:-1]*pU[1:-1]**2/C[1:-1]))
    zetaE = source["zeta_E"]+float(weights @ (r[1:-1]*pU[1:-1]**2/(2*rU[1:-1])))
    Q = M*p[-1]**2/12+r[-1]**2*p[-1]*(rU[-1]*pV[-1]+rV[-1]*pU[-1])/(3*C[-1])
    zetaJ = zetaE+p[-1]**2/12-r[-1]*p[-1]*pU[-1]/(6*rU[-1])
    return {"t": t, "U": upper, "V": V, "x": float(r[-1]),
            "mu_E": muE, "mu_J": muE+float(Q), "frame_mass": float(Q),
            "zeta_E": zetaE, "zeta_J": float(zetaJ),
            "p": float(p[-1]), "p_U": float(pU[-1]),
            "p_V": float(pV[-1]), "r_U": float(rU[-1]),
            "chi_U": float(sf["chi_U"][-1]),
            "chi_U_lower_limit_bound": float(np.max(sf["chi_U_lower_limit_bound"])),
            "radial_nodes": len(points), "source_mu_E": source["mu_E"],
            "source_zeta_E": source["zeta_E"]}


def clock_rate(x, rho, mu_plus, nu_minus, zeta_J):
    """Return d(V1/V0)/ds, including the background V displacement.

    V1 is the coefficient at fixed proper wall time. Dropping the V1/V0
    term in V1_s before forming this ratio gives a different equation.
    """
    response = wall_response_coefficients(x, rho, mu_plus, nu_minus)
    g = reference_coefficients(x)
    d = g["v"]-g["beta_plus"]
    zeta_J = _number(zeta_J, "Jordan areal-null lapse coefficient")
    E = zeta_J+(response["displacement_rate"]-response["delta_beta_plus"])/d
    return float(E/(2*A_RADIUS*d))


def _wall_trace(wedge, geometry, t, field_order, budget):
    V = min(geometry["Vb"], max(geometry["Va"],
            geometry["Va"]+(geometry["Vb"]-geometry["Va"])*t))
    U = float(wedge.geometry["wall"](V))
    bg = background(U, V)
    raw_x = float(bg["r"])
    if not .45-2e-13 <= raw_x <= .55+2e-13:
        raise ValueError("wall radius exits the admitted reference strip")
    x = min(.55, max(.45, raw_x))
    g = reference_coefficients(x)
    v, beta = float(g["v"]), float(g["beta_plus"])
    d = v-beta
    Vdot = -V/(2*A_RADIUS*d)
    if not Vdot > 0:
        raise ValueError("exterior proper-time clock lost its positive normalization")
    scalar = _field(wedge, np.array([U]), t, field_order, budget)
    if scalar["chi_U"][0] != 0:
        raise ValueError("the wall trace violates the no-return condition chi_U=0")
    w, j = float(scalar["chi"][0]), float(scalar["chi_V"][0])*Vdot
    material = {key: float(value) for key, value in field_mass_ingredients(x, w, j).items()}
    return {"U": U, "V": V, "x": x, "v": v, "beta_plus": beta,
            "beta_minus": float(g["beta_minus"]), "Vdot": Vdot,
            "ds_dt": (geometry["Vb"]-geometry["Va"])/Vdot,
            "w": w, "j": j, **material,
            "wall_physical_p_U_defect": float(scalar["p_U"][0]+w*bg["r_U"]/x**2),
            "reference_clock_normalization_defect":
                float(bg["C"]*wedge.geometry["wall"](V, 1)*Vdot**2-1)}


def _gauss_antiderivatives():
    nodes, _ = np.polynomial.legendre.leggauss(4)
    basis = []
    for i, point in enumerate(nodes):
        polynomial = Polynomial([1.])
        for j, other in enumerate(nodes):
            if i != j:
                polynomial *= Polynomial([-other, 1.])/(point-other)
        integral = polynomial.integ()
        basis.append(integral-integral(-1.))
    return nodes, tuple(basis)


class WallSolution:
    """Dense wall coefficients with a separately quadrature-controlled clock."""
    def __init__(self, wedge, geometry, cuts, segments, budget, settings):
        self.wedge, self.geometry = wedge, dict(geometry)
        self.cuts, self.segments, self.budget = cuts, segments, budget
        self.settings = settings
        self.radial_cache = {}
        self.maxima = dict.fromkeys(("mass_path_defect", "frame_mass_defect",
            "wall_lapse_frame_defect", "wall_physical_p_U_defect",
            "reference_clock_normalization_defect", "angular_jump_defect",
            "two_displacement_equations_defect", "clock_normalization_defect",
            "chi_U_lower_limit_bound"), 0.)

    def _wall(self, t):
        t = _time(t)
        self.budget.check()
        index = min(np.searchsorted(self.cuts, t, side="right")-1, len(self.segments)-1)
        Ip, Im, ratio = self.segments[index](t)
        row = _wall_trace(self.wedge, self.geometry, t, self.settings["field_order"], self.budget)
        mu, nu = Ip+row["Q_plus"], Im+row["Q_minus"]
        rho = row["v"]*ratio
        response = wall_response_coefficients(row["x"], rho, mu, nu)
        return {"t": t, **row, "mu_plus": float(mu), "nu_minus": float(nu),
                "rho": float(rho), "rho_over_v": float(ratio),
                "integral_L_plus": float(Ip), "integral_L_minus": float(Im),
                "rho_s": float(response["displacement_rate"]),
                **{key: float(value) for key, value in response.items()}}

    def _radial(self, t):
        t = float(t)
        if t not in self.radial_cache:
            keys = ("radial_order", "field_order", "source_order", "source_scale")
            self.radial_cache[t] = radial_response(self.wedge, self.geometry, t,
                budget=self.budget, **{key: self.settings[key] for key in keys})
        return self.radial_cache[t]

    def _diagnostics(self, wall, radial, eta=0.):
        x, v, bp = wall["x"], wall["v"], wall["beta_plus"]
        eta_s = clock_rate(x, wall["rho"], wall["mu_plus"], wall["nu_minus"], radial["zeta_J"])
        delta_d = wall["displacement_rate"]-wall["delta_beta_plus"]
        E = radial["zeta_J"]+delta_d/(v-bp)
        V1 = wall["V"]*eta
        V1_s = wall["Vdot"]*eta+wall["V"]*eta_s
        defects = {
            "mass_path_defect": radial["mu_J"]-wall["mu_plus"],
            "frame_mass_defect": radial["frame_mass"]-wall["Q_plus"],
            "wall_lapse_frame_defect": radial["zeta_J"]-radial["zeta_E"]-wall["p"]**2/4,
            "wall_physical_p_U_defect": wall["wall_physical_p_U_defect"],
            "reference_clock_normalization_defect": wall["reference_clock_normalization_defect"],
            "angular_jump_defect": wall["delta_beta_minus"]-wall["delta_beta_plus"]-K_ELL*wall["rho"],
            "two_displacement_equations_defect": wall["displacement_rate"]-wall["cap_side_displacement_rate"],
            "clock_normalization_defect": V1_s/wall["Vdot"]-eta+E,
            "chi_U_lower_limit_bound": radial["chi_U_lower_limit_bound"],
        }
        for key, value in defects.items():
            self.maxima[key] = max(self.maxima[key], abs(float(value)))
        return {**defects, "eta": float(eta), "eta_s": eta_s,
                "V1": float(V1), "V1_s": float(V1_s)}

    def sample(self, t):
        t = _time(t)
        wall, radial = self._wall(t), self._radial(t)
        panel = min(int(t*self.settings["time_panels"]), self.settings["time_panels"]-1)
        local = 2*(t-self.time_cuts[panel])*self.settings["time_panels"]-1
        weights = np.array([polynomial(local) for polynomial in self.clock_basis])
        eta = self.clock_prefix[panel]+float(weights @ self.clock_values[panel])/(2*self.settings["time_panels"])
        return {**wall, "radial_mu_E": radial["mu_E"], "radial_mu_J": radial["mu_J"],
                "zeta_E": radial["zeta_E"], "zeta_J": radial["zeta_J"],
                **self._diagnostics(wall, radial, eta)}

    def summary(self, times=tuple(np.linspace(0, 1, 9))):
        times = tuple(times)
        if not 1 <= len(times) <= 129:
            raise ValueError("display at most 129 wall time cuts")
        keys = ("t", "U", "V", "x", "beta_plus", "mu_plus", "nu_minus",
                "cap_total_mass", "rho", "rho_s", "displacement_rate", "delta_beta_plus",
                "delta_normal_acceleration", "radial_mu_E", "radial_mu_J",
                "zeta_E", "zeta_J", "eta", "eta_s", "V1", "V1_s", "mass_path_defect")
        rows = [{key: row[key] for key in keys} for row in (self.sample(t) for t in times)]
        return {"schema": SCHEMA,
                "domain": {key: self.geometry[key] for key in ("U0", "Ua", "Ub", "Va", "Vb")},
                **self.settings, "clock_time_order": 4,
                "maximum_absolute_defects": dict(self.maxima), "wall_samples": rows,
                "radial_evaluations": len(self.radial_cache),
                "rhs_calls": self.budget.rhs_calls, "metric_values": self.budget.metric_values,
                "source_nodes": self.budget.source_nodes,
                "process_cpu_seconds": time.process_time()-self.budget.started,
                "scope": "formal moving-wall and exterior V-clock coefficients; cap metric and central observer remain open"}


def solve_wall_evolution(wedge_field, geometry, *, radial_order=2, field_order=32,
                         time_panels=8, source_order=64, source_scale=1.,
                         tolerance=2e-9, budget=None):
    """Integrate quiet-anchored masses/displacement, then the exterior clock.

    Four-point Gauss panels integrate the clock. Its partial-panel values use
    the integrated cubic through those same rates. Vary time_panels separately
    from the field's first-derivative quadrature order and radial rule.
    """
    _domain(wedge_field, geometry)
    radial_order = _choice(radial_order, RADIAL_ORDERS, "radial order")
    field_order = _choice(field_order, (32, 64), "field derivative order")
    source_order = _choice(source_order, (32, 64), "source order")
    time_panels = _choice(time_panels, TIME_PANELS, "time panels")
    source_scale = _number(source_scale, "source scale")
    if source_scale not in (0., 1.):
        raise ValueError("source scale is the original 1 or zero-source control 0")
    tolerance = _number(tolerance, "wall tolerance")
    if not 1e-12 <= tolerance <= 1e-7:
        raise ValueError("wall tolerance must lie between 1e-12 and 1e-7")
    budget = _budget(wedge_field, budget)
    settings = dict(radial_order=radial_order, field_order=field_order,
                    time_panels=time_panels, source_order=source_order,
                    source_scale=source_scale, tolerance=tolerance)
    # Probe quiet/source compatibility before advancing any wall state.
    for t in (0., .5, 1.):
        radial_response(wedge_field, geometry, t, upper_U=geometry["U0"],
            budget=budget, **{key: settings[key] for key in
                            ("radial_order", "field_order", "source_order", "source_scale")})
    first = _wall_trace(wedge_field, geometry, 0., field_order, budget)
    if first["w"] != 0 or first["j"] != 0:
        raise ValueError("wall evolution requires the fixed quiet initial corner")

    def rate(t, state):
        budget.check(rhs=True)
        row = _wall_trace(wedge_field, geometry, float(t), field_order, budget)
        mu, nu = state[0]+row["Q_plus"], state[1]+row["Q_minus"]
        forcing = mu/row["x"]+row["beta_plus"]*(mu-nu)/(K_ELL*row["x"]**2)
        return row["ds_dt"]*np.array([row["L_plus"], row["L_minus"], forcing/row["v"]**2])

    cuts = np.unique(np.r_[0., .25, .75, 1., np.asarray(getattr(wedge_field, "cuts", []))])
    if cuts[0] != 0 or cuts[-1] != 1 or len(cuts) > 2*MAX_QUERY:
        raise ValueError("wall evolution received an invalid dense time partition")
    state, segments = np.zeros(3), []
    for lower, upper in zip(cuts[:-1], cuts[1:]):
        solved = solve_ivp(rate, (lower, upper), state, method="DOP853", rtol=tolerance,
            atol=1e-12, first_step=upper-lower, max_step=min(1/32, upper-lower), dense_output=True)
        if not solved.success or solved.t[-1] != upper:
            raise RuntimeError("wall integration did not reach its admitted endpoint")
        state = solved.y[:, -1]
        segments.append(solved.sol)
    result = WallSolution(wedge_field, geometry, cuts, segments, budget, settings)
    node, result.clock_basis = _gauss_antiderivatives()
    result.time_cuts = np.linspace(0., 1., time_panels+1)
    result.clock_values = np.zeros((time_panels, 4))
    for panel, lower in enumerate(result.time_cuts[:-1]):
        for index, point in enumerate(lower+(node+1)/(2*time_panels)):
            wall, radial = result._wall(float(point)), result._radial(float(point))
            diagnostic = result._diagnostics(wall, radial)
            result.clock_values[panel, index] = diagnostic["eta_s"]*wall["ds_dt"]
    weights = np.array([polynomial(1.) for polynomial in result.clock_basis])/(2*time_panels)
    result.clock_prefix = np.r_[0., np.cumsum(result.clock_values @ weights)]
    budget.check()
    return result
