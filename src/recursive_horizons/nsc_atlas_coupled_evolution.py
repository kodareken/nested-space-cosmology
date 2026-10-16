"""One common characteristic-exterior, moving-wall and polar-cap evolution.

The exterior stores the increment from its supplied source, on eta=Ubar/W(x).
The cap stores the weighted finite-volume variables of nsc_atlas_polar_cap.
First-order spatial upwinding and fixed-step SSPRK2 are a numerical prototype;
their residuals and refinements do not certify a continuum solution.
"""
from collections.abc import Mapping
import copy
import math

import numpy as np

from .nsc_atlas_characteristic_initialization import Budget
from .nsc_atlas_exterior_characteristic_slice import exterior_characteristic_slice
from .nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from .nsc_atlas_polar_cap import cap_metric, cap_rhs, cap_diagnostics


FIELDS = ("r", "L", "c")
WALL_FIELDS = (*FIELDS, "T", "s")


def _finite(value, name, shape=None):
    raw = np.asarray(value)
    if (raw.dtype.kind not in "iuf" or
            (isinstance(value, (list, tuple)) and any(
                isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat))):
        raise ValueError(name+" must be real and finite")
    result = np.array(raw, dtype=float, copy=True)
    if (not np.all(np.isfinite(result)) or
            (shape is not None and result.shape != shape)):
        raise ValueError(name+" has invalid shape or nonfinite entries")
    return result


def _grid(value, name, minimum, maximum):
    grid = _finite(value, name)
    if (grid.ndim != 1 or not minimum <= len(grid)-1 <= maximum
            or grid[0] != 0. or grid[-1] != 1.
            or not np.allclose(np.diff(grid), 1/(len(grid)-1), rtol=1e-12, atol=1e-14)):
        raise ValueError(name+" must be a uniform closed unit grid")
    return grid


def _state(state, ne, nc):
    if not isinstance(state, Mapping) or set(state) != {"increment", "wall", "cap"}:
        raise ValueError("state must contain increment, wall and cap")
    result = {}
    for group, names, shape in (("increment", FIELDS, (ne-1,)),
                                ("wall", WALL_FIELDS, ()),
                                ("cap", ("f", "Pi", "Q"), (nc,))):
        if not isinstance(state[group], Mapping) or set(state[group]) != set(names):
            raise ValueError("invalid "+group+" state fields")
        result[group] = {key: _finite(state[group][key], group+"."+key, shape) for key in names}
    return result


def canonical_field(c, r, alpha):
    """Normalized canonical field, with its exact zero-amplitude limit."""
    u = np.asarray(c)/np.asarray(r)
    z = np.sqrt(alpha/6)*u
    ratio = np.ones_like(z, dtype=float)
    np.divide(np.arcsinh(z), z, out=ratio, where=z != 0)
    return u*ratio


def quiet_initial_state(bottom, source, *, cap_cells, T=0.):
    """Assemble already evaluated quiet data; no source/geometry integration.

    bottom must include both endpoints on the initial uniform exterior grid.
    The initial wall clock T is supplied; s is elapsed Einstein proper time.
    """
    if type(cap_cells) is not int or not 3 <= cap_cells <= 128:
        raise ValueError("cap_cells must be an integer in [3,128]")
    r = _finite(bottom["r"], "bottom.r")
    if r.ndim != 1 or not 3 <= len(r) <= 513:
        raise ValueError("bottom requires 3..513 nodes")
    rows = {key: _finite(bottom["log_C" if key == "L" else key], key, r.shape) for key in FIELDS}
    if np.any(r <= 0) or np.any(rows["c"] != 0):
        raise ValueError("initial bottom must have positive radius and exactly quiet scalar")
    source_values = {key: float(source["log_C" if key == "L" else key]) for key in FIELDS}
    if any(abs(rows[key][0]-source_values[key]) > 5e-13 for key in FIELDS):
        raise ValueError("source and quiet bottom do not share their corner")
    return {"increment": {key: rows[key][1:-1]-source_values[key] for key in FIELDS},
            "wall": {**{key: float(rows[key][-1]) for key in FIELDS},
                     "T": float(_finite(T, "T", ())), "s": 0.},
            "cap": {key: np.zeros(cap_cells) for key in ("f", "Pi", "Q")}}


def coupled_stage(state, *, eta, y_faces, source, gauge, alpha, tube, budget,
                  exterior_quadrature_order=2, cap_quadrature_order=4,
                  diagnostics=False):
    """Evaluate all regional rates from one simultaneous state.

    source supplies r,log_C,c and their x derivatives at this same x.
    gauge supplies W,W_x,eta_plus,eta_plus_x for the fixed reference wall
    label W(x). The wall lapse is evolved; it is not prescribed again.
    Independent numerical incoming jets are compared with the hybrid map,
    never overwritten to make the coupling residual vanish.
    """
    if not isinstance(budget, Budget):
        raise ValueError("a shared caller Budget is required")
    budget.check()
    eta = _grid(eta, "eta", 2, 512)
    y = _grid(y_faces, "y_faces", 3, 128)
    ne, nc = len(eta)-1, len(y)-1
    state = _state(state, ne, nc)
    alpha = float(_finite(alpha, "alpha", ()))
    if alpha < 0:
        raise ValueError("alpha must be nonnegative")
    W, Wx, boost, boost_x = (float(_finite(gauge[k], "gauge."+k, ()))
                            for k in ("W", "W_x", "eta_plus", "eta_plus_x"))
    if W <= 0 or Wx <= 0 or not math.isclose(boost, -.5*math.log(Wx), rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("require positive fixed wall label W,W_x and its compatible boost")
    source_values = {key: float(_finite(source["log_C" if key == "L" else key], key, ())) for key in FIELDS}
    source_rates = {key: float(_finite(source[key], key, ())) for key in ("r_x", "log_C_x", "c_x")}
    wall = state["wall"]
    inc = {key: np.r_[0., state["increment"][key], wall[key]-source_values[key]] for key in FIELDS}
    X = {key: inc[key]+source_values[key] for key in FIELDS}
    # Negative moving-grid characteristic speed: information enters from
    # the wall, hence forward upwinding. Last backward jet is diagnostic.
    gradients = {key: np.r_[np.diff(inc[key]), inc[key][-1]-inc[key][-2]]/np.diff(eta)[0]
                 for key in FIELDS}
    jets = {"r_Ubar": gradients["r"]/W, "log_C_Ubar": gradients["L"]/W,
            "c_Ubar": gradients["c"]/W}
    ext = exterior_characteristic_slice(eta*W, X, jets, source_rates,
        alpha=alpha, quadrature_order=exterior_quadrature_order, tube=tube, budget=budget)
    metric = cap_metric(y, float(wall["r"]), state["cap"], alpha=alpha,
                        quadrature_order=cap_quadrature_order, budget=budget)
    outgoing = {"r_V": float(ext["r_x"][-1]), "c_V": float(ext["c_x"][-1]),
                "log_C_V": float(ext["log_C_x"][-1])}
    common = {"r": wall["r"], "c": wall["c"], "eta_plus": boost}
    cap_wall = {key: metric["wall"][key] for key in ("A", "delta", "w_minus")}
    C = math.exp(float(wall["L"]))
    first = hybrid_wall_boundary(common, outgoing, cap_wall,
                                {"C_wall": C, "log_C_s": 0.}, alpha=alpha)
    xs = float(first["Vdot_plus"])
    Lx = 2*boost_x-2*float(first["k_plus"])/xs+2*outgoing["log_C_V"]
    hybrid = hybrid_wall_boundary(common, outgoing, cap_wall,
                                 {"C_wall": C, "log_C_s": Lx*xs}, alpha=alpha)
    wall_rates = {"r": float(hybrid["r_s"])/xs, "c": float(hybrid["c_s"])/xs,
                  "L": Lx, "T": float(hybrid["T_s"])/xs, "s": 1/xs}
    wall_f = float(canonical_field(wall["c"], wall["r"], alpha))
    cap = cap_rhs(metric, T_x=wall_rates["T"], radius_x=wall_rates["r"],
                  incoming_w_plus=float(hybrid["w_plus"]), wall_f=wall_f, budget=budget)
    longitudinal = {"r": ext["r_x"], "L": ext["log_C_x"], "c": ext["c_x"]}
    mapped_source = {"r": source_rates["r_x"], "L": source_rates["log_C_x"], "c": source_rates["c_x"]}
    rates = {"increment": {key: (longitudinal[key]-mapped_source[key]
                              +eta*(Wx/W)*gradients[key])[1:-1] for key in FIELDS},
             "wall": wall_rates, "cap": cap["rates"]}
    wall_residual = {key: wall_rates[key]-(longitudinal[key][-1]+Wx*gradients[key][-1]/W)
                     for key in FIELDS}
    budget.check(2)
    result = {"rates": rates, "exterior": X, "exterior_first_jets": jets,
              "exterior_longitudinal": longitudinal, "cap_metric": metric,
              "hybrid": hybrid, "wall_f": wall_f,
              "wall_kinematic_defect": wall_residual,
              "wall_value_cell_difference": cap["wall_value_cell_difference"],
              "courant_rate": max(float(Wx/W/np.diff(eta)[0]), cap["characteristic_courant_rate"]),
              "source_constraint": (None if "source_constraint" not in source
                                    else float(source["source_constraint"]))}
    if diagnostics:
        result["cap_diagnostics"] = cap_diagnostics(metric, budget=budget)
    return result


def _combine(a, b, wa, wb):
    return {group: {key: wa*a[group][key]+wb*b[group][key] for key in a[group]} for group in a}


class CoupledEvolutionError(RuntimeError):
    """A stopped fixed-grid calculation with its original stage context."""
    def __init__(self, diagnostic):
        self.diagnostic = diagnostic
        super().__init__(diagnostic["reason"])


def evolve_coupled(initial, x_nodes, *, eta, y_faces, source, gauge, alpha,
                   tube, budget, save_every=1, maximum_courant=.4,
                   exterior_quadrature_order=2, cap_quadrature_order=4):
    """Fixed-step SSPRK2 with no retry, step change or source replacement.

    source(x) and gauge(x) return simultaneous rows and must share the
    caller's accounting. A caller may supply authenticated/cached rows;
    caching and its construction costs must be explicit in its record.
    Snapshots save accepted states, not internal predictor stages.
    """
    if not isinstance(budget, Budget):
        raise ValueError("a shared caller Budget is required")
    budget.check()
    x = _finite(x_nodes, "x_nodes")
    if (x.ndim != 1 or not 2 <= len(x) <= 2049 or np.any(np.diff(x) <= 0)
            or not np.allclose(np.diff(x), np.diff(x)[0], rtol=1e-10, atol=1e-14)):
        raise ValueError("use 1..2048 fixed uniform positive time steps")
    eta = _grid(eta, "eta", 2, 512)
    y_faces = _grid(y_faces, "y_faces", 3, 128)
    state = _state(initial, len(eta)-1, len(y_faces)-1)
    if type(save_every) is not int or save_every < 1 or (len(x)-1) % save_every:
        raise ValueError("save_every must divide the number of steps")
    if not math.isfinite(maximum_courant) or not 0 < maximum_courant <= .5:
        raise ValueError("require 0<maximum_courant<=.5")
    stages, max_courant = 0, 0.
    saved = [{"x": float(x[0]), "state": copy.deepcopy(state)}]
    phase, step = "initial", 0
    attempt_state = state

    def evaluate(value, at, dt, diagnostics=False):
        nonlocal stages, max_courant, attempt_state
        attempt_state = value
        stage = coupled_stage(value, eta=eta, y_faces=y_faces, source=source(float(at)),
            gauge=gauge(float(at)), alpha=alpha, tube=tube, budget=budget,
            exterior_quadrature_order=exterior_quadrature_order,
            cap_quadrature_order=cap_quadrature_order, diagnostics=diagnostics)
        stages += 1
        courant = dt*stage["courant_rate"]
        max_courant = max(max_courant, courant)
        if courant > maximum_courant:
            raise ValueError("declared characteristic Courant limit exceeded")
        return stage

    try:
        for step in range(len(x)-1):
            dt = float(x[step+1]-x[step])
            phase = "first_stage"
            first = evaluate(state, x[step], dt)
            predictor = _combine(state, first["rates"], 1., dt)
            phase = "predictor_stage"
            second = evaluate(predictor, x[step+1], dt)
            corrected = _combine(predictor, second["rates"], 1., dt)
            state = _combine(state, corrected, .5, .5)
            if (step+1) % save_every == 0:
                saved.append({"x": float(x[step+1]), "state": copy.deepcopy(state)})
        phase = "final_diagnostics"
        final = evaluate(state, x[-1], 0., diagnostics=True)
    except (ValueError, RuntimeError, FloatingPointError, OverflowError) as error:
        raise CoupledEvolutionError({"phase": phase, "step": step,
            "x_left": float(x[step]), "x_right": float(x[min(step+1, len(x)-1)]),
            "reason": str(error), "completed_stages": stages,
            "charged_evaluations": budget.evaluations, "max_courant": max_courant,
            "accepted_state": copy.deepcopy(state),
            "attempt_state": copy.deepcopy(attempt_state)}) from error
    budget.check()
    return {"snapshots": saved, "final_stage": final, "stages": stages,
            "steps": len(x)-1, "max_courant": max_courant,
            "spatial_scheme": "first-order exterior upwind and regular cap finite volume",
            "time_scheme": "fixed-step SSPRK2", "retries": 0}
