"""Independent pure constraints from saved departures and exact vacuum jets.

Only centred differences of saved departures are numerical derivatives.
Analytic Schwarzschild jets and the supplied source's actual jets are added
afterwards. No mixed evolution equation, prefix or incoming wall law is used.
The resulting sampled defects are not continuum error enclosures.
"""
from collections.abc import Mapping

import numpy as np

from .nsc_atlas_coupled_diagnostics import (
    FIELDS, SOURCE_KEYS, GAUGE_KEYS, _array, _table, _uniform, _summary,
)
from .nsc_atlas_vacuum_departure import vacuum_slice

MAX_TIME_NODES = 129
MAX_ETA_NODES = 513


def departure_history_diagnostics(snapshots, *, eta, gauge, source, alpha,
                                  tube, budget, source_fronts=()):
    """Evaluate saved-history constraints at fixed Ubar, excluding all edges.

    Each snapshot supplies x, state.departure interior arrays r,L,c and
    physical.wall values r,L,c. The latter recover the terminal departure;
    the source endpoint has zero departure. The actual field is
    X=X_vacuum+(X_source-X_vacuum_at_source)+Z.

    gauge supplies W,W_x,W_xx arrays; source supplies r,r_x,r_xx,log_C,
    log_C_x,c,c_x arrays from one consistent source representation. Its
    curvature is not replaced by a source ODE or inferred from its samples.
    Uniform histories have 4..129 times and 4..513 eta nodes. Each analytic
    vacuum query is one time row of at most 513 points.

    Charges are 4*Nt*Ne+9*(Nt-2)*(Ne-2)+Nt: one analytic vacuum point and
    three field assembly values per node, nine centred derivative points
    per interior node, and one source check per time. All work shares the
    supplied live budget. Boundary stencils are excluded, not extrapolated.
    """
    if not callable(getattr(budget, "check", None)):
        raise ValueError("a caller-owned cumulative budget is required")
    budget.check()
    if not isinstance(snapshots, (list, tuple)) or not 4 <= len(snapshots) <= MAX_TIME_NODES:
        raise ValueError("snapshots must contain 4..129 saved rows")
    try:
        x = _uniform([row["x"] for row in snapshots], "saved x", 4, MAX_TIME_NODES)
    except (KeyError, TypeError) as error:
        raise ValueError("each snapshot must supply x, state and physical wall") from error
    eta = _uniform(eta, "eta", 4, MAX_ETA_NODES)
    if eta[0] != 0 or eta[-1] != 1:
        raise ValueError("eta must include endpoints 0 and 1")
    nt, ne = len(x), len(eta)
    alpha = float(_array(alpha, "alpha", ()))
    if alpha < 0:
        raise ValueError("alpha must be nonnegative")
    gauge = _table(gauge, GAUGE_KEYS, "gauge", (nt,))
    source = _table(source, SOURCE_KEYS, "source", (nt,))
    if np.any(gauge["W"] <= 0) or np.any(gauge["W_x"] <= 0):
        raise ValueError("the original wall gauge requires W,W_x>0")
    if not isinstance(tube, Mapping) or set(tube) != set(FIELDS):
        raise ValueError("tube must contain exactly r,L,c bounds")
    bounds = {key: _array(tube[key], "tube."+key, (2,)) for key in FIELDS}
    if any(not low < high for low, high in bounds.values()) or bounds["r"][0] <= 0:
        raise ValueError("tube bounds must increase and radius floor must be positive")
    with np.errstate(over="ignore", under="ignore"):
        lapse_bounds = np.exp(bounds["L"])
    if np.any(lapse_bounds <= 0) or not np.all(np.isfinite(lapse_bounds)):
        raise ValueError("tube lapse must remain positive and representable")
    fronts = np.asarray(source_fronts)
    if fronts.size:
        fronts = _array(source_fronts, "source_fronts")
        if fronts.ndim != 1 or len(fronts) > 2 or np.any(np.diff(fronts) <= 0):
            raise ValueError("source_fronts must contain at most two ordered fixed-x labels")
    else:
        fronts = np.empty(0)

    source_values = {"r": source["r"], "L": source["log_C"], "c": source["c"]}
    Z = {key: np.zeros((nt, ne)) for key in FIELDS}
    walls = {key: np.empty(nt) for key in FIELDS}
    for i, row in enumerate(snapshots):
        try:
            interiors = _table(row["state"]["departure"], FIELDS, "departure", (ne-2,))
            wall = _table(row["physical"]["wall"], FIELDS, "physical.wall", ())
        except (KeyError, TypeError) as error:
            raise ValueError("snapshot must supply state.departure and physical.wall") from error
        for key in FIELDS:
            Z[key][i, 1:-1] = interiors[key]
            walls[key][i] = wall[key]
    vacuum_keys = ("r", "L", "r_Ubar", "r_x", "r_Ubar_Ubar", "r_xx",
                   "log_C_Ubar", "log_C_x")
    vacuum = {key: np.empty((nt, ne)) for key in vacuum_keys}
    for i, at in enumerate(x):
        row = vacuum_slice(eta*gauge["W"][i], float(at), budget=budget)
        for key in vacuum:
            vacuum[key][i] = row[key]
    background = {"r": vacuum["r"], "L": vacuum["L"], "c": np.zeros((nt, ne))}
    offsets = {"r": source["r"]-x, "L": source["log_C"], "c": source["c"]}
    X = {}
    with np.errstate(over="ignore", invalid="ignore"):
        for key in FIELDS:
            Z[key][:, -1] = walls[key]-background[key][:, -1]-offsets[key]
            X[key] = background[key]+offsets[key][:, None]+Z[key]
            X[key][:, 0], X[key][:, -1] = source_values[key], walls[key]
    if not all(np.all(np.isfinite(v)) for v in (*X.values(), *Z.values())):
        raise ValueError("nonfinite saved departure reconstruction")
    margins = {}
    for key in FIELDS:
        low, high = bounds[key]
        minimum, maximum = float(np.min(X[key])), float(np.max(X[key]))
        if minimum < low or maximum > high:
            raise ValueError("saved state exits the declared tube: "+key)
        margins[key] = {"min": minimum, "max": maximum,
                        "lower_margin": minimum-low, "upper_margin": high-maximum}
    budget.check(3*nt*ne)

    take = np.s_[1:-1, 1:-1]
    dx, de = x[1]-x[0], eta[1]-eta[0]
    W = gauge["W"][1:-1, None]
    b = gauge["W_x"][1:-1, None]/W
    transport = eta[None, 1:-1]*b
    second_map = eta[None, 1:-1]*(2*b*b-gauge["W_xx"][1:-1, None]/W)
    first_offsets = {"r": source["r_x"]-1., "L": source["log_C_x"], "c": source["c_x"]}
    derivatives = {}
    budget.check(9*(nt-2)*(ne-2)+nt)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            for key in FIELDS:
                field = Z[key]
                d_eta = (field[1:-1, 2:]-field[1:-1, :-2])/(2*de)
                d_x = (field[2:, 1:-1]-field[:-2, 1:-1])/(2*dx)
                vu = vacuum[("log_C" if key == "L" else key)+"_Ubar"][take] if key != "c" else 0.
                vx = vacuum[("log_C" if key == "L" else key)+"_x"][take] if key != "c" else 0.
                derivatives[key+"_Ubar"] = vu+d_eta/W
                derivatives[key+"_x"] = vx+first_offsets[key][1:-1, None]+d_x-transport*d_eta
                if key == "r":
                    d_etaeta = (field[1:-1, 2:]-2*field[1:-1, 1:-1]+field[1:-1, :-2])/de**2
                    d_xx = (field[2:, 1:-1]-2*field[1:-1, 1:-1]+field[:-2, 1:-1])/dx**2
                    d_xeta = (field[2:, 2:]-field[2:, :-2]-field[:-2, 2:]+field[:-2, :-2])/(4*dx*de)
                    derivatives["r_Ubar_Ubar"] = vacuum["r_Ubar_Ubar"][take]+d_etaeta/W**2
                    derivatives["r_xx"] = (vacuum["r_xx"][take]+source["r_xx"][1:-1, None]
                        +d_xx-2*transport*d_xeta+transport**2*d_etaeta+second_map*d_eta)
            r, L, c = (X[key][take] for key in FIELDS)
            u, C = c/r, np.exp(L)
            D = 1+alpha*u*u/6
            kinetic = {}
            for coordinate in ("Ubar", "x"):
                u_a = (derivatives["c_"+coordinate]-u*derivatives["r_"+coordinate])/r
                kinetic[coordinate] = alpha*r*u_a*u_a/(2*D)
            E_U = derivatives["r_Ubar_Ubar"]-derivatives["L_Ubar"]*derivatives["r_Ubar"]+kinetic["Ubar"]
            E_x = derivatives["r_xx"]-derivatives["L_x"]*derivatives["r_x"]+kinetic["x"]
            mass = r*(1+4*derivatives["r_Ubar"]*derivatives["r_x"]/C)/2
            us = source["c"]/source["r"]
            us_x = (source["c_x"]-us*source["r_x"])/source["r"]
            source_constraint = (source["r_xx"]-source["log_C_x"]*source["r_x"]
                +alpha*source["r"]*us_x**2/(2*(1+alpha*us**2/6)))
    except FloatingPointError as error:
        raise ValueError("nonfinite saved-departure diagnostic arithmetic") from error
    arrays = {"E_Ubar": E_U, "E_x": E_x, "mass": mass}
    if not all(np.all(np.isfinite(v)) for v in (*arrays.values(), *derivatives.values(), source_constraint)):
        raise ValueError("nonfinite saved-departure diagnostics")
    near_rows = np.zeros(nt-2, dtype=bool)
    for front in fronts:
        near_rows |= (x[:-2] <= front) & (front <= x[2:])
    near = np.broadcast_to(near_rows[:, None], E_U.shape).copy()
    masks = {"all_interior": np.ones(E_U.shape, dtype=bool), "front_near": near, "away_from_front": ~near}
    summaries = {key: {name: _summary(value, mask) for name, mask in masks.items()}
                 for key, value in arrays.items()}
    budget.check()
    return {"x": x, "eta": eta, "gauge": gauge, "source_fronts": fronts,
        "states": X, "departure": Z, "interior_x": x[1:-1], "interior_eta": eta[1:-1],
        "interior_Ubar": W*eta[None, 1:-1], **arrays, "derivatives": derivatives,
        "masks": masks, "summaries": summaries, "source_constraint": source_constraint,
        "source_constraint_max_abs": float(np.max(abs(source_constraint))), "tube_margins": margins,
        "settings": {"difference_order": 2, "boundary_stencils": "excluded",
            "time_nodes": nt, "eta_nodes": ne, "analytic_vacuum_jets": True,
            "charged_evaluations": 4*nt*ne+9*(nt-2)*(ne-2)+nt},
        "scope": {"source_integrations": 0, "time_steps": 0, "mixed_equations_used": False,
            "cap_constraint_evaluated": False, "source_error_removed": False,
            "residual_error_enclosed": False, "vacuum_background_differenced": False}}
