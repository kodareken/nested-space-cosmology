"""Independent pure-null constraints from saved moving-grid histories.

Only centred finite differences of saved increments and supplied source
jets enter these diagnostics.  No evolution RHS, prefix integral, hybrid
incoming jet, or solved mixed equation is used.  Values are sampled defects,
not a uniform residual bound or a constraint-propagation certificate.
"""
from collections.abc import Mapping

import numpy as np

FIELDS = ("r", "L", "c")
SOURCE_KEYS = ("r", "r_x", "r_xx", "log_C", "log_C_x", "c", "c_x")
GAUGE_KEYS = ("W", "W_x", "W_xx")
MAX_VALUES = 131073
OVERLAP_TIME_TOLERANCE = 1e-14


def _array(value, name, shape=None):
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf" or raw.size == 0:
            raise ValueError
        if isinstance(value, (list, tuple)) and any(
                isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat):
            raise ValueError
        result = np.array(raw, dtype=float, copy=True)
        if not np.all(np.isfinite(result)) or (shape is not None and result.shape != shape):
            raise ValueError
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(name + " must contain finite real values with the required shape") from error
    return result


def _table(value, keys, name, shape):
    if not isinstance(value, Mapping) or not set(keys) <= set(value):
        raise ValueError(name + " lacks required fields: " + ", ".join(keys))
    return {key: _array(value[key], name+"."+key, shape) for key in keys}


def _uniform(value, name, minimum, maximum):
    grid = _array(value, name)
    if (grid.ndim != 1 or not minimum <= len(grid) <= maximum
            or np.any(np.diff(grid) <= 0)
            or not np.allclose(np.diff(grid), np.diff(grid)[0], rtol=1e-10, atol=1e-14)):
        raise ValueError(name + " must be a bounded uniform increasing grid")
    return grid


def _summary(values, mask):
    selected = values[mask]
    if not len(selected):
        return {"count": 0, "min": None, "max": None, "max_abs": None, "rms": None}
    # Scale the square to avoid overflow when a finite diagnostic is large.
    largest = float(np.max(abs(selected)))
    rms = largest*float(np.sqrt(np.mean((selected/largest)**2))) if largest else 0.
    return {"count": int(len(selected)), "min": float(np.min(selected)),
            "max": float(np.max(selected)), "max_abs": largest, "rms": rms}


def coupled_history_diagnostics(snapshots, *, eta, gauge, source, alpha,
                                tube, budget, source_fronts=()):
    """Differentiate saved exterior states at fixed Ubar, excluding all edges.

    snapshots is the coupled evolution's [{x,state}] sequence.  Its saved
    state has interior `increment` arrays and scalar `wall` values r,L,c.
    Cap fields/clocks are retained by the caller but not used here.  eta
    includes both endpoints 0,1.  Both grids need at least four nodes.

    gauge supplies arrays W,W_x,W_xx at the saved x; source supplies arrays
    r,r_x,r_xx,log_C,log_C_x,c,c_x from one consistent source representation.
    They are supplied data, not inferred from the regional finite
    differences.  The source's own pure constraint is returned separately.

    With Y=X-X_source, b=W_x/W and a=eta*b, the fixed-Ubar derivative is
    Dx-a*Deta.  Its second application includes
    eta*(2*b*b-W_xx/W)*Y_eta.  All partial differences are central second
    order, including the four-corner mixed difference.  No stencil touches
    beyond the saved grid.  `front_near` includes a time stencil intersecting
    either original constant-x source front; boundaries are never included.

    tube gives engineering bounds r,L,c, not a proved invariant set.  The
    caller's live cumulative budget is mandatory; 3*Nx*Neta assembly values,
    9*(Nx-2)*(Neta-2) derivative points, and Nx source checks are charged.
    """
    if not callable(getattr(budget, "check", None)):
        raise ValueError("a caller-owned cumulative budget is required")
    budget.check()
    if not isinstance(snapshots, (list, tuple)) or not 4 <= len(snapshots) <= 2049:
        raise ValueError("snapshots must contain 4..2049 saved rows")
    try:
        x = _uniform([row["x"] for row in snapshots], "saved x", 4, 2049)
    except (KeyError, TypeError) as error:
        raise ValueError("each snapshot must supply x and state") from error
    eta = _uniform(eta, "eta", 4, 513)
    if eta[0] != 0 or eta[-1] != 1:
        raise ValueError("eta must include endpoints 0 and 1")
    nt, ne = len(x), len(eta)
    if nt*ne > MAX_VALUES:
        raise ValueError("saved history exceeds the bounded grid-value allowance")
    alpha = float(_array(alpha, "alpha", ()))
    if alpha < 0:
        raise ValueError("alpha must be nonnegative")
    gauge = _table(gauge, GAUGE_KEYS, "gauge", (nt,))
    if np.any(gauge["W"] <= 0) or np.any(gauge["W_x"] <= 0):
        raise ValueError("the original wall gauge requires W,W_x>0")
    source = _table(source, SOURCE_KEYS, "source", (nt,))
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
    values = {"r": source["r"], "L": source["log_C"], "c": source["c"]}
    Y = {key: np.zeros((nt, ne)) for key in FIELDS}
    for i, row in enumerate(snapshots):
        try:
            state = row["state"]
            if not isinstance(state, Mapping):
                raise ValueError
            increments = _table(state["increment"], FIELDS, "increment", (ne-2,))
            wall = _table(state["wall"], FIELDS, "wall", ())
        except (KeyError, TypeError) as error:
            raise ValueError("snapshot state must supply increment and wall") from error
        for key in FIELDS:
            Y[key][i, 1:-1] = increments[key]
            Y[key][i, -1] = wall[key]-values[key][i]
    X = {key: Y[key]+values[key][:, None] for key in FIELDS}
    margins = {}
    for key in FIELDS:
        low, high = bounds[key]
        minimum, maximum = float(np.min(X[key])), float(np.max(X[key]))
        if minimum < low or maximum > high:
            raise ValueError("saved state exits the declared tube: "+key)
        margins[key] = {"min": minimum, "max": maximum,
                        "lower_margin": minimum-low, "upper_margin": high-maximum}
    budget.check(3*nt*ne)
    dx, de = x[1]-x[0], eta[1]-eta[0]
    W = gauge["W"][1:-1, None]
    b = gauge["W_x"][1:-1, None]/W
    at = eta[None, 1:-1]*b
    second_map = eta[None, 1:-1]*(2*b*b-gauge["W_xx"][1:-1, None]/W)
    derivatives = {}
    first_source = {"r": source["r_x"], "L": source["log_C_x"], "c": source["c_x"]}
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            for key in FIELDS:
                field = Y[key]
                d_eta = (field[1:-1, 2:]-field[1:-1, :-2])/(2*de)
                d_x = (field[2:, 1:-1]-field[:-2, 1:-1])/(2*dx)
                derivatives[key+"_Ubar"] = d_eta/W
                derivatives[key+"_x"] = first_source[key][1:-1, None]+d_x-at*d_eta
                if key == "r":
                    d_etaeta = (field[1:-1, 2:]-2*field[1:-1, 1:-1]+field[1:-1, :-2])/de**2
                    d_xx = (field[2:, 1:-1]-2*field[1:-1, 1:-1]+field[:-2, 1:-1])/dx**2
                    d_xeta = (field[2:, 2:]-field[2:, :-2]-field[:-2, 2:]+field[:-2, :-2])/(4*dx*de)
                    derivatives["r_Ubar_Ubar"] = d_etaeta/W**2
                    derivatives["r_xx"] = (source["r_xx"][1:-1, None]+d_xx-2*at*d_xeta
                                             +at*at*d_etaeta+second_map*d_eta)
            r, L, c = (X[key][1:-1, 1:-1] for key in FIELDS)
            u, C = c/r, np.exp(L)
            D = 1+alpha*u*u/6
            kinetic = {}
            for coordinate in ("Ubar", "x"):
                u_a = (derivatives["c_"+coordinate]-u*derivatives["r_"+coordinate])/r
                kinetic[coordinate] = alpha*r*u_a*u_a/(2*D)
            E_U = (derivatives["r_Ubar_Ubar"]-derivatives["L_Ubar"]*derivatives["r_Ubar"]
                   +kinetic["Ubar"])
            E_x = derivatives["r_xx"]-derivatives["L_x"]*derivatives["r_x"]+kinetic["x"]
            mass = r*(1+4*derivatives["r_Ubar"]*derivatives["r_x"]/C)/2
            us = source["c"]/source["r"]
            source_u_x = (source["c_x"]-us*source["r_x"])/source["r"]
            source_constraint = (source["r_xx"]-source["log_C_x"]*source["r_x"]
                                 +alpha*source["r"]*source_u_x**2/(2*(1+alpha*us**2/6)))
    except FloatingPointError as error:
        raise ValueError("nonfinite saved-history diagnostic arithmetic") from error
    arrays = {"E_Ubar": E_U, "E_x": E_x, "mass": mass}
    if not all(np.all(np.isfinite(value)) for value in (*arrays.values(), *derivatives.values(), source_constraint)):
        raise ValueError("nonfinite saved-history diagnostics")
    near_rows = np.zeros(nt-2, dtype=bool)
    for front in fronts:
        near_rows |= (x[:-2] <= front) & (front <= x[2:])
    near = np.broadcast_to(near_rows[:, None], E_U.shape).copy()
    masks = {"all_interior": np.ones(E_U.shape, dtype=bool),
             "front_near": near, "away_from_front": ~near}
    budget.check(9*(nt-2)*(ne-2)+nt)
    summaries = {key: {name: _summary(value, mask) for name, mask in masks.items()}
                 for key, value in arrays.items()}
    return {
        "x": x, "eta": eta, "gauge": gauge, "source_fronts": fronts,
        "states": X, "increment": Y, "interior_x": x[1:-1],
        "interior_eta": eta[1:-1], "interior_Ubar": W*eta[None, 1:-1],
        **arrays, "derivatives": derivatives, "masks": masks,
        "summaries": summaries, "source_constraint": source_constraint,
        "source_constraint_max_abs": float(np.max(abs(source_constraint))),
        "tube_margins": margins,
        "settings": {"difference_order": 2, "boundary_stencils": "excluded",
                     "time_nodes": nt, "eta_nodes": ne,
                     "charged_evaluations": 3*nt*ne+9*(nt-2)*(ne-2)+nt},
        "scope": {"source_integrations": 0, "time_steps": 0,
                  "mixed_equations_used": False, "cap_constraint_evaluated": False,
                  "source_error_removed": False, "residual_error_enclosed": False},
    }


def compare_rectangle_overlap(states, x, eta, W, rectangle_x, rectangle_Ubar,
                              rectangle_states, *, budget):
    """Compare saved fields on the old rectangle with no time interpolation.

    Current state arrays have axes (x,eta); old arrays and returned
    differences have axes (Ubar,x).  Each old time must match one unique
    saved time, exactly or within the fixed absolute 1e-14 coordinate
    roundoff tolerance.  The actual mismatch is returned.  No time
    interpolation, source call, or regional equation is evaluated.

    At each match, Ubar/W is interpolated linearly on eta and separately
    with a three-node Lagrange polynomial centred on its nearest node.
    Both rules interpolate X directly; subtracting a common source value
    gives the same rules on Y because their weights sum to one.  The
    difference of the rules measures interpolation sensitivity, not its
    error.  Coordinates at most eight machine epsilons outside a closed
    eta endpoint are clamped as roundoff, with the displacement reported;
    all larger excursions fail.  No physical extrapolation is allowed.

    Grids have at most 513 nodes on each axis.  The caller's cumulative
    budget is charged 15 values per target (three fields, two linear and
    three quadratic samples), before any interpolation is performed.
    """
    if not callable(getattr(budget, "check", None)):
        raise ValueError("a caller-owned cumulative budget is required")
    budget.check()
    x = _uniform(x, "current x", 3, 513)
    eta = _uniform(eta, "current eta", 3, 513)
    old_x = _uniform(rectangle_x, "rectangle x", 3, 513)
    old_U = _uniform(rectangle_Ubar, "rectangle Ubar", 3, 513)
    if eta[0] != 0 or eta[-1] != 1 or old_U[0] < 0:
        raise ValueError("require eta endpoints 0,1 and nonnegative rectangle Ubar")
    W = _array(W, "W", (len(x),))
    if np.any(W <= 0):
        raise ValueError("W must remain positive")
    current = _table(states, FIELDS, "states", (len(x), len(eta)))
    old = _table(rectangle_states, FIELDS, "rectangle_states", (len(old_U), len(old_x)))
    if np.any(current["r"] <= 0) or np.any(old["r"] <= 0):
        raise ValueError("overlap radii must remain positive")
    indices = []
    for value in old_x:
        right = int(np.searchsorted(x, value))
        choices = [j for j in (right-1, right) if 0 <= j < len(x)]
        matched = [j for j in choices if abs(x[j]-value) <= OVERLAP_TIME_TOLERANCE]
        if len(matched) != 1:
            raise ValueError("rectangle time has no unique saved match; time interpolation is forbidden")
        indices.append(matched[0])
    indices = np.array(indices, dtype=int)
    if np.any(np.diff(indices) <= 0):
        raise ValueError("rectangle times must match distinct increasing current times")
    eta_requested = old_U[:, None]/W[indices][None, :]
    epsilon = 8*np.finfo(float).eps
    if np.any(eta_requested < -epsilon) or np.any(eta_requested > 1+epsilon):
        raise ValueError("rectangle exits the current retained domain; extrapolation is forbidden")
    eta_used = np.clip(eta_requested, 0., 1.)
    count = len(old_U)*len(old_x)
    budget.check(15*count)
    linear = {key: np.empty((len(old_U), len(old_x))) for key in FIELDS}
    quadratic = {key: np.empty_like(linear[key]) for key in FIELDS}
    stencil_starts = np.empty(eta_used.shape, dtype=int)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            for j, index in enumerate(indices):
                q = eta_used[:, j]
                right = np.clip(np.searchsorted(eta, q), 1, len(eta)-1)
                left = right-1
                nearest = np.where(q-eta[left] <= eta[right]-q, left, right)
                starts = np.clip(nearest-1, 0, len(eta)-3)
                stencil_starts[:, j] = starts
                stencil = starts[:, None]+np.arange(3)
                nodes = eta[stencil]
                weights = np.ones_like(nodes)
                for k in range(3):
                    for m in range(3):
                        if m != k:
                            weights[:, k] *= (q-nodes[:, m])/(nodes[:, k]-nodes[:, m])
                for key in FIELDS:
                    linear[key][:, j] = np.interp(q, eta, current[key][index])
                    quadratic[key][:, j] = np.sum(weights*current[key][index, stencil], axis=1)
            dl = {key: linear[key]-old[key] for key in FIELDS}
            dq = {key: quadratic[key]-old[key] for key in FIELDS}
            sensitivity = {key: quadratic[key]-linear[key] for key in FIELDS}
    except FloatingPointError as error:
        raise ValueError("nonfinite overlap interpolation arithmetic") from error
    if not all(np.all(np.isfinite(value)) for group in (dl, dq, sensitivity) for value in group.values()):
        raise ValueError("nonfinite overlap differences")
    full = np.ones(eta_used.shape, dtype=bool)
    interior = np.zeros_like(full)
    interior[1:-1, 1:-1] = True
    summaries = {key: {
        method: {"full": _summary(values[key], full), "interior": _summary(values[key], interior)}
        for method, values in (("linear", dl), ("quadratic", dq),
                               ("quadratic_minus_linear", sensitivity))} for key in FIELDS}
    budget.check()
    return {
        "rectangle_x": old_x, "rectangle_Ubar": old_U,
        "matched_time_indices": indices, "time_coordinate_difference": x[indices]-old_x,
        "max_abs_time_coordinate_difference": float(np.max(abs(x[indices]-old_x))),
        "eta_requested": eta_requested, "eta_used": eta_used,
        "max_eta_endpoint_roundoff_adjustment": float(np.max(abs(eta_requested-eta_used))),
        "quadratic_stencil_starts": stencil_starts,
        "difference_linear": dl, "difference_quadratic": dq,
        "interpolation_sensitivity": sensitivity, "summaries": summaries,
        "settings": {"time_coordinate_tolerance": OVERLAP_TIME_TOLERANCE,
                     "eta_endpoint_tolerance": epsilon, "charged_evaluations": 15*count,
                     "output_axes": ["rectangle_Ubar", "rectangle_x"]},
        "scope": {"time_interpolations": 0, "source_integrations": 0,
                  "time_steps": 0, "physical_extrapolation": False,
                  "interpolation_error_enclosed": False},
    }
