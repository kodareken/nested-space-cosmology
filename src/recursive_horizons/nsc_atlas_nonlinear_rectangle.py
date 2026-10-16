"""Finite exterior characteristic march and independent nodal diagnostics.

Axes are (Ubar,x). Only supplied r,L,c boundary values enter the march;
this module constructs no source, changes no gauge and selects no amplitude.
Cell algebraic residuals, finite-difference constraints and continuum error
are distinct quantities. None of the reported differences is an enclosure.
"""
from collections.abc import Mapping
import copy

import numpy as np

from .nsc_atlas_characteristic_cell import CellSolveError, solve_characteristic_cell

NAMES = ("r", "L", "c")
JET_NAMES = ("r_Ubar", "r_x", "log_C_Ubar", "log_C_x", "c_Ubar", "c_x")
MAX_AXIS_POINTS = 65
CORNER_TOLERANCE = 5e-13
CELL_OPTIONS = {"atol", "rtol", "reference_scales", "max_updates", "max_backtracks", "maximum_condition"}


class RectangleSolveError(RuntimeError):
    """A stopped cell with its supplied corners and cumulative counts."""

    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        self.reason = diagnostics["reason"]
        super().__init__(f"rectangle cell {diagnostics['indices']}: {self.reason}")


def _array(value, label, shape=None):
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf" or raw.size == 0:
            raise ValueError
        if isinstance(value, (list, tuple)) and any(
                isinstance(v, (bool, np.bool_)) for v in np.asarray(value, dtype=object).flat):
            raise ValueError
        out = np.array(raw, dtype=float, copy=True)
        if not np.all(np.isfinite(out)) or (shape is not None and out.shape != shape):
            raise ValueError
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(label+" must contain finite real values with the required shape") from error
    return out


def _axis(value, label, minimum=2):
    axis = _array(value, label)
    if axis.ndim != 1 or not minimum <= len(axis) <= MAX_AXIS_POINTS:
        raise ValueError(f"{label} must have {minimum}..{MAX_AXIS_POINTS} points")
    steps = np.diff(axis)
    if np.any(steps <= 0) or not np.allclose(steps, steps[0], rtol=1e-10, atol=0):
        raise ValueError(label+" must be strictly increasing and uniform")
    return axis


def _fields(value, shape, label):
    if not isinstance(value, Mapping) or set(value) != set(NAMES):
        raise ValueError(label+" must contain exactly r,L,c")
    result = {key: _array(value[key], label+"."+key, shape) for key in NAMES}
    if np.any(result["r"] <= 0):
        raise ValueError(label+" radius must be positive")
    with np.errstate(over="ignore", under="ignore"):
        C = np.exp(result["L"])
    if np.any(C <= 0) or not np.all(np.isfinite(C)):
        raise ValueError(label+" lapse must be positive and representable")
    return result


def _alpha(value):
    result = _array(value, "alpha", ())
    if result < 0:
        raise ValueError("alpha must be nonnegative")
    return float(result)


def _budget_count(budget):
    if not callable(getattr(budget, "check", None)):
        raise ValueError("budget must expose check(count) and evaluations")
    count = getattr(budget, "evaluations", None)
    if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)) or count < 0:
        raise ValueError("budget evaluations must be a nonnegative integer")
    return int(count)


def _validated_tube(tube, boundaries):
    if not isinstance(tube, Mapping) or set(tube) != set(NAMES):
        raise ValueError("tube must contain exactly r,L,c")
    result = {}
    for name in NAMES:
        pair = _array(tube[name], "tube."+name, (2,))
        if pair[0] >= pair[1] or (name == "r" and pair[0] <= 0):
            raise ValueError("tube bounds must be increasing with positive radius floor")
        if name == "L":
            with np.errstate(over="ignore", under="ignore"):
                lapses = np.exp(pair)
            if np.any(lapses <= 0) or not np.all(np.isfinite(lapses)):
                raise ValueError("tube lapse bounds must be representable")
        if any(np.any(values[name] < pair[0]) or np.any(values[name] > pair[1]) for values in boundaries):
            raise ValueError("supplied boundary exits the declared tube")
        result[name] = tuple(map(float, pair))
    return result


def solve_rectangle(Ubar, x, source, bottom, *, alpha, tube, budget, cell_options=None):
    """March immutable characteristic values with one shared finite budget.

    source is X(0,x), bottom is X(Ubar,x_a), each containing arrays r,L,c.
    A corner difference <=5e-13 is admitted as roundoff and recorded; the
    source value is the stored common corner. Other differences are errors.
    Each cell uses SW=(i-1,j-1), NW=(i-1,j), SE=(i,j-1). No retry,
    subdivision, amplitude change or replacement boundary is performed.
    The caller owns the original physical domain and the declared tube.
    """
    Ubar, x = _axis(Ubar, "Ubar"), _axis(x, "x")
    alpha = _alpha(alpha)
    source = _fields(source, (len(x),), "source")
    bottom = _fields(bottom, (len(Ubar),), "bottom")
    tube = _validated_tube(tube, (source, bottom))
    if cell_options is not None and (not isinstance(cell_options, Mapping)
                                    or not set(cell_options) <= CELL_OPTIONS):
        raise ValueError("cell_options contains an unsupported option")
    options = {} if cell_options is None else copy.deepcopy(dict(cell_options))
    corner_difference = {name: float(source[name][0]-bottom[name][0]) for name in NAMES}
    if max(abs(value) for value in corner_difference.values()) > CORNER_TOLERANCE:
        raise ValueError("source and bottom corner values are incompatible")
    start_count = _budget_count(budget)
    states = {name: np.empty((len(Ubar), len(x))) for name in NAMES}
    for name in NAMES:
        states[name][:, 0] = bottom[name]
        states[name][0, :] = source[name]
    integers = ("updates", "backtrack_trials", "rhs_evaluations", "jacobian_evaluations")
    reals = ("residual_norm", "correction_norm", "condition_number")
    cells = {name: np.zeros((len(Ubar)-1, len(x)-1), dtype=int if name in integers else float)
             for name in integers+reals}
    cells["seed_kind"] = np.full((len(Ubar)-1, len(x)-1), "", dtype="U24")
    completed = 0
    for i in range(1, len(Ubar)):
        for j in range(1, len(x)):
            corners = {label: {name: float(states[name][ii, jj]) for name in NAMES}
                       for label, ii, jj in (("sw", i-1, j-1), ("nw", i-1, j), ("se", i, j-1))}
            try:
                row = solve_characteristic_cell(corners["sw"], corners["nw"], corners["se"],
                    du=float(Ubar[i]-Ubar[i-1]), dv=float(x[j]-x[j-1]), alpha=alpha,
                    region="exterior", tube=tube, budget=budget, **options)
            except (CellSolveError, RuntimeError, ValueError, FloatingPointError) as error:
                detail = error.diagnostics if isinstance(error, CellSolveError) else {"detail": str(error)}
                reason = error.reason if isinstance(error, CellSolveError) else (
                    "budget_exhausted" if isinstance(error, RuntimeError) and "budget" in str(error)
                    else "cell_exception")
                raise RectangleSolveError({"indices": [i, j],
                    "coordinates": {"Ubar": Ubar[i-1:i+1].tolist(), "x": x[j-1:j+1].tolist()},
                    "input_corners": corners, "reason": reason, "cell": detail,
                    "completed_cells": completed, "budget_start": start_count,
                    "budget_end": _budget_count(budget),
                    "completed_rhs_evaluations": int(np.sum(cells["rhs_evaluations"])),
                    "completed_jacobian_evaluations": int(np.sum(cells["jacobian_evaluations"]))}) from error
            for name in NAMES:
                states[name][i, j] = row["ne"][name]
            for name in cells:
                cells[name][i-1, j-1] = row[name]
            completed += 1
    budget.check()
    summary = {"cell_count": completed, "budget_start": start_count, "budget_end": _budget_count(budget),
               **{name: int(np.sum(cells[name])) for name in integers if name != "updates"},
               "max_updates": int(np.max(cells["updates"])),
               "seed_counts": {str(name): int(count) for name, count in zip(*np.unique(cells["seed_kind"], return_counts=True))},
               **{"max_"+name: float(np.max(cells[name])) for name in reals}}
    return {"Ubar": Ubar, "x": x, "states": states, "cells": cells, "summary": summary,
            "source": source, "bottom": bottom,
            "settings": {"alpha": alpha, "tube": tube, "cell_options": options,
                         "axis_order": ["Ubar", "x"], "corner_tolerance": CORNER_TOLERANCE,
                         "corner_choice": "source", "corner_difference": corner_difference}}


def _derivatives(values, spacing, axis):
    """Second-order centred interiors and independent one-sided endpoints.

    Low endpoint coefficients are (-3,4,-1)/(2h) for the first derivative
    and (2,-5,4,-1)/h² for the second. High endpoint coefficients reflect
    these with the first-derivative sign reversed.
    """
    f = np.moveaxis(values, axis, 0)
    first, second = np.empty_like(f), np.empty_like(f)
    first[1:-1] = (f[2:]-f[:-2])/(2*spacing)
    second[1:-1] = (f[2:]-2*f[1:-1]+f[:-2])/spacing**2
    first[0], first[-1] = (-3*f[0]+4*f[1]-f[2])/(2*spacing), (3*f[-1]-4*f[-2]+f[-3])/(2*spacing)
    second[0] = (2*f[0]-5*f[1]+4*f[2]-f[3])/spacing**2
    second[-1] = (2*f[-1]-5*f[-2]+4*f[-3]-f[-4])/spacing**2
    return np.moveaxis(first, 0, axis), np.moveaxis(second, 0, axis)


def _stats(values, mask):
    selected = values[mask]
    if not selected.size:
        return {"count": 0, "max_abs": None, "rms": None}
    return {"count": int(selected.size), "max_abs": float(np.max(abs(selected))),
            "rms": float(np.sqrt(np.mean(selected*selected)))}


def _summaries(fields, masks):
    return {region: {name: _stats(values, mask) for name, values in fields.items()}
            for region, mask in masks.items()}


def grid_diagnostics(Ubar, x, states, *, alpha, source_fronts=(), source_jets=None, bottom_jets=None):
    """Pure nodal diagnostics independent of every solved mixed-cell RHS.

    No source/ODE/PDE is called. Pure second derivatives are finite
    differences, never reconstructed from Einstein equations. Pulse fronts
    are original support_x labels; their crossing stencils stay in the
    all-interior summary and also get a separate front-near summary.
    A C3, piecewise W4,infinity pulse is not assumed globally C4. Observed
    derivative refinement must be assessed independently of state refinement.
    Optional boundary jets use the six JET_NAMES and may include other keys.
    """
    Ubar, x = _axis(Ubar, "Ubar", 4), _axis(x, "x", 4)
    states = _fields(states, (len(Ubar), len(x)), "states")
    alpha = _alpha(alpha)
    if not isinstance(source_fronts, (tuple, list, np.ndarray)):
        raise ValueError("source_fronts must be empty or two original support labels")
    fronts = np.asarray(source_fronts)
    if fronts.size:
        fronts = _array(source_fronts, "source_fronts", (2,))
        if not x[0] < fronts[0] < fronts[1] < x[-1]:
            raise ValueError("source fronts must lie strictly inside x")
    elif fronts.shape != (0,):
        raise ValueError("source_fronts must be empty or two original support labels")
    derivatives = {}
    for name in NAMES:
        for axis, axis_name, grid in ((0, "Ubar", Ubar), (1, "x", x)):
            first, second = _derivatives(states[name], grid[1]-grid[0], axis)
            prefix = "log_C" if name == "L" else name
            derivatives[prefix+"_"+axis_name] = first
            if name == "r":
                derivatives["r_"+axis_name+"_"+axis_name] = second
    r, L, c = (states[name] for name in NAMES)
    u = c/r
    D = 1+alpha*u*u/6
    fields = {}
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        for axis, name in (("Ubar", "E_U"), ("x", "E_x")):
            ra = derivatives["r_"+axis]
            ua = (derivatives["c_"+axis]-u*ra)/r
            fields[name] = derivatives["r_"+axis+"_"+axis]-derivatives["log_C_"+axis]*ra+alpha*r*ua*ua/(2*D)
        fields["m"] = r*(1+4*derivatives["r_Ubar"]*derivatives["r_x"]*np.exp(-L))/2
    if not all(np.all(np.isfinite(v)) for v in (*fields.values(), *derivatives.values())):
        raise ValueError("nonfinite finite-difference diagnostics")
    interior = np.zeros(r.shape, dtype=bool)
    interior[1:-1, 1:-1] = True
    front_columns = np.zeros(len(x), dtype=bool)
    for front in fronts:
        front_columns[1:-1] |= (x[:-2] <= front) & (front <= x[2:])
    near = interior & front_columns[None, :]
    masks = {"interior": interior, "front_near": near, "away": interior & ~near}
    boundary_defects, boundary_summary = {}, {}
    pulse = np.zeros(len(x), dtype=bool)
    postpulse = np.zeros(len(x), dtype=bool)
    if len(fronts):
        pulse = (x > fronts[0]) & (x < fronts[1])
        postpulse = x >= fronts[1]
    source_masks = {"all": np.ones(len(x), dtype=bool), "pulse": pulse, "postpulse": postpulse}
    for label, supplied, size, selection in (("source", source_jets, len(x), (0, slice(None))),
                                            ("bottom", bottom_jets, len(Ubar), (slice(None), 0))):
        if supplied is None:
            continue
        if not isinstance(supplied, Mapping) or not set(JET_NAMES) <= set(supplied):
            raise ValueError(label+" jets must include all six first derivatives")
        errors = {key: derivatives[key][selection]-_array(supplied[key], label+"."+key, (size,))
                  for key in JET_NAMES}
        boundary_defects[label] = errors
        selections = source_masks if label == "source" else {"all": np.ones(size, dtype=bool)}
        boundary_summary[label] = _summaries(errors, selections)
    return {"Ubar": Ubar, "x": x, "alpha": alpha, "source_fronts": fronts.copy(),
            **fields, "derivatives": derivatives, "masks": masks,
            "summary": _summaries(fields, masks), "source_masks": source_masks,
            "boundary_defects": boundary_defects, "boundary_summary": boundary_summary}


def compare_common_grids(coarse, fine, *, coarse_zero=None, fine_zero=None):
    """Common interior diagnostics on nested two-to-one grids.

    Returned differences are fine-minus-coarse. Optional zero-amplitude
    diagnostics additionally compare (finite-zero) between resolutions;
    these paired quantities never replace either raw constraint residual.
    The front mask is the union of both actual stencil-crossing masks.
    """
    for axis in ("Ubar", "x"):
        c, f = _axis(coarse[axis], axis, 4), _axis(fine[axis], axis, 4)
        if len(f) != 2*len(c)-1 or not np.allclose(c, f[::2], rtol=0, atol=2e-14):
            raise ValueError("comparison requires nested two-to-one common coordinates")
    if coarse["alpha"] != fine["alpha"] or not np.array_equal(coarse["source_fronts"], fine["source_fronts"]):
        raise ValueError("comparison must preserve amplitude and original source fronts")
    shape = (len(coarse["Ubar"]), len(coarse["x"]))
    selection = (slice(1, -1), slice(1, -1))
    downsample = (slice(2, -2, 2), slice(2, -2, 2))
    raw_coarse = {key: _array(coarse[key], "coarse."+key, shape)[selection] for key in ("E_U", "E_x", "m")}
    fine_shape = (len(fine["Ubar"]), len(fine["x"]))
    raw_fine = {key: _array(fine[key], "fine."+key, fine_shape)[downsample] for key in raw_coarse}
    near = np.asarray(coarse["masks"]["front_near"], dtype=bool)[selection]
    near = near | np.asarray(fine["masks"]["front_near"], dtype=bool)[downsample]
    masks = {"interior": np.ones(near.shape, dtype=bool), "front_near": near, "away": ~near}
    differences = {key: raw_fine[key]-raw_coarse[key] for key in raw_coarse}
    result = {"Ubar": np.asarray(coarse["Ubar"])[1:-1].copy(),
              "x": np.asarray(coarse["x"])[1:-1].copy(), "masks": masks,
              "coarse": raw_coarse, "fine": raw_fine, "differences": differences,
              "coarse_summary": _summaries(raw_coarse, masks),
              "fine_summary": _summaries(raw_fine, masks), "summary": _summaries(differences, masks)}
    if (coarse_zero is None) != (fine_zero is None):
        raise ValueError("paired comparison needs both zero-amplitude baselines")
    if coarse_zero is not None:
        for baseline, actual in ((coarse_zero, coarse), (fine_zero, fine)):
            if baseline["alpha"] != 0 or not all(np.array_equal(baseline[axis], actual[axis]) for axis in ("Ubar", "x", "source_fronts")):
                raise ValueError("zero baselines must share the exact grids and source fronts")
        zero_comparison = compare_common_grids(coarse_zero, fine_zero)
        paired_coarse = {key: raw_coarse[key]-zero_comparison["coarse"][key] for key in raw_coarse}
        paired_fine = {key: raw_fine[key]-zero_comparison["fine"][key] for key in raw_fine}
        paired_difference = {key: paired_fine[key]-paired_coarse[key] for key in raw_coarse}
        result["paired"] = {"coarse": paired_coarse, "fine": paired_fine,
                            "differences": paired_difference, "summary": _summaries(paired_difference, masks)}
    return result
