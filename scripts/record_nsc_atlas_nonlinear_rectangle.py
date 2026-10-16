#!/usr/bin/env python3
"""Create or authenticate the fixed four-grid nonlinear exterior experiment.

Checking reconstructs finite-difference diagnostics from saved nodal arrays;
it never integrates a source or solves a characteristic cell. The experiment
declaration below is frozen with the producer, rather than binding live prose.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import time
import uuid

import numpy as np
from validation_source_successors import verify_saved_bindings, verify_saved_source
import scipy
from scipy.special import lambertw

from recursive_horizons.nsc_atlas_characteristic_initialization import (
    Budget, source_geometry, source_jet, validate_geometry,
)
from recursive_horizons.nsc_atlas_validation import source_jet_roundoff_tolerances
from recursive_horizons.nsc_atlas_nonlinear_source import ExactSource, source_amplitude_bound


ROOT = Path(__file__).resolve().parents[1]
RELATIVE_OUTPUT = "results/development/nsc-atlas-nonlinear-rectangle-v1"
SCHEMA = "nsc-atlas-nonlinear-rectangle-v1"
STATUS = "FINITE_AMPLITUDE_EXTERIOR_RECTANGLE"
SCOPE = ("original exterior rectangle only; sampled grid sensitivity and discrete null constraints; "
         "no wall-born cells, cap evolution, central readout or uniform error enclosure")
TUBE = {"r": [.2, .8], "L": [-2., 2.], "c": [-4., 4.]}
CELL_OPTIONS = {"atol": [1e-12]*3, "rtol": 1e-12, "reference_scales": [1.]*3,
                "max_updates": 12, "max_backtracks": 8, "maximum_condition": 1e12}
SETTINGS = {"geometry_order": 64, "cells": [16, 32], "source_nodes": 33,
            "source_rtol": 2e-10, "source_atol": 1e-12, "source_max_step_divisor": 32,
            "amplitudes": ["zero", "half_source_threshold"], "region": "exterior",
            "global_cpu_seconds": 120., "global_evaluations": 60000,
            "source_cpu_seconds": 30., "source_evaluations": 20000,
            "cell_options": CELL_OPTIONS, "tube": TUBE}
CASES = (("zero", 16), ("zero", 32), ("finite", 16), ("finite", 32))
SOURCE_KEYS = {"x", "rJ", "rJ_x", "h", "h_x", "r", "C", "log_C", "c",
               "r_Ubar", "r_x", "c_Ubar", "c_x", "log_C_Ubar", "log_C_x", "m", "r_xx",
               "source_constraint", "source_constraint_jordan"}
BOTTOM_KEYS = (SOURCE_KEYS-{"r_xx", "source_constraint", "source_constraint_jordan"}) | {"Ubar", "U"}
SOURCES = (
    "scripts/record_nsc_atlas_nonlinear_rectangle.py",
    "tests/test_record_nsc_atlas_nonlinear_rectangle.py",
    "src/recursive_horizons/nsc_atlas_nonlinear_rectangle.py",
    "tests/test_nsc_atlas_nonlinear_rectangle.py",
    "src/recursive_horizons/nsc_atlas_nonlinear_source.py",
    "tests/test_nsc_atlas_nonlinear_source.py",
    "src/recursive_horizons/nsc_atlas_characteristic_cell.py",
    "tests/test_nsc_atlas_characteristic_cell.py",
    "src/recursive_horizons/nsc_atlas_nonlinear_bulk.py",
    "tests/test_nsc_atlas_nonlinear_bulk.py",
    "src/recursive_horizons/nsc_atlas_characteristic_initialization.py",
    "src/recursive_horizons/nsc_atlas_field_preparation.py",
    "src/recursive_horizons/nsc_atlas_field_transfer.py",
    "src/recursive_horizons/nsc_atlas_causal_reference.py",
    "src/recursive_horizons/nsc_atlas_causal_observers.py",
    "src/recursive_horizons/nsc_atlas_connection_slice.py",
    "src/recursive_horizons/nsc_atlas_moving_boundaries.py",
    "src/recursive_horizons/__init__.py", "pyproject.toml", "scripts/run_python.py",
)


def plain(value):
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def encoded(value):
    return (json.dumps(plain(value), indent=2, sort_keys=True, allow_nan=False)+"\n").encode()


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def runtime_versions():
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


class SourceBudget(Budget):
    """A source's own bounded lifetime, with every charge also on its parent."""
    def __init__(self, parent):
        super().__init__(maximum_evaluations=20000, cpu_seconds=30.)
        self.parent = parent

    def check(self, count=0):
        if type(count) is not int or count < 0:
            raise ValueError("source budget charges must be nonnegative integers")
        self.evaluations += count
        self.parent.check(count)
        if self.evaluations > self.maximum_evaluations or time.process_time()-self.started > self.cpu_seconds:
            raise RuntimeError("rectangle source child budget exhausted")


def _boundary_states(rows, step=1):
    return {key: np.asarray(rows[source_key])[::step]
            for key, source_key in (("r", "r"), ("L", "log_C"), ("c", "c"))}


def _boundary_jets(rows, step=1):
    return {key: np.asarray(value)[::step] for key, value in rows.items()}


def schwarzschild(Ubar, x):
    """Independent exact vacuum geometry in the experiment's affine gauge."""
    Ubar, x = np.asarray(Ubar)[:, None], np.asarray(x)[None, :]
    q = 1-Ubar/3.2
    product = q*(1-x/.8)*np.exp(x/.8)
    radius = .8*(1+lambertw(-product/math.e).real)
    lapse = x/radius*np.exp((x-radius)/.8)
    return {"r": radius, "L": np.log(x/radius)+(x-radius)/.8, "C": lapse}


def statistics(values):
    values = np.asarray(values, dtype=float)
    return {"maximum_absolute": float(np.max(np.abs(values))),
            "rms": float(np.sqrt(np.mean(values*values)))}


def _state_difference(left, right):
    arrays = {key: np.asarray(left[key])-np.asarray(right[key]) for key in ("r", "L", "c")}
    return {"arrays": arrays, "all_nodes": {key: statistics(value) for key, value in arrays.items()},
            "interior": {key: statistics(value[1:-1, 1:-1]) for key, value in arrays.items()}}


def diagnostics(record):
    """Pure saved-array assessment; no source integration or cell solving."""
    from recursive_horizons.nsc_atlas_nonlinear_rectangle import grid_diagnostics, compare_common_grids
    rows, vacuum, quiet = [], [], []
    for grid in record["grids"]:
        label, cells = grid["amplitude"], grid["cell_count_per_axis"]
        data = record["initial_data"][label]
        step = 32//cells
        diagnostic = grid_diagnostics(np.asarray(grid["Ubar"]), np.asarray(grid["x"]),
            {key: np.asarray(value) for key, value in grid["states"].items()}, alpha=data["alpha"],
            source_fronts=record["geometry"]["support_x"],
            source_jets=_boundary_jets(data["source"], step),
            bottom_jets=_boundary_jets(data["bottom"], step))
        rows.append(diagnostic)
        quiet.append({"amplitude": label, "cells": cells,
            "maximum_absolute_c_before_source_front": float(np.max(np.abs(
                np.asarray(grid["states"]["c"])[:, np.asarray(grid["x"]) <= record["geometry"]["support_x"][0]])))})
        if label == "zero":
            exact = schwarzschild(grid["Ubar"], grid["x"])
            vacuum.append({"cells": cells,
                "metric_error": {key: statistics(np.asarray(grid["states"][key])-exact[key])
                                 for key in ("r", "L")},
                "lapse_error": statistics(np.exp(grid["states"]["L"])-exact["C"]),
                "finite_difference_mass_error": statistics(np.asarray(diagnostic["m"])-.4)})
    refinements = {}
    for label, coarse_index in (("zero", 0), ("finite", 2)):
        coarse, fine = record["grids"][coarse_index:coarse_index+2]
        common = {key: np.asarray(value)[::2, ::2] for key, value in fine["states"].items()}
        refinements[label] = _state_difference(common, coarse["states"])
    paired = [_state_difference(record["grids"][j+2]["states"], record["grids"][j]["states"])
              for j in (0, 1)]
    paired_refinement = _state_difference(
        {key: np.asarray(value)[::2, ::2] for key, value in paired[1]["arrays"].items()}, paired[0]["arrays"])
    return plain({"grids": rows, "vacuum_controls": vacuum, "quiet_field_controls": quiet,
        "state_refinement": refinements,
        "finite_minus_zero": paired, "finite_minus_zero_refinement": paired_refinement,
        "constraint_refinement_zero": compare_common_grids(rows[0], rows[1]),
        "constraint_refinement_finite": compare_common_grids(rows[2], rows[3],
            coarse_zero=rows[0], fine_zero=rows[1]),
        "difference_convention": "fine common nodes minus coarse; paired states are finite minus zero"})


def calculate():
    from recursive_horizons.nsc_atlas_nonlinear_rectangle import solve_rectangle
    budget = Budget(maximum_evaluations=60000, cpu_seconds=120.)
    phase = {"name": "geometry"}
    try:
        geometry = source_geometry(64, budget)
        geometry_evaluations = budget.evaluations
        CJ0 = -3.2/geometry["U0"]
        Ubar = np.linspace(0., CJ0*(geometry["Ua"]-geometry["U0"]), 33)
        x = np.linspace(.4, geometry["x_last"], 33)
        initial = {}
        bounds = None
        for label in ("zero", "finite"):
            phase = {"name": "source", "amplitude": label}
            child = SourceBudget(budget)
            if bounds is None:
                bounds = source_amplitude_bound(geometry, budget=child)
            alpha = 0. if label == "zero" else bounds["alpha_star"]/2
            source = ExactSource(geometry, alpha=alpha, budget=child,
                                 rtol=2e-10, atol=1e-12)
            source_rows, bottom_rows = source.sample(x), source.quiet_bottom(Ubar)
            child.check()
            initial[label] = {"alpha": alpha, "domain": source.domain,
                "source": plain(source_rows), "bottom": plain(bottom_rows),
                "max_step": source.max_step, "nfev": source.nfev,
                "execution": {"evaluations": child.evaluations,
                              "cpu_seconds": time.process_time()-child.started}}
        record = {"schema": SCHEMA, "status": STATUS, "scope": SCOPE,
            "runtime": runtime_versions(), "settings": SETTINGS, "geometry": geometry,
            "source_bounds": bounds, "initial_data": initial, "grids": [],
            "refinement_is_rigorous_error_enclosure": False}
        for label, cells in CASES:
            phase = {"name": "march", "amplitude": label, "cells": cells}
            step, data = 32//cells, initial[label]
            grid = solve_rectangle(Ubar[::step], x[::step],
                _boundary_states(data["source"], step), _boundary_states(data["bottom"], step),
                alpha=data["alpha"], tube=TUBE, budget=budget, cell_options=CELL_OPTIONS)
            record["grids"].append(plain(dict(grid, amplitude=label, cell_count_per_axis=cells)))
            budget.check()
        phase = {"name": "diagnostics"}
        record["diagnostics"] = diagnostics(record)
        budget.check()
        record["execution"] = {"cpu_seconds": time.process_time()-budget.started,
            "evaluations": budget.evaluations, "geometry_evaluations": geometry_evaluations,
            "source_integrations": 2, "source_segments": 6, "rectangle_marches": 4,
            "successful_cells": 2*(16*16+32*32), "leading_field_solves": 0,
            "wall_solves": 0, "cap_solves": 0}
        return plain(record)
    except Exception as error:
        error.rectangle_failure = {"phase": phase, "exception": type(error).__name__,
            "message": str(error), "diagnostics": getattr(error, "diagnostics", {}),
            "execution": {"cpu_seconds": time.process_time()-budget.started,
                          "evaluations": budget.evaluations}}
        raise


def _array(value, shape, label):
    if any(isinstance(v, (bool, np.bool_)) for v in np.asarray(value, dtype=object).flat):
        raise ValueError(label+" contains a boolean")
    array = np.asarray(value)
    if array.shape != shape or array.dtype.kind not in "iuf" or not np.all(np.isfinite(array)):
        raise ValueError(label+" has invalid finite numeric shape")
    return array.astype(float)


def _count(value, maximum, label):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError(label+" is not a bounded integer")


def _real(value, label, low=0., high=math.inf):
    if type(value) not in (int, float) or not math.isfinite(value) or not low <= value <= high:
        raise ValueError(label+" is outside its finite domain")


def _close(actual, expected, label, atol=2e-13):
    if not np.allclose(actual, expected, rtol=2e-13, atol=atol):
        raise ValueError(label+" identity mismatch")


def _same_recomputed(actual, expected):
    """Compare recomputed arithmetic, retaining exact structure and counts."""
    if isinstance(expected, dict):
        return isinstance(actual, dict) and set(actual) == set(expected) and all(
            _same_recomputed(actual[key], value) for key, value in expected.items())
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(
            _same_recomputed(a, b) for a, b in zip(actual, expected))
    if type(expected) is float:
        return (type(actual) in (int, float) and math.isfinite(actual) and math.isfinite(expected)
                and math.isclose(actual, expected, rel_tol=2e-12, abs_tol=2e-13))
    return type(actual) is type(expected) and actual == expected


def _validate_initial(record):
    geometry = record["geometry"]
    if (set(geometry) != {"order", "U0", "Ua", "Ub", "Va", "Vb", "x_corner", "x_last", "support_x", "support_V"}
            or type(geometry["order"]) is not int or geometry["order"] != 64):
        raise ValueError("original geometry inventory/order mismatch")
    validate_geometry(geometry)
    bounds = source_amplitude_bound(geometry)
    if (not _same_recomputed(record["source_bounds"], bounds)
            or record["source_bounds"]["test_alpha"] != record["source_bounds"]["alpha_star"]/2
            or set(record["initial_data"]) != {"zero", "finite"}):
        raise ValueError("source threshold or amplitude inventory mismatch")
    CJ0 = -3.2/geometry["U0"]
    Ubar = np.linspace(0., CJ0*(geometry["Ua"]-geometry["U0"]), 33)
    x = np.linspace(.4, geometry["x_last"], 33)
    expected_domain = {key: geometry[key] for key in ("U0", "Ua", "Ub", "Va", "Vb")}
    expected_domain.update(x_a=.4, x_b=geometry["x_last"], Ubar_a=0., Ubar_b=Ubar[-1], CJ0=CJ0)
    p, px, _ = source_jet(x, geometry)
    g_roundoff, gx_roundoff = source_jet_roundoff_tolerances(x, geometry)
    quiet = (x <= geometry["support_x"][0]) | (x >= geometry["support_x"][1])
    g, gx = np.where(quiet, 0., x*p), np.where(quiet, 0., p+x*px)
    for label in ("zero", "finite"):
        data = record["initial_data"][label]
        # This is the retained model input, already authenticated against the
        # recomputed bound above; do not replace it with platform-roundoff data.
        alpha = 0. if label == "zero" else record["source_bounds"]["test_alpha"]
        if (type(data["alpha"]) not in (int, float) or data["alpha"] != alpha
                or data["domain"] != expected_domain or data["max_step"] != bounds["length"]/32):
            raise ValueError("source amplitude, domain or step mismatch")
        _count(data["nfev"], 20000, "source RHS calls")
        _count(data["execution"]["evaluations"], 20000, "source evaluations")
        _real(data["execution"]["cpu_seconds"], "source CPU", high=30)
        # Every RHS/profile call, both 33-point queries and initialization checks are charged.
        expected_charges = data["nfev"]+66+3+(3 if label == "zero" else 0)
        if data["execution"]["evaluations"] != expected_charges:
            raise ValueError("source evaluation accounting mismatch")
        for kind, keys, coordinates in (("source", SOURCE_KEYS, x), ("bottom", BOTTOM_KEYS, Ubar)):
            row = data[kind]
            if set(row) != keys:
                raise ValueError(kind+" first-jet inventory mismatch")
            arrays = {key: _array(value, (33,), kind+"."+key) for key, value in row.items()}
            if np.any(arrays["r"] <= 0) or np.any(arrays["C"] <= 0):
                raise ValueError("initial radius/lapse is not positive")
            _close(np.exp(arrays["log_C"]), arrays["C"], "initial lapse")
            _close(arrays["m"], arrays["r"]*(1+4*arrays["r_Ubar"]*arrays["r_x"]/arrays["C"])/2,
                   "initial mass")
            if kind == "source":
                _close(arrays["x"], coordinates, "source coordinates", atol=0)
                _close(arrays["c"], g, "original compact source", atol=2e-13+g_roundoff)
                _close(arrays["c_x"], gx, "original compact source derivative", atol=2e-13+gx_roundoff)
                _close(arrays["rJ"], x+alpha*arrays["h"], "source scaled radius")
                _close(arrays["rJ_x"], 1+alpha*arrays["h_x"], "source scaled slope")
                _close(arrays["C"], 1-alpha*g*g/(6*arrays["rJ"]**2), "source conformal factor")
                _close(arrays["r"], arrays["rJ"]*np.sqrt(arrays["C"]), "source Einstein radius")
                _close(arrays["r_x"], (arrays["rJ"]*arrays["rJ_x"]-alpha*g*gx/6)/arrays["r"],
                       "source Einstein radial jet")
                pJ, pJx = g/arrays["rJ"], (gx-g*arrays["rJ_x"]/arrays["rJ"])/arrays["rJ"]
                _close(arrays["log_C_x"], -alpha*pJ*pJx/(3*arrays["C"]), "source lapse jet")
            else:
                _close(arrays["Ubar"], coordinates, "bottom coordinates", atol=0)
                _close(arrays["U"], geometry["U0"]+Ubar/CJ0, "original bottom U")
                _close(arrays["x"], np.full(33, .4), "bottom x", atol=0)
                vacuum = schwarzschild(Ubar, np.array([.4]))
                for key, actual_key in (("r", "r"), ("L", "log_C"), ("C", "C")):
                    _close(arrays[actual_key], vacuum[key][:, 0], "quiet bottom "+key)
                for key in ("c", "c_Ubar", "c_x", "h", "h_x"):
                    if np.any(arrays[key] != 0):
                        raise ValueError("quiet bottom field mismatch")
    return Ubar, x


def validate_record(record):
    try:
        if (record["schema"] != SCHEMA or record["status"] != STATUS or record["scope"] != SCOPE
                or encoded(record["settings"]) != encoded(SETTINGS)
                or record["refinement_is_rigorous_error_enclosure"] is not False):
            raise ValueError("rectangle declaration mismatch")
        runtime = record["runtime"]
        if (set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(value, str) or not value or not value[0].isdigit() or "." not in value
                       for value in runtime.values())):
            raise ValueError("runtime inventory mismatch")
        Ubar, x = _validate_initial(record)
        if len(record["grids"]) != 4:
            raise ValueError("four fixed grids required")
        cell_charges = 0
        previous_end = 960+sum(item["execution"]["evaluations"] for item in record["initial_data"].values())
        for grid, (label, cells) in zip(record["grids"], CASES):
            if (grid["amplitude"] != label or type(grid["cell_count_per_axis"]) is not int
                    or grid["cell_count_per_axis"] != cells):
                raise ValueError("grid identity mismatch")
            step, data = 32//cells, record["initial_data"][label]
            source, bottom = _boundary_states(data["source"], step), _boundary_states(data["bottom"], step)
            corner = {key: float(source[key][0]-bottom[key][0]) for key in source}
            expected_settings = {"alpha": data["alpha"], "tube": TUBE, "cell_options": CELL_OPTIONS,
                "axis_order": ["Ubar", "x"], "corner_tolerance": 5e-13,
                "corner_choice": "source", "corner_difference": corner}
            if (encoded(grid["settings"]) != encoded(expected_settings)
                    or encoded(grid["source"]) != encoded(source) or encoded(grid["bottom"]) != encoded(bottom)
                    or max(abs(value) for value in corner.values()) > 5e-13):
                raise ValueError("march settings or supplied boundaries mismatch")
            if not np.array_equal(grid["Ubar"], Ubar[::step]) or not np.array_equal(grid["x"], x[::step]):
                raise ValueError("grid coordinates mismatch")
            if set(grid["states"]) != {"r", "L", "c"}:
                raise ValueError("grid state inventory mismatch")
            for key, values in grid["states"].items():
                array = _array(values, (cells+1, cells+1), "state "+key)
                if np.any(array < TUBE[key][0]) or np.any(array > TUBE[key][1]):
                    raise ValueError("state exits fixed engineering tube")
                if (not np.array_equal(array[0], source[key])
                        or not np.array_equal(array[1:, 0], bottom[key][1:])):
                    raise ValueError("march boundary differs from saved original input")
            if set(grid["cells"]) != {"updates", "backtrack_trials", "rhs_evaluations", "jacobian_evaluations",
                                      "residual_norm", "correction_norm", "condition_number", "seed_kind"}:
                raise ValueError("cell diagnostic inventory mismatch")
            seeds = np.asarray(grid["cells"]["seed_kind"])
            if seeds.shape != (cells, cells) or not np.all(np.isin(seeds, ("affine_predictor", "convex_corner_seed"))):
                raise ValueError("cell seed inventory mismatch")
            affine_inside = np.ones((cells, cells), dtype=bool)
            for key, values in grid["states"].items():
                values = np.asarray(values)
                affine = values[:-1, 1:]+values[1:, :-1]-values[:-1, :-1]
                affine_inside &= (affine >= TUBE[key][0]) & (affine <= TUBE[key][1])
            if not np.array_equal(seeds, np.where(affine_inside, "affine_predictor", "convex_corner_seed")):
                raise ValueError("cell seed differs from fixed predictor/tube rule")
            for key in ("rhs_evaluations", "jacobian_evaluations"):
                values = _array(grid["cells"][key], (cells, cells), "cell "+key)
                if np.asarray(grid["cells"][key]).dtype.kind not in "iu" or np.any(values < 1) or np.any(values > 200):
                    raise ValueError("cell charges outside bounded integer range")
                cell_charges += int(values.sum())
            for key, upper in (("updates", 12), ("backtrack_trials", 96), ("residual_norm", 1.), ("correction_norm", 1.),
                               ("condition_number", 1e12)):
                values = _array(grid["cells"][key], (cells, cells), "cell "+key)
                if np.any(values < 0) or np.any(values > upper):
                    raise ValueError("cell convergence/iteration declaration invalid")
                if key in ("updates", "backtrack_trials") and np.asarray(grid["cells"][key]).dtype.kind not in "iu":
                    raise ValueError("cell iteration counts must be integers")
            updates, trials, rhs, jac = (np.asarray(grid["cells"][key]) for key in
                                        ("updates", "backtrack_trials", "rhs_evaluations", "jacobian_evaluations"))
            if (np.any(jac != updates+1) or np.any(trials < updates) or np.any(trials > 8*updates)
                    or np.any(rhs < jac+updates) or np.any(rhs > jac+trials)
                    or np.any(np.asarray(grid["cells"]["condition_number"]) < 1)):
                raise ValueError("cell Newton/RHS/Jacobian accounting mismatch")
            costs = sum(int(np.sum(grid["cells"][key])) for key in ("rhs_evaluations", "jacobian_evaluations"))
            expected_summary = {"cell_count": cells*cells, "budget_start": previous_end, "budget_end": previous_end+costs,
                "seed_counts": {str(name): int(count) for name, count in zip(*np.unique(seeds, return_counts=True))},
                **{key: int(np.sum(grid["cells"][key])) for key in
                   ("backtrack_trials", "rhs_evaluations", "jacobian_evaluations")},
                "max_updates": int(np.max(grid["cells"]["updates"])),
                **{"max_"+key: float(np.max(grid["cells"][key])) for key in
                   ("residual_norm", "correction_norm", "condition_number")}}
            if not _same_recomputed(grid["summary"], expected_summary):
                raise ValueError("grid summary or cumulative transition mismatch")
            previous_end += costs
        if not _same_recomputed(record["diagnostics"], plain(diagnostics(record))):
            raise ValueError("saved-array diagnostics disagree with nodal states")
        execution = record["execution"]
        _real(execution["cpu_seconds"], "global CPU", high=120)
        _count(execution["evaluations"], 60000, "global evaluations")
        expected = {"geometry_evaluations": 960, "source_integrations": 2, "source_segments": 6,
                    "rectangle_marches": 4, "successful_cells": 2560, "leading_field_solves": 0,
                    "wall_solves": 0, "cap_solves": 0}
        if any(type(execution[key]) is not int or execution[key] != value for key, value in expected.items()):
            raise ValueError("experiment solve/count declaration mismatch")
        if execution["evaluations"] != 960+cell_charges+sum(
                item["execution"]["evaluations"] for item in record["initial_data"].values()):
            raise ValueError("cumulative budget accounting mismatch")
    except (KeyError, TypeError, IndexError, OverflowError) as error:
        raise ValueError("incomplete rectangle record") from error
    return True


def _failure_safe(value):
    value = plain(value)
    if isinstance(value, dict):
        return {key: _failure_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_failure_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return {"nonfinite": repr(value)}
    return value


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    bindings, runtime = source_bindings(), runtime_versions()
    started = time.process_time()
    try:
        record = calculate()
        if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
            raise RuntimeError("rectangle inputs or runtime changed during calculation")
        validate_record(record)
        elapsed = time.process_time()-started
        if elapsed > 120:
            raise RuntimeError("rectangle creation and diagnostic CPU budget exhausted")
        record["execution"]["cpu_seconds"] = elapsed
    except Exception as error:
        try:
            after = source_bindings()
        except OSError as binding_error:
            after = {"unavailable": str(binding_error)}
        receipt = {"schema": SCHEMA+"-failed-attempt", "status": "FAILED_NOT_COMPLETED",
                   "settings": SETTINGS, "source_bindings_before": bindings,
                   "source_bindings_after": after, "runtime": runtime,
                   "failure": getattr(error, "rectangle_failure", {"exception": type(error).__name__,
                                                                     "message": str(error)})}
        path = ROOT/".build/nonlinear-rectangle-failures"/(uuid.uuid4().hex+".json")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(encoded(_failure_safe(receipt)))
        error.failure_receipt = str(path)
        raise
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings,
                "outputs": {"rectangle.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"rectangle.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"],
            "vacuum_controls": record["diagnostics"]["vacuum_controls"],
            "state_refinement": {key: value["interior"] for key, value in
                                 record["diagnostics"]["state_refinement"].items()}}


def check(directory):
    directory = Path(directory)
    if ({path.name for path in directory.iterdir()} != {"rectangle.json", "manifest.json"}
            or any(path.is_symlink() for path in directory.iterdir())):
        raise ValueError("rectangle output inventory mismatch")
    manifest = json.loads((directory/"manifest.json").read_bytes())
    if manifest.get("schema") != SCHEMA+"-manifest" or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES):
        raise ValueError("rectangle source binding or manifest schema mismatch")
    raw = (directory/"rectangle.json").read_bytes()
    if manifest.get("outputs") != {"rectangle.json": {"sha256": digest(raw), "bytes": len(raw)}}:
        raise ValueError("rectangle numerical payload hash/size mismatch")
    validate_record(json.loads(raw))
    return {"checked": str(directory), "source_integrations": 0, "rectangle_marches": 0, "cell_solves": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    try:
        result = check(args.output) if args.check else create(args.output)
    except Exception as error:
        if hasattr(error, "failure_receipt"):
            print(json.dumps({"failed": str(error), "diagnostic": error.failure_receipt}))
        raise
    print(json.dumps(result, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
