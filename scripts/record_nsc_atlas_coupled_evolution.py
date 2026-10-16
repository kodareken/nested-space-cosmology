#!/usr/bin/env python3
"""Create/check the first fixed five-case cap-wall-exterior calculation.

Authentication uses retained arrays and reconstructs the final radial cap
constraints/prefix stage. It never integrates a source, reference gauge or
time evolution. The original exterior rectangle is authenticated and reused
as an independent saved-data comparator; it is never modified or evolved.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import time
import uuid

import numpy as np
from validation_source_successors import verify_saved_bindings, verify_saved_source
from scipy.special import lambertw

import record_nsc_atlas_nonlinear_rectangle as previous
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget, source_geometry, validate_geometry
from recursive_horizons.nsc_atlas_nonlinear_source import ExactSource, source_amplitude_bound
from recursive_horizons.nsc_atlas_characteristic_trace import SourceSnapshot
from recursive_horizons.nsc_atlas_validation import source_snapshot_roundoff_tolerances
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_coupled_evolution import quiet_initial_state, coupled_stage, evolve_coupled
from recursive_horizons.nsc_atlas_coupled_diagnostics import (
    coupled_history_diagnostics, compare_rectangle_overlap,
)

ROOT = Path(__file__).resolve().parents[1]
RELATIVE_OUTPUT = "results/development/nsc-atlas-coupled-evolution-v1"
SCHEMA = "nsc-atlas-coupled-evolution-v1"
STATUS = "FINITE_AMPLITUDE_COUPLED_WINDOW"
SCOPE = ("original incoming-label window with evolved exterior, moving wall and polar cap; "
         "sampled discretization/constraint diagnostics, not a uniform error enclosure, "
         "stability result or post-wall central-arrival evolution")
TUBE = previous.TUBE
CASES = (("zero", 16, 256), ("zero", 32, 256), ("finite", 16, 256),
         ("finite", 32, 256), ("finite", 32, 512))
SETTINGS = dict(geometry_order=64, gauge_order=64, gauge_nodes=513,
    source_exact_nodes=257, source_snapshot_knots=129, source_cache_nodes=513,
    source_rtol=2e-10, source_atol=1e-12, source_max_step_divisor=32,
    saved_time_nodes=129, exterior_quadrature_order=2, cap_quadrature_order=4,
    maximum_courant=.4, tube=TUBE, global_cpu_seconds=180., global_evaluations=8000000,
    source_cpu_seconds=30., source_evaluations=20000,
    cases=[dict(amplitude=a, cells=n, steps=k) for a, n, k in CASES],
    clock_origin="T=s=0 at the quiet wall event", source_continuation="original exact compact reduced field",
    lookup="exact stored x keys; no time interpolation",
    overlap="authenticated original 32-cell rectangle at common Ubar,x; linear/quadratic eta sensitivity",
    spatial_scheme="first-order exterior upwind and regular cap finite volume", time_scheme="fixed-step SSPRK2")
SOURCES = tuple(dict.fromkeys((*previous.SOURCES,
    previous.RELATIVE_OUTPUT+"/rectangle.json", previous.RELATIVE_OUTPUT+"/manifest.json",
    "scripts/record_nsc_atlas_coupled_evolution.py", "tests/test_record_nsc_atlas_coupled_evolution.py",
    *[f"{base}/{name}.py" for name in (
        "nsc_atlas_characteristic_trace", "nsc_atlas_reference_wall_gauge",
        "nsc_atlas_exterior_characteristic_slice", "nsc_atlas_nonlinear_boundary",
        "nsc_atlas_polar_cap", "nsc_atlas_coupled_evolution", "nsc_atlas_coupled_diagnostics")
      for base in ("src/recursive_horizons",)],
    *[f"tests/test_{name}.py" for name in (
        "nsc_atlas_characteristic_trace", "nsc_atlas_reference_wall_gauge",
        "nsc_atlas_exterior_characteristic_slice", "nsc_atlas_nonlinear_boundary",
        "nsc_atlas_hybrid_boundary", "nsc_atlas_polar_cap", "nsc_atlas_coupled_evolution",
        "nsc_atlas_coupled_diagnostics")],
)))
plain, encoded, digest = previous.plain, previous.encoded, previous.digest
statistics, runtime_versions = previous.statistics, previous.runtime_versions
SourceBudget = previous.SourceBudget
_array, _close, _count, _real = previous._array, previous._close, previous._count, previous._real


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def subset(rows, indices):
    return {key: np.asarray(value)[indices] for key, value in rows.items()}


class StoredRows:
    """Exact-key reads of already charged inputs; never interpolate/fallback."""
    def __init__(self, rows):
        self.rows = {key: np.asarray(value) for key, value in rows.items()}
        self.indices = {float(value): i for i, value in enumerate(self.rows["x"])}
        self.calls = 0

    def __call__(self, x):
        if float(x) not in self.indices:
            raise ValueError("requested time is absent from the immutable input cache")
        self.calls += 1
        i = self.indices[float(x)]
        return {key: value[i] for key, value in self.rows.items()}


def final_summary(stage):
    metric = stage["cap_metric"]
    return plain(dict(wall_rates=stage["rates"]["wall"], hybrid=stage["hybrid"],
        wall_f=stage["wall_f"], wall_kinematic_defect=stage["wall_kinematic_defect"],
        wall_value_cell_difference=stage["wall_value_cell_difference"],
        courant_rate=stage["courant_rate"], source_constraint=stage["source_constraint"],
        cap_center=metric["center"], cap_wall=metric["wall"],
        cap_faces={key: metric[key] for key in ("r_faces", "m_faces", "n_faces", "A_faces", "delta_faces", "B_faces")},
        cap_diagnostics=stage["cap_diagnostics"]))


def compact_history(row):
    return {key: row[key] for key in ("E_Ubar", "E_x", "mass", "summaries", "source_constraint",
        "source_constraint_max_abs", "tube_margins", "settings", "scope")}


def _histories(record, case, budget, stride=1):
    samples = case["snapshots"][::stride]
    indices = np.arange(0, 513, 4*stride)
    return coupled_history_diagnostics(samples, eta=np.linspace(0., 1., case["cells"]+1),
        gauge=subset(record["gauge"]["rows"], indices),
        source=subset(record["initial_data"][case["amplitude"]]["snapshot_rows"], indices),
        alpha=case["alpha"], tube=TUBE, budget=budget, source_fronts=record["geometry"]["support_x"])


def _state_tables(case, history):
    return {"exterior": history["states"],
        "wall": {key: np.array([s["state"]["wall"][key] for s in case["snapshots"]])
                 for key in ("r", "L", "c", "T", "s")},
        "cap": {key: np.array([s["state"]["cap"][key] for s in case["snapshots"]])
                for key in ("f", "Pi", "Q")}}


def _restricted(table):
    result = {"exterior": {k: v[:, ::2] for k, v in table["exterior"].items()}, "wall": table["wall"], "cap": {}}
    n = table["cap"]["f"].shape[1]
    y = np.linspace(0., 1., n+1)
    for key, power in (("f", 1), ("Pi", 3), ("Q", 2)):
        w = np.diff(y**power)
        v = table["cap"][key]
        result["cap"][key] = (v[:, ::2]*w[::2]+v[:, 1::2]*w[1::2])/(w[::2]+w[1::2])
    return result


def _difference(a, b):
    return {group: {key: np.asarray(a[group][key])-np.asarray(value) for key, value in b[group].items()}
            for group in b}


def _difference_summary(a, b):
    values = _difference(a, b)
    return {"statistics": {group: {key: statistics(value) for key, value in fields.items()}
                           for group, fields in values.items()}, "wall_arrays": values["wall"]}


def diagnostics(record, budget):
    """Only saved-array derivatives, algebra and analytic vacuum comparisons."""
    original_path = ROOT/previous.RELATIVE_OUTPUT/"rectangle.json"
    original_raw = original_path.read_bytes()
    original = json.loads(original_raw)
    _equal(record["geometry"], original["geometry"], "rectangle overlap geometry")
    histories, tables, controls, overlaps = [], [], [], []
    for case in record["cases"]:
        full, coarsened = _histories(record, case, budget), _histories(record, case, budget, 2)
        histories.append({"retained129": compact_history(full), "coarsened65": compact_history(coarsened)})
        table = _state_tables(case, full)
        tables.append(table)
        x = np.array([row["x"] for row in case["snapshots"]])
        old = original["grids"][1 if case["amplitude"] == "zero" else 3]
        _equal(case["alpha"], original["initial_data"][case["amplitude"]]["alpha"],
               "rectangle overlap amplitude")
        overlap = compare_rectangle_overlap(table["exterior"], x,
            np.linspace(0., 1., case["cells"]+1),
            np.asarray(record["gauge"]["rows"]["W"])[::4], old["x"], old["Ubar"],
            old["states"], budget=budget)
        overlaps.append(overlap)
        quiet = x <= record["geometry"]["support_x"][0]
        control = {"quiet_scalar_maximum": max(float(np.max(abs(table["exterior"]["c"][quiet]))),
            *[float(np.max(abs(value[quiet]))) for value in table["cap"].values()])}
        if case["amplitude"] == "zero":
            gauge = subset(record["gauge"]["rows"], slice(None, None, 4))
            U = gauge["W"][:, None]*np.linspace(0., 1., case["cells"]+1)[None, :]
            product = (1-U/3.2)*(1-x[:, None]/.8)*np.exp(x[:, None]/.8)
            radius = .8*(1+lambertw(-product/math.e).real)
            lapse = np.log(x[:, None]/radius)+(x[:, None]-radius)/.8
            control["vacuum_exterior"] = {key: statistics(table["exterior"][key]-value)
                for key, value in (("r", radius), ("L", lapse))}
            control["vacuum_wall"] = {key: statistics(table["wall"][key]-gauge[key]) for key in ("r", "L", "T", "s")}
        controls.append(control)
    paired = [_difference(tables[2], tables[0]), _difference(tables[3], tables[1])]
    source_defects = {}
    for label, data in record["initial_data"].items():
        exact, snapshot = data["exact_rows"], subset(data["snapshot_rows"], slice(None, None, 2))
        source_defects[label] = {key: statistics((np.asarray(snapshot[key])-np.asarray(exact[key]))[1::2])
            for key in ("h", "h_x", "r", "r_x", "r_xx", "log_C", "log_C_x", "r_Ubar", "log_C_Ubar", "c_Ubar")}
        source_defects[label]["actual_snapshot_constraint_maximum"] = float(np.max(abs(np.asarray(data["snapshot_rows"]["source_constraint"]))))
    return plain(dict(histories=histories, controls=controls,
        original_rectangle_overlap={"reference": previous.RELATIVE_OUTPUT+"/rectangle.json",
            "sha256": digest(original_raw), "cases": overlaps,
            "interpretation": "independent discretization comparison; interpolation sensitivity is not an error enclosure"},
        source_interpolation_withheld_257_odd=source_defects,
        spatial_zero=_difference_summary(_restricted(tables[1]), tables[0]),
        spatial_finite=_difference_summary(_restricted(tables[3]), tables[2]),
        finite_minus_zero=[_difference_summary(tables[j+2], tables[j]) for j in (0, 1)],
        paired_spatial=_difference_summary(_restricted(paired[1]), paired[0]),
        time_512_minus_256=_difference_summary(tables[4], tables[3]),
        comparison_coordinates="common x, exterior eta and cap y; fine cap averages restricted with their declared weights",
        difference_convention="fine minus coarse; paired responses are separately evolved finite minus zero",
        time_history_coarsening="129 versus65 stored slices only; no additional evolution"))


def calculate(budget=None):
    budget = Budget(maximum_evaluations=SETTINGS["global_evaluations"], cpu_seconds=SETTINGS["global_cpu_seconds"]) if budget is None else budget
    phase, record = {"name": "geometry"}, {}
    try:
        geometry = source_geometry(64, budget)
        # This pure three-point check is charged to the global allowance;
        # its optional budget argument is restricted to a source child's
        # 20k/30s ceiling and must not receive the larger parent budget.
        budget.check(3)
        bounds = source_amplitude_bound(geometry)
        x = np.linspace(.4, geometry["x_last"], 513)
        W0 = -3.2/geometry["U0"]*(geometry["Ua"]-geometry["U0"])
        initial = {}
        for label in ("zero", "finite"):
            phase = {"name": "source", "amplitude": label}
            alpha = 0. if label == "zero" else bounds["test_alpha"]
            child = SourceBudget(budget)
            live = ExactSource(geometry, alpha=alpha, budget=child, rtol=2e-10, atol=1e-12)
            exact = live.sample(x[::2])
            bottom = live.quiet_bottom(np.linspace(0., W0, 33))
            child.check()
            child_execution = dict(evaluations=child.evaluations, cpu_seconds=time.process_time()-child.started)
            snapshot = SourceSnapshot(geometry, alpha=alpha, x_nodes=x[::4], source_rows=subset(exact, slice(None, None, 2)), budget=budget)
            cached = snapshot.sample(x)
            initial[label] = plain(dict(alpha=alpha, domain=live.domain, exact_rows=exact,
                bottom=bottom, snapshot_rows=cached, snapshot_summary=snapshot.summary(),
                source_nfev=live.nfev, source_max_step=live.max_step, source_execution=child_execution))
        phase = {"name": "reference_gauge"}
        reference = ReferenceWallGauge(geometry, order=64, budget=budget)
        gauge_rows = reference.sample(x)
        record = dict(schema=SCHEMA, status=STATUS, scope=SCOPE, settings=SETTINGS,
            runtime=runtime_versions(), geometry=geometry, source_bounds=bounds, initial_data=initial,
            gauge=plain(dict(rows=gauge_rows, summary=reference.summary())), cases=[],
            refinement_is_rigorous_error_enclosure=False)
        construction_end = budget.evaluations
        for label, cells, steps in CASES:
            phase = {"name": "coupled_evolution", "amplitude": label, "cells": cells, "steps": steps}
            data = initial[label]
            source, gauge = StoredRows(data["snapshot_rows"]), StoredRows(gauge_rows)
            start = budget.evaluations
            init = quiet_initial_state(subset(data["bottom"], slice(None, None, 32//cells)), source(.4), cap_cells=cells)
            solution = evolve_coupled(init, x[::512//steps], eta=np.linspace(0., 1., cells+1),
                y_faces=np.linspace(0., 1., cells+1), source=source, gauge=gauge,
                alpha=data["alpha"], tube=TUBE, budget=budget, save_every=steps//128,
                maximum_courant=.4, exterior_quadrature_order=2, cap_quadrature_order=4)
            record["cases"].append(plain(dict(amplitude=label, cells=cells, steps=steps, alpha=data["alpha"],
                snapshots=solution["snapshots"], final=final_summary(solution["final_stage"]),
                execution=dict(budget_start=start, budget_end=budget.evaluations,
                    stages=solution["stages"], source_cache_reads=source.calls, gauge_cache_reads=gauge.calls,
                    max_courant=solution["max_courant"], retries=solution["retries"]))))
        phase = {"name": "saved_array_diagnostics"}
        diagnostic_start = budget.evaluations
        record["diagnostics"] = diagnostics(record, budget)
        budget.check()
        record["execution"] = dict(cpu_seconds=time.process_time()-budget.started, evaluations=budget.evaluations,
            construction_evaluations=construction_end, diagnostic_evaluations=budget.evaluations-diagnostic_start,
            source_integrations=2, source_segments=6, gauge_constructions=1, gauge_queries=513,
            snapshot_constructions=2, snapshot_queries=1026, coupled_evolutions=5,
            time_steps=sum(k for _, _, k in CASES), stage_evaluations=sum(2*k+1 for _, _, k in CASES),
            historical_record_replays=0, old_rectangle_solves=0)
        return plain(record)
    except Exception as error:
        error.coupled_failure = plain(dict(phase=phase, exception=type(error).__name__, message=str(error),
            stage=getattr(error, "diagnostic", None), evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started, completed_cases=record.get("cases", []),
            geometry=record.get("geometry"), initial_data=record.get("initial_data"), gauge=record.get("gauge")))
        raise


def _equal(actual, expected, label):
    # Saved data hashes remain exact. Re-evaluated arithmetic is portable at
    # binary64 roundoff, while keys, shapes, integer counters and flags remain
    # exact. Descendant checkers share this validation-only comparison.
    if not previous._same_recomputed(plain(actual), plain(expected)):
        raise ValueError(label+" mismatch")


def _source_snapshot_equal(actual, expected, geometry, alpha, *, source_knots=None):
    """Condition rebuilt diagnostic arithmetic; retain exact query coordinates."""
    label = "retained source interpolant"
    if not isinstance(actual, dict) or set(actual) != set(expected):
        raise ValueError(label+" inventory mismatch")
    tolerances = source_snapshot_roundoff_tolerances(expected, geometry, alpha, source_knots=source_knots)
    for key, value in expected.items():
        a, b = np.asarray(actual[key]), np.asarray(value)
        if (a.dtype.kind not in "iuf" or a.shape != b.shape
                or any(isinstance(v, (bool, np.bool_)) for v in np.asarray(actual[key], dtype=object).flat)
                or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b))):
            raise ValueError(label+" shape/type/nonfinite mismatch: "+key)
        if key == "x":
            if not np.array_equal(a.astype(float).view(np.uint64), b.astype(float).view(np.uint64)):
                raise ValueError(label+" exact coordinates mismatch")
            continue
        allowance = 2e-13+2e-12*abs(b)+tolerances.get(key, 0.)
        error = abs(a-b)
        quiet = (b == 0) & (tolerances[key] == 0) if key in ("c", "c_x", "c_xx") else np.zeros(b.shape, dtype=bool)
        failed = (error > allowance) | (quiet & (a != 0))
        if np.any(failed):
            worst = int(np.argmax(np.where(failed, error/allowance, -1)))
            raise ValueError(label+" mismatch: "+key+
                f"; max_abs_difference={float(np.max(error)):.17g}; "
                f"worst_index={worst}; difference={float(error.flat[worst]):.17g}; "
                f"allowance={float(np.broadcast_to(allowance, b.shape).flat[worst]):.17g}")


def _table(rows, length, label):
    if not isinstance(rows, dict) or not rows:
        raise ValueError(label+" must be a nonempty table")
    return {key: _array(value, (length,), label+"."+key) for key, value in rows.items()}


def _validate_gauge(record):
    rows = _table(record["gauge"]["rows"], 513, "gauge")
    keys = {"x", "U", "V", "W", "W_x", "W_xx", "eta_plus", "eta_plus_x", "r",
        "r_x_total", "r_xx_total", "C", "L", "L_x_total", "r_Ubar", "r_x", "log_C_Ubar",
        "log_C_x", "c", "c_Ubar", "c_x", "s", "s_x", "x_s", "Ubar_s", "T", "T_x",
        "T_s", "r_s", "beta_minus", "beta_plus", "inverse_iterations", "inverse_residual", "normalization_defect"}
    if set(rows) != keys:
        raise ValueError("reference row inventory mismatch")
    x = np.linspace(.4, record["geometry"]["x_last"], 513)
    _equal(rows["x"], x, "gauge coordinates")
    R = rows["r"]
    if R[0] != .45 or R[-1] != .55 or np.any(np.diff(R) <= 0):
        raise ValueError("reference radius domain mismatch")
    g = record["geometry"]
    CJ = -3.2/g["U0"]
    V = (1-x/.8)*np.exp(x/.8)/g["U0"]
    U = (1-R/.8)*np.exp(R/.8)/V
    _close(rows["V"], V, "reference V")
    _close(rows["U"], U, "reference U")
    _close(rows["W"], CJ*(rows["U"]-g["U0"]), "reference W")
    bm, bp = 1/(5*R**2)+3*R/4, 1/(5*R**2)-5*R/4
    v = np.sqrt(bm*bm-1+R*R)
    dv = (25*R/8-3/(10*R**2)-4/(25*R**5))/(2*v)
    dbp = -2/(5*R**3)-5/4
    lam = -x/(.8*(.8-x))
    a = -1/(v*(v-bp))
    Rx = 1.6*lam/a
    Rxx = Rx*((-1/(.8-x)**2)/lam-(-dv/v-(dv-dbp)/(v-bp))*Rx)
    Z = -R*Rx/(.8*(.8-R))-lam
    Zx = -Rx*Rx/(.8-R)**2-R*Rxx/(.8*(.8-R))+1/(.8-x)**2
    Wx, Wxx = CJ*rows["U"]*Z, CJ*rows["U"]*(Z*Z+Zx)
    for key, expected in (("r_x_total", Rx), ("r_xx_total", Rxx), ("W_x", Wx), ("W_xx", Wxx),
        ("eta_plus", -.5*np.log(Wx)), ("eta_plus_x", -.5*Wxx/Wx),
        ("L", np.log(x/R)+(x-R)/.8), ("s_x", Rx/v), ("T_x", bm*Rx/((1-R*R)*v))):
        _close(rows[key], expected, "reference "+key, atol=5e-12)
    _close(rows["C"], np.exp(rows["L"]), "reference lapse")
    rU, rx = rows["C"]*(.8-x)/(4*x), (1-rows["W"]/3.2)*rows["C"]
    for key, expected in (("r_Ubar", rU), ("r_x", rx), ("log_C_Ubar", -(1/R+1/.8)*rU),
        ("log_C_x", 1/x+1/.8-(1/R+1/.8)*rx), ("L_x_total", 1/x+1/.8-(1/R+1/.8)*Rx),
        ("T_s", bm/(1-R*R)), ("x_s", v/Rx), ("Ubar_s", Wx*v/Rx),
        ("beta_minus", bm), ("beta_plus", bp), ("r_s", v)):
        _close(rows[key], expected, "reference "+key, atol=5e-12)
    for key in ("c", "c_Ubar", "c_x"):
        if np.any(rows[key] != 0):
            raise ValueError("reference scalar must be zero")
    _close(rows["normalization_defect"], rows["C"]*Wx/rows["s_x"]**2-1, "reference normalization")
    for key in ("T", "s"):
        if rows[key][0] != 0 or np.any(np.diff(rows[key]) <= 0):
            raise ValueError("reference clock anchor/orientation mismatch")
    if np.max(abs(rows["inverse_residual"])) > 3e-13:
        raise ValueError("reference inverse residual exceeds tolerance")
    if (np.any(rows["inverse_iterations"] != np.floor(rows["inverse_iterations"]))
            or np.any(rows["inverse_iterations"] < 0) or np.any(rows["inverse_iterations"] > 32)):
        raise ValueError("reference inverse iteration count mismatch")
    summary = record["gauge"]["summary"]
    if (summary["query_points"] != 513 or summary["cache_hits"] != 0 or summary["inverse_calls"] != 511
            or summary["quadrature_order"] != 64 or summary["source_integrations"] != 0
            or summary["regional_marches"] != 0 or summary["comparison_is_rigorous_enclosure"] is not False):
        raise ValueError("reference construction inventory mismatch")
    _count(summary["quadrature_nodes"], SETTINGS["global_evaluations"], "reference nodes")
    if summary["quadrature_nodes"] % 64 or summary["evaluations"] != 3+513+summary["quadrature_nodes"]:
        raise ValueError("reference charge mismatch")
    _close(summary["maximum_inverse_residual"], np.max(abs(rows["inverse_residual"])), "reference inverse maximum", atol=0)
    if summary["maximum_inverse_iterations"] != int(np.max(rows["inverse_iterations"])):
        raise ValueError("reference inverse iteration maximum mismatch")
    _equal(summary["domain"], dict(x_a=.4, x_b=x[-1], r_a=.45, r_b=.55,
        W_a=CJ*(g["Ua"]-g["U0"]), W_b=CJ*(g["Ub"]-g["U0"])), "reference domain")
    if summary["schema"] != "nsc-atlas-reference-wall-gauge-v1" or abs(summary["endpoint_integral_defect"]) > 3e-13:
        raise ValueError("reference schema/endpoint mismatch")
    return rows


def _validate_inputs(record, budget):
    # Authentication of the independent comparator uses only its saved
    # arrays. Its file identities are also bound in this record's manifest.
    previous.check(ROOT/previous.RELATIVE_OUTPUT)
    budget.check()
    geometry = record["geometry"]
    validate_geometry(geometry)
    if geometry["order"] != 64 or set(record["initial_data"]) != {"zero", "finite"}:
        raise ValueError("original input inventory mismatch")
    _equal(record["source_bounds"], source_amplitude_bound(geometry), "source-only amplitude bound")
    x = np.linspace(.4, geometry["x_last"], 513)
    W0 = -3.2/geometry["U0"]*(geometry["Ua"]-geometry["U0"])
    for label in ("zero", "finite"):
        data = record["initial_data"][label]
        alpha = 0. if label == "zero" else record["source_bounds"]["test_alpha"]
        if data["alpha"] != alpha or data["source_max_step"] != (x[-1]-x[0])/32:
            raise ValueError("source amplitude or max step mismatch")
        exact = _table(data["exact_rows"], 257, "exact source")
        if set(exact) != previous.SOURCE_KEYS:
            raise ValueError("exact source first-jet inventory mismatch")
        _equal(exact["x"], x[::2], "exact source coordinates")
        # Check every retained exact row's frame against its own h,h_x,
        # without solving the source. This separate snapshot is not fitted
        # to, or substituted for, the declared129-knot source interpolant.
        check_source = SourceSnapshot(geometry, alpha=alpha, x_nodes=x[::2], source_rows=exact, budget=budget)
        budget.check(2*257)
        actual = check_source._frame(x[::2], exact["h"], exact["h_x"], np.zeros(257))
        Q = actual["source_ode_h_xx"]
        original = check_source._frame(x[::2], exact["h"], exact["h_x"], Q)
        for key in ("r_xx", "source_constraint", "source_constraint_jordan"):
            _close(exact[key], original[key], "exact source "+key)
        _close(exact["m"], exact["r"]*(1+4*exact["r_Ubar"]*exact["r_x"]/exact["C"])/2, "source geometric mass")
        snapshot = SourceSnapshot(geometry, alpha=alpha, x_nodes=x[::4], source_rows=subset(exact, slice(None, None, 2)), budget=budget)
        _source_snapshot_equal(data["snapshot_rows"], snapshot.sample(x), geometry, alpha,
                               source_knots=subset(exact, slice(None, None, 2)))
        _equal(data["snapshot_summary"], snapshot.summary(), "snapshot construction summary")
        _equal(data["domain"], snapshot.domain, "source domain")
        bottom = _table(data["bottom"], 33, "quiet bottom")
        if set(bottom) != previous.BOTTOM_KEYS:
            raise ValueError("quiet bottom inventory mismatch")
        ub = np.linspace(0., W0, 33)
        _equal(bottom["Ubar"], ub, "bottom coordinates")
        vacuum = previous.schwarzschild(ub, [.4])
        for key, expected in (("r", vacuum["r"][:, 0]), ("log_C", vacuum["L"][:, 0]), ("C", vacuum["C"][:, 0])):
            _close(bottom[key], expected, "quiet bottom "+key)
        for key in ("c", "c_Ubar", "c_x", "h", "h_x"):
            if np.any(bottom[key] != 0):
                raise ValueError("nonquiet initial bottom")
        _count(data["source_nfev"], 20000, "source RHS calls")
        if data["source_execution"]["evaluations"] != data["source_nfev"]+3+257+33:
            raise ValueError("source child charge mismatch")
        _real(data["source_execution"]["cpu_seconds"], "source CPU", high=30.)
    return _validate_gauge(record)


def validate_record(record, budget=None, *, require_creation=True):
    """Authenticate retained data with an independent bounded check budget."""
    budget = Budget(maximum_evaluations=8000000, cpu_seconds=180.) if budget is None else budget
    authentication_start = budget.evaluations
    try:
        if (record["schema"] != SCHEMA or record["status"] != STATUS or record["scope"] != SCOPE
                or record["refinement_is_rigorous_error_enclosure"] is not False):
            raise ValueError("coupled declaration mismatch")
        _equal(record["settings"], SETTINGS, "fixed calculation settings")
        if set(record["runtime"]) != {"python", "numpy", "scipy"}:
            raise ValueError("runtime inventory mismatch")
        gauge = _validate_inputs(record, budget)
        construction = 960+3+record["gauge"]["summary"]["evaluations"]
        construction += sum(d["source_execution"]["evaluations"]+d["snapshot_summary"]["evaluations"] for d in record["initial_data"].values())
        if len(record["cases"]) != 5:
            raise ValueError("all five ordered cases are required")
        previous_end = construction
        x = np.asarray(gauge["x"])
        for case, (label, n, steps) in zip(record["cases"], CASES):
            data = record["initial_data"][label]
            if (case["amplitude"], case["cells"], case["steps"], case["alpha"]) != (label, n, steps, data["alpha"]):
                raise ValueError("case identity mismatch")
            if len(case["snapshots"]) != 129:
                raise ValueError("each case must retain129 time slices")
            for sample, at in zip(case["snapshots"], x[::4]):
                if sample["x"] != at or set(sample["state"]) != {"increment", "wall", "cap"}:
                    raise ValueError("saved state/time inventory mismatch")
                for group, names, shape in (("increment", ("r", "L", "c"), (n-1,)),
                        ("wall", ("r", "L", "c", "T", "s"), ()), ("cap", ("f", "Pi", "Q"), (n,))):
                    if set(sample["state"][group]) != set(names):
                        raise ValueError("state field inventory mismatch")
                    for key in names:
                        _array(sample["state"][group][key], shape, group+"."+key)
            for key in ("T", "s"):
                values = np.array([sample["state"]["wall"][key] for sample in case["snapshots"]])
                if values[0] != 0 or np.any(np.diff(values) <= 0):
                    raise ValueError("evolved wall clock must increase from quiet zero")
            initial = quiet_initial_state(subset(data["bottom"], slice(None, None, 32//n)),
                                         subset(data["snapshot_rows"], 0), cap_cells=n)
            _equal(case["snapshots"][0]["state"], initial, "initial simultaneous state")
            stage = coupled_stage(case["snapshots"][-1]["state"], eta=np.linspace(0., 1., n+1),
                y_faces=np.linspace(0., 1., n+1), source=subset(data["snapshot_rows"], -1),
                gauge=subset(gauge, -1), alpha=data["alpha"], tube=TUBE, budget=budget,
                exterior_quadrature_order=2, cap_quadrature_order=4, diagnostics=True)
            _equal(case["final"], final_summary(stage), "final algebra/radial reconstruction")
            execution = case["execution"]
            stages = 2*steps+1
            charge = stages*(32*n+7)+3*n
            if (execution["budget_start"] != previous_end or execution["budget_end"] != previous_end+charge
                    or execution["stages"] != stages or execution["source_cache_reads"] != stages+1
                    or execution["gauge_cache_reads"] != stages or execution["retries"] != 0):
                raise ValueError("case cumulative accounting mismatch")
            _real(execution["max_courant"], "maximum Courant", high=.4)
            previous_end += charge
        before = budget.evaluations
        _equal(record["diagnostics"], diagnostics(record, budget), "saved-array diagnostics")
        diagnostic_charges = budget.evaluations-before
        execution = record["execution"]
        expected = dict(evaluations=previous_end+diagnostic_charges, construction_evaluations=construction,
            diagnostic_evaluations=diagnostic_charges, source_integrations=2, source_segments=6,
            gauge_constructions=1, gauge_queries=513, snapshot_constructions=2, snapshot_queries=1026,
            coupled_evolutions=5, time_steps=sum(k for _, _, k in CASES),
            stage_evaluations=sum(2*k+1 for _, _, k in CASES), historical_record_replays=0, old_rectangle_solves=0)
        _equal({key: value for key, value in execution.items() if key != "cpu_seconds"}, expected, "global execution counts")
        _count(execution["evaluations"], 8000000, "global evaluations")
        _real(execution["cpu_seconds"], "global CPU", high=180.)
        budget.check()
        authentication_charges = budget.evaluations-authentication_start
        if require_creation:
            creation = record["creation_execution"]
            expected_counts = dict(calculation_evaluations=execution["evaluations"],
                validation_evaluations=authentication_charges,
                total_evaluations=execution["evaluations"]+authentication_charges)
            _equal({key: value for key, value in creation.items() if key != "cpu_seconds"}, expected_counts,
                   "creation cumulative allowance")
            _count(creation["total_evaluations"], 8000000, "creation evaluations")
            _real(creation["cpu_seconds"], "creation CPU", low=execution["cpu_seconds"], high=180.)
    except (KeyError, TypeError, IndexError, OverflowError) as error:
        raise ValueError("incomplete coupled record") from error
    return {"source_integrations": 0, "time_evolutions": 0, "reference_quadrature": 0,
            "original_rectangle_authentications": 1,
            "final_cap_radial_reconstructions": 5, "final_prefix_wall_reconstructions": 5,
            "authentication_evaluations": budget.evaluations-authentication_start}


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: use a fresh directory or --check")
    bindings, runtime, started = source_bindings(), runtime_versions(), time.process_time()
    budget = Budget(maximum_evaluations=8000000, cpu_seconds=180., started=started)
    try:
        record = calculate(budget)
        if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
            raise RuntimeError("coupled sources or runtime changed during calculation")
        authentication = validate_record(record, budget, require_creation=False)
        elapsed = time.process_time()-started
        if elapsed > 180:
            raise RuntimeError("coupled creation/authentication CPU budget exhausted")
        record["creation_execution"] = dict(calculation_evaluations=record["execution"]["evaluations"],
            validation_evaluations=authentication["authentication_evaluations"],
            total_evaluations=budget.evaluations, cpu_seconds=elapsed)
        budget.check()
    except Exception as error:
        try:
            after = source_bindings()
        except OSError as secondary:
            after = {"unavailable": str(secondary)}
        receipt = dict(schema=SCHEMA+"-failed-attempt", status="FAILED_NOT_COMPLETED", settings=SETTINGS,
            source_bindings_before=bindings, source_bindings_after=after, runtime=runtime,
            failure=getattr(error, "coupled_failure", dict(exception=type(error).__name__, message=str(error))))
        path = ROOT/".build/coupled-evolution-failures"/(uuid.uuid4().hex+".json")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(encoded(previous._failure_safe(receipt)))
        error.failure_receipt = str(path)
        raise
    raw = encoded(record)
    manifest = dict(schema=SCHEMA+"-manifest", source_bindings=bindings,
                    outputs={"coupled.json": dict(sha256=digest(raw), bytes=len(raw))})
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"coupled.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"],
            "creation_execution": record["creation_execution"], "controls": record["diagnostics"]["controls"]}


def check(directory):
    directory = Path(directory)
    if ({path.name for path in directory.iterdir()} != {"coupled.json", "manifest.json"}
            or any(path.is_symlink() for path in directory.iterdir())):
        raise ValueError("coupled output inventory mismatch")
    manifest = json.loads((directory/"manifest.json").read_bytes())
    if manifest.get("schema") != SCHEMA+"-manifest" or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES):
        raise ValueError("coupled source binding/manifest mismatch")
    raw = (directory/"coupled.json").read_bytes()
    if manifest.get("outputs") != {"coupled.json": dict(sha256=digest(raw), bytes=len(raw))}:
        raise ValueError("coupled payload hash/size mismatch")
    return {"checked": str(directory), **validate_record(json.loads(raw))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
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
