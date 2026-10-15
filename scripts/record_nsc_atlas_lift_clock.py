#!/usr/bin/env python3
"""Compare wall clock panels while reusing authenticated outgoing cap jets.

One leading field and two wall solutions share the finite comparison budget.
No cap moment or outgoing ray is recomputed. Authentication rebuilds only
algebra and the declared boundary polynomials from saved inputs.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

import record_nsc_atlas_local_lift as predecessor
from record_nsc_atlas_scalar_inputs import BASE, KNOT_KEYS, digest, encoded, runtime_versions, _plain
from recursive_horizons.nsc_atlas_adaptive_lift import AdaptiveBoundaryLift
from recursive_horizons.nsc_atlas_scalar_lift import wall_boundary_derivative
from recursive_horizons.nsc_atlas_wall_evolution import _gauss_antiderivatives

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-lift-clock-comparison-v1"
SCOPE = "wall clock refinement of the formal boundary lift; no evolved scalar or uniform error enclosure"
INPUT = "results/development/nsc-atlas-local-lift-v1"
PANELS = (16, 32)
SOURCES = tuple(dict.fromkeys((*predecessor.SOURCES, INPUT+"/lift.json", INPUT+"/manifest.json",
    "docs/nsc-atlas-scalar-evolution.md", "scripts/record_nsc_atlas_lift_clock.py",
    "tests/test_record_nsc_atlas_lift_clock.py")))
WALL_KEYS = ("t", "U", "V", "x", "v", "beta_plus", "beta_minus", "Vdot", "w", "j", "p", "A",
             "rho", "rho_s", "V1", "V1_s", "eta", "eta_s", "ds_dt", "reference_q_plus")
BOUNDARY_KEYS = tuple(dict.fromkeys((*KNOT_KEYS, "x", "V0_s", "b_s", "r_U", "r_U_s", "c1_r", "c1_r_s",
    "e", "e_s", "displacement_at_fixed_V", "displacement_at_fixed_V_s", "c0_rr", "c0_rr_s",
    "V0_ss", "backscatter_term", "displacement_term", "d", "d_s", "q_minus_s", "q_minus_ss", "q_plus_s",
    "J_minus", "J_minus_s", "q_minus", "q_plus", "wall_profile")))
INVARIANTS = ("e", "e_s", "r_U", "r_U_s", "c0_rr", "c0_rr_s", "d", "d_s", "q_minus_s", "q_minus_ss", "q_plus_s")


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def input_snapshot():
    path = ROOT/INPUT
    predecessor.check(path)
    raw = {name: (path/name).read_bytes() for name in ("lift.json", "manifest.json")}
    binding = {"path": INPUT, "files": {name: {"sha256": digest(data), "bytes": len(data)} for name, data in raw.items()},
               "source_validation": "record_nsc_atlas_local_lift.check", "physical_solves": 0}
    return json.loads(raw["lift.json"]), binding


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("invalid finite clock datum: "+name)
    return float(value)


def restore_lift(domain, rows, *, geometry=None, budget=None, validation=True):
    """Restore only saved polynomials; optional live geometry/budget supports evolution."""
    if [row["t"] for row in rows] != list(predecessor.DIRECT_TIMES):
        raise ValueError("clock comparison requires the unchanged 81 direct times")
    by_time = {row["t"]: dict(row) for row in rows}
    if geometry is None:
        by_V = {row["V"]: row["U"] for row in rows}
        def wall(V):
            values = np.asarray(V)
            return np.array([by_V[float(v)] for v in values.ravel()]).reshape(values.shape)
        geometry = {"wall": wall}
    if budget is None:
        budget = SimpleNamespace(check=lambda: None)
    jets = SimpleNamespace(domain=domain, budget=budget, geometry=geometry, cache=by_time,
                           sample=lambda t: dict(by_time[t]))
    return AdaptiveBoundaryLift(jets, knot_times=predecessor.FINE_TIMES, nodes=5,
        validation_times=predecessor.WITHHELD if validation else ())


def packed_fit(domain, rows):
    result = _plain(restore_lift(domain, rows).summary())
    result["knot_samples"] = [{key: row[key] for key in KNOT_KEYS} for row in result["knot_samples"]]
    return result


def _clock_data(clock, panels):
    values = np.asarray(clock["values"], dtype=float)
    prefix = np.asarray(clock["prefix"], dtype=float)
    if (set(clock) != {"values", "prefix"} or values.shape != (panels, 4) or prefix.shape != (panels+1,)
            or not np.all(np.isfinite(values)) or not np.all(np.isfinite(prefix)) or prefix[0] != 0):
        raise ValueError("invalid saved four-node wall clock panels")
    _, basis = _gauss_antiderivatives()
    weights = np.array([polynomial(1.) for polynomial in basis])/(2*panels)
    expected = np.r_[0., np.cumsum(values @ weights)]
    if not np.allclose(prefix, expected, rtol=2e-13, atol=2e-14):
        raise ValueError("clock prefix does not integrate its saved panel rates")
    return values, prefix, basis


def rebuild_solution(old, panels, walls, clock):
    """Recompute boundary pullbacks and clock diagnostics, without physical queries."""
    if type(panels) is not int or panels not in PANELS:
        raise ValueError("clock comparison permits only panels16/32")
    if [row["t"] for row in walls] != list(predecessor.DIRECT_TIMES):
        raise ValueError("clock samples must retain all81 original times")
    values, prefix, basis = _clock_data(clock, panels)
    reference = {row["t"]: row for row in old["direct_samples"]}
    rows, diagnostics = [], []
    compatibility = dict.fromkeys(("U", "V", "x", "Vdot", "w", "q_plus", "q_minus", *INVARIANTS), 0.)
    for supplied in walls:
        if set(supplied) != set(WALL_KEYS):
            raise ValueError("wall input inventory mismatch")
        wall = {key: _finite(value, key) for key, value in supplied.items()}
        t, saved = wall["t"], reference[wall["t"]]
        if wall["Vdot"] <= 0 or wall["ds_dt"] <= 0:
            raise ValueError("clock comparison needs positive reference proper-time rates")
        q = wall["reference_q_plus"]
        current = {"U": wall["U"], "V": wall["V"], "x": wall["x"], "Vdot": wall["Vdot"],
                   "w": wall["w"], "q_plus": q, "q_minus": q-2*math.atanh(wall["x"])}
        expected = {**saved, "Vdot": saved["V0_s"], "w": saved["wall_profile"]}
        derivative = wall_boundary_derivative(wall, saved["J_minus"], saved["J_minus_s"])
        current.update({key: derivative[key] for key in INVARIANTS})
        for key, value in current.items():
            compatibility[key] = max(compatibility[key], abs(value-expected[key]))
            if not math.isclose(value, expected[key], rel_tol=2e-11, abs_tol=2e-12):
                raise ValueError("clock refinement changed reference/leading data: "+key)
        if not math.isclose(wall["ds_dt"]*wall["Vdot"], old["domain"]["Vb"]-old["domain"]["Va"], rel_tol=2e-13, abs_tol=2e-14):
            raise ValueError("clock sample changed normalized-time conversion")
        row = {**saved, **derivative, "V0_s": wall["Vdot"],
               "b_U": derivative["b_V"]/saved["U_V"]}
        rows.append({key: row[key] for key in BOUNDARY_KEYS})
        panel = min(int(t*panels), panels-1)
        local = 2*(t-panel/panels)*panels-1
        eta = prefix[panel]+sum(p(local)*r for p, r in zip(basis, values[panel]))/(2*panels)
        eta_t = sum(p.deriv()(local)*r for p, r in zip(basis, values[panel]))
        eta_s = eta_t/wall["ds_dt"]
        V1s = wall["Vdot"]*eta+wall["V"]*eta_s
        if (not math.isclose(eta, wall["eta"], rel_tol=2e-12, abs_tol=2e-13)
                or not math.isclose(wall["V1"], wall["V"]*eta, rel_tol=2e-12, abs_tol=2e-13)
                or not math.isclose(wall["V1_s"], wall["Vdot"]*wall["eta"]+wall["V"]*wall["eta_s"], rel_tol=2e-12, abs_tol=2e-13)):
            raise ValueError("clock value or continuum slope violates its saved definition")
        polynomial_derivative = wall_boundary_derivative({**wall, "V1_s": float(V1s)}, saved["J_minus"], saved["J_minus_s"])
        actual_gap = (derivative["b_s"]-polynomial_derivative["b_s"])/(wall["Vdot"]*saved["U_V"])
        predicted_gap = derivative["r_U"]*derivative["c0_rr"]*wall["v"]/wall["Vdot"]*wall["V"]*(wall["eta_s"]-eta_s)/(wall["Vdot"]*saved["U_V"])
        diagnostics.append({"t": t, "eta_value_defect": float(eta-wall["eta"]),
            "eta_t_polynomial": float(eta_t), "eta_t_continuum": wall["eta_s"]*wall["ds_dt"],
            "eta_t_gap": float(eta_t-wall["eta_s"]*wall["ds_dt"]),
            "V1_s_polynomial": float(V1s), "V1_s_gap": float(V1s-wall["V1_s"]),
            "b_U_clock_gap": float(actual_gap), "predicted_b_U_clock_gap": float(predicted_gap),
            "clock_slope_identity_defect": float(actual_gap-predicted_gap)})
    fit = packed_fit(old["domain"], rows)
    return {"wall_time_panels": panels, "wall_samples": walls, "clock_panels": clock,
        "boundary_samples": rows, "lift": fit, "clock_diagnostics": diagnostics,
        "maximum_absolute_reference_differences": compatibility,
        "maximum_absolute_clock_diagnostics": {key: max(abs(row[key]) for row in diagnostics)
            for key in ("eta_value_defect", "eta_t_gap", "V1_s_gap", "b_U_clock_gap", "clock_slope_identity_defect")}}


def comparison_diagnostics(old, solutions):
    target = old["diagnostics"]["interpolation_target"]
    old_fine = old["solutions"][1]
    arrays = [solution["boundary_samples"] for solution in solutions]
    reference = old["direct_samples"]
    keys = ("b", "b_U", "c1_r", "displacement_at_fixed_V")
    def difference(a, b):
        return {key: max(abs(left[key]-right[key]) for left, right in zip(a, b)) for key in keys}
    return {"interpolation_target": target,
        "sampled_b_target_met": [s["lift"]["maximum_absolute_defects"]["b_defect"] <= target["absolute_b_target"] for s in solutions],
        "boundary_clock16_32_differences": difference(*arrays),
        "boundary_differences_from_clock8": [difference(rows, reference) for rows in arrays],
        "clock8_maximum_absolute_defects": old_fine["maximum_absolute_defects"],
        "maximum_lift_difference": max(abs(solutions[0]["lift"]["terminal_lift"]-solutions[1]["lift"]["terminal_lift"]),
            *(abs(a["B_lift"]-b["B_lift"]) for a, b in zip(solutions[0]["lift"]["validation_samples"], solutions[1]["lift"]["validation_samples"]))),
        "maximum_lift_differences_from_clock8": [max(abs(s["lift"]["terminal_lift"]-old_fine["terminal_lift"]),
            *(abs(a["B_lift"]-b["B_lift"]) for a, b in zip(s["lift"]["validation_samples"], old_fine["validation_samples"]))) for s in solutions]}


def _prepare_shared():
    from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, solve_prepared_field
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
    from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
    from recursive_horizons.nsc_atlas_wedge_field import WedgeField
    shared = PreparedFieldBudget()
    metrics = BulkBudget(started=shared.started, external=shared)
    geometry = source_geometry()
    field = solve_prepared_field(BASE["field_nodes"], BASE["geometry_knots"], budget=shared)
    return shared, metrics, geometry, field, WedgeField(field)


def compare():
    from recursive_horizons.nsc_atlas_wall_evolution import solve_wall_evolution
    old, binding = input_snapshot()
    shared, metrics, geometry, field, wedge = _prepare_shared()
    solutions = []
    for panels in PANELS:
        settings = {key.removeprefix("wall_"): value for key, value in BASE.items() if key.startswith("wall_")}
        settings["time_panels"] = panels
        wall = solve_wall_evolution(wedge, geometry, budget=metrics, source_scale=1., **settings)
        rows = []
        for t in predecessor.DIRECT_TIMES:
            raw = wall.sample(t)
            raw["reference_q_plus"] = float(field.geometry["arrival"](raw["V"]))
            rows.append({key: float(raw[key]) for key in WALL_KEYS})
        clock = {"values": wall.clock_values.tolist(), "prefix": wall.clock_prefix.tolist()}
        solutions.append(rebuild_solution(old, panels, rows, clock))
        shared.check()
        metrics.check()
    return _plain({"schema": SCHEMA, "scope": SCOPE, "input": binding, "domain": old["domain"],
        "runtime": runtime_versions(), "comparison_base": dict(BASE), "varied_quantity": "wall_time_panels",
        "settings": [{**BASE, "wall_time_panels": panels} for panels in PANELS],
        "solutions": solutions, "diagnostics": comparison_diagnostics(old, solutions),
        "reused_ray_quantities": ["J_minus", "J_minus_s"], "uniform_error_enclosure": False,
        "execution": {"cpu_seconds": time.process_time()-shared.started, "cpu_seconds_max": 120.,
            "scalar_kernel_values": shared.kernel_values, "scalar_kernel_values_max": 30_000_000,
            "metric_values": metrics.metric_values, "source_quadrature_nodes": metrics.source_nodes,
            "field_rhs_calls": field.rhs_calls, "wall_rhs_calls": metrics.rhs_calls,
            "field_solves": 1, "wall_solves": 2, "reused_boundary_jets": 81,
            "new_cap_queries": 0, "new_ray_quadrature_nodes": 0}})


def validate_record(record):
    old, binding = input_snapshot()
    if (record.get("schema") != SCHEMA or record.get("scope") != SCOPE or record.get("input") != binding
            or record.get("domain") != old["domain"] or record.get("comparison_base") != BASE
            or record.get("varied_quantity") != "wall_time_panels"
            or record.get("settings") != [{**BASE, "wall_time_panels": p} for p in PANELS]
            or record.get("reused_ray_quantities") != ["J_minus", "J_minus_s"]
            or record.get("uniform_error_enclosure") is not False):
        raise ValueError("clock comparison input, scope or settings mismatch")
    solutions = record["solutions"]
    if len(solutions) != 2:
        raise ValueError("clock comparison requires exactly two solutions")
    for panels, solution in zip(PANELS, solutions):
        expected = rebuild_solution(old, panels, solution["wall_samples"], solution["clock_panels"])
        if not predecessor._same_saved(solution, expected):
            raise ValueError("clock boundary algebra or saved polynomial mismatch")
    if not predecessor._same_saved(record["diagnostics"], comparison_diagnostics(old, solutions)):
        raise ValueError("clock comparison diagnostics mismatch")
    cost = record["execution"]
    required = dict(cpu_seconds_max=120., scalar_kernel_values_max=30_000_000,
        field_solves=1, wall_solves=2, reused_boundary_jets=81, new_cap_queries=0, new_ray_quadrature_nodes=0)
    if any(isinstance(cost.get(key), bool) or cost.get(key) != value for key, value in required.items()):
        raise ValueError("clock comparison solve counts or ceilings changed")
    if not 0 <= _finite(cost["cpu_seconds"], "CPU cost") <= 120:
        raise ValueError("clock comparison CPU budget exceeded")
    for key, limit in (("scalar_kernel_values", 30_000_000), ("metric_values", 30_000_000),
                       ("source_quadrature_nodes", 30_000_000), ("field_rhs_calls", 6000), ("wall_rhs_calls", 12000)):
        if type(cost.get(key)) is not int or not 0 <= cost[key] <= limit:
            raise ValueError("clock comparison counter exceeded: "+key)
    if cost["metric_values"]+cost["source_quadrature_nodes"] > 30_000_000:
        raise ValueError("clock comparison cumulative metric budget exceeded")
    runtime = record.get("runtime")
    if (not isinstance(runtime, dict) or set(runtime) != {"python", "numpy", "scipy"}
            or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
        raise ValueError("clock comparison runtime metadata mismatch")
    return True


def _summary(record):
    return {"withheld_maxima": [s["lift"]["maximum_absolute_defects"] for s in record["solutions"]],
        "sampled_b_target_met": record["diagnostics"]["sampled_b_target_met"],
        "clock_slope_gaps": [s["maximum_absolute_clock_diagnostics"] for s in record["solutions"]],
        "maximum_lift_difference": record["diagnostics"]["maximum_lift_difference"]}


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only clock record: use a fresh directory or --check")
    sources, runtime = source_bindings(), runtime_versions()
    record = compare()
    if sources != source_bindings() or runtime != runtime_versions() or runtime != record["runtime"]:
        raise RuntimeError("clock comparison sources or runtime changed")
    validate_record(record)
    raw = encoded(record)
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"clock.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded({"schema": SCHEMA+"-manifest", "source_bindings": sources,
        "outputs": {"clock.json": {"sha256": digest(raw), "bytes": len(raw)}}}))
    return {"created": str(directory), "execution": record["execution"], "summary": _summary(record)}


def check(directory):
    directory = Path(directory)
    if {p.name for p in directory.iterdir()} != {"clock.json", "manifest.json"}:
        raise ValueError("clock record inventory mismatch")
    if any((directory/name).is_symlink() for name in ("clock.json", "manifest.json")):
        raise ValueError("clock evidence must not use symlinks")
    manifest = json.loads((directory/"manifest.json").read_bytes())
    raw = (directory/"clock.json").read_bytes()
    if (manifest.get("schema") != SCHEMA+"-manifest" or manifest.get("source_bindings") != source_bindings()
            or manifest.get("outputs") != {"clock.json": {"sha256": digest(raw), "bytes": len(raw)}}):
        raise ValueError("clock record source or output binding mismatch")
    record = json.loads(raw)
    validate_record(record)
    return {"checked": str(directory), "physical_solves": 0, "summary": _summary(record)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
