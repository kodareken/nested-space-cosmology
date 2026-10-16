#!/usr/bin/env python3
"""Record four bounded scalar evolutions using authenticated saved inputs.

The check path restores coefficient polynomials and checks retained discrete
identities. It never evolves a field or silently regenerates missing inputs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
from validation_source_successors import verify_saved_bindings, verify_saved_source

import record_nsc_atlas_forcing_cache as cache_record
import record_nsc_atlas_lift_clock as clock_record
from record_nsc_atlas_scalar_inputs import _plain, _real, _count, _equal, digest, encoded, runtime_versions

ROOT = Path(__file__).resolve().parents[1]
CACHE = "results/development/nsc-atlas-forcing-cache-v1"
CLOCK = "results/development/nsc-atlas-lift-clock-v1"
SCHEMA = "nsc-atlas-scalar-evolution-comparison-v1"
STATUS = "FORMAL_EXTERIOR_SCALAR_CORRECTION_REFINEMENT"
SCOPE = "formal exterior cubic scalar coefficient; sampled sensitivity, no central waveform or uniform remainder"
TIMES = tuple(i/16 for i in range(17))
CASES = ((17, 129), (33, 129), (33, 65), (33, 129))
ODE = dict(tolerance=2e-6, atol=1e-8, max_step=1/64, rhs_calls_max=6000)
ODE_CASES = (ODE, ODE, ODE, {**ODE, "tolerance": 2e-7, "atol": 1e-9})
Q = .13882099234357653
ROW_KEYS = ("t", "V", "U_wall", "U_wall_V", "g1", "g1_V", "B_wall", "b_wall",
    "z_wall", "y_wall", "z_wall_V", "y_wall_V", "source_defect", "source_rate_defect",
    "wall_flux", "inactive_state_spread")
READOUTS = ("z_wall", "y_wall", "z_wall_V", "y_wall_V")
SOURCES = tuple(dict.fromkeys((*cache_record.SOURCES, *clock_record.SOURCES,
    CACHE+"/cache.json", CACHE+"/manifest.json", CLOCK+"/clock.json", CLOCK+"/manifest.json",
    "src/recursive_horizons/nsc_atlas_scalar_evolution.py", "tests/test_nsc_atlas_scalar_evolution.py",
    "scripts/record_nsc_atlas_scalar_evolution.py", "tests/test_record_nsc_atlas_scalar_evolution.py")))


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def input_snapshot():
    cache_record.check(ROOT/CACHE)
    clock_record.check(ROOT/CLOCK)
    cache = json.loads((ROOT/CACHE/"cache.json").read_bytes())
    clock = json.loads((ROOT/CLOCK/"clock.json").read_bytes())
    if (clock["solutions"][1]["wall_time_panels"] != 32
            or clock["diagnostics"]["sampled_b_target_met"][1] is not True):
        raise ValueError("scalar propagation requires the clock32 lift to pass its unchanged sampled target")
    if any(cache["domain"][key] != clock["domain"][key] for key in cache_record.DOMAIN_KEYS):
        raise ValueError("retained forcing and clock lift domains disagree")
    files = (CACHE+"/cache.json", CACHE+"/manifest.json", CLOCK+"/clock.json", CLOCK+"/manifest.json")
    return cache, clock, {"files": {name: digest((ROOT/name).read_bytes()) for name in files},
        "clock_time_panels": 32, "clock_sampled_b_target_met": True,
        "restored_forcing_slices": 65, "restored_boundary_jets": 81,
        "new_leading_field_solves": 0, "new_wall_solves": 0}


def restore_inputs():
    from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, _geometry
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
    from recursive_horizons.nsc_atlas_forcing_cache import ForcingSliceProvider, ForcingTimeCache
    shared = PreparedFieldBudget()
    metrics = BulkBudget(started=shared.started, external=shared)
    cache, clock, receipt = input_snapshot()
    geometry = _geometry(257)
    d = cache["domain"]
    for key in ("Ua", "Ub", "Va", "Vb"):
        if abs(geometry[key]-d[key]) > 8*np.finfo(float).eps*max(1., abs(d[key])):
            raise ValueError("restored geometry endpoints differ from retained domain")
    defects = dict(wall_coordinate=0., wall_derivative=0., arrival_coordinate=0.)
    for row in cache["direct_samples"]:
        defects["wall_coordinate"] = max(defects["wall_coordinate"], abs(float(geometry["wall"](row["V"]))-row["U_wall"]))
    boundary = clock["solutions"][1]["boundary_samples"]
    for row in boundary:
        defects["wall_coordinate"] = max(defects["wall_coordinate"], abs(float(geometry["wall"](row["V"]))-row["U"]))
        defects["wall_derivative"] = max(defects["wall_derivative"], abs(float(geometry["wall"](row["V"], 1))-row["U_V"]))
        defects["arrival_coordinate"] = max(defects["arrival_coordinate"], abs(float(geometry["arrival"](row["V"]))-row["q_plus"]))
    if any(value > 2e-13 for value in defects.values()):
        raise ValueError("full reference geometry differs from retained input coordinates")
    wedge = SimpleNamespace(domain=dict(d), geometry=geometry, budget=shared)
    provider = ForcingSliceProvider.from_records(wedge, cache["source_geometry"], cache["direct_samples"], budget=metrics)
    forcing = {count: ForcingTimeCache(wedge, cache["source_geometry"], np.linspace(0, 1, 129),
        np.linspace(0, 1, count), budget=metrics, slice_provider=provider) for count in (17, 33)}
    lift = clock_record.restore_lift(d, boundary, geometry=geometry, budget=metrics, validation=False)
    shared.check()
    return SimpleNamespace(shared=shared, budget=metrics, domain=d, geometry=geometry, forcing=forcing,
        lift=lift, cache=cache, clock=clock, receipt=receipt, geometry_defects=defects)


def input_indicators(cache):
    d = cache["domain"]
    du, dv = d["Ub"]-d["U0"], d["Vb"]-d["Va"]
    values = []
    for row in cache["diagnostics"]["maximum_absolute_defects"]:
        f = max(row["temporal_F_defect"], row["combined_F_defect"])
        g = max(row["temporal_G_defect"], row["combined_G_defect"])
        values.append((row["g1_defect"]+2*dv*f+du*dv*g)/(1-Q))
    pair = dict(F=0., G=0., g1=0.)
    for left, right in zip(*(solution["validation_samples"] for solution in cache["solutions"])):
        for key in pair:
            pair[key] = max(pair[key], float(np.max(abs(np.asarray(left["interpolated"][key])-
                                                      np.asarray(right["interpolated"][key])))))
    return {"Q": Q, "formula": "(d_g+2*Delta_V*d_F+Delta_U*Delta_V*d_G)/(1-Q)",
        "source": "docs/nsc-atlas-scalar-forcing.md", "lift_held_fixed": True,
        "cache_to_direct_sampled_y_indicators": values, "cache_pair_sampled_input_differences": pair,
        "cache_pair_sampled_y_indicator": (pair["g1"]+2*dv*pair["F"]+du*dv*pair["G"])/(1-Q),
        "uniform_bound": False, "scope": "sampled reduced-field sensitivity; not a wall-derivative bound"}


def comparison_diagnostics(solutions):
    axes = []
    for name, first, second in (("cache_time", 0, 1), ("correction_grid", 2, 1), ("ODE_tolerance", 1, 3)):
        left, right = (solutions[index]["wall_samples"] for index in (first, second))
        differences = {key: max(abs(a[key]-b[key]) for a, b in zip(left, right)) for key in READOUTS}
        peak = max(abs(row["y_wall"]) for row in right)
        axes.append({"axis": name, "solution_indices": [first, second], "maximum_absolute_differences": differences,
            "fine_sampled_peak_y": peak, "sampled_y_difference_over_peak": differences["y_wall"]/peak if peak else None})
    cache_difference = axes[0]["maximum_absolute_differences"]["y_wall"]
    ode_difference = axes[2]["maximum_absolute_differences"]["y_wall"]
    return {"axes": axes, "sampled_readout_times": list(TIMES), "uniform_error_enclosure": False,
            "source_or_lift_changes_between_axes": False, "ODE_tolerance_refined": True,
            "ODE_y_difference_over_cache_y_difference": ode_difference/cache_difference if cache_difference else None}


def compare():
    from recursive_horizons.nsc_atlas_scalar_evolution import solve_scalar_correction
    restored = restore_inputs()
    solutions = []
    for (count, nodes), ode in zip(CASES, ODE_CASES):
        solution = solve_scalar_correction(restored.forcing[count], restored.lift,
            np.linspace(restored.domain["U0"], restored.domain["Ub"], nodes), **ode)
        saved = solution.summary(TIMES)
        saved["cache_times"] = count
        saved["interior_samples"] = [_plain(solution.bulk(np.minimum(solution.grid, row["U_wall"]), row["t"]))
                                     for row in saved["wall_samples"]]
        solutions.append(_plain(saved))
    restored.shared.check()
    restored.budget.check()
    return {"schema": SCHEMA, "status": STATUS, "scope": SCOPE, "runtime": runtime_versions(),
        "domain": restored.domain, "source_scale": 1., "geometry_knots": 257, "clock_time_panels": 32,
        "inputs": restored.receipt, "geometry_restoration_defects": restored.geometry_defects,
        "solutions": solutions, "diagnostics": comparison_diagnostics(solutions),
        "sampled_input_indicators": input_indicators(restored.cache),
        "refinement_is_rigorous_error_enclosure": False,
        "execution": {"cpu_seconds": time.process_time()-restored.shared.started,
            "cpu_seconds_max": restored.shared.cpu_seconds, "scalar_kernel_values": restored.shared.kernel_values,
            "scalar_kernel_values_max": restored.shared.kernel_values_max,
            "rhs_calls": sum(row["rhs_calls"] for row in solutions),
            "readout_rhs_calls": sum(row["readout_rhs_calls"] for row in solutions),
            "correction_solves": 4, "leading_field_solves": 0, "wall_solves": 0,
            "new_forcing_slices": 0, "new_boundary_jets": 0}}


def _validate_solution(saved, restored, count, nodes, ode):
    from recursive_horizons.nsc_atlas_scalar_evolution import SCHEMA as SOLVER_SCHEMA, SCOPE as SOLVER_SCOPE, reference_kernel
    if saved["schema"] != SOLVER_SCHEMA or saved["scope"] != SOLVER_SCOPE or saved["cache_times"] != count:
        raise ValueError("scalar evolution solution declaration mismatch")
    cache_record._same(saved["settings"], dict(grid_nodes=nodes, **ode), "solver settings")
    cache_record._same(saved["domain"], restored.domain, "solution domain")
    rows, interiors = saved["wall_samples"], saved["interior_samples"]
    if (not isinstance(rows, list) or len(rows) != 17 or [row["t"] for row in rows] != list(TIMES)
            or not isinstance(interiors, list) or len(interiors) != 17):
        raise ValueError("scalar evolution readout inventory mismatch")
    d = restored.domain
    grid = np.linspace(d["U0"], d["Ub"], nodes)
    forcing, lift = restored.forcing[count], restored.lift
    for row, interior in zip(rows, interiors):
        if set(row) != set(ROW_KEYS) or set(interior) != {"U", "V", "z", "y"}:
            raise ValueError("scalar evolution readout fields mismatch")
        for key, value in row.items():
            _real(value, key)
        for key, values in interior.items():
            if not isinstance(values, list) or len(values) != nodes:
                raise ValueError("scalar evolution interior shape mismatch")
            for value in values:
                _real(value, key)
        t = row["t"]
        V = min(d["Vb"], max(d["Va"], d["Va"]+(d["Vb"]-d["Va"])*t))
        wall = min(d["Ub"], max(d["Ua"], float(restored.geometry["wall"](V))))
        wall_V = float(restored.geometry["wall"](V, 1))
        U = np.minimum(grid, wall)
        cache_record._same(interior["U"], U.tolist(), "retained interior coordinates")
        cache_record._same(interior["V"], np.full(nodes, V).tolist(), "retained interior time")
        supplied = forcing.eval(U, t)
        z, y = np.asarray(interior["z"]), np.asarray(interior["y"])
        B = np.asarray(lift.B(U))
        cache_record._same(y.tolist(), (z+supplied["g1"]+B).tolist(), "retained lifted field")
        if (z[0] != 0. or any(value != z[-1] for u, value in zip(U, z) if u == wall)
                or row["source_defect"] != 0. or row["source_rate_defect"] != 0.
                or row["inactive_state_spread"] < 0):
            raise ValueError("scalar evolution source or clipped continuation mismatch")
        integrand = reference_kernel(U, V)*y+supplied["G"]
        # Reconstruct the finite-grid equation from saved values, without
        # integrating an ODE or differentiating a saved waveform.
        prefix = float(np.sum((integrand[:-1]+integrand[1:])*np.diff(U)/2))
        zV = float(supplied["F"][-1]-supplied["Fsource"]-prefix)
        b = float(lift.b_at_U(wall))
        expected = dict(V=V, U_wall=wall, U_wall_V=wall_V, g1=supplied["g1"], g1_V=supplied["g1_V"],
            B_wall=float(B[-1]), b_wall=b, z_wall=float(z[-1]), y_wall=float(y[-1]), z_wall_V=zV,
            y_wall_V=float(supplied["g1_V"]+zV+wall_V*b), wall_flux=float(supplied["F"][-1]))
        if any(not _equal(row[key], value) for key, value in expected.items()):
            raise ValueError("scalar evolution saved source, lift or discrete wall-rate identity mismatch")
        if t <= .25 and max(abs(z).max(), abs(y).max(), *(abs(row[key]) for key in READOUTS)) > 1e-24:
            raise ValueError("scalar evolution lost the quiet interval")
    for key, row_key in (("maximum_absolute_source_defect", "source_defect"),
                         ("maximum_absolute_source_rate_defect", "source_rate_defect"),
                         ("maximum_absolute_wall_flux", "wall_flux")):
        if saved[key] != max(abs(row[row_key]) for row in rows):
            raise ValueError("scalar evolution defect maximum mismatch")
    _count(saved["rhs_calls"], "RHS calls", 6000)
    if saved["rhs_calls"] < 1 or type(saved["readout_rhs_calls"]) is not int or saved["readout_rhs_calls"] != 17:
        raise ValueError("scalar evolution readout or solve count mismatch")
    _count(saved["kernel_values"], "kernel values", 30_000_000)
    if saved["kernel_values"] != nodes*(saved["rhs_calls"]+17):
        raise ValueError("scalar evolution kernel accounting mismatch")
    cuts = saved["split_times"]
    if not isinstance(cuts, list) or not 2 <= len(cuts) <= 513:
        raise ValueError("scalar evolution split inventory mismatch")
    for value in cuts:
        _real(value, "split time")
    if cuts[0] != 0 or cuts[-1] != 1 or any(b <= a for a, b in zip(cuts, cuts[1:])):
        raise ValueError("scalar evolution splits do not cover the interval")
    _count(saved["segments"], "segments", 512)
    if saved["segments"] != len(cuts)-1:
        raise ValueError("scalar evolution segment count mismatch")
    for t in (*np.linspace(0, 1, count), *lift.knot_times, .25, .75):
        if min(abs(cut-t) for cut in cuts) > 8*np.finfo(float).eps:
            raise ValueError("scalar evolution omitted a retained interpolation knot")
    wall_cuts = restored.geometry["wall"](d["Va"]+(d["Vb"]-d["Va"])*np.asarray(cuts))
    for U in grid[(grid > d["Ua"]) & (grid < d["Ub"])]:
        if min(abs(wall_cuts-U)) > 2e-13:
            raise ValueError("scalar evolution omitted a moving-wall cell birth")
    _real(saved["maximum_split_merge_width"], "split merge width", True)
    if saved["maximum_split_merge_width"] > 8*np.finfo(float).eps:
        raise ValueError("scalar evolution merged distinct split times")


def validate_record(record):
    try:
        if (record.get("schema") != SCHEMA or record.get("status") != STATUS or record.get("scope") != SCOPE
                or record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)
                or type(record.get("geometry_knots")) is not int or record["geometry_knots"] != 257
                or type(record.get("clock_time_panels")) is not int or record["clock_time_panels"] != 32):
            raise ValueError("scalar evolution declaration mismatch")
        runtime = record["runtime"]
        if (set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("scalar evolution runtime mismatch")
        restored = restore_inputs()
        cache_record._same(record["domain"], restored.domain, "restored domain")
        cache_record._same(record["inputs"], restored.receipt, "input authentication")
        cache_record._same(record["geometry_restoration_defects"], restored.geometry_defects, "restored geometry")
        if not isinstance(record["solutions"], list) or len(record["solutions"]) != 4:
            raise ValueError("scalar evolution requires four fixed-axis solutions")
        for saved, (count, nodes), ode in zip(record["solutions"], CASES, ODE_CASES):
            _validate_solution(saved, restored, count, nodes, ode)
        cache_record._same(record["sampled_input_indicators"], input_indicators(restored.cache), "sampled input indicator")
        cache_record._same(record["diagnostics"], comparison_diagnostics(record["solutions"]), "axis comparison")
        cost = record["execution"]
        _real(cost["cpu_seconds"], "CPU time", True)
        _real(cost["cpu_seconds_max"], "CPU ceiling")
        _count(cost["scalar_kernel_values_max"], "kernel ceiling")
        if cost["cpu_seconds"] > 120 or cost["cpu_seconds_max"] != 120 or cost["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("scalar evolution budget changed or exceeded")
        for key, maximum in (("scalar_kernel_values", 30_000_000), ("rhs_calls", 24000), ("readout_rhs_calls", 68)):
            _count(cost[key], key, maximum)
        for key in ("rhs_calls", "readout_rhs_calls"):
            if cost[key] != sum(solution[key] for solution in record["solutions"]):
                raise ValueError("scalar evolution cumulative RHS accounting mismatch")
        if cost["scalar_kernel_values"] != sum(solution["kernel_values"] for solution in record["solutions"]):
            raise ValueError("scalar evolution cumulative kernel accounting mismatch")
        for key, expected in (("correction_solves", 4), ("leading_field_solves", 0), ("wall_solves", 0),
                              ("new_forcing_slices", 0), ("new_boundary_jets", 0)):
            if type(cost[key]) is not int or cost[key] != expected:
                raise ValueError("scalar evolution restored-input accounting mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete scalar evolution record") from error
    return True


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    bindings, runtime = source_bindings(), runtime_versions()
    record = compare()
    if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
        raise RuntimeError("calculation inputs or runtime changed during comparison")
    validate_record(record)
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings,
        "outputs": {"evolution.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"evolution.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"], "diagnostics": record["diagnostics"]}


def check(directory):
    directory = Path(directory)
    path = directory/"manifest.json"
    manifest = json.loads(path.read_bytes())
    if (path.is_symlink() or manifest.get("schema") != SCHEMA+"-manifest"
            or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)):
        raise ValueError("scalar evolution source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "evolution.json"}
            or set(manifest.get("outputs", {})) != {"evolution.json"}):
        raise ValueError("scalar evolution inventory mismatch")
    path = directory/"evolution.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["evolution.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("scalar evolution output mismatch")
    validate_record(json.loads(raw))
    return {"checked": str(directory), "physical_solves": 0, "propagated_scalar_solves": 0,
            "saved_discrete_readouts_checked": 68}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
