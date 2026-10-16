#!/usr/bin/env python3
"""Create/check a five-case successor using authenticated, unchanged v1 inputs.

Creation reuses source/gauge caches and counts inherited authentication inside
the same allowance. Checking reconstructs saved fields, final radial/prefix
stages and diagnostics, but never integrates a source, gauge or time evolution.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import time
import uuid

import numpy as np
from recursive_horizons.nsc_atlas_validation import collect_departure_replay
from validation_source_successors import verify_saved_bindings, verify_saved_source

import record_nsc_atlas_coupled_evolution as previous
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_compatible_evolution import (
    quiet_departures, compatible_stage, evolve_compatible,
)
from recursive_horizons.nsc_atlas_coupled_diagnostics import (
    coupled_history_diagnostics, compare_rectangle_overlap,
)
from recursive_horizons.nsc_atlas_departure_diagnostics import departure_history_diagnostics
from recursive_horizons.nsc_atlas_reduced_cap import (
    reduced_cap_diagnostics, _reconstruct, _wall_conversion,
)
from recursive_horizons.nsc_atlas_vacuum_departure import vacuum_slice


ROOT = Path(__file__).resolve().parents[1]
RELATIVE_OUTPUT = "results/development/nsc-atlas-compatible-evolution-v1"
SCHEMA = "nsc-atlas-compatible-evolution-v1"
STATUS = "FINITE_AMPLITUDE_COMPATIBLE_WINDOW"
SCOPE = ("original incoming-label window with vacuum departures and a conservative reduced cap; "
         "sampled discretization and constraint diagnostics, not a uniform error enclosure, "
         "stability result or post-wall central-arrival evolution")
CASES, TUBE = previous.CASES, previous.TUBE
SETTINGS = dict(cases=[dict(amplitude=a, cells=n, steps=k) for a, n, k in CASES],
    input=previous.RELATIVE_OUTPUT+"/coupled.json", input_cache_nodes=513,
    saved_time_nodes=129, exterior_quadrature_order=2, cap_quadrature_order=4,
    maximum_courant=.4, tube=TUBE, global_cpu_seconds=180., global_evaluations=8000000,
    source_continuation=previous.SETTINGS["source_continuation"],
    clock_origin=previous.SETTINGS["clock_origin"], lookup=previous.SETTINGS["lookup"],
    spatial_scheme="analytic vacuum plus upwind departures; conservative reduced cap characteristics",
    time_scheme="fixed-step SSPRK2", new_source_integrations=0, new_reference_quadratures=0,
    cap_restriction="sum adjacent physical dr integrals Ip and Iq at common y cells",
    history_diagnostics="raw physical finite differences and analytic-vacuum departure differences separately")
SOURCES = tuple(dict.fromkeys((*previous.SOURCES,
    previous.RELATIVE_OUTPUT+"/coupled.json", previous.RELATIVE_OUTPUT+"/manifest.json",
    "scripts/record_nsc_atlas_compatible_evolution.py", "tests/test_record_nsc_atlas_compatible_evolution.py",
    *[f"src/recursive_horizons/{name}.py" for name in (
        "nsc_atlas_vacuum_departure", "nsc_atlas_reduced_cap", "nsc_atlas_compatible_evolution",
        "nsc_atlas_departure_diagnostics")],
    *[f"tests/test_{name}.py" for name in (
        "nsc_atlas_vacuum_departure", "nsc_atlas_reduced_cap", "nsc_atlas_compatible_evolution",
        "nsc_atlas_departure_diagnostics")],
)))
plain, encoded, digest = previous.plain, previous.encoded, previous.digest
statistics, runtime_versions = previous.statistics, previous.runtime_versions
subset, StoredRows = previous.subset, previous.StoredRows
_equal, _array, _count, _real, _close = (previous._equal, previous._array, previous._count,
                                      previous._real, previous._close)


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def load_parent(budget):
    """Authenticate v1 with this caller's live allowance, not a fresh budget."""
    directory = ROOT/previous.RELATIVE_OUTPUT
    if ({p.name for p in directory.iterdir()} != {"coupled.json", "manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("inherited coupled inventory mismatch")
    raw, manifest_raw = (directory/"coupled.json").read_bytes(), (directory/"manifest.json").read_bytes()
    manifest = json.loads(manifest_raw)
    if (manifest.get("schema") != previous.SCHEMA+"-manifest"
            or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=previous.SOURCES)
            or manifest.get("outputs") != {"coupled.json": dict(sha256=digest(raw), bytes=len(raw))}):
        raise ValueError("inherited coupled hash/source binding mismatch")
    parent = json.loads(raw)
    before = budget.evaluations
    authentication = previous.validate_record(parent, budget)
    charges = budget.evaluations-before
    if charges != authentication["authentication_evaluations"]:
        raise ValueError("inherited authentication accounting mismatch")
    identity = dict(path=previous.RELATIVE_OUTPUT+"/coupled.json", sha256=digest(raw),
        manifest_sha256=digest(manifest_raw), authentication=authentication)
    return parent, identity, charges


def final_summary(stage, budget):
    metric, prefix, rhs = stage["cap_metric"], stage["exterior_prefix"], stage["cap_rhs"]
    return plain(dict(physical_wall=stage["physical_wall"], wall_rates=stage["wall_rates"],
        wall_departure_rates=stage["rates"]["wall"], hybrid=stage["hybrid"],
        wall_kinematic_defect=stage["wall_kinematic_defect"],
        canonical_wall_flux_defect=stage["canonical_wall_flux_defect"],
        conformal_wall_flux_defect=stage["conformal_wall_flux_defect"],
        vacuum_wall_anchor_defect=stage["vacuum_wall_anchor_defect"],
        courant_rate=stage["courant_rate"], cap_center=metric["center"], cap_wall=metric["wall"],
        cap_faces={k: metric[k] for k in ("r_faces", "m_faces", "n_faces", "A_faces", "delta_faces", "B_faces", "chi_faces")},
        cap_diagnostics=reduced_cap_diagnostics(metric, budget=budget),
        cap_flux={k: rhs[k] for k in ("p_flux", "q_flux", "reaction_integral", "wall_chi_x",
                                      "wall_f_x", "wall_c_x", "constructed_wall_flux_defect")},
        removed_vacuum_quadrature={k: prefix[k] for k in (
            "removed_radial_quadrature_error", "removed_lapse_quadrature_error")}))


def physical_view(state, x, source, gauge, eta, alpha, budget):
    """Reconstruct saved physical values; no metric constraint or wall solve."""
    base = vacuum_slice(eta*float(gauge["W"]), x, budget=budget)
    wall = {k: float(gauge[k])+float(state["wall"][k]) for k in ("r", "L", "T", "s")}
    if wall["r"] <= 0:
        raise ValueError("saved physical wall radius must be positive")
    chi = float(np.cumsum(np.asarray(state["cap"]["Iq"]))[-1])
    wall["c"] = wall["r"]*_wall_conversion(chi/wall["r"], alpha)
    increment = {}
    for key in ("r", "L", "c"):
        src = float(source["log_C" if key == "L" else key])
        vacuum = base[key] if key != "c" else np.zeros_like(eta)
        offset = src-vacuum[0]
        increment[key] = vacuum[1:-1]+offset+np.asarray(state["departure"][key])-src
    return {"increment": increment, "wall": wall}


def _histories(parent, case, budget, stride=1):
    snapshots = case["snapshots"][::stride]
    indices = np.arange(0, 513, 4*stride)
    inputs = dict(eta=np.linspace(0., 1., case["cells"]+1),
        gauge=subset(parent["gauge"]["rows"], indices),
        source=subset(parent["initial_data"][case["amplitude"]]["snapshot_rows"], indices),
        alpha=case["alpha"], tube=TUBE, budget=budget, source_fronts=parent["geometry"]["support_x"])
    raw = coupled_history_diagnostics([dict(x=s["x"], state=s["physical"]) for s in snapshots], **inputs)
    departures = departure_history_diagnostics(snapshots, **inputs)
    return raw, departures


def _tables(case, history):
    return {"exterior": history["states"],
        "wall": {k: np.array([s["physical"]["wall"][k] for s in case["snapshots"]])
                 for k in ("r", "L", "c", "T", "s")},
        "cap": {k: np.array([s["state"]["cap"][k] for s in case["snapshots"]]) for k in ("Ip", "Iq")}}


def _old_tables(parent, case, budget):
    rows = parent["initial_data"][case["amplitude"]]["snapshot_rows"]
    arrays = {k: [] for k in ("r", "L", "c")}
    for i, sample in enumerate(case["snapshots"]):
        for key in arrays:
            src = rows["log_C" if key == "L" else key][4*i]
            arrays[key].append(np.r_[src, np.asarray(sample["state"]["increment"][key])+src,
                                     sample["state"]["wall"][key]])
    budget.check(3*129*(case["cells"]+1))
    return dict(exterior={k: np.array(v) for k, v in arrays.items()},
        wall={k: np.array([s["state"]["wall"][k] for s in case["snapshots"]]) for k in ("r", "L", "c", "T", "s")})


def _restricted(table):
    return dict(exterior={k: v[:, ::2] for k, v in table["exterior"].items()}, wall=table["wall"],
                cap={k: v[:, ::2]+v[:, 1::2] for k, v in table["cap"].items()})


def _cap_characteristics(case, budget):
    output = {key: [] for key in ("a_minus", "canonical_wall_chi", "canonical_center_f")}
    for sample in case["snapshots"]:
        n, R = case["cells"], sample["physical"]["wall"]["r"]
        budget.check(4*n+1)
        rec = _reconstruct(R*np.linspace(0., 1., n+1), sample["state"]["cap"])
        output["a_minus"].append(float(rec["minus"][-1]+rec["minus_slope"][-1]*(R-rec["midpoint"][-1])))
        output["canonical_wall_chi"].append(float(rec["chi"][-1]))
        output["canonical_center_f"].append(float(rec["Q0"]))
    return dict(arrays=output, statistics={k: statistics(v) for k, v in output.items()},
        interpretation="sampled stated reconstruction; no radial metric solve or propagation error enclosure")


def diagnostics(record, parent, budget, *, replay_allowances=None):
    histories, tables, controls, overlaps, old_differences, caps = [], [], [], [], [], []
    old_rectangle = json.loads((ROOT/previous.previous.RELATIVE_OUTPUT/"rectangle.json").read_bytes())
    gauge129 = subset(parent["gauge"]["rows"], slice(None, None, 4))
    for case, old in zip(record["cases"], parent["cases"]):
        full, departed = _histories(parent, case, budget)
        coarse, coarse_departed = _histories(parent, case, budget, 2)
        histories.append(dict(raw129=previous.compact_history(full), raw65=previous.compact_history(coarse),
            departure129=previous.compact_history(departed), departure65=previous.compact_history(coarse_departed)))
        collect_departure_replay(replay_allowances, len(histories)-1,
            parent["initial_data"][case["amplitude"]]["snapshot_rows"], case["alpha"],
            departure129=departed, departure65=coarse_departed)
        table = _tables(case, full)
        tables.append(table)
        old_table = _old_tables(parent, old, budget)
        old_differences.append(previous._difference_summary(
            {k: table[k] for k in ("exterior", "wall")}, old_table))
        caps.append(_cap_characteristics(case, budget))
        x = np.array([s["x"] for s in case["snapshots"]])
        quiet = x <= parent["geometry"]["support_x"][0]
        controls.append(dict(quiet_scalar_maximum=max(float(np.max(abs(table["exterior"]["c"][quiet]))),
            *[float(np.max(abs(value[quiet]))) for value in table["cap"].values()]),
            vacuum_departure_maximum={key: max(float(np.max(abs(np.asarray(s["state"]["departure"][key]))))
                for s in case["snapshots"]) for key in ("r", "L")},
            vacuum_wall_departure={key: statistics([s["state"]["wall"][key] for s in case["snapshots"]])
                for key in ("r", "L", "T", "s")},
            interpretation="vacuum preservation control only for alpha=0; finite rows retain physical response"))
        old_grid = old_rectangle["grids"][1 if case["amplitude"] == "zero" else 3]
        overlaps.append(compare_rectangle_overlap(table["exterior"], x,
            np.linspace(0., 1., case["cells"]+1), gauge129["W"], old_grid["x"], old_grid["Ubar"],
            old_grid["states"], budget=budget))
    paired = [previous._difference(tables[2], tables[0]), previous._difference(tables[3], tables[1])]
    return plain(dict(histories=histories, controls=controls, cap_characteristics=caps,
        successor_minus_v1=old_differences,
        original_rectangle_overlap=dict(reference=previous.previous.RELATIVE_OUTPUT+"/rectangle.json", cases=overlaps,
            interpretation="independent saved discretization and eta interpolation sensitivity; not an error enclosure"),
        source_interpolation_withheld_257_odd=copy.deepcopy(parent["diagnostics"]["source_interpolation_withheld_257_odd"]),
        spatial_zero=previous._difference_summary(_restricted(tables[1]), tables[0]),
        spatial_finite=previous._difference_summary(_restricted(tables[3]), tables[2]),
        finite_minus_zero=[previous._difference_summary(tables[j+2], tables[j]) for j in (0, 1)],
        paired_spatial=previous._difference_summary(_restricted(paired[1]), paired[0]),
        time_512_minus_256=previous._difference_summary(tables[4], tables[3]),
        cap_restriction=SETTINGS["cap_restriction"],
        comparison_coordinates="common x, exterior eta and cap y labels; physical radii can differ between cases",
        difference_convention="successor minus v1; fine minus coarse; finite minus separately evolved zero",
        constructed_identities_are_accuracy_tests=False,
        time_history_coarsening="129 versus65 saved slices; no additional evolution"))


def calculate(budget=None):
    budget = Budget(maximum_evaluations=8000000, cpu_seconds=180.) if budget is None else budget
    phase, record = {"name": "inherited_authentication"}, {}
    try:
        parent, identity, auth_charges = load_parent(budget)
        x = np.asarray(parent["gauge"]["rows"]["x"])
        record = dict(schema=SCHEMA, status=STATUS, scope=SCOPE, settings=SETTINGS,
            runtime=runtime_versions(), parent_input=identity, geometry=parent["geometry"],
            source_bounds=parent["source_bounds"], cases=[], refinement_is_rigorous_error_enclosure=False)
        for label, cells, steps in CASES:
            phase = dict(name="compatible_evolution", amplitude=label, cells=cells, steps=steps)
            source = StoredRows(parent["initial_data"][label]["snapshot_rows"])
            gauge = StoredRows(parent["gauge"]["rows"])
            start = budget.evaluations
            solution = evolve_compatible(quiet_departures(exterior_cells=cells, cap_cells=cells),
                x[::512//steps], eta=np.linspace(0., 1., cells+1), y_faces=np.linspace(0., 1., cells+1),
                source=source, gauge=gauge, alpha=parent["initial_data"][label]["alpha"], tube=TUBE,
                budget=budget, save_every=steps//128, maximum_courant=.4,
                exterior_quadrature_order=2, cap_quadrature_order=4)
            final = final_summary(solution["final_stage"], budget)
            record["cases"].append(plain(dict(amplitude=label, cells=cells, steps=steps,
                alpha=parent["initial_data"][label]["alpha"], snapshots=solution["snapshots"], final=final,
                execution=dict(budget_start=start, budget_end=budget.evaluations,
                    stages=solution["stages"], source_cache_reads=source.calls, gauge_cache_reads=gauge.calls,
                    max_courant=solution["max_courant"], retries=solution["retries"])) ))
        phase = {"name": "saved_array_diagnostics"}
        before = budget.evaluations
        record["diagnostics"] = diagnostics(record, parent, budget)
        budget.check()
        record["execution"] = dict(evaluations=budget.evaluations, cpu_seconds=time.process_time()-budget.started,
            inherited_authentication_evaluations=auth_charges, diagnostic_evaluations=budget.evaluations-before,
            source_integrations=0, evolution_source_snapshot_constructions=0,
            inherited_source_snapshot_reconstructions=4, gauge_quadratures=0,
            coupled_evolutions=5, time_steps=sum(k for _, _, k in CASES),
            stage_evaluations=sum(2*k+1 for _, _, k in CASES), inherited_authentications=1, old_rectangle_solves=0)
        return plain(record)
    except Exception as error:
        error.compatible_failure = plain(dict(phase=phase, exception=type(error).__name__, message=str(error),
            stage=getattr(error, "diagnostic", None), evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started, completed_cases=record.get("cases", []),
            parent_input=record.get("parent_input")))
        raise


def validate_record(record, budget=None, *, require_creation=True):
    budget = Budget(maximum_evaluations=8000000, cpu_seconds=180.) if budget is None else budget
    start = budget.evaluations
    try:
        if (record["schema"] != SCHEMA or record["status"] != STATUS or record["scope"] != SCOPE
                or record["refinement_is_rigorous_error_enclosure"] is not False):
            raise ValueError("compatible declaration mismatch")
        _equal(record["settings"], SETTINGS, "settings")
        if (set(record["runtime"]) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v for v in record["runtime"].values())):
            raise ValueError("runtime inventory mismatch")
        parent, identity, auth_charges = load_parent(budget)
        _equal(record["parent_input"], identity, "inherited input")
        _equal(record["geometry"], parent["geometry"], "inherited geometry")
        _equal(record["source_bounds"], parent["source_bounds"], "inherited source bounds")
        if len(record["cases"]) != len(CASES):
            raise ValueError("all five ordered cases required")
        x = np.asarray(parent["gauge"]["rows"]["x"])
        end = auth_charges
        for case, (label, n, steps) in zip(record["cases"], CASES):
            data = parent["initial_data"][label]
            if (case["amplitude"], case["cells"], case["steps"], case["alpha"]) != (label, n, steps, data["alpha"]):
                raise ValueError("case identity mismatch")
            if len(case["snapshots"]) != 129:
                raise ValueError("require129 saved time slices")
            eta = np.linspace(0., 1., n+1)
            for j, (sample, at) in enumerate(zip(case["snapshots"], x[::4])):
                if sample["x"] != at or set(sample) != {"x", "state", "physical"}:
                    raise ValueError("saved row/time inventory mismatch")
                state = sample["state"]
                if set(state) != {"departure", "wall", "cap"}:
                    raise ValueError("state inventory mismatch")
                for group, names, shape in (("departure", ("r", "L", "c"), (n-1,)),
                        ("wall", ("r", "L", "T", "s"), ()), ("cap", ("Ip", "Iq"), (n,))):
                    if set(state[group]) != set(names):
                        raise ValueError("state field inventory mismatch")
                    for key in names:
                        _array(state[group][key], shape, group+"."+key)
                expected = physical_view(state, at, subset(data["snapshot_rows"], 4*j),
                    subset(parent["gauge"]["rows"], 4*j), eta, data["alpha"], budget)
                _equal(sample["physical"], expected, "raw-to-physical reconstruction")
            _equal(case["snapshots"][0]["state"], quiet_departures(exterior_cells=n, cap_cells=n), "quiet initial departures")
            for key in ("T", "s"):
                clocks = np.array([s["physical"]["wall"][key] for s in case["snapshots"]])
                if clocks[0] != 0 or np.any(np.diff(clocks) <= 0):
                    raise ValueError("physical clocks must increase from quiet zero")
            stage = compatible_stage(case["snapshots"][-1]["state"], x=x[-1], eta=eta, y_faces=eta,
                source=subset(data["snapshot_rows"], -1), gauge=subset(parent["gauge"]["rows"], -1),
                alpha=data["alpha"], tube=TUBE, budget=budget,
                exterior_quadrature_order=2, cap_quadrature_order=4)
            _equal(case["final"], final_summary(stage, budget), "final radial/prefix reconstruction")
            execution, stages = case["execution"], 2*steps+1
            charge = stages*(36*n+11)+3*n
            expected_counts = dict(budget_start=end, budget_end=end+charge, stages=stages,
                source_cache_reads=stages, gauge_cache_reads=stages, retries=0)
            _equal({k: v for k, v in execution.items() if k != "max_courant"}, expected_counts, "case accounting")
            _real(execution["max_courant"], "Courant", high=.4)
            end += charge
        before = budget.evaluations
        replay_allowances = {}
        rebuilt = diagnostics(record, parent, budget, replay_allowances=replay_allowances)
        _equal(record["diagnostics"], rebuilt, "saved diagnostics", roundoff=replay_allowances)
        diagnostic_charges = budget.evaluations-before
        execution = record["execution"]
        expected = dict(evaluations=end+diagnostic_charges, inherited_authentication_evaluations=auth_charges,
            diagnostic_evaluations=diagnostic_charges, source_integrations=0,
            evolution_source_snapshot_constructions=0, inherited_source_snapshot_reconstructions=4,
            gauge_quadratures=0, coupled_evolutions=5, time_steps=sum(k for _, _, k in CASES),
            stage_evaluations=sum(2*k+1 for _, _, k in CASES), inherited_authentications=1, old_rectangle_solves=0)
        _equal({k: v for k, v in execution.items() if k != "cpu_seconds"}, expected, "execution accounting")
        _count(execution["evaluations"], 8000000, "calculation charges")
        _real(execution["cpu_seconds"], "calculation CPU", high=180.)
        budget.check()
        authentication_charges = budget.evaluations-start
        if require_creation:
            creation = record["creation_execution"]
            expected_creation = dict(calculation_evaluations=execution["evaluations"],
                validation_evaluations=authentication_charges,
                total_evaluations=execution["evaluations"]+authentication_charges)
            _equal({k: v for k, v in creation.items() if k != "cpu_seconds"}, expected_creation, "creation allowance")
            _count(creation["total_evaluations"], 8000000, "creation charges")
            _real(creation["cpu_seconds"], "creation CPU", low=execution["cpu_seconds"], high=180.)
    except (KeyError, TypeError, IndexError, OverflowError) as error:
        raise ValueError("incomplete compatible record") from error
    return dict(source_integrations=0, time_evolutions=0, reference_quadratures=0,
        inherited_authentications=1, inherited_final_cap_radial_reconstructions=5,
        final_cap_radial_reconstructions=5, final_prefix_wall_reconstructions=5,
        authentication_evaluations=budget.evaluations-start)


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: use a fresh directory or --check")
    bindings, runtime, started = source_bindings(), runtime_versions(), time.process_time()
    budget = Budget(maximum_evaluations=8000000, cpu_seconds=180., started=started)
    try:
        record = calculate(budget)
        if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
            raise RuntimeError("compatible sources or runtime changed during calculation")
        authentication = validate_record(record, budget, require_creation=False)
        if bindings != source_bindings() or runtime != runtime_versions():
            raise RuntimeError("compatible sources or runtime changed during authentication")
        record["creation_execution"] = dict(calculation_evaluations=record["execution"]["evaluations"],
            validation_evaluations=authentication["authentication_evaluations"], total_evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-started)
        budget.check()
    except Exception as error:
        try:
            after = source_bindings()
        except OSError as secondary:
            after = {"unavailable": str(secondary)}
        receipt = dict(schema=SCHEMA+"-failed-attempt", status="FAILED_NOT_COMPLETED", settings=SETTINGS,
            source_bindings_before=bindings, source_bindings_after=after, runtime=runtime,
            failure=getattr(error, "compatible_failure", dict(exception=type(error).__name__, message=str(error))))
        path = ROOT/".build/compatible-evolution-failures"/(uuid.uuid4().hex+".json")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(encoded(previous.previous._failure_safe(receipt)))
        error.failure_receipt = str(path)
        raise
    raw = encoded(record)
    manifest = dict(schema=SCHEMA+"-manifest", source_bindings=bindings,
                    outputs={"compatible.json": dict(sha256=digest(raw), bytes=len(raw))})
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"compatible.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return dict(created=str(directory), execution=record["execution"],
        creation_execution=record["creation_execution"], controls=record["diagnostics"]["controls"])


def check(directory):
    directory = Path(directory)
    if ({p.name for p in directory.iterdir()} != {"compatible.json", "manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("compatible output inventory mismatch")
    manifest = json.loads((directory/"manifest.json").read_bytes())
    if manifest.get("schema") != SCHEMA+"-manifest" or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES):
        raise ValueError("compatible source binding/manifest mismatch")
    raw = (directory/"compatible.json").read_bytes()
    if manifest.get("outputs") != {"compatible.json": dict(sha256=digest(raw), bytes=len(raw))}:
        raise ValueError("compatible payload hash/size mismatch")
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
