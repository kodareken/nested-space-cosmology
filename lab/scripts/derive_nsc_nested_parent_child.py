#!/usr/bin/env python3
"""One bounded live nested-pair campaign, with physical two-way controls.

All cases start from their own source-consistent initial data. No source or
geometry is reset during a trajectory. Only this driver launches the study.
Published inputs are preserved; successor outputs are never overwritten.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
from scipy.integrate import simpson
from recursive_horizons import nsc_nested_parent_child as model
from recursive_horizons import nsc_nested_parent_child_response as response
from recursive_horizons import nsc_regional_energy_exchange as regional
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-nested-parent-child-v1.json"
PAYLOAD = OUT.with_suffix(".npz")
CPU_LIMIT = 600.0
PAYLOAD_LIMIT = 64 * 1024**2
TARGET = .3
FRAME_COUNT = 31
COARSE_NF = 128
FINE_NF = 256
WEIGHTS = {
    "baseline": (.75, .75, .5, .5, .25, .25),
    "parent_only": (.85, .85, .5, .5, .15, .15),
    "child_only": (.75, .75, .6, .6, .25, .25),
}
DEPENDENCIES = (
    Path(__file__),
    LAB / "src/recursive_horizons/nsc_nested_parent_child.py",
    LAB / "src/recursive_horizons/nsc_nested_parent_child_response.py",
    LAB / "src/recursive_horizons/nsc_spherical_coupling.py",
    LAB / "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    LAB / "src/recursive_horizons/nsc_spherical_feedback_action.py",
    LAB / "src/recursive_horizons/nsc_conformal_adm_source.py",
    LAB / "src/recursive_horizons/nsc_regional_energy_exchange.py",
    LAB / "src/recursive_horizons/nsc_evolving_reduction.py",
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def jsonable(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imag": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(v) for v in value]
    return value


def state_columns(state):
    return np.vstack((state.phi0, state.phi1))


def state_hash(state):
    hasher = hashlib.sha256()
    for name in model.STATE_NAMES:
        array = np.ascontiguousarray(getattr(state, name))
        hasher.update(name.encode())
        hasher.update(str(array.dtype).encode())
        hasher.update(str(array.shape).encode())
        hasher.update(array.tobytes())
    return hasher.hexdigest()


def lifted_rate(pair, rate, bundle):
    nodal = bundle["nodal_rate"]
    fields = [galerkin.prolong_geometry(pair.grid, getattr(nodal, name))
              for name in model.GEOMETRY_NAMES + model.MOMENTUM_NAMES]
    fields += [pair.grid.U_f @ nodal.phi0, pair.grid.U_f @ nodal.phi1]
    return coupling.CauchyRate(*fields, rate.force_L, rate.force_Q,
                              rate.force_beta, rate.fieldwork_power)


def normal_energy_slope(system, state, rate, source):
    """Polarize the actual column source, then differentiate FL/r."""
    p0, p1 = system.momentum @ state.phi0, system.momentum @ state.phi1
    dp0, dp1 = system.momentum @ rate.phi0, system.momentum @ rate.phi1
    weights = system.occupations[None, :]
    kdot = np.sum(weights * (
        rate.phi0.conj() * (-1j * p1) + state.phi0.conj() * (-1j * dp1)
        + rate.phi1.conj() * (1j * p0) + state.phi1.conj() * (1j * dp0)), axis=1).real
    sdot = np.sum(weights * (
        rate.phi0.conj() * state.phi1 + state.phi0.conj() * rate.phi1
        + rate.phi1.conj() * state.phi0 + state.phi1.conj() * rate.phi0), axis=1).real
    fldot = system.multiplicity * (kdot / state.Q - source["K"] * rate.Q / state.Q**2
                                   + system.kappa * sdot)
    return fldot / state.r - source["force_L"] * rate.r / state.r**2


def observe(pair, state, mark, hierarchy):
    rate, bundle = model.rates(pair, state, return_bundle=True)
    fine, system = bundle["fine_state"], bundle["fine_system"]
    lifted = lifted_rate(pair, rate, bundle)
    ledger = regional.matter_ledger(system, fine)
    terms = regional.proper_balance_terms(system, fine, lifted, ledger)
    slope = normal_energy_slope(system, fine, lifted, bundle["source"])
    defect = slope + system.derivative @ terms["flux_nodal"] - terms["source_nodal"]
    probability = np.sum((abs(fine.phi0)**2 + abs(fine.phi1)**2)
                         * pair.weights[None, :], axis=1)
    intervals = {"child": pair.child_interval, "parent": pair.parent_interval,
                 "left_annulus": (0., 1.), "right_annulus": (3., 4.)}
    windows = {}
    for name, interval in intervals.items():
        flux_edges = model._periodic_values(pair.grid, terms["flux_nodal"], interval)
        integrate = lambda values: model.interval_integral(pair.grid, values / system.dx, interval)
        windows[name] = {
            "interval": list(interval),
            "normal_shell": integrate(ledger["normal_energy_nodal"]),
            "boundary_inflow": float((flux_edges[0] - flux_edges[1]) / system.dx),
            "pressure_work": integrate(terms["proper_pressure_work"]),
            "lapse_gradient_work": integrate(terms["momentum_lapse_work"]),
            "balance_defect": integrate(defect),
            "probability": integrate(probability),
        }
    h = model.hamiltonian(pair, state)
    accounting = response.block_accounting(h, state_columns(state), pair.weights,
                                          hierarchy, multiplicity=system.multiplicity)
    blocks = response.nested_blocks(h, hierarchy)
    constraints = galerkin.conformal_constraint_transport(
        pair.grid, bundle["nodal_state"], bundle["nodal_rate"], bundle)
    metrics = model.metrics(pair, state, rate)
    column_gram = state_columns(state).conj().T @ state_columns(state)
    roots = np.sqrt(pair.weights)
    eigenvalues = np.linalg.eigvalsh(roots[:, None] * column_gram * roots[None, :])
    coefficients = {name: {
        "parent_norm": float(np.linalg.norm(getattr(rate, name)[pair.geometry_parent_indices])),
        "child_norm": float(np.linalg.norm(getattr(rate, name)[pair.geometry_child_indices])),
    } for name in model.GEOMETRY_NAMES + model.MOMENTUM_NAMES}
    return {"time": float(mark), "metrics": metrics, "windows": windows,
            "energy": model.energy(pair, state),
            "fieldwork_power": rate.fieldwork_power,
            "field_energy": accounting["source_energy"],
            "field_cross_energy": float(np.sum(accounting["interaction_energy_upper"])),
            "field_energy_closure": accounting["energy_closure_gap"],
            "child_occupation": np.diag(accounting["child_covariance"]).real,
            "parent_detail_occupation": np.diag(accounting["detail_covariance"]).real,
            "link_frobenius": float(np.linalg.norm(blocks["child_detail"])),
            "gram_gap": float(np.max(abs(column_gram - np.eye(6)))),
            "covariance_eigenvalue_min": float(eigenvalues[0]),
            "covariance_eigenvalue_max": float(eigenvalues[-1]),
            "canonical_rates": coefficients, "constraint_transport": constraints}


def integral(times, values):
    return {"simpson": float(simpson(values, x=times)),
            "trapezoid": float(np.trapezoid(values, x=times)),
            "indicator": float(abs(simpson(values, x=times) - np.trapezoid(values, x=times)))}


def summary(pair, initial, final, rows):
    times = np.array([row["time"] for row in rows])
    windows = {}
    for name in rows[0]["windows"]:
        change = rows[-1]["windows"][name]["normal_shell"] - rows[0]["windows"][name]["normal_shell"]
        terms = {key: integral(times, [row["windows"][name][key] for row in rows])
                 for key in ("boundary_inflow", "pressure_work", "lapse_gradient_work", "balance_defect")}
        balance = sum(terms[key]["simpson"] for key in
                      ("boundary_inflow", "pressure_work", "lapse_gradient_work"))
        windows[name] = {"shell_change": change, "terms": terms,
                         "accounted_change": balance, "closure_gap": change - balance,
                         "closure_after_measured_projection_defect": change - balance - terms["balance_defect"]["simpson"]}
    clocks = np.array([row["metrics"]["clock_rates"] for row in rows])
    clock_integrals = np.vstack((np.zeros(3), np.cumsum(.5 * np.diff(times)[:, None]
                               * (clocks[:-1] + clocks[1:]), axis=0)))
    forcing = integral(times, [row["constraint_transport"]["forcing_norm"] for row in rows])
    constraint_budget = rows[0]["constraint_transport"]["constraint_norm"] + forcing["trapezoid"]
    return {"target_reached": abs(times[-1] - TARGET) < 1e-12,
            "initial_state_sha256": state_hash(initial), "final_state_sha256": state_hash(final),
            "initial": rows[0], "final": rows[-1], "windows": windows,
            "clock_protocol": "continuous normal clocks at x=1,2,3; common coordinate T comparisons",
            "normal_clocks": clock_integrals.tolist(),
            "total_energy_change": rows[-1]["energy"] - rows[0]["energy"],
            "field_energy_change": rows[-1]["field_energy"] - rows[0]["field_energy"],
            "coordinate_metric_work": integral(times, [row["fieldwork_power"] for row in rows]),
            "minimum_r": min(row["metrics"]["r_min"] for row in rows),
            "minimum_Q": min(row["metrics"]["Q_min"] for row in rows),
            "maximum_gram_gap": max(row["gram_gap"] for row in rows),
            "minimum_covariance_eigenvalue": min(row["covariance_eigenvalue_min"] for row in rows),
            "maximum_covariance_eigenvalue": max(row["covariance_eigenvalue_max"] for row in rows),
            "maximum_field_energy_closure": max(abs(row["field_energy_closure"]) for row in rows),
            "sampled_constraint_forcing_integral": forcing,
            "sampled_constraint_end_budget": constraint_budget,
            "sampled_constraint_end_margin": constraint_budget - rows[-1]["constraint_transport"]["constraint_norm"],
            "constraint_time_integral_certified": False,
            "canonical_groups_evolve": {
                group: max(row["canonical_rates"][name][group + "_norm"]
                           for row in rows for name in model.STATE_NAMES[:6])
                for group in ("parent", "child")}}


def compare_cases(records):
    assessments = {}
    observables = {
        "parent_to_child_state": ("parent_only", "child_probability"),
        "child_to_parent_state": ("child_only", "parent_annulus_probability"),
        "parent_to_child_modal_state": ("parent_only", "child_occupation"),
        "child_to_parent_detail_modal_state": ("child_only", "parent_detail_occupation"),
        "parent_to_child_geometry": ("parent_only", "child_r_proper_mean"),
        "child_to_parent_annulus_geometry": ("child_only", "parent_annulus_r_proper_mean"),
    }
    def curve(nf, cap, control, observable):
        rows = records[f"nf{nf}_{control}_dt{cap:g}"]["rows"]
        if observable == "child_probability":
            values = np.array([row["windows"]["child"]["probability"] for row in rows])
        elif observable == "parent_annulus_probability":
            values = np.array([row["windows"]["left_annulus"]["probability"]
                               + row["windows"]["right_annulus"]["probability"] for row in rows])
        elif observable.endswith("occupation"):
            values = np.array([row[observable] for row in rows])
        else:
            values = np.array([row["metrics"][observable] for row in rows])
        return values - values[0]
    for name, (control, observable) in observables.items():
        effects = {(nf, cap): curve(nf, cap, control, observable) - curve(nf, cap, "baseline", observable)
                   for nf in (COARSE_NF, FINE_NF) for cap in (.001, .0005)}
        effect = float(np.max(abs(effects[(FINE_NF, .0005)])))
        time_gap = float(np.max(abs(effects[(FINE_NF, .001)] - effects[(FINE_NF, .0005)])))
        space_gap = float(np.max(abs(effects[(COARSE_NF, .0005)] - effects[(FINE_NF, .0005)])))
        assessments[name] = {"control": control, "observable": observable,
            "effect_max": effect, "time_indicator": time_gap, "space_indicator": space_gap,
            "time_fraction": time_gap / effect if effect else None,
            "space_fraction": space_gap / effect if effect else None,
            "resolved_at_one_percent": bool(effect > 0 and max(time_gap, space_gap) < .01 * effect),
            "initial_offset_removed": True, "matched_proper_time_claimed": False}
    return assessments


def response_on_trajectory(pair, states, times, budget):
    start = len(times) // 3
    selected_times = np.asarray(times[start:], dtype=float)
    selected = states[start:]
    @lru_cache(maxsize=4)
    def actual_hamiltonian(mark):
        upper = min(max(np.searchsorted(times, mark, side="right"), 1), len(times) - 1)
        lower = upper - 1
        fraction = (mark - times[lower]) / (times[upper] - times[lower])
        state = model.NestedState(*(getattr(states[lower], name) * (1 - fraction)
                                  + getattr(states[upper], name) * fraction for name in model.STATE_NAMES))
        return model.hamiltonian(pair, state)
    hierarchy = response.nested_frame(pair.reference_columns)
    return response.retained_response(actual_hamiltonian, selected_times,
        state_columns(selected[0]), pair.weights, hierarchy,
        autonomous_columns=np.array([state_columns(state) for state in selected]),
        cpu_limit_s=budget)


def save(record, arrays):
    if OUT.exists() or PAYLOAD.exists():
        raise FileExistsError("successor evidence already exists; use --check")
    temporary = PAYLOAD.with_suffix(".npz.tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    if temporary.stat().st_size > PAYLOAD_LIMIT:
        temporary.unlink()
        raise RuntimeError("successor payload exceeds 64 MiB")
    os.replace(temporary, PAYLOAD)
    record["payload_sha256"] = digest(PAYLOAD)
    record["payload_bytes"] = PAYLOAD.stat().st_size
    record["cpu_seconds"] = time.process_time() - record.pop("_started_cpu")
    OUT.write_text(json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n")


def run():
    if OUT.exists() or PAYLOAD.exists():
        raise FileExistsError("refusing to overwrite existing scientific evidence")
    started = time.process_time()
    bindings = {path.relative_to(LAB).as_posix(): digest(path) for path in DEPENDENCIES}
    record = {"schema": "NSC-NESTED-PARENT-CHILD-v1", "status": "INCOMPLETE",
        "_started_cpu": started, "target": TARGET, "cpu_budget_seconds": CPU_LIMIT,
        "payload_cap_bytes": PAYLOAD_LIMIT, "source_bindings": bindings,
        "one_global_action": True, "parent_interval": [0., 4.], "child_interval": [1., 3.],
        "coordinate_scale_ratio": 2., "physical_scale_ratio_imposed": False,
        "source_layout": "separated", "weights": WEIGHTS, "results": {}, "forecasts": {},
        "comparison_clock": "same global coordinate T; each proper clock recorded separately",
        "source_changed_during_evolution": False, "geometry_prescribed": False,
        "historical_paper_rewritten": False, "tests_are_not_trajectory_evidence": True}
    arrays, baseline = {}, None
    try:
        for nf in (COARSE_NF, FINE_NF):
            for control, weights in WEIGHTS.items():
                pair = model.build_pair(nf, occupations=weights)
                initial, preparation = model.initial_state(pair)
                hierarchy = response.nested_frame(pair.reference_columns)
                for cap in (.001, .0005):
                    name = f"nf{nf}_{control}_dt{cap:g}"
                    dt, omega, _ = model.stable_timestep(pair, initial, cap)
                    pilot_cpu = time.process_time()
                    model.rk4_step(pair, initial, dt)
                    step_cost = time.process_time() - pilot_cpu
                    frame_cpu = time.process_time()
                    first = observe(pair, initial, 0., hierarchy)
                    frame_cost = time.process_time() - frame_cpu
                    forecast = 1.5 * (np.ceil(TARGET / dt) * step_cost + FRAME_COUNT * frame_cost)
                    record["forecasts"][name] = {"step_cpu_seconds": step_cost,
                        "frame_cpu_seconds": frame_cost, "omega": omega,
                        "forecast_cpu_seconds": float(forecast)}
                    if time.process_time() - started + forecast + 30 > CPU_LIMIT:
                        raise RuntimeError("next case does not fit declared shared CPU budget")
                    print("admitted", name, "CPU forecast", round(forecast, 3), flush=True)
                    state = initial.copy()
                    times = np.linspace(0., TARGET, FRAME_COUNT)
                    states, rows = [initial.copy()], [first]
                    mark = 0.
                    for frame_time in times[1:]:
                        while mark < frame_time - 1e-13:
                            step, _, _ = model.stable_timestep(pair, state, cap)
                            step = min(step, frame_time - mark)
                            state = model.rk4_step(pair, state, step)
                            mark += step
                            if time.process_time() - started > CPU_LIMIT - 10:
                                raise RuntimeError("shared CPU budget reached; checkpoint preserved")
                        states.append(state.copy())
                        rows.append(observe(pair, state, mark, hierarchy))
                        if rows[-1]["gram_gap"] > 1e-8:
                            raise RuntimeError("column Gram lost admissibility at the declared threshold")
                    data = {"nf": nf, "step_cap": cap, "control": control,
                        "preparation": preparation, "source": pair.source_metadata,
                        "geometry_basis": pair.geometry_metadata, "rows": rows,
                        "summary": summary(pair, initial, state, rows)}
                    record["results"][name] = data
                    for field in model.STATE_NAMES:
                        arrays[name + "_" + field] = np.array([getattr(s, field) for s in states])
                    arrays[name + "_times"] = times
                    arrays[name + "_reference_columns"] = pair.reference_columns
                    if nf == FINE_NF and control == "baseline" and cap == .0005:
                        baseline = pair, states, times
                    print("completed", name, "energy drift", data["summary"]["total_energy_change"], flush=True)
        record["two_way_effects"] = compare_cases(record["results"])
        pair, states, times = baseline
        reduction_budget = min(90., CPU_LIMIT - (time.process_time() - started) - 5)
        if reduction_budget <= 0:
            raise RuntimeError("no admitted CPU budget remains for local response")
        result = response_on_trajectory(pair, states, times, reduction_budget)
        record["local_response"] = result
        all_targets = all(case["summary"]["target_reached"] for case in record["results"].values())
        record["status"] = "MEASURED_NESTED_PAIR" if all_targets else "INCOMPLETE"
        primary_effects = ("parent_to_child_state", "child_to_parent_state",
                          "parent_to_child_geometry", "child_to_parent_annulus_geometry")
        record["primary_two_way_effects"] = list(primary_effects)
        record["all_two_way_effects_resolved"] = all(record["two_way_effects"][name]["resolved_at_one_percent"]
                                                  for name in primary_effects)
        record["nested_schur"] = response.nested_schur(
            model.hamiltonian(pair, states[-1]), response.nested_frame(pair.reference_columns))
    except (RuntimeError, ValueError, coupling.PositiveChartExit) as error:
        record["stop_reason"] = str(error)
        if "state" in locals():
            for field in model.STATE_NAMES:
                arrays["interrupted_checkpoint_" + field] = np.array(getattr(state, field), copy=True)
            record["interrupted_time"] = mark if "mark" in locals() else 0.
    if bindings != {path.relative_to(LAB).as_posix(): digest(path) for path in DEPENDENCIES}:
        raise RuntimeError("scientific source changed during campaign")
    save(record, arrays)
    return record


def verify():
    record = json.loads(OUT.read_text())
    if record.get("schema") != "NSC-NESTED-PARENT-CHILD-v1" or digest(PAYLOAD) != record["payload_sha256"]:
        raise ValueError("invalid successor schema or payload")
    for name, expected in record["source_bindings"].items():
        if digest(LAB / name) != expected:
            raise ValueError("source binding changed: " + name)
    with np.load(PAYLOAD, allow_pickle=False) as stored:
        for name, case in record["results"].items():
            states = [model.NestedState(*(stored[name + "_" + field][index] for field in model.STATE_NAMES))
                      for index in (0, -1)]
            if state_hash(states[0]) != case["summary"]["initial_state_sha256"] or state_hash(states[1]) != case["summary"]["final_state_sha256"]:
                raise ValueError("saved Cauchy state mismatch: " + name)
    if "two_way_effects" in record and compare_cases(record["results"]) != record["two_way_effects"]:
        raise ValueError("two-way effects do not replay")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = verify() if args.check else run()
    print(result["status"], "CPU", result["cpu_seconds"],
          "resolved two-way effects", result.get("all_two_way_effects_resolved"), flush=True)
