#!/usr/bin/env python3
"""Bounded source-weight controls in the same conformal spherical action.

The six saved T0 columns and total occupation3 stay fixed. Each α control
recomputes the source and solves the existing initial radius/shift ansatz.
Four nf256, dtcap.0005 trajectories run from0 toward.3 under600CPU total.
Their finite/coarse error domain is explicit; uniform is not assumed to
annihilate every current of a covariance with an empty complement.
"""
from __future__ import annotations
import json
import math
import os
from pathlib import Path
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_frames as frames
import derive_nsc_spherical_conformal_continuation as continuation
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import nsc_regeneration_controls as controls
from recursive_horizons.nsc_spherical_coupling import CauchyState, PositiveChartExit

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-spherical-conformal-source-controls-v1.json"
NPZ = LAB / "results/development/nsc-spherical-conformal-source-controls-v1.npz"
DRIVER = Path(__file__).resolve()
TARGET = .3
FRAME = .005
DT_CAP = .0005
CPU_BUDGET = 600.
SAFETY = 1.25
RUN_PLAN = (("uniform_alpha0", 0.), ("reversed_alpha_minus1", -1.),
            ("imbalance_alpha_0_95", .95), ("imbalance_alpha_1_05", 1.05))
SERIES_KEYS = ("time", "field_energy", "gravity_energy", "energy", "fieldwork_power", "constraint_norm", "forcing_norm",
               "quadrature_forcing_norm", "projection_forcing_norm", "full_C_l2", "full_D_l2", "full_hamilton_max", "full_momentum_max",
               "proper_min", "proper_max", "r_min", "r_max", "Q_min", "Q_max", "gram_gap", "omega", "clock_rate", "proper_clock")
WINDOW_KEYS = ("normal_shell", "coordinate_matter_energy", "coordinate_gravity_energy", "normal_boundary_flux", "normal_left_outward_flux",
               "normal_right_outward_flux", "proper_pressure_work", "momentum_lapse_work", "coordinate_metric_work", "normal_balance_defect",
               "probability", "probability_share", "probability_boundary_flux", "probability_balance_defect", "normal_internal_flux_l2")


def source_weights(alpha):
    alpha = float(alpha)
    weights = np.array((.5 + .25 * alpha, .5 + .25 * alpha, .5, .5, .5 - .25 * alpha, .5 - .25 * alpha))
    if not np.isfinite(weights).all() or np.min(weights) < 0 or np.max(weights) > 1:
        raise ValueError("source weights outside CAR interval")
    if abs(np.sum(weights) - 3) > 1e-14:
        raise ValueError("total occupation changed")
    return weights


def strip_sample(row):
    return {key: row[key] for key in SERIES_KEYS if key in row} | {"occupation": row["occupation"]}


def signed_surface_summary(rows):
    times = [row["time"] for row in rows]
    result = []
    for index in range(4):
        per_surface = {}
        for channel in ("normal_left_outward_flux", "normal_right_outward_flux", "normal_boundary_flux"):
            values = np.array([row["windows"][index][channel] for row in rows])
            intervals = np.diff(times)
            per_surface[channel] = {
                "signed_integral": float(np.sum(intervals * .5 * (values[:-1] + values[1:]))),
                "absolute_integral": float(np.sum(intervals * .5 * (abs(values[:-1]) + abs(values[1:])))),
                "endpoint": float(values[-1]),
                "quadrature": episode.integrate(values, FRAME),
            }
        result.append({"interval": list(frames.BOUNDS[index:index + 2]), **per_surface})
    return result


def summarize(samples, physical):
    initial, final = samples[0], samples[-1]
    integrals = {name: continuation.scalar_integral(samples, name) for name in
                 ("forcing_norm", "quadrature_forcing_norm", "projection_forcing_norm", "fieldwork_power", "clock_rate")}
    r0 = initial["constraint_norm"]
    end_budget = r0 + integrals["forcing_norm"]["trapezoid"]
    physical_summary = frames.summarize(physical) if len(physical) > 1 else None
    return {"initial": initial, "final": final, "time": final["time"], "integrals": integrals,
            "field_energy_change": final["field_energy"] - initial["field_energy"],
            "gravity_energy_change": final["gravity_energy"] - initial["gravity_energy"],
            "energy_balance_change": final["energy"] - initial["energy"],
            "fieldwork_difference": final["field_energy"] - initial["field_energy"] - integrals["fieldwork_power"]["trapezoid"],
            "constraint_end_sampled_budget": end_budget, "constraint_end_margin": end_budget - final["constraint_norm"],
            "positive_chart": all(row["r_min"] > 0 and row["Q_min"] > 0 for row in samples),
            "gram_gap_max": max(row["gram_gap"] for row in samples),
            "occupation_change": np.asarray(final["occupation"]) - np.asarray(initial["occupation"]),
            "physical": physical_summary, "surfaces": signed_surface_summary(physical),
            "continuum_constraint_certified": False, "observable_error_certified": False, "renewal": False}


def run():
    start = time.process_time()
    paths = {"v1_source_json": episode.OUT, "v1_source_npz": episode.NPZ, "v2_source_json": continuation.OUT,
             "v2_source_npz": continuation.NPZ, "driver": DRIVER, "galerkin": galerkin._MODULE_PATH,
             "source_constructor": Path(controls.__file__), "coupling": galerkin._COUPLING_PATH,
             "action": galerkin._ACTION_PATH, "conformal_source": galerkin._SOURCE_PATH}
    before = {name: episode.sha256(path) for name, path in paths.items()}
    grid = galerkin.build_grid(256, quadrature=1024, gauge="conformal")
    with np.load(episode.NPZ, allow_pickle=False) as v1:
        original = CauchyState(**{name: v1["nf256_initial_" + name].copy() for name in episode.STATE_NAMES})
    observer = np.vstack((original.phi0, original.phi1))[:, :2].copy()
    arrays = {"original_T0_phi0": original.phi0, "original_T0_phi1": original.phi1, "original_observer": observer}
    record = {"schema": "NSC-SPHERICAL-CONFORMAL-SOURCE-CONTROLS-v1", "gauge": "conformal", "nf": 256, "nq": 1024,
              "target": TARGET, "dt_cap": DT_CAP, "cpu_budget_seconds": CPU_BUDGET, "forecast_safety": SAFETY,
              "run_plan": RUN_PLAN, "source_parameter": "(.5+.25alpha,.5+.25alpha,.5,.5,.5-.25alpha,.5-.25alpha)",
              "source_inventory_changed": False, "prepared_columns_changed": False, "total_occupation": 3.,
              "mean_current_deleted": False, "midtrajectory_reset": False, "radius_imposed": False,
              "constructor": "nsc_regeneration_controls.prepare_population: source radius Newton plus projected periodic shift antiderivative",
              "admission": "positive/source-derived finite chart, antiderivative accepted, unchanged columns; historical initial tolerance flags are reported, not retuned",
              "uniform_annihilation_assumed": False, "baseline_alpha1_sealed": True,
              "source_bindings_before": before, "results": {}, "series_schema": SERIES_KEYS, "window_schema": WINDOW_KEYS,
              "error_domain": "exploratory nf256 controls; each initial full/held-out residual and actual forcing budget reported. Baseline nf256/nf512 and dt indicators are reference scales, not a control-specific continuum certificate"}
    prepared = []
    for label, alpha in RUN_PLAN:
        owned, state, constructor = controls.prepare_population(grid, original, source_weights(alpha))
        admitted = constructor["chart_admissible"] and constructor["source_derived"] and constructor["columns_unchanged"] and constructor["momentum"]["antiderivative_accepted"]
        if not admitted:
            raise RuntimeError(f"source preparation failed for {label}: {constructor.get('blocker')}")
        sample = strip_sample(episode.observe(owned, state, state, observer, 0.))
        sample["proper_clock"] = 0.
        prepared.append((label, alpha, owned, state, constructor, sample))
        for name in episode.STATE_NAMES:
            arrays[label + "_initial_" + name] = getattr(state, name)
        arrays[label + "_weights"] = owned.fine.occupations
        print("prepared", label, "fullC", constructor["full_hamilton_max"], "meanj", constructor["current_mean"], "CPU", time.process_time() - start, flush=True)
    pilot_start = time.process_time()
    _, _, owned, state, _, _ = prepared[0]
    trial = galerkin.rk4_step(owned, state, DT_CAP)
    episode.observe(owned, trial, state, observer, DT_CAP)
    step_cost = time.process_time() - pilot_start
    frame_start = time.process_time()
    frames.frame_ledger(owned, state, 0.)
    frame_cost = time.process_time() - frame_start
    remaining = CPU_BUDGET - (time.process_time() - start) - 5.
    forecast = len(prepared) * (math.ceil(TARGET / DT_CAP) * step_cost + math.ceil(TARGET / FRAME) * frame_cost)
    admitted_target = SAFETY * forecast <= remaining
    duration = TARGET if admitted_target else max(0., math.floor(TARGET * remaining / (SAFETY * forecast) / FRAME) * FRAME)
    record["forecast"] = {"measured_full_step_cpu": step_cost, "measured_frame_cpu": frame_cost,
                          "forecast_four_cases_cpu": forecast, "full_target_admitted": admitted_target,
                          "selected_common_endpoint": duration, "preparation_cpu": time.process_time() - start}
    print("forecast", record["forecast"], flush=True)
    if duration < FRAME:
        raise RuntimeError("measured preparation/step cost leaves no discriminating frame within600CPU")
    for label, alpha, owned, initial, constructor, first in prepared:
        current = initial.copy()
        samples, physical = [first], []
        physical.append(frames.frame_ledger(owned, current, 0.)[0])
        time_value, next_frame = 0., FRAME
        stop_reason = None
        run_start = time.process_time()
        cost = step_cost
        while time_value < duration - 1e-12:
            if time.process_time() - start + max(1., SAFETY * cost) + 5. >= CPU_BUDGET:
                stop_reason = "EXPLORATORY_CPU_BUDGET"
                break
            omega = galerkin.subspace_frequency(owned, current)[0]
            dt = min(DT_CAP, 1.4 / max(omega, 1.), next_frame - time_value, duration - time_value)
            if dt <= 1e-12:
                next_frame += FRAME
                continue
            step_start = time.process_time()
            try:
                candidate = galerkin.rk4_step(owned, current, dt)
                following = strip_sample(episode.observe(owned, candidate, initial, observer, time_value + dt))
                if following["gram_gap"] > 1e-8 or dt * following["omega"] > 2 * np.sqrt(2):
                    stop_reason = "NUMERICAL_ADMISSIBILITY"
                    break
                following["proper_clock"] = samples[-1]["proper_clock"] + .5 * dt * (samples[-1]["clock_rate"] + following["clock_rate"])
                current = candidate
                time_value += dt
                samples.append(following)
                if abs(time_value - next_frame) < 1e-10 or abs(time_value - duration) < 1e-10:
                    physical.append(frames.frame_ledger(owned, current, time_value)[0])
                    next_frame += FRAME
            except PositiveChartExit as error:
                stop_reason = error.reason
                break
            cost = max(cost, time.process_time() - step_start)
            if len(samples) % 200 == 0:
                print("progress", label, time_value, "CPU", time.process_time() - start, flush=True)
        if abs(physical[-1]["time"] - time_value) > 1e-12:
            physical.append(frames.frame_ledger(owned, current, time_value)[0])
        summary = summarize(samples, physical)
        summary.update(alpha=alpha, weights=owned.fine.occupations, constructor=constructor,
                       stop_reason=stop_reason, completed=stop_reason is None and abs(time_value - duration) < 1e-12,
                       cpu_seconds=time.process_time() - run_start)
        record["results"][label] = summary
        arrays[label + "_series"] = np.array([[row[key] for key in SERIES_KEYS] for row in samples])
        arrays[label + "_occupation"] = np.array([row["occupation"] for row in samples])
        arrays[label + "_physical_times"] = np.array([row["time"] for row in physical])
        arrays[label + "_window_series"] = np.array([[[window[key] for key in WINDOW_KEYS] for window in row["windows"]] for row in physical])
        for name in episode.STATE_NAMES:
            arrays[label + "_final_" + name] = getattr(current, name)
        print("complete", label, "T", time_value, "CPU", summary["cpu_seconds"], "shellleader", summary["physical"]["final"]["shell_leader"], "flowx0", summary["surfaces"][0]["normal_left_outward_flux"]["signed_integral"], flush=True)
    record["cpu_seconds"] = time.process_time() - start
    record["source_bindings_after"] = {name: episode.sha256(path) for name, path in paths.items()}
    record["sealed_sources_unchanged"] = record["source_bindings_after"] == before
    record["all_target_reached"] = duration == TARGET and all(row["completed"] for row in record["results"].values())
    record["verdict"] = "MEASURED_COARSE_SOURCE_CONTROLS" if record["all_target_reached"] else "PARTIAL_SOURCE_CONTROL_CHECKPOINTS"
    record["continuum_constraint_certified"] = False
    record["observable_error_certified"] = False
    with NPZ.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    record["payload_sha256"] = episode.sha256(NPZ)
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    record["payload_bytes"] = OUT.stat().st_size + NPZ.stat().st_size
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    print("saved", OUT, "CPU", record["cpu_seconds"], "payload", record["payload_bytes"], flush=True)
    return record


if __name__ == "__main__":
    run()
