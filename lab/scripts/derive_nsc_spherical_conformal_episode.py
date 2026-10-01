#!/usr/bin/env python3
"""Bounded same-action conformal episode from the unchanged saved T=0 data.

The legacy episode/replay owners and payloads are never written. Four runs
compare nf256/nf512, dt .001/.0005. A measured full-step cost selects a
common .005 frame endpoint within 300 CPU seconds. The original fixed
region-0 pair is read directly from the saved columns. A frozen-geometry
control uses the exact constant-Q conformal Dirac group, not a new source.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np

import derive_nsc_spherical_feedback_episode as legacy
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_spherical_coupling import CauchyState, PositiveChartExit, field_energy, gravity_energy

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-spherical-conformal-episode-v1.json"
NPZ = LAB / "results/development/nsc-spherical-conformal-episode-v1.npz"
NOTE = LAB / "docs/nsc-spherical-conformal-gauge.md"
DRIVER = Path(__file__).resolve()
SCHEMA = "NSC-SPHERICAL-CONFORMAL-EPISODE-v1"
TARGET = .05
FRAME = .005
CPU_BUDGET = 300.
PAYLOAD_LIMIT = 64 * 1024 * 1024
FORECAST_SAFETY = 1.25
RUN_PLAN = tuple({"name": f"nf{nf}_dt_{dt:.4f}", "nf": nf, "dt": dt, "gauge": "conformal"}
                 for nf in (256, 512) for dt in (.001, .0005))
STATE_NAMES = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frozen_columns(state, length, kappa, time_value):
    """Exact exp(-it H0) for the saved constant-Q initial conformal metric."""
    if np.ptp(state.Q) > 1e-12:
        raise ValueError("frozen control requires the saved constant initial Q")
    count = state.phi0.shape[0]
    coordinate = np.arange(count) * length / count
    phase = np.exp(1j * np.pi * coordinate / length)[:, None]
    c0 = np.fft.fft(state.phi0 / phase, axis=0)
    c1 = np.fft.fft(state.phi1 / phase, axis=0)
    wave = (np.fft.fftfreq(count) * count + .5) * 2 * np.pi / length
    mass = kappa * float(state.Q[0])
    omega = np.sqrt(wave ** 2 + mass ** 2)
    cosine = np.cos(time_value * omega)[:, None]
    sine = (np.sin(time_value * omega) / omega)[:, None]
    out0 = cosine * c0 - 1j * sine * (mass - 1j * wave)[:, None] * c1
    out1 = cosine * c1 - 1j * sine * (mass + 1j * wave)[:, None] * c0
    return phase * np.fft.ifft(out0, axis=0), phase * np.fft.ifft(out1, axis=0)


def occupation(observer, phi0, phi1, weights):
    projected = observer.conj().T @ np.vstack((phi0, phi1))
    return np.sum(abs(projected) ** 2 * weights[None, :], axis=1)


def observe(grid, state, initial, observer, time_value):
    rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    fine, system, source = bundle["fine_state"], bundle["fine_system"], bundle["source"]
    transport = galerkin.conformal_constraint_transport(grid, state, rate, bundle, include_vectors=True)
    h, d = transport["h"], transport["d"]
    c = h / fine.Q
    constraints = galerkin.constraint_diagnostics(grid, state, source, fine)
    r_dot = galerkin.prolong_geometry(grid, rate.r)
    proper = r_dot / (fine.r * fine.Q)
    gram = state.phi0.conj().T @ state.phi0 + state.phi1.conj().T @ state.phi1
    frozen0, frozen1 = frozen_columns(initial, grid.length, system.kappa, time_value)
    center_index = int(round(grid.nq / grid.length))  # x=1, center of region 0
    field = field_energy(system, fine, (source["image0"], source["image1"]))
    gravity = gravity_energy(system, fine)
    report = {
        "time": float(time_value), "field_energy": field, "gravity_energy": gravity, "energy": field + gravity,
        "fieldwork_power": rate.fieldwork_power,
        "constraint_norm": transport["constraint_norm"], "forcing_norm": transport["forcing_norm"],
        "quadrature_forcing_norm": transport["quadrature_forcing_norm"],
        "projection_forcing_norm": transport["projection_forcing_norm"],
        "constraint_energy_rate": transport["constraint_energy_rate"], "sbp_identity_error": transport["sbp_identity_error"],
        "full_C_l2": float(np.sqrt(grid.dx_q * np.sum(c ** 2))),
        "full_D_l2": float(np.sqrt(grid.dx_q * np.sum(d ** 2))),
        "full_C_D_l2": float(np.sqrt(grid.dx_q * np.sum(c ** 2 + d ** 2))),
        "full_hamilton_max": constraints["full_hamilton_max"], "full_momentum_max": constraints["full_momentum_max"],
        "proper_min": float(np.min(proper)), "proper_max": float(np.max(proper)), "proper_mean": float(np.mean(proper)),
        "r_min": float(np.min(fine.r)), "r_max": float(np.max(fine.r)),
        "Q_min": float(np.min(fine.Q)), "Q_max": float(np.max(fine.Q)),
        "gram_gap": float(np.max(abs(gram - np.eye(6)))),
        "omega": galerkin.subspace_frequency(grid, state)[0],
        "clock_rate": float(fine.r[center_index] * fine.Q[center_index]),
        "occupation": occupation(observer, state.phi0, state.phi1, system.occupations),
        "frozen_occupation": occupation(observer, frozen0, frozen1, system.occupations),
    }
    return report


def integrate(values, dt):
    """Independent trapezoid/paired Simpson indicators on uniform samples."""
    values = np.asarray(values, dtype=float)
    trapezoid = float(dt * (np.sum(values[1:-1]) + .5 * (values[0] + values[-1])))
    if (len(values) - 1) % 2:
        return {"trapezoid": trapezoid, "simpson": None, "indicator": None}
    simpson = float(dt / 3 * (values[0] + values[-1] + 4 * np.sum(values[1:-1:2]) + 2 * np.sum(values[2:-1:2])))
    return {"trapezoid": trapezoid, "simpson": simpson, "indicator": abs(simpson - trapezoid)}


def select_common_duration(elapsed, costs, target=TARGET):
    forecast = sum(math.ceil(target / spec["dt"]) * costs[spec["nf"]] for spec in RUN_PLAN)
    remainder = max(0., CPU_BUDGET - elapsed - 5.)
    duration = target if forecast * FORECAST_SAFETY <= remainder else math.floor(target * remainder / (forecast * FORECAST_SAFETY) / FRAME) * FRAME
    return max(0., min(target, duration)), forecast


def summarize(samples, dt, completed):
    first, last = samples[0], samples[-1]
    integrals = {name: integrate([row[name] for row in samples], dt) for name in
                 ("forcing_norm", "quadrature_forcing_norm", "projection_forcing_norm", "fieldwork_power", "clock_rate")}
    step_force = .5 * dt * (np.array([row["forcing_norm"] for row in samples[1:]]) + np.array([row["forcing_norm"] for row in samples[:-1]]))
    cumulative = np.concatenate(([0.], np.cumsum(step_force)))
    clock_values = np.array([row["clock_rate"] for row in samples])
    clocks = np.concatenate(([0.], np.cumsum(.5 * dt * (clock_values[:-1] + clock_values[1:]))))
    for row, clock in zip(samples, clocks):
        row["proper_clock"] = float(clock)
    constraint_end_budget = first["constraint_norm"] + integrals["forcing_norm"]["trapezoid"]
    field_change = last["field_energy"] - first["field_energy"]
    gravity_change = last["gravity_energy"] - first["gravity_energy"]
    frozen_clock_rate = first["clock_rate"]
    comparison_clock = min(clocks[-1], samples[-1]["time"] * frozen_clock_rate)
    occupations = np.array([row["occupation"] for row in samples])
    frozen_occ = np.array([row["frozen_occupation"] for row in samples])
    frozen_clocks = np.array([row["time"] * frozen_clock_rate for row in samples])
    matched_coupled = np.array([np.interp(comparison_clock, clocks, occupations[:, column]) for column in range(2)])
    matched_frozen = np.array([np.interp(comparison_clock, frozen_clocks, frozen_occ[:, column]) for column in range(2)])
    return {
        "completed": completed, "time": last["time"], "dt": dt, "steps": len(samples) - 1,
        "initial": first, "final": last, "integrals": integrals,
        "field_energy_change": field_change, "gravity_energy_change": gravity_change,
        "energy_balance_change": field_change + gravity_change,
        "fieldwork_difference": field_change - integrals["fieldwork_power"]["trapezoid"],
        "constraint_end_sampled_budget": constraint_end_budget,
        "constraint_end_margin": constraint_end_budget - last["constraint_norm"],
        "max_sampled_budget_violation": float(max(row["constraint_norm"] - first["constraint_norm"] - forcing for row, forcing in zip(samples, cumulative))),
        "positive_chart": all(row["r_min"] > 0 and row["Q_min"] > 0 for row in samples),
        "r_min_window": min(row["r_min"] for row in samples), "Q_min_window": min(row["Q_min"] for row in samples),
        "gram_gap_max": max(row["gram_gap"] for row in samples),
        "energy_drift_max": max(abs(row["energy"] - first["energy"]) for row in samples),
        "sbp_identity_error_max": max(row["sbp_identity_error"] for row in samples),
        "occupation_change": last["occupation"] - first["occupation"],
        "coupled_minus_frozen_at_equal_coordinate_time": last["occupation"] - last["frozen_occupation"],
        "proper_clock_comparison": {"time": float(comparison_clock), "coupled": matched_coupled,
                                    "frozen": matched_frozen, "difference": matched_coupled - matched_frozen,
                                    "sampling_interpolation_certified": False},
        "continuum_constraint_certified": False, "forcing_time_integral_certified": False,
        "renewal": False,
    }


def comparison_rows(results):
    keys = ("field_energy_change", "gravity_energy_change", "constraint_end_margin", "energy_balance_change")
    rows = []
    for nf in (256, 512):
        a, b = results[f"nf{nf}_dt_0.0010"], results[f"nf{nf}_dt_0.0005"]
        rows.append({"comparison": f"time_nf{nf}", **{key: abs(a[key] - b[key]) for key in keys},
                     "occupation_max": float(np.max(abs(a["final"]["occupation"] - b["final"]["occupation"]))),
                     "constraint_norm_difference": abs(a["final"]["constraint_norm"] - b["final"]["constraint_norm"]),
                     "forcing_integral_difference": abs(a["integrals"]["forcing_norm"]["trapezoid"] - b["integrals"]["forcing_norm"]["trapezoid"])})
    a, b = results["nf256_dt_0.0005"], results["nf512_dt_0.0005"]
    rows.append({"comparison": "space_dt_0.0005", **{key: abs(a[key] - b[key]) for key in keys},
                 "occupation_max": float(np.max(abs(a["final"]["occupation"] - b["final"]["occupation"]))),
                 "constraint_norm_difference": abs(a["final"]["constraint_norm"] - b["final"]["constraint_norm"])})
    return rows


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def save(record, arrays):
    raw_bytes = sum(np.asarray(value).nbytes for value in arrays.values())
    if raw_bytes > PAYLOAD_LIMIT:
        raise ValueError("uncompressed arrays exceed payload budget")
    with NPZ.with_suffix(".npz.tmp").open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    os.replace(NPZ.with_suffix(".npz.tmp"), NPZ)
    record["payload_sha256"] = sha256(NPZ)
    record["payload_uncompressed_array_bytes"] = raw_bytes
    for _ in range(2):
        OUT.write_text(json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n")
        record["payload_bytes"] = OUT.stat().st_size + NPZ.stat().st_size
    OUT.write_text(json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n")
    if OUT.stat().st_size + NPZ.stat().st_size > PAYLOAD_LIMIT:
        raise ValueError("saved payload exceeds 64MiB")


def run():
    started = time.process_time()
    bound_paths = {"v5_json": legacy.V5_JSON, "v5_npz": legacy.V5_NPZ,
                   "legacy_episode_json": legacy.OUT, "legacy_episode_npz": legacy.NPZ,
                   "driver": DRIVER, "galerkin": galerkin._MODULE_PATH,
                   "coupling": galerkin._COUPLING_PATH, "action": galerkin._ACTION_PATH,
                   "source": galerkin._SOURCE_PATH, "proof_note": NOTE}
    before = {name: sha256(path) for name, path in bound_paths.items()}
    arrays, loaded, costs = {}, {}, {}
    record = {"schema": SCHEMA, "gauge": "conformal", "run_plan": RUN_PLAN,
              "target_time": TARGET, "cpu_budget_seconds": CPU_BUDGET, "payload_limit_bytes": PAYLOAD_LIMIT,
              "constraints_varied_before_gauge": True, "lapse": "L=Q", "shift": "beta=0",
              "initial_data": "unchanged saved T=0 v5 Cauchy fields and columns", "initial_reset": False,
              "source_changed": False, "radius_model_changed": False, "renewal": False,
              "constraint_transport": "hc_t=D_x+f_h; D_t=hc_x+f_D",
              "constraint_bound": "R(t)<=R(0)+integral||f||_2 dt",
              "spatial_quadrature_indicator": "unprojected fine-grid constraint transport defect; not a continuum enclosure",
              "time_integral_indicator": "paired Simpson minus trapezoid and independent dt refinement; no validated integral",
              "clock_protocol": "normal worldline at fixed x=1; d tau=r(x=1) Q(x=1) dt; frozen control uses same initial clock",
              "observer": "original fixed T=0 region-0 pair; weighted occupation without isotropic copy factor",
              "frozen_control": "exact constant-Q initial conformal Dirac group, same Gaussian columns and initial geometry",
              "work": "sum(FQ Qdot+FL Ldot), Ldot=lifted Qdot, beta_dot=0",
              "source_hashes_before": before, "results": {}, "series": {}}
    for nf in (256, 512):
        build_start = time.process_time()
        grid, state, metadata = legacy.load_v5_state(f"nf{nf}")
        grid.gauge = "conformal"  # New private grid, no shared calibration arrays are changed.
        observer = np.vstack((state.phi0, state.phi1))[:, :2].copy()
        initial_sample = observe(grid, state, state, observer, 0.)
        pilot_start = time.process_time()
        trial = galerkin.rk4_step(grid, state, .0005)
        observe(grid, trial, state, observer, .0005)
        costs[nf] = time.process_time() - pilot_start
        loaded[nf] = (grid, state, observer, initial_sample)
        record[f"setup_nf{nf}"] = {"gauge": "conformal", "build_and_pilot_cpu": time.process_time() - build_start,
                                    "full_diagnostic_step_cpu": costs[nf], "v5_sha256": metadata["v5_sha256"],
                                    "occupation_gap": metadata["occupation_gap"], "quadrature_radius_gap": metadata["quadrature_radius_gap"],
                                    "initial_full_hamilton_max": initial_sample["full_hamilton_max"],
                                    "initial_full_momentum_max": initial_sample["full_momentum_max"]}
        for name in STATE_NAMES:
            arrays[f"nf{nf}_initial_{name}"] = getattr(state, name)
        print("measured_step", nf, costs[nf], "cpu_total", time.process_time() - started, flush=True)
    elapsed = time.process_time() - started
    duration, forecast = select_common_duration(elapsed, costs)
    record["forecast"] = {"elapsed_cpu": elapsed, "requested_four_case_cpu": forecast, "safety_factor": FORECAST_SAFETY,
                          "selected_common_time": duration, "budget_margin_forecast": CPU_BUDGET - elapsed - FORECAST_SAFETY * forecast}
    print("forecast", record["forecast"], flush=True)
    if duration < FRAME:
        record["verdict"] = "BUDGET_BLOCKED_AFTER_MEASURED_FORECAST"
        record["cpu_seconds"] = time.process_time() - started
        save(record, arrays)
        return record
    for spec in RUN_PLAN:
        grid, initial, observer, initial_sample = loaded[spec["nf"]]
        current = initial.copy()
        samples = [dict(initial_sample)]
        frames = {name: [getattr(current, name).copy()] for name in STATE_NAMES}
        frame_times = [0.]
        actual_times = np.linspace(0., duration, int(round(duration / spec["dt"])) + 1)
        stopped = None
        cost = costs[spec["nf"]]
        run_start = time.process_time()
        for time_value in actual_times[1:]:
            if time.process_time() - started + max(1., FORECAST_SAFETY * cost) + 3. >= CPU_BUDGET:
                stopped = "CPU_BUDGET"
                break
            step_start = time.process_time()
            try:
                candidate = galerkin.rk4_step(grid, current, spec["dt"])
                sample = observe(grid, candidate, initial, observer, time_value)
                if sample["r_min"] <= 0 or sample["Q_min"] <= 0:
                    raise PositiveChartExit("positive_chart", time_value, candidate)
                if spec["dt"] * sample["omega"] > 2 * np.sqrt(2):
                    stopped = "ABSOLUTE_RK4_IMAGINARY_LIMIT"
                    break
                if sample["gram_gap"] > 1e-6:
                    stopped = "GRAM_DEFECT_EXCEEDS_1e-6"
                    break
                current = candidate
                samples.append(sample)
                if abs(time_value / FRAME - round(time_value / FRAME)) < 1e-8:
                    frame_times.append(float(time_value))
                    for name in STATE_NAMES:
                        frames[name].append(getattr(current, name).copy())
            except PositiveChartExit as error:
                stopped = error.reason
                break
            cost = max(cost, time.process_time() - step_start)
        summary = summarize(samples, spec["dt"], stopped is None and abs(samples[-1]["time"] - duration) < 1e-12)
        summary["stop_reason"] = stopped
        summary["cpu_seconds"] = time.process_time() - run_start
        record["results"][spec["name"]] = summary
        record["series"][spec["name"]] = samples
        arrays[spec["name"] + "_times"] = np.array([row["time"] for row in samples])
        arrays[spec["name"] + "_frame_times"] = np.array(frame_times)
        for name in STATE_NAMES:
            arrays[spec["name"] + "_frames_" + name] = np.array(frames[name])
            arrays[spec["name"] + "_final_" + name] = getattr(current, name)
        print("run", spec["name"], "T", summary["time"], "CPU", summary["cpu_seconds"], "exchange", summary["field_energy_change"], "R", summary["final"]["constraint_norm"], "margin", summary["constraint_end_margin"], flush=True)
    all_completed = all(row["completed"] for row in record["results"].values())
    record["verdict"] = "MEASURED_SAME_ACTION_CONSTRAINT_FORCING" if all_completed else "PARTIAL_DIAGNOSTIC_CHECKPOINT"
    record["comparisons"] = comparison_rows(record["results"]) if all_completed else []
    record["full_target_reached"] = all_completed and duration == TARGET
    record["continuum_constraint_certified"] = False
    record["observable_error_certified"] = False
    record["cpu_seconds"] = time.process_time() - started
    record["source_hashes_after"] = {name: sha256(path) for name, path in bound_paths.items()}
    record["bound_sources_unchanged_during_run"] = record["source_hashes_after"] == before
    save(record, arrays)
    print("saved", OUT, "cpu", record["cpu_seconds"], "payload", record["payload_bytes"], flush=True)
    return record


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise SystemExit("No arguments: writes only the named new v1 diagnostic outputs")
    run()
