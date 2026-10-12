#!/usr/bin/env python3
"""Population, frozen-geometry, and one continuation of the saved T=0.05 episode.

The episode driver is not modified. Rates, the column source, and the saved
trajectory stay with their owners. Coarse work stops at 600 CPU seconds.
A finer run is started only when its result can change the renewal conclusion.
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

try:
    from threadpoolctl import threadpool_limits

    threadpool_limits(limits=1).__enter__()
    THREADPOOL_LIMIT = 1
except Exception as error:
    THREADPOOL_LIMIT = "not applied: " + type(error).__name__

import derive_nsc_spherical_feedback_episode as episode
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_regeneration_controls import (
    BASELINE_OCCUPATIONS,
    CAR_GAP_MAX,
    CONTENT_EFFECT,
    CONTROL_FLOORS,
    CRITERION,
    PERTURBATIONS_FROM_PROXY,
    REVERSED_OCCUPATIONS,
    SHARE_STOP,
    UNIFORM_OCCUPATIONS,
    assess_windows,
    bracket_g,
    car_report,
    chart_from_fields,
    dynamic_stop_reason,
    frozen_rk4_step,
    galerkin_half_identity_report,
    geometry_bytes,
    imbalance_perturbation,
    join_series,
    load_episode_final,
    mode_energies,
    prepare_population,
    proper_clock_rates,
    separate_ledgers,
    shell_views,
    source_arrays,
    state_sha256,
    uniform_complement_localization,
)
from recursive_horizons.nsc_spherical_coupling import CauchyState, PositiveChartExit, chart_failure

LAB = episode.LAB
OUT = LAB / "results" / "development" / "nsc-regeneration-controls-v1.json"
NPZ = LAB / "results" / "development" / "nsc-regeneration-controls-v1.npz"
OUT_V2 = LAB / "results" / "development" / "nsc-regeneration-controls-v2.json"
NOTE = LAB / "docs" / "nsc-regeneration-controls.md"
SCHEMA = "NSC-REGENERATION-CONTROLS-v1"
SCHEMA_V2 = "NSC-REGENERATION-CONTROLS-v2"
REVIEW = "01787780-8ed3-4fea-837d-ec21c899939a"
COARSE_BUDGET = 600.0
FINE_BUDGET = 1800.0
DT = 0.0005
DURATION = 0.05
JOIN = 0.05
FRAME = 0.005
STEP_COST = {512: 47.301563 / 100.0, 256: 11.520004 / 100.0}
RUN_NAME = {512: "nf512_dt_0_0005", 256: "nf256_dt_0_0005"}
SCALAR_KEYS = (
    "energy", "field_energy", "gravity_energy", "coordinate_work", "proper_work",
    "diagnosed_work", "diagnosed_proper_work", "proper_min", "proper_max", "proper_mean",
    "Q_min", "r_min", "r_max", "chi_min", "chi_max", "full_h", "projected_h", "held_h",
    "full_d", "projected_d", "held_d", "mean_rho", "mean_current", "G_min",
    "proper_clock", "leader_clock", "leader", "leader_share", "leader_content",
    "packet_flux", "reservoir_flux", "realized_r", "realized_Q", "realized_chi",
    "unitarity", "number",
)


CLOCK = None


class Clock:
    def __init__(self, already=0.0):
        self.start = time.process_time()
        self.already = float(already)

    def __call__(self):
        return self.already + time.process_time() - self.start


def estimate_cpu(fermions, duration=DURATION, dt=DT):
    steps = int(np.ceil(float(duration) / float(dt) - 1e-12))
    return float(STEP_COST[int(fermions)] * steps)


def checkpoint(record, arrays):
    record["cpu_seconds"] = float(CLOCK())
    episode.write_json(OUT, record)
    if arrays:
        episode.write_npz(NPZ, arrays)


def declare(record, arrays, text, estimated, kind):
    record.setdefault("declarations", []).append({
        "when_cpu": float(CLOCK()),
        "kind": kind,
        "estimated_cpu": float(estimated),
        "text": text,
    })
    print("DECLARE", kind, round(estimated, 2), text, flush=True)
    checkpoint(record, arrays)


def input_hashes():
    return {
        "episode_json": episode.sha256(episode.OUT),
        "episode_npz": episode.sha256(episode.NPZ),
        "v5_json": episode.sha256(episode.V5_JSON),
        "v5_npz": episode.sha256(episode.V5_NPZ),
    }


def saved_episode():
    return json.loads(episode.OUT.read_text())


def saved_exchange():
    balance = saved_episode()["energy_balance"]["nf512_dt_0_0005"]
    return abs(float(balance["field_energy_change"]))


def saved_series(fermions):
    run = RUN_NAME[int(fermions)]
    with np.load(episode.NPZ, allow_pickle=False) as data:
        return {
            "time": np.array(data[run + "_time"], dtype=float),
            "content": np.array(data[run + "_window_normal"], dtype=float),
            "flux": np.array(data[run + "_window_proper_flux"], dtype=float),
            "field_energy": np.array(data[run + "_field_energy"], dtype=float),
            "energy": np.array(data[run + "_energy"], dtype=float),
            "chi_max": np.array(data[run + "_chi_max"], dtype=float),
            "proper_max": np.array(data[run + "_proper_max"], dtype=float),
            "proper_min": np.array(data[run + "_proper_min"], dtype=float),
            "Q_min": np.array(data[run + "_Q_min"], dtype=float),
            "r_min": np.array(data[run + "_r_min"], dtype=float),
            "r_max": np.array(data[run + "_r_max"], dtype=float),
            "mean_current": np.array(data[run + "_mean_current"], dtype=float),
            "mean_rho": np.array(data[run + "_mean_rho"], dtype=float),
            "full_h": np.array(data[run + "_full_hamilton_max"], dtype=float),
            "coordinate_work": np.array(data[run + "_coordinate_work_integral"], dtype=float),
        }


def baseline_bracket(fermions, grid):
    run = RUN_NAME[int(fermions)]
    with np.load(episode.NPZ, allow_pickle=False) as data:
        rho = np.array(data[run + "_frame_rho"][0], dtype=float)
        conformal = np.array(data[run + "_frame_Q"][0], dtype=float)
        radius = np.array(data[run + "_frame_r"][0], dtype=float)
    gee = bracket_g(conformal, rho, grid)
    return {
        "G_min": float(np.min(gee)),
        "G_max": float(np.max(gee)),
        "positive_G": bool(np.min(gee) > 0.0),
        "r_min": float(np.min(radius)),
        "Q_min": float(np.min(conformal)),
        "positive_r": bool(np.min(radius) > 0.0),
        "positive_Q": bool(np.min(conformal) > 0.0),
        "source": "saved episode frame at T=0",
    }


def observation_summary(grid, state):
    windows, partition = episode.windows_for(grid)
    sample, _rate, fields = episode.observe(grid, state, windows, partition)
    view = shell_views(sample)
    chart = chart_from_fields(grid, fields)
    occupations, totals, leader_energy = mode_energies(sample, view["leader"])
    global_rate, leader_rate = proper_clock_rates(grid, fields, windows, view["leader"])
    return {
        "energy": float(sample["energy"]),
        "field_energy": float(sample["field_energy"]),
        "gravity_energy": float(sample["gravity_energy"]),
        "full_h": float(sample["full_hamilton_max"]),
        "projected_h": float(sample["projected_hamilton_max"]),
        "held_h": float(sample["held_out_hamilton_max"]),
        "full_d": float(sample["full_momentum_max"]),
        "projected_d": float(sample["projected_momentum_max"]),
        "held_d": float(sample["held_out_momentum_max"]),
        "mean_rho": float(sample["mean_rho"]),
        "mean_current": float(sample["mean_current"]),
        "mean_subtracted": False,
        "Q_min": float(sample["Q_min"]),
        "r_min": float(sample["r_min"]),
        "chi_max": float(sample["chi_max"]),
        "proper_min": float(sample["proper_min"]),
        "proper_max": float(sample["proper_max"]),
        "number": float(sample["number"]),
        "unitarity": float(sample["unitarity"]),
        "window_normal": view["content"].tolist(),
        "window_proper_flux": view["flux"].tolist(),
        "shares": view["shares"].tolist(),
        "leader": int(view["leader"]),
        "leader_share": float(view["leader_share"]),
        "packet_flux": float(view["packet_flux"]),
        "reservoir_flux": float(view["reservoir_flux"]),
        "G_min": chart["G_min"],
        "positive_G": chart["positive_G"],
        "positive_r": chart["positive_r"],
        "positive_Q": chart["positive_Q"],
        "mode_occupations": occupations.tolist(),
        "mode_energy": totals.tolist(),
        "mode_leader_energy": leader_energy.tolist(),
        "proper_clock_rate": global_rate,
        "leader_clock_rate": leader_rate,
        "car": car_report(state),
    }


def _frame(absolute):
    if absolute <= 1e-15:
        return True
    nearest = round(float(absolute) / FRAME) * FRAME
    return abs(float(absolute) - nearest) <= 1e-10


def integrate(grid, state, duration, dt, budget, elapsed, *, t0=0.0, frozen=False,
              stop_on_delocalization=False, label="run"):
    """Owned RK4, or the same column rates with geometry velocities removed."""
    steps = int(np.ceil(float(duration) / float(dt) - 1e-12))
    dt_used = float(duration) / steps
    windows, partition = episode.windows_for(grid)
    current = state.copy()
    start_geometry = geometry_bytes(current)
    occupation_bytes = np.ascontiguousarray(grid.fine.occupations).tobytes()
    radius0 = np.array(current.r, copy=True)
    conformal0 = np.array(current.Q, copy=True)
    chi0 = np.array(current.chi, copy=True)
    rows = []
    checkpoints = []
    coordinate_work = 0.0
    proper_work = 0.0
    diagnosed_work = 0.0
    diagnosed_proper_work = 0.0
    proper_clock = 0.0
    leader_clock = 0.0
    previous_power = None
    previous_proper = None
    previous_global = None
    previous_leader = None
    refused_max = 0.0
    stop_reason = None
    chart_stop = False
    completed_steps = 0
    step_cost = 0.0
    cpu_start = time.process_time()
    error_text = None
    start_sha = state_sha256(current)
    pending_q_change = 0.0
    gram_gaps = []
    last_checkpoint = None

    def remember(sample, fields, absolute, applied_coordinate, applied_proper):
        nonlocal proper_clock, leader_clock, last_checkpoint
        view = shell_views(sample)
        chart = chart_from_fields(grid, fields)
        gram = car_report(current)
        gram_gaps.append(float(gram["gap"]))
        global_rate, leader_rate = proper_clock_rates(grid, fields, windows, view["leader"])
        row = {
            "time": float(absolute),
            "energy": float(sample["energy"]),
            "field_energy": float(sample["field_energy"]),
            "gravity_energy": float(sample["gravity_energy"]),
            "coordinate_work": float(applied_coordinate),
            "proper_work": float(applied_proper),
            "diagnosed_work": float(diagnosed_work),
            "diagnosed_proper_work": float(diagnosed_proper_work),
            "proper_min": float(sample["proper_min"]),
            "proper_max": float(sample["proper_max"]),
            "proper_mean": float(sample["proper_mean"]),
            "Q_min": float(sample["Q_min"]),
            "r_min": float(sample["r_min"]),
            "r_max": float(sample["r_max"]),
            "chi_min": float(sample["chi_min"]),
            "chi_max": float(sample["chi_max"]),
            "full_h": float(sample["full_hamilton_max"]),
            "projected_h": float(sample["projected_hamilton_max"]),
            "held_h": float(sample["held_out_hamilton_max"]),
            "full_d": float(sample["full_momentum_max"]),
            "projected_d": float(sample["projected_momentum_max"]),
            "held_d": float(sample["held_out_momentum_max"]),
            "mean_rho": float(sample["mean_rho"]),
            "mean_current": float(sample["mean_current"]),
            "G_min": float(chart["G_min"]),
            "positive_G": bool(chart["positive_G"]),
            "proper_clock": float(proper_clock),
            "leader_clock": float(leader_clock),
            "leader": int(view["leader"]),
            "leader_share": float(view["leader_share"]),
            "leader_content": float(view["content"][view["leader"]]),
            "packet_flux": float(view["packet_flux"]),
            "reservoir_flux": float(view["reservoir_flux"]),
            "realized_r": float(np.max(np.abs(current.r - radius0))),
            "realized_Q": float(np.max(np.abs(current.Q - conformal0))),
            "realized_chi": float(np.max(np.abs(current.chi - chi0))),
            "unitarity": float(sample["unitarity"]),
            "number": float(sample["number"]),
            "gram_gap": float(gram["gap"]),
            "Q_change": float(pending_q_change),
            "positive_L": bool(chart["positive_L"]),
            "content": np.array(view["content"], dtype=float),
            "flux": np.array(view["flux"], dtype=float),
            "global_rate": float(global_rate),
            "leader_rate": float(leader_rate),
            "lifted_fieldwork": float(sample["lifted_fieldwork"]),
            "proper_pressure_power": float(sample["proper_pressure_power"]),
        }
        if _frame(absolute) or abs(absolute - (t0 + duration)) <= 1e-10:
            occupations, totals, leader_energy = mode_energies(sample, view["leader"])
            checkpoints.append({
                key: row[key] for key in row if key not in ("content", "flux", "global_rate", "leader_rate",
                                                             "lifted_fieldwork", "proper_pressure_power")
            })
            checkpoints[-1]["window_normal"] = view["content"].tolist()
            checkpoints[-1]["window_proper_flux"] = view["flux"].tolist()
            checkpoints[-1]["shares"] = view["shares"].tolist()
            checkpoints[-1]["mode_occupations"] = occupations.tolist()
            checkpoints[-1]["mode_energy"] = totals.tolist()
            checkpoints[-1]["mode_leader_energy"] = leader_energy.tolist()
            checkpoints[-1]["positive_G"] = bool(chart["positive_G"])
            checkpoints[-1]["positive_r"] = bool(chart["positive_r"])
            checkpoints[-1]["positive_Q"] = bool(chart["positive_Q"])
            checkpoints[-1]["positive_L"] = bool(chart["positive_L"])
            checkpoints[-1]["exit_reason"] = None
            checkpoints[-1]["kept"] = True
            last_checkpoint = checkpoints[-1]
        else:
            last_checkpoint = {
                "time": float(absolute),
                "positive_G": bool(chart["positive_G"]),
                "positive_r": bool(chart["positive_r"]),
                "positive_Q": bool(chart["positive_Q"]),
                "positive_L": bool(chart["positive_L"]),
                "G_min": float(chart["G_min"]),
                "Q_min": float(sample["Q_min"]),
                "r_min": float(sample["r_min"]),
                "gram_gap": float(gram["gap"]),
                "Q_change": float(pending_q_change),
                "leader_share": float(view["leader_share"]),
                "exit_reason": None,
                "kept": True,
            }
        rows.append(row)
        return chart

    def mark_exit(reason):
        if last_checkpoint is None:
            checkpoints.append({
                "time": float(t0),
                "exit_reason": reason,
                "kept": False,
                "chart_stop": True,
            })
            return
        if not checkpoints or abs(float(checkpoints[-1].get("time", -1.0)) - float(last_checkpoint["time"])) > 1e-12:
            checkpoints.append(dict(last_checkpoint))
        checkpoints[-1]["exit_reason"] = reason
        checkpoints[-1]["kept"] = True

    def open_reason(sample, share):
        gram = car_report(current)
        if float(sample["unitarity"]) > CAR_GAP_MAX:
            gram = dict(gram)
            gram["admissible"] = False
        delocalized = bool(stop_on_delocalization and share < SHARE_STOP)
        return dynamic_stop_reason(None, gram, delocalized=delocalized)

    try:
        try:
            sample, _rate, fields = episode.observe(grid, current, windows, partition)
        except PositiveChartExit as exit_chart:
            stop_reason = str(exit_chart.reason)
            chart_stop = True
            mark_exit(stop_reason)
        else:
            pending_q_change = 0.0
            chart = remember(sample, fields, t0, 0.0, 0.0)
            previous_power = rows[-1]["lifted_fieldwork"]
            previous_proper = rows[-1]["proper_pressure_power"]
            previous_global = rows[-1]["global_rate"]
            previous_leader = rows[-1]["leader_rate"]
            reason, is_chart = open_reason(sample, rows[-1]["leader_share"])
            if reason:
                stop_reason = reason
                chart_stop = is_chart
                mark_exit(reason)
            for step_index in range(steps):
                if stop_reason:
                    break
                if not episode.budget_allows(elapsed(), step_cost if step_index else 0.0, budget):
                    stop_reason = "CPU_BUDGET"
                    mark_exit(stop_reason)
                    break
                started = time.process_time()
                try:
                    if frozen:
                        candidate, refused = frozen_rk4_step(grid, current, dt_used)
                        refused_max = max(refused_max, abs(refused))
                    else:
                        candidate = galerkin.rk4_step(grid, current, dt_used)
                except PositiveChartExit as exit_chart:
                    stop_reason = str(exit_chart.reason)
                    chart_stop = True
                    mark_exit(stop_reason)
                    break
                failure = chart_failure(grid.fine, galerkin.prolong_state(grid, candidate))
                if failure:
                    stop_reason = str(failure)
                    chart_stop = True
                    mark_exit(stop_reason)
                    break
                pending_q_change = float(np.max(np.abs(np.asarray(candidate.Q) - np.asarray(current.Q))))
                try:
                    sample, _rate, fields = episode.observe(grid, candidate, windows, partition)
                except PositiveChartExit as exit_chart:
                    stop_reason = str(exit_chart.reason)
                    chart_stop = True
                    mark_exit(stop_reason)
                    break
                current = candidate
                if frozen and geometry_bytes(current) != start_geometry:
                    stop_reason = "FROZEN_GEOMETRY_DRIFT"
                    mark_exit(stop_reason)
                    break
                if np.ascontiguousarray(grid.fine.occupations).tobytes() != occupation_bytes:
                    stop_reason = "OCCUPATIONS_CHANGED"
                    mark_exit(stop_reason)
                    break
                absolute = t0 + (step_index + 1) * dt_used
                diagnosed_work += 0.5 * dt_used * (previous_power + float(sample["lifted_fieldwork"]))
                diagnosed_proper_work += 0.5 * dt_used * (previous_proper + float(sample["proper_pressure_power"]))
                previous_power = float(sample["lifted_fieldwork"])
                previous_proper = float(sample["proper_pressure_power"])
                if not frozen:
                    coordinate_work = diagnosed_work
                    proper_work = diagnosed_proper_work
                view = shell_views(sample)
                global_rate, leader_rate = proper_clock_rates(grid, fields, windows, view["leader"])
                proper_clock += 0.5 * dt_used * (previous_global + global_rate)
                leader_clock += 0.5 * dt_used * (previous_leader + leader_rate)
                previous_global = global_rate
                previous_leader = leader_rate
                remember(sample, fields, absolute, coordinate_work, proper_work)
                completed_steps = step_index + 1
                step_cost = time.process_time() - started
                reason, is_chart = open_reason(sample, rows[-1]["leader_share"])
                if reason:
                    stop_reason = reason
                    chart_stop = chart_stop or is_chart
                    mark_exit(reason)
                if step_index == 0 or completed_steps % 10 == 0 or completed_steps == steps or stop_reason:
                    print(
                        label, "step", completed_steps, "T", round(absolute, 6),
                        "cpu", round(time.process_time() - cpu_start, 3),
                        "share", round(rows[-1]["leader_share"], 4),
                        "chi", round(rows[-1]["chi_max"], 4),
                        "G", rows[-1]["G_min"],
                        "Q", rows[-1]["Q_min"],
                        flush=True,
                    )
    except PositiveChartExit as exit_chart:
        stop_reason = stop_reason or str(exit_chart.reason)
        chart_stop = True
        mark_exit(stop_reason)
    except Exception:
        error_text = traceback.format_exc()
        stop_reason = stop_reason or "CODE_EXCEPTION"
        print(error_text, flush=True)
    stacked = {
        "time": np.asarray([row["time"] for row in rows], dtype=float),
        "content": np.vstack([row["content"] for row in rows]) if rows else np.zeros((0, 4)),
        "flux": np.vstack([row["flux"] for row in rows]) if rows else np.zeros((0, 4)),
    }
    for key in SCALAR_KEYS:
        stacked[key] = np.asarray([row[key] for row in rows], dtype=float)
    attained = float(stacked["time"][-1]) if rows else float(t0)
    return {
        "series": stacked,
        "checkpoints": checkpoints,
        "state": current,
        "attained_time": attained,
        "requested_time": float(t0 + duration),
        "t0": float(t0),
        "dt": float(dt_used),
        "steps_requested": steps,
        "steps_completed": completed_steps,
        "completed": bool(completed_steps == steps and stop_reason is None),
        "stop_reason": stop_reason,
        "chart_stop": chart_stop,
        "cpu_seconds": float(time.process_time() - cpu_start),
        "step_cost": float(step_cost),
        "frozen": bool(frozen),
        "geometry_unchanged": bool(frozen and geometry_bytes(current) == start_geometry),
        "occupations_unchanged": bool(np.ascontiguousarray(grid.fine.occupations).tobytes() == occupation_bytes),
        "applied_fieldwork_is_zero": bool(frozen),
        "refused_stage_power_max": float(refused_max),
        "positive_G": bool(rows and all(row["positive_G"] for row in rows)),
        "positive_r": bool(rows and np.min(stacked["r_min"]) > 0.0),
        "positive_Q": bool(rows and np.min(stacked["Q_min"]) > 0.0),
        "positive_L": bool(rows and all(row["positive_L"] for row in rows)),
        "g_used_as_dynamic_stop": False,
        "gram_gap_max": None if not gram_gaps else float(max(gram_gaps)),
        "Q_change_max": None if not rows else float(max(row["Q_change"] for row in rows)),
        "start_sha256": start_sha,
        "end_sha256": state_sha256(current),
        "error": error_text,
        "mean_source_subtracted": False,
    }


def store_series(arrays, name, result):
    series = result["series"]
    for key, value in series.items():
        arrays[name + "_" + key] = np.asarray(value)
    arrays[name + "_final_Q"] = np.asarray(result["state"].Q)
    arrays[name + "_final_r"] = np.asarray(result["state"].r)
    arrays[name + "_final_chi"] = np.asarray(result["state"].chi)
    arrays[name + "_final_p_Q"] = np.asarray(result["state"].p_Q)
    arrays[name + "_final_p_r"] = np.asarray(result["state"].p_r)
    arrays[name + "_final_p_chi"] = np.asarray(result["state"].p_chi)
    arrays[name + "_final_phi0"] = np.asarray(result["state"].phi0)
    arrays[name + "_final_phi1"] = np.asarray(result["state"].phi1)


def balances(series, reference_exchange):
    if len(series["time"]) < 1:
        return {"empty": True}
    energy_change = float(series["energy"][-1] - series["energy"][0])
    field_change = float(series["field_energy"][-1] - series["field_energy"][0])
    gravity_change = float(series["gravity_energy"][-1] - series["gravity_energy"][0])
    work = float(series["coordinate_work"][-1])
    scale = abs(float(reference_exchange))
    return {
        "energy_change": energy_change,
        "field_energy_change": field_change,
        "gravity_energy_change": gravity_change,
        "applied_coordinate_work": work,
        "proper_pressure_work": float(series["proper_work"][-1]),
        "diagnosed_coordinate_work": float(series["diagnosed_work"][-1]),
        "versus_field_change": episode.balance_against_exchange(energy_change, field_change),
        "versus_applied_work": episode.balance_against_exchange(energy_change, work),
        "saved_exchange_scale": scale,
        "energy_change_over_saved_exchange": None if scale == 0.0 else abs(energy_change) / scale,
        "field_change_over_saved_exchange": None if scale == 0.0 else abs(field_change) / scale,
    }


def public_result(result, reference_exchange):
    series = result["series"]
    if len(series["time"]) < 1:
        return {"completed": False, "stop_reason": result.get("stop_reason"), "empty": True}
    structure = assess_windows(series["content"], series["flux"]) if len(series["time"]) >= 2 else None
    return {
        "attained_time": result["attained_time"],
        "requested_time": result["requested_time"],
        "completed": result["completed"],
        "stop_reason": result["stop_reason"],
        "chart_stop": result["chart_stop"],
        "cpu_seconds": result["cpu_seconds"],
        "dt": result["dt"],
        "steps_completed": result["steps_completed"],
        "steps_requested": result["steps_requested"],
        "frozen": result["frozen"],
        "geometry_unchanged": result["geometry_unchanged"] if result["frozen"] else None,
        "occupations_unchanged": result["occupations_unchanged"],
        "mean_source_subtracted": False,
        "applied_fieldwork_is_zero": result["applied_fieldwork_is_zero"],
        "refused_stage_power_max": result["refused_stage_power_max"],
        "positive_G": result["positive_G"],
        "positive_r": result["positive_r"],
        "positive_Q": result["positive_Q"],
        "start_sha256": result["start_sha256"],
        "end_sha256": result["end_sha256"],
        "balance": balances(series, reference_exchange),
        "structure": structure,
        "leader_content_change": None if structure is None else structure["leader_content_change"],
        "field_energy_change": float(series["field_energy"][-1] - series["field_energy"][0]),
        "chi_max_change": float(series["chi_max"][-1] - series["chi_max"][0]),
        "proper_max_change": float(series["proper_max"][-1] - series["proper_max"][0]),
        "realized_radius_end": float(series["realized_r"][-1]),
        "realized_chi_end": float(series["realized_chi"][-1]),
        "proper_clock_end": float(series["proper_clock"][-1]),
        "leader_clock_end": float(series["leader_clock"][-1]),
        "packet_flux_end": float(series["packet_flux"][-1]),
        "reservoir_flux_end": float(series["reservoir_flux"][-1]),
        "mean_current_start": float(series["mean_current"][0]),
        "mean_current_end": float(series["mean_current"][-1]),
        "mean_rho_start": float(series["mean_rho"][0]),
        "mean_rho_end": float(series["mean_rho"][-1]),
        "full_h_end": float(series["full_h"][-1]),
        "projected_h_end": float(series["projected_h"][-1]),
        "held_h_end": float(series["held_h"][-1]),
        "full_d_end": float(series["full_d"][-1]),
        "G_min_end": float(series["G_min"][-1]),
        "Q_min_end": float(series["Q_min"][-1]),
        "checkpoints": result["checkpoints"],
        "error": result["error"],
    }


def movement_row(name, primary, other, floor):
    primary = float(primary)
    other = float(other)
    movement = abs(primary - other)
    return {
        "name": name,
        "primary": primary,
        "other": other,
        "effect_scale": abs(primary),
        "movement": movement,
        "floor": float(floor),
        "status": episode.classify_movement(primary, movement, floor),
    }


def room_for(estimated, budget, elapsed):
    return bool(float(elapsed()) + float(estimated) <= float(budget))


def run_evolution(record, arrays, name, grid, state, *, budget, elapsed, frozen, t0, stop_share, reference):
    estimated = estimate_cpu(grid.nf, DURATION, DT)
    if not room_for(estimated, budget, elapsed):
        record["stages"][name] = {
            "completed": False,
            "stop_reason": "CPU_BUDGET_BEFORE_START",
            "estimated_cpu": estimated,
            "cpu_at_refusal": float(elapsed()),
        }
        checkpoint(record, arrays)
        return None
    start_note = ""
    if t0 > 0.0:
        start_note = (
            f" Starts bitwise from the saved final state {state_sha256(state)}. "
            "No reset, pulse, or extra force."
        )
    declare(
        record, arrays,
        (
            f"{name}: duration {DURATION} from T={t0}, dt {DT}, nf {grid.nf}, "
            f"frozen_geometry {frozen}. Localization stop share {SHARE_STOP if stop_share else 'off'}. "
            f"Effect floors {CONTROL_FLOORS}. Estimated CPU {estimated:.1f}.{start_note}"
        ),
        estimated,
        "continuation" if t0 > 0.0 else "evolution",
    )
    result = integrate(
        grid, state, DURATION, DT, budget, elapsed,
        t0=t0, frozen=frozen, stop_on_delocalization=stop_share, label=name,
    )
    store_series(arrays, name, result)
    checkpoint(record, arrays)
    record["stages"][name] = public_result(result, reference)
    record["stages"][name]["completed_flag"] = result["completed"]
    checkpoint(record, arrays)
    return result


def compare_pair(saved, evolved, prefix):
    """Geometry or population effect against the saved coupled episode."""
    rows = []
    if evolved is None or len(evolved["series"]["time"]) < 2:
        return rows
    series = evolved["series"]
    pairs = (
        ("leader_shell", float(saved["content"][-1, 0] - saved["content"][0, 0]),
         float(series["content"][-1, 0] - series["content"][0, 0]),
         CONTROL_FLOORS["leader_shell"]),
        ("field_energy", float(saved["field_energy"][-1] - saved["field_energy"][0]),
         float(series["field_energy"][-1] - series["field_energy"][0]),
         CONTROL_FLOORS["field_energy"]),
        ("chi_max", float(saved["chi_max"][-1] - saved["chi_max"][0]),
         float(series["chi_max"][-1] - series["chi_max"][0]),
         CONTROL_FLOORS["chi_max"]),
        ("proper_velocity_max", float(saved["proper_max"][-1] - saved["proper_max"][0]),
         float(series["proper_max"][-1] - series["proper_max"][0]),
         CONTROL_FLOORS["proper_velocity"]),
    )
    for name, coupled, control, floor in pairs:
        effect = coupled - control
        rows.append({
            "name": prefix + name,
            "coupled_change": coupled,
            "control_change": control,
            "effect_requiring_the_difference": effect,
            "floor": floor,
            "above_floor": bool(abs(effect) >= floor),
        })
    return rows


def continuation_assessment(fermions, result):
    saved = saved_series(fermions)
    series = result["series"]
    _time, content = join_series(saved["time"], saved["content"], series["time"], series["content"])
    _time, flux = join_series(saved["time"], saved["flux"], series["time"], series["flux"])
    joined = assess_windows(content, flux)
    second = assess_windows(series["content"], series["flux"]) if len(series["time"]) >= 2 else None
    accounting = {
        "field_energy_change": float(series["field_energy"][-1] - series["field_energy"][0]),
        "energy_change": float(series["energy"][-1] - series["energy"][0]),
        "reported_with_windows": False,
    }
    return {
        "joined_from_T0": joined,
        "second_episode": second,
        "ledgers": None if second is None else separate_ledgers(series["content"], series["flux"], accounting),
        "maintained_structure": bool(result["completed"] and second and second["maintained"]),
        "renewed_structure": bool(joined["renewed"]),
        "one_drift": bool(joined["one_drift"]),
        "candidate_regime": bool(joined["candidate_regime"]),
        "candidate_regime_is_programme_requirement": False,
        "splitting_one_drift_is_renewal": False,
        "splitting_one_drift_is_complete_renewal": False,
    }


def render_note(record):
    conclusion = record.get("conclusion") or {}
    lines = [
        "# Regeneration controls on the saved spherical episode",
        "",
        "Thin consumer of the owned Galerkin step and the saved T=0.05 episode.",
        "No new radius force and no mean subtraction. The episode driver was not edited.",
        "",
        "## Criterion, declared before the evolutions",
        "",
        CRITERION["localization"],
        "",
        CRITERION["maintained"],
        "",
        CRITERION["renewed"],
        "",
        CRITERION["evolution_gate"],
        "",
        f"Coarse CPU cap {COARSE_BUDGET:.0f}. Fine cap {FINE_BUDGET:.0f}, used only if a candidate or an unresolved comparison can change the conclusion.",
        "",
        "## What the runs did",
        "",
        conclusion.get("statement") or "No conclusion was stored.",
        "",
        f"CPU seconds: {record.get('cpu_seconds')}.",
        f"Maintained structure: {conclusion.get('maintained_structure')}.",
        f"Renewed structure: {conclusion.get('renewed_structure')}.",
        f"One drift: {conclusion.get('one_drift')}.",
        f"Candidate regime proxy: {conclusion.get('candidate_regime')}.",
        "That proxy is not a programme requirement.",
        f"Perturbations executed: {conclusion.get('perturbations_executed')}.",
        "G > 0 is the initial convex-radius hypothesis. The dynamic chart is r, Q, L.",
        "Successor assessment: `lab/results/development/nsc-regeneration-controls-v2.json`.",
        "",
        "Initial-data, control, and continuation numbers are in",
        "`lab/results/development/nsc-regeneration-controls-v1.json`.",
        "The episode JSON and NPZ, and the v5 files, were read and not rewritten.",
        "",
        "## Requested interface, not applied to the shared driver",
        "",
        "The episode integrator has no frozen-geometry switch, no absolute start time,",
        "and no localization stop. This driver calls `rk4_step` and `rates` directly.",
        "A later owner can add `evolve_episode(..., t0, geometry='coupled'|'frozen', stop=...)`.",
        "",
        "```sh",
        "python scripts/lab.py scripts/derive_nsc_regeneration_controls.py --check",
        "python scripts/lab.py -m pytest tests/test_nsc_regeneration_controls.py -q",
        "```",
        "",
    ]
    NOTE.write_text("\n".join(lines))


def verify_saved():
    if not OUT.is_file() or not NPZ.is_file():
        raise FileNotFoundError("regeneration control record is missing")
    record = json.loads(OUT.read_text())
    if record.get("schema") != SCHEMA:
        raise AssertionError("schema mismatch")
    current = input_hashes()
    for name, digest in current.items():
        if record.get("inputs_after", record.get("inputs_before", {})).get(name) != digest:
            raise AssertionError(f"input hash drifted: {name}")
    if record.get("inputs_before") != record.get("inputs_after"):
        raise AssertionError("inputs changed during the campaign")
    continuation = record.get("stages", {}).get("continuation_nf512")
    if continuation is None:
        raise AssertionError("primary continuation is missing")
    state = load_episode_final(episode.NPZ, "nf512_dt_0_0005")
    if state_sha256(state) != continuation["start_sha256"]:
        raise AssertionError("continuation did not start from the saved final state")
    if continuation.get("mean_source_subtracted") is not False:
        raise AssertionError("mean source was subtracted")
    frozen = record.get("stages", {}).get("frozen_nf512")
    if frozen is None or frozen.get("geometry_unchanged") is not True:
        raise AssertionError("frozen geometry control did not keep the geometry")
    conclusion = record["conclusion"]
    if conclusion.get("splitting_one_drift_is_renewal") is not False:
        raise AssertionError("drift splitting was treated as renewal")
    if conclusion.get("perturbations_executed") and not conclusion.get("candidate_regime"):
        raise AssertionError("perturbations ran without a candidate")
    return record


def campaign():
    global CLOCK
    if OUT.exists() and NPZ.exists():
        existing = json.loads(OUT.read_text())
        if existing.get("status") == "MEASURED":
            raise FileExistsError(f"refusing to overwrite {OUT.name}")
    started = time.process_time()
    CLOCK = Clock(0.0)
    reference = saved_exchange()
    record = {
        "schema": SCHEMA,
        "status": "RUNNING",
        "question": (
            "Do uniform or reversed occupations, or frozen geometry, change the "
            "saved coupled episode, and does continuing that episode from T=0.05 "
            "to T=0.10 maintain or renew a localized structure?"
        ),
        "criterion": CRITERION,
        "control_floors": CONTROL_FLOORS,
        "renewal": False,
        "mean_source_subtracted": False,
        "mean_current_subtracted": False,
        "new_dynamical_framework": False,
        "episode_driver_modified": False,
        "incoming_gate_prerequisite": False,
        "thread_limits": {
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
            "NUMEXPR_NUM_THREADS": os.environ.get("NUMEXPR_NUM_THREADS"),
            "VECLIB_MAXIMUM_THREADS": os.environ.get("VECLIB_MAXIMUM_THREADS"),
            "threadpoolctl": THREADPOOL_LIMIT,
        },
        "inputs_before": input_hashes(),
        "stages": {},
        "initial_data": {},
        "comparisons": {},
        "declarations": [],
    }
    arrays = {"schema": np.array(SCHEMA)}
    declare(
        record, arrays,
        (
            "Before any evolution: uniform [0.5]*6 and reversed "
            "[0.25,0.25,0.5,0.5,0.75,0.75] reuse the saved mode columns. "
            "The source is recomputed. Radius and shift momentum use the owned "
            "action. A Newton stall is a solver floor. The initial convex-radius "
            "hypothesis asks for positive G, r, and Q. Dynamic steps stop on "
            f"chart_failure for r, Q, and L. Projected residual gate {1e-4}, full residual "
            f"gate {1e-2}. Mean current is not subtracted. Estimated CPU under 120."
        ),
        120.0,
        "initial-data",
    )
    loaded = {}
    prepared_cache = {}
    for label, fermions in (("nf256", 256), ("nf512", 512)):
        grid, state, metadata = episode.load_v5_state(label)
        loaded[fermions] = (grid, state, metadata)
        record["initial_data"][label + "_baseline"] = {
            "occupations": [float(value) for value in metadata["occupations"]],
            "car": car_report(state),
            "bracket_from_saved_frame": baseline_bracket(fermions, grid),
            "saved_initial_mean_current": float(saved_series(fermions)["mean_current"][0]),
            "saved_final_mean_current": float(saved_series(fermions)["mean_current"][-1]),
            "source_re_solved": False,
            "trajectory_rerun": False,
        }
        _source, rho_base, current_base = source_arrays(grid, state)
        for name, occupations in (("uniform", UNIFORM_OCCUPATIONS), ("reversed", REVERSED_OCCUPATIONS)):
            owned, solved, report = prepare_population(grid, state, occupations)
            _other, rho_new, current_new = source_arrays(owned, solved)
            report["rho_gap_from_baseline"] = float(np.max(np.abs(rho_new - rho_base)))
            report["current_gap_from_baseline"] = float(np.max(np.abs(current_new - current_base)))
            report["observation"] = observation_summary(owned, solved)
            report["grid_occupations_left_baseline"] = bool(np.allclose(grid.fine.occupations, metadata["occupations"]))
            record["initial_data"][f"{label}_{name}"] = report
            prepared_cache[(fermions, name)] = (owned, solved)
            arrays[f"{label}_{name}_r"] = np.asarray(solved.r)
            arrays[f"{label}_{name}_p_Q"] = np.asarray(solved.p_Q)
            arrays[f"{label}_{name}_phi0"] = np.asarray(solved.phi0)
            print(
                "prepared", label, name, "evolve", report["evolve"],
                "floor", report["solver_floor"], "proj", report["projected_hamilton_max"],
                "full", report["full_hamilton_max"], "G", report["G_min"],
                "meanj", report["current_mean"], "cpu", round(CLOCK(), 2),
                flush=True,
            )
            checkpoint(record, arrays)
    reference = saved_exchange()
    results = {}
    for fermions in (512, 256):
        grid, state, _metadata = loaded[fermions]
        name = f"frozen_nf{fermions}"
        results[name] = run_evolution(
            record, arrays, name, grid, state.copy(),
            budget=COARSE_BUDGET, elapsed=CLOCK, frozen=True, t0=0.0,
            stop_share=False, reference=reference,
        )
    for fermions in (512, 256):
        grid, _state, _metadata = loaded[fermions]
        for population in ("uniform", "reversed"):
            report = record["initial_data"][f"nf{fermions}_{population}"]
            name = f"{population}_nf{fermions}"
            if not report["evolve"]:
                record["stages"][name] = {
                    "completed": False,
                    "stop_reason": "INITIAL_DATA_NOT_EVOLVED",
                    "blocker": report.get("blocker"),
                    "projected_hamilton_max": report.get("projected_hamilton_max"),
                    "full_hamilton_max": report.get("full_hamilton_max"),
                    "solver_floor": report.get("solver_floor"),
                    "physical_failure": report.get("physical_failure"),
                }
                checkpoint(record, arrays)
                continue
            owned, solved = prepared_cache[(fermions, population)]
            results[name] = run_evolution(
                record, arrays, name, owned, solved,
                budget=COARSE_BUDGET, elapsed=CLOCK, frozen=False, t0=0.0,
                stop_share=False, reference=reference,
            )
    for fermions in (512, 256):
        grid, _state, _metadata = loaded[fermions]
        final = load_episode_final(episode.NPZ, RUN_NAME[fermions])
        name = f"continuation_nf{fermions}"
        results[name] = run_evolution(
            record, arrays, name, grid, final,
            budget=COARSE_BUDGET, elapsed=CLOCK, frozen=False, t0=JOIN,
            stop_share=True, reference=reference,
        )
        if results[name] is not None and results[name]["start_sha256"] != state_sha256(final):
            raise RuntimeError("continuation start hash drifted before the first step")
    record["comparisons"]["frozen_nf512"] = compare_pair(saved_series(512), results.get("frozen_nf512"), "frozen_")
    record["comparisons"]["frozen_nf256"] = compare_pair(saved_series(256), results.get("frozen_nf256"), "frozen_")
    for population in ("uniform", "reversed"):
        record["comparisons"][population + "_nf512"] = compare_pair(
            saved_series(512), results.get(population + "_nf512"), population + "_"
        )
        record["comparisons"][population + "_nf256"] = compare_pair(
            saved_series(256), results.get(population + "_nf256"), population + "_"
        )
    assessments = {}
    candidate = False
    for fermions in (512, 256):
        result = results.get(f"continuation_nf{fermions}")
        if result is None or len(result["series"]["time"]) < 2:
            continue
        assessments[fermions] = continuation_assessment(fermions, result)
        record["stages"][f"continuation_nf{fermions}"]["renewal_assessment"] = assessments[fermions]
        candidate = candidate or assessments[fermions]["candidate_regime"]
    for name, result in results.items():
        if result is None or name.startswith("continuation"):
            continue
        structure = record["stages"].get(name, {}).get("structure")
        if structure and structure.get("candidate_regime"):
            candidate = True
            record["stages"][name]["candidate_regime"] = True
    primary = assessments.get(512)
    other = assessments.get(256)
    fine_used = 0.0
    fine_origin = CLOCK()
    disagreement = False
    if primary and other:
        disagreement = bool(primary["renewed_structure"] != other["renewed_structure"]
                            or primary["maintained_structure"] != other["maintained_structure"])
    if disagreement and primary:
        fine_elapsed = lambda: CLOCK() - fine_origin
        declare(
            record, arrays,
            "nf512 and nf256 disagree on maintenance or renewal. One half-step continuation can change that conclusion.",
            estimate_cpu(512, DURATION, DT / 2),
            "fine",
        )
        grid, _state, _metadata = loaded[512]
        final = load_episode_final(episode.NPZ, RUN_NAME[512])
        fine = integrate(
            grid, final, DURATION, DT / 2, FINE_BUDGET, fine_elapsed,
            t0=JOIN, frozen=False, stop_on_delocalization=True, label="continuation_nf512_half_dt",
        )
        fine_used = fine["cpu_seconds"]
        store_series(arrays, "continuation_nf512_half_dt", fine)
        record["stages"]["continuation_nf512_half_dt"] = public_result(fine, reference)
        record["stages"]["continuation_nf512_half_dt"]["renewal_assessment"] = continuation_assessment(512, fine)
        candidate = candidate or record["stages"]["continuation_nf512_half_dt"]["renewal_assessment"]["candidate_regime"]
    perturbations_executed = False
    if PERTURBATIONS_FROM_PROXY and candidate:
        fine_elapsed = lambda: CLOCK() - fine_origin
        declare(
            record, arrays,
            "A candidate regime was found. Plus and minus five percent imbalance perturbations re-solve the source and evolve nf512 once each.",
            2.0 * estimate_cpu(512),
            "perturbation",
        )
        grid, state, _metadata = loaded[512]
        for fraction in (0.05, -0.05):
            if fine_elapsed() + estimate_cpu(512) > FINE_BUDGET:
                record["stages"][f"perturbation_{fraction}"] = {"completed": False, "stop_reason": "FINE_BUDGET"}
                break
            occupations = imbalance_perturbation(BASELINE_OCCUPATIONS, fraction)
            owned, solved, report = prepare_population(grid, state, occupations)
            record["initial_data"][f"perturbation_{fraction}"] = report
            if not report["evolve"]:
                record["stages"][f"perturbation_{fraction}"] = {
                    "completed": False,
                    "stop_reason": "INITIAL_DATA_NOT_EVOLVED",
                    "projected_hamilton_max": report["projected_hamilton_max"],
                    "full_hamilton_max": report["full_hamilton_max"],
                }
                continue
            perturbed = integrate(
                owned, solved, DURATION, DT, FINE_BUDGET, fine_elapsed,
                t0=0.0, frozen=False, stop_on_delocalization=False,
                label=f"perturbation_{fraction}",
            )
            fine_used += perturbed["cpu_seconds"]
            store_series(arrays, f"perturbation_{fraction}", perturbed)
            record["stages"][f"perturbation_{fraction}"] = public_result(perturbed, reference)
            perturbations_executed = True
    maintained = bool(primary and primary["maintained_structure"])
    renewed = bool(primary and primary["renewed_structure"])
    one_drift = bool(primary and primary["one_drift"])
    if primary is None:
        statement = "The primary continuation did not produce a series. No renewal claim is made."
    elif renewed:
        statement = (
            "The joined shell series from T=0 through the second episode reverses "
            "or moves the localized leader by at least the declared 0.10 content scale, "
            "and the end state stays localized. That is renewed structure."
        )
    elif maintained and one_drift:
        statement = (
            "The second episode keeps the same localized leader, and the joined "
            "series is one drift. Splitting that drift at T=0.05 is not renewal."
        )
    elif maintained:
        statement = "The localized leader persists through the second episode. Renewal was not demonstrated."
    else:
        reason = None if results.get("continuation_nf512") is None else results["continuation_nf512"]["stop_reason"]
        statement = (
            "Maintained localized structure was not demonstrated through T=0.10. "
            f"Stop reason: {reason}."
        )
    record["conclusion"] = {
        "maintained_structure": maintained,
        "renewed_structure": renewed,
        "one_drift": one_drift,
        "candidate_regime": bool(candidate),
        "perturbations_executed": bool(perturbations_executed),
        "splitting_one_drift_is_renewal": False,
        "resolution_disagreement": bool(disagreement),
        "fine_cpu_seconds": float(fine_used),
        "statement": statement,
        "saved_field_exchange_scale": reference,
    }
    record["renewal"] = renewed
    record["inputs_after"] = input_hashes()
    record["inputs_preserved"] = record["inputs_before"] == record["inputs_after"]
    record["cpu_seconds"] = float(CLOCK())
    record["wall_note"] = "cpu_seconds is process time and includes preparation"
    record["process_cpu_this_run"] = float(time.process_time() - started)
    record["status"] = "MEASURED" if primary is not None else "PARTIAL"
    render_note(record)
    checkpoint(record, arrays)
    print("STATUS", record["status"], "cpu", record["cpu_seconds"], flush=True)
    print(statement, flush=True)
    return record


def assess_successor(write=False):
    """Read-only successor. v1 JSON and NPZ stay byte-identical."""
    if not OUT.is_file() or not NPZ.is_file():
        raise FileNotFoundError("v1 regeneration record is missing")
    json_hash = episode.sha256(OUT)
    npz_hash = episode.sha256(NPZ)
    record = json.loads(OUT.read_text())
    if record.get("schema") != SCHEMA or record.get("status") != "MEASURED":
        raise AssertionError("v1 is not the measured immutable record")
    fine = record["stages"]["continuation_nf512"]
    coarse = record["stages"]["continuation_nf256"]
    uniform = record["stages"]["uniform_nf512"]["structure"]
    with np.load(NPZ, allow_pickle=False) as data:
        content = np.array(data["continuation_nf512_content"], dtype=float)
        flux = np.array(data["continuation_nf512_flux"], dtype=float)
        times = np.array(data["continuation_nf512_time"], dtype=float)
        q_min = np.array(data["continuation_nf512_Q_min"], dtype=float)
        g_min = np.array(data["continuation_nf512_G_min"], dtype=float)
        chi_max = np.array(data["continuation_nf512_chi_max"], dtype=float)
        field = np.array(data["continuation_nf512_field_energy"], dtype=float)
        energy = np.array(data["continuation_nf512_energy"], dtype=float)
        unitarity = np.array(data["continuation_nf512_unitarity"], dtype=float)
        final = CauchyState(
            Q=np.array(data["continuation_nf512_final_Q"], dtype=float, copy=True),
            r=np.array(data["continuation_nf512_final_r"], dtype=float, copy=True),
            chi=np.array(data["continuation_nf512_final_chi"], dtype=float, copy=True),
            p_Q=np.array(data["continuation_nf512_final_p_Q"], dtype=float, copy=True),
            p_r=np.array(data["continuation_nf512_final_p_r"], dtype=float, copy=True),
            p_chi=np.array(data["continuation_nf512_final_p_chi"], dtype=float, copy=True),
            phi0=np.array(data["continuation_nf512_final_phi0"], copy=True),
            phi1=np.array(data["continuation_nf512_final_phi1"], copy=True),
        )
    if episode.sha256(OUT) != json_hash or episode.sha256(NPZ) != npz_hash:
        raise AssertionError("v1 bytes changed while the successor was reading them")
    grid = galerkin.build_grid(512)
    omega, quadrature_omega = galerkin.subspace_frequency(grid, final)
    stable_step, omega_again, _quadrature_again = galerkin.stable_timestep(grid, final)
    prolonged = galerkin.prolong_state(grid, final)
    failure = chart_failure(grid.fine, prolonged)
    gram = car_report(final)
    lapse_peak = float(np.max(grid.fine.length_density))
    shift_peak = float(np.max(np.abs(grid.fine.shift)))
    wavenumber = float(np.max(np.abs(grid.modes_f)) * 2.0 * np.pi / grid.length)
    mass_peak = float(np.max(np.abs(grid.fine.kappa * grid.fine.length_density)))
    omega_hat = (lapse_peak / q_min + shift_peak) * wavenumber + mass_peak
    dt = float(fine["dt"])
    cap = 1.4
    crossed = np.flatnonzero(dt * omega_hat >= cap)
    index = None if crossed.size == 0 else int(crossed[0])
    accounting = {
        "field_energy_change": float(field[-1] - field[0]),
        "energy_change": float(energy[-1] - energy[0]),
        "balance_over_exchange": fine["balance"]["versus_field_change"]["balance_over_exchange"],
        "applied_coordinate_work": fine["balance"]["applied_coordinate_work"],
        "reported_with_windows": False,
    }
    ledgers = separate_ledgers(content, flux, accounting)
    if index is None:
        crossing_ledgers = None
    else:
        crossing_ledgers = separate_ledgers(
            content[: index + 1], flux[: index + 1],
            {
                "field_energy_change": float(field[index] - field[0]),
                "energy_change": float(energy[index] - energy[0]),
                "reported_with_windows": False,
            },
        )
    algebra = galerkin_half_identity_report(galerkin.build_grid(14, quadrature=64))
    step_cpu = float(fine["cpu_seconds"] / fine["steps_completed"])
    successor = {
        "schema": SCHEMA_V2,
        "status": "ASSESSED",
        "review": REVIEW,
        "dynamics_rerun": False,
        "new_action": False,
        "new_forces": False,
        "v1_overwritten": False,
        "predecessor": {
            "json": "lab/results/development/nsc-regeneration-controls-v1.json",
            "npz": "lab/results/development/nsc-regeneration-controls-v1.npz",
            "json_sha256": json_hash,
            "npz_sha256": npz_hash,
            "schema": record["schema"],
            "immutable": True,
        },
        "preserved": {
            "attained_time": fine["attained_time"],
            "Q_min_end": fine["Q_min_end"],
            "G_min_end": fine["G_min_end"],
            "chi_max_start": fine["checkpoints"][0]["chi_max"],
            "chi_max_end": fine["checkpoints"][-1]["chi_max"],
            "leader_content_change": fine["leader_content_change"],
            "packet_flux_end": fine["packet_flux_end"],
            "reservoir_flux_end": fine["reservoir_flux_end"],
            "flux_sum_end": fine["structure"]["flux_sum_end"],
            "field_energy_change": fine["field_energy_change"],
            "balance_over_exchange": fine["balance"]["versus_field_change"]["balance_over_exchange"],
            "stop_reason": fine["stop_reason"],
            "chart_stop": fine["chart_stop"],
            "positive_G_diagnostic": fine["positive_G"],
            "dt": fine["dt"],
            "cpu_seconds": fine["cpu_seconds"],
            "steps_completed": fine["steps_completed"],
            "nf256_Q_min_end": coarse["Q_min_end"],
            "nf256_field_energy_change": coarse["field_energy_change"],
            "nf256_leader_content_change": coarse["leader_content_change"],
            "uniform_start_share": uniform["start_share"],
            "uniform_end_share": uniform["end_share"],
        },
        "ledgers": ledgers,
        "uniform_six": {
            **uniform_complement_localization(),
            "window_share_start": uniform["start_share"],
            "window_share_end": uniform["end_share"],
            "window_proxy_localized": uniform["localized_end"],
        },
        "galerkin_half_identity": algebra,
        "stability": {
            "endpoint_omega": float(omega),
            "endpoint_omega_repeat": float(omega_again),
            "quadrature_omega": float(quadrature_omega),
            "owned_stable_timestep": float(stable_step),
            "record_dt": dt,
            "dt_times_omega": float(dt * omega),
            "owned_cap": cap,
            "absolute_rk4_limit": float(2.0 * np.sqrt(2.0)),
            "above_owned_cap": bool(dt * omega > cap),
            "below_absolute_rk4_limit": bool(dt * omega < 2.0 * np.sqrt(2.0)),
            "chart_failure": failure,
            "gram_gap": gram["gap"],
            "gram_admissible": gram["admissible"],
            "positive_r": bool(float(np.min(prolonged.r)) > 0.0),
            "positive_Q": bool(float(np.min(prolonged.Q)) > 0.0),
            "positive_L": bool(float(np.min(grid.fine.length_density)) > 0.0),
            "G_min_end": float(g_min[-1]),
            "series_diagnostic": "L_peak/Q_min with the gauge lapse peak and shift held fixed",
            "series_endpoint_relative_gap": float(abs(omega_hat[-1] - omega) / omega),
            "cap_crossing_index": index,
            "cap_crossing_time": None if index is None else float(times[index]),
            "cap_crossing_Q_min": None if index is None else float(q_min[index]),
            "cap_crossing_ledgers": crossing_ledgers,
            "unitarity_max": float(np.max(unitarity)),
            "chi_max_end": float(chi_max[-1]),
        },
        "next_event": {
            "executed": False,
            "name": "one owned RK4 step from the stored T=0.10 state at stable_timestep",
            "dt": float(stable_step),
            "measured_cpu_seconds_per_recorded_step": step_cpu,
            "estimated_cpu_seconds": step_cpu,
            "full_slice_estimated_cpu_seconds": float(estimate_cpu(512)),
            "why": (
                "The stored T=0.10 state remains on the r, Q, L chart with the column "
                "Gram inside 1e-8. dt times the owned subspace frequency is already above "
                "1.4 and below 2√2. The lapse-peak reconstruction crosses 1.4 at a stored "
                "sample whose content, flux, reversal, and energy exchange keep the same "
                "reading as T=0.10. One step at the owned stable timestep is the smallest "
                "dynamical probe that can change a ledger. It is not run here."
            ),
        },
        "policy": {
            "g_is_initial_convex_radius_hypothesis": True,
            "initial_domain": "constant Q, chi = p_r = p_chi = 0, r = y^2",
            "g_is_dynamic_chart_stop": False,
            "dynamic_chart": "chart_failure on r, Q, and L",
            "candidate_regime_is_programme_requirement": False,
            "maintained_with_throughflow_allowed": True,
            "splitting_one_drift_is_complete_renewal": False,
            "fine_identity_is_algebraic_control": True,
        },
    }
    successor = episode.jsonable(successor)
    if write:
        episode.write_json(OUT_V2, successor)
        if episode.sha256(OUT) != json_hash or episode.sha256(NPZ) != npz_hash:
            raise AssertionError("writing the successor changed v1")
    return successor


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if "--assess" in arguments:
        successor = assess_successor(write=True)
        print(json.dumps({
            "schema": successor["schema"],
            "json_sha256": successor["predecessor"]["json_sha256"],
            "npz_sha256": successor["predecessor"]["npz_sha256"],
            "attained_time": successor["preserved"]["attained_time"],
            "Q_min_end": successor["preserved"]["Q_min_end"],
            "dt_times_omega": successor["stability"]["dt_times_omega"],
            "next_event_executed": successor["next_event"]["executed"],
        }, indent=2))
        return 0
    if "--check" in arguments:
        record = verify_saved()
        print(json.dumps({
            "status": record["status"],
            "cpu_seconds": record["cpu_seconds"],
            "maintained_structure": record["conclusion"]["maintained_structure"],
            "renewed_structure": record["conclusion"]["renewed_structure"],
            "one_drift": record["conclusion"]["one_drift"],
        }, indent=2))
        return 0
    campaign()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
