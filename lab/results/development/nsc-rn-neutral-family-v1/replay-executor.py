#!/usr/bin/env python3
"""Execute the neutral RN family with the frozen parent APIs.

This process owns the campaign directory, the prepare marker, and the six
case processes. It does not replace the PG step, the source, or the writer.
"""
from __future__ import annotations

import os

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"
os.environ["PYTHONUNBUFFERED"] = "1"

import gc
import json
import math
import platform
import resource
import sys
import time
import traceback
from multiprocessing import get_context
from pathlib import Path

import numpy as np

import derive_nsc_rn_parent as driver
from recursive_horizons import nsc_rn_observables as observables
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_source as source

COMMIT = "ab03859ae553d197f99430146985598fdeee3382"
CAMPAIGN = driver.LAB / "results" / "development" / "nsc-rn-neutral-family-v1"
CPU_BUDGET = float(driver.CPU_BUDGET_SECONDS)
WALL_BUDGET = 7200.0
MEMORY_BUDGET = int(driver.MEMORY_BYTES)
SNAPSHOT_OVER_RM = 0.5
SCALAR_STEPS = 10
COORDINATE_STATIONS = (5.0, 10.0, 20.0)
TARGETS = (20.0, 10.0, 5.0)
RECORDED_SPEEDUPS = (
    {"points": 869, "speedup": 2.055},
    {"points": 1750, "speedup": 10.744},
    {"points": 3921, "speedup": 30.308},
)
THREAD_NAMES = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
)


def emit(payload):
    print(json.dumps(payload, sort_keys=True, default=str), flush=True)


def to_jsonable(value):
    if isinstance(value, np.ndarray):
        array = np.asarray(value)
        if array.size > 8:
            summary = {
                "shape": [int(item) for item in array.shape],
                "dtype": str(array.dtype),
                "sha256": driver.array_sha256(array),
            }
            if not np.iscomplexobj(array) and array.size and np.isfinite(array).all():
                summary["max_abs"] = float(np.max(np.abs(array)))
            return summary
        return driver._plain(array)
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def exclusive_bytes(path, data):
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o444)
    with os.fdopen(descriptor, "wb") as handle:
        handle.write(data)
    return len(data)


def git_head():
    import subprocess
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=driver.ROOT, capture_output=True, text=True, check=True, timeout=5,
    )
    return result.stdout.strip()


def rss_tree_bytes(root_pid):
    import subprocess
    result = subprocess.run(
        ["ps", "-ax", "-o", "pid=,ppid=,rss="],
        capture_output=True, text=True, check=True, timeout=5,
    )
    rows = []
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        rows.append((int(parts[0]), int(parts[1]), int(parts[2]) * 1024))
    children = {pid: [] for pid, _parent, _rss in rows}
    own = {}
    for pid, parent, rss in rows:
        own[pid] = rss
        children.setdefault(parent, []).append(pid)
    total = 0
    stack = [root_pid]
    seen = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += own.get(pid, 0)
        stack.extend(children.get(pid, ()))
    return total


def bulk_mask(size):
    mask = np.ones(size, dtype=bool)
    if size > 20:
        mask[:10] = False
        mask[-10:] = False
    return mask


def curvature_margins(radius, r4, ricci, vacuum):
    if r4 is None or ricci is None:
        return {"available": False}
    scalar = np.asarray(r4, dtype=float)
    square = np.asarray(ricci, dtype=float)
    mask = bulk_mask(scalar.size)
    departure = np.abs(square - np.asarray(vacuum, dtype=float))
    return {
        "available": True,
        "R4_full_max_abs": float(np.max(np.abs(scalar))),
        "R4_bulk_10_max_abs": float(np.max(np.abs(scalar[mask]))),
        "Ricci2_full_max": float(np.max(square)),
        "Ricci2_bulk_10_max": float(np.max(square[mask])),
        "Ricci2_vacuum_departure_full_max": float(np.max(departure)),
        "Ricci2_vacuum_departure_bulk_10_max": float(np.max(departure[mask])),
        "nodes": int(scalar.size),
        "bulk_nodes_excluded_each_end": 10,
    }


def stress_report(phi, occupations, lapse, shift, radius, weights, derivative):
    evaluated = source.source_evaluate(
        phi, occupations, lapse, shift, np.ones_like(radius), radius, weights, derivative,
    )
    rho = np.asarray(evaluated["physical"]["rho"], dtype=float)
    radial = np.asarray(evaluated["physical"]["radial_pressure"], dtype=float)
    angular = np.asarray(evaluated["physical"]["angular_pressure"], dtype=float)
    radial_null = rho + radial
    angular_null = rho + angular
    trace = np.asarray(evaluated["massless_trace_residual"], dtype=float)
    return evaluated, {
        "available": True,
        "source": "nsc_rn_source.source_evaluate physical stress",
        "null_combination": "rho+pressure from the existing orthonormal stress",
        "null_energy_theorem": False,
        "rho_min": float(np.min(rho)),
        "rho_max": float(np.max(rho)),
        "radial_pressure_min": float(np.min(radial)),
        "angular_pressure_min": float(np.min(angular)),
        "radial_null_min": float(np.min(radial_null)),
        "angular_null_min": float(np.min(angular_null)),
        "negative_radial_null_fraction": float(np.mean(radial_null < 0.0)),
        "negative_angular_null_fraction": float(np.mean(angular_null < 0.0)),
        "massless_trace_max_abs": float(np.max(np.abs(trace))),
    }, rho, radial, angular


def mean_radius(phi, radius, weights):
    density = np.sum(np.abs(np.asarray(phi)) ** 2, axis=0)
    total = float(np.sum(weights * density))
    if total <= 0.0:
        return None
    return float(np.sum(weights * radius * density) / total)


def fixed_diag(state, grid, metric, initial_norm):
    phi_rate, _sat = pg.dirac_rates(
        grid, state["phi"], metric["lapse"], metric["shift"], kappa=state["kappa"],
    )
    jets = pg.metric_jets(grid["radius"], metric["lapse"], metric["shift"])
    zeros = np.zeros(grid["points"])
    account = source.energy_account(
        state["phi"], state["occupations"], metric["lapse"], metric["shift"],
        grid["radius"], grid["weights"], grid["derivative"],
        phi_rate=phi_rate, lapse_t=zeros, shift_t=zeros,
        lapse_r=jets["lapse_r"], shift_r=jets["shift_r"],
    )
    car = pg.current_occupation_car(state["phi"], grid["weights"], state["occupations"])
    located = observables.trapping_horizon(grid["radius"], metric["lapse"], metric["shift"])
    ledger = driver.sat_norm_ledger(
        state["stocks"]["sat_debit"], state["probability_stocks"], initial_norm,
    )
    return {
        "probability_stocks": dict(state["probability_stocks"]),
        "inner_mass": float(state["inner_mass"]),
        "clocks": dict(state["clocks"]),
        "sat_norm_stock": ledger,
        "E_coordinate": account["E_coordinate"],
        "E_normal": account["E_normal"],
        "metric_power": account["metric_power"],
        "trapping_horizon": located["trapping_horizon"],
        "car_min_eigenvalue": car["min_eigenvalue"],
        "car_max_eigenvalue": car["max_eigenvalue"],
        "car_psd": car["psd"],
        "occupations_renormalized": False,
        "metric_reconstructed": False,
        "control": driver.FIXED_CONTROL,
        "second_einstein_solution": False,
    }


def coupled_snapshot(coupled, grid, rates, packet, initial_norm):
    observation = driver.station_observation(
        coupled, grid, rates, packet=packet, initial_norm=initial_norm,
    )
    vacuum = np.asarray(grid["model"].invariants(grid["radius"])["Ricci2"], dtype=float)
    margins = curvature_margins(
        grid["radius"], observation.get("R4"), observation.get("Ricci2"), vacuum,
    )
    stress, stress_summary, rho, radial, angular = stress_report(
        coupled["phi"], coupled["occupations"], rates["lapse"], rates["shift"],
        grid["radius"], grid["weights"], grid["derivative"],
    )
    return observation, margins, stress_summary, rho, radial, angular, vacuum


def compact_sample(t_over, coupled_obs, fixed, margins, stress, defects, balance_max):
    energy = coupled_obs.get("energy") or {}
    car = (coupled_obs.get("car") or {}).get("current") or {}
    stocks = coupled_obs.get("probability_stocks") or {}
    fixed_stocks = fixed.get("probability_stocks") or {}
    trapping = coupled_obs.get("trapping_horizon")
    fixed_trapping = fixed.get("trapping_horizon")
    closure = (coupled_obs.get("sat_norm_stock") or {}).get("probability_closure_residual")
    fixed_closure = (fixed.get("sat_norm_stock") or {}).get("probability_closure_residual")

    def diff(left, right):
        if left is None or right is None:
            return None
        return float(left) - float(right)

    return {
        "t_over_rm": float(t_over),
        "coupled_remaining": stocks.get("remaining"),
        "coupled_excision": stocks.get("excision_outflow"),
        "coupled_outer": stocks.get("outer_outflow"),
        "fixed_remaining": fixed_stocks.get("remaining"),
        "fixed_excision": fixed_stocks.get("excision_outflow"),
        "fixed_outer": fixed_stocks.get("outer_outflow"),
        "coupled_inner_mass": coupled_obs.get("charged_mass_inner"),
        "fixed_inner_mass": fixed.get("inner_mass"),
        "coupled_trapping": trapping,
        "fixed_trapping": fixed_trapping,
        "coupled_E_coordinate": energy.get("E_coordinate"),
        "fixed_E_coordinate": fixed.get("E_coordinate"),
        "coupled_E_normal": energy.get("E_normal"),
        "fixed_E_normal": fixed.get("E_normal"),
        "metric_power": energy.get("metric_power"),
        "fixed_metric_power": fixed.get("metric_power"),
        "car_min": car.get("min_eigenvalue"),
        "car_psd": car.get("psd"),
        "fixed_car_min": fixed.get("car_min_eigenvalue"),
        "fixed_car_psd": fixed.get("car_psd"),
        "closure_residual": closure,
        "fixed_closure_residual": fixed_closure,
        "mass_slope_residual": energy.get("mass_slope_residual"),
        "normal_continuity_residual": energy.get("normal_continuity_residual"),
        "balance_residual_max_abs": balance_max,
        "projection_defects": defects,
        "curvature": margins,
        "stress": stress,
        "delta_remaining": diff(stocks.get("remaining"), fixed_stocks.get("remaining")),
        "delta_excision": diff(stocks.get("excision_outflow"), fixed_stocks.get("excision_outflow")),
        "delta_outer": diff(stocks.get("outer_outflow"), fixed_stocks.get("outer_outflow")),
        "delta_inner_mass": diff(coupled_obs.get("charged_mass_inner"), fixed.get("inner_mass")),
        "delta_E_coordinate": diff(energy.get("E_coordinate"), fixed.get("E_coordinate")),
        "delta_trapping": diff(trapping, fixed_trapping),
        "clock_t": coupled_obs.get("clock_t"),
        "fixed_clocks": fixed.get("clocks"),
    }


def events_from_series(series, initial_norm):
    t = np.asarray(series["t_over_rm"], dtype=float)
    excision = np.asarray(series["coupled_excision"], dtype=float)
    outer = np.asarray(series["coupled_outer"], dtype=float)
    events = []
    if t.size < 2 or initial_norm <= 0.0:
        return events
    capture_rate = np.diff(excision)
    reflect_rate = np.diff(outer)
    capture_index = int(np.argmax(capture_rate))
    reflect_index = int(np.argmax(reflect_rate))
    events.append({
        "name": "peak_capture_rate",
        "t_over_rm_left": float(t[capture_index]),
        "t_over_rm_right": float(t[capture_index + 1]),
        "delta_excision": float(capture_rate[capture_index]),
    })
    events.append({
        "name": "peak_reflection_rate",
        "t_over_rm_left": float(t[reflect_index]),
        "t_over_rm_right": float(t[reflect_index + 1]),
        "delta_outer": float(reflect_rate[reflect_index]),
    })
    for level, name, values in (
        (0.5, "capture_fraction_0.5", excision),
        (0.01, "outer_fraction_0.01", outer),
    ):
        hits = np.nonzero(values >= level * float(initial_norm))[0]
        if hits.size:
            events.append({
                "name": name,
                "occurred": True,
                "t_over_rm": float(t[int(hits[0])]),
                "value": float(values[int(hits[0])]),
                "threshold_times_initial_norm": level,
            })
        else:
            events.append({
                "name": name,
                "occurred": False,
                "threshold_times_initial_norm": level,
            })
    return events


def packet_gate(packet, grid):
    admitted = driver.admit_prepared_packet(packet)
    frequencies = [float(item) for item in packet["spectrum"]["requested"]]
    reasons = []
    if len(frequencies) != 9 or any(item <= 0.0 for item in frequencies):
        reasons.append("packet is not nine positive exterior frequencies")
    if admitted["car_rank"] != 1:
        reasons.append("prepared CAR rank is not 1")
    if admitted["inward_at_center"] is not True:
        reasons.append("packet current is not inward at the centre")
    phi = np.asarray(packet["phi"])
    if phi.shape != (2, grid["points"], 1):
        reasons.append("phi is not the rank-one column on this grid")
    if packet.get("absorbing_matrix_projection") is not False:
        reasons.append("absorbing matrix projection is not the real constraint path")
    if "gram" not in packet.get("car", {}):
        reasons.append("weighted Gram is missing")
    occupations = np.asarray(packet["occupations"], dtype=float)
    if occupations.shape != (1,) or not (0.0 < float(occupations[0]) < 1.0):
        reasons.append("filling left the open interval (0, 1)")
    return admitted, frequencies, reasons


def backend_record(grid):
    import scipy
    return {
        "derivative_layout": grid["derivative_layout"],
        "derivative_family": grid["derivative_family"],
        "production_operator": bool(grid["production_operator"]),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "python": sys.version,
        "platform": platform.platform(),
        "threads": {name: os.environ.get(name) for name in THREAD_NAMES},
    }


def dense_reference(mu):
    geometry = driver.case_geometry(mu)
    points = int(geometry["route"]["points"])
    outer = float(geometry["route"]["outer"])
    mass = float(geometry["mass"])
    r_m = float(geometry["r_m"])
    started = time.process_time()
    dense = pg.build_grid(mass, r_m=r_m, points=points, r_out=outer, layout="dense")
    csr = pg.build_grid(mass, r_m=r_m, points=points, r_out=outer, layout="csr")
    coefficient_gap = float(np.max(np.abs(
        csr["derivative"].toarray() - np.asarray(dense["derivative"], dtype=float)
    )))
    phi = np.zeros((2, points, 1), dtype=complex)
    phi[0, :, 0] = 1.0e-4
    phi[1, :, 0] = 5.0e-5
    occupation = np.array([1.0e-4])
    probe = pg.make_state(csr, phi, occupations=occupation, inner_mass=mass, killing_frequency=1.0)
    rates = pg.stage_rates(probe, csr)
    metric = driver.freeze_sourced_metric(rates)
    step = driver.matched_step(rates["cfl_dt"], r_m)
    dense_rates = pg.stage_rates(
        pg.make_state(dense, phi, occupations=occupation, inner_mass=mass, killing_frequency=1.0),
        dense,
    )
    rate_gap = float(np.max(np.abs(rates["phi_rate"] - dense_rates["phi_rate"])))
    lapse_gap = float(np.max(np.abs(rates["lapse"] - dense_rates["lapse"])))
    shift_gap = float(np.max(np.abs(rates["shift"] - dense_rates["shift"])))

    def once(grid):
        coupled = pg.make_state(grid, phi, occupations=occupation, inner_mass=mass, killing_frequency=1.0)
        fixed = driver._clone_state(coupled)
        clock = time.process_time()
        marched = driver.march_pair(coupled, fixed, grid, metric, dt=step["dt"], steps=1)
        return marched, time.process_time() - clock

    dense_march, dense_cpu = once(dense)
    csr_march, csr_cpu = once(csr)
    phi_gap = float(np.max(np.abs(dense_march["coupled"]["phi"] - csr_march["coupled"]["phi"])))
    fixed_gap = float(np.max(np.abs(dense_march["fixed"]["phi"] - csr_march["fixed"]["phi"])))
    elapsed = time.process_time() - started
    del dense, dense_rates, dense_march
    gc.collect()
    return {
        "mass_over_rm": float(mu),
        "points": points,
        "coefficient_max_abs": coefficient_gap,
        "source_phi_rate_max_abs": rate_gap,
        "lapse_max_abs": lapse_gap,
        "shift_max_abs": shift_gap,
        "rk4_phi_max_abs": phi_gap,
        "fixed_rk4_phi_max_abs": fixed_gap,
        "dense_pair_cpu_seconds": dense_cpu,
        "csr_pair_cpu_seconds": csr_cpu,
        "measured_pair_speedup": None if csr_cpu <= 0.0 else dense_cpu / csr_cpu,
        "measurement_cpu_seconds": elapsed,
        "independent_of_saved_speedup_table": True,
    }


class CpuMeter:
    def __init__(self, shared):
        self.shared = shared
        self.last = time.process_time()

    def add(self):
        now = time.process_time()
        delta = now - self.last
        self.last = now
        with self.shared.get_lock():
            self.shared.value += float(delta)
        return float(delta)

    def used(self):
        with self.shared.get_lock():
            return float(self.shared.value)


def append_sample(series, coupled, fixed, radius, weights, r_m):
    series["t_over_rm"].append(float(coupled["clocks"]["t"]) / float(r_m))
    series["coupled_remaining"].append(float(coupled["probability_stocks"]["remaining"]))
    series["coupled_excision"].append(float(coupled["probability_stocks"]["excision_outflow"]))
    series["coupled_outer"].append(float(coupled["probability_stocks"]["outer_outflow"]))
    series["coupled_inner_mass"].append(float(coupled["inner_mass"]))
    series["coupled_sat"].append(float(coupled["stocks"]["sat_debit"]))
    series["coupled_killing"].append(float(coupled["clocks"]["killing"]))
    series["coupled_normal"].append(float(coupled["clocks"]["normal"]))
    series["fixed_remaining"].append(float(fixed["probability_stocks"]["remaining"]))
    series["fixed_excision"].append(float(fixed["probability_stocks"]["excision_outflow"]))
    series["fixed_outer"].append(float(fixed["probability_stocks"]["outer_outflow"]))
    series["fixed_inner_mass"].append(float(fixed["inner_mass"]))
    series["fixed_sat"].append(float(fixed["stocks"]["sat_debit"]))
    series["fixed_killing"].append(float(fixed["clocks"]["killing"]))
    series["fixed_normal"].append(float(fixed["clocks"]["normal"]))
    centre = mean_radius(coupled["phi"], radius, weights)
    series["coupled_mean_radius"].append(None if centre is None else float(centre))


def empty_series():
    return {key: [] for key in (
        "t_over_rm", "coupled_remaining", "coupled_excision", "coupled_outer",
        "coupled_inner_mass", "coupled_sat", "coupled_killing", "coupled_normal",
        "fixed_remaining", "fixed_excision", "fixed_outer", "fixed_inner_mass",
        "fixed_sat", "fixed_killing", "fixed_normal", "coupled_mean_radius",
    )}


def slice_series(series, start):
    arrays = {}
    for key, values in series.items():
        chunk = values[start:]
        if key == "coupled_mean_radius":
            arrays["scalar_" + key] = np.array(
                [math.nan if item is None else float(item) for item in chunk], dtype=float,
            )
        else:
            arrays["scalar_" + key] = np.asarray(chunk, dtype=float)
    return arrays


def stage_name(t_over):
    milli = int(round(float(t_over) * 10.0))
    if milli % 10 == 0:
        return f"station-{milli // 10}"
    return f"station-{milli // 10}.{milli % 10}"


def write_failure(directory, commit, reason, coupled, fixed, metric, extra):
    record = {
        "producing_commit": commit,
        "blocker": reason,
        "fixed_control": driver.FIXED_CONTROL,
        "second_einstein_solution": False,
        "marched_by_api": True,
        "clock_t": None if coupled is None else coupled["clocks"]["t"],
        **extra,
    }
    arrays = None if coupled is None or fixed is None or metric is None else driver.state_arrays(
        coupled, fixed, metric,
    )
    try:
        return driver.write_stage(directory, "station-failure", to_jsonable(record), arrays, blob_match=True)
    except driver.CampaignBlocker as blocked:
        return {"blocker": blocked.payload.get("blocker"), "bytes_written": 0}


def observe_and_write(directory, coupled, fixed, grid, metric, packet, initial_norm,
                      initial_occupations, series, series_start, t_over, step, commit,
                      mu, eps, include_packet):
    rates = pg.stage_rates(coupled, grid)
    observation, margins, stress_summary, rho, radial, angular, _vacuum = coupled_snapshot(
        coupled, grid, rates, packet if include_packet else None, initial_norm,
    )
    fixed = fixed_diag(fixed, grid, metric, initial_norm)
    defects = to_jsonable(rates["dynamic_contract"])
    balance = float(np.max(np.abs(rates["balance_residual"])))
    sample = compact_sample(
        t_over, observation, fixed, margins, stress_summary, defects, balance,
    )
    sample["occupations_unchanged"] = bool(np.array_equal(coupled["occupations"], initial_occupations))
    sample["phi_sha256"] = driver.array_sha256(coupled["phi"])
    sample["fixed_phi_sha256"] = driver.array_sha256(fixed_phi_holder(fixed))
    events = events_from_series(series, initial_norm)
    record = {
        "station_over_rm": float(t_over),
        "mass_over_rm": float(mu),
        "energy_fraction": float(eps),
        "fixed_control": driver.FIXED_CONTROL,
        "second_einstein_solution": False,
        "producing_commit": commit,
        "step": step,
        "backend": backend_record(grid),
        "marched_by_api": True,
        "observation": to_jsonable(observation),
        "fixed_diagnostic": to_jsonable(fixed),
        "sample": sample,
        "events_so_far": events,
        "current_car_recomputed": True,
        "cached_car_not_used_as_current": True,
        "phi_not_renormalized": True,
    }
    arrays = driver.state_arrays(coupled, fixed_state_from_diag(fixed, metric, coupled), metric)
    # The helper above is wrong if I pass a diagnostic dict. Rebuild from the live fixed state
    # before this call. This placeholder is replaced by the caller-supplied arrays path.
    raise RuntimeError("internal snapshot writer requires the live fixed state")


def fixed_phi_holder(state):
    return state["phi"]


def run_case(spec):
    os.environ.update({name: "1" for name in THREAD_NAMES})
    meter = CpuMeter(spec["cpu_used"])
    name = spec["name"]
    directory = Path(spec["campaign"]) / name
    commit = spec["commit"]
    mu = float(spec["mu"])
    eps = float(spec["eps"])
    result = {
        "case": name,
        "mass_over_rm": mu,
        "energy_fraction": eps,
        "ok": False,
        "evolved": False,
        "blocker": None,
        "stations_reached": [],
        "stopped": None,
    }
    coupled = fixed = metric = None
    try:
        record = json.loads(Path(spec["calibration_path"]).read_text())
        decision = driver.admission_decision(record)
        if not decision["admitted"]:
            raise driver.CampaignBlocker(decision["blocker"])
        if not driver.frozen_blob_match(commit, driver.producers()):
            raise driver.CampaignBlocker("producer bytes are not the frozen root commit")
        geometry = driver.case_geometry(mu)
        grid = pg.build_grid(
            geometry["mass"], r_m=geometry["r_m"], points=int(geometry["route"]["points"]),
            r_out=float(geometry["route"]["outer"]),
        )
        if grid["derivative_layout"] != "csr" or grid["production_operator"] is not True:
            raise driver.CampaignBlocker("production grid is not the CSR D(2,1) operator")
        prep_clock = time.process_time()
        packet = source.prepare_radial_packet(
            mass_over_rm=mu, energy_target=eps, rank=driver.PILOT_RANK, model=grid["model"],
            radii=grid["radius"], weights=grid["weights"], derivative=grid["derivative"],
        )
        prep_seconds = time.process_time() - prep_clock
        admitted, frequencies, reasons = packet_gate(packet, grid)
        if reasons:
            raise driver.CampaignBlocker("; ".join(reasons))
        phi = np.asarray(packet["phi"])
        coupled = pg.make_state(
            grid, phi, inner_mass=geometry["mass"],
            killing_frequency=float(packet["radial_moments"]["requested_omega"]),
            occupations=np.asarray(packet["occupations"], dtype=float),
        )
        initial_rates = pg.stage_rates(coupled, grid)
        metric = driver.freeze_sourced_metric(initial_rates)
        fixed = driver._clone_state(coupled)
        if not np.array_equal(coupled["phi"], fixed["phi"]):
            raise driver.CampaignBlocker("fixed control did not start from the same phi")
        if coupled["occupations"].tolist() != fixed["occupations"].tolist():
            raise driver.CampaignBlocker("fixed control did not start from the same filling")
        initial_norm = float(coupled["probability_stocks"]["remaining"])
        initial_occupations = np.array(coupled["occupations"], copy=True)
        step = driver.matched_step(initial_rates["cfl_dt"], geometry["r_m"])
        obs_clock = time.process_time()
        _observation, _margins, _stress, _rho, _radial, _angular, _vacuum = coupled_snapshot(
            coupled, grid, initial_rates, packet, initial_norm,
        )
        _fixed_view = fixed_diag(fixed, grid, metric, initial_norm)
        obs_seconds = time.process_time() - obs_clock
        probe_coupled = driver._clone_state(coupled)
        probe_fixed = driver._clone_state(fixed)
        step_clock = time.process_time()
        driver.march_pair(probe_coupled, probe_fixed, grid, metric, dt=step["dt"], steps=1)
        pair_seconds = time.process_time() - step_clock
        del probe_coupled, probe_fixed
        meter.add()
        forecast_body = {
            "prep_seconds": prep_seconds,
            "pair_step_seconds": pair_seconds,
            "observation_seconds": obs_seconds,
            "dt": step["dt"],
            "points": int(grid["points"]),
            "cfl_dt": float(initial_rates["cfl_dt"]),
            "limited_by": step["limited_by"],
        }
        spec["queue"].put({
            "case": name,
            "ok": True,
            "mu": mu,
            "eps": eps,
            "forecast": forecast_body,
            "cpu_seconds": meter.used(),
        })
        remaining_wait = max(1.0, float(spec["deadline"]) - time.time())
        if not spec["released"].wait(remaining_wait):
            written = write_failure(
                directory, commit, "admission wait expired before the march",
                coupled, fixed, metric, {"forecast": forecast_body},
            )
            result.update({
                "blocker": "admission wait expired before the march",
                "stopped": "admission_wait",
                "bytes_written": written.get("bytes_written", 0),
                "forecast": forecast_body,
            })
            return result
        choice = dict(spec["decision"].get(name, {"admitted": False, "blocker": "no admission record"}))
        if not choice.get("admitted"):
            written = write_failure(
                directory, commit, choice.get("blocker") or "timing admission refused the case",
                coupled, fixed, metric, {"forecast": forecast_body, "choice": choice},
            )
            result.update({
                "blocker": choice.get("blocker") or "timing admission refused the case",
                "stopped": "admission",
                "bytes_written": written.get("bytes_written", 0),
                "forecast": forecast_body,
                "packet": packet_public(packet, admitted, frequencies),
            })
            return result
        target = float(choice["target_over_rm"])
        cap = float(choice["cpu_cap_seconds"])
        march_started = time.process_time()
        series = empty_series()
        append_sample(series, coupled, fixed, grid["radius"], grid["weights"], geometry["r_m"])
        written = write_snapshot(
            directory, coupled, fixed, grid, metric, packet, initial_norm, initial_occupations,
            series, 0, 0.0, step, commit, mu, eps, True, initial_rates,
        )
        bytes_written = int(written["bytes_written"])
        stations = [0.0]
        samples = {0.0: written["sample"]}
        series_start = 0
        next_snapshot = SNAPSHOT_OVER_RM
        stopped = "station"
        blocker = None
        while coupled["clocks"]["t"] < target * geometry["r_m"] - 1e-12:
            if spec["stop_event"].is_set() or time.time() > float(spec["deadline"]) - 20.0:
                stopped = "wall_or_memory"
                blocker = "wall, memory, or parent stop reached before the next chunk"
                break
            spent_march = time.process_time() - march_started
            chunk_cap = cap - driver.ADMISSION_FACTOR * spent_march
            if chunk_cap <= 0.0:
                stopped = "cpu_budget"
                blocker = "case CPU cap reached before the next chunk"
                break
            used = meter.used()
            if driver.ADMISSION_FACTOR * (used + max(pair_seconds, 1.0e-6)) > CPU_BUDGET:
                stopped = "cpu_budget"
                blocker = "aggregate CPU projection reached the 21600 s cap"
                break
            snap_time = min(next_snapshot, target) * geometry["r_m"]
            scalar_time = coupled["clocks"]["t"] + SCALAR_STEPS * step["dt"]
            piece = min(snap_time, scalar_time)
            marched = driver.march_pair(
                coupled, fixed, grid, metric, dt=step["dt"], station_time=piece,
                cpu_seconds=chunk_cap,
            )
            coupled = marched["coupled"]
            fixed = marched["fixed"]
            meter.add()
            append_sample(series, coupled, fixed, grid["radius"], grid["weights"], geometry["r_m"])
            if marched["stopped"] == "cpu_budget":
                stopped = "cpu_budget"
                blocker = "march_pair stopped before the chunk exceeded the CPU projection"
                break
            if coupled["clocks"]["t"] + 1e-10 >= next_snapshot * geometry["r_m"] or (
                coupled["clocks"]["t"] + 1e-10 >= target * geometry["r_m"]
                and next_snapshot >= target
            ):
                landed = coupled["clocks"]["t"] / geometry["r_m"]
                snapshot = write_snapshot(
                    directory, coupled, fixed, grid, metric, packet, initial_norm,
                    initial_occupations, series, series_start, landed, step, commit, mu, eps,
                    False, None,
                )
                bytes_written += int(snapshot["bytes_written"])
                stations.append(float(landed))
                key = round(float(landed), 4)
                if any(abs(key - item) < 5.0e-3 for item in (0.0,) + COORDINATE_STATIONS) or abs(key - target) < 5.0e-3:
                    samples[key] = snapshot["sample"]
                series_start = len(series["t_over_rm"])
                emit({
                    "phase": "snapshot",
                    "case": name,
                    "t_over_rm": landed,
                    "remaining": snapshot["sample"]["coupled_remaining"],
                    "excision": snapshot["sample"]["coupled_excision"],
                    "outer": snapshot["sample"]["coupled_outer"],
                })
                next_snapshot = round(next_snapshot + SNAPSHOT_OVER_RM, 4)
                if next_snapshot > target + 1.0e-9:
                    stopped = "station"
                    break
        if stopped != "station":
            failure = write_failure(
                directory, commit, blocker or stopped, coupled, fixed, metric,
                {
                    "stopped": stopped,
                    "t_over_rm": coupled["clocks"]["t"] / geometry["r_m"],
                    "events_so_far": events_from_series(series, initial_norm),
                    "target_over_rm": target,
                },
            )
            bytes_written += int(failure.get("bytes_written", 0))
        meter.add()
        result.update({
            "ok": stopped == "station",
            "evolved": True,
            "blocker": blocker,
            "stopped": stopped,
            "stations_reached": stations,
            "target_over_rm": target,
            "cpu_cap_seconds": cap,
            "bytes_written": bytes_written,
            "forecast": forecast_body,
            "packet": packet_public(packet, admitted, frequencies),
            "geometry": geometry_public(geometry, grid),
            "samples": {str(key): value for key, value in samples.items()},
            "events": events_from_series(series, initial_norm),
            "final_t_over_rm": coupled["clocks"]["t"] / geometry["r_m"],
            "occupations_unchanged": bool(np.array_equal(coupled["occupations"], initial_occupations)),
            "max_rss_bytes": peak_rss(),
            "backend": backend_record(grid),
            "closure_max_abs": max(
                abs(float(item.get("closure_residual") or 0.0)) for item in samples.values()
            ) if samples else None,
        })
        return result
    except Exception as error:
        meter.add()
        reason = str(error)
        if isinstance(error, driver.CampaignBlocker):
            reason = error.payload.get("blocker") or reason
        written = {"bytes_written": 0}
        if directory.exists():
            written = write_failure(
                directory, commit, reason, coupled, fixed, metric,
                {"traceback": traceback.format_exc(limit=20)},
            )
        try:
            spec["queue"].put({
                "case": name, "ok": False, "mu": mu, "eps": eps,
                "blocker": reason, "cpu_seconds": meter.used(),
            })
        except Exception:
            pass
        result.update({
            "blocker": reason,
            "stopped": "exception",
            "bytes_written": written.get("bytes_written", 0),
            "traceback": traceback.format_exc(limit=20),
            "max_rss_bytes": peak_rss(),
        })
        return result


def fixed_state_from_diag(diagnostic, metric, template):
    del diagnostic, metric, template
    raise RuntimeError("unused")


def write_snapshot(directory, coupled, fixed, grid, metric, packet, initial_norm,
                   initial_occupations, series, series_start, t_over, step, commit,
                   mu, eps, include_packet, rates):
    if rates is None:
        rates = pg.stage_rates(coupled, grid)
    observation, margins, stress_summary, rho, radial, angular, _vacuum = coupled_snapshot(
        coupled, grid, rates, packet if include_packet else None, initial_norm,
    )
    fixed_view = fixed_diag(fixed, grid, metric, initial_norm)
    defects = to_jsonable(rates["dynamic_contract"])
    balance = float(np.max(np.abs(rates["balance_residual"])))
    sample = compact_sample(t_over, observation, fixed_view, margins, stress_summary, defects, balance)
    sample["occupations_unchanged"] = bool(np.array_equal(coupled["occupations"], initial_occupations))
    sample["phi_sha256"] = driver.array_sha256(coupled["phi"])
    sample["fixed_phi_sha256"] = driver.array_sha256(fixed["phi"])
    record = {
        "station_over_rm": float(t_over),
        "mass_over_rm": float(mu),
        "energy_fraction": float(eps),
        "fixed_control": driver.FIXED_CONTROL,
        "second_einstein_solution": False,
        "producing_commit": commit,
        "step": step,
        "backend": backend_record(grid),
        "marched_by_api": True,
        "observation": to_jsonable(observation),
        "fixed_diagnostic": to_jsonable(fixed_view),
        "sample": sample,
        "events_so_far": events_from_series(series, initial_norm),
        "current_car_recomputed": True,
        "cached_car_not_used_as_current": True,
        "phi_not_renormalized": sample["occupations_unchanged"],
        "G_N": float(grid["G_N"]),
        "C_F": float(grid["C_F"]),
        "P2": float(grid["magnetic_r2"]),
        "flux": float(grid["flux"]),
        "A": float(grid["A"]),
    }
    arrays = driver.state_arrays(coupled, fixed, metric)
    arrays.update({
        "radius": np.asarray(grid["radius"], dtype=float),
        "coupled_mass": np.asarray(rates["mass"], dtype=float),
        "R4": np.asarray(observation["R4"], dtype=float),
        "Ricci2": np.asarray(observation["Ricci2"], dtype=float),
        "rho": np.asarray(rho, dtype=float),
        "radial_pressure": np.asarray(radial, dtype=float),
        "angular_pressure": np.asarray(angular, dtype=float),
    })
    arrays.update(slice_series(series, series_start))
    written = driver.write_stage(
        directory, stage_name(t_over), to_jsonable(record), arrays, blob_match=True,
    )
    written["sample"] = sample
    return written


def packet_public(packet, admitted, frequencies):
    moments = packet["radial_moments"]
    return {
        "frequencies": frequencies,
        "car_rank": admitted["car_rank"],
        "inward_at_center": admitted["inward_at_center"],
        "center_probability_current": admitted["center_probability_current"],
        "mean_radius": admitted["mean_radius"],
        "width_rms": admitted["width_rms"],
        "requested_center": admitted["requested_center"],
        "requested_width": admitted["requested_width"],
        "packet_mode_energy": admitted["packet_mode_energy"],
        "killing_energy": admitted["killing_energy"],
        "filling": float(packet["packet_filling"]),
        "actual_energy_target": float(packet["actual_energy_target"]),
        "energy_target": float(packet["energy_target"]),
        "G_N": float(packet["G_N"]),
        "phi_sha256": driver.array_sha256(packet["phi"]),
        "occupation_sha256": driver.array_sha256(packet["occupations"]),
        "phi_shape": admitted["phi_shape"],
        "width_is_fit": False,
        "absorbing_matrix_projection": False,
    }


def geometry_public(geometry, grid):
    return {
        "mass": float(geometry["mass"]),
        "r_m": float(geometry["r_m"]),
        "r_minus": float(geometry["r_minus"]),
        "r_plus": float(geometry["r_plus"]),
        "excision": float(geometry["excision"]),
        "points": int(grid["points"]),
        "spacing": float(grid["spacing"]),
        "outer": float(grid["r_out"]),
        "G_N": float(grid["G_N"]),
        "C_F": float(grid["C_F"]),
        "A": float(grid["A"]),
        "P2": float(grid["magnetic_r2"]),
        "flux": float(grid["flux"]),
        "route_points_requested": int(geometry["route"]["points"]),
        "route_spacing": float(geometry["route"]["spacing"]),
    }


def peak_rss():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss if sys.platform == "darwin" else rss * 1024)


def case_name(mu, eps):
    return f"mu-{mu:.3f}-eps-{eps:.3f}"


def forecast_march(body, target):
    steps = int(math.ceil(float(target) / float(body["dt"])))
    observations = int(math.floor(float(target) / SNAPSHOT_OVER_RM))
    expected = body["pair_step_seconds"] * steps + body["observation_seconds"] * observations
    return {
        "target_over_rm": float(target),
        "steps": steps,
        "observations_after_zero": observations,
        "expected_cpu_seconds": expected,
        "forecast_cpu_seconds": driver.ADMISSION_FACTOR * expected,
    }


def choose_target(messages, cpu_already, wall_left):
    runnable = [item for item in messages if item.get("ok") and item.get("forecast")]
    refused = [item for item in messages if not item.get("ok")]
    remaining = CPU_BUDGET - cpu_already
    selected = None
    table = {}
    for target in TARGETS:
        rows = [forecast_march(item["forecast"], target) for item in runnable]
        total = sum(row["forecast_cpu_seconds"] for row in rows)
        slowest = max((row["expected_cpu_seconds"] for row in rows), default=0.0)
        table[str(target)] = {
            "forecast_cpu_seconds": total,
            "slowest_expected_seconds": slowest,
            "fits_cpu": total <= remaining,
            "fits_wall": slowest <= wall_left,
        }
        if selected is None and total <= remaining and slowest <= wall_left:
            selected = target
    caps = {}
    if selected is not None and runnable:
        rows = [forecast_march(item["forecast"], selected) for item in runnable]
        weight = sum(row["expected_cpu_seconds"] for row in rows)
        for item, row in zip(runnable, rows):
            share = remaining if weight <= 0.0 else remaining * row["expected_cpu_seconds"] / weight
            caps[item["case"]] = {
                "admitted": True,
                "target_over_rm": selected,
                "cpu_cap_seconds": share,
                "forecast": row,
                "blocker": None,
            }
    for item in runnable:
        if item["case"] not in caps:
            caps[item["case"]] = {
                "admitted": False,
                "target_over_rm": None,
                "cpu_cap_seconds": 0.0,
                "blocker": "timing admission refused the family at stations 5, 10, and 20",
                "table": table,
            }
    for item in refused:
        caps[item["case"]] = {
            "admitted": False,
            "blocker": item.get("blocker") or "case forecast failed",
            "target_over_rm": None,
            "cpu_cap_seconds": 0.0,
        }
    return selected, caps, table, remaining


def file_hashes(root):
    found = {}
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "summary.json":
            found[str(path.relative_to(root))] = driver.file_sha256(path)
    return found


def monitor(stop_event, deadline):
    root = os.getpid()
    while not stop_event.is_set():
        if time.time() > deadline:
            stop_event.set()
            emit({"phase": "stop", "reason": "wall"})
            return
        try:
            total = rss_tree_bytes(root)
        except Exception:
            total = 0
        if total > MEMORY_BUDGET:
            stop_event.set()
            emit({"phase": "stop", "reason": "memory", "rss_bytes": total})
            return
        time.sleep(2.0)


def main():
    wall_start = time.time()
    cpu_start = time.process_time()
    deadline = wall_start + WALL_BUDGET
    head = git_head()
    if not driver.frozen_blob_match(COMMIT, driver.producers()):
        emit({
            "phase": "refused",
            "blocker": "producer bytes are not the frozen root commit",
            "head": head,
            "frozen_commit": COMMIT,
        })
        return 2
    if CAMPAIGN.exists():
        emit({"phase": "refused", "blocker": "campaign directory already exists", "path": str(CAMPAIGN)})
        return 2
    emit({"phase": "calibration_start", "commit": COMMIT})
    calibration = driver.calibrate()
    decision = driver.admission_decision(calibration)
    emit({
        "phase": "calibration",
        "admitted": decision["admitted"],
        "blocker": decision["blocker"],
        "residuals": calibration.get("measured_residuals"),
        "measurement_sha256": calibration.get("measurement_sha256"),
        "executed": calibration.get("calibration_executed"),
    })
    if not decision["admitted"]:
        return 2
    dense_rows = []
    for mu in driver.MASS_OVER_RM:
        row = dense_reference(mu)
        dense_rows.append(row)
        emit({"phase": "dense_pair", **row})
    temp_calibration = Path("/tmp/nsc-rn-neutral-family-calibration.json")
    encoded = (json.dumps(to_jsonable(calibration), indent=2, sort_keys=True) + "\n").encode()
    temp_calibration.write_bytes(encoded)
    try:
        prepared = driver.prepare(
            CAMPAIGN, producer_commit=COMMIT, calibration_record=calibration,
            inputs=(str(temp_calibration),),
        )
    except driver.CampaignBlocker as blocked:
        emit({"phase": "refused", "blocker": blocked.payload.get("blocker"), "mode": "prepare"})
        return 2
    exclusive_bytes(CAMPAIGN / "calibration.json", encoded)
    script_bytes = Path(__file__).resolve().read_bytes()
    exclusive_bytes(CAMPAIGN / "replay-executor.py", script_bytes)
    durable_calibration = CAMPAIGN / "calibration.json"
    cases = [(mu, eps) for mu in driver.MASS_OVER_RM for eps in driver.ENERGY_FRACTIONS]
    for mu, eps in cases:
        driver.prepare(
            CAMPAIGN / case_name(mu, eps), producer_commit=COMMIT,
            calibration_record=calibration, inputs=(str(durable_calibration),),
        )
    emit({"phase": "prepared", "cases": [case_name(mu, eps) for mu, eps in cases]})
    ctx = get_context("spawn")
    manager = ctx.Manager()
    queue = manager.Queue()
    released = manager.Event()
    decision_box = manager.dict()
    stop_event = manager.Event()
    cpu_used = ctx.Value("d", time.process_time() - cpu_start)
    import threading
    watcher = threading.Thread(target=monitor, args=(stop_event, deadline), daemon=True)
    watcher.start()
    specs = []
    for mu, eps in cases:
        specs.append({
            "name": case_name(mu, eps),
            "mu": mu,
            "eps": eps,
            "campaign": str(CAMPAIGN),
            "commit": COMMIT,
            "calibration_path": str(durable_calibration),
            "queue": queue,
            "released": released,
            "decision": decision_box,
            "stop_event": stop_event,
            "cpu_used": cpu_used,
            "deadline": deadline,
        })
    pool = ctx.Pool(processes=6)
    async_results = [pool.apply_async(run_case, (spec,)) for spec in specs]
    messages = []
    try:
        while len(messages) < 6:
            remaining = max(1.0, deadline - time.time())
            messages.append(queue.get(timeout=remaining))
            emit({"phase": "forecast_message", "message": messages[-1]})
    except Exception as error:
        stop_event.set()
        released.set()
        emit({"phase": "forecast_incomplete", "blocker": str(error), "received": len(messages)})
        messages.extend(
            {"case": spec["name"], "ok": False, "mu": spec["mu"], "eps": spec["eps"],
             "blocker": "forecast message missing"}
            for spec in specs if spec["name"] not in {item.get("case") for item in messages}
        )
    selected, caps, table, remaining_cpu = choose_target(
        messages, float(cpu_used.value), max(0.0, deadline - time.time() - 30.0),
    )
    for name, choice in caps.items():
        decision_box[name] = choice
    emit({
        "phase": "admission",
        "target_over_rm": selected,
        "remaining_cpu_seconds": remaining_cpu,
        "table": table,
        "caps": {key: to_jsonable(value) for key, value in caps.items()},
    })
    released.set()
    results = []
    for item in async_results:
        try:
            results.append(item.get(timeout=max(30.0, deadline - time.time() + 120.0)))
        except Exception as error:
            results.append({"ok": False, "blocker": str(error), "stopped": "worker_exception"})
    pool.close()
    pool.join()
    stop_event.set()
    try:
        driver.write_dependency_manifest(
            CAMPAIGN, producer_commit=COMMIT,
            inputs=(str(durable_calibration), str(CAMPAIGN / "replay-executor.py")),
        )
    except driver.CampaignBlocker as blocked:
        emit({"phase": "dependency_blocker", "blocker": blocked.payload.get("blocker")})
    hashes = file_hashes(CAMPAIGN)
    summary = {
        "schema": driver.SCHEMA,
        "mode": "neutral-family-v1",
        "producing_commit": COMMIT,
        "git_head_at_execution": head,
        "campaign": str(CAMPAIGN),
        "calibration_admitted": True,
        "calibration_measurement_sha256": calibration.get("measurement_sha256"),
        "calibration_residuals": calibration.get("measured_residuals"),
        "calibration_source_binding": calibration.get("source_binding"),
        "dense_reference": dense_rows,
        "recorded_speedups": RECORDED_SPEEDUPS,
        "recorded_speedups_are_prior_measurements": True,
        "timing_admission": {"target_over_rm": selected, "table": table, "remaining_cpu_seconds": remaining_cpu},
        "cases": results,
        "prepare_bytes": prepared.get("bytes_written"),
        "cpu_seconds_aggregate": float(cpu_used.value),
        "wall_seconds": time.time() - wall_start,
        "memory_budget_bytes": MEMORY_BUDGET,
        "peak_rss_tree_bytes": None,
        "concurrency": {"processes": 6, "threads_per_process": 1, "pool": "spawn"},
        "output_sha256": hashes,
        "producer_hashes": driver.producers(),
        "fixed_control": driver.FIXED_CONTROL,
        "second_einstein_solution": False,
        "stations_requested": list(COORDINATE_STATIONS),
        "station_40": False,
        "event_horizon_claimed": False,
        "renewal_claimed": False,
        "baryogenesis_claimed": False,
        "cosmological_constant_claimed": False,
        "scope": (
            "Finite family comparison of transmission, trapping-horizon and field "
            "feedback, late remaining probability, and numerical margins. "
            "The frozen-metric arm is a control."
        ),
    }
    try:
        summary["peak_rss_tree_bytes"] = rss_tree_bytes(os.getpid())
    except Exception:
        summary["peak_rss_tree_bytes"] = None
    try:
        written = driver.write_assessment(
            CAMPAIGN / "summary.json", to_jsonable(summary), producer_commit=COMMIT,
        )
        summary_hash = written["sha256"]
        summary_bytes = written["bytes_written"]
    except driver.CampaignBlocker as blocked:
        emit({"phase": "summary_blocker", "blocker": blocked.payload.get("blocker")})
        return 3
    emit({
        "phase": "done",
        "summary_sha256": summary_hash,
        "summary_bytes": summary_bytes,
        "cpu_seconds_aggregate": summary["cpu_seconds_aggregate"],
        "wall_seconds": summary["wall_seconds"],
        "target_over_rm": selected,
        "cases": [
            {
                "case": item.get("case"),
                "ok": item.get("ok"),
                "stopped": item.get("stopped"),
                "blocker": item.get("blocker"),
                "final_t_over_rm": item.get("final_t_over_rm"),
            }
            for item in results
        ],
    })
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except driver.CampaignBlocker as blocked:
        emit({"phase": "refused", "blocker": blocked.payload.get("blocker")})
        sys.exit(2)
