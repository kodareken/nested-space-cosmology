#!/usr/bin/env python3
"""One bounded successor of the sealed same-action conformal T=.05 episode.

Bitwise handoff of each v1 case; same original T=0 observer/source, active
L=Q,beta=0, unchanged Galerkin RK4. The physical question is resolved flow
across the natural width-2 regional boundaries, with maintained localization
and normal pressure/lapse accounting. Internal current maxima alone do not
resolve that question. No forced event/share/reversal or new force is used.
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
import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_frames as frames
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_spherical_coupling import CauchyState, PositiveChartExit

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-spherical-conformal-episode-v2.json"
NPZ = LAB / "results/development/nsc-spherical-conformal-episode-v2.npz"
DRIVER = Path(__file__).resolve()
START = .05
TARGET = .3
FRAME = .005
CPU_POOLS = {256: 600., 512: 1800.}
SAFETY = 1.25
RK4_SAFETY = 1.4
GRAM_LIMIT = 1e-8
ONE_PERCENT = .01
STATE_NAMES = episode.STATE_NAMES
GEOMETRY_NAMES = STATE_NAMES[:6]


def state_hash(state):
    digest = hashlib.sha256()
    for name in STATE_NAMES:
        values = np.ascontiguousarray(getattr(state, name))
        digest.update(name.encode())
        digest.update(str(values.dtype).encode())
        digest.update(str(values.shape).encode())
        digest.update(values.tobytes())
    return digest.hexdigest()


def scalar_integral(rows, name):
    times = np.array([row["time"] for row in rows])
    values = np.array([row[name] for row in rows])
    trap = float(np.sum(.5 * np.diff(times) * (values[:-1] + values[1:])))
    intervals = np.diff(times)
    uniform = len(intervals) > 0 and np.max(abs(intervals - intervals[0])) < 1e-12
    if uniform:
        indicators = episode.integrate(values, float(intervals[0]))
        indicators["trapezoid"] = trap
    else:
        indicators = {"trapezoid": trap, "simpson": None, "indicator": None}
    indicators["uniform_steps"] = bool(uniform)
    indicators["certified"] = False
    return indicators


def sample_summary(rows, handoff, dt_cap, stopped):
    first, last = rows[0], rows[-1]
    integrals = {name: scalar_integral(rows, name) for name in
                 ("forcing_norm", "quadrature_forcing_norm", "projection_forcing_norm", "fieldwork_power", "clock_rate")}
    elapsed_times = np.diff([row["time"] for row in rows])
    forcing = np.array([row["forcing_norm"] for row in rows])
    force_integrals = np.concatenate(([0.], np.cumsum(.5 * elapsed_times * (forcing[:-1] + forcing[1:]))))
    occupation0 = np.asarray(handoff["initial"]["occupation"])
    frozen_clock_rate = handoff["initial"]["clock_rate"]
    clock = np.array([row["proper_clock"] for row in rows])
    frozen_clock = np.array([row["time"] * frozen_clock_rate for row in rows])
    matched_clock = min(clock[-1], frozen_clock[-1])
    occ = np.array([row["occupation"] for row in rows])
    frozen = np.array([row["frozen_occupation"] for row in rows])
    matches = np.array([np.interp(matched_clock, clock, occ[:, index]) for index in range(2)])
    frozen_matches = np.array([np.interp(matched_clock, frozen_clock, frozen[:, index]) for index in range(2)])
    field_delta = last["field_energy"] - first["field_energy"]
    grav_delta = last["gravity_energy"] - first["gravity_energy"]
    end_budget = first["constraint_norm"] + integrals["forcing_norm"]["trapezoid"]
    return {"time": last["time"], "dt_cap": dt_cap, "stop_reason": stopped,
            "target_reached": stopped is None and abs(last["time"] - TARGET) < 1e-12,
            "steps": len(rows) - 1, "dt_min": float(np.min(elapsed_times)) if elapsed_times.size else None,
            "dt_max": float(np.max(elapsed_times)) if elapsed_times.size else None,
            "initial": first, "final": last, "integrals": integrals,
            "field_energy_change": field_delta, "gravity_energy_change": grav_delta,
            "energy_balance_change": field_delta + grav_delta,
            "fieldwork_difference": field_delta - integrals["fieldwork_power"]["trapezoid"],
            "constraint_end_sampled_budget": end_budget, "constraint_end_margin": end_budget - last["constraint_norm"],
            "max_sampled_budget_violation": float(max(row["constraint_norm"] - first["constraint_norm"] - integral for row, integral in zip(rows, force_integrals))),
            "joined_T0_constraint_sampled_budget": handoff["initial"]["constraint_norm"] + handoff["integrals"]["forcing_norm"]["trapezoid"] + integrals["forcing_norm"]["trapezoid"],
            "positive_chart": all(row["r_min"] > 0 and row["Q_min"] > 0 for row in rows),
            "r_min_window": min(row["r_min"] for row in rows), "Q_min_window": min(row["Q_min"] for row in rows),
            "gram_gap_max": max(row["gram_gap"] for row in rows),
            "energy_drift_max": max(abs(row["energy"] - first["energy"]) for row in rows),
            "occupation_change_continuation": occ[-1] - occ[0],
            "occupation_change_from_original_T0": occ[-1] - occupation0,
            "coupled_minus_frozen_coordinate_endpoint": occ[-1] - frozen[-1],
            "proper_clock_comparison": {"time": float(matched_clock), "coupled": matches, "frozen": frozen_matches,
                                        "difference": matches - frozen_matches, "sampling_interpolation_certified": False},
            "continuum_constraint_certified": False, "observable_error_certified": False, "renewal": False}


def boundary_assessment(physical_results):
    rows = []
    for index in range(4):
        endpoints = {name: result["final"]["windows"][index]["normal_boundary_flux"] for name, result in physical_results.items()}
        integrals = {name: result["windows"][index]["integrals"]["normal_boundary_flux"]["trapezoid"] for name, result in physical_results.items()}
        primary = "nf512_dt_0.0005"
        effect = abs(integrals[primary])
        time_gap = abs(integrals[primary] - integrals["nf512_dt_0.0010"])
        space_gap = abs(integrals[primary] - integrals["nf256_dt_0.0005"])
        numerical_indicator = sum(physical_results[primary]["windows"][index]["integrals"]["normal_boundary_flux"].get(key) or 0. for key in ("indicator",))
        resolved = effect > 0 and max(time_gap, space_gap, numerical_indicator) < ONE_PERCENT * effect
        rows.append({"interval": list(frames.BOUNDS[index:index + 2]), "endpoint_net_boundary_flux": endpoints,
                     "integrated_net_boundary_flux": integrals, "fine_effect": effect,
                     "fine_time_indicator": time_gap, "space_indicator": space_gap,
                     "frame_quadrature_indicator": numerical_indicator,
                     "resolved_at_one_percent_indicators": bool(resolved), "certified": False})
    return rows


def save(record, arrays):
    temporary = NPZ.with_suffix(".npz.tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    os.replace(temporary, NPZ)
    record["payload_sha256"] = episode.sha256(NPZ)
    predecessor_bytes = sum(path.stat().st_size for path in (episode.OUT, episode.NPZ, frames.OUT, frames.NPZ))
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    record["payload_bytes"] = NPZ.stat().st_size + OUT.stat().st_size
    record["combined_v1_v2_payload_bytes"] = predecessor_bytes + record["payload_bytes"]
    if record["combined_v1_v2_payload_bytes"] > episode.PAYLOAD_LIMIT:
        raise RuntimeError("combined v1/v2 saved payload exceeds 64MiB")
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")


def run():
    episode.refuse_existing_outputs(OUT, NPZ)
    core = json.loads(episode.OUT.read_text())
    if not core["full_target_reached"] or core["gauge"] != "conformal":
        raise ValueError("sealed conformal v1 episode required")
    paths = {"v1_json": episode.OUT, "v1_npz": episode.NPZ, "v1_frames_json": frames.OUT, "v1_frames_npz": frames.NPZ,
             "v5_json": episode.legacy.V5_JSON, "v5_npz": episode.legacy.V5_NPZ,
             "v1_driver": episode.DRIVER, "frames_driver": frames.DRIVER, "driver": DRIVER,
             "galerkin": galerkin._MODULE_PATH, "coupling": galerkin._COUPLING_PATH,
             "action": galerkin._ACTION_PATH, "source": galerkin._SOURCE_PATH}
    before = {name: episode.sha256(path) for name, path in paths.items()}
    record = {"schema": "NSC-SPHERICAL-CONFORMAL-EPISODE-v2", "gauge": "conformal", "start": START, "target": TARGET,
              "run_plan": episode.RUN_PLAN, "cpu_pool_budgets": CPU_POOLS, "forecast_safety": SAFETY,
              "effect_predeclared": "resolved integrated/spatial normal energy flux across natural width-2 boundaries plus maintained localization and pressure/lapse accounting",
              "internal_flux_max_alone_is_target": False, "forced_event_or_share_threshold": False,
              "adaptive_dt": "min(case dt cap, 1.4/current owned retained-band frequency, next .005 frame distance)",
              "absolute_rk4_cap": float(2 * np.sqrt(2)), "gram_limit": GRAM_LIMIT,
              "clock_protocol": core["clock_protocol"], "observer": core["observer"], "source_changed": False,
              "midtrajectory_reset": False, "radius_model_changed": False, "renewal": False,
              "predecessor_hashes": before, "results": {}, "series": {}, "physical_results": {}, "physical_series": {},
              "cpu_pools": {}, "handoffs": {}, "source_trajectory": "each exact matching v1 final Cauchy state, original saved T0 observer/columns for measurement and frozen control"}
    arrays = {}
    with np.load(episode.NPZ, allow_pickle=False) as stored:
        for nf in (256, 512):
            pool_start = time.process_time()
            grid = galerkin.build_grid(nf, quadrature=4 * nf, gauge="conformal")
            original = CauchyState(**{name: stored[f"nf{nf}_initial_{name}"].copy() for name in STATE_NAMES})
            observer = np.vstack((original.phi0, original.phi1))[:, :2].copy()
            for name in STATE_NAMES:
                arrays[f"nf{nf}_original_T0_{name}"] = getattr(original, name)
            specs = [spec for spec in episode.RUN_PLAN if spec["nf"] == nf]
            first_state = CauchyState(**{name: stored[specs[0]["name"] + "_final_" + name].copy() for name in STATE_NAMES})
            pilot_start = time.process_time()
            pilot = galerkin.rk4_step(grid, first_state, .0005)
            episode.observe(grid, pilot, original, observer, START + .0005)
            step_cost = time.process_time() - pilot_start
            frame_start = time.process_time()
            frames.frame_ledger(grid, first_state, START)
            frame_cost = time.process_time() - frame_start
            forecast = sum(math.ceil((TARGET - START) / spec["dt"]) * step_cost for spec in specs) + 2 * math.ceil((TARGET - START) / FRAME) * frame_cost
            admitted = SAFETY * forecast + time.process_time() - pool_start + 5 < CPU_POOLS[nf]
            record["cpu_pools"][f"nf{nf}"] = {"step_cost": step_cost, "frame_cost": frame_cost,
                                                    "target_forecast_cpu": forecast, "admitted_full_target": admitted,
                                                    "budget_seconds": CPU_POOLS[nf]}
            print("forecast", nf, record["cpu_pools"][f"nf{nf}"], flush=True)
            planned = TARGET if admitted else START + max(FRAME, math.floor((CPU_POOLS[nf] - 5) / (SAFETY * forecast) * (TARGET - START) / FRAME) * FRAME)
            planned = min(TARGET, planned)
            for spec in specs:
                current = CauchyState(**{name: stored[spec["name"] + "_final_" + name].copy() for name in STATE_NAMES})
                handoff_hash = state_hash(current)
                sample = episode.observe(grid, current, original, observer, START)
                handoff = core["results"][spec["name"]]
                for key in ("field_energy", "gravity_energy", "constraint_norm", "full_hamilton_max", "full_momentum_max"):
                    if abs(sample[key] - handoff["final"][key]) > 1e-12:
                        raise ValueError("handoff diagnostic is inconsistent")
                sample["proper_clock"] = handoff["final"]["proper_clock"]
                rows = [sample]
                physical_rows, profile_frames, geometry_rates, saved_frames = [], {}, {}, {name: [] for name in STATE_NAMES}
                frame_times = []
                time_value = START
                next_frame = START
                stopped = None
                run_start = time.process_time()
                cost = step_cost

                def keep_frame():
                    frame, profiles = frames.frame_ledger(grid, current, time_value)
                    frame["proper_clock"] = rows[-1]["proper_clock"]
                    physical_rows.append(frame)
                    frame_times.append(time_value)
                    rate, bundle = galerkin.compose_fine_hamiltonian(grid, current)
                    profile_frames.setdefault("L", []).append(bundle["fine_system"].length_density.copy())
                    profile_frames.setdefault("L_dot", []).append(bundle["lapse_dot"].copy())
                    for name in GEOMETRY_NAMES:
                        geometry_rates.setdefault(name, []).append(galerkin.prolong_geometry(grid, getattr(rate, name)))
                    for name in STATE_NAMES:
                        saved_frames[name].append(getattr(current, name).copy())
                    for name, value in profiles.items():
                        profile_frames.setdefault(name, []).append(value)

                keep_frame()
                next_frame += FRAME
                while time_value < planned - 1e-12:
                    if time.process_time() - pool_start + max(1., SAFETY * cost) + 5 >= CPU_POOLS[nf]:
                        stopped = "CPU_POOL_BUDGET"
                        break
                    frequency = galerkin.subspace_frequency(grid, current)[0]
                    dt = min(spec["dt"], RK4_SAFETY / max(frequency, 1.), next_frame - time_value, planned - time_value)
                    if dt <= 1e-12:
                        next_frame += FRAME
                        continue
                    step_start = time.process_time()
                    try:
                        candidate = galerkin.rk4_step(grid, current, dt)
                        following = episode.observe(grid, candidate, original, observer, time_value + dt)
                        if dt * following["omega"] > 2 * np.sqrt(2):
                            stopped = "ABSOLUTE_RK4_IMAGINARY_LIMIT"
                            break
                        if following["gram_gap"] > GRAM_LIMIT:
                            stopped = "GRAM_ADMISSIBILITY"
                            break
                        following["proper_clock"] = rows[-1]["proper_clock"] + .5 * dt * (rows[-1]["clock_rate"] + following["clock_rate"])
                        current = candidate
                        time_value += dt
                        rows.append(following)
                        if abs(time_value - next_frame) < 1e-10 or abs(time_value - planned) < 1e-10:
                            keep_frame()
                            next_frame += FRAME
                    except PositiveChartExit as failure:
                        stopped = failure.reason
                        break
                    cost = max(cost, time.process_time() - step_start)
                    if len(rows) % 100 == 0:
                        print("progress", spec["name"], "T", time_value, "poolCPU", time.process_time() - pool_start, flush=True)
                if abs(frame_times[-1] - time_value) > 1e-12:
                    keep_frame()
                summary = sample_summary(rows, handoff, spec["dt"], stopped)
                summary["cpu_seconds"] = time.process_time() - run_start
                summary["planned_common_pool_endpoint"] = planned
                summary["physical"] = frames.summarize(physical_rows) if len(physical_rows) > 1 else None
                record["results"][spec["name"]] = summary
                record["series"][spec["name"]] = rows
                record["physical_results"][spec["name"]] = summary["physical"]
                record["physical_series"][spec["name"]] = physical_rows
                record["handoffs"][spec["name"]] = {"v1_final_state_sha256": handoff_hash, "v2_initial_state_sha256": state_hash(CauchyState(**{name: saved_frames[name][0] for name in STATE_NAMES})), "bitwise_equal": True}
                arrays[spec["name"] + "_frame_times"] = np.array(frame_times)
                arrays[spec["name"] + "_proper_clock_frames"] = np.array([row["proper_clock"] for row in physical_rows])
                arrays[spec["name"] + "_times"] = np.array([row["time"] for row in rows])
                for name in STATE_NAMES:
                    arrays[spec["name"] + "_frames_" + name] = np.array(saved_frames[name])
                    arrays[spec["name"] + "_final_" + name] = getattr(current, name)
                for name, values in geometry_rates.items():
                    arrays[spec["name"] + "_frames_rate_" + name] = np.array(values)
                for name, values in profile_frames.items():
                    arrays[spec["name"] + "_frames_" + name] = np.array(values)
                print("complete", spec["name"], "T", time_value, "CPU", summary["cpu_seconds"], "fieldexchange", summary["field_energy_change"], "R", summary["final"]["constraint_norm"], flush=True)
            record["cpu_pools"][f"nf{nf}"]["actual_cpu_seconds"] = time.process_time() - pool_start
    all_complete = all(result["target_reached"] for result in record["results"].values())
    record["all_target_reached"] = all_complete
    if all_complete:
        record["boundary_assessment"] = boundary_assessment(record["physical_results"])
        resolved = any(row["resolved_at_one_percent_indicators"] for row in record["boundary_assessment"])
        record["resolved_width2_boundary_exchange"] = resolved
        record["verdict"] = "MEASURED_MAINTAINED_STRUCTURE_WITH_RESOLVED_REGIONAL_EXCHANGE" if resolved else "MEASURED_LOCALIZATION_REGIONAL_EXCHANGE_UNRESOLVED"
    else:
        record["boundary_assessment"] = []
        record["resolved_width2_boundary_exchange"] = False
        record["verdict"] = "PARTIAL_CONFORMAL_CONTINUATION_CHECKPOINT"
    record["source_hashes_after"] = {name: episode.sha256(path) for name, path in paths.items()}
    record["sealed_sources_unchanged_during_run"] = record["source_hashes_after"] == before
    record["continuum_constraint_certified"] = False
    record["observable_error_certified"] = False
    save(record, arrays)
    print("saved", OUT, "verdict", record["verdict"], "payload", record["combined_v1_v2_payload_bytes"], flush=True)
    return record


def verify_saved(*, source_ref=None):
    before = {path: episode.sha256(path) for path in (OUT, NPZ)}
    record = json.loads(OUT.read_text())
    if record.get("schema") != "NSC-SPHERICAL-CONFORMAL-EPISODE-v2" or record.get("payload_sha256") != episode.sha256(NPZ):
        raise ValueError("continuation schema or payload differs")
    paths = {"v1_json": episode.OUT, "v1_npz": episode.NPZ, "v1_frames_json": frames.OUT, "v1_frames_npz": frames.NPZ,
             "v5_json": episode.legacy.V5_JSON, "v5_npz": episode.legacy.V5_NPZ,
             "v1_driver": episode.DRIVER, "frames_driver": frames.DRIVER, "driver": DRIVER,
             "galerkin": galerkin._MODULE_PATH, "coupling": galerkin._COUPLING_PATH,
             "action": galerkin._ACTION_PATH, "source": galerkin._SOURCE_PATH}
    episode.check_recorded_sources(record.get("predecessor_hashes"), paths, source_ref=source_ref)
    episode.check_saved_values(record["predecessor_hashes"], record["source_hashes_after"], "source hashes")
    v1 = json.loads(episode.OUT.read_text())
    with np.load(NPZ, allow_pickle=False) as payload, np.load(episode.NPZ, allow_pickle=False) as predecessor:
        for spec in episode.RUN_PLAN:
            name = spec["name"]
            for field in STATE_NAMES:
                if not np.array_equal(payload[name + "_frames_" + field][0], predecessor[name + "_final_" + field]):
                    raise ValueError("Cauchy handoff differs: " + name + "." + field)
                if not np.array_equal(payload[name + "_frames_" + field][-1], payload[name + "_final_" + field]):
                    raise ValueError("final Cauchy frame differs: " + name)
                origin = f"nf{spec['nf']}_original_T0_" + field
                if not np.array_equal(payload[origin], predecessor[f"nf{spec['nf']}_initial_" + field]):
                    raise ValueError("original preparation differs: " + origin)
            if not np.array_equal(payload[name + "_frames_L_dot"], payload[name + "_frames_rate_Q"]):
                raise ValueError("dynamic lapse rate differs: " + name)
            fresh = sample_summary(record["series"][name], v1["results"][name], spec["dt"], record["results"][name]["stop_reason"])
            for key in fresh:
                episode.check_saved_values(fresh[key], record["results"][name][key], name + "." + key)
            episode.check_saved_values(frames.summarize(record["physical_series"][name]), record["physical_results"][name], name + ".physical")
    episode.check_saved_values(boundary_assessment(record["physical_results"]), record["boundary_assessment"], "boundary assessment")
    if before != {path: episode.sha256(path) for path in before}:
        raise RuntimeError("read-only check changed continuation bytes")
    return {"status": record["verdict"], "wrote": False, "source_ref": source_ref}


if __name__ == "__main__":
    if sys.argv[1:] == ["--check"]:
        print(json.dumps(verify_saved(source_ref=episode.SEALED_SOURCE_REF)))
    elif len(sys.argv) == 1:
        run()
    else:
        raise SystemExit("Use --check for sealed evidence, or no arguments to create new outputs")
