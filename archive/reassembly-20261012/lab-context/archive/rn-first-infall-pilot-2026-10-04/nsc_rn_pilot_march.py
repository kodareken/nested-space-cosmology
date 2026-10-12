#!/usr/bin/env python3
"""Scratch matched-pair march. Not a sealed record and not a frozen-commit run."""
import hashlib
import json
import os
import time
from pathlib import Path

for name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[name] = "1"

import numpy as np

import derive_nsc_rn_parent as driver
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_source as source

ROOT = Path("/Users/admin/Documents/BlackHoles-Infinity")
SCRATCH = Path("/tmp/nsc-rn-pilot-rank9")
SCRATCH.mkdir(parents=True, exist_ok=True)
CPU_CAP = 600.0
WALL_CAP = 1000.0
OWNERS = (
    "lab/src/recursive_horizons/nsc_rn_pg.py",
    "lab/src/recursive_horizons/nsc_rn_source.py",
    "lab/src/recursive_horizons/nsc_rn_observables.py",
    "lab/src/recursive_horizons/nsc_rn_reference.py",
)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(name, payload):
    path = SCRATCH / name
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("WROTE", path, flush=True)


def moments_of(phi, grid):
    return source.radial_moments(phi[:, :, 0], grid["radius"], grid["weights"])


def partition(phi, grid, center=12.0):
    density = np.sum(np.abs(phi[:, :, 0]) ** 2, axis=0)
    weights = grid["weights"]
    total = float(np.sum(weights * density))
    radius = grid["radius"]
    def share(mask):
        return float(np.sum(weights[mask] * density[mask]) / total)
    return {
        "norm": total,
        "prob_abs_r_minus_12_le_2": share(np.abs(radius - center) <= 2.0),
        "prob_r_gt_20": share(radius > 20.0),
        "prob_r_lt_8": share(radius < 8.0),
    }


def accounting(coupled, fixed, grid, metric, rates, packet_admission):
    observed = driver.station_observation(coupled, grid, rates, packet=None)
    coupled_moments = moments_of(coupled["phi"], grid)
    fixed_moments = moments_of(fixed["phi"], grid)
    closure = (
        coupled["probability_stocks"]["remaining"]
        + coupled["probability_stocks"]["excision_outflow"]
        + coupled["probability_stocks"]["outer_outflow"]
    )
    fixed_closure = (
        fixed["probability_stocks"]["remaining"]
        + fixed["probability_stocks"]["excision_outflow"]
        + fixed["probability_stocks"]["outer_outflow"]
    )
    return {
        "t": coupled["clocks"]["t"],
        "t_over_rm": coupled["clocks"]["t"] / grid["r_m"],
        "fixed_t": fixed["clocks"]["t"],
        "coupled_clocks": coupled["clocks"],
        "fixed_clocks": fixed["clocks"],
        "coupled_probability": coupled["probability_stocks"],
        "fixed_probability": fixed["probability_stocks"],
        "coupled_mass_stocks": coupled["stocks"],
        "fixed_mass_stocks": fixed["stocks"],
        "coupled_inner_mass": coupled["inner_mass"],
        "fixed_inner_mass": fixed["inner_mass"],
        "coupled_flux_sum": closure,
        "fixed_flux_sum": fixed_closure,
        "coupled_moments": coupled_moments,
        "fixed_moments": fixed_moments,
        "coupled_partition": partition(coupled["phi"], grid),
        "fixed_partition": partition(fixed["phi"], grid),
        "lapse_difference_max": float(np.max(np.abs(rates["lapse"] - metric["lapse"]))),
        "shift_difference_max": float(np.max(np.abs(rates["shift"] - metric["shift"]))),
        "lapse_rate_max": float(np.max(np.abs(rates["lapse_rate"]))),
        "shift_rate_max": float(np.max(np.abs(rates["shift_rate"]))),
        "trapping_horizon": observed["trapping_horizon"],
        "charged_mass_inner": observed["charged_mass_inner"],
        "charged_mass_outer": observed["charged_mass_outer"],
        "curvature_status": observed["curvature_status"],
        "R4": None,
        "Ricci2": None,
        "required_time_jet_primitive": observed["required_time_jet_primitive"],
        "analytic_curvature_injected": False,
        "initial_packet": packet_admission,
    }


def save_checkpoint(tag, coupled, fixed, metric, grid, report):
    arrays = driver.state_arrays(coupled, fixed, metric)
    arrays["radius"] = np.asarray(grid["radius"])
    path = SCRATCH / f"{tag}.npz"
    np.savez_compressed(path, **arrays)
    sidecar = {
        "tag": tag,
        "payload_sha256": sha256(path),
        "payload_bytes": path.stat().st_size,
        "accounting": report,
    }
    dump(f"{tag}.json", sidecar)
    return str(path), sidecar["payload_sha256"]


def main():
    before = {name: sha256(ROOT / name) for name in OWNERS}
    dump("producer-hashes-before.json", before)
    geometry = driver.case_geometry(1.04)
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
        r_out=geometry["route"]["outer"],
    )
    started = time.perf_counter()
    packet = source.prepare_radial_packet(
        mass_over_rm=1.04, energy_target=1.0e-3, rank=driver.PILOT_RANK,
        model=grid["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"],
    )
    prepare_seconds = time.perf_counter() - started
    admission = driver.admit_prepared_packet(packet)
    admission["width_is_fit"] = False
    admission["requested_width"] = float(packet["radial_moments"]["requested_width"])
    admission["measured_width_rms"] = float(packet["radial_moments"]["width_rms"])
    admission["prepare_seconds"] = prepare_seconds
    admission["grid_points"] = int(grid["points"])
    admission["initial_partition"] = partition(packet["phi"], grid)
    dump("initial-packet.json", admission)
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
        raise SystemExit("prepared phi diverged before the march")
    if not np.array_equal(metric["lapse"], initial_rates["lapse"]):
        raise SystemExit("frozen lapse is not the sourced initial lapse")
    step = driver.matched_step(initial_rates["cfl_dt"], geometry["r_m"])
    if step["dt"] > 0.001 * geometry["r_m"] + 1e-15 or step["dt"] > step["cfl_dt"]:
        raise SystemExit("step exceeds the cap or the full CFL")
    dump("initial-step.json", step)
    save_checkpoint(
        "t0", coupled, fixed, metric, grid,
        accounting(coupled, fixed, grid, metric, initial_rates, admission),
    )
    wall0 = time.perf_counter()
    cpu = 0.0
    targets = (5.0, 10.0)
    history = []
    reached = []
    stop_reason = "not_started"
    for target in targets:
        if cpu >= CPU_CAP or (time.perf_counter() - wall0) >= WALL_CAP:
            stop_reason = "cap_before_target"
            break
        while coupled["clocks"]["t"] < target * geometry["r_m"] - 1e-12:
            if cpu >= CPU_CAP or (time.perf_counter() - wall0) >= WALL_CAP:
                stop_reason = "cap_during_target"
                break
            remaining = CPU_CAP - cpu
            chunk_time = min(target * geometry["r_m"], coupled["clocks"]["t"] + 0.25)
            marched = driver.march_pair(
                coupled, fixed, grid, metric, dt=step["dt"],
                station_time=chunk_time,
                cpu_seconds=remaining * driver.ADMISSION_FACTOR,
            )
            coupled = marched["coupled"]
            fixed = marched["fixed"]
            cpu += float(marched["cpu_seconds"])
            rates = pg.stage_rates(coupled, grid)
            report = accounting(coupled, fixed, grid, metric, rates, admission)
            report["chunk_stopped"] = marched["stopped"]
            report["chunk_steps"] = marched["steps"]
            report["cpu_seconds"] = cpu
            report["wall_seconds"] = time.perf_counter() - wall0
            tag = f"t{coupled['clocks']['t']:.4f}"
            path, digest = save_checkpoint(tag, coupled, fixed, metric, grid, report)
            history.append({
                "tag": tag, "path": path, "sha256": digest,
                "t_over_rm": report["t_over_rm"], "cpu_seconds": cpu,
                "stopped": marched["stopped"],
            })
            print(
                f"t/rm={report['t_over_rm']:.4f} cpu={cpu:.1f}s "
                f"stop={marched['stopped']} curv={report['curvature_status']}",
                flush=True,
            )
            if marched["stopped"] == "cpu_budget":
                stop_reason = "cpu_budget"
                break
        else:
            reached.append(target)
            stop_reason = "station"
            continue
        break
    after = {name: sha256(ROOT / name) for name in OWNERS}
    final = {
        "reached_stations": reached,
        "stop_reason": stop_reason,
        "cpu_seconds": cpu,
        "wall_seconds": time.perf_counter() - wall0,
        "prepare_seconds": prepare_seconds,
        "dt": step["dt"],
        "step": step,
        "points": int(grid["points"]),
        "rank": driver.PILOT_RANK,
        "hashes_before": before,
        "hashes_after": after,
        "hashes_unchanged": before == after,
        "history": history,
        "fixed_control": driver.FIXED_CONTROL,
        "sealed_lab_results_written": False,
    }
    dump("summary.json", final)
    print("DONE", json.dumps({
        "reached": reached, "stop": stop_reason, "cpu": cpu,
        "wall": final["wall_seconds"], "unchanged": before == after,
    }), flush=True)


if __name__ == "__main__":
    main()
