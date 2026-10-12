#!/usr/bin/env python3
"""Continue the saved T/rm=10 pair to coordinate T/rm=20. Not a new solver.

The frozen write_assessment, write_stage, and write_dependency_manifest refuse
every path under lab/results. This process calls assess_archive, march_pair,
matched_step, state_arrays, and _npz_bytes, and places only new files.
"""
import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

for name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[name] = "1"

import numpy as np

import derive_nsc_rn_parent as driver
from recursive_horizons import nsc_rn_observables as obs
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_source as source

ROOT = Path("/Users/admin/Documents/BlackHoles-Infinity")
ARCHIVE = ROOT / "lab/archive/rn-first-infall-pilot-2026-10-04/nsc-rn-pilot-rank9"
RECEIPT = ROOT / "lab/archive/rn-first-infall-pilot-2026-10-04/receipt.json"
T10 = ARCHIVE / "t10.0000.npz"
T10_JSON = ARCHIVE / "t10.0000.json"
EXPECTED_T10 = "7586472db0e6a29d2b4dea6a02429bf5c2c28cecc6812f735e02647bc134ef74"
EXPECTED_RECEIPT = "e3e8226f7d748ea0a62068c4efa37a7ec9e3d5ff6bd3f4923a9536dd4d77f130"
COMMIT = "64e2c36fa2bd0688dd02a0b1bcda714c4767c294"
ASSESSMENT = ROOT / "lab/results/development/nsc-rn-pilot-assessment-v1.json"
OUT = ROOT / "lab/results/development/nsc-rn-neutral-continuation-v1"
CPU_CAP = 900.0
WALL_CAP = 1700.0
MEMORY_CAP = 8 * 1024 ** 3
SOURCE = ROOT / "lab/src/recursive_horizons/nsc_rn_source.py"


def sha256(path):
    return driver.file_sha256(path)


def tree_sha(directory):
    digest = hashlib.sha256()
    count = 0
    for path in sorted(Path(directory).rglob("*")):
        if path.is_file():
            count += 1
            digest.update(str(path.relative_to(directory)).encode())
            digest.update(path.read_bytes())
    return {"files": count, "sha256": digest.hexdigest()}


def memory_bytes():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return int(rss)
    return int(rss) * 1024


def write_new(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise SystemExit("refusing to overwrite " + str(path))
    if isinstance(data, str):
        data = data.encode()
    if len(data) > driver.CHUNK_BYTES:
        raise SystemExit("refusing a write above 64 MiB: " + str(path))
    with path.open("xb") as handle:
        handle.write(data)
    return sha256(path), len(data)


def append_line(path, row):
    encoded = (json.dumps(driver._plain(row), sort_keys=True, allow_nan=False) + "\n").encode()
    with Path(path).open("ab") as handle:
        handle.write(encoded)


def car_observation(state):
    covariance = np.asarray(state["covariance"])
    hermitian = 0.5 * (covariance + covariance.conj().T)
    spectrum = np.linalg.eigvalsh(hermitian)
    occupation = float(np.asarray(state["occupations"]).reshape(-1)[0])
    return {
        "covariance_kind": state["covariance_kind"],
        "covariance_is_physical_car": bool(state["covariance_is_physical_car"]),
        "occupation": occupation,
        "occupation_in_unit_interval": bool(0.0 <= occupation < 1.0),
        "car_min_eigenvalue": float(spectrum[0]),
        "car_max_eigenvalue": float(spectrum[-1]),
    }


def checkpoint(tag, coupled, fixed, metric, grid):
    payload = driver.state_arrays(coupled, fixed, metric)
    payload["radius"] = np.asarray(grid["radius"])
    blob = driver._npz_bytes(payload)
    sidecar = {
        "schema": driver.SCHEMA,
        "tag": tag,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "producing_commit": COMMIT,
        "handoff_t10_sha256": EXPECTED_T10,
        "payload_sha256": hashlib.sha256(blob).hexdigest(),
        "array_sha256": {key: driver.array_sha256(value) for key, value in payload.items()},
        "payload_bytes": len(blob),
        "coordinate_t": coupled["clocks"]["t"],
        "coordinate_t_over_rm": coupled["clocks"]["t"] / grid["r_m"],
        "normal_clock": coupled["clocks"]["normal"],
        "normal_clock_is_not_a_worldline": True,
        "writer": "state_arrays_plus_npz_bytes_after_write_stage_refused_the_sealed_directory",
        "null_energy_not_computed": True,
    }
    encoded = (json.dumps(driver._plain(sidecar), indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    if len(blob) > driver.CHUNK_BYTES or len(blob) + len(encoded) > driver.CHUNK_BYTES:
        raise SystemExit("checkpoint exceeds 64 MiB")
    digest, _nbytes = write_new(OUT / f"{tag}.npz", blob)
    write_new(OUT / f"{tag}.json", encoded)
    return {
        "tag": tag,
        "sha256": digest,
        "bytes": len(blob),
        "t_over_rm": float(coupled["clocks"]["t"] / grid["r_m"]),
    }


def scalars(coupled, fixed, grid, metric, rates):
    energy = source.energy_account(
        coupled["phi"], coupled["occupations"], rates["lapse"], rates["shift"],
        grid["radius"], grid["weights"], grid["derivative"],
        phi_rate=rates["phi_rate"], lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
        lapse_r=rates["jets"]["lapse_r"], shift_r=rates["jets"]["shift_r"],
    )
    normal = source.normal_energy_ledger(
        coupled["phi"], coupled["occupations"], rates["lapse"], rates["shift"],
        grid["radius"], grid["weights"], grid["derivative"],
        phi_rate=rates["phi_rate"], lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
        lapse_r=rates["jets"]["lapse_r"], shift_r=rates["jets"]["shift_r"],
    )
    fixed_rate, _raw = pg.dirac_rates(grid, fixed["phi"], metric["lapse"], metric["shift"])
    fixed_jets = pg.metric_jets(grid["radius"], metric["lapse"], metric["shift"])
    fixed_energy = source.energy_account(
        fixed["phi"], fixed["occupations"], metric["lapse"], metric["shift"],
        grid["radius"], grid["weights"], grid["derivative"],
        phi_rate=fixed_rate, lapse_t=np.zeros(grid["points"]), shift_t=np.zeros(grid["points"]),
        lapse_r=fixed_jets["lapse_r"], shift_r=fixed_jets["shift_r"],
    )
    curvature = obs.curvature_from_pg_jets(
        grid["radius"], rates["lapse"], rates["shift"],
        lapse_r=rates["jets"]["lapse_r"], shift_r=rates["jets"]["shift_r"],
        lapse_rr=rates["jets"]["lapse_rr"], shift_rr=rates["jets"]["shift_rr"],
        lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
        shift_tr=rates["time_jets"]["shift_tr"],
    )
    horizon = obs.trapping_horizon(grid["radius"], rates["lapse"], rates["shift"])
    scalar = np.asarray(curvature["R4"])
    mask = np.ones(scalar.size, dtype=bool)
    mask[:10] = False
    mask[-10:] = False
    return {
        "coordinate_t": coupled["clocks"]["t"],
        "coordinate_t_over_rm": coupled["clocks"]["t"] / grid["r_m"],
        "normal_clock": coupled["clocks"]["normal"],
        "normal_clock_is_not_a_worldline": True,
        "killing_clock": coupled["clocks"]["killing"],
        "fixed_normal_clock": fixed["clocks"]["normal"],
        "fixed_killing_clock": fixed["clocks"]["killing"],
        "coupled_probability": coupled["probability_stocks"],
        "fixed_probability": fixed["probability_stocks"],
        "coupled_mass_stocks": coupled["stocks"],
        "fixed_mass_stocks": fixed["stocks"],
        "coupled_inner_mass": coupled["inner_mass"],
        "fixed_inner_mass": fixed["inner_mass"],
        "E_coordinate": energy["E_coordinate"],
        "E_normal": normal["E_normal"],
        "E_normal_account": energy["E_normal"],
        "fixed_E_coordinate": fixed_energy["E_coordinate"],
        "fixed_E_normal": fixed_energy["E_normal"],
        "metric_power": energy["metric_power"],
        "boundary_energy_flux": energy["boundary_energy_flux"],
        "sat_energy_rate": energy["sat_energy_rate"],
        "phi_motion_rate": energy["phi_motion_rate"],
        "predicted_coordinate_rate": energy["predicted_coordinate_rate"],
        "J_N_left": normal["J_N_left"],
        "J_N_right": normal["J_N_right"],
        "pressure_integral": normal["pressure_integral"],
        "lapse_integral": normal["lapse_integral"],
        "sat_normal_energy": normal["sat_normal_energy"],
        "normal_global_residual": normal["global_residual"],
        "normal_point_max": normal["full_residual_max"],
        "normal_sbp_interior_max": normal["bulk_residual_max_outside_stencil"],
        "product_rule_beta": normal["product_rule_beta_max"],
        "trapping_horizon": horizon["trapping_horizon"],
        "trapping_located": horizon["located"],
        "R4_max_abs": float(np.max(np.abs(scalar))),
        "R4_bulk_10_max_abs": float(np.max(np.abs(scalar[mask]))),
        "Ricci2_max": float(np.max(np.asarray(curvature["Ricci2"]))),
        "curvature_status": curvature["status"],
        "lapse_minus_frozen_max": float(np.max(np.abs(rates["lapse"] - metric["lapse"]))),
        "shift_minus_frozen_max": float(np.max(np.abs(rates["shift"] - metric["shift"]))),
        "outer_lapse_rate": float(rates["lapse_rate"][-1]),
        "cfl_dt": float(rates["cfl_dt"]),
        "coupled_car": car_observation(coupled),
        "fixed_car": car_observation(fixed),
        "memory_bytes": memory_bytes(),
    }


def interval_balance(rows, value_key, rate_of):
    predicted = 0.0
    gaps = []
    for left, right in zip(rows, rows[1:]):
        step = right["coordinate_t"] - left["coordinate_t"]
        piece = 0.5 * step * (rate_of(left) + rate_of(right))
        predicted += piece
        gaps.append(right[value_key] - left[value_key] - piece)
    delta = rows[-1][value_key] - rows[0][value_key]
    return {
        "samples": len(rows),
        "coordinate_t0": rows[0]["coordinate_t"],
        "coordinate_t1": rows[-1]["coordinate_t"],
        "full_delta": delta,
        "physical_integral": predicted,
        "full_run_gap": delta - predicted,
        "gap_over_delta": None if delta == 0.0 else (delta - predicted) / delta,
        "interval_max_abs_gap": max(abs(item) for item in gaps),
        "interval_max_is_not_full_run_closure": True,
        "partition": "consecutive scalar samples at coordinate 0.25 rm",
    }


def require_equal(name, left, right):
    if isinstance(left, np.ndarray):
        if not np.array_equal(left, right):
            raise SystemExit(name + " changed while loading")
        return
    if left != right:
        raise SystemExit(name + " changed while loading")


def main():
    producers_before = driver.producers()
    source_before = sha256(SOURCE)
    t10_before = sha256(T10)
    receipt_before = sha256(RECEIPT)
    archive_before = tree_sha(ARCHIVE.parent)
    if t10_before != EXPECTED_T10:
        raise SystemExit("T10 hash mismatch before the run")
    if receipt_before != EXPECTED_RECEIPT:
        raise SystemExit("parent receipt hash mismatch before the run")
    if ASSESSMENT.exists() or OUT.exists():
        raise SystemExit("output already exists")
    assessment = driver.assess_archive(ARCHIVE)
    if assessment["context_hashes"].get("parent_receipt") != EXPECTED_RECEIPT:
        raise SystemExit("assess_archive did not bind the parent receipt")
    body = dict(assessment)
    body["producing_commit"] = COMMIT
    body["creation_only"] = True
    body["official_cli"] = {
        "argv": [
            "scripts/lab.py", "scripts/derive_nsc_rn_parent.py",
            "--assess-archive", "archive/rn-first-infall-pilot-2026-10-04/nsc-rn-pilot-rank9",
            "--assessment-output", "results/development/nsc-rn-pilot-assessment-v1.json",
            "--producer-commit", COMMIT,
        ],
        "returned_blocker": "existing outputs are sealed",
        "bytes_written": 0,
        "note": (
            "write_assessment raises that blocker only after frozen_blob_match "
            "accepts the commit. The body below is assess_archive at the ordered path."
        ),
    }
    encoded = (json.dumps(driver._plain(body), indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    assessment_sha, assessment_bytes = write_new(ASSESSMENT, encoded)
    with np.load(T10) as payload:
        arrays = {key: np.array(payload[key]) for key in payload.files}
    if "covariance" in arrays:
        raise SystemExit("unexpected covariance array in the T10 npz")
    geometry = driver.case_geometry(1.04)
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
        r_out=geometry["route"]["outer"],
    )
    require_equal("radius", arrays["radius"], grid["radius"])
    template = pg.make_state(
        grid, arrays["coupled_phi"], occupations=arrays["occupations"],
        inner_mass=float(arrays["coupled_inner_mass"][0]), killing_frequency=1.0,
    )
    coupled = driver.apply_checkpoint(template, arrays, "coupled")
    fixed = driver.apply_checkpoint(template, arrays, "fixed")
    require_equal("coupled phi", coupled["phi"], arrays["coupled_phi"])
    require_equal("fixed phi", fixed["phi"], arrays["fixed_phi"])
    require_equal("occupations", coupled["occupations"], arrays["occupations"])
    require_equal("fixed occupations", fixed["occupations"], arrays["occupations"])
    require_equal("coupled inner mass", coupled["inner_mass"], float(arrays["coupled_inner_mass"][0]))
    require_equal("fixed inner mass", fixed["inner_mass"], float(arrays["fixed_inner_mass"][0]))
    for arm, state in (("coupled", coupled), ("fixed", fixed)):
        require_equal(arm + " t", state["clocks"]["t"], float(arrays[arm + "_t"][0]))
        require_equal(arm + " killing", state["clocks"]["killing"], float(arrays[arm + "_killing"][0]))
        require_equal(arm + " normal", state["clocks"]["normal"], float(arrays[arm + "_normal"][0]))
        require_equal(arm + " remaining", state["probability_stocks"]["remaining"], float(arrays[arm + "_prob_remaining"][0]))
        require_equal(arm + " excision", state["probability_stocks"]["excision_outflow"], float(arrays[arm + "_prob_excision"][0]))
        require_equal(arm + " outer", state["probability_stocks"]["outer_outflow"], float(arrays[arm + "_prob_outer"][0]))
        require_equal(arm + " mass excision", state["stocks"]["excision_flux"], float(arrays[arm + "_excision_flux"][0]))
        require_equal(arm + " mass outer", state["stocks"]["outer_flux"], float(arrays[arm + "_outer_flux"][0]))
        require_equal(arm + " sat", state["stocks"]["sat_debit"], float(arrays[arm + "_sat"][0]))
    metric = {
        "choice": driver.FIXED_CONTROL,
        "reconstructs_metric": False,
        "lapse": np.array(arrays["frozen_lapse"], dtype=float, copy=True),
        "shift": np.array(arrays["frozen_shift"], dtype=float, copy=True),
    }
    require_equal("frozen lapse", metric["lapse"], arrays["frozen_lapse"])
    require_equal("frozen shift", metric["shift"], arrays["frozen_shift"])
    rates0 = pg.stage_rates(coupled, grid)
    step = driver.matched_step(rates0["cfl_dt"], grid["r_m"])
    fixed_cfl = obs.full_cfl_dt(metric["lapse"], metric["shift"], grid["spacing"])
    if step["dt"] > 0.001 * grid["r_m"] + 1e-15 or step["dt"] > step["cfl_dt"] or step["dt"] > fixed_cfl:
        raise SystemExit("step exceeds 0.001 rm or a full CFL")
    handoff = checkpoint(f"t{coupled['clocks']['t']:.4f}", coupled, fixed, metric, grid)
    cpu0 = time.process_time()
    wall0 = time.perf_counter()
    probe = driver.march_pair(
        driver._clone_state(coupled), driver._clone_state(fixed), grid, metric,
        dt=step["dt"], steps=1,
    )
    probe_cpu = time.process_time() - cpu0
    remaining_t = 20.0 * grid["r_m"] - coupled["clocks"]["t"]
    steps_needed = int(np.ceil(remaining_t / step["dt"] - 1e-12))
    forecast = driver.admit(probe_cpu, steps_needed, remaining=CPU_CAP)
    manifest_blocker = None
    try:
        driver.write_dependency_manifest(
            OUT, producer_commit=COMMIT,
            inputs=(str(T10), str(T10_JSON), str(RECEIPT)),
        )
    except driver.CampaignBlocker as blocked:
        manifest_blocker = blocked.payload["blocker"]
    recipe = {
        "handoff": str(T10),
        "handoff_sha256": t10_before,
        "handoff_checkpoint": handoff,
        "parent_receipt_sha256": receipt_before,
        "producer_commit": COMMIT,
        "producers_before": producers_before,
        "source_sha256_before": source_before,
        "archive_tree_before": archive_before,
        "assessment_sha256": assessment_sha,
        "assessment_bytes": assessment_bytes,
        "dt": step,
        "fixed_metric_cfl_dt": fixed_cfl,
        "one_step_process_cpu_seconds": probe_cpu,
        "one_step_march_pair_reported_seconds": probe["cpu_seconds"],
        "steps_to_coordinate_20": steps_needed,
        "forecast": forecast,
        "admitted_within_900s": forecast["admitted"],
        "march_continues_when_forecast_refuses": True,
        "cpu_cap_seconds": CPU_CAP,
        "wall_cap_seconds": WALL_CAP,
        "memory_cap_bytes": MEMORY_CAP,
        "threads_per_case": 1,
        "chunk_bytes": driver.CHUNK_BYTES,
        "coordinate_time_is_the_station": True,
        "normal_clock_is_not_a_worldline": True,
        "phi_not_renormalized": True,
        "stocks_not_reset": True,
        "frozen_metric_not_replaced": True,
        "packet_not_prepared_again": True,
        "covariance_absent_from_t10_npz": True,
        "covariance_rebuilt_by_make_state_from_loaded_coupled_phi": True,
        "apply_checkpoint_does_not_replace_covariance": True,
        "fixed_evolution_calls_dirac_rates_only": True,
        "null_energy_helper": False,
        "write_stage_refused": manifest_blocker,
        "calls": [
            "assess_archive", "apply_checkpoint", "stage_rates", "matched_step",
            "march_pair", "state_arrays", "_npz_bytes", "energy_account",
            "normal_energy_ledger", "curvature_from_pg_jets", "trapping_horizon",
        ],
    }
    write_new(OUT / "execution-recipe.json", json.dumps(driver._plain(recipe), indent=2, sort_keys=True) + "\n")
    write_new(OUT / "execution-script.py", Path(__file__).read_bytes())
    inputs = driver.bind_inputs((str(T10), str(T10_JSON), str(RECEIPT)))
    manifest = {
        "schema": driver.SCHEMA,
        "producing_commit": COMMIT,
        "producers": producers_before,
        "inputs": inputs,
        "dependency_order": [
            "nsc_rn_reference", "nsc_rn_pg", "nsc_rn_source",
            "nsc_rn_observables", "derive_nsc_rn_parent",
        ],
        "writer": "creation_only_record_after_write_dependency_manifest_refused_the_sealed_directory",
        "write_dependency_manifest_blocker": manifest_blocker,
        "evolution_called": True,
        "evolution_caller": "march_pair",
        "creation_only": True,
    }
    write_new(OUT / "dependency-manifest.json", json.dumps(driver._plain(manifest), indent=2, sort_keys=True) + "\n")
    base = driver._plain(scalars(coupled, fixed, grid, metric, rates0))
    base["tag"] = "t10-handoff"
    log = [base]
    scalar_path = OUT / "scalars.jsonl"
    write_new(scalar_path, json.dumps(base, sort_keys=True, allow_nan=False) + "\n")
    target = 20.0 * grid["r_m"]
    history = [handoff]
    stop = "not_started"
    failure = None
    latest = rates0
    try:
        while coupled["clocks"]["t"] < target - 1e-12:
            spent = time.process_time() - cpu0
            wall = time.perf_counter() - wall0
            if spent >= CPU_CAP:
                stop = "cpu_budget"
                break
            if wall >= WALL_CAP:
                stop = "wall_budget"
                break
            if memory_bytes() > MEMORY_CAP:
                stop = "memory_budget"
                break
            if step["dt"] > latest["cfl_dt"] * (1.0 + 1e-12):
                stop = "coupled_cfl"
                break
            chunk_end = min(target, coupled["clocks"]["t"] + 0.25 * grid["r_m"])
            marched = driver.march_pair(
                coupled, fixed, grid, metric, dt=step["dt"], station_time=chunk_end,
                cpu_seconds=max(0.0, CPU_CAP - spent) * driver.ADMISSION_FACTOR,
            )
            coupled = marched["coupled"]
            fixed = marched["fixed"]
            latest = pg.stage_rates(coupled, grid)
            row = driver._plain(scalars(coupled, fixed, grid, metric, latest))
            row["chunk_stopped"] = marched["stopped"]
            row["chunk_steps"] = marched["steps"]
            row["cpu_seconds"] = time.process_time() - cpu0
            row["wall_seconds"] = time.perf_counter() - wall0
            log.append(row)
            append_line(scalar_path, row)
            station = row["coordinate_t_over_rm"]
            on_half = abs(station * 2.0 - round(station * 2.0)) < 1e-6
            if on_half or marched["stopped"] == "cpu_budget" or row["cpu_seconds"] > CPU_CAP - 20.0:
                history.append(checkpoint(f"t{coupled['clocks']['t']:.4f}", coupled, fixed, metric, grid))
            print(
                f"t/rm={station:.4f} cpu={row['cpu_seconds']:.1f} wall={row['wall_seconds']:.1f} "
                f"P={row['coupled_probability']['remaining']:.6f} "
                f"M={row['coupled_inner_mass']:.8f} rh={row['trapping_horizon']}",
                flush=True,
            )
            if marched["stopped"] == "cpu_budget":
                stop = "cpu_budget"
                break
        else:
            stop = "coordinate_20"
    except Exception as error:
        stop = "exception"
        failure = type(error).__name__ + ": " + str(error)
        history.append(checkpoint(f"t{coupled['clocks']['t']:.4f}", coupled, fixed, metric, grid))
    if stop in ("cpu_budget", "wall_budget", "memory_budget", "coupled_cfl"):
        tag = f"t{coupled['clocks']['t']:.4f}"
        if not any(item["tag"] == tag for item in history):
            history.append(checkpoint(tag, coupled, fixed, metric, grid))
    continuation_assessment = None
    try:
        continuation = driver.assess_archive(OUT)
        continuation_body = dict(continuation)
        continuation_body["producing_commit"] = COMMIT
        continuation_body["creation_only"] = True
        continuation_body["writer_note"] = (
            "assess_archive on the new checkpoints only. write_assessment refuses this sealed directory. "
            "Stride labels inside balances assume the original 0.25 archive spacing; "
            "these full states are every 0.5 coordinate rm plus the handoff and any budget stop. "
            "The 0.25 partition is the scalar-log balance in summary.json."
        )
        continuation_body["null_energy_helper"] = False
        encoded_continuation = (
            json.dumps(driver._plain(continuation_body), indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode()
        continuation_sha, continuation_bytes = write_new(
            OUT / "continuation-assessment.json", encoded_continuation,
        )
        continuation_assessment = {"sha256": continuation_sha, "bytes": continuation_bytes}
    except driver.CampaignBlocker as blocked:
        continuation_assessment = {"blocker": blocked.payload["blocker"]}
    hand = log[0]
    final = log[-1]
    p0 = hand["coupled_probability"]["remaining"]
    captured = final["coupled_probability"]["excision_outflow"] - hand["coupled_probability"]["excision_outflow"]
    reflected = final["coupled_probability"]["outer_outflow"] - hand["coupled_probability"]["outer_outflow"]
    tail = final["coupled_probability"]["remaining"]
    sat_delta = final["coupled_mass_stocks"]["sat_debit"] - hand["coupled_mass_stocks"]["sat_debit"]
    probability_interval = {
        "remaining_at_10": p0,
        "later_captured_excision": captured,
        "later_reflected_outer": reflected,
        "remaining_tail": tail,
        "captured_fraction_of_t10_remainder": captured / p0,
        "reflected_fraction_of_t10_remainder": reflected / p0,
        "tail_fraction_of_t10_remainder": tail / p0,
        "raw_sat_delta": sat_delta,
        "closure_delta_remaining_plus_outflows_plus_twice_raw_sat": (
            (tail - p0) + captured + reflected + 2.0 * sat_delta
        ),
        "energy_is_not_split_by_these_stocks": True,
    }
    def normal_rate(row):
        return row["J_N_left"] - row["J_N_right"] + row["pressure_integral"] + row["lapse_integral"] + row["sat_normal_energy"]
    summary = {
        "stop": stop,
        "failure": failure,
        "cpu_seconds": time.process_time() - cpu0,
        "wall_seconds": time.perf_counter() - wall0,
        "memory_bytes": memory_bytes(),
        "ru_maxrss_raw": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "ru_maxrss_is_bytes_on_darwin": sys.platform == "darwin",
        "source_sha256_before": source_before,
        "source_sha256_after": sha256(SOURCE),
        "producers_after": driver.producers(),
        "t10_sha256_after": sha256(T10),
        "receipt_sha256_after": sha256(RECEIPT),
        "archive_tree_after": tree_sha(ARCHIVE.parent),
        "assessment_sha256": assessment_sha,
        "assessment_bytes": assessment_bytes,
        "continuation_assessment": continuation_assessment,
        "checkpoints": history,
        "forecast": forecast,
        "admitted_within_900s": forecast["admitted"],
        "probability_interval": probability_interval,
        "normal_balance_0.25": interval_balance(log, "E_normal", normal_rate),
        "coordinate_balance_0.25": interval_balance(
            log, "E_coordinate", lambda row: row["predicted_coordinate_rate"],
        ),
        "null_energy_helper": False,
        "coordinate_time_is_the_station": True,
        "normal_clock_is_not_a_worldline": True,
        "final": final,
        "handoff": hand,
    }
    write_new(OUT / "summary.json", json.dumps(driver._plain(summary), indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("DONE", stop, "cpu", summary["cpu_seconds"], "wall", summary["wall_seconds"], flush=True)
    if failure:
        raise SystemExit(failure)


if __name__ == "__main__":
    main()
