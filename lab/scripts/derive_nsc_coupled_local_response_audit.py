#!/usr/bin/env python3
"""Independent column-control audit on the sealed local-response v2 window.

The geometry is the same stored piecewise-linear Q, and the observer is the
original T=0 pair. This driver does not call the consumer's temporal solver,
cross split, projection, or omission helpers. It reuses the checked Galerkin
matrix representation, transports columns, and never forms a propagator W.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(variable, "1")

import numpy as np
from scipy.sparse.linalg import expm_multiply

from recursive_horizons.nsc_coupled_local_response import build_case_grid, galerkin_matrix

LAB = Path(__file__).resolve().parents[1]
DEV = LAB / "results" / "development"
ORIGIN = DEV / "nsc-spherical-feedback-episode-v1.npz"
EPISODE = DEV / "nsc-regeneration-episode-v1.npz"
RESPONSE = DEV / "nsc-coupled-local-response-v2.npz"
RESPONSE_JSON = RESPONSE.with_suffix(".json")
OUTPUT = DEV / "nsc-coupled-local-response-audit-v1.json"
OUTPUT_NPZ = OUTPUT.with_suffix(".npz")
CASE = "nf512_dtmax_0_00025"
HANDOFF = "nf512_dt_0_0005"
WEIGHTS = np.array([0.75, 0.75, 0.5, 0.5, 0.25, 0.25])
SEALED = {
    ORIGIN: "2864d3a8d3ba413e840c395961f14a3d71eabee702c27fd9f47c7d333eb24b77",
    EPISODE: "17364df24efead5814c3ab2a2ab5f82d9e09bf412c8c82e68d57a3b38498d7a1",
    RESPONSE: "84c09302413d48c2f5dfaba32090b8169d0a4c1dbbe9dfc7a971debaaa525d14",
    RESPONSE_JSON: "5330ccef1ba2fc4d80e0e0bdc0a98b1117e851f3873a46a0e057d63aae7457bb",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def covariance(amplitudes, weights=WEIGHTS):
    return np.einsum("tak,k,tbk->tab", amplitudes, weights, amplitudes.conj())


def occupation(matrix):
    return np.diagonal(matrix, axis1=1, axis2=2).real


def max_abs(values):
    return float(np.max(np.abs(values)))


def array_id(values):
    contiguous = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(contiguous.dtype).encode("ascii"))
    digest.update(str(contiguous.shape).encode("ascii"))
    digest.update(contiguous.view(np.uint8))
    return digest.hexdigest()


def load_inputs():
    for path, expected in SEALED.items():
        if sha256(path) != expected:
            raise RuntimeError(f"sealed source changed: {path.name}")
    with np.load(ORIGIN, allow_pickle=False) as origin, np.load(EPISODE, allow_pickle=False) as episode:
        original = np.vstack((origin["nf512_initial_phi0"], origin["nf512_initial_phi1"]))
        observer = np.array(original[:, :2], copy=True)
        columns = np.vstack((episode[f"{CASE}_frame_phi0"][0], episode[f"{CASE}_frame_phi1"][0]))
        handoff = np.vstack((origin[f"{HANDOFF}_final_phi0"], origin[f"{HANDOFF}_final_phi1"]))
        weights = np.array(origin["nf512_occupations"], copy=True)
        if not np.array_equal(columns, handoff) or not np.array_equal(weights, WEIGHTS):
            raise RuntimeError("the original weighted handoff is required")
        times = np.array(episode[f"{CASE}_frame_time"], copy=True)
        geometry = np.array(episode[f"{CASE}_frame_quad_Q"], copy=True)
        coarse = np.array(episode[f"{CASE}_frame_coarse_Q"][0], copy=True)
        autonomous = np.concatenate((episode[f"{CASE}_frame_phi0"], episode[f"{CASE}_frame_phi1"]), axis=1)
        primary_phi0 = np.array(original[: original.shape[0] // 2], copy=True)
    with np.load(RESPONSE, allow_pickle=False) as response:
        arrays = {name: np.array(response[name], copy=True) for name in response.files}
    record = json.loads(RESPONSE_JSON.read_text())
    return {
        "observer": observer, "columns": columns, "weights": weights, "times": times,
        "geometry": geometry, "coarse": coarse, "phi0": primary_phi0,
        "autonomous": autonomous, "saved": arrays, "record": record,
    }


def saved_consistency_errors(inputs):
    """Independent checks of numerical quantities, including the full cross diagonal."""
    errors = []
    saved, record = inputs["saved"], inputs["record"]
    observer, columns = inputs["observer"], inputs["columns"]
    for source_key, request_key in (("observer", "observer"), ("columns", "initial_state"),
                                    ("geometry", "geometry"), ("weights", "weights")):
        if array_id(inputs[source_key]) != record["request"][request_key]:
            errors.append("binding_" + request_key)
    if not np.array_equal(inputs["weights"], WEIGHTS):
        errors.append("original_weights")
    if not np.array_equal(saved["times"], inputs["times"]):
        errors.append("full_times")
    if not np.allclose(observer.conj().T @ observer, np.eye(2), rtol=0.0, atol=1e-13):
        errors.append("observer_gram")
    projected = np.einsum("ia,tik->tak", observer.conj(), inputs["autonomous"])
    if max_abs(occupation(covariance(projected)) - saved["occupation_autonomous"]) > 1e-12:
        errors.append("fixed_original_observer")
    initial = observer.conj().T @ columns
    exterior = columns - observer @ initial
    cross_norm = float(np.linalg.norm((initial * WEIGHTS) @ exterior.conj().T))
    if abs(cross_norm - record["initial_blocks"]["C_AE_frobenius"]) > 1e-12:
        errors.append("actual_initial_cross")
    panel = record["phases"]["full_controls"]["omissions"]
    mixed = saved["covariance_retained"] - saved["covariance_without_cross"]
    metrics = {
        "occupation_diagonal_movement": max_abs(occupation(mixed)),
        "separation_max_frobenius": float(np.max(np.linalg.norm(mixed, axis=(1, 2)))),
    }
    for key, measured in metrics.items():
        if abs(measured - panel["initial_cross"][key]) > 1e-12:
            errors.append("full_cross_" + key)
    if max_abs(occupation(saved["covariance_retained"]) - saved["occupation_retained"]) > 1e-12:
        errors.append("cross_baseline")
    probe = record["phases"]["probe_controls"]["omissions"]["initial_cross"]
    if abs(max_abs(saved["probe_cross_occupation_movement"]) - probe["occupation_diagonal_movement"]) > 1e-12:
        errors.append("probe_cross_diagonal")
    for name in ("memory", "outside_drive"):
        measured = max_abs(saved["occupation_" + name + "_off"] - saved["occupation_full"])
        if abs(measured - panel[name]["occupation_error_against_reference"]["max_abs"]) > 1e-12:
            errors.append(name + "_omission")
    return errors


def replay(inputs, stop_nodes=None, deadline=None):
    """Independent full split and memoryless trapezoid, on source columns only.

    P_E H P_E is represented in the original space, with two zero local
    directions. Its exponentials act on six exterior columns. This avoids a
    dense exterior frame and every time-indexed propagator.
    """
    times = inputs["times"] if stop_nodes is None else inputs["times"][:stop_nodes]
    observer, columns = inputs["observer"], inputs["columns"]
    grid, factors = build_case_grid({"phi0": inputs["phi0"], "frame_Q": inputs["geometry"], "coarse_Q": inputs["coarse"]})
    del grid

    def hamiltonian(mark):
        index = min(int(np.searchsorted(inputs["times"], mark, side="right") - 1), inputs["times"].size - 2)
        fraction = (mark - inputs["times"][index]) / (inputs["times"][index + 1] - inputs["times"][index])
        radial = (1.0 - fraction) * inputs["geometry"][index] + fraction * inputs["geometry"][index + 1]
        return galerkin_matrix(factors, radial)

    retained_initial = observer @ (observer.conj().T @ columns)
    exterior_initial = columns - retained_initial
    full = np.concatenate((retained_initial, exterior_initial), axis=1)
    free_exterior = np.array(exterior_initial, copy=True)
    memoryless = observer.conj().T @ columns
    local_full = np.empty((times.size, 2, 12), dtype=complex)
    local_memoryless = np.empty((times.size, 2, 6), dtype=complex)
    local_full[0] = observer.conj().T @ full
    local_memoryless[0] = memoryless
    h_now = hamiltonian(float(times[0]))
    exterior_projection_defect = 0.0
    for index, width in enumerate(np.diff(times)):
        if deadline is not None and time.process_time() > deadline:
            raise RuntimeError("independent control replay exhausted its CPU budget")
        a_now = observer.conj().T @ h_now @ observer
        drive_now = observer.conj().T @ h_now @ free_exterior
        for part in range(2):
            midpoint = float(times[index] + (part + 0.5) * width / 2.0)
            h_mid = hamiltonian(midpoint)
            coefficient = -1j * width / 2.0
            full = expm_multiply(coefficient * h_mid, full, traceA=coefficient * np.trace(h_mid))
            hv = h_mid @ observer
            local = observer.conj().T @ hv
            exterior = h_mid - observer @ hv.conj().T - hv @ observer.conj().T + observer @ local @ observer.conj().T
            free_exterior = expm_multiply(coefficient * exterior, free_exterior, traceA=coefficient * (np.trace(h_mid) - np.trace(local)))
            exterior_projection_defect = max(exterior_projection_defect, max_abs(observer.conj().T @ free_exterior))
        h_new = hamiltonian(float(times[index + 1]))
        a_new = observer.conj().T @ h_new @ observer
        drive_new = observer.conj().T @ h_new @ free_exterior
        right = memoryless - 1j * width / 2.0 * (a_now @ memoryless + drive_now + drive_new)
        memoryless = np.linalg.solve(np.eye(2) + 1j * width / 2.0 * a_new, right)
        local_full[index + 1] = observer.conj().T @ full
        local_memoryless[index + 1] = memoryless
        h_now = h_new
    a, e = local_full[:, :, :6], local_full[:, :, 6:]
    covariance_full = covariance(a + e)
    covariance_without_cross = covariance(a) + covariance(e)
    return {
        "times": times, "covariance_full": covariance_full,
        "covariance_without_cross": covariance_without_cross,
        "covariance_cross": covariance_full - covariance_without_cross,
        "occupation_full": occupation(covariance_full),
        "occupation_drive_off": occupation(covariance(a)),
        "occupation_memory_off": occupation(covariance(local_memoryless)),
        "exterior_projection_defect": exterior_projection_defect,
        "column_width": 12, "free_exterior_width": 6,
        "dense_propagator_history_bytes": 0,
        "full_column_workspace_bytes": int(full.nbytes),
    }


def execute(budget_s=300.0):
    started = time.process_time()
    inputs = load_inputs()
    errors = saved_consistency_errors(inputs)
    if errors:
        raise RuntimeError("saved control consistency: " + ", ".join(errors))
    pilot_started = time.process_time()
    replay(inputs, stop_nodes=3, deadline=started + budget_s)
    pilot_cpu = time.process_time() - pilot_started
    # Include a second grid build and give the measured two-step cost linear headroom.
    forecast = pilot_cpu * 4.5
    if time.process_time() - started + forecast > budget_s:
        raise RuntimeError(f"forecast {forecast:g}s exceeds remaining CPU budget")
    print(f"independent pilot {pilot_cpu:.6f}s, full forecast {forecast:.6f}s", flush=True)
    full = replay(inputs, deadline=started + budget_s)
    saved = inputs["saved"]
    saved_cross = saved["covariance_retained"] - saved["covariance_without_cross"]
    metrics = {
        "independent_full_versus_saved_full_occupation": max_abs(full["occupation_full"] - saved["occupation_full"]),
        "streamed_cross_versus_independent_full_cross_occupation": max_abs(occupation(saved_cross) - occupation(full["covariance_cross"])),
        "streamed_cross_versus_independent_full_cross_frobenius": float(np.max(np.linalg.norm(saved_cross - full["covariance_cross"], axis=(1, 2)))),
        "streamed_drive_off_versus_independent_full_drive_off": max_abs(saved["occupation_outside_drive_off"] - full["occupation_drive_off"]),
        "independent_memoryless_versus_saved_memoryless": max_abs(full["occupation_memory_off"] - saved["occupation_memory_off"]),
        "independent_full_cross_occupation_movement": max_abs(occupation(full["covariance_cross"])),
        "saved_full_cross_occupation_movement": max_abs(occupation(saved_cross)),
        "saved_probe_cross_occupation_movement": max_abs(saved["probe_cross_occupation_movement"]),
        "independent_full_drive_occupation_movement": max_abs(full["occupation_drive_off"] - full["occupation_full"]),
        "independent_memory_occupation_movement": max_abs(full["occupation_memory_off"] - full["occupation_full"]),
    }
    arrays = {key: value for key, value in full.items() if isinstance(value, np.ndarray)}
    arrays["saved_full_cross_occupation_movement"] = occupation(saved_cross)
    record = {
        "schema": "NSC-COUPLED-LOCAL-RESPONSE-AUDIT-v1", "status": "MEASURED_INDEPENDENT_CONTROLS",
        "case": CASE, "time_window": full["times"][[0, -1]].tolist(), "stored_frames": int(full["times"].size),
        "observer": "original T=0 region-0 pair, columns 0 and 1, no QR or phase reset",
        "preparation": "bitwise transported Phi(0.05); original Gaussian weights", "weights": WEIGHTS.tolist(),
        "geometry": "same stored piecewise-linear quadrature Q with static lapse and shift",
        "midpoint_substeps": 2, "memoryless_method": "independent projected-space exterior transport and local trapezoid",
        "reused_owner": "checked Fourier-Galerkin matrix representation only",
        "temporal_consumer_helpers_called": False, "production_rk4_rerun": False,
        "source_bindings": {str(path.relative_to(LAB)): digest for path, digest in SEALED.items()},
        "saved_consistency_errors": errors, "metrics": metrics,
        "error_scope": "cross superposition roundoff is an algebraic check; cross reduction error is measured against independent full split in occupation or covariance units",
        "pilot_cpu_seconds": pilot_cpu, "full_forecast_seconds": forecast,
        "cpu_seconds": time.process_time() - started, "cpu_budget_seconds": budget_s,
        "allocation": {key: value for key, value in full.items() if not isinstance(value, np.ndarray)},
        "code_sha256": sha256(__file__), "renewal_claimed": False, "stress_claimed": False,
    }
    np.savez_compressed(OUTPUT_NPZ, **arrays)
    record["payload_npz_sha256"] = sha256(OUTPUT_NPZ)
    OUTPUT.write_text(json.dumps(record, indent=2) + "\n")
    for _ in range(3):
        record["payload_bytes"] = OUTPUT.stat().st_size + OUTPUT_NPZ.stat().st_size
        OUTPUT.write_text(json.dumps(record, indent=2) + "\n")
    if record["payload_bytes"] > 64 * 1024 * 1024:
        raise RuntimeError("supplemental payload exceeds 64 MiB")
    for path, expected in SEALED.items():
        if sha256(path) != expected:
            raise RuntimeError(f"audit changed sealed source: {path.name}")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-s", type=float, default=300.0)
    args = parser.parse_args()
    result = execute(args.budget_s)
    print(json.dumps(result["metrics"], indent=2))
    print(f"CPU {result['cpu_seconds']:.6f}s; payload {result['payload_bytes']} bytes")
