"""Causal photon response of one paired-phase tangent on the original conformal ODE.

The baseline state and its initial tangent advance together to T=0.3.
Redshift derivatives use the fixed null events of the longest saved ray.
This does not rerun a root solve, change the action, or identify H0.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

from recursive_horizons.nsc_conformal_preparation_tangent import preparation_tangent
from recursive_horizons.nsc_conformal_tangent_evolution import (
    FIELDS, conformal_jvp, conformal_tangent_constraints, variational_rk4_step,
)
from recursive_horizons.nsc_spherical_coupling import (
    PERIOD, TOL_FD_RELATIVE, CauchyState, PositiveChartExit,
)
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    build_grid, compose_fine_hamiltonian, conformal_constraint_transport,
    constraint_diagnostics,
)
from measure_nsc_conformal_geometric_average import (
    CURVATURE, EPISODE, EXPECTED_PAYLOAD, sha256,
)
from measure_nsc_conformal_preparation_tangent import (
    GEOMETRY_NAMES, PREFIX, STATE_NAMES, _archive_state,
)


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-conformal-causal-response-summary-v1"
COMPLETED = "completed_causal_photon_response"
PARTIAL = "partial_causal_photon_response"
SETUP_CPU = 140.0
COARSE_CPU = 1000.0
FINE_CPU = 2000.0
TOTAL_CPU = 3600.0
BYTE_LIMIT = 24_000
RAW_BYTE_LIMIT = 2_000_000
SETUP_THREADS = 8
WORKER_THREADS = 4
WORKERS = 2
PREFLIGHT_H = 1e-4
BASELINE_TOLERANCE = 1e-8
LIGHTCONE_Z_TOLERANCE = 1e-8
TARGET_TIME = 0.3
DIAGNOSTIC_TIMES = (0.0, 0.05, 0.175, 0.3)
RECEPTION_TIME = 0.3
RECEPTION_X = 1.0
EMISSION_TIME = 0.175
LOOKBACK = 0.125
NODE_TOLERANCE = 1e-12
REFERENCE_FRACTION = -0.27203771374589336
CASES = (
    {"name": "nf512_dt_0.0010", "dt": 0.001, "cpu_seconds": COARSE_CPU},
    {"name": "nf512_dt_0.0005", "dt": 0.0005, "cpu_seconds": FINE_CPU},
)
SOURCE_PATHS = (
    "scripts/measure_nsc_conformal_causal_response.py",
    "tests/test_measure_nsc_conformal_causal_response.py",
    "scripts/measure_nsc_conformal_preparation_tangent.py",
    "scripts/measure_nsc_conformal_geometric_average.py",
    "src/recursive_horizons/nsc_conformal_tangent_evolution.py",
    "src/recursive_horizons/nsc_conformal_preparation_tangent.py",
    "src/recursive_horizons/nsc_spherical_coupling.py",
    "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "src/recursive_horizons/nsc_spherical_feedback_action.py",
    "src/recursive_horizons/nsc_covariant_operator.py",
    "src/recursive_horizons/nsc_conformal_adm_source.py",
    "src/recursive_horizons/nsc_atlas_cosmological_distance.py",
)
LIGHTCONE_SUMMARY = ROOT / "results" / "summaries" / "nsc-conformal-lightcone-v1.json"
DEFAULT_GEOMETRY = ROOT / ".build" / "data" / "nsc-spherical-conformal-episode-v2.npz"
DEFAULT_RAW = ROOT / ".build" / "results" / "nsc-conformal-causal-response-v1.npz"
THREAD_VARIABLES = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
)


class CpuAllocationExceeded(RuntimeError):
    """A stage used more process CPU than its allocation. There is no retry."""

    def __init__(self, stage, used, limit):
        super().__init__(
            f"{stage} exceeded {limit:g} process-CPU seconds ({used:.3f}); no automatic retry"
        )
        self.stage = stage
        self.used = float(used)
        self.limit = float(limit)


class CpuClock:
    """process_time checkpoint with the same role as the preparation measure budget."""

    def __init__(self, seconds, started=None):
        self.seconds = float(seconds)
        self.started = time.process_time() if started is None else float(started)

    @property
    def used(self):
        return time.process_time() - self.started

    def check(self, stage):
        used = self.used
        if used > self.seconds:
            raise CpuAllocationExceeded(stage, used, self.seconds)
        return used


def require_fresh(path, kind):
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"fresh {kind} path already exists: {path}")
    if not path.parent.is_dir():
        raise ValueError(f"{kind} parent must exist: {path.parent}")
    return path


def quadrature_index(xi_q, position):
    """Node of a uniform quadrature grid. This is not a coarse geometry index."""
    nodes = np.asarray(xi_q, dtype=float)
    position = float(position)
    if nodes.ndim != 1 or nodes.size < 2 or not np.isfinite(nodes).all():
        raise ValueError("quadrature nodes must be a finite one-dimensional grid")
    index = int(np.argmin(np.abs(nodes - position)))
    if abs(float(nodes[index]) - position) > NODE_TOLERANCE:
        raise ValueError(f"x={position} is not a quadrature node")
    return index


def prolonged_value(row, values):
    """One row of A_g applied to a coarse geometry field."""
    return float(np.asarray(row, dtype=float) @ np.asarray(values, dtype=float))


def production_nodes(xi_q):
    receiver = quadrature_index(xi_q, RECEPTION_X)
    emitters = {
        -1: quadrature_index(xi_q, RECEPTION_X - LOOKBACK),
        1: quadrature_index(xi_q, RECEPTION_X + LOOKBACK),
    }
    if len(xi_q) == 2048 and (receiver, emitters[-1], emitters[1]) != (256, 224, 288):
        raise ValueError("production quadrature nodes are not 256, 224 and 288")
    return {"receiver": receiver, "emitters": emitters}


def log_redshift_derivative(receiver, emitter):
    """d log(1+z) = d log N_receiver - d log N_emitter at fixed events."""
    for sample in (receiver, emitter):
        if sample["r"] <= 0.0 or sample["Q"] <= 0.0:
            raise ValueError("redshift sample left the positive chart")
    receiver_dlog = receiver["ur"] / receiver["r"] + receiver["uQ"] / receiver["Q"]
    emitter_dlog = emitter["ur"] / emitter["r"] + emitter["uQ"] / emitter["Q"]
    return receiver_dlog - emitter_dlog


def event_sample(grid, state, tangent, index, time_value, position):
    radius = prolonged_value(grid.A_g[index], state.r)
    scale = prolonged_value(grid.A_g[index], state.Q)
    ur = prolonged_value(grid.A_g[index], tangent.r)
    uq = prolonged_value(grid.A_g[index], tangent.Q)
    if radius <= 0.0 or scale <= 0.0:
        raise PositiveChartExit("event_lapse_left_positive_chart", time_value, state.copy())
    return {
        "time": float(time_value),
        "x": float(position),
        "quadrature_index": int(index),
        "r": radius,
        "Q": scale,
        "N": radius * scale,
        "ur": ur,
        "uQ": uq,
        "dlogN": ur / radius + uq / scale,
    }


def ray_row(direction, receiver, emitter):
    derivative = log_redshift_derivative(receiver, emitter)
    one_plus = receiver["N"] / emitter["N"]
    return {
        "direction": int(direction),
        "coordinate_lookback": LOOKBACK,
        "reception_time": receiver["time"],
        "reception_x": receiver["x"],
        "emission_time": emitter["time"],
        "emission_x": emitter["x"],
        "receiver_N": receiver["N"],
        "emitter_N": emitter["N"],
        "one_plus_redshift": one_plus,
        "redshift": one_plus - 1.0,
        "receiver_r": receiver["r"],
        "receiver_Q": receiver["Q"],
        "receiver_ur": receiver["ur"],
        "receiver_uQ": receiver["uQ"],
        "emitter_r": emitter["r"],
        "emitter_Q": emitter["Q"],
        "emitter_ur": emitter["ur"],
        "emitter_uQ": emitter["uQ"],
        "receiver_dlogN": receiver["dlogN"],
        "emitter_dlogN": emitter["dlogN"],
        "log_redshift_derivative": derivative,
        "events_fixed": True,
        "tangent_interpolation": False,
    }


def rays_from_events(receiver, emitters):
    if receiver is None:
        return []
    rows = []
    for direction in (-1, 1):
        emitter = emitters.get(direction)
        if emitter is not None:
            rows.append(ray_row(direction, receiver, emitter))
    return rows


def forecast_projection(elapsed, step_cost, steps_done, steps_required):
    """CPU to finish the remaining steps at the measured first-step cost."""
    if steps_done > steps_required or step_cost < 0.0 or elapsed < 0.0:
        raise ValueError("forecast counts must be ordered and nonnegative")
    return float(elapsed + step_cost * (steps_required - steps_done))


def forecast_accepts(elapsed, step_cost, steps_done, steps_required, limit):
    projected = forecast_projection(elapsed, step_cost, steps_done, steps_required)
    return projected <= float(limit), projected


def energy_variation(grid, rate, tangent):
    """Omega(X, u): coarse dx_g on geometry and 2 M c Im on the columns."""
    value = grid.dx_g * sum(
        np.sum(getattr(rate, coordinate) * getattr(tangent, momentum)
               - getattr(rate, momentum) * getattr(tangent, coordinate))
        for coordinate, momentum in (("Q", "p_Q"), ("r", "p_r"), ("chi", "p_chi"))
    )
    weights = grid.fine.occupations[None, :]
    for name in ("phi0", "phi1"):
        value += 2.0 * grid.fine.multiplicity * np.sum(
            weights * np.imag(getattr(rate, name).conj() * getattr(tangent, name))
        )
    return float(value)


def column_gram(phi0, phi1):
    return phi0.conj().T @ phi0 + phi1.conj().T @ phi1


def tangent_gram(state, tangent):
    return (
        tangent.phi0.conj().T @ state.phi0 + state.phi0.conj().T @ tangent.phi0
        + tangent.phi1.conj().T @ state.phi1 + state.phi1.conj().T @ tangent.phi1
    )


def diagnostic_record(grid, state, tangent, time_value, steps, energy_reference):
    rate, bundle = compose_fine_hamiltonian(grid, state)
    transport = conformal_constraint_transport(grid, state, rate, bundle)
    tangent_report = conformal_tangent_constraints(
        grid, state, tangent, rate=rate, bundle=bundle,
    )
    full = constraint_diagnostics(grid, state, bundle["source"], bundle["fine_state"])
    variation = energy_variation(grid, rate, tangent)
    reference = variation if energy_reference is None else energy_reference
    return {
        "time": float(time_value),
        "steps": int(steps),
        "full_hamilton_max": float(full["full_hamilton_max"]),
        "full_momentum_max": float(full["full_momentum_max"]),
        "constraint_norm": float(transport["constraint_norm"]),
        "forcing_norm": float(transport["forcing_norm"]),
        "quadrature_forcing_norm": float(transport["quadrature_forcing_norm"]),
        "projection_forcing_norm": float(transport["projection_forcing_norm"]),
        "constraint_tangent_norm": float(tangent_report["constraint_tangent_norm"]),
        "forcing_tangent_norm": float(tangent_report["forcing_tangent_norm"]),
        "hamilton_tangent_max_abs": float(tangent_report["hamilton_tangent_max_abs"]),
        "momentum_tangent_max_abs": float(tangent_report["momentum_tangent_max_abs"]),
        "shift_source_tangent_max_abs": float(tangent_report["shift_source_tangent_max_abs"]),
        "gram_error_max_abs": float(np.max(np.abs(column_gram(state.phi0, state.phi1) - np.eye(6)))),
        "tangent_gram_error_max_abs": float(np.max(np.abs(tangent_gram(state, tangent)))),
        "energy_variation": variation,
        "energy_variation_drift": variation - reference,
        "u_max_abs": [
            float(np.max(np.abs(getattr(tangent, name)))) for name in FIELDS
        ],
    }, reference


def offset_state(state, tangent, step):
    return CauchyState(**{
        name: getattr(state, name) + step * getattr(tangent, name) for name in FIELDS
    })


def block_relative(exact, approximate):
    gap = float(np.max(np.abs(np.asarray(exact) - np.asarray(approximate))))
    scale = max(1.0, float(np.max(np.abs(exact))), float(np.max(np.abs(approximate))))
    return {"max_abs": gap, "scale": scale, "relative": gap / scale}


def preflight_jvp(grid, state, tangent, h=PREFLIGHT_H):
    """One central difference of the full rate. Two compose calls, fixed h."""
    step = float(h)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("preflight step must be finite and positive")
    plus, _bundle_plus = compose_fine_hamiltonian(grid, offset_state(state, tangent, step))
    minus, _bundle_minus = compose_fine_hamiltonian(grid, offset_state(state, tangent, -step))
    exact = conformal_jvp(grid, state, tangent)
    blocks = []
    for name in FIELDS:
        finite_difference = (getattr(plus, name) - getattr(minus, name)) / (2.0 * step)
        item = block_relative(getattr(exact, name), finite_difference)
        item["field"] = name
        blocks.append(item)
    return {
        "h": step,
        "compose_calls": 2,
        "gate": float(TOL_FD_RELATIVE),
        "accepted": all(item["relative"] <= TOL_FD_RELATIVE for item in blocks),
        "blocks": blocks,
        "repeated_step_search": False,
        "dense_jacobian": False,
    }


def pack_state(state):
    return {name: np.array(getattr(state, name), copy=True) for name in STATE_NAMES}


def evolve_case(grid, state, tangent, *, dt, target, cpu_seconds, name,
                diagnostic_times=DIAGNOSTIC_TIMES, sample_events=True,
                checkpoint_times=(0.05, 0.3), started=None):
    """Advance z and u with variational RK4. Stop rather than shorten the target."""
    clock = CpuClock(cpu_seconds, started)
    step = float(dt)
    target = float(target)
    required = int(round(target / step))
    if required < 1 or abs(required * step - target) > 1e-9:
        raise ValueError("target is not an integer number of accepted time steps")
    diagnostic_at = {}
    for mark in diagnostic_times:
        count = int(round(float(mark) / step))
        if count < 0 or count > required or abs(count * step - float(mark)) > 1e-9:
            raise ValueError("diagnostic time is not on this step sequence")
        diagnostic_at[count] = float(mark)
    nodes = production_nodes(grid.xi_q) if sample_events else None
    checkpoints = {float(mark) for mark in checkpoint_times}
    diagnostics = []
    events = {"receiver": None, "emitters": {}}
    stored = {}
    energy_reference = None
    forecast = None
    reason = None
    steps = 0
    time_value = 0.0
    last_step = None

    def record(mark):
        nonlocal energy_reference
        row, energy_reference = diagnostic_record(
            grid, state, tangent, mark, steps, energy_reference,
        )
        diagnostics.append(row)

    def remember(mark):
        if mark in checkpoints:
            stored[f"{mark:.3f}"] = pack_state(state)
        if nodes is None:
            return
        if abs(mark - EMISSION_TIME) <= 1e-12:
            for direction, index in nodes["emitters"].items():
                position = RECEPTION_X + direction * LOOKBACK
                events["emitters"][direction] = event_sample(
                    grid, state, tangent, index, mark, position,
                )
        if abs(mark - RECEPTION_TIME) <= 1e-12:
            events["receiver"] = event_sample(
                grid, state, tangent, nodes["receiver"], mark, RECEPTION_X,
            )

    try:
        clock.check("before initial diagnostic")
        if 0 in diagnostic_at:
            record(diagnostic_at[0])
            remember(diagnostic_at[0])
        clock.check("after initial diagnostic")
        while steps < required:
            if last_step is not None and clock.used + last_step > clock.seconds:
                reason = "cpu_allocation_exhausted"
                break
            clock.check("before step")
            step_started = time.process_time()
            state, tangent = variational_rk4_step(
                grid, state, tangent, step, checkpoint=lambda: clock.check("rk4 stage"),
            )
            last_step = time.process_time() - step_started
            steps += 1
            time_value = steps * step
            if steps in diagnostic_at:
                record(diagnostic_at[steps])
                remember(diagnostic_at[steps])
            if steps == 1:
                accepted, projected = forecast_accepts(
                    clock.used, last_step, 1, required, clock.seconds,
                )
                forecast = {
                    "first_step_cpu": float(last_step),
                    "elapsed_after_first_step": float(clock.used),
                    "projected_cpu": projected,
                    "required_steps": required,
                    "accepted": bool(accepted),
                    "endpoint_shortened": False,
                }
                if not accepted:
                    reason = "forecast_cannot_reach_full_endpoint"
                    break
    except CpuAllocationExceeded:
        reason = "cpu_allocation_exhausted"
    except PositiveChartExit as error:
        reason = "positive_chart:" + str(error.reason)
    reached = reason is None and steps == required and abs(time_value - target) <= 1e-9
    rays = rays_from_events(events["receiver"], events["emitters"]) if sample_events else []
    report = {
        "name": name,
        "dt": step,
        "nf": int(grid.nf),
        "ng": int(grid.ng),
        "nq": int(grid.nq),
        "target": target,
        "time": float(time_value),
        "steps": int(steps),
        "required_steps": required,
        "target_reached": bool(reached),
        "reason": reason,
        "cpu_seconds": float(clock.used),
        "forecast": forecast,
        "diagnostics": diagnostics,
        "rays": rays,
        "photon_result": bool(reached and len(rays) == 2),
        "energy_variation_is_hard_gate": False,
        "continuum_constraint_certified": False,
        "uniform_stage_constraint_claim": False,
        "u_max_abs_fields": list(FIELDS),
        "constraint_projection": False,
        "chirality_projection": False,
        "endpoint_shortened": False,
    }
    cache = {
        "name": name,
        "time": float(time_value),
        "target_reached": bool(reached),
        "state": pack_state(state),
        "tangent": pack_state(tangent),
        "checkpoints": stored,
    }
    return plain(report), cache


def field_scale_gap(actual, reference, tolerance=BASELINE_TOLERANCE):
    actual = np.asarray(actual)
    reference = np.asarray(reference)
    if actual.shape != reference.shape:
        raise ValueError("baseline fields do not have the same shape")
    gap = float(np.max(np.abs(actual - reference)))
    scale = max(1.0, float(np.max(np.abs(actual))), float(np.max(np.abs(reference))))
    return {"max_abs": gap, "scale": scale, "within_tolerance": bool(gap <= tolerance * scale)}


def archived_frame(data, case, time_value):
    times = np.asarray(data[case + "_frame_times"], dtype=float)
    hits = np.flatnonzero(np.abs(times - float(time_value)) <= 1e-9)
    if hits.size != 1:
        raise ValueError(f"{case} has no unique saved frame at {time_value}")
    index = int(hits[0])
    return {
        name: np.array(data[case + "_frames_" + name][index], copy=True)
        for name in STATE_NAMES
    }


def compare_checkpoint(actual, reference, tolerance=BASELINE_TOLERANCE):
    fields = {
        name: field_scale_gap(actual[name], reference[name], tolerance) for name in STATE_NAMES
    }
    failing = [name for name, item in fields.items() if not item["within_tolerance"]]
    return {
        "max_abs": {name: item["max_abs"] for name, item in fields.items()},
        "scale": {name: item["scale"] for name, item in fields.items()},
        "failing_fields": failing,
        "within_tolerance": not failing,
        "tolerance": tolerance,
        "bitwise_gate": False,
    }


def lightcone_redshifts(path):
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    found = {}
    for row in record.get("rows", []):
        observables = row.get("observables", {})
        if (row.get("case") in {case["name"] for case in CASES}
                and abs(float(observables.get("coordinate_lookback", -1)) - LOOKBACK) < 1e-12):
            found[(row["case"], int(observables["direction"]))] = float(observables["redshift"])
    return found


def redshift_metadata(rays, case, references, tolerance=LIGHTCONE_Z_TOLERANCE):
    rows = []
    for ray in rays:
        key = (case, int(ray["direction"]))
        if key not in references:
            rows.append({"direction": ray["direction"], "available": False})
            continue
        gap = float(ray["redshift"] - references[key])
        rows.append({
            "direction": int(ray["direction"]),
            "available": True,
            "saved_redshift": references[key],
            "computed_redshift": float(ray["redshift"]),
            "difference": gap,
            "within_tolerance": bool(abs(gap) <= tolerance),
            "bounds_tangent": False,
        })
    return {"tolerance": tolerance, "rows": rows}


def temporal_envelope(reports):
    by_name = {report["name"]: report for report in reports}
    coarse = by_name.get("nf512_dt_0.0010")
    fine = by_name.get("nf512_dt_0.0005")
    if coarse is None or fine is None or not coarse.get("photon_result") or not fine.get("photon_result"):
        return None
    rows = []
    for direction in (-1, 1):
        coarse_ray = next(ray for ray in coarse["rays"] if ray["direction"] == direction)
        fine_ray = next(ray for ray in fine["rays"] if ray["direction"] == direction)
        fine_value = float(fine_ray["log_redshift_derivative"])
        coarse_value = float(coarse_ray["log_redshift_derivative"])
        rows.append({
            "direction": direction,
            "fine_coefficient": fine_value,
            "coarse_coefficient": coarse_value,
            "absolute_difference": abs(fine_value - coarse_value),
        })
    return {
        "quantity": "derivative of log(1+z) per epsilon",
        "fine_dt": 0.0005,
        "coarse_dt": 0.001,
        "rows": rows,
        "empirical_temporal_envelope": True,
        "continuum_bound": False,
        "spatial_bound": False,
        "physical_error_bound": False,
        "convergence_claimed": False,
    }


def case_succeeded(report):
    return bool(report.get("target_reached") and report.get("photon_result") and len(report.get("rays", [])) == 2)


def summary_status(reports):
    named = {report.get("name") for report in reports}
    if named == {case["name"] for case in CASES} and all(case_succeeded(report) for report in reports):
        return COMPLETED
    return PARTIAL


def plain(value):
    if isinstance(value, dict):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if isinstance(value, np.ndarray):
        raise TypeError("arrays belong in the raw cache, not the summary")
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def dump_summary(payload):
    raw = json.dumps(plain(payload), indent=2, sort_keys=True, allow_nan=False) + "\n"
    size = len(raw.encode("utf-8"))
    if size > BYTE_LIMIT:
        raise ValueError(f"summary exceeds {BYTE_LIMIT} bytes ({size})")
    return raw


def observed_threads():
    return [item.get("num_threads") for item in threadpool_info()]


def run_case(payload):
    """One numerical process. Rebuild the grid here; receive z and u as arrays."""
    threads = int(payload["threads"])
    for variable in THREAD_VARIABLES:
        os.environ[variable] = str(threads)
    started = time.process_time()
    try:
        with threadpool_limits(limits=threads):
            threads_seen = observed_threads()
            grid = build_grid(int(payload["nf"]), quadrature=int(payload["nq"]), gauge="conformal")
            state = CauchyState(**{
                name: np.array(payload["state"][name], copy=True) for name in STATE_NAMES
            })
            tangent = CauchyState(**{
                name: np.array(payload["tangent"][name], copy=True) for name in STATE_NAMES
            })
            report, cache = evolve_case(
                grid, state, tangent,
                dt=float(payload["dt"]),
                target=float(payload["target"]),
                cpu_seconds=float(payload["cpu_seconds"]),
                name=payload["name"],
                diagnostic_times=tuple(payload["diagnostic_times"]),
                sample_events=bool(payload["sample_events"]),
                checkpoint_times=tuple(payload["checkpoint_times"]),
                started=started,
            )
            report["observed_num_threads"] = threads_seen
            report["requested_threads"] = threads
            return {"report": report, "cache": cache}
    except (Exception, PositiveChartExit) as error:
        return {"report": {
            "name": payload.get("name"),
            "dt": payload.get("dt"),
            "target": payload.get("target"),
            "time": 0.0,
            "steps": 0,
            "target_reached": False,
            "reason": f"{type(error).__name__}: {error}",
            "cpu_seconds": float(time.process_time() - started),
            "diagnostics": [],
            "events": {"receiver": None, "emitters": []},
            "rays": [],
            "photon_result": False,
            "forecast": None,
        }, "cache": None}


def collect_cases(payloads, workers=WORKERS):
    if len(payloads) != workers:
        raise ValueError("this producer launches one worker per retained time step")
    context = multiprocessing.get_context("spawn")
    results = []
    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        futures = {pool.submit(run_case, payload): payload for payload in payloads}
        for future, payload in futures.items():
            try:
                results.append(future.result())
            except Exception as error:
                results.append({"report": {
                    "name": payload.get("name"),
                    "target_reached": False,
                    "reason": f"{type(error).__name__}: {error}",
                    "rays": [],
                    "photon_result": False,
                    "time": 0.0,
                    "steps": 0,
                }, "cache": None})
    return results


def initial_tangent(prepared):
    """Geometry tangents other than the solved radius response stay zero at T=0."""
    return CauchyState(
        Q=np.zeros_like(prepared.delta_r),
        r=np.array(prepared.delta_r, dtype=float, copy=True),
        chi=np.zeros_like(prepared.delta_r),
        p_Q=np.zeros_like(prepared.delta_r),
        p_r=np.zeros_like(prepared.delta_r),
        p_chi=np.zeros_like(prepared.delta_r),
        phi0=np.array(prepared.phi0_tangent, dtype=complex, copy=True),
        phi1=np.array(prepared.phi1_tangent, dtype=complex, copy=True),
    )


def frozen_sources():
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
    ).strip()
    hashes = {}
    for name in SOURCE_PATHS:
        raw = (ROOT / name).read_bytes()
        try:
            committed = subprocess.check_output(["git", "show", f"{commit}:{name}"], cwd=ROOT)
        except subprocess.CalledProcessError as error:
            raise ValueError("causal response must be frozen before application: " + name) from error
        if raw != committed:
            raise ValueError("causal response must be frozen before application: " + name)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    return commit, hashes


def authenticate(geometry):
    commit, bindings = frozen_sources()
    old = json.loads(CURVATURE.read_text(encoding="utf-8"))["source_bindings"]
    episode = json.loads(EPISODE.read_text(encoding="utf-8"))
    if (
        sha256(geometry) != EXPECTED_PAYLOAD
        or sha256(EPISODE) != old["episode_v2"]
        or episode["payload_sha256"] != EXPECTED_PAYLOAD
        or episode["gauge"] != "conformal"
        or episode["all_target_reached"] is not True
        or episode["source_changed"] is not False
        or episode["midtrajectory_reset"] is not False
    ):
        raise ValueError("original conformal episode is not authenticated")
    for key, path in (
        ("galerkin", "src/recursive_horizons/nsc_spherical_galerkin_coupling.py"),
        ("action", "src/recursive_horizons/nsc_spherical_feedback_action.py"),
    ):
        if bindings[path] != old[key]:
            raise ValueError("original action or projection changed")
    return commit, bindings, episode


def write_raw(path, caches):
    payload = {}
    for cache in caches:
        if not cache or cache.get("state") is None or cache.get("tangent") is None:
            continue
        prefix = cache["name"]
        payload[prefix + "_time"] = np.float64(cache["time"])
        payload[prefix + "_target_reached"] = np.bool_(cache["target_reached"])
        for name in STATE_NAMES:
            payload[prefix + "_" + name] = cache["state"][name]
            payload[prefix + "_u_" + name] = cache["tangent"][name]
    if not payload:
        return None
    with path.open("wb") as handle:
        np.savez_compressed(handle, **payload)
    size = path.stat().st_size
    if size > RAW_BYTE_LIMIT:
        path.unlink()
        raise ValueError(f"raw cache exceeds {RAW_BYTE_LIMIT} bytes ({size})")
    return {"path": str(path), "sha256": sha256(path), "bytes": int(size),
            "contents": "final reached baseline and tangent for each case that returned a state"}


def baseline_reports(geometry, results):
    reports = []
    with np.load(geometry, allow_pickle=False) as data:
        for result in results:
            cache = result.get("cache") or {}
            comparison = {"times": {}, "within_tolerance": True}
            for label, state in cache.get("checkpoints", {}).items():
                try:
                    reference = archived_frame(data, result["report"]["name"], float(label))
                    gap = compare_checkpoint(state, reference)
                except (KeyError, ValueError, OSError) as error:
                    gap = {"within_tolerance": False, "error": f"{type(error).__name__}: {error}"}
                comparison["times"][label] = gap
                comparison["within_tolerance"] = comparison["within_tolerance"] and bool(
                    gap.get("within_tolerance")
                )
            if not cache.get("checkpoints"):
                comparison["within_tolerance"] = False
                comparison["error"] = "no baseline checkpoint was reached"
            reports.append(comparison)
    return reports


def prepare_application(geometry, clock):
    clock.check("before grid")
    grid = build_grid(512, quadrature=2048, gauge="conformal")
    if (grid.ng, grid.nf, grid.nq, grid.length, grid.gauge) != (511, 512, 2048, PERIOD, "conformal"):
        raise ValueError("this application is the one conformal 511/512/2048 response")
    production_nodes(grid.xi_q)
    clock.check("after grid")
    with np.load(geometry, allow_pickle=False) as data:
        state = _archive_state(data, grid)
    gram_defect = float(np.max(np.abs(column_gram(state.phi0, state.phi1) - np.eye(6))))
    if gram_defect > 1e-8:
        raise ValueError(f"original columns are not orthonormal: gram defect {gram_defect}")
    clock.check("before preparation tangent")
    prepared = preparation_tangent(grid, state, checkpoint=lambda: clock.check("preparation tangent"))
    tangent = initial_tangent(prepared)
    clock.check("before preflight")
    preflight = preflight_jvp(grid, state, tangent)
    clock.check("after preflight")
    couplings = {
        "A": float(grid.fine.A),
        "C_W": float(grid.fine.C_W),
        "C_F": float(grid.fine.C_F),
        "flux": float(grid.fine.flux),
        "kappa": int(grid.fine.kappa),
        "multiplicity": int(grid.fine.multiplicity),
        "occupations": [float(value) for value in grid.fine.occupations],
        "changed": False,
    }
    return state, tangent, prepared.report, preflight, couplings


def application_payload(state, tangent, spec):
    return {
        "name": spec["name"],
        "dt": spec["dt"],
        "target": TARGET_TIME,
        "cpu_seconds": spec["cpu_seconds"],
        "nf": 512,
        "nq": 2048,
        "threads": WORKER_THREADS,
        "diagnostic_times": DIAGNOSTIC_TIMES,
        "sample_events": True,
        "checkpoint_times": (0.05, 0.3),
        "state": pack_state(state),
        "tangent": pack_state(tangent),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    args = parser.parse_args()
    output = require_fresh(args.output, "summary")
    raw_path = require_fresh(args.raw, "raw cache")
    started = time.process_time()
    clock = CpuClock(SETUP_CPU, started)
    commit, bindings, _episode = authenticate(args.geometry)
    clock.check("after authentication")
    try:
        with threadpool_limits(limits=SETUP_THREADS):
            state, tangent, preparation, preflight, couplings = prepare_application(
                args.geometry, clock,
            )
    except CpuAllocationExceeded as error:
        output.write_text(dump_summary({
            "schema": SCHEMA,
            "status": PARTIAL,
            "reason": str(error),
            "source_commit": commit,
            "cases": [],
            "temporal_envelope": None,
            "photon_result": False,
        }), encoding="utf-8", newline="\n")
        return 1
    setup_cpu = clock.used
    background = preparation["background"]
    observer = preparation["observer"]
    result = {
        "schema": SCHEMA,
        "status": PARTIAL,
        "source_commit": commit,
        "source_sha256": bindings,
        "input": {
            "payload_sha256": EXPECTED_PAYLOAD,
            "keys": PREFIX + "<field>",
            "episode_sha256": sha256(EPISODE),
        },
        "initial_constraint": {
            "full_hamilton_max_abs": background["full_hamilton_max_abs"],
            "reference_full_hamilton_max_abs": 9.978453729075909e-06,
            "original_newton_threshold_met": background["original_newton_threshold_met"],
            "initial_defect_removed": background["initial_defect_removed"],
            "zero_root_assumed": False,
            "preserved": True,
        },
        "initial_fractional_clock": {
            "rebuilt": observer["delta_N_over_N_per_epsilon"],
            "reference": REFERENCE_FRACTION,
            "difference": observer["delta_N_over_N_per_epsilon"] - REFERENCE_FRACTION,
            "fit_target": False,
        },
        "preflight": preflight,
        "limits": {
            "setup_cpu_seconds": SETUP_CPU,
            "coarse_cpu_seconds": COARSE_CPU,
            "fine_cpu_seconds": FINE_CPU,
            "total_cpu_seconds": TOTAL_CPU,
            "summary_bytes": BYTE_LIMIT,
            "raw_bytes": RAW_BYTE_LIMIT,
            "workers": WORKERS,
            "setup_threads": SETUP_THREADS,
            "blas_threads_per_worker": WORKER_THREADS,
        },
        "setup_cpu_seconds": setup_cpu,
        "couplings": couplings,
        "counts": {
            "dense_factorizations": 1,
            "root_solves": 0,
            "preflight_compose_calls": 2,
            "parameter_scans": 0,
            "retarded_kernels": 0,
        },
        "scope": (
            "Two original time steps of the same conformal Hamiltonian, nf512 and nq2048. "
            "The paired-phase tangent is differentiated along the evolving state. "
            "Photon rows are fixed-event normal-observer derivatives, not an H0 or sky fit. "
            "The temporal difference is an empirical envelope, not a continuum, spatial, "
            "or physical error bound. No Lambda, formation apparatus, finite-epsilon "
            "observation, retarded kernel, or constraint projection is claimed."
        ),
    }
    if not preflight["accepted"]:
        result["reason"] = "preflight_rejected"
        result["cases"] = []
        result["temporal_envelope"] = None
        result["photon_result"] = False
        result["runtime_cpu_seconds"] = time.process_time() - started
        output.write_text(dump_summary(result), encoding="utf-8", newline="\n")
        print(f"Wrote {output}; preflight rejected; no evolution was launched.")
        return 1
    numerical = setup_cpu
    if numerical > TOTAL_CPU:
        raise RuntimeError("setup already exceeds the total CPU envelope")
    payloads = [application_payload(state, tangent, spec) for spec in CASES]
    collected = collect_cases(payloads, WORKERS)
    by_name = {item["report"]["name"]: item for item in collected}
    ordered = [by_name[spec["name"]] for spec in CASES if spec["name"] in by_name]
    for item in collected:
        if item not in ordered:
            ordered.append(item)
    try:
        comparisons = baseline_reports(args.geometry, ordered)
    except (OSError, ValueError, KeyError) as error:
        comparisons = [{"within_tolerance": False, "error": f"{type(error).__name__}: {error}"}
                       for _item in ordered]
    references = {}
    if LIGHTCONE_SUMMARY.is_file():
        try:
            references = lightcone_redshifts(LIGHTCONE_SUMMARY)
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            references = {}
    cases = []
    numerical_cases = 0.0
    for item, comparison in zip(ordered, comparisons):
        report = dict(item["report"])
        report["baseline"] = comparison
        report["lightcone_redshift"] = redshift_metadata(
            report.get("rays", []), report.get("name"), references,
        )
        numerical_cases += float(report.get("cpu_seconds") or 0.0)
        cases.append(report)
    envelope = temporal_envelope(cases)
    result.update({
        "status": summary_status(cases),
        "cases": cases,
        "temporal_envelope": envelope,
        "baseline_within_tolerance": all(item.get("baseline", {}).get("within_tolerance") for item in cases),
        "numerical_cpu_seconds": numerical + numerical_cases,
        "runtime_cpu_seconds": time.process_time() - started,
        "reproduce": (
            "python scripts/run_python.py scripts/measure_nsc_conformal_causal_response.py "
            f"--geometry {args.geometry} --output {output} --raw {raw_path}"
        ),
    })
    if result["numerical_cpu_seconds"] > TOTAL_CPU:
        result["status"] = PARTIAL
        result["reason"] = "total_cpu_envelope_exceeded"
    try:
        raw_record = write_raw(raw_path, [item.get("cache") for item in ordered])
    except ValueError as error:
        raw_record = {"error": str(error)}
        result["status"] = PARTIAL
        result["reason"] = "raw_cache_rejected"
    result["raw"] = raw_record
    output.write_text(dump_summary(result), encoding="utf-8", newline="\n")
    print(
        f"Wrote {output}; status {result['status']}; "
        f"numerical CPU {result['numerical_cpu_seconds']:.3f}s."
    )
    return 0 if result["status"] == COMPLETED else 1


if __name__ == "__main__":
    raise SystemExit(main())
