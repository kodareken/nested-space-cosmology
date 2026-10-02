#!/usr/bin/env python3
"""Time dense and FFT carriers at 1, 2, and 4 threads. No production writes.

The speed gate uses repeated RK4 steps plus one physical observation at
nf=256, quadrature 4*nf, on one thread. That is the thread count pinned by
the conformal episode drivers. Build time is reported separately.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

_THREAD_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def _preset_threads(argv):
    if "--all-threads" in argv or "--threads" not in argv:
        return
    count = int(argv[argv.index("--threads") + 1])
    for name in _THREAD_VARIABLES:
        os.environ[name] = str(count)


if __name__ == "__main__":
    _preset_threads(sys.argv)


import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


def _rss_bytes():
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value
    return value * 1024


def _time_steps(grid, state, dt, repeats):
    current = state.copy()
    galerkin.rk4_step(grid, current, dt)
    start_wall = time.perf_counter()
    start_cpu = time.process_time()
    for _step in range(repeats):
        current = galerkin.rk4_step(grid, current, dt)
        galerkin.observation(grid, current)
    return {
        "wall_seconds": time.perf_counter() - start_wall,
        "cpu_seconds": time.process_time() - start_cpu,
        "steps": int(repeats),
        "final": current,
    }


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if value is None:
        return None
    return value


def _note(message):
    print(message, file=sys.stderr, flush=True)


def run_measurement(threads, repeats):
    _note(f"threads={threads} building dense nf={backend.PRODUCTION_RESOLUTION_NF}")
    with threadpool_limits(limits=int(threads)), backend.fft_thread_limit(int(threads)):
        build_start = time.perf_counter()
        dense = galerkin.build_grid(
            backend.PRODUCTION_RESOLUTION_NF,
            quadrature=backend.PRODUCTION_QUADRATURE_FACTOR * backend.PRODUCTION_RESOLUTION_NF,
            gauge="conformal",
        )
        dense_build = time.perf_counter() - build_start
        fft_start = time.perf_counter()
        fft_grid = backend.make_fft_grid(dense)
        fft_build = time.perf_counter() - fft_start
        _note(f"threads={threads} dense build {dense_build:.3f}s fft build {fft_build:.3f}s; comparing stages")
        manufactured = backend.manufactured_state(dense, seed=3)
        saved, saved_time = backend.load_saved_conformal_state()
        dt, omega, quadrature_omega = galerkin.stable_timestep(dense, manufactured)
        saved_dt, saved_omega, _saved_quadrature = galerkin.stable_timestep(dense, saved)
        manufactured_oracle = backend.evaluate_stage(dense, manufactured, dt)
        manufactured_fft = backend.evaluate_stage(fft_grid, manufactured, dt)
        manufactured_report = backend.compare_stages(manufactured_oracle, manufactured_fft)
        saved_oracle = backend.evaluate_stage(dense, saved, saved_dt)
        saved_fft = backend.evaluate_stage(fft_grid, saved, saved_dt)
        saved_report = backend.compare_stages(saved_oracle, saved_fft)
        _note(f"threads={threads} timing {repeats} RK4+observe steps")
        dense_timing = _time_steps(dense, manufactured, dt, repeats)
        fft_timing = _time_steps(fft_grid, manufactured, dt, repeats)
        _note(
            f"threads={threads} dense {dense_timing['wall_seconds']:.3f}s "
            f"fft {fft_timing['wall_seconds']:.3f}s"
        )
        speedup = dense_timing["wall_seconds"] / fft_timing["wall_seconds"]
        equivalent = bool(manufactured_report["equivalent"] and saved_report["equivalent"])
        gap_names = (
            "energy", "rate_force_L", "rate_force_Q", "rate_force_beta",
            "force_L_over_dx", "force_Q_over_dx", "force_beta_over_dx",
            "jet_R_h", "jet_Q_dot", "jet_Q_ddot", "step_Q", "step_r", "step_phi0", "step_phi1",
        )
        return {
            "threads": int(threads),
            "resolution_nf": backend.PRODUCTION_RESOLUTION_NF,
            "quadrature": dense.nq,
            "repeats": int(repeats),
            "manufactured_dt": dt,
            "manufactured_omega": omega,
            "quadrature_omega": quadrature_omega,
            "saved_time": saved_time,
            "saved_dt": saved_dt,
            "saved_omega": saved_omega,
            "dense_build_seconds": dense_build,
            "fft_build_seconds": fft_build,
            "dense_wall_seconds": dense_timing["wall_seconds"],
            "fft_wall_seconds": fft_timing["wall_seconds"],
            "dense_cpu_seconds": dense_timing["cpu_seconds"],
            "fft_cpu_seconds": fft_timing["cpu_seconds"],
            "dense_seconds_per_step": dense_timing["wall_seconds"] / repeats,
            "fft_seconds_per_step": fft_timing["wall_seconds"] / repeats,
            "wall_speedup_dense_over_fft": speedup,
            "equivalent": equivalent,
            "peak_rss_bytes": _rss_bytes(),
            "manufactured_one_step_movement": backend.step_movement(manufactured, manufactured_oracle["stepped"]),
            "saved_one_step_movement": backend.step_movement(saved, saved_oracle["stepped"]),
            "repeated_step_backend_gap": backend.step_movement(dense_timing["final"], fft_timing["final"]),
            "manufactured_backend_gap": {name: manufactured_report[name] for name in gap_names},
            "saved_backend_gap": {name: saved_report[name] for name in gap_names},
            "normalization": manufactured_report["oracle_normalization"],
            "gate_min_speedup": backend.PRODUCTION_MIN_SPEEDUP,
            "production_thread": threads == backend.PRODUCTION_THREAD_COUNT,
            "select_fft_at_this_thread_count": bool(
                equivalent and speedup >= backend.PRODUCTION_MIN_SPEEDUP
            ),
        }


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, choices=(1, 2, 4), default=1)
    parser.add_argument("--all-threads", action="store_true")
    parser.add_argument("--repeats", type=int, default=3)
    return parser


def _run_child(threads, repeats):
    env = os.environ.copy()
    for name in _THREAD_VARIABLES:
        env[name] = str(threads)
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--threads",
        str(threads),
        "--repeats",
        str(repeats),
    ]
    completed = subprocess.run(
        command, env=env, cwd=Path.cwd(), text=True, stdout=subprocess.PIPE, stderr=None, check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"benchmark thread count {threads} failed")
    return json.loads(completed.stdout)


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.repeats < 1:
        raise ValueError("repeats must be positive")
    if args.all_threads:
        reports = [_run_child(threads, args.repeats) for threads in (1, 2, 4)]
        production = next(item for item in reports if item["threads"] == backend.PRODUCTION_THREAD_COUNT)
        summary = {
            "threads": reports,
            "production_thread": backend.PRODUCTION_THREAD_COUNT,
            "production_speedup": production["wall_speedup_dense_over_fft"],
            "production_equivalent": production["equivalent"],
            "fft_selected": bool(production["select_fft_at_this_thread_count"]),
            "writes": "none",
        }
        print(json.dumps(_jsonable(summary), indent=2, allow_nan=False))
        return 0
    report = run_measurement(args.threads, args.repeats)
    print(json.dumps(_jsonable(report), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
