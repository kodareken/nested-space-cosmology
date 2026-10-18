"""One paired-preparation initial tangent of the authenticated conformal slice.

No root solve, time step, parameter scan, or retarded kernel. The archived
application is a separate step and is not run by the manufactured tests.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

from recursive_horizons.nsc_conformal_preparation_tangent import preparation_tangent
from recursive_horizons.nsc_spherical_coupling import PERIOD, CauchyState
from recursive_horizons.nsc_spherical_galerkin_coupling import build_grid
from measure_nsc_conformal_geometric_average import (
    CURVATURE, EPISODE, EXPECTED_PAYLOAD, sha256,
)


ROOT = Path(__file__).resolve().parents[1]
CPU_LIMIT = 140.0
BYTE_LIMIT = 24_000
THREAD_LIMIT = 8
STATE_NAMES = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
GEOMETRY_NAMES = STATE_NAMES[:6]
PREFIX = "nf512_original_T0_"
SOURCE_PATHS = (
    "scripts/measure_nsc_conformal_preparation_tangent.py",
    "src/recursive_horizons/nsc_conformal_preparation_tangent.py",
    "tests/test_nsc_conformal_preparation_tangent.py",
    "tests/test_nsc_conformal_preparation_fourier_control.py",
    "scripts/measure_nsc_conformal_geometric_average.py",
    "src/recursive_horizons/nsc_spherical_coupling.py",
    "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "src/recursive_horizons/nsc_spherical_feedback_action.py",
    "src/recursive_horizons/nsc_covariant_operator.py",
    "src/recursive_horizons/nsc_conformal_adm_source.py",
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
            raise ValueError(
                "preparation tangent must be frozen before application: " + name
            ) from error
        if raw != committed:
            raise ValueError("preparation tangent must be frozen before application: " + name)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    return commit, hashes


def _archive_state(data, grid):
    missing = [PREFIX + name for name in STATE_NAMES if PREFIX + name not in data.files]
    if missing:
        raise ValueError("original T0 keys are absent: " + ", ".join(missing))
    arrays = {}
    for name in GEOMETRY_NAMES:
        value = np.array(data[PREFIX + name], copy=True)
        if value.shape != (grid.ng,) or value.dtype != np.float64:
            raise ValueError(f"{PREFIX + name} is not float64 of shape {(grid.ng,)}")
        arrays[name] = value
    for name in ("phi0", "phi1"):
        value = np.array(data[PREFIX + name], copy=True)
        if value.shape != (grid.nf, 6) or value.dtype != np.complex128:
            raise ValueError(f"{PREFIX + name} is not complex128 of shape {(grid.nf, 6)}")
        arrays[name] = value
    return CauchyState(**arrays)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--geometry", type=Path,
        default=ROOT / ".build/data/nsc-spherical-conformal-episode-v2.npz",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("use a fresh summary path")
    if not args.output.parent.is_dir():
        raise ValueError("output parent must exist")
    started = time.process_time()

    def budget():
        if time.process_time() - started > CPU_LIMIT:
            raise RuntimeError(
                "preparation tangent exceeded 140 process-CPU seconds; no automatic retry"
            )

    commit, bindings = frozen_sources()
    old = json.loads(CURVATURE.read_text())["source_bindings"]
    episode = json.loads(EPISODE.read_text())
    budget()
    if (
        sha256(args.geometry) != EXPECTED_PAYLOAD
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
    with threadpool_limits(limits=THREAD_LIMIT):
        grid = build_grid(512, quadrature=2048, gauge="conformal")
        if (grid.ng, grid.nf, grid.nq, grid.length, grid.gauge) != (511, 512, 2048, PERIOD, "conformal"):
            raise ValueError("this application is the one conformal 511/512/2048 tangent")
        budget()
        with np.load(args.geometry, allow_pickle=False) as data:
            state = _archive_state(data, grid)
        gram = state.phi0.conj().T @ state.phi0 + state.phi1.conj().T @ state.phi1
        gram_defect = float(np.max(np.abs(gram - np.eye(6))))
        if gram_defect > 1e-8:
            raise ValueError(f"original columns are not orthonormal: gram defect {gram_defect}")
        budget()
        tangent = preparation_tangent(grid, state, checkpoint=budget)
        budget()
        pools = [
            {key: item.get(key) for key in ("internal_api", "num_threads", "prefix", "user_api")}
            for item in threadpool_info()
        ]

    result = {
        "schema": "nsc-conformal-preparation-tangent-summary-v1",
        "status": "initial_constraint_level_response",
        "source_commit": commit,
        "source_sha256": bindings,
        "input": {
            "episode_sha256": sha256(EPISODE),
            "payload_sha256": EXPECTED_PAYLOAD,
            "payload_bytes": args.geometry.stat().st_size,
            "keys": PREFIX + "<field>",
            "frames_loaded": False,
            "authenticated_before_factorization": True,
        },
        "counts": {
            "ng": grid.ng,
            "nf": grid.nf,
            "nq": grid.nq,
            "columns": 6,
            "pairs": 3,
            "dense_factorizations": 1,
            "root_solves": 0,
            "time_steps": 0,
            "parameter_scans": 0,
            "retarded_kernels": 0,
        },
        "threads": {"requested": THREAD_LIMIT, "pools": pools},
        "tangent": tangent.report,
        "limits": {
            "cpu_seconds": CPU_LIMIT,
            "summary_bytes": BYTE_LIMIT,
            "dense_factorizations": 1,
        },
        "runtime_cpu_seconds": time.process_time() - started,
        "scope": (
            "One initial chiral constraint-level tangent on the original T=0 slice. "
            "Derivative per epsilon at fixed initial Q. Not an observed length, "
            "Lambda, formation history, light-cone shift, or retarded kernel. "
            "The projected solve residual does not certify invertibility. "
            "The full quadrature remainder stays in the summary and is not removed."
            " The original nonzero constraint defect is preserved, not solved."
        ),
    }
    raw = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if len(raw.encode()) > BYTE_LIMIT:
        raise ValueError("preparation summary exceeds 24KB")
    args.output.write_text(raw, encoding="utf-8", newline="\n")
    observer = tangent.report["observer"]
    held = tangent.report["directional_check"]["held_out_max_abs"]
    print(
        f"Wrote {args.output}; d(r)/d(eps) at x=1 is {observer['delta_r_per_epsilon']:.12g}; "
        f"d(r)/d(eps)/r = dN/N = {observer['delta_N_over_N_per_epsilon']:.12g}; "
        f"held-out max {held:.12g}; "
        f"{result['runtime_cpu_seconds']:.3f} process-CPU seconds; no evolution."
    )


if __name__ == "__main__":
    main()
