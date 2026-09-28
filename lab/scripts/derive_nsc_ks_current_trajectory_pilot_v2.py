#!/usr/bin/env python3
"""Capture one current-history DOP853 cell for a bounded continuous-proof pilot."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import signal
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.path.insert(0, str(ROOT/"scripts"))
import derive_nsc_ks_source_control_v2 as C
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid, usual_axial_support
from recursive_horizons.nsc_ks_streaming_trajectory import stream_ks_trajectory
from recursive_horizons.nsc_ks_family_pool import acquire_dirac_lock, release_dirac_lock
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = ROOT/"results/development/nsc-ks-current-trajectory-pilot-v2.json"
PAYLOAD = ROOT/"results/development/artifacts/nsc-ks-current-trajectory-pilot-v2.npz"
LOCK = ROOT/"results/development/nsc-ks-gate-dirac.lock"
OPTIONS = dict(rtol=5e-14, atol=5e-19, max_step=1/8192,
               tangents="zero", step_control="joint")
CELL = 122
CPU_CAP = 120.0
OWNERS = (
    "scripts/derive_nsc_ks_current_trajectory_pilot_v2.py",
    "scripts/derive_nsc_ks_source_control_v2.py",
    "src/recursive_horizons/nsc_ks_streaming_trajectory.py",
    "src/recursive_horizons/nsc_ks_trajectory.py",
    "src/recursive_horizons/nsc_ks_difference_envelope.py",
    "src/recursive_horizons/nsc_ks_source_envelope.py",
    "src/recursive_horizons/nsc_local_incoming_family.py",
)


def descriptor(path):
    return {"path": str(path.relative_to(ROOT)),
            "sha256": sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size}


def run():
    if OUTPUT.exists() or PAYLOAD.exists():
        raise FileExistsError("pilot capture exists; use --check")
    history = json.loads((ROOT/C.HISTORY).read_text())
    family = LocalIncomingFamily(np.array(history["history"]["coefficients"]))
    identity = profile_identity(family, include_normal_window=True)
    if identity != history["profile_identity"]:
        raise ValueError("history identity mismatch")
    archive = RetainedUpstreamArchive(ROOT)
    source, _negative, initial, first, channel, selections = C.select_source(archive, (14, 1))
    grid, target = computational_z_grid(1024, .4), np.array(family.collocation_nodes(129))
    saved, count = {}, 0
    def emit(segment):
        nonlocal count
        if count == CELL:
            saved.update(rho_start=segment.rho_start, rho_end=segment.rho_end,
                         start=np.array(segment.start), end=np.array(segment.end),
                         dense_corrections=np.array(segment.coefficients))
        count += 1
    def stop(*_):
        raise TimeoutError("current-history trajectory pilot CPU cap")
    acquire_dirac_lock(LOCK, Path(__file__).name, os.getpid())
    start, wall = time.process_time(), time.perf_counter()
    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CPU_CAP)
    try:
        prepared, summary = stream_ks_trajectory(
            source, initial, family, grid, target, first.mass, first.angular, first.rho_up,
            axial_support=usual_axial_support(), on_segment=emit, **OPTIONS)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
        release_dirac_lock(LOCK, os.getpid())
    if not saved:
        raise ValueError("requested proof cell was not reached")
    arrays = {**{k: np.asarray(v) for k, v in saved.items()},
              "source_covariance": source.covariance,
              "source_weights": source.column_weights, "source_energies": source.energies,
              "initial_columns": initial, "computational_z": grid, "target_z": target}
    raw = deterministic_npz_bytes(arrays)
    publish_exclusive_file(ROOT, str(PAYLOAD.relative_to(ROOT)), raw)
    result = {
        "schema": "NSC-KS-CURRENT-TRAJECTORY-PILOT-v2",
        "status": "OPEN: one saved cell, no continuous field certificate",
        "profile_identity": identity, "history_path": C.HISTORY,
        "family": [14, 1], "source_rows": selections,
        "source_preparation_digest": prepared.fixed_preparation_digest,
        "mass": first.mass, "angular": first.angular, "rho_up": first.rho_up,
        "multiplicity": channel["copy_count"]*channel["degeneracy"],
        "source_columns": len(source.energies), "grid_nodes": 1024, "period_length": .4,
        "cell_index": CELL, "accepted_steps": summary.accepted_steps,
        "cell_rho_interval": [saved["rho_end"], saved["rho_start"]],
        "interval": list(family.interval), "options": OPTIONS,
        "payload": descriptor(PAYLOAD),
        "source_hashes": {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
        "input_hashes": {**archive.input_hashes, C.HISTORY: sha256((ROOT/C.HISTORY).read_bytes()).hexdigest()},
        "runtime": {"CPU_seconds": time.process_time()-start,
                    "wall_seconds": time.perf_counter()-wall, "CPU_cap": CPU_CAP},
        "scope": {"saved_cells": 1, "full_source_covered": False,
                  "continuous_field_bound": None, "source_error_bound": None,
                  "physical_EXISTENCE_certificate": False,
                  "step_control": "zero-tangent DOP853 dense-output pilot; not a production cache"},
    }
    publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)),
        (json.dumps(result, indent=2, sort_keys=True)+"\n").encode())
    return result


def check():
    result = json.loads(OUTPUT.read_text())
    for name, expected in {**result["source_hashes"], **result["input_hashes"]}.items():
        if sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError("pilot dependency changed: " + name)
    if descriptor(PAYLOAD) != result["payload"]:
        raise ValueError("pilot payload changed")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--run", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    r = run() if args.run else check()
    print(json.dumps({k: r[k] for k in ("status", "payload", "runtime", "cell_rho_interval")}, indent=2))
