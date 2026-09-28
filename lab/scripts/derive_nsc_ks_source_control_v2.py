#!/usr/bin/env python3
"""Four unchanged source-energy rows per family: a bounded field control."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import block_diag

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_family_pool import acquire_dirac_lock, release_dirac_lock
from recursive_horizons.nsc_ks_local_constraints import _family_matter
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_signed_state import negative_angular_partner, s3_conjugate
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid, usual_axial_support
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_evolved_incoming_constraints import coefficients_from_records

HISTORY = "results/development/nsc-ks-gate-history-lm-broyden.json"
FIELDS = ("columns", "axial_columns", "column_tangents", "axial_tangents")
OWNERS = (
    "scripts/derive_nsc_ks_source_control_v2.py",
    "src/recursive_horizons/nsc_ks_difference_envelope.py",
    "src/recursive_horizons/nsc_ks_primal_step_control.py",
    "src/recursive_horizons/nsc_ks_source_envelope.py",
    "src/recursive_horizons/nsc_ks_local_constraints.py",
    "src/recursive_horizons/nsc_ks_matter_difference.py",
    "src/recursive_horizons/nsc_ks_signed_state.py",
    "src/recursive_horizons/nsc_ks_retained_upstream_archive.py",
    "src/recursive_horizons/nsc_local_incoming_family.py",
)
LOCK = ROOT / "results/development/nsc-ks-gate-dirac.lock"


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def output_paths(name):
    if not re.fullmatch("[a-z0-9-]+", name):
        raise ValueError("simple lowercase control name required")
    return (ROOT / "results/development" / f"nsc-ks-source-control-v2-{name}.json",
            ROOT / "results/development/artifacts" / f"nsc-ks-source-control-v2-{name}.npz")


def select_source(archive, key, count=4):
    """Retain complete original coherent 3-column blocks without refitting."""
    if type(count) is not int or count < 1:
        raise ValueError("positive number of original rows required")
    entries = archive.family_entries(key)
    negative = {batch.panel_name: batch for batch, _ in entries if batch.energy_sign < 0}
    pool = []
    for batch, channel in entries:
        if batch.energy_sign < 0:
            continue
        energies = batch.source.energies
        if len(energies) % 3:
            raise ValueError("three source columns per original energy required")
        for offset in range(0, len(energies), 3):
            if not np.all(energies[offset:offset+3] == energies[offset]):
                raise ValueError("source-energy coherence block changed")
            pool.append((float(energies[offset]), batch, channel, offset))
    pool.sort(key=lambda row: row[0])
    if len(pool) < count:
        raise ValueError("not enough original energy rows")
    chosen = [pool[i] for i in np.linspace(0, len(pool)-1, count, dtype=int)]
    covariance, neg_covariance, weights, energies, initial, selections = [], [], [], [], [], []
    for energy, batch, channel, offset in chosen:
        sl = slice(offset, offset+3)
        neg = negative[batch.panel_name + "/negative-partner"]
        if (not np.array_equal(neg.source.energies[sl], -batch.source.energies[sl])
                or not np.array_equal(neg.source.column_weights[sl], batch.source.column_weights[sl])
                or not np.array_equal(neg.initial_columns[:, sl],
                                      s3_conjugate(batch.initial_columns[:, sl]))):
            raise ValueError("original signed preparation relation changed")
        first = chosen[0][1]
        if (batch.mass, batch.angular, batch.rho_up) != (first.mass, first.angular, first.rho_up):
            raise ValueError("source rows belong to different channels")
        for source in (batch.source, neg.source):
            off = np.array(source.covariance[sl], copy=True)
            off[:, sl] = 0
            if np.count_nonzero(off):
                raise ValueError("row restriction would discard cross-energy coherence")
        covariance.append(batch.source.covariance[sl, sl])
        neg_covariance.append(neg.source.covariance[sl, sl])
        weights.extend(batch.source.column_weights[sl])
        energies.extend(batch.source.energies[sl])
        initial.append(batch.initial_columns[:, sl])
        selections.append({"panel": batch.panel_name, "row_offset": offset, "energy": energy})
    positive = FixedSourcePreparation(block_diag(*covariance), np.array(weights), np.array(energies))
    neg_source = FixedSourcePreparation(block_diag(*neg_covariance), np.array(weights), -np.array(energies))
    first, channel = chosen[0][1:3]
    return positive, neg_source, np.concatenate(initial, axis=1), first, channel, selections


def run(args):
    output, payload = output_paths(args.name)
    if output.exists() or payload.exists():
        raise FileExistsError("control output exists; use --check or a new name")
    history_path = ROOT / HISTORY
    history = json.loads(history_path.read_text())
    family = LocalIncomingFamily(np.array(history["history"]["coefficients"]))
    if profile_identity(family, include_normal_window=True) != history["profile_identity"]:
        raise ValueError("history identity differs")
    grid = computational_z_grid(args.grid, args.length)
    target = np.array(family.collocation_nodes(args.nodes))
    archive = RetainedUpstreamArchive(ROOT)
    cpath = ROOT / "results/development/nsc-incoming-surface-coefficients.json"
    bpath = ROOT / "results/development/nsc-incoming-surface-regular-branch.json"
    coefficients = coefficients_from_records(json.loads(cpath.read_text()), json.loads(bpath.read_text()))
    settings = {
        "grid_nodes": args.grid, "length": args.length, "target_nodes": target.tolist(),
        "integrator": args.integrator, "max_step": args.max_step,
        "rtol": args.rtol, "atol": args.atol, "step_control": "primal",
        "tangents": args.tangents,
    }
    arrays, rows = {}, []
    signatures = {path: digest(ROOT / path) for path in OWNERS}
    code_commit = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    lock = acquire_dirac_lock(LOCK, Path(__file__).name, os.getpid())
    cpu, wall = time.process_time(), time.perf_counter()

    def exhausted(*_):
        raise TimeoutError("original-source control CPU budget exhausted")

    previous = signal.signal(signal.SIGPROF, exhausted)
    signal.setitimer(signal.ITIMER_PROF, args.cpu_budget)
    try:
        for key in ((14, 1), (32, -1)):
            positive, negative, initial, batch, channel, selection = select_source(archive, key)
            before = time.process_time()
            state = evolve_ks_difference_envelope(
                positive, initial, family, grid, target, batch.mass, batch.angular, batch.rho_up,
                axial_support=usual_axial_support(), integrator=args.integrator,
                max_step=args.max_step, rtol=args.rtol, atol=args.atol, tangents=args.tangents,
                step_control="primal",
                tangent_rtol=args.rtol if args.tangents == "all" else None,
                tangent_atol=args.atol if args.tangents == "all" else None)
            partner = negative_angular_partner(state, negative)
            for suffix, prepared, signed in (
                    ("positive", state, channel),
                    ("negative", partner, {**channel, "angular_sign": -channel["angular_sign"]})):
                prefix = f"{key[0]}_{key[1]}/{suffix}/"
                for field in FIELDS:
                    arrays[prefix + field] = np.asarray(getattr(prepared, field))
                normalized = _channel_record("_", {"_": signed})
                matter = _family_matter(
                    prepared, normalized, coefficients, len(prepared.column_tangents), target)
                arrays[prefix + "matter_change"] = np.asarray(matter["correction"])
                arrays[prefix + "matter_tangent"] = np.asarray(matter["tangent"])
            rows.append({
                "family": list(key), "selections": selection,
                "preparation_digest": state.fixed_preparation_digest,
                "negative_preparation_digest": partner.fixed_preparation_digest,
                "CPU_seconds": time.process_time() - before,
                "diagnostics": dict(state.diagnostics),
            })
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
        release_dirac_lock(lock, os.getpid())
    if signatures != {path: digest(ROOT / path) for path in OWNERS}:
        raise ValueError("propagation code changed during the control")
    record = {
        "schema": "NSC-KS-ORIGINAL-SOURCE-CONTROL-v2",
        "status": "OPEN: finite-source field control; no continuum or full-source enclosure",
        "profile_identity": history["profile_identity"], "history": HISTORY,
        "solver": settings, "families": rows,
        "input_hashes": {**archive.input_hashes, HISTORY: digest(history_path),
                        str(cpath.relative_to(ROOT)): digest(cpath), str(bpath.relative_to(ROOT)): digest(bpath)},
        "source_hashes": signatures, "source_interpolated": False,
        "code_commit": code_commit,
        "source_runs": 0, "full_source_covered": False,
        "physical_EXISTENCE_certificate": False, "physical_NONEXISTENCE_certificate": False,
        "runtime": {"CPU_seconds": time.process_time()-cpu, "wall_seconds": time.perf_counter()-wall},
    }
    raw = deterministic_npz_bytes(arrays)
    record["payload"] = {"path": str(payload.relative_to(ROOT)), "sha256": sha256(raw).hexdigest(), "bytes": len(raw)}
    payload.parent.mkdir(parents=True, exist_ok=True)
    publish_exclusive_file(ROOT, str(payload.relative_to(ROOT)), raw)
    publish_exclusive_file(ROOT, str(output.relative_to(ROOT)),
                           (json.dumps(record, sort_keys=True, indent=2, allow_nan=False)+"\n").encode())
    return record


def load(name):
    output, payload = output_paths(name)
    record = json.loads(output.read_text())
    if record["schema"] != "NSC-KS-ORIGINAL-SOURCE-CONTROL-v2":
        raise ValueError("wrong source-control schema")
    for path, expected in {**record["source_hashes"], **record["input_hashes"]}.items():
        if digest(ROOT / path) != expected:
            raise ValueError("source-control dependency changed: " + path)
    if digest(payload) != record["payload"]["sha256"] or payload.stat().st_size != record["payload"]["bytes"]:
        raise ValueError("source-control payload changed")
    with np.load(payload, allow_pickle=False) as handle:
        arrays = {key: np.array(handle[key]) for key in handle.files}
    return record, arrays


def compare(left_name, right_name):
    left, x = load(left_name)
    right, y = load(right_name)
    if (left["profile_identity"] != right["profile_identity"]
            or left["input_hashes"] != right["input_hashes"]
            or left["source_hashes"] != right["source_hashes"]
            or left["solver"]["target_nodes"] != right["solver"]["target_nodes"]):
        raise ValueError("comparison changes source/history/code or physical sample points")
    for a, b in zip(left["families"], right["families"]):
        for key in ("family", "selections", "preparation_digest", "negative_preparation_digest"):
            if a[key] != b[key]:
                raise ValueError("source-control family changed")
    if x.keys() != y.keys():
        raise ValueError("field coverage differs")
    differences = {}
    for key in x:
        if x[key].shape != y[key].shape:
            raise ValueError("field shapes differ")
        delta = x[key]-y[key]
        differences[key] = (np.max(np.abs(delta), axis=0).tolist()
                            if key.endswith("matter_change") else float(np.max(np.abs(delta), initial=0)))
    return {
        "left": left_name, "right": right_name, "differences": differences,
        "indicator_not_a_bound": True, "full_source_covered": False,
        "physical_EXISTENCE_certificate": False,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--compare", nargs=2)
    p.add_argument("--name")
    p.add_argument("--grid", type=int, default=256)
    p.add_argument("--length", type=float, default=1.6)
    p.add_argument("--nodes", type=int, default=129)
    p.add_argument("--integrator", choices=("cf4", "dop853"), default="cf4")
    p.add_argument("--max-step", type=float, default=0.0005)
    p.add_argument("--rtol", type=float, default=2e-11)
    p.add_argument("--atol", type=float, default=2e-13)
    p.add_argument("--tangents", choices=("zero", "all"), default="zero")
    p.add_argument("--cpu-budget", type=float, default=600)
    args = p.parse_args()
    if args.compare:
        result = compare(*args.compare)
    elif not args.name:
        p.error("--name required")
    elif args.check:
        result = load(args.name)[0]
    else:
        if not np.isfinite(args.cpu_budget) or args.cpu_budget <= 0:
            p.error("positive CPU budget required")
        result = run(args)
    print(json.dumps(result if args.compare else {
        "status": result["status"], "runtime": result["runtime"],
        "families": result["families"]}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
