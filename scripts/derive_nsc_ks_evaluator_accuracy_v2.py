#!/usr/bin/env python3
"""Bounded, source-fixed family calibration; indicators are never enclosures."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import derive_nsc_ks_gate_value as V
from recursive_horizons.nsc_ks_energy_propagator import KSEnergyPropagator
from recursive_horizons.nsc_ks_family_pool import acquire_dirac_lock, release_dirac_lock
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeBinding, _sample_axial_profiles, usual_axial_support)
from recursive_horizons.nsc_ks_signed_state import negative_angular_partner
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

SCHEMA = "NSC-KS-EVALUATOR-ACCURACY-v2"
OWNERS = (
    "scripts/derive_nsc_ks_evaluator_accuracy_v2.py",
    "src/recursive_horizons/nsc_ks_difference_envelope.py",
    "src/recursive_horizons/nsc_ks_energy_propagator.py",
    "src/recursive_horizons/nsc_ks_primal_step_control.py",
    "src/recursive_horizons/nsc_ks_source_envelope.py",
    "src/recursive_horizons/nsc_ks_signed_state.py",
    "src/recursive_horizons/nsc_ks_local_constraints.py",
)
HISTORY = "results/development/nsc-ks-gate-history-lm-broyden.json"
FIELDS = ("columns", "axial_columns", "column_tangents", "axial_tangents")


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def paths(name):
    V.require_label(name)
    return (
        ROOT / "results/development" / f"nsc-ks-evaluator-accuracy-v2-{name}.json",
        ROOT / "results/development/artifacts" / f"nsc-ks-evaluator-accuracy-v2-{name}",
    )


def write_once(path, data):
    """Publish a complete checkpoint without replacing scientific bytes."""
    raw = data if isinstance(data, bytes) else (
        json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    created = False
    try:
        with temporary.open("xb") as stream:
            created = True
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        # A hard-link publication is atomic and refuses an existing destination.
        os.link(temporary, path)
    finally:
        if created and temporary.exists():
            temporary.unlink()


def descriptor(path):
    path = Path(path)
    return {"path": str(path.relative_to(ROOT)), "sha256": digest(path),
            "bytes": path.stat().st_size}


def checked_file(item):
    relative = Path(item["path"])
    path = ROOT / relative
    if relative.is_absolute() or ".." in relative.parts or not path.resolve().is_relative_to(ROOT):
        raise ValueError("calibration dependency escapes repository")
    if path.stat().st_size != item["bytes"] or digest(path) != item["sha256"]:
        raise ValueError("calibration dependency changed: " + str(relative))
    return path


def context(history, nodes, overrides):
    ctx = V.context(nodes, overrides, history)
    # Include the implementation even when an older context owner omits it.
    ctx["hashes"] = {**ctx["hashes"], **{name: digest(ROOT / name) for name in OWNERS}}
    return ctx


def run(name, history, nodes, overrides, families, cpu_budget):
    output, directory = paths(name)
    if output.exists():
        raise FileExistsError("calibration already recorded; use --check")
    if not np.isfinite(cpu_budget) or cpu_budget <= 0:
        raise ValueError("CPU budget must be positive")
    ctx = context(history, nodes, overrides)
    keys = tuple(sorted(set(families)))
    if len(keys) != len(families) or any(key not in ctx["archived"] for key in keys):
        raise ValueError("distinct original positive-family keys required")
    binding = {
        "profile_identity": ctx["identity"], "history_path": ctx["history_path"],
        "target": ctx["target"].tolist(), "solver": ctx["solver"],
        "family_keys": [list(key) for key in keys], "source_hashes": ctx["hashes"],
    }
    binding_path = directory / "binding.json"
    if binding_path.exists():
        if json.loads(binding_path.read_text()) != binding:
            raise ValueError("checkpoint belongs to another source/history/numerics")
    else:
        write_once(binding_path, binding)

    lock = acquire_dirac_lock(V.LOCK, Path(__file__).name, os.getpid())
    started, wall = time.process_time(), time.perf_counter()
    rows, complete = [], True

    def exhausted(*_):
        raise TimeoutError("calibration CPU budget exhausted")

    previous = signal.signal(signal.SIGPROF, exhausted)
    try:
        for key in keys:
            stem = f"family-{key[0]:02d}-{key[1]:+d}"
            row_path, payload = directory / (stem + ".json"), directory / (stem + ".npz")
            if row_path.exists():
                row = json.loads(row_path.read_text())
                checked_file(row["payload"])
                if row["binding_sha256"] != digest(binding_path):
                    raise ValueError("family checkpoint binding changed")
                rows.append(row)
                continue
            if payload.exists():
                raise ValueError("orphan calibration payload: " + str(payload))
            remaining = cpu_budget - (time.process_time() - started)
            if remaining <= 0:
                complete = False
                break
            begin = time.process_time()
            signal.setitimer(signal.ITIMER_PROF, remaining)
            try:
                matter, signed, operator = V.evolve_family(ctx, key)
            except TimeoutError:
                complete = False
                break
            finally:
                signal.setitimer(signal.ITIMER_PROF, 0)
            arrays = {field: np.asarray(getattr(operator, field)) for field in V.KEYS}
            arrays["matter_change"] = np.asarray(matter)
            meta = {
                "interval": operator.interval, "mass": operator.mass,
                "angular": operator.angular, "rho_up": operator.rho_up,
                "diagnostics": dict(operator.diagnostics),
                "family_records": signed,
            }
            arrays["metadata_json"] = np.frombuffer(
                json.dumps(meta, sort_keys=True).encode(), dtype=np.uint8)
            write_once(payload, deterministic_npz_bytes(arrays))
            row = {
                "family": list(key), "payload": descriptor(payload),
                "binding_sha256": digest(binding_path),
                "matter_change_maxima": np.max(np.abs(matter), axis=0).tolist(),
                "CPU_seconds": time.process_time() - begin,
                "operator_diagnostics": dict(operator.diagnostics),
            }
            write_once(row_path, row)
            rows.append(row)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
        release_dirac_lock(lock, os.getpid())
    result = {
        "schema": SCHEMA, "status": "OPEN: family calibration indicators only",
        "indicator_not_a_bound": True, "complete": complete,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        **binding, "binding": descriptor(binding_path), "families": rows,
        "runtime": {"CPU_seconds": time.process_time() - started,
                    "wall_seconds": time.perf_counter() - wall, "CPU_budget": cpu_budget},
    }
    if complete:
        write_once(output, result)
    return result


def load(name, *, verify_sources=True):
    output, _directory = paths(name)
    result = json.loads(output.read_text())
    if result["schema"] != SCHEMA or not result["complete"]:
        raise ValueError("complete calibration record required")
    binding_path = checked_file(result["binding"])
    binding = json.loads(binding_path.read_text())
    for field, value in binding.items():
        if result.get(field) != value:
            raise ValueError("calibration header/binding differs: " + field)
    if verify_sources:
        for path, expected in result["source_hashes"].items():
            relative = Path(path)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("unsafe source path")
            if digest(ROOT / relative) != expected:
                raise ValueError("calibration source changed: " + path)
    if [row["family"] for row in result["families"]] != result["family_keys"]:
        raise ValueError("calibration family coverage differs")
    for row in result["families"]:
        if row["binding_sha256"] != result["binding"]["sha256"]:
            raise ValueError("family binding differs")
        path = checked_file(row["payload"])
        with np.load(path, allow_pickle=False) as values:
            matter = np.asarray(values["matter_change"])
            if matter.shape != (len(result["target"]), 2) or not np.isfinite(matter).all():
                raise ValueError("invalid family constraint array")
            if np.max(np.abs(matter), axis=0).tolist() != row["matter_change_maxima"]:
                raise ValueError("family maxima do not reconstruct")
    return result


def restore_operator(ctx, row):
    path = checked_file(row["payload"])
    with np.load(path, allow_pickle=False) as handle:
        arrays = {name: np.array(handle[name]) for name in V.KEYS}
        matter = np.array(handle["matter_change"])
        meta = json.loads(handle["metadata_json"].tobytes())
    metric, solver = ctx["family"].metric(), ctx["solver"]
    w, u = _sample_axial_profiles(metric.directions, ctx["grid"])
    binding = KSEnvelopeBinding(
        ctx["grid"], w, u, tuple(metric.amplitudes),
        tuple(d.inner_radius for d in metric.directions),
        tuple(d.outer_radius for d in metric.directions), usual_axial_support(),
        meta["rho_up"], 1.0, solver["rtol"], solver["atol"], solver["max_step"])
    operator = KSEnergyPropagator(
        meta["interval"], *(arrays[name] for name in V.KEYS),
        meta["mass"], meta["angular"], meta["rho_up"], binding, meta["diagnostics"])
    return operator, matter


def physical_states(ctx, key, operator):
    """Reconstruct the actual retained signed sources, not identity columns."""
    entries = []
    for batch, channel in ctx["archive"].family_entries(key):
        normalized = _channel_record("_", {"_": channel})
        tag = f"{batch.group}_{batch.angular_sign}_E{batch.energy_sign:+d}:{batch.panel_name}"
        entries.append((tag, batch, normalized))
    groups = V.H.group_signed_operator_families(entries)
    if len(groups) != 1:
        raise ValueError("one operator family required")
    states = {}
    for batch_id, batch, _channel in groups[0]["applies"]:
        states[batch_id] = operator.apply(batch.source, batch.initial_columns)
    for pos_id, neg_id, batch, _channel in groups[0]["maps"]:
        states[neg_id] = negative_angular_partner(states[pos_id], batch.source)
    return states


def compare(left_name, right_name):
    left, right = load(left_name), load(right_name)
    for field in ("profile_identity", "history_path", "target", "family_keys", "source_hashes"):
        if left[field] != right[field]:
            raise ValueError("comparison changes physical input or code: " + field)
    contexts = [context(row["history_path"], len(row["target"]), row["solver"])
                for row in (left, right)]
    comparisons = []
    for lrow, rrow in zip(left["families"], right["families"]):
        lop, lmatter = restore_operator(contexts[0], lrow)
        rop, rmatter = restore_operator(contexts[1], rrow)
        key = tuple(lrow["family"])
        lstate, rstate = (physical_states(ctx, key, op)
                          for ctx, op in zip(contexts, (lop, rop)))
        if lstate.keys() != rstate.keys():
            raise ValueError("signed source coverage changed")
        differences = dict.fromkeys(FIELDS, 0.0)
        for batch_id in lstate:
            a, b = lstate[batch_id], rstate[batch_id]
            if a.fixed_preparation_digest != b.fixed_preparation_digest:
                raise ValueError("comparison changes upstream preparation")
            for field in FIELDS:
                delta = np.asarray(getattr(a, field)) - np.asarray(getattr(b, field))
                differences[field] = max(differences[field], float(np.max(np.abs(delta), initial=0)))
        comparisons.append({
            "family": list(key), "physical_field_max_differences": differences,
            "tangent_directions": int(lop.tangent.shape[1]),
            "matter_change_max_difference": np.max(np.abs(lmatter-rmatter), axis=0).tolist(),
        })
    return {
        "schema": "NSC-KS-EVALUATOR-COMPARISON-v2",
        "left": descriptor(paths(left_name)[0]), "right": descriptor(paths(right_name)[0]),
        "left_solver": left["solver"], "right_solver": right["solver"],
        "families": comparisons, "indicator_not_a_bound": True,
        "full_source_covered": False,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--run", action="store_true")
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--compare", nargs=2, metavar=("LEFT", "RIGHT"))
    parser.add_argument("--name")
    parser.add_argument("--history", default=HISTORY)
    parser.add_argument("--nodes", type=int, default=129)
    parser.add_argument("--families", nargs="+", default=["14:1", "32:-1"])
    parser.add_argument("--cpu-budget", type=float, default=1200)
    parser.add_argument("--grid", type=int, default=256)
    parser.add_argument("--length", type=float, default=1.6)
    parser.add_argument("--degree", type=int, default=32)
    parser.add_argument("--integrator", choices=("cf4", "dop853"), default="cf4")
    parser.add_argument("--max-step", type=float, default=0.0005)
    parser.add_argument("--rtol", type=float, default=2e-11)
    parser.add_argument("--atol", type=float, default=2e-13)
    args = parser.parse_args()
    if args.compare:
        result = compare(*args.compare)
    elif not args.name:
        parser.error("--name is required for --run/--check")
    elif args.check:
        result = load(args.name)
    else:
        settings = {"grid_nodes": args.grid, "length": args.length, "degree": args.degree,
                    "integrator": args.integrator, "max_step": args.max_step,
                    "rtol": args.rtol, "atol": args.atol}
        families = [tuple(map(int, item.split(":"))) for item in args.families]
        result = run(args.name, args.history, args.nodes, settings, families, args.cpu_budget)
    if args.compare:
        print(json.dumps(result, indent=2, allow_nan=False))
    else:
        print(json.dumps({key: result[key] for key in (
            "status", "complete", "families", "runtime", "physical_EXISTENCE_certificate")},
            indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
