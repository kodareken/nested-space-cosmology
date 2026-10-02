"""Numerics-only successor of the frozen leading-Einstein diagnostic.

Physically scaled full local reaction bound; no action, RHS, source or state
change. Existing leading runner/pool/checkpoint machinery is reused.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import hashlib
import json
import math
from pathlib import Path
import resource
from types import SimpleNamespace
import time

import numpy as np
from threadpoolctl import threadpool_limits

from . import nsc_discovery_leading_einstein as leading

SCHEMA = "NSC-DISCOVERY-LEADING-STEP-CONTROL-v2"
PARENT_SCHEMA = "NSC-DISCOVERY-LEADING-EINSTEIN-v1"
PARENT_COMMIT = "e6e784b32213ea2368b7fc863c37c72fba6aa4e3"
NUMERICAL_MODE = "full28_physically_scaled_reaction_with_absolute_moving_scale_rate"
ROOT = leading.ROOT
SOURCE = leading.OUTPUT
OUTPUT = leading.LAB / "results/development/nsc-discovery-leading-einstein-v2"
MAX_NEW_CPU = 165.
MAX_TOTAL_CPU = 600.
RAW_RESTRICTION = leading.step_restriction
SQRT_FIELD_COMPONENTS = 24


class Zero:
    def __matmul__(self, values):
        return np.zeros_like(values)


def local_system(system):
    return SimpleNamespace(A=system.A, C_F=system.C_F, flux=system.flux,
        derivative=Zero(), momentum=Zero(), kappa=system.kappa,
        occupations=system.occupations, multiplicity=system.multiplicity,
        dx=system.dx, length_density=system.length_density, shift=system.shift)


def physical_scales(grid, fine, system):
    a, _Z, _mag = leading.coefficients(system)
    kg = leading.episode.geometry_band_wavenumber(grid)
    kf = float(2*np.pi*np.max(abs(grid.modes_f))/grid.length)
    frequency = max(kg, kf, 1.)
    scales = np.column_stack((fine.Q, fine.r, 2*abs(a)*fine.r**2*frequency/fine.Q,
                              2*abs(a)*fine.r*frequency))
    if not np.isfinite(scales).all() or np.min(scales) <= 0:
        raise ValueError("physical local scales must be finite and positive")
    return scales, frequency, kg, kf


def scaled_local_reaction(grid, fine, system):
    """Infinity row norm of S^-1 J0 S for all four geometry +24 real field slots."""
    scales, frequency, kg, kf = physical_scales(grid, fine, system)
    rows = np.zeros((grid.nq, 28))
    sqrt_dx = np.sqrt(grid.dx_q)
    local = local_system(system)
    for component in range(28):
        tangent = leading.State(*(np.zeros_like(getattr(fine, name)) for name in leading.FIELDS))
        if component < 4:
            setattr(tangent, leading.REAL_FIELDS[component], scales[:, component].copy())
        else:
            index = component-4
            field = tangent.phi0 if index < 12 else tangent.phi1
            index %= 12
            field[:, index % 6] = sqrt_dx*(1 if index < 6 else 1j)
        image = leading.fine_jvp(local, fine, tangent)
        geometric = np.column_stack([getattr(image, name) for name in leading.REAL_FIELDS])/scales
        fields = np.column_stack((image.phi0.real, image.phi0.imag, image.phi1.real, image.phi1.imag))/sqrt_dx
        rows += abs(np.column_stack((geometric, fields)))
    return rows, scales, frequency, kg, kf


def step_restriction(pair, state, step_cap):
    if not math.isfinite(step_cap) or step_cap <= 0:
        raise ValueError("step cap must be finite and positive")
    rate, bundle = leading.rates(pair, state, return_bundle=True)
    fine, system, grid = bundle["fine_state"], bundle["fine_system"], pair.grid
    rows, scales, frequency, kg, kf = scaled_local_reaction(grid, fine, system)
    qt = leading.galerkin.prolong_geometry(grid, pair.geometry_map@rate.Q)
    rt = leading.galerkin.prolong_geometry(grid, pair.geometry_map@rate.r)
    g = np.column_stack((qt/fine.Q, rt/fine.r, 2*rt/fine.r-qt/fine.Q, rt/fine.r))
    frozen_norm = float(np.max(rows))
    moving_norm = float(np.max(abs(g)))
    bound = frozen_norm+moving_norm
    gradient = 2*max(float(np.max(abs(system.derivative@fine.Q/fine.Q))),
                     float(np.max(abs(system.derivative@fine.r/fine.r))))
    kquad = leading.episode.quadrature_principal_wavenumber(grid)
    omega = max(kquad, kf)+bound+gradient
    dt = min(float(step_cap), leading.episode.RK4_HALF_STABILITY/max(omega, 1.))
    return dt, {"dt": dt, "step_cap": float(step_cap), "combined_omega": omega,
        "numerical_mode": NUMERICAL_MODE, "local_real_variables": 28,
        "frozen_scaled_reaction_norm": frozen_norm, "absolute_moving_scale_rate": moving_norm,
        "reaction_majorant": bound, "gradient_majorant": gradient,
        "actual_geometry_band_k": kg, "field_k": kf, "quadrature_k": kquad, "scale_frequency": frequency,
        "scales": "Q,r,2|a|r^2 omega*/Q,2|a|r omega*; psi=Phi/sqrt(dx_q),24 real field slots",
        "moving_scale_sign_cancellation_used": False, "physical_growth_removed": False,
        "principal_coordinate_speed": 1., "field_cfl_alone": False, "stability_certificate": False,
        "scope": "full local weighted norm plus moving-scale triangle bound, principal and gradient terms; not a global stability certificate"}


def worker_initializer():
    leading.SCHEMA = SCHEMA
    leading.step_restriction = step_restriction
    leading.worker_initializer()


def worker_identity():
    """Read-only initializer identity used by the bounded spawn smoke check."""
    return {"leading_schema": leading.SCHEMA, "episode_schema": leading.episode.SCHEMA,
            "scaled_step_installed": leading.step_restriction is step_restriction
                                     and leading.episode.step_restriction is step_restriction,
            "original_RHS_module": leading.rates.__module__, "trajectory_executed": False}


def pool(max_workers):
    return ProcessPoolExecutor(max_workers=max_workers, initializer=worker_initializer)


@contextmanager
def numerical_adapter():
    """Three scoped runtime hooks; frozen action/RHS files and functions are untouched."""
    saved = leading.SCHEMA, leading.step_restriction, leading.pool
    leading.SCHEMA, leading.step_restriction, leading.pool = SCHEMA, step_restriction, pool
    try:
        with threadpool_limits(limits=1), leading.backend.fft_thread_limit(1):
            yield
    finally:
        leading.SCHEMA, leading.step_restriction, leading.pool = saved


def source_hashes(parent):
    result = dict(parent["source_hashes"])
    for directory, name in (
        ("src/recursive_horizons", "nsc_discovery_leading_step_control.py"),
        ("scripts", "derive_nsc_discovery_leading_step_control.py"),
        ("tests", "test_nsc_discovery_leading_step_control.py"),
        ("docs", "nsc-discovery-leading-step-control.md")):
        relative = "lab/"+directory+"/"+name
        result[relative] = leading.episode.file_sha256(ROOT/relative)
    return result


def authenticate_parent(source=SOURCE):
    directory = Path(source).resolve()
    manifest = json.loads((directory/"manifest.json").read_text())
    if manifest.get("schema") != PARENT_SCHEMA or manifest.get("producing_commit") != PARENT_COMMIT:
        raise ValueError("successor requires the frozen e6e784 leading-v1 parent")
    leading.verify_sources(manifest)
    ledger = json.loads((directory/"cpu-ledger.json").read_text())
    if ledger.get("reserved", 0) != 0:
        raise ValueError("parent still has active CPU reservations")
    return manifest, float(ledger["spent"])


def copy_case(source, destination, case_id, *, producer_commit=None):
    """Copy exact latest NPZ bytes and observation prefix; only JSON numerics metadata changes."""
    source, destination = Path(source).resolve(), Path(destination).resolve()
    record, arrays = leading.episode.load_checkpoint(source, case_id)
    ordinal = int(record["ordinal"])
    npz = source/f"{case_id}-{ordinal:06d}.npz"
    original_json = source/f"{case_id}-{ordinal:06d}.json"
    payload = npz.read_bytes()
    old_obs = source/f"{case_id}-observations.jsonl"
    prefix = old_obs.read_bytes() if old_obs.exists() else b""
    prefix_hash = hashlib.sha256(prefix).hexdigest()
    completed = record.get("status") in ("station_reached", "chart_exit")
    updated = dict(record, schema=SCHEMA, ordinal=0, numerical_mode=NUMERICAL_MODE,
                   producing_commit=producer_commit, source_producing_commit=PARENT_COMMIT,
                   predecessor_status=record["status"], predecessor_stop=record.get("stop"),
                   status=record["status"] if completed else "RESUME_PREPARED",
                   stop=record.get("stop") if completed else None,
                   snapshot_kind="exact_numerical_handoff", evolved=False,
                   already_complete=completed, no_resolve_or_reprepare=True,
                   old_numerics_domain={"mode": "e6e784_raw_local_reaction_norm",
                       "through_checkpoint_time": record["coordinate_time"],
                       "last_observation_time": record.get("channel_time"),
                       "prefix_bytes": len(prefix), "prefix_sha256": prefix_hash,
                       "prefix_not_relabelled_as_new_numerics": True},
                   numerical_handoff={"directory": str(source.relative_to(ROOT)),
                       "case_id": case_id, "ordinal": ordinal,
                       "json_sha256": leading.episode.file_sha256(original_json),
                       "npz_sha256": leading.episode.file_sha256(npz),
                       "all_array_hashes": dict(record["array_sha256"])})
    encoded = (json.dumps(updated, indent=2, allow_nan=False)+"\n").encode()
    if len(payload)+len(encoded) > leading.MAX_CHUNK:
        raise ValueError("successor handoff exceeds the existing64-MiB chunk bound")
    leading.episode._commit_bytes(destination, case_id, 0, ".npz", payload, leading.MAX_CHUNK)
    leading.episode._commit_bytes(destination, case_id, 0, ".json", encoded, leading.MAX_CHUNK)
    with (destination/f"{case_id}-observations.jsonl").open("xb") as stream:
        stream.write(prefix)
    reloaded, _ = leading.episode.load_checkpoint(destination, case_id, 0)
    if reloaded["array_sha256"] != record["array_sha256"] or reloaded["arrays_sha256"] != record["arrays_sha256"]:
        raise ValueError("copied handoff lost exact source/state bytes")
    return {"case_id": case_id, "ordinal": 0, "npz": f"{case_id}-000000.npz",
            "json": f"{case_id}-000000.json", "coordinate_time": record["coordinate_time"],
            "steps": record["steps"], "already_complete": completed,
            "prefix_bytes": len(prefix), "prefix_sha256": prefix_hash}


def cpu_usage():
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+child.ru_utime+child.ru_stime


def prepare(source=SOURCE, output=OUTPUT, *, execute=False, producer_commit=None, cpu_budget=MAX_NEW_CPU):
    if not 0 < cpu_budget <= MAX_NEW_CPU:
        raise ValueError("successor CPU budget must be in (0,165]")
    report = {"schema": SCHEMA, "numerical_mode": NUMERICAL_MODE, "status": "PREFLIGHT",
        "old_action_RHS_unchanged": True, "no_interpolation_resolve_reeigenframe_or_reset": True,
        "fresh_CPU_budget": cpu_budget, "total_parent_successor_limit": MAX_TOTAL_CPU,
        "one_existing_pool": True, "workers": 4, "numerical_threads_per_worker": 1,
        "stability_certificate": False, "regeneration_proof": False}
    if not execute:
        return report
    if not producer_commit:
        raise ValueError("successor preparation requires a frozen producer commit")
    started = cpu_usage()
    parent, spent = authenticate_parent(source)
    if spent+cpu_budget > MAX_TOTAL_CPU+1e-9:
        raise ValueError("parent and successor budget would exceed600 CPU seconds")
    destination = leading.episode.assert_campaign_output(output)
    if destination.exists():
        raise FileExistsError("successor already exists")
    sources = source_hashes(parent)
    leading.preparation._git_hashes(PARENT_COMMIT, parent["source_hashes"])
    commit = leading.preparation._git_hashes(producer_commit, sources)
    inputs = {str(path.resolve().relative_to(ROOT)): leading.episode.file_sha256(path)
              for path in Path(source).iterdir() if path.suffix in (".json", ".npz", ".jsonl")}
    destination.mkdir(parents=True, exist_ok=False)
    cases = [copy_case(source, destination, case["case_id"], producer_commit=commit) for case in parent["cases"]]
    used = cpu_usage()-started
    if used >= cpu_budget:
        raise RuntimeError("successor handoff preparation exhausted its fresh CPU budget")
    for name, digest in inputs.items():
        if leading.episode.file_sha256(ROOT/name) != digest:
            raise ValueError("parent changed while successor handoffs were copied")
    manifest = dict(parent)
    manifest.update(report, status="PREPARED_RESUME", stage=0, cases=cases,
        producing_commit=commit, source_hashes=sources, input_hashes=inputs,
        source_directory=str(Path(source).resolve().relative_to(ROOT)),
        source_producing_commit=PARENT_COMMIT, parent_CPU_seconds=spent,
        cpu_budget_seconds=cpu_budget, preparation_cpu_seconds=used,
        child_cpu_seconds=0., evolved=False, results=[],
        numerical_admission_review={"scope": "coordinator-relayed independent constrained review",
            "raw_local_norm_range": [1e5, 1e6], "scaled_local_norm_range": [228., 438.],
            "full_local_geometry_field_eigen_magnitude_below": 123.,
            "true_local_growth_range": [29., 38.], "growth_not_removed": True,
            "whole_finite_weighted_norm_reference": [3435., 7799.],
            "global_norm_computed_per_step": False, "global_certificate": False})
    leading.episode._write_json(destination/"manifest.json", manifest)
    leading.episode._ledger_update(destination, pilot=used, budget=cpu_budget)
    return manifest


def verify_prefixes(output, manifest):
    for case in manifest["cases"]:
        original, _ = leading.episode.load_checkpoint(output, case["case_id"], 0)
        prefix = original["old_numerics_domain"]
        content = (Path(output)/f"{case['case_id']}-observations.jsonl").read_bytes()[:prefix["prefix_bytes"]]
        if hashlib.sha256(content).hexdigest() != prefix["prefix_sha256"]:
            raise ValueError("old observation byte prefix changed")


def run(output=OUTPUT, *, workers=4, cpu_budget=MAX_NEW_CPU, max_steps=None):
    manifest = json.loads((Path(output)/"manifest.json").read_text())
    if manifest.get("schema") != SCHEMA or manifest.get("numerical_mode") != NUMERICAL_MODE:
        raise ValueError("wrong successor numerical schema/mode")
    if not 0 < cpu_budget <= min(MAX_NEW_CPU, manifest["cpu_budget_seconds"]):
        raise ValueError("run cannot enlarge the fresh successor budget")
    if manifest["parent_CPU_seconds"]+cpu_budget > MAX_TOTAL_CPU+1e-9:
        raise ValueError("parent plus successor exceeds600 CPU seconds")
    leading.verify_sources(manifest)
    verify_prefixes(output, manifest)
    started = cpu_usage()
    with numerical_adapter():
        result = leading.run(output, workers=workers, cpu_budget=cpu_budget, max_steps=max_steps)
    ledger_path = Path(output)/"cpu-ledger.json"
    ledger = json.loads(ledger_path.read_text())
    # Include process/reader/worker-init overhead beyond the existing worker
    # timings. Never silently enlarge the budget if accounting reaches it.
    accounted_run = sum(row["child_cpu_seconds"] for row in result.get("results", []))
    extra = max(0., cpu_usage()-started-accounted_run)
    if extra:
        _, ledger = leading.episode._ledger_update(output, pilot=extra, budget=cpu_budget)
    result.update(numerical_mode=NUMERICAL_MODE, successor_CPU_seconds=ledger["spent"],
        parent_and_successor_CPU_seconds=manifest["parent_CPU_seconds"]+ledger["spent"],
        CPU_accounting_includes_reader_and_process_overhead=True,
        status="BUDGET_STOP" if ledger["spent"] >= cpu_budget else result["status"],
        primitive_bottleneck_removed_is_regeneration_proof=False)
    verify_prefixes(output, result)
    leading.episode._write_json(Path(output)/"manifest.json", result)
    return result


def check(output=OUTPUT):
    manifest = json.loads((Path(output)/"manifest.json").read_text())
    if manifest.get("schema") != SCHEMA or manifest.get("numerical_mode") != NUMERICAL_MODE:
        raise ValueError("wrong successor numerical schema/mode")
    leading.verify_sources(manifest)
    verify_prefixes(output, manifest)
    for case in manifest["cases"]:
        handoff, _ = leading.episode.load_checkpoint(output, case["case_id"], 0)
        old = handoff["numerical_handoff"]
        original, _ = leading.episode.load_checkpoint(ROOT/old["directory"], old["case_id"], old["ordinal"])
        if original["array_sha256"] != handoff["array_sha256"] or original["arrays_sha256"] != handoff["arrays_sha256"]:
            raise ValueError("successor handoff differs from exact parent state")
        if case["already_complete"] and leading.episode._next_ordinal(output, case["case_id"]) != 1:
            raise ValueError("an already-completed frozen arm was evolved again")
    with numerical_adapter():
        result = leading.check(output)
    result.update(numerical_mode=NUMERICAL_MODE, parent_handoffs_exact=True, old_prefix_bytes_unchanged=True)
    return result


def paired_step_check(source=SOURCE, case_ids=None):
    """Bounded local full-step/two-half-step comparison; never resumes a case."""
    if case_ids is None:
        case_ids = ("nf256_coherent_dt0.001", "nf256_coherent_dt0.0005")
    parent, _spent = authenticate_parent(source)
    before = {name: leading.episode.file_sha256(ROOT/name) for name in parent["source_hashes"]}
    rows = []
    with threadpool_limits(limits=1), leading.backend.fft_thread_limit(1):
        for case_id in case_ids:
            record, arrays = leading.episode.load_checkpoint(source, case_id)
            pair = leading.backend.make_fft_pair(leading.episode.pair_from_arrays(arrays, record))
            state = leading.state_from_arrays(arrays)
            old_rate = leading.rates(pair, state)
            old_dt, _ = RAW_RESTRICTION(pair, state, record["step_cap"])
            dt, restriction = step_restriction(pair, state, record["step_cap"])
            with numerical_adapter():
                new_rate = leading.rates(pair, state)
                full = leading.rk4_step(pair, state, dt)
                half = leading.rk4_step(pair, leading.rk4_step(pair, state, dt/2), dt/2)
            fine = leading.fine_state(pair, state)
            ff, fh = leading.fine_state(pair, full), leading.fine_state(pair, half)
            scales, _frequency, _kg, _kf = physical_scales(pair.grid, fine, leading.active_system(pair.grid, fine))
            geo = np.column_stack([getattr(ff, name)-getattr(fh, name) for name in leading.REAL_FIELDS])/scales
            phi_diff = np.vstack((full.phi0-half.phi0, full.phi1-half.phi1))
            weights = np.sqrt(pair.weights)
            relative_field = float(np.linalg.norm(phi_diff*weights)) / max(
                float(np.linalg.norm(np.vstack((half.phi0, half.phi1))*weights)), 1e-300)
            identical = all(np.array_equal(getattr(old_rate, name), getattr(new_rate, name)) for name in leading.FIELDS)
            if not identical:
                raise ValueError("numerical adapter changed the leading RHS")
            rows.append({"case_id": case_id, "saved_time": record["coordinate_time"],
                "old_dt": old_dt, "new_dt": dt, "restriction": restriction,
                "one_step_two_half_geometric_scaled_max_gap": float(np.max(abs(geo))),
                "one_step_two_half_weighted_field_relative_gap": relative_field,
                "RHS_bit_identical": identical, "new_trajectory_or_checkpoint_written": False})
    after = {name: leading.episode.file_sha256(ROOT/name) for name in before}
    if before != after:
        raise ValueError("frozen action/source files changed during paired step check")
    return {"checks": rows, "frozen_source_hashes_unchanged": True, "bounded_local_steps_only": True,
            "global_stability_certificate": False}
