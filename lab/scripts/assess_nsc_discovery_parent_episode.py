#!/usr/bin/env python3
"""Read completed parent episodes. Snapshot NPZ SHA-256 and dtype+shape array hashes.

Historical producers are the frozen commit via Git, not the live parent-episode file.
Saved canonical pi are not signed again. Raw Hamiltonian jets stay diagnostic.
Stdout is the default. --output creates one new readonly file of at most 64 MiB.
"""
from __future__ import annotations

import argparse, hashlib, json, os, sys, time
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import numpy as np
from threadpoolctl import threadpool_limits
from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_extent as extent
from recursive_horizons import nsc_discovery_leading_einstein as leading
from recursive_horizons import nsc_discovery_parent_observer as observer
from recursive_horizons import nsc_discovery_tidal as tidal
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import provenance

SCHEMA, COMMIT, LIMIT = "NSC-DISCOVERY-PARENT-EPISODE-ASSESSMENT-v1", "8c47d0c17fa5962f1d2210f55db3a96c6457a82e", 64*1024*1024
REGIONS = observer.REGION_ORDER
DEFAULTS = ("results/development/nsc-discovery-parent-episode-magnetic-t3-v1",
            "results/development/nsc-discovery-parent-episode-magnetic-nf256-balanced-v1")
_SEEN = {}

def file_hash(path):
    return episode.file_sha256(path)

def summary(values):
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.isfinite(array).all():
        raise ValueError("nonfinite or empty diagnostic")
    scale = float(np.max(np.abs(array)))
    rms = 0. if scale == 0 else scale*float(np.sqrt(np.mean((array/scale)**2)))
    return {"min": float(np.min(array)), "max": float(np.max(array)), "max_abs": scale, "rms": rms}

def authenticate_producers(record):
    commit, hashes = record.get("producing_commit"), record.get("source_hashes") or {}
    if commit != COMMIT or not hashes:
        raise ValueError("snapshot is not bound to the frozen parent producer commit")
    for path, digest in hashes.items():
        key = (commit, path, digest)
        if key not in _SEEN:
            provenance.resolve_pinned_source_bytes(episode.REPO, path, digest, commit=commit)
            _SEEN[key] = True
    return {"producing_commit": commit, "producer_count": len(hashes),
            "historical_producers_authenticated": True, "current_parent_episode_bytes_used": False}

def authenticate_snapshot(directory, json_path):
    json_path = Path(json_path)
    record = json.loads(json_path.read_text())
    npz_path = json_path.with_suffix(".npz")
    for path in (json_path, npz_path):
        if not path.is_file() or path.stat().st_mode & 0o222 or path.stat().st_size > LIMIT:
            raise ValueError("snapshot chunk must be one immutable file within 64 MiB: "+path.name)
    if file_hash(npz_path) != record.get("arrays_sha256"):
        raise ValueError("snapshot NPZ SHA-256 mismatch: "+npz_path.name)
    with np.load(npz_path, allow_pickle=False) as saved:
        arrays = {name: np.array(saved[name], copy=True) for name in saved.files}
    if set(arrays) != set(record.get("array_sha256") or {}):
        raise ValueError("snapshot array inventory mismatch: "+npz_path.name)
    for name, digest in record["array_sha256"].items():
        if episode.array_sha256(arrays[name]) != digest:
            raise ValueError("snapshot array hash mismatch: "+name)
    if int(record.get("column_rank", 0)) != 2 or record.get("momenta_already_signed") is not True:
        raise ValueError("refusing a snapshot that is not stored rank-2 already-signed canonical pi")
    if record.get("projector_source_reset") or record.get("source_reset"):
        raise ValueError("refusing a projector source reset")
    for path, digest in (record.get("input_hashes") or {}).items():
        if file_hash(episode.LAB/path) != digest:
            raise ValueError("historical parent input changed: "+path)
    return record, arrays

def reconstruct(record, arrays):
    """Saved W, occupations, canonical pi, Phi, and clocks. No source solve."""
    nf, weights = int(record["nf"]), np.array(arrays["source_weights"], dtype=float, copy=True)
    if weights.shape != (2,):
        raise ValueError("stored occupation rank is not 2")
    carrier = backend.make_fft_grid(galerkin.build_grid(nf, gauge="conformal"))
    grid = replace(carrier, fine=replace(carrier.fine, occupations=weights.copy()))
    state = leading.State(arrays["Q"], arrays["r"], arrays["pi_Q"], arrays["pi_r"], arrays["phi0"], arrays["phi1"])
    if state.phi0.shape != (nf, 2) or arrays["W"].shape != (grid.ng, grid.ng):
        raise ValueError("stored field or frame does not match the carrier")
    pair = SimpleNamespace(
        grid=grid, geometry_map=np.array(arrays["W"], dtype=float, copy=True), weights=weights,
        reference_columns=np.array(arrays["observer_columns"], copy=True),
        source_phi0=np.array(arrays["source_phi0"], copy=True), source_phi1=np.array(arrays["source_phi1"], copy=True),
        source_columns=np.vstack((arrays["source_phi0"], arrays["source_phi1"])),
        source_metadata=dict(record.get("source_metadata") or {}), geometry_metadata=dict(record.get("geometry_metadata") or {}),
        child_interval=tuple(record["child_interval"]), parent_interval=tuple(record["parent_interval"]),
        clock_locations=tuple(float(value) for value in record["clock_locations"]),
        protected_collar=tuple(record["protected_collar"]),
        parent_annulus=tuple(tuple(float(x) for x in piece) for piece in record["parent_annulus"]),
        column_rank=2, common_k=float(record["common_k"]))
    return pair, state

def measure(record, arrays, producers):
    pair, state = reconstruct(record, arrays)
    mode = record.get("control_mode", "coupled")
    rate, bundle = leading.rates(pair, state, return_bundle=True, control_mode=mode)
    curvature, velocity, acceleration = leading.metric_jets(pair, state, rate, bundle, control_mode=mode)
    row = observer.observe(pair, state, float(record["coordinate_time"]), control_mode=mode,
                           bundle=dict(bundle, leading_rate=rate, metric_jets=(curvature, velocity, acceleration)))
    fine, system, grid = bundle["fine_state"], bundle["fine_system"], pair.grid
    proper, partition, masks = fine.r*fine.Q, observer.region_partition(grid.length), observer.region_masks(grid)
    probability, total = row["channels"]["probability"], float(row["channels"]["probability"]["total"])
    locations = np.asarray(pair.clock_locations, dtype=float)
    centre_index = int(np.argmin(np.abs(locations-0.5*float(grid.length))))
    rates = leading.clock_rates(pair, state)
    gram = state.phi0.conj().T@state.phi0+state.phi1.conj().T@state.phi1
    car = np.linalg.eigvalsh(np.sqrt(pair.weights)[:, None]*gram*np.sqrt(pair.weights)[None, :]).real
    mass = (np.abs(fine.phi0)**2+np.abs(fine.phi1)**2)*pair.weights[None, :]/grid.dx_q
    child = [float(extent.real_interval_integral(grid, mass[:, index], pair.child_interval)) for index in range(2)]
    actual = float(extent.real_interval_integral(grid, mass.sum(axis=1), pair.child_interval))
    dx = lambda values: tidal.metric.spectral_dx(values, grid.length)
    raw_velocity = leading.State(*bundle["unprojected_rates"])
    raw_acceleration = leading.fine_jvp(system, fine, raw_velocity)
    raw = tidal.curvature_from_local_jets(
        fine.Q, fine.r, raw_velocity.Q, dx(fine.Q), raw_acceleration.Q, dx(raw_velocity.Q), dx(dx(fine.Q)),
        raw_velocity.r, dx(fine.r), raw_acceleration.r, dx(raw_velocity.r), dx(dx(fine.r)))
    euler = -8*np.pi*float(system.A)*fine.r**3*fine.Q**2*curvature["tides"]["R4"]
    channels = row["channels"]
    return {"case_id": record["case_id"], "nf": int(record["nf"]), "population_id": record["population_id"],
        "sign_name": record["sign_name"], "time": float(record["coordinate_time"]), "ordinal": int(record["ordinal"]),
        "snapshot_kind": record["snapshot_kind"], "producers": producers,
        "child_fraction": probability["child"]/total,
        "disjoint_fractions": {name: probability[name]/total for name in REGIONS},
        "disjoint_probability": {name: float(probability[name]) for name in REGIONS}, "probability_total": total,
        "proper_lengths": {name: observer._piece_integral(grid, proper, partition[name]) for name in REGIONS},
        "areal_r": {name: summary(fine.r[masks[name]]) for name in REGIONS},
        "N": {name: summary(proper[masks[name]]) for name in REGIONS}, "Nmin": float(np.min(proper)),
        "clocks": {"locations": locations.tolist(),
            "integrated_tau": np.asarray(arrays["normal_clocks"], dtype=float).tolist(),
            "saved_rates": np.asarray(arrays["clock_rates"], dtype=float).tolist(), "recomputed_rates": rates.tolist(),
            "saved_rate_gap_max": float(np.max(np.abs(rates-arrays["clock_rates"]))),
            "centre_index": centre_index, "centre_tau": float(arrays["normal_clocks"][centre_index]),
            "definition": "dτ/dT=rQ; integrated τ is the saved clock"},
        "tides": {name: summary(curvature["tides"][name]) for name in observer.TIDE_CHANNELS},
        "euler_trace_residual": summary(euler), "euler_definition": "-8*pi*A*r^3*Q^2*R4 on projected leading jets",
        "raw_same_point_hamiltonian_R4": summary(raw["tides"]["R4"]), "raw_jets_replace_projected_metric": False,
        "ancestry": {"column_indices": [0, 1], "weights": [float(value) for value in pair.weights],
            "child_probability": child, "sum": float(sum(child)), "actual_child_probability": actual,
            "sum_minus_actual": float(sum(child)-actual), "observer_rebased": False, "particles": False,
            "tag": "original_source_column_order; evolved Phi with original occupations"},
        "energy": leading.energy(pair, state), "gram_distance_from_identity": float(np.max(np.abs(gram-np.eye(2)))),
        "CAR_eigenvalues": car.tolist(), "child_normal_accounting": {
            "normal_energy_rate": channels["normal_energy_rate"]["child"], "pressure": channels["pressure"]["child"],
            "lapse": channels["lapse"]["child"], "work_ledger": {key: record["work_ledger"][key] for key in (
                "coordinate_fieldwork", "pressure_work", "lapse_work", "boundary_child", "boundary_parent", "quadrature")},
            "quadrature_is_a_bound": False},
        "continuum_certified": False, "physical_R4_claimed": False, "renewal_asserted": False}

def series_report(directory, case_id, population):
    path = directory/f"{case_id}-observations.jsonl"
    rows = []
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            result = ((row.get("external_observer") or {}).get("result") or {})
            probability = ((result.get("channels") or {}).get("probability") or {})
            total = probability.get("total")
            if not total:
                continue
            rate = ((result.get("channels") or {}).get("normal_energy_rate") or {}).get("child")
            rows.append({"time": float(row["time"]), "child_fraction": float(probability["child"])/float(total),
                         "child_normal_energy_rate": None if rate is None else float(rate)})
    if not rows:
        return {"available": False}
    return {"available": True, "population_id": population, "samples": len(rows),
        "minimum": min(rows, key=lambda item: item["child_fraction"]),
        "maximum": max(rows, key=lambda item: item["child_fraction"]),
        "nearest": {str(mark): min(rows, key=lambda item: abs(item["time"]-mark)) for mark in (0., 1., 1.25, 2.25, 3.)},
        "column_ancestry_at_series_only_times": False, "full_snapshot_times": [0., 1., 3.],
        "ancestry_absent_because": "full snapshots exist only at T=0,1,3; T=1.25 and T=2.25 are scalar samples"}

def campaign(directory):
    directory = Path(directory).expanduser().resolve()
    manifest = json.loads((directory/"manifest.json").read_text())
    producers, sidecars, measured, series = authenticate_producers(manifest), {str(directory/"manifest.json"): file_hash(directory/"manifest.json")}, [], {}
    for case in manifest["cases"]:
        case_id, chunks = case["case_id"], sorted(directory.glob(case["case_id"]+"-*.json"))
        if not chunks:
            raise ValueError("completed case has no snapshot: "+case_id)
        for json_path in chunks:
            record, arrays = authenticate_snapshot(directory, json_path)
            if record["case_id"] != case_id:
                raise ValueError("snapshot case does not match the manifest")
            authenticate_producers(record)
            measured.append(measure(record, arrays, producers))
        stream = directory/f"{case_id}-observations.jsonl"
        if stream.is_file():
            sidecars[str(stream)] = file_hash(stream)
        series[case_id] = series_report(directory, case_id, measured[-1]["population_id"])
    if any(file_hash(path) != digest for path, digest in sidecars.items()):
        raise ValueError("manifest or scalar stream changed during the read")
    present = {row["case_id"] for row in measured}
    absent = [item["case_id"] for item in manifest.get("confirmation_cases") or [] if item["case_id"] not in present]
    return {"directory": str(directory), "manifest_status": manifest.get("status"), "measured_cases": sorted(present),
            "confirmation_cases_without_snapshots": absent, "stations": measured, "scalar_series": series, "sidecar_sha256": sidecars}

def comparison(stations):
    def pick(row, path):
        for key in path:
            row = row[key]
        return row
    left = next(row for row in stations if row["case_id"] == "nf128_plus_balanced_dt0.001" and abs(row["time"]-3) < 1e-8)
    right = next(row for row in stations if row["case_id"] == "nf256_plus_balanced_dt0.0005" and abs(row["time"]-3) < 1e-8)
    paths = {"child_fraction": ("child_fraction",), "child_proper_length": ("proper_lengths", "child"),
             "centre_tau": ("clocks", "centre_tau"), "R4_global_max_abs": ("tides", "R4", "max_abs"),
             "euler_trace_residual_max_abs": ("euler_trace_residual", "max_abs")}
    drifts = []
    for case_id in sorted({row["case_id"] for row in stations}):
        group = [row for row in stations if row["case_id"] == case_id]
        first = np.asarray(group[0]["CAR_eigenvalues"])
        drifts.append({"case_id": case_id,
            "energy_total_span": max(row["energy"]["total"] for row in group)-min(row["energy"]["total"] for row in group),
            "gram_distance_span": max(row["gram_distance_from_identity"] for row in group)-min(row["gram_distance_from_identity"] for row in group),
            "CAR_span_max": float(np.max(np.abs(np.asarray([row["CAR_eigenvalues"] for row in group])-first))),
            "field_energy_scale": max(abs(row["energy"]["field"]) for row in group)})
    return {"nf_pair": [128, 256], "population_id": "balanced", "sign_name": "plus", "time": 3.0,
            "values_nf128_then_nf256": {name: [pick(left, path), pick(right, path)] for name, path in paths.items()},
            "nf512_measured": False, "global_conservation_spans": drifts,
            "late_curvature_is_a_continuum_value": False, "autonomous_renewal": False}

def report(directories):
    started = time.process_time()
    modules = (backend, episode, extent, leading, observer, tidal, galerkin, provenance)
    sources = {str(Path(module.__file__).resolve()): file_hash(module.__file__) for module in modules}
    sources[str(Path(__file__).resolve())] = file_hash(__file__)
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        campaigns = [campaign(path if Path(path).is_absolute() else episode.LAB/path) for path in directories]
    stations = [row for item in campaigns for row in item["stations"]]
    body = {"schema": SCHEMA, "campaigns": campaigns, "comparison": comparison(stations),
        "measurement_source_hashes": sources, "evolved": False, "reprepared": False, "resigned": False, "healed": False,
        "trajectory_rewritten": False, "input_source_clock_bytes_unchanged": True,
        "unmeasured": {"nf512": "no completed NF512 parent episode in these directories",
            "nf256_other_populations": "only nf256 plus balanced reached T=3",
            "column_ancestry_outside_0_1_3": "scalar samples, including T=1.25 and T=2.25, have no full state",
            "raw_hamiltonian_jets": "diagnostic only; projected leading jets remain the metric",
            "work_ledger_quadrature": "trapezoid coordinate-time indicator, not an error bound",
            "continuum_or_physical_R4": "not claimed", "autonomous_renewal": "not claimed"}}
    encoded = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    body["content_sha256"] = hashlib.sha256(encoded).hexdigest()
    body["cpu_seconds"] = time.process_time()-started
    if any(episode.file_sha256(path) != digest for path, digest in sources.items()):
        raise ValueError("evaluator changed during the read")
    return body

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", nargs="+", type=Path, default=[episode.LAB/path for path in DEFAULTS])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.output and args.output.exists():
            raise FileExistsError("refusing to overwrite assessment output")
        payload = json.dumps(episode.jsonable(report(args.records)), indent=2, allow_nan=False)+"\n"
        if len(payload.encode()) > LIMIT:
            raise ValueError("assessment output exceeds 64 MiB")
        if args.output:
            resolved = args.output.expanduser().resolve()
            if any(resolved.is_relative_to(Path(path).expanduser().resolve()) for path in args.records):
                raise ValueError("assessment output must stay outside the episode directories")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x") as stream:
                stream.write(payload)
            args.output.chmod(0o444)
        print(payload, end="")
        return 0
    except (ValueError, RuntimeError, FileNotFoundError, FileExistsError, KeyError, OSError) as error:
        print(type(error).__name__+": "+str(error), file=sys.stderr)
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
