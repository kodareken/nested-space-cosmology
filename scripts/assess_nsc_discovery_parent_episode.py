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

SCHEMA = "NSC-DISCOVERY-PARENT-EPISODE-ASSESSMENT-v2"
V1_SCHEMA = "NSC-DISCOVERY-PARENT-EPISODE-ASSESSMENT-v1"
V1_DIGEST = "79dbc6a7bca6010c160a92784dd5d46f57e56ba29b79e6e9d77eb5b79f992bdf"
V1_RECORD = "results/development/nsc-discovery-parent-episode-assessment-v1.json"
COMMIT, LIMIT = "8c47d0c17fa5962f1d2210f55db3a96c6457a82e", 64*1024*1024
RESPONSIVE_SCHEMA = "NSC-DISCOVERY-PARENT-RESPONSIVE-CONTROL-ASSESSMENT-v1"
REGIONS = observer.REGION_ORDER
DEFAULTS = ("results/development/nsc-discovery-parent-episode-magnetic-t3-v1",
            "results/development/nsc-discovery-parent-episode-magnetic-nf256-balanced-v1")
RESPONSIVE = ("results/development/nsc-discovery-parent-strong-t1-v1",
              "results/development/nsc-discovery-parent-strong-t3-v1",
              "results/development/nsc-discovery-parent-frozen-parent-heavy-v1",
              "results/development/nsc-discovery-parent-empty-source-v1")
MAGNETIC = "results/development/nsc-discovery-parent-episode-magnetic-t3-v1"
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

def authenticate_producers(record, pinned=COMMIT):
    commit, hashes = record.get("producing_commit"), record.get("source_hashes") or {}
    if pinned is not None and commit != pinned or not hashes:
        raise ValueError("snapshot is not bound to the frozen parent producer commit")
    if not isinstance(commit, str) or len(commit) != 40:
        raise ValueError("full producer commit required")
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

def _share(part, total):
    return 0.0 if total == 0.0 else part/total

def measure(record, arrays, producers, detail=False):
    pair, state = reconstruct(record, arrays)
    mode = record.get("control_mode", "coupled")
    # Source-free columns already carry the empty source. Their geometry rate stays the coupled rate.
    rate_mode = "coupled" if mode in ("source_free", "source-free") else mode
    rate, bundle = leading.rates(pair, state, return_bundle=True, control_mode=rate_mode)
    curvature, velocity, acceleration = leading.metric_jets(pair, state, rate, bundle, control_mode=rate_mode)
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
    result = {"case_id": record["case_id"], "nf": int(record["nf"]), "population_id": record["population_id"],
        "sign_name": record["sign_name"], "time": float(record["coordinate_time"]), "ordinal": int(record["ordinal"]),
        "snapshot_kind": record["snapshot_kind"], "producers": producers,
        "child_fraction": _share(probability["child"], total),
        "disjoint_fractions": {name: _share(probability[name], total) for name in REGIONS},
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
    if detail:
        header = record.get("observation") if isinstance(record.get("observation"), dict) else {}
        header_time = header.get("time")
        instant = float(record["coordinate_time"])
        stale = header_time is not None and abs(float(header_time)-instant) > 1e-8
        peak = row["child_peak"]
        result["child_peak"] = {"value": peak["value"], "x": peak["x"], "signed_distance": peak["signed_distance"],
            "is_global_peak": peak["is_global_peak"], "global_peak_region": row["flags"]["global_peak_region"]}
        result["common_k"] = float(record["common_k"])
        result["parent_k"] = None if record.get("parent_k") is None else float(record["parent_k"])
        result["control_mode"] = mode
        result["checkpoint_header_audit"] = {
            "present": bool(header), "time": None if header_time is None else float(header_time), "stale": stale,
            "used_as_measurement": False,
            "stored_child_proper_length": header.get("child_proper_length"),
            "stored_child_probability": header.get("child_probability"),
            "recomputed_child_proper_length": result["proper_lengths"]["child"],
            "length_header_minus_recomputed": None if header.get("child_proper_length") is None else
                float(header["child_proper_length"])-result["proper_lengths"]["child"]}
    return result

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

def _canon(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)

def _measurement_sha256(campaigns, compared):
    return hashlib.sha256(_canon({"campaigns": campaigns, "comparison": compared}).encode()).hexdigest()

def _sealed_v1():
    path = episode.LAB/V1_RECORD
    sealed = json.loads(path.read_text())
    predecessor = {key: value for key, value in sealed.items() if key not in ("content_sha256", "cpu_seconds")}
    if sealed.get("schema") != V1_SCHEMA or hashlib.sha256(_canon(predecessor).encode()).hexdigest() != sealed.get("content_sha256"):
        raise ValueError("sealed v1 assessment record failed its own content digest")
    if sealed["content_sha256"] != V1_DIGEST:
        raise ValueError("sealed v1 assessment digest is not the accepted predecessor")
    return path, sealed

def report(directories):
    started = time.process_time()
    modules = (backend, episode, extent, leading, observer, tidal, galerkin, provenance)
    sources = {str(Path(module.__file__).resolve()): file_hash(module.__file__) for module in modules}
    script = str(Path(__file__).resolve())
    sources[script] = file_hash(script)
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        campaigns = [campaign(path if Path(path).is_absolute() else episode.LAB/path) for path in directories]
    stations = [row for item in campaigns for row in item["stations"]]
    compared = comparison(stations)
    measurement = _measurement_sha256(campaigns, compared)
    sealed_path, sealed = _sealed_v1()
    predecessor_measurement = _measurement_sha256(sealed["campaigns"], sealed["comparison"])
    if measurement != predecessor_measurement:
        raise ValueError("default measurement diverged from the sealed v1 campaigns and comparison")
    body = {"schema": SCHEMA, "predecessor": {
            "schema": V1_SCHEMA, "content_sha256": V1_DIGEST,
            "sealed_record": V1_RECORD, "sealed_record_sha256": file_hash(sealed_path),
            "measurement_sha256": predecessor_measurement},
        "campaigns": campaigns, "comparison": compared, "measurement_sha256": measurement,
        "measurement_unchanged_from_predecessor": True,
        "measurement_source_hashes": sources, "evolved": False, "reprepared": False, "resigned": False, "healed": False,
        "trajectory_rewritten": False, "input_source_clock_bytes_unchanged": True,
        "unmeasured": {"nf512": "no completed NF512 parent episode in these directories",
            "nf256_other_populations": "only nf256 plus balanced reached T=3",
            "column_ancestry_outside_0_1_3": "scalar samples, including T=1.25 and T=2.25, have no full state",
            "raw_hamiltonian_jets": "diagnostic only; projected leading jets remain the metric",
            "work_ledger_quadrature": "trapezoid coordinate-time indicator, not an error bound",
            "continuum_or_physical_R4": "not claimed", "autonomous_renewal": "not claimed"}}
    encoded = _canon(body).encode()
    body["content_sha256"] = hashlib.sha256(encoded).hexdigest()
    body["cpu_seconds"] = time.process_time()-started
    for path, digest in sources.items():
        if episode.file_sha256(path) != digest:
            raise ValueError("evaluator changed during the read")
    return body

def _gap(left, right):
    return float(np.max(np.abs(np.asarray(left)-np.asarray(right))))

def _reference_arrays(directory, case_id):
    directory = episode.LAB/directory if not Path(directory).is_absolute() else Path(directory)
    record, arrays = authenticate_snapshot(directory, directory/f"{case_id}-000000.json")
    authenticate_producers(record)
    return record, arrays

def responsive_report(directories):
    """Recompute control stations from saved arrays. Checkpoint observation headers are not the measurement."""
    started = time.process_time()
    modules = (backend, episode, extent, leading, observer, tidal, galerkin, provenance)
    sources = {str(Path(module.__file__).resolve()): file_hash(module.__file__) for module in modules}
    script = str(Path(__file__).resolve())
    sources[script] = file_hash(script)
    magnetic_balanced, magnetic_balanced_arrays = _reference_arrays(MAGNETIC, "nf128_plus_balanced_dt0.001")
    magnetic_parent, magnetic_parent_arrays = _reference_arrays(MAGNETIC, "nf128_plus_parent_heavy_dt0.001")
    campaigns, seen = [], {}
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        for directory in directories:
            directory = Path(directory).expanduser().resolve() if Path(directory).is_absolute() else (episode.LAB/directory).resolve()
            manifest = json.loads((directory/"manifest.json").read_text())
            producers = authenticate_producers(manifest, pinned=None)
            if manifest.get("predecessor_source_hashes") and manifest.get("predecessor_producing_commit"):
                authenticate_producers({"producing_commit": manifest["predecessor_producing_commit"],
                                        "source_hashes": manifest["predecessor_source_hashes"]}, pinned=None)
            measured, audits = [], []
            for case in manifest["cases"]:
                case_id = case["case_id"]
                chunks = sorted(directory.glob(case_id+"-*.json"))
                if not chunks:
                    raise ValueError("completed case has no snapshot: "+case_id)
                handoff = None
                for json_path in chunks:
                    record, arrays = authenticate_snapshot(directory, json_path)
                    authenticate_producers(record, pinned=None)
                    digest = record["arrays_sha256"]
                    prior = seen.get((case_id, digest))
                    header = record.get("observation") if isinstance(record.get("observation"), dict) else {}
                    header_time = None if header.get("time") is None else float(header["time"])
                    audit = {"case_id": case_id, "ordinal": int(record["ordinal"]), "snapshot_kind": record["snapshot_kind"],
                             "coordinate_time": float(record["coordinate_time"]), "header_time": header_time,
                             "header_stale": header_time is not None and abs(header_time-float(record["coordinate_time"])) > 1e-8,
                             "header_used_as_measurement": False, "duplicate_of": prior}
                    audits.append(audit)
                    if handoff is None:
                        handoff = {name: np.array(arrays[name], copy=True) for name in ("Q", "r")}
                    if prior is not None:
                        continue
                    row = measure(record, arrays, producers, detail=True)
                    row["geometry_change_from_case_handoff"] = {"Q": _gap(arrays["Q"], handoff["Q"]), "r": _gap(arrays["r"], handoff["r"])}
                    measured.append(row)
                    seen[(case_id, digest)] = {"directory": directory.name, "ordinal": int(record["ordinal"])}
            campaigns.append({"directory": directory.relative_to(episode.LAB).as_posix(), "manifest_status": manifest.get("status"),
                "producing_commit": manifest.get("producing_commit"), "control_mode": manifest.get("control_mode"),
                "stations_declared": manifest.get("stations"), "stations": measured, "header_audits": audits})
    strong = [row for item in campaigns for row in item["stations"] if item["directory"].endswith("strong-t1-v1") or item["directory"].endswith("strong-t3-v1")]
    strong_handoff = next(row for row in strong if row["time"] == 0.0 and row["sign_name"] == "plus")
    frozen = [row for item in campaigns for row in item["stations"] if "frozen-parent-heavy" in item["directory"]]
    empty = [row for item in campaigns for row in item["stations"] if "empty-source" in item["directory"]]
    selected_times = (1.0, 1.25, 1.5, 2.25, 3.0)
    selected = [row for row in strong if any(abs(row["time"]-mark) < 1e-8 for mark in selected_times) and row["snapshot_kind"] in ("station", "handoff")]
    # T=1 is the completed strong-source station. Later duplicate handoffs are array copies and stay in the audit.
    selected = [row for row in selected if not (abs(row["time"]-1.0) < 1e-8 and row["snapshot_kind"] != "station")]
    selected = sorted(selected, key=lambda row: (row["sign_name"], row["time"]))
    def anchor(row):
        return {"case_id": row["case_id"], "time": row["time"], "sign_name": row["sign_name"], "control_mode": row["control_mode"],
                "parent_k": row["parent_k"], "legacy_common_k": row["common_k"],
                "child_fraction": row["child_fraction"], "child_proper_length": row["proper_lengths"]["child"],
                "child_peak": row["child_peak"], "centre_tau": row["clocks"]["centre_tau"], "Nmin": row["Nmin"],
                "projected_R4_max_abs": row["tides"]["R4"]["max_abs"], "raw_R4_max_abs": row["raw_same_point_hamiltonian_R4"]["max_abs"],
                "euler_max_abs": row["euler_trace_residual"]["max_abs"], "energy_total": row["energy"]["total"],
                "CAR": row["CAR_eigenvalues"], "ancestry_sum_minus_actual": row["ancestry"]["sum_minus_actual"],
                "header_stale": row["checkpoint_header_audit"]["stale"]}
    weight_ratio = (np.asarray(strong_handoff["ancestry"]["weights"])/magnetic_balanced_arrays["source_weights"]).tolist()
    strong_dir = next(Path(item) for item in directories if "strong-t1-v1" in str(item))
    strong_dir = strong_dir if strong_dir.is_absolute() else episode.LAB/strong_dir
    strong_record, strong_arrays = authenticate_snapshot(strong_dir, strong_dir/"nf128_plus_balanced_dt0.001-000000.json")
    attribution = {
        "strong_versus_magnetic_plus_balanced_handoff": {
            "weight_ratio": weight_ratio,
            "k": {"authoritative": "parent_k",
                  "weak": float(magnetic_balanced["parent_k"]), "strong": float(strong_record["parent_k"]),
                  "legacy_common_k": {"authoritative": False,
                                      "weak_episode": float(magnetic_balanced["common_k"]),
                                      "strong_episode": float(strong_record["common_k"]),
                                      "qualification": "weak episode common_k is a legacy label and is not the parent-record k"}},
            "initial_Q_gap": _gap(strong_arrays["Q"], magnetic_balanced_arrays["Q"]),
            "initial_r_gap": _gap(strong_arrays["r"], magnetic_balanced_arrays["r"]),
            "initial_phi_gap": max(_gap(strong_arrays["phi0"], magnetic_balanced_arrays["phi0"]),
                                   _gap(strong_arrays["phi1"], magnetic_balanced_arrays["phi1"])),
            "initial_pi_Q_gap": _gap(strong_arrays["pi_Q"], magnetic_balanced_arrays["pi_Q"]),
            "initial_pi_r_gap": _gap(strong_arrays["pi_r"], magnetic_balanced_arrays["pi_r"]),
            "initial_weight_gap": _gap(strong_arrays["source_weights"], magnetic_balanced_arrays["source_weights"]),
            "source_force_alone": False,
            "reason": "strong occupations and the solved canonical momenta both differ from the magnetic balanced handoff; authoritative parent_k also differs. The legacy weak common_k is not that comparison"},
        "frozen_versus_magnetic_plus_parent_heavy_handoff": {
            "same_canonical_handoff": True, "control": "frozen_geometry",
            "isolates_prescribed_geometry": True},
        "empty_versus_magnetic_plus_balanced_handoff": {
            "phi_deleted": True, "weights_retained": True, "momenta_re_solved": True,
            "force_deleted_from_the_coupled_trajectory": False}}
    # Fill numeric gaps from the actual handoff arrays stored on the first measured strong/frozen/empty rows by reloading is already done.
    body = {"schema": RESPONSIVE_SCHEMA, "selected_strong_observations": [anchor(row) for row in selected],
            "campaigns": campaigns, "attribution": attribution,
            "magnetic_reference_commits": magnetic_balanced["producing_commit"],
            "measurement_source_hashes": sources, "evolved": False, "reprepared": False, "resigned": False,
            "header_copied_into_measurement": False, "instability_claimed": False, "new_threshold": None,
            "late_curvature_is_an_instability": False, "physical_R4_claimed": False, "renewal_asserted": False,
            "unmeasured": {"source_force_isolated_from_k": "no completed episode changes only the source weights or only k",
                           "nf512": "not in these directories",
                           "stale_headers": "strong T3 intermediate observation headers remain the T=1 text; measurements use the arrays",
                           "ledger_quadrature": "trapezoid indicator, not a bound"}}
    return _responsive_finish(body, sources, script, started, magnetic_balanced, magnetic_balanced_arrays,
                              magnetic_parent, magnetic_parent_arrays, strong, frozen, empty, directories)

def _responsive_finish(body, sources, script, started, magnetic_balanced, magnetic_balanced_arrays,
                       magnetic_parent, magnetic_parent_arrays, strong, frozen, empty, directories):
    frozen_dir = next(Path(d) for d in directories if "frozen-parent-heavy" in str(d))
    frozen_dir = frozen_dir if frozen_dir.is_absolute() else episode.LAB/frozen_dir
    _, frozen_arrays = authenticate_snapshot(frozen_dir, next(frozen_dir.glob("*-000000.json")))
    empty_dir = next(Path(d) for d in directories if "empty-source" in str(d))
    empty_dir = empty_dir if empty_dir.is_absolute() else episode.LAB/empty_dir
    _, empty_arrays = authenticate_snapshot(empty_dir, next(empty_dir.glob("*-000000.json")))
    frozen_gap = max(_gap(frozen_arrays[name], magnetic_parent_arrays[name]) for name in ("Q", "r", "pi_Q", "pi_r", "phi0", "phi1", "source_weights"))
    body["attribution"]["frozen_versus_magnetic_plus_parent_heavy_handoff"].update(
        same_canonical_handoff=frozen_gap == 0.0, initial_state_gap=frozen_gap,
        geometry_change_at_T3=next(row["geometry_change_from_case_handoff"] for row in frozen if abs(row["time"]-3) < 1e-8))
    body["attribution"]["empty_versus_magnetic_plus_balanced_handoff"].update(
        initial_Q_gap=_gap(empty_arrays["Q"], magnetic_balanced_arrays["Q"]),
        initial_weight_gap=_gap(empty_arrays["source_weights"], magnetic_balanced_arrays["source_weights"]),
        initial_phi_max=float(max(np.max(np.abs(empty_arrays["phi0"])), np.max(np.abs(empty_arrays["phi1"])))),
        initial_pi_r_gap=_gap(empty_arrays["pi_r"], magnetic_balanced_arrays["pi_r"]))
    spans = []
    for item in body["campaigns"]:
        for case_id in sorted({row["case_id"] for row in item["stations"]}):
            group = [row for row in item["stations"] if row["case_id"] == case_id]
            first = np.asarray(group[0]["CAR_eigenvalues"])
            spans.append({"directory": item["directory"], "case_id": case_id,
                "energy_total_span": max(row["energy"]["total"] for row in group)-min(row["energy"]["total"] for row in group),
                "gram_distance_span": max(row["gram_distance_from_identity"] for row in group)-min(row["gram_distance_from_identity"] for row in group),
                "CAR_span_max": float(np.max(np.abs(np.asarray([row["CAR_eigenvalues"] for row in group])-first))),
                "field_energy_scale": max(abs(row["energy"]["field"]) for row in group)})
    body["conservation"] = spans
    encoded = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    body["content_sha256"] = hashlib.sha256(encoded).hexdigest()
    body["cpu_seconds"] = time.process_time()-started
    for path, digest in sources.items():
        if episode.file_sha256(path) != digest:
            raise ValueError("evaluator changed during the read")
    return body

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", nargs="+", type=Path)
    parser.add_argument("--responsive-controls", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.output and args.output.exists():
            raise FileExistsError("refusing to overwrite assessment output")
        if args.responsive_controls:
            directories = args.records or [episode.LAB/path for path in RESPONSIVE]
            result = responsive_report(directories)
        else:
            directories = args.records or [episode.LAB/path for path in DEFAULTS]
            result = report(directories)
        payload = json.dumps(episode.jsonable(result), indent=2, allow_nan=False)+"\n"
        if len(payload.encode()) > LIMIT:
            raise ValueError("assessment output exceeds 64 MiB")
        if args.output:
            resolved = args.output.expanduser().resolve()
            if any(resolved.is_relative_to(Path(path).expanduser().resolve()) for path in directories):
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
