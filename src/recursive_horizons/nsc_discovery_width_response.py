"""Finite physical source-width response using the existing action and Jv.

The continuous frame adapter extends the owned packet/Lowdin construction;
geometry, observer, angular inputs and occupations are fixed. Stage files are
exclusive immutable commits. Forecast precedes the independent held-out arm.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import io
import json
from pathlib import Path
import time

import numpy as np

from . import nsc_discovery_episode as episode
from . import nsc_discovery_family as family
from . import nsc_discovery_prediction as prediction
from . import nsc_discovery_response as response
from . import nsc_nested_parent_child as model
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

SCHEMA = "NSC-DISCOVERY-WIDTH-RESPONSE-v1"
NONLINEAR_SCHEMA = "NSC-DISCOVERY-WIDTH-RESPONSE-NONLINEAR-v2"
NONLINEAR_HELD_OUT_WIDTH = 1.03
ABSOLUTE_TAU_TARGET = 0.3578554631682531
CURVATURE_STAGES = ("curvature-prepare", "curvature-prediction", "curvature-measurement")
CONFIRMATION_SCHEMA = "NSC-DISCOVERY-WIDTH-NUMERICAL-CONFIRMATION-v1"
CONFIRMATION_STAGES = ("confirmation-lock", "confirmation-measurement")
CONFIRMATION_CPU_CAP = 120.0
LEGACY_LINEAR_OUTPUT = episode.LAB / "results" / "development" / "nsc-discovery-width-response-v1"
BASE_WIDTH = 1.0
HELD_OUT_WIDTH = 1.05
DERIVATIVE_H = (0.001, 0.0005)
CPU_CAP = 300.0
FORECAST_FACTOR = 1.5
CHUNK_LIMIT = 64 * 1024**2
FIELDS = model.STATE_NAMES
STAGES = ("prepare", "prediction", "measurement")


def _digest(value):
    return hashlib.sha256(json.dumps(episode.jsonable(value), sort_keys=True,
                                    separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def inputs(nf=128, duration=0.01, step_cap=episode.MATCHED_STEP_CAP,
           production=False, cpu_budget=CPU_CAP):
    nf = int(nf)
    decision = prediction.accepted_duration(duration, production=production)
    if decision["production"] and float(duration) != 0.3:
        raise ValueError("width response production target is T=0.3")
    if nf not in (128, 256):
        raise ValueError("campaign uses frozen NF128 or NF256; smaller grids are test helpers only")
    if not np.isfinite(step_cap) or float(step_cap) <= 0:
        raise ValueError("step cap must be positive and finite")
    if not np.isfinite(cpu_budget) or not 0 < float(cpu_budget) <= CPU_CAP:
        raise ValueError("aggregate CPU allowance must be positive and at most 300 seconds")
    value = {"nf": nf, "duration": float(duration), "step_cap": float(step_cap),
             "production": bool(decision["production"]), "cpu_budget_seconds": float(cpu_budget),
             "forecast_factor": FORECAST_FACTOR, "base_width": BASE_WIDTH,
             "held_out_width": HELD_OUT_WIDTH, "derivative_h": list(DERIVATIVE_H),
             "weights": family.occupation_weights(1.0).tolist(),
             "float_hex": {"duration": float(duration).hex(), "step_cap": float(step_cap).hex(),
                           "held_out_width": HELD_OUT_WIDTH.hex(),
                           "derivative_h": [h.hex() for h in DERIVATIVE_H]},
             "primary": "all-six-column probability in the fixed child (1,3)",
             "observer_reference_varied": False, "W_varied": False,
             "Q0_varied": False, "angular_inputs_varied": False,
             "momentum_representation": episode.CANONICAL_PI,
             "coordinate_pullback": False, "universal_scale_power": None,
             "fixed_inputs": {"locked_coefficients": coupling.locked_coefficients(),
                              "kappa": int(coupling.KAPPA), "multiplicity": 4*int(coupling.KAPPA),
                              "carrier_k": float(coupling.CARRIER_K),
                              "Q0": float(coupling.CALIBRATION["b0"]/coupling.CALIBRATION["a0"])}}
    value["inputs_sha256"] = _digest(value)
    return value


def source_columns(grid, width):
    """Full-precision width, fixed columns/order/phases; no catalogue matching."""
    width = float(width)
    if not np.isfinite(width) or not 0.0 < width < 2.0:
        raise ValueError("source width must lie strictly between 0 and 2")
    supports = ((0., 2.-width), (2.-width, 2.+width), (2.+width, 4.))
    norms = coupling._lobe_norms()
    phase = np.exp(1j * coupling.CALIBRATION["phase"])
    raw = np.column_stack([column for left, right in supports
                           for column in family._packet_pair(grid.xi_f, grid.dx_f, left, right, norms, phase)])
    reference0, reference1, _ = coupling.prepare_rank6(grid.xi_f, grid.dx_f)
    child = (np.vstack((reference0[:, 2:4], reference1[:, 2:4])).copy()
             if width == BASE_WIDTH else raw[:, 2:4].copy())
    child_gap = float(np.max(np.abs(child.conj().T @ child - np.eye(2))))
    child_repaired = child_gap > 1e-9
    if child_repaired:
        if width == BASE_WIDTH:
            raise ValueError("refusing to repair the owned nominal child frame")
        child, _ = coupling.lowdin(child)
    outer = np.column_stack((raw[:, :2], raw[:, 4:]))
    overlap_before = child.conj().T @ outer
    projected = outer - child @ overlap_before
    outer, outer_gram = coupling.lowdin(projected)
    columns = np.column_stack((outer[:, :2], child, outer[:, 2:]))
    phi0, phi1 = columns[:grid.nf].copy(), columns[grid.nf:].copy()
    gram = float(np.max(np.abs(columns.conj().T @ columns - np.eye(6))))
    if gram > 1e-9:
        raise ValueError("continuous source Gram defect exceeded the owned frame tolerance")
    weights = family.occupation_weights(1.)
    tails = family._tails(grid, phi0, phi1, supports)
    reference_child = np.vstack((reference0[:, 2:4], reference1[:, 2:4]))
    eigenvalues = prediction.assert_admissible_source(phi0, phi1, weights)
    meta = {"width_scale": width, "half_width": width, "width_hex": width.hex(),
            "supports": [list(x) for x in supports], "source_support_closures": [list(x) for x in supports],
            "gram_max": gram, "child_pair_lowdin_explicit": child_repaired,
            "child_frame_origin": "owned_original_middle" if width == BASE_WIDTH else "declared_continuous_width_lobes",
            "child_support": list(supports[1]), "outer_annuli": [list(supports[0]), list(supports[2])],
            "child_frame_altered": False, "child_frame_silently_altered": False,
            "middle_columns_equal_owned_original": bool(np.array_equal(child, reference_child)),
            "child_source_equals_reference": bool(np.array_equal(child, reference_child)),
            "source_equals_reference": bool(np.array_equal(columns, np.vstack((reference0, reference1)))),
            "overlap_before_max": float(np.max(np.abs(overlap_before))),
            "overlap_after_max": float(np.max(np.abs(child.conj().T @ outer))),
            "outer_gram_before_lowdin_max": float(np.max(np.abs(outer_gram-np.eye(4)))),
            "outer_gram_min": float(np.min(np.linalg.eigvalsh(outer_gram))),
            "child_included_in_outer_lowdin": False, "column_order_changed": False,
            "unequal_weight_alignment_rotations": False,
            "explicit_orthogonalization": family.EXPLICIT_ORTHOGONALIZATION,
            "tails": tails, "tail_max": float(max(tails)),
            "fine_outside_support_power_fraction": tails, "exact_spatial_support": False,
            "mean_current": family._probe_current_mean(grid, phi0, phi1), "mean_current_deleted": False,
            "CAR_min": float(np.min(eigenvalues)), "CAR_max": float(np.max(eigenvalues)),
            "occupation_trace": float(np.sum(weights))}
    meta["frame_sha256"] = episode.array_sha256(columns)
    return phi0, phi1, meta


def pair_at_width(base_pair, width):
    phi0, phi1, meta = source_columns(base_pair.grid, width)
    columns = np.vstack((phi0, phi1))
    metadata = dict(base_pair.source_metadata, **meta)
    return replace(base_pair, source_phi0=family._readonly(phi0), source_phi1=family._readonly(phi1),
                   source_columns=family._readonly(columns), source_metadata=metadata)


def frame_direction(pair, state, h, *, base_width=BASE_WIDTH):
    h = float(h)
    if h not in DERIVATIVE_H:
        raise ValueError("frame direction uses predeclared h=0.001 or 0.0005")
    plus0, plus1, plus_meta = source_columns(pair.grid, float(base_width)+h)
    minus0, minus1, minus_meta = source_columns(pair.grid, float(base_width)-h)
    tangent = response.zero_tangent(pair.grid)
    tangent.phi0 = (plus0-minus0)/(2*h)
    tangent.phi1 = (plus1-minus1)/(2*h)
    prepared = response.prepared_nested_radius_tangent(pair, state, tangent)
    if not prepared.available or prepared.delta_r is None:
        raise ValueError("implicit initial radius direction unavailable: " + str(prepared.missing_primitive))
    tangent.r = np.array(prepared.delta_r, copy=True)
    source = np.asarray(pair.source_columns)
    derivative = np.vstack((tangent.phi0, tangent.phi1))
    covariance_derivative = ((derivative * pair.weights) @ source.conj().T
                             + (source * pair.weights) @ derivative.conj().T)
    report = {"base_width": float(base_width), "h": h, "h_hex": h.hex(), "plus": plus_meta, "minus": minus_meta,
              "frame_direction_sha256": episode.array_sha256(derivative),
              "frame_direction_norm": float(np.linalg.norm(derivative)),
              "tangency_max": float(np.max(np.abs(source.conj().T @ derivative + derivative.conj().T @ source))),
              "covariance_derivative_frobenius": float(np.linalg.norm(covariance_derivative)),
              "physical_covariance_change_measured": bool(np.linalg.norm(covariance_derivative) > 0),
              "implicit_radius_linear_residual": prepared.linear_residual_max,
              "radius_jacobian_owner": prepared.jacobian_owner,
              "radius_direction_norm": float(np.linalg.norm(tangent.r)),
              "newton_used_for_tangent": prepared.newton_used,
              "frame_difference_is_analytic_width_derivative": False,
              "full_state_Jv_transport_is_analytic": True,
              "initial_delta_phi_retained": True, "delta_occupations": tangent.occupations.tolist()}
    arrays = {"plus_columns": np.vstack((plus0, plus1)), "minus_columns": np.vstack((minus0, minus1))}
    return tangent, report, arrays


def proper_geometry(pair, state, width, tangent=None, *, width_derivative=1.0):
    fine, _system = model._active(pair, state)
    density = fine.r * fine.Q
    support = (2.-float(width), 2.+float(width))
    parent = model.interval_integral(pair.grid, density, model.PARENT_INTERVAL)
    child = model.interval_integral(pair.grid, density, model.CHILD_INTERVAL)
    source = model.interval_integral(pair.grid, density, support)
    if min(parent, child, source) <= 0:
        raise ValueError("proper geometry left the positive domain")
    result = {"parent_proper_length": parent, "fixed_child_proper_length": child,
              "source_support_proper_length": source,
              "parent_to_source_support_ratio": parent/source,
              "parent_to_fixed_child_ratio": parent/child,
              "source_support_to_fixed_child_ratio": source/child,
              "source_support": list(support), "support_is_exact_packet_wall": False}
    if tangent is not None:
        nodal = prediction._nodal_tangent(pair, tangent)
        lifted = response._prolong_tangent(pair.grid, nodal)
        delta_density = lifted.r * fine.Q + fine.r * lifted.Q
        dp = model.interval_integral(pair.grid, delta_density, model.PARENT_INTERVAL)
        dc = model.interval_integral(pair.grid, delta_density, model.CHILD_INTERVAL)
        interior = model.interval_integral(pair.grid, delta_density, support)
        endpoints = float(width_derivative) * float(np.sum(model._periodic_values(pair.grid, density, support)))
        ds = interior + endpoints
        result["width_derivative_at_fixed_t"] = {
            "source_support_interior": interior, "source_support_endpoint_term": endpoints,
            "source_support_length": ds, "parent_length": dp, "fixed_child_length": dc,
            "parent_to_source_support_ratio": (dp-source**-1*parent*ds)/source,
            "parent_to_fixed_child_ratio": (dp-child**-1*parent*dc)/child}
    return result


def gaussian_report(pair, state):
    columns = np.vstack((state.phi0, state.phi1))
    values = prediction.assert_admissible_source(state.phi0, state.phi1, pair.weights)
    return {"CAR_min": float(np.min(values)), "CAR_max": float(np.max(values)),
            "column_gram_max": float(np.max(np.abs(columns.conj().T @ columns-np.eye(6)))),
            "actual_covariance_trace": float(np.sum(np.abs(columns)**2 * pair.weights[None, :])),
            "source_occupation_trace": float(np.sum(pair.weights)),
            "CAR_complement_roundoff_tolerance": 1e-8}


def _state_arrays(state, prefix):
    return {prefix+name: np.array(getattr(state, name), copy=True) for name in FIELDS}


def _state_from(arrays, prefix):
    return model.NestedState(*(np.array(arrays[prefix+name], copy=True) for name in FIELDS))


def _tangent_arrays(tangent, prefix):
    return {prefix+name: np.array(getattr(tangent, name), copy=True) for name in FIELDS+("occupations",)}


def _tangent_from(arrays, prefix):
    return response.StateTangent(*(np.array(arrays[prefix+name], copy=True) for name in FIELDS+("occupations",)))


def _pair_arrays(pair):
    return {"W": pair.geometry_map, "source_phi0": pair.source_phi0, "source_phi1": pair.source_phi1,
            "observer_columns": pair.reference_columns, "source_weights": pair.weights}


def _pair_record(pair):
    return {"nf": int(pair.grid.nf), "coarse_indices": pair.geometry_coarse_indices.tolist(),
            "child_indices": pair.geometry_child_indices.tolist(), "parent_indices": pair.geometry_parent_indices.tolist(),
            "source_metadata": pair.source_metadata, "geometry_metadata": pair.geometry_metadata,
            "clock_locations": list(pair.clock_locations)}


def _pin_identity():
    identity = family.producer_identity()
    for path in (Path(__file__).resolve(), episode.LAB / "scripts" / "derive_nsc_discovery_width_response.py"):
        relative = path.relative_to(episode.REPO).as_posix()
        digest = episode.file_sha256(path)
        commit = identity["working_tree_commit"]
        try:
            family.resolve_pinned_source_bytes(episode.REPO, relative, digest, commit=commit)
        except (ValueError, RuntimeError):
            commit = None
        identity["files"].append({"path": relative, "sha256": digest, "commit": commit})
    identity["all_sources_commit_pinned"] = all(item["commit"] for item in identity["files"])
    identity["implementation_revision"] = SCHEMA
    return identity


def _output(output):
    path = episode.assert_campaign_output(output)
    if path == family.LEGACY_DIRECTORY.resolve() or family.LEGACY_DIRECTORY.resolve() in path.parents:
        raise PermissionError("sealed family evidence is not a width-response destination")
    return path


def _new_stage(directory, stage):
    if any((Path(directory)/(stage+extension)).exists() for extension in (".json", ".npz")):
        raise FileExistsError("immutable stage already exists: " + stage)


def seal(output, stage, record, arrays):
    directory = _output(output)
    directory.mkdir(parents=True, exist_ok=True)
    json_path, npz_path = directory/(stage+".json"), directory/(stage+".npz")
    _new_stage(directory, stage)
    for name, value in arrays.items():
        if not np.isfinite(np.asarray(value)).all():
            raise ValueError("nonfinite checkpoint array: " + name)
    stream = io.BytesIO()
    np.savez_compressed(stream, **arrays)
    payload = stream.getvalue()
    record = dict(record, schema=record.get("schema", SCHEMA), stage=stage, npz=npz_path.name,
                  payload_sha256=hashlib.sha256(payload).hexdigest(),
                  array_sha256={name: episode.array_sha256(value) for name, value in arrays.items()},
                  immutable=True, chunk_limit_bytes=CHUNK_LIMIT)
    text = json.dumps(episode.jsonable(record), indent=2, allow_nan=False).encode()+b"\n"
    if max(len(payload), len(text)) > CHUNK_LIMIT:
        raise ValueError("stage exceeds 64 MiB chunk limit")
    with npz_path.open("xb") as handle:
        handle.write(payload)
    npz_path.chmod(0o444)
    with json_path.open("xb") as handle:
        handle.write(text)
    json_path.chmod(0o444)
    return record


def load(output, stage):
    directory = _output(output)
    path = directory/(stage+".json")
    record = json.loads(path.read_text())
    payload_path = directory/record["npz"]
    expected_schema = (CONFIRMATION_SCHEMA if stage in CONFIRMATION_STAGES else
                       NONLINEAR_SCHEMA if stage in CURVATURE_STAGES else SCHEMA)
    if payload_path != directory/(stage+".npz") or record["schema"] != expected_schema or record["stage"] != stage:
        raise ValueError("unexpected stage identity")
    if episode.file_sha256(payload_path) != record["payload_sha256"]:
        raise ValueError("stage payload hash mismatch")
    for source in (path, payload_path):
        if source.stat().st_size > CHUNK_LIMIT or source.stat().st_mode & 0o222:
            raise ValueError("stage is writable or exceeds chunk limit")
    with np.load(payload_path, allow_pickle=False) as stored:
        arrays = {name: np.array(stored[name], copy=True) for name in stored.files}
    if set(arrays) != set(record["array_sha256"]):
        raise ValueError("checkpoint array names differ")
    for name, value in arrays.items():
        if episode.array_sha256(value) != record["array_sha256"][name]:
            raise ValueError("checkpoint array identity mismatch: " + name)
    if _digest({k:v for k,v in record["inputs"].items() if k != "inputs_sha256"}) != record["inputs"]["inputs_sha256"]:
        raise ValueError("full-precision numerical input digest mismatch")
    return record, arrays


def _spent(record):
    return float(record.get("aggregate_cpu_seconds", 0.0))


def _budget(config, prior, start, cost=0.0, remaining_steps=0):
    spent = prior + time.process_time()-start
    projected = spent + FORECAST_FACTOR * float(cost) * int(remaining_steps)
    return spent >= config["cpu_budget_seconds"] or projected > config["cpu_budget_seconds"], spent, projected


def prepare(output, *, nf=128, duration=0.01, step_cap=episode.MATCHED_STEP_CAP,
            production=False, cpu_budget=CPU_CAP):
    directory = _output(output)
    if any((directory/(stage+".json")).exists() or (directory/(stage+".npz")).exists() for stage in STAGES):
        raise FileExistsError("width-response successor already has a stage")
    config = inputs(nf, duration, step_cap, production, cpu_budget)
    started = time.process_time()
    pair, basis_pins = family.assemble_pair(nf, 1., 1.)
    arrays = _pair_arrays(pair)
    directions, reports = [], []
    try:
        state, solver = family.solve_prepared(pair)
        arrays.update(_state_arrays(state, "initial_"))
        for index, h in enumerate(DERIVATIVE_H):
            halted, _cpu, _forecast = _budget(config, 0., started)
            if halted:
                raise TimeoutError("CPU allowance exhausted during preparation")
            tangent, report, frames = frame_direction(pair, state, h)
            directions.append(tangent); reports.append(report)
            arrays.update(_tangent_arrays(tangent, f"d{index}_"))
            arrays.update({f"d{index}_"+name: value for name, value in frames.items()})
    except (ValueError, RuntimeError, TimeoutError, FloatingPointError) as error:
        failed = {"inputs": config, "pair": _pair_record(pair), "basis_pins": basis_pins,
                  "producer_identity": _pin_identity(),
                  "status": "budget_stop" if isinstance(error, TimeoutError) else "initial_preparation_failed",
                  "error": type(error).__name__+": "+str(error), "directions": reports,
                  "non_existence_claimed": False, "physical_conclusion_forced": False,
                  "held_out_source_not_constructed": True, "aggregate_cpu_seconds": float(time.process_time()-started)}
        return seal(directory, "prepare", failed, arrays)
    frame_gap = np.linalg.norm(np.vstack((directions[0].phi0-directions[1].phi0,
                                          directions[0].phi1-directions[1].phi1)))
    elapsed = time.process_time()-started
    record = {"inputs": config, "pair": _pair_record(pair), "basis_pins": basis_pins,
              "producer_identity": _pin_identity(), "status": "prepared" if elapsed < cpu_budget else "budget_stop",
              "aggregate_cpu_seconds": float(time.process_time()-started),
              "solver": solver, "directions": reports, "frame_direction_h_difference_norm": float(frame_gap),
              "initial_geometry": [proper_geometry(pair, state, BASE_WIDTH, tangent) for tangent in directions],
              "fixed_Q0_preparation_limitation": "Q0=b0/a0; chi and canonical momenta initially zero; only source columns and implicit radius vary",
              "initial_state_sha256": episode.state_sha256(state),
              "held_out_source_not_constructed": True, "held_out_measurement_performed": False,
              "forecast_factor": FORECAST_FACTOR, "pool_count": 0,
              "budget_limit": "CPU checked between operations; dense preparation cannot be interrupted mid-solve"}
    record["aggregate_cpu_seconds"] = float(time.process_time()-started)
    if record["aggregate_cpu_seconds"] >= config["cpu_budget_seconds"]:
        record["status"] = "budget_stop"
    return seal(directory, "prepare", record, arrays)


def predict(output):
    directory = _output(output)
    _new_stage(directory, "prediction")
    started = time.process_time()
    prepared, saved = load(directory, "prepare")
    family.authenticate_producer({"producer_identity": prepared["producer_identity"]})
    if prepared["status"] != "prepared":
        raise ValueError("preparation stopped; use a new successor")
    config = prepared["inputs"]
    pair = episode.pair_from_arrays(saved, prepared["pair"])
    state = _state_from(saved, "initial_")
    tangents = [_tangent_from(saved, f"d{i}_") for i in range(2)]
    pair, state, backend = prediction._resolve(pair, state, "fft")
    mark, tau, delta_tau, steps = 0., 0., [0., 0.], 0
    cost = 0.; status = "predicted"; stop = None; projected = _spent(prepared)
    while mark < config["duration"]-1e-13:
        limit, _omega, _quadrature = model.stable_timestep(pair, state, config["step_cap"])
        dt = min(float(limit), config["duration"]-mark)
        halted, _cpu, projected = _budget(config, _spent(prepared), started, cost,
                                         int(np.ceil((config["duration"]-mark)/dt)))
        if halted:
            status, stop = "budget_stop", {"forecast_cpu_seconds": projected}; break
        batch_start = time.process_time()
        try:
            results = [prediction.coupled_rk4_step(pair, state, tangent, dt) for tangent in tangents]
        except (coupling.PositiveChartExit, FloatingPointError, ValueError) as error:
            status, stop = "chart_or_numerical_stop", {"error": type(error).__name__+": "+str(error)}; break
        if any(not np.isfinite(getattr(result[kind], name)).all()
               for result in results for kind in ("state", "tangent") for name in FIELDS):
            status, stop = "chart_or_numerical_stop", {"error": "nonfinite full-state or Jv output"}
            break
        if any(not np.array_equal(getattr(results[0]["state"], name), getattr(results[1]["state"], name)) for name in FIELDS):
            raise RuntimeError("two Jv directions changed the shared baseline differently")
        state = results[0]["state"]
        tau += results[0]["tau"]
        for i, result in enumerate(results):
            delta_tau[i] += result["delta_tau"]
            tangents[i] = result["tangent"]
        mark += dt; steps += 1
        cost = max(cost, time.process_time()-batch_start)
    coefficients = []
    for index, tangent in enumerate(tangents):
        readout = prediction.matched_readout(pair, state, tangent, tau=tau,
                                             delta_tau=delta_tau[index], coordinate_time=mark)
        coefficient = readout["delta_child_regional_content_tau"]
        coefficients.append({"h": DERIVATIVE_H[index], "coefficient_at_equal_tau": coefficient,
                             "readout": readout,
                             "held_out_prediction": readout["child_regional_content"] + (HELD_OUT_WIDTH-BASE_WIDTH)*coefficient,
                             "proper_geometry": proper_geometry(pair, state, BASE_WIDTH, tangent)})
    eigenvalues = prediction.assert_admissible_source(state.phi0, state.phi1, pair.weights)
    arrays = _pair_arrays(pair); arrays.update(_state_arrays(state, "baseline_"))
    arrays.update(_state_arrays(_state_from(saved, "initial_"), "initial_"))
    for i, tangent in enumerate(tangents):
        arrays.update(_tangent_arrays(tangent, f"d{i}_"))
    record = {"inputs": config, "pair": prepared["pair"], "producer_identity": prepared["producer_identity"],
              "prepare_json_sha256": episode.file_sha256(directory/"prepare.json"),
              "status": status, "stop": stop, "coordinate_time": mark, "target_tau": tau,
              "steps": steps, "coefficients": coefficients,
              "coefficient_h_difference": abs(coefficients[0]["coefficient_at_equal_tau"]-coefficients[1]["coefficient_at_equal_tau"]),
              "held_out_frame_difference_indicator": abs(coefficients[0]["held_out_prediction"]-coefficients[1]["held_out_prediction"]),
              "aggregate_cpu_seconds": _spent(prepared)+float(time.process_time()-started),
              "last_batch_cpu_seconds": cost, "forecast_cpu_seconds": projected,
              "backend": backend, "CAR_min": float(np.min(eigenvalues)), "CAR_max": float(np.max(eigenvalues)),
              "gaussian_state": gaussian_report(pair, state),
              "W_sha256": episode.array_sha256(pair.geometry_map),
              "observer_sha256": episode.array_sha256(pair.reference_columns),
              "initial_delta_phi_retained": True, "delta_weights_zero": True,
              "forecast_locked_before_measurement": True, "held_out_measurement_performed": False,
              "continuum_certified": False, "universal_scale_power": None,
              "space_time_evolution_error_bound": None,
              "method": "existing coupled_rk4_step and analytic full-state Jv; numerical source-frame width directions"}
    record["aggregate_cpu_seconds"] = _spent(prepared)+float(time.process_time()-started)
    if record["aggregate_cpu_seconds"] >= config["cpu_budget_seconds"]:
        record["status"] = "budget_stop"
    return seal(directory, "prediction", record, arrays)


def measure(output):
    directory = _output(output)
    _new_stage(directory, "measurement")
    started = time.process_time()
    forecast, saved = load(directory, "prediction")
    prepared, _prepared_arrays = load(directory, "prepare")
    if forecast["prepare_json_sha256"] != episode.file_sha256(directory/"prepare.json"):
        raise ValueError("prediction preparation identity changed")
    family.authenticate_producer({"producer_identity": forecast["producer_identity"]})
    if forecast["status"] != "predicted":
        raise ValueError("forecast did not reach target; held-out arm not launched")
    config = forecast["inputs"]
    initial_cost = float(prepared["solver"].get("initial_cpu_seconds", 0.0))
    halted, _cpu, projected = _budget(config, _spent(forecast), started, initial_cost, 1)
    if halted:
        record = {"inputs": config, "producer_identity": forecast["producer_identity"],
                  "prediction_json_sha256": episode.file_sha256(directory/"prediction.json"),
                  "status": "budget_stop", "held_out_constructor_called": False,
                  "forecast_cpu_seconds": projected, "physical_conclusion_forced": False,
                  "aggregate_cpu_seconds": _spent(forecast)+float(time.process_time()-started)}
        return seal(directory, "measurement", record, _state_arrays(_state_from(saved, "baseline_"), "baseline_"))
    dense = episode.pair_from_arrays(saved, forecast["pair"])
    changed = pair_at_width(dense, HELD_OUT_WIDTH)
    # The independent source constructor is called only after prediction.json is sealed.
    try:
        state, solver = family.solve_prepared(changed)
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        record = {"inputs": config, "producer_identity": forecast["producer_identity"],
                  "prediction_json_sha256": episode.file_sha256(directory/"prediction.json"),
                  "status": "initial_preparation_failed", "error": type(error).__name__+": "+str(error),
                  "non_existence_claimed": False, "physical_conclusion_forced": False,
                  "aggregate_cpu_seconds": _spent(forecast)+float(time.process_time()-started)}
        return seal(directory, "measurement", record, {"held_source_columns": changed.source_columns})
    initial = state.copy()
    initial_geometry = proper_geometry(changed, state, HELD_OUT_WIDTH)
    pair, state, backend = prediction._resolve(changed, state, "fft")
    mark, tau, steps, cost = 0., 0., 0, 0.
    target, ceiling = forecast["target_tau"], 1.25*config["duration"]
    if not config["production"]:
        ceiling = min(ceiling, prediction.NONPRODUCTION_DURATION)
    status, stop, bracket, event = "tau_unbracketed", None, None, None
    while tau < target and mark < ceiling-1e-13:
        limit, _omega, _quadrature = model.stable_timestep(pair, state, config["step_cap"])
        dt = min(float(limit), ceiling-mark)
        halted, _cpu, projected = _budget(config, _spent(forecast), started, cost,
                                         int(np.ceil((ceiling-mark)/dt))+prediction.ROOT_BISECTIONS)
        if halted:
            status, stop = "budget_stop", {"forecast_cpu_seconds": projected}; break
        batch_start = time.process_time()
        previous = state.copy(); previous_tau = tau; previous_mark = mark
        try:
            result = prediction.nonlinear_rk4_step(pair, state, dt)
        except (coupling.PositiveChartExit, FloatingPointError, ValueError) as error:
            status, stop = "chart_or_numerical_stop", {"error": type(error).__name__+": "+str(error)}; break
        if any(not np.isfinite(getattr(result["state"], name)).all() for name in FIELDS):
            status, stop = "chart_or_numerical_stop", {"error": "nonfinite held-out state output"}
            break
        state = result["state"]; tau += result["tau"]; mark += dt; steps += 1
        cost = max(cost, time.process_time()-batch_start)
        if previous_tau < target <= tau:
            bracket = {"state": previous, "tau": previous_tau, "time": previous_mark,
                       "dt_hi": dt, "tau_hi": tau}
            halted, _cpu, projected = _budget(config, _spent(forecast), started, cost, prediction.ROOT_BISECTIONS)
            if halted:
                status, stop = "budget_stop", {"forecast_cpu_seconds": projected}; break
            try:
                event = prediction._root_step_to_tau(pair, bracket, target)
            except (coupling.PositiveChartExit, ValueError, FloatingPointError) as error:
                status, stop = "clock_root_failed", {"error": type(error).__name__+": "+str(error)}
                break
            state, tau, mark = event["state"], event["tau"], event["time"]
            status = "measured_at_equal_tau"
            break
    geometry = proper_geometry(pair, state, HELD_OUT_WIDTH)
    observed = response.child_regional_content(pair.grid, model.reconstruct_state(pair, state))
    base = forecast["coefficients"][-1]["readout"]["child_regional_content"]
    predicted = forecast["coefficients"][-1]["held_out_prediction"]
    clock_indicator = None; interpolation_indicator = None
    if event is not None:
        zero = response.zero_tangent(pair.grid)
        final_readout = prediction.matched_readout(pair, state, zero, tau=tau, delta_tau=0., coordinate_time=mark)
        clock_indicator = abs(final_readout["child_regional_content_dot"] * event["tau_residual"] / final_readout["tau_dot"])
        low = prediction._observables(pair, bracket["state"])["child_regional_content"]
        high_state = prediction.nonlinear_rk4_step(pair, bracket["state"], bracket["dt_hi"])["state"]
        high = prediction._observables(pair, high_state)["child_regional_content"]
        fraction = (target-bracket["tau"])/(bracket["tau_hi"]-bracket["tau"])
        interpolation_indicator = abs(low+fraction*(high-low)-observed)
    physical_effect = observed-base if event is not None else None
    available_indicator = forecast["held_out_frame_difference_indicator"]
    if clock_indicator is not None:
        available_indicator += clock_indicator + interpolation_indicator
    arrays = _pair_arrays(pair); arrays.update(_state_arrays(initial, "initial_")); arrays.update(_state_arrays(state, "held_"))
    if bracket is not None:
        arrays.update(_state_arrays(bracket["state"], "bracket_"))
    eigenvalues = prediction.assert_admissible_source(state.phi0, state.phi1, pair.weights)
    record = {"inputs": config, "pair": _pair_record(changed), "producer_identity": forecast["producer_identity"],
              "prediction_json_sha256": episode.file_sha256(directory/"prediction.json"),
              "prepare_json_sha256": episode.file_sha256(directory/"prepare.json"),
              "status": status, "stop": stop, "coordinate_time": mark, "proper_time": tau,
              "target_tau": target, "tau_residual": tau-target, "event_ceiling": ceiling,
              "event": None if event is None else {k:v for k,v in event.items() if k != "state"},
              "steps": steps, "held_out_solver": solver, "initial_geometry": initial_geometry,
              "final_geometry": geometry, "observed_child_fixed_window_content": observed,
              "predicted_child_fixed_window_content": predicted,
              "observed_effect": physical_effect, "predicted_effect": predicted-base,
              "signed_nonlinear_remainder": None if event is None else observed-predicted,
              "frame_difference_indicator": forecast["held_out_frame_difference_indicator"],
              "clock_event_indicator": clock_indicator, "bracket_interpolation_indicator": interpolation_indicator,
              "available_numerical_indicator_sum": available_indicator,
              "effect_exceeds_available_indicators": None if physical_effect is None else abs(physical_effect)>available_indicator,
              "uncertainty_status": "space/time evolution refinement not measured; available indicators are not an error bound",
              "initial_radius_copied": False, "source_field_reset": False,
              "W_unchanged": episode.array_sha256(pair.geometry_map)==forecast["W_sha256"],
              "observer_unchanged": episode.array_sha256(pair.reference_columns)==forecast["observer_sha256"],
              "CAR_min": float(np.min(eigenvalues)), "CAR_max": float(np.max(eigenvalues)),
              "gaussian_state": gaussian_report(pair, state),
              "aggregate_cpu_seconds": _spent(forecast)+float(time.process_time()-started),
              "physical_conclusion_forced": False, "universal_percent_gate": None,
              "universal_scale_power": None, "continuum_certified": False, "backend": backend}
    record["aggregate_cpu_seconds"] = _spent(forecast)+float(time.process_time()-started)
    if record["aggregate_cpu_seconds"] >= config["cpu_budget_seconds"]:
        record["status"] = "budget_stop"
    return seal(directory, "measurement", record, arrays)


def check(output):
    directory = _output(output)
    if any((directory/(stage+".json")).exists() for stage in CONFIRMATION_STAGES):
        return check_confirmation(directory)
    if any((directory/(stage+".json")).exists() for stage in CURVATURE_STAGES):
        return check_curvature(directory)
    checked = []; previous = None; reports = {}
    for stage in STAGES:
        if not (directory/(stage+".json")).exists():
            if (directory/(stage+".npz")).exists():
                raise ValueError("uncommitted stage payload: " + stage)
            continue
        if stage != "prepare" and any(not (directory/(earlier+".json")).exists() for earlier in STAGES[:STAGES.index(stage)]):
            raise ValueError("stage predecessor missing: " + stage)
        record, arrays = load(directory, stage)
        identity = family.authenticate_producer({"producer_identity": record["producer_identity"]})
        if previous is not None:
            key = "prepare_json_sha256" if stage == "prediction" else "prediction_json_sha256"
            if record[key] != episode.file_sha256(directory/(previous+".json")):
                raise ValueError("stage predecessor identity changed")
            if record["inputs"] != reports[previous]["inputs"]:
                raise ValueError("stage numerical inputs changed")
        if record.get("aggregate_cpu_seconds", 0) > record["inputs"]["cpu_budget_seconds"] and record["status"] not in ("budget_stop", "initial_preparation_failed"):
            raise ValueError("completed stage exceeded the aggregate CPU allowance")
        if "W" in arrays:
            if episode.array_sha256(arrays["W"]) != episode.array_sha256(load(directory, "prepare")[1]["W"]):
                raise ValueError("stage changed W")
            original = load(directory, "prepare")[1]
            for name in ("observer_columns", "source_weights"):
                if not np.array_equal(arrays[name], original[name]):
                    raise ValueError("stage changed fixed " + name)
            frozen = family.load_frozen_basis(record["inputs"]["nf"])
            if not np.array_equal(arrays["W"], frozen["W"]):
                raise ValueError("stage W differs from frozen portable basis")
            prediction.assert_admissible_source(arrays["source_phi0"], arrays["source_phi1"], arrays["source_weights"])
            if stage == "prepare" and "initial_phi0" in arrays:
                if not np.array_equal(arrays["initial_phi0"], arrays["source_phi0"]) or not np.array_equal(arrays["initial_phi1"], arrays["source_phi1"]):
                    raise ValueError("prepared initial field differs from its source")
                if "initial_state_sha256" in record and episode.state_sha256(_state_from(arrays, "initial_")) != record["initial_state_sha256"]:
                    raise ValueError("initial physical state hash differs")
            if stage == "measurement":
                for name in ("phi0", "phi1"):
                    if not np.array_equal(arrays["initial_"+name], arrays["source_"+name]):
                        raise ValueError("held-out initial field differs from its own source")
        checked.append({"stage": stage, "status": record["status"], "producer": identity["status"]})
        previous = stage; reports[stage] = record
    if not checked:
        raise FileNotFoundError("no committed width-response stages")
    return {"schema": SCHEMA, "ok": True, "checked": checked, "evolved": False, "repaired": False}


def _nonlinear_output(output):
    directory = _output(output)
    if directory == LEGACY_LINEAR_OUTPUT.resolve() or LEGACY_LINEAR_OUTPUT.resolve() in directory.parents:
        raise PermissionError("linear-v1 evidence is immutable; nonlinear response needs a new successor")
    if directory.is_relative_to(episode.REPO) and not any(part.startswith("nsc-discovery-width-response-nonlinear-") for part in directory.parts):
        raise PermissionError("nonlinear evidence needs its own nsc-discovery-width-response-nonlinear- successor")
    return directory


def cached_baseline(directory, *, production=False, expected_nf=None):
    """Authenticate only old preparation/prediction; never read its held-out data."""
    prepared, _initial = load(directory, "prepare")
    forecast, arrays = load(directory, "prediction")
    if forecast["prepare_json_sha256"] != episode.file_sha256(Path(directory)/"prepare.json"):
        raise ValueError("cached baseline predecessor changed")
    family.authenticate_producer({"producer_identity": forecast["producer_identity"]})
    if forecast["status"] != "predicted":
        raise ValueError("cached baseline did not complete")
    nf = int(forecast["inputs"]["nf"])
    if expected_nf is not None and int(expected_nf) != nf:
        raise ValueError("cached baseline resolution differs")
    if production and abs(float(forecast["target_tau"])-ABSOLUTE_TAU_TARGET) > 1e-12:
        raise ValueError("cached baseline is not at the declared absolute proper time; redo baseline")
    if not np.array_equal(arrays["W"], family.load_frozen_basis(nf)["W"]):
        raise ValueError("cached baseline changed the frozen complete W")
    current_fixed = inputs(nf)["fixed_inputs"]
    if forecast["inputs"].get("fixed_inputs") != current_fixed:
        raise ValueError("cached action/angular/initial-Q inputs differ")
    if not np.array_equal(arrays["source_weights"], family.occupation_weights(1.)):
        raise ValueError("cached baseline source weights differ")
    pair = episode.pair_from_arrays(arrays, forecast["pair"])
    actual = np.vstack(source_columns(pair.grid, BASE_WIDTH)[:2])
    if not np.array_equal(actual, np.vstack((arrays["source_phi0"], arrays["source_phi1"]))):
        raise ValueError("cached nominal source differs from the continuous frame owner")
    return prepared, forecast, arrays


def nonlinear_preflight(baseline):
    prepared, forecast, _arrays = cached_baseline(baseline, production=True)
    maximum = float(forecast.get("last_batch_cpu_seconds", 0.))/2.
    steps = int(forecast["steps"])
    average = max(0., _spent(forecast)-_spent(prepared))/max(2*steps, 1)
    initial = float(prepared["solver"].get("initial_cpu_seconds", 0.))
    neighbor = FORECAST_FACTOR*4*steps*average
    # One full Jv step is a conservative proxy for the cheaper nonlinear arm.
    held = FORECAST_FACTOR*steps*average
    return {"schema": NONLINEAR_SCHEMA, "readonly": True, "bytes_written": 0,
            "neighbor_baselines": [BASE_WIDTH+sign*h for h in DERIVATIVE_H for sign in (-1, 1)],
            "absolute_tau_target": ABSOLUTE_TAU_TARGET, "new_held_out_width": NONLINEAR_HELD_OUT_WIDTH,
            "cached_average_one_direction_step_cpu_indicator": average,
            "cached_maximum_one_direction_step_cpu_indicator": maximum,
            "four_neighbor_forecast_cpu_indicator": neighbor,
            "held_arm_conservative_cpu_indicator": held,
            "four_dense_solves_cpu_indicator": FORECAST_FACTOR*4*initial,
            "fresh_batch_forecast_cpu_indicator": neighbor+held+FORECAST_FACTOR*4*initial,
            "nominal_coordinate_steps": steps,
            "extended_event_ceiling_maximum_cost_indicator": FORECAST_FACTOR*4*int(np.ceil(.375/forecast["inputs"]["step_cap"]))*maximum,
            "forecast_is_indicator_not_production_admission": True,
            "aggregate_cpu_cap": CPU_CAP, "old_width1p05_measurement_read": False}


def root_retarded_event(pair, bracket, target, *, stepper=None, admission=None):
    """Root the clock using the same RK step for state, tangent and clock tangent.

    This event adapter changes no RHS. Every trial starts with the bracket's
    actual state/tangent and integrated clocks, never a reset source frame.
    """
    stepper = prediction.coupled_rk4_step if stepper is None else stepper
    lo, hi = 0., float(bracket["dt_hi"])
    best = None; count = 0
    for _index in range(prediction.ROOT_BISECTIONS):
        if admission is not None and not admission(prediction.ROOT_BISECTIONS-_index):
            raise TimeoutError("CPU forecast disallows retarded clock root")
        mid = 0.5*(lo+hi)
        trial = stepper(pair, bracket["state"], bracket["tangent"], mid)
        tau = float(bracket["tau"])+float(trial["tau"])
        count += 1
        if tau < float(target):
            lo = mid
        else:
            hi = mid
            best = {"state": trial["state"], "tangent": trial["tangent"],
                    "tau": tau, "delta_tau": float(bracket["delta_tau"])+float(trial["delta_tau"]),
                    "time": float(bracket["time"])+mid, "dt": mid}
        if hi-lo <= float(bracket["dt_hi"])*1e-12:
            break
    if best is None:
        # The upper endpoint is a computed crossing, not a fabricated event.
        trial = stepper(pair, bracket["state"], bracket["tangent"], float(bracket["dt_hi"]))
        tau = float(bracket["tau"])+float(trial["tau"])
        if tau < target:
            raise ValueError("retarded clock bracket has no crossing")
        best = {"state": trial["state"], "tangent": trial["tangent"], "tau": tau,
                "delta_tau": float(bracket["delta_tau"])+float(trial["delta_tau"]),
                "time": float(bracket["time"])+float(bracket["dt_hi"]), "dt": float(bracket["dt_hi"])}
    best.update(tau_residual=best["tau"]-float(target), bracket_width=hi-lo, bisections=count,
                tangent_carried_on_same_rooted_step=True, uses_linear_clock_correction_as_event=False)
    return best


def march_retarded(pair, state, tangent, config, *, target_tau, prior_cpu, stage_start, with_tangent=True):
    pair, state, info = prediction._resolve(pair, state, "fft")
    mark, tau, delta_tau, steps, cost = 0., 0., 0., 0, 0.
    ceiling = 1.25*config["duration"] if target_tau is not None else config["duration"]
    if not config["production"]:
        ceiling = min(ceiling, prediction.NONPRODUCTION_DURATION)
    status, stop, event = "tau_unbracketed", None, None
    def stepper(pair, state, tangent, dt):
        if with_tangent:
            return prediction.coupled_rk4_step(pair, state, tangent, dt)
        result = prediction.nonlinear_rk4_step(pair, state, dt)
        return dict(result, tangent=tangent)
    while mark < ceiling-1e-13:
        limit, _omega, _quad = model.stable_timestep(pair, state, config["step_cap"])
        dt = min(float(limit), ceiling-mark)
        roots = prediction.ROOT_BISECTIONS if target_tau is not None else 0
        halted, _spent_cpu, projected = _budget(config, prior_cpu, stage_start, cost,
                                               int(np.ceil((ceiling-mark)/float(limit)))+roots)
        if halted:
            status, stop = "budget_stop", {"forecast_cpu_seconds": projected}; break
        bracket = {"state": state, "tangent": tangent, "tau": tau, "delta_tau": delta_tau,
                   "time": mark, "dt_hi": dt}
        before = time.process_time()
        try:
            advanced = stepper(pair, state, tangent, dt)
            if any(not np.isfinite(getattr(advanced[kind], name)).all() for kind in ("state", "tangent") for name in FIELDS):
                raise FloatingPointError("nonfinite retarded step")
        except (ValueError, FloatingPointError, coupling.PositiveChartExit) as error:
            status, stop = "chart_or_numerical_stop", {"error": type(error).__name__+": "+str(error)}; break
        cost = max(cost, time.process_time()-before)
        upper_tau = tau+advanced["tau"]
        if target_tau is not None and tau < target_tau <= upper_tau:
            def admission(remaining):
                return not _budget(config, prior_cpu, stage_start, cost, remaining)[0]
            try:
                event = root_retarded_event(pair, bracket, target_tau, stepper=stepper, admission=admission)
            except (TimeoutError, ValueError, coupling.PositiveChartExit) as error:
                status, stop = ("budget_stop" if isinstance(error, TimeoutError) else "clock_root_failed"), {"error": str(error)}
                break
            state, tangent = event["state"], event["tangent"]
            mark, tau, delta_tau = event["time"], event["tau"], event["delta_tau"]
            steps += 1; status = "event_reached"; break
        state, tangent = advanced["state"], advanced["tangent"]
        mark += dt; tau = upper_tau; delta_tau += advanced["delta_tau"]; steps += 1
    if target_tau is None and mark >= config["duration"]-1e-13:
        status = "coordinate_baseline_reached"
    readout = prediction.matched_readout(pair, state, tangent, tau=tau, delta_tau=delta_tau, coordinate_time=mark)
    return {"pair": pair, "state": state, "tangent": tangent, "time": mark, "tau": tau,
            "delta_tau": delta_tau, "status": status, "stop": stop, "event": event,
            "readout": readout, "steps": steps, "backend": info, "max_step_cpu": cost,
            "coordinate_ceiling": ceiling}


def curvature_from_slopes(slopes, h):
    """Centered derivative of full first retarded coefficients, not a Hessian."""
    h = float(h)
    return (float(slopes[1])-float(slopes[0]))/(2*h)


def prepare_curvature(output, *, baseline=None, nf=None, duration=0.01,
                      step_cap=episode.MATCHED_STEP_CAP, production=False, cpu_budget=CPU_CAP,
                      redo_baseline=False):
    directory = _nonlinear_output(output)
    if any((directory/(stage+suffix)).exists() for stage in CURVATURE_STAGES for suffix in (".json", ".npz")):
        raise FileExistsError("nonlinear successor already has a stage")
    start = time.process_time(); cached = None
    if baseline is not None and not redo_baseline:
        cached = cached_baseline(baseline, production=production, expected_nf=nf)
        nf = cached[1]["inputs"]["nf"]
    nf = int(128 if nf is None else nf)
    config = inputs(nf, duration, step_cap, production, cpu_budget)
    target = ABSOLUTE_TAU_TARGET if production else (None if cached is None else cached[1]["target_tau"])
    config.update(experiment="nonlinear-v2", held_out_width=NONLINEAR_HELD_OUT_WIDTH,
                  outer_derivative_h=list(DERIVATIVE_H), inner_frame_h=DERIVATIVE_H[-1],
                  common_absolute_tau_target=target,
                  old_width1p05_used_for_fit=False, second_response_method="finite differences of full retarded slopes at a common clock event")
    config["float_hex"]["held_out_width"] = NONLINEAR_HELD_OUT_WIDTH.hex()
    config["inputs_sha256"] = _digest({k:v for k,v in config.items() if k != "inputs_sha256"})
    pair, pins = family.assemble_pair(nf, 1., 1.)
    arrays = _pair_arrays(pair); nominal = {"cached": cached is not None}
    if cached is not None:
        _old_prepare, old, saved = cached
        for name in FIELDS:
            arrays["nominal_"+name] = saved["baseline_"+name]
        for name in FIELDS+("occupations",):
            arrays["nominal_d_"+name] = saved["d1_"+name]
        nominal.update(readout=old["coefficients"][-1]["readout"],
                       first_coefficient=old["coefficients"][-1]["coefficient_at_equal_tau"],
                       first_coefficient_h_indicator=old["coefficient_h_difference"],
                       target_tau=old["target_tau"], coordinate_time=old["coordinate_time"],
                       baseline_directory=str(Path(baseline).resolve()),
                       prediction_json_sha256=episode.file_sha256(Path(baseline)/"prediction.json"),
                       prediction_npz_sha256=old["payload_sha256"],
                       producing_identity=old["producer_identity"], old_measurement_read=False)
    else:
        try:
            state, solver = family.solve_prepared(pair)
            arrays.update(_state_arrays(state, "nominal_initial_"))
            tangent, report, _frames = frame_direction(pair, state, DERIVATIVE_H[-1])
            arrays.update(_tangent_arrays(tangent, "nominal_d_"))
            nominal.update(solver=solver, initial_direction=report, first_coefficient_h_indicator=None)
        except (ValueError, RuntimeError) as error:
            record = {"schema": NONLINEAR_SCHEMA, "inputs": config, "pair": _pair_record(pair),
                      "nominal": nominal, "status": "initial_preparation_failed",
                      "error": type(error).__name__+": "+str(error), "producer_identity": _pin_identity(),
                      "non_existence_claimed": False, "new_held_out_source_constructed": False,
                      "aggregate_cpu_seconds": float(time.process_time()-start)}
            return seal(directory, CURVATURE_STAGES[0], record, arrays)
    branches = []; failures = []
    for index, (outer_h, sign) in enumerate((h, s) for h in DERIVATIVE_H for s in (-1, 1)):
        center = BASE_WIDTH+sign*outer_h
        try:
            if _budget(config, 0., start)[0]:
                raise TimeoutError("CPU allowance exhausted during neighbor preparation")
            changed = pair_at_width(pair, center)
            state, solver = family.solve_prepared(changed)
            tangent, report, frames = frame_direction(changed, state, DERIVATIVE_H[-1], base_width=center)
            coarse0, coarse1, _m0 = source_columns(pair.grid, center+DERIVATIVE_H[0])
            low0, low1, _m1 = source_columns(pair.grid, center-DERIVATIVE_H[0])
            coarse = np.vstack(((coarse0-low0)/(2*DERIVATIVE_H[0]), (coarse1-low1)/(2*DERIVATIVE_H[0])))
            fine = np.vstack((tangent.phi0, tangent.phi1))
            prefix = f"b{index}_"
            arrays.update(_state_arrays(state, prefix)); arrays.update(_tangent_arrays(tangent, prefix+"d_"))
            arrays[prefix+"source_phi0"] = changed.source_phi0; arrays[prefix+"source_phi1"] = changed.source_phi1
            arrays.update({prefix+key: value for key,value in frames.items()})
            branches.append({"index": index, "width": center, "width_hex": center.hex(), "sign": sign,
                             "outer_h": outer_h, "solver": solver, "pair": _pair_record(changed),
                             "initial_direction": report, "frame_h_comparison_norm": float(np.linalg.norm(coarse-fine)),
                             "initial_state_sha256": episode.state_sha256(state),
                             "initial_geometry": proper_geometry(changed, state, center, tangent)})
        except (ValueError, RuntimeError, TimeoutError) as error:
            failures.append({"index": index, "width": center, "error": type(error).__name__+": "+str(error),
                             "non_existence_claimed": False})
    record = {"schema": NONLINEAR_SCHEMA, "inputs": config, "pair": _pair_record(pair), "basis_pins": pins,
              "nominal": nominal, "branches": branches, "failures": failures,
              "status": "prepared" if len(branches) == 4 else "preparation_incomplete",
              "producer_identity": _pin_identity(), "aggregate_cpu_seconds": float(time.process_time()-start),
              "new_held_out_source_constructed": False, "old_width1p05_used_for_fit": False,
              "fixed_Q0_preparation_limitation": "Q0 fixed; chi/momenta initially zero; each neighbor has its own solved radius"}
    if record["aggregate_cpu_seconds"] >= config["cpu_budget_seconds"]:
        record["status"] = "budget_stop"
    return seal(directory, CURVATURE_STAGES[0], record, arrays)


def predict_curvature(output):
    directory = _nonlinear_output(output); _new_stage(directory, CURVATURE_STAGES[1]); start = time.process_time()
    prepared, saved = load(directory, CURVATURE_STAGES[0])
    family.authenticate_producer({"producer_identity": prepared["producer_identity"]})
    if prepared["status"] != "prepared":
        raise ValueError("nonlinear preparation incomplete; no forecast")
    config = prepared["inputs"]; pair = episode.pair_from_arrays(saved, prepared["pair"])
    target = config["common_absolute_tau_target"]; nominal = dict(prepared["nominal"]); arrays = _pair_arrays(pair)
    status = "predicted"; rows = []
    if not nominal["cached"]:
        state = _state_from(saved, "nominal_initial_"); tangent = _tangent_from(saved, "nominal_d_")
        result = march_retarded(pair, state, tangent, config, target_tau=target,
                                prior_cpu=_spent(prepared), stage_start=start)
        if result["status"] not in ("event_reached", "coordinate_baseline_reached"):
            status = result["status"]
        target = result["tau"] if target is None else target
        nominal.update(readout=result["readout"], first_coefficient=result["readout"]["delta_child_regional_content_tau"],
                       target_tau=target, coordinate_time=result["time"], status=result["status"])
        arrays.update(_state_arrays(result["state"], "nominal_")); arrays.update(_tangent_arrays(result["tangent"], "nominal_d_"))
    else:
        arrays.update({name: value for name,value in saved.items() if name.startswith("nominal_")})
    if status == "predicted":
        for branch in prepared["branches"]:
            index = branch["index"]; prefix = f"b{index}_"
            branch_arrays = dict(saved, source_phi0=saved[prefix+"source_phi0"], source_phi1=saved[prefix+"source_phi1"])
            changed = episode.pair_from_arrays(branch_arrays, branch["pair"])
            result = march_retarded(changed, _state_from(saved, prefix), _tangent_from(saved, prefix+"d_"), config,
                                    target_tau=target, prior_cpu=_spent(prepared), stage_start=start)
            arrays.update(_state_arrays(result["state"], prefix)); arrays.update(_tangent_arrays(result["tangent"], prefix+"d_"))
            arrays[prefix+"source_phi0"] = changed.source_phi0; arrays[prefix+"source_phi1"] = changed.source_phi1
            rows.append({"index": index, "width": branch["width"], "sign": branch["sign"], "outer_h": branch["outer_h"],
                         "status": result["status"], "stop": result["stop"], "time": result["time"],
                         "tau": result["tau"], "delta_tau": result["delta_tau"], "target_tau": target,
                         "readout": result["readout"], "first_retarded_coefficient": result["readout"]["delta_child_regional_content_tau"],
                         "event": None if result["event"] is None else {k:v for k,v in result["event"].items() if k not in ("state", "tangent")},
                         "proper_geometry": proper_geometry(result["pair"], result["state"], branch["width"], result["tangent"]),
                         "gaussian_state": gaussian_report(result["pair"], result["state"]),
                         "source_reset": False, "tangent_carried_to_common_clock_event": True})
            if result["status"] != "event_reached":
                status = result["status"]; break
    coefficients = []; delta = NONLINEAR_HELD_OUT_WIDTH-BASE_WIDTH
    if status == "predicted":
        for h in DERIVATIVE_H:
            selected = sorted([row for row in rows if row["outer_h"] == h], key=lambda row: row["sign"])
            second = curvature_from_slopes([row["first_retarded_coefficient"] for row in selected], h)
            base_value = nominal["readout"]["child_regional_content"]
            linear = base_value+delta*nominal["first_coefficient"]
            coefficients.append({"outer_h": h, "second_proper_clock_coefficient": second,
                                 "linear_forecast": linear, "quadratic_correction": .5*delta**2*second,
                                 "quadratic_forecast": linear+.5*delta**2*second})
    record = {"schema": NONLINEAR_SCHEMA, "inputs": config, "pair": prepared["pair"], "status": status,
              "producer_identity": prepared["producer_identity"], "prepare_curvature_json_sha256": episode.file_sha256(directory/(CURVATURE_STAGES[0]+".json")),
              "nominal": nominal, "target_tau": target, "neighbors": rows, "curvature_coefficients": coefficients,
              "curvature_h_difference": None if len(coefficients) != 2 else abs(coefficients[0]["second_proper_clock_coefficient"]-coefficients[1]["second_proper_clock_coefficient"]),
              "quadratic_forecast_h_indicator": None if len(coefficients) != 2 else abs(coefficients[0]["quadratic_forecast"]-coefficients[1]["quadratic_forecast"]),
              "forecast_locked_before_new_measurement": True, "new_held_out_source_constructed": False,
              "old_width1p05_used_for_fit": False, "hessian_integrated": False,
              "aggregate_cpu_seconds": _spent(prepared)+float(time.process_time()-start),
              "space_time_evolution_error_bound": None, "all_green_threshold": None,
              "universal_scale_power": None}
    if record["aggregate_cpu_seconds"] >= config["cpu_budget_seconds"]:
        record["status"] = "budget_stop"
    return seal(directory, CURVATURE_STAGES[1], record, arrays)


def measure_curvature(output):
    directory = _nonlinear_output(output); _new_stage(directory, CURVATURE_STAGES[2]); start = time.process_time()
    forecast, saved = load(directory, CURVATURE_STAGES[1])
    if forecast["status"] != "predicted":
        raise ValueError("quadratic forecast incomplete; new held-out source not launched")
    family.authenticate_producer({"producer_identity": forecast["producer_identity"]})
    config = forecast["inputs"]; pair = episode.pair_from_arrays(saved, forecast["pair"])
    changed = pair_at_width(pair, NONLINEAR_HELD_OUT_WIDTH)
    try:
        if _budget(config, _spent(forecast), start)[0]:
            raise TimeoutError("aggregate allowance exhausted before new held-out preparation")
        state, solver = family.solve_prepared(changed)
    except (ValueError, RuntimeError, TimeoutError) as error:
        record = {"schema": NONLINEAR_SCHEMA, "inputs": config, "producer_identity": forecast["producer_identity"],
                  "prediction_curvature_json_sha256": episode.file_sha256(directory/(CURVATURE_STAGES[1]+".json")),
                  "status": "budget_stop" if isinstance(error, TimeoutError) else "initial_preparation_failed",
                  "error": type(error).__name__+": "+str(error), "non_existence_claimed": False,
                  "aggregate_cpu_seconds": _spent(forecast)+float(time.process_time()-start)}
        return seal(directory, CURVATURE_STAGES[2], record, {"held_source_columns": changed.source_columns})
    initial = state.copy(); initial_geometry = proper_geometry(changed, state, NONLINEAR_HELD_OUT_WIDTH)
    result = march_retarded(changed, state, response.zero_tangent(changed.grid), config,
                            target_tau=forecast["target_tau"], prior_cpu=_spent(forecast), stage_start=start,
                            with_tangent=False)
    observed = result["readout"]["child_regional_content"]
    nominal = forecast["nominal"]["readout"]["child_regional_content"]
    estimates = forecast["curvature_coefficients"][-1]
    measured = result["status"] == "event_reached"
    event_indicator = None
    if measured:
        event_indicator = abs(result["readout"]["child_regional_content_dot"]*result["event"]["tau_residual"]/result["readout"]["tau_dot"])
    nominal_indicator = forecast["nominal"].get("first_coefficient_h_indicator")
    indicator = forecast["quadratic_forecast_h_indicator"]+(0. if nominal_indicator is None else
                abs(NONLINEAR_HELD_OUT_WIDTH-BASE_WIDTH)*nominal_indicator)+(event_indicator or 0.)
    effect = observed-nominal if measured else None
    remainder = observed-estimates["quadratic_forecast"] if measured else None
    arrays = _pair_arrays(result["pair"]); arrays.update(_state_arrays(initial, "initial_")); arrays.update(_state_arrays(result["state"], "held_"))
    record = {"schema": NONLINEAR_SCHEMA, "inputs": config, "pair": _pair_record(changed),
              "producer_identity": forecast["producer_identity"],
              "prediction_curvature_json_sha256": episode.file_sha256(directory/(CURVATURE_STAGES[1]+".json")),
              "status": "measured_at_common_absolute_tau" if measured else result["status"],
              "stop": result["stop"], "target_tau": forecast["target_tau"], "proper_time": result["tau"],
              "coordinate_time": result["time"], "delta_tau": result["delta_tau"],
              "event": None if result["event"] is None else {k:v for k,v in result["event"].items() if k not in ("state", "tangent")},
              "observed_child_fixed_window_content": observed,
              "linear_forecast": estimates["linear_forecast"], "quadratic_forecast": estimates["quadratic_forecast"],
              "observed_effect": effect, "linear_signed_remainder": None if not measured else observed-estimates["linear_forecast"],
              "quadratic_signed_remainder": remainder,
              "quadratic_remainder_fraction_of_effect": None if effect is None or effect == 0 else abs(remainder/effect),
              "available_h_and_event_indicator": indicator, "event_content_indicator": event_indicator,
              "effect_exceeds_available_indicators": None if effect is None else abs(effect)>indicator,
              "space_time_evolution_error_bound": None, "all_green_threshold": None,
              "old_width1p05_used_for_fit": False, "held_out_width_predeclared": NONLINEAR_HELD_OUT_WIDTH,
              "own_initial_solver": solver, "initial_radius_copied": False, "source_field_reset": False,
              "initial_geometry": initial_geometry, "final_geometry": proper_geometry(result["pair"], result["state"], NONLINEAR_HELD_OUT_WIDTH),
              "gaussian_state": gaussian_report(result["pair"], result["state"]),
              "W_unchanged": np.array_equal(result["pair"].geometry_map, saved["W"]),
              "observer_unchanged": np.array_equal(result["pair"].reference_columns, saved["observer_columns"]),
              "aggregate_cpu_seconds": _spent(forecast)+float(time.process_time()-start),
              "universal_scale_power": None, "physical_conclusion_forced": False}
    if record["aggregate_cpu_seconds"] >= config["cpu_budget_seconds"]:
        record["status"] = "budget_stop"
    return seal(directory, CURVATURE_STAGES[2], record, arrays)


def check_curvature(output):
    directory = _nonlinear_output(output); checked = []; previous = None; config = None; reference = None
    for index, stage in enumerate(CURVATURE_STAGES):
        if not (directory/(stage+".json")).exists():
            if (directory/(stage+".npz")).exists():
                raise ValueError("uncommitted nonlinear payload")
            continue
        if index and not (directory/(CURVATURE_STAGES[index-1]+".json")).exists():
            raise ValueError("nonlinear predecessor missing")
        record, arrays = load(directory, stage)
        if record["schema"] != NONLINEAR_SCHEMA:
            raise ValueError("nonlinear stage has another experiment schema")
        family.authenticate_producer({"producer_identity": record["producer_identity"]})
        if config is not None and record["inputs"] != config:
            raise ValueError("nonlinear numerical inputs changed")
        config = record["inputs"]
        if config["held_out_width"] != NONLINEAR_HELD_OUT_WIDTH or config["old_width1p05_used_for_fit"] is not False:
            raise ValueError("nonlinear held-out contract changed")
        if record.get("aggregate_cpu_seconds", 0.) > config["cpu_budget_seconds"] and record["status"] not in ("budget_stop", "preparation_incomplete", "initial_preparation_failed"):
            raise ValueError("completed nonlinear stage exceeded the aggregate allowance")
        if config["production"] and config["common_absolute_tau_target"] != ABSOLUTE_TAU_TARGET:
            raise ValueError("nonlinear proper-clock target changed")
        if index:
            key = "prepare_curvature_json_sha256" if index == 1 else "prediction_curvature_json_sha256"
            if record[key] != episode.file_sha256(directory/(CURVATURE_STAGES[index-1]+".json")):
                raise ValueError("nonlinear predecessor identity changed")
        if "W" in arrays:
            if reference is None:
                reference = arrays
            for name in ("W", "observer_columns", "source_weights"):
                if not np.array_equal(reference[name], arrays[name]):
                    raise ValueError("nonlinear stage changed fixed "+name)
            if not np.array_equal(arrays["W"], family.load_frozen_basis(config["nf"])["W"]):
                raise ValueError("nonlinear W differs from the frozen domain")
        if index == 0 and record["nominal"]["cached"]:
            path = Path(record["nominal"]["baseline_directory"])
            _prepared, old, _cached = cached_baseline(path, production=config["production"], expected_nf=config["nf"])
            if episode.file_sha256(path/"prediction.json") != record["nominal"]["prediction_json_sha256"] or old["payload_sha256"] != record["nominal"]["prediction_npz_sha256"]:
                raise ValueError("cached historical nominal identity changed")
        if index == 1 and record["status"] == "predicted":
            if len(record["neighbors"]) != 4 or len(record["curvature_coefficients"]) != 2:
                raise ValueError("nonlinear slope stencil incomplete")
            for row in record["neighbors"]:
                if row["status"] != "event_reached" or row["target_tau"] != record["target_tau"]:
                    raise ValueError("neighbor coefficient is not at the common clock event")
                if not row["event"]["tangent_carried_on_same_rooted_step"]:
                    raise ValueError("neighbor event dropped its tangent")
        checked.append({"stage": stage, "status": record["status"]}); previous = record
    if not checked:
        raise FileNotFoundError("no nonlinear width stages")
    return {"schema": NONLINEAR_SCHEMA, "ok": True, "checked": checked, "repaired": False, "evolved": False}


def _confirmation_output(output):
    directory = _output(output)
    if directory.is_relative_to(episode.REPO) and not any(part.startswith("nsc-discovery-width-response-confirmation-") for part in directory.parts):
        raise PermissionError("numerical confirmation needs a new nsc-discovery-width-response-confirmation- output")
    return directory


def confirmation_cases():
    return [{"nf": nf, "step_cap": cap, "width": w}
            for nf, cap in ((256, .00025), (128, .0005))
            for w in (BASE_WIDTH, NONLINEAR_HELD_OUT_WIDTH)]


def lock_curvature_confirmation(output, reference, *, production=False, cpu_budget=CONFIRMATION_CPU_CAP):
    directory = _confirmation_output(output)
    if directory == Path(reference).resolve() or Path(reference).resolve() in directory.parents:
        raise PermissionError("confirmation output cannot extend the reference evidence")
    _new_stage(directory, CONFIRMATION_STAGES[0])
    _new_stage(directory, CONFIRMATION_STAGES[1])
    start = time.process_time(); forecast, saved = load(reference, CURVATURE_STAGES[1])
    family.authenticate_producer({"producer_identity": forecast["producer_identity"]})
    if forecast["status"] != "predicted":
        raise ValueError("reference quadratic forecast is incomplete")
    if forecast["inputs"]["production"] and not production:
        raise PermissionError("numerical confirmation is a root production batch")
    if not np.isfinite(cpu_budget) or not 0 < cpu_budget <= CONFIRMATION_CPU_CAP:
        raise ValueError("confirmation aggregate allowance must lie in (0,120]")
    config = inputs(256, forecast["inputs"]["duration"], .00025, production, cpu_budget)
    config.update(experiment="single_numerical_confirmation", held_out_width=NONLINEAR_HELD_OUT_WIDTH,
                  cases=confirmation_cases(), common_absolute_tau_target=forecast["target_tau"],
                  no_new_response_coefficients=True, original_forecast_recalibrated=False,
                  maximum_confirmation_batches=1)
    config["float_hex"]["held_out_width"] = NONLINEAR_HELD_OUT_WIDTH.hex()
    config["inputs_sha256"] = _digest({k:v for k,v in config.items() if k != "inputs_sha256"})
    if config["production"] and forecast["target_tau"] != ABSOLUTE_TAU_TARGET:
        raise ValueError("reference proper-clock target differs from the original fine target")
    original = forecast["curvature_coefficients"][-1]
    locked = {"nominal_content": forecast["nominal"]["readout"]["child_regional_content"],
              "first_coefficient": forecast["nominal"]["first_coefficient"],
              "second_coefficient": original["second_proper_clock_coefficient"],
              "linear_forecast": original["linear_forecast"], "quadratic_forecast": original["quadratic_forecast"],
              "curvature_h_indicator": forecast["quadratic_forecast_h_indicator"],
              "all_green_threshold": None}
    record = {"schema": CONFIRMATION_SCHEMA, "inputs": config, "status": "forecast_locked",
              "reference_directory": str(Path(reference).resolve()),
              "reference_prediction_json_sha256": episode.file_sha256(Path(reference)/(CURVATURE_STAGES[1]+".json")),
              "reference_prediction_npz_sha256": forecast["payload_sha256"],
              "reference_producer_identity": forecast["producer_identity"],
              "locked_original_forecast": locked, "locked_forecast_sha256": _digest(locked),
              "producer_identity": _pin_identity(), "own_control_sources_constructed": False,
              "aggregate_cpu_seconds": float(time.process_time()-start)}
    arrays = {name: saved[name] for name in ("W", "observer_columns", "source_weights")}
    record["producing_commit"] = record["producer_identity"]["working_tree_commit"]
    record["producing_commit_authenticated"] = record["producer_identity"]["all_sources_commit_pinned"]
    return seal(directory, CONFIRMATION_STAGES[0], record, arrays)


def measure_confirmation(output):
    directory = _confirmation_output(output); _new_stage(directory, CONFIRMATION_STAGES[1]); start = time.process_time()
    locked, _saved = load(directory, CONFIRMATION_STAGES[0])
    family.authenticate_producer({"producer_identity": locked["producer_identity"]})
    if _digest(locked["locked_original_forecast"]) != locked["locked_forecast_sha256"]:
        raise ValueError("locked original coefficients changed")
    config = locked["inputs"]; rows = []; arrays = {}; failures = []
    for index, case in enumerate(config["cases"]):
        prefix = f"a{index}_"
        try:
            if _budget(config, _spent(locked), start)[0]:
                raise TimeoutError("confirmation aggregate allowance exhausted")
            pair, _pins = family.assemble_pair(case["nf"], 1., 1.)
            changed = pair_at_width(pair, case["width"])
            state, solver = family.solve_prepared(changed)
            initial = state.copy()
            arm_config = dict(config, nf=case["nf"], step_cap=case["step_cap"])
            result = march_retarded(changed, state, response.zero_tangent(changed.grid), arm_config,
                target_tau=config["common_absolute_tau_target"], prior_cpu=_spent(locked), stage_start=start,
                with_tangent=False)
            arrays.update({prefix+name: value for name,value in _pair_arrays(changed).items()})
            arrays.update(_state_arrays(initial, prefix+"initial_")); arrays.update(_state_arrays(result["state"], prefix+"event_"))
            content = result["readout"]["child_regional_content"]
            event = result["event"]
            clock_indicator = None if event is None else abs(result["readout"]["child_regional_content_dot"]*event["tau_residual"]/result["readout"]["tau_dot"])
            row = {"index": index, **case, "pair": _pair_record(changed), "status": result["status"],
                   "stop": result["stop"], "time": result["time"], "proper_time": result["tau"],
                   "target_tau": config["common_absolute_tau_target"], "absolute_child_content": content,
                   "absolute_discrepancy_from_original_forecast": content-(locked["locked_original_forecast"]["nominal_content"]
                        if case["width"] == BASE_WIDTH else locked["locked_original_forecast"]["quadratic_forecast"]),
                   "clock_event_indicator": clock_indicator,
                   "event": None if event is None else {k:v for k,v in event.items() if k not in ("state", "tangent")},
                   "own_initial_solver": solver, "initial_radius_copied": False, "source_field_reset": False,
                   "initial_geometry": proper_geometry(changed, initial, case["width"]),
                   "event_geometry": proper_geometry(result["pair"], result["state"], case["width"]),
                   "gaussian_state": gaussian_report(result["pair"], result["state"]),
                   "control_baseline_is_a_new_sealed_forecast": False}
            rows.append(row)
            if result["status"] != "event_reached":
                break
        except (ValueError, RuntimeError, TimeoutError) as error:
            failures.append({"index": index, **case, "error": type(error).__name__+": "+str(error),
                             "non_existence_claimed": False}); break
    groups = []
    original_effect = locked["locked_original_forecast"]["quadratic_forecast"]-locked["locked_original_forecast"]["nominal_content"]
    for nf in (256, 128):
        pair_rows = [row for row in rows if row["nf"] == nf and row["status"] == "event_reached"]
        if len(pair_rows) == 2:
            baseline, held = sorted(pair_rows, key=lambda row: row["width"])
            effect = held["absolute_child_content"]-baseline["absolute_child_content"]
            groups.append({"nf": nf, "step_cap": baseline["step_cap"], "measured_effect": effect,
                           "baseline_absolute_content": baseline["absolute_child_content"],
                           "held_absolute_content": held["absolute_child_content"],
                           "original_locked_predicted_effect": original_effect,
                           "effect_discrepancy_from_original_forecast": effect-original_effect,
                           "baseline_not_used_to_recalibrate_forecast": True})
    # Only after the original forecast lock and controls: read old observations as a numerical comparator.
    reference_measurement, _reference_arrays = load(locked["reference_directory"], CURVATURE_STAGES[2])
    family.authenticate_producer({"producer_identity": reference_measurement["producer_identity"]})
    if reference_measurement["prediction_curvature_json_sha256"] != locked["reference_prediction_json_sha256"]:
        raise ValueError("reference measurement used another forecast")
    reference_effect = reference_measurement.get("observed_effect")
    fine = next((row for row in groups if row["nf"] == 256), None)
    coarse = next((row for row in groups if row["nf"] == 128), None)
    indicators = {"time_effect_indicator": None if fine is None or reference_effect is None else abs(fine["measured_effect"]-reference_effect),
                  "space_effect_indicator": None if fine is None or coarse is None else abs(reference_effect-coarse["measured_effect"]) if reference_effect is not None else None,
                  "time_comparison": "NF256 cap0.00025 controls minus original NF256 cap0.0005 effect",
                  "space_comparison": "NF128 cap0.0005 controls minus original NF256 cap0.0005 effect",
                  "not_a_certified_error_bound": True, "additional_refinement_automatically_requested": False}
    original_gap = None if reference_effect is None else reference_effect-original_effect
    movements = (None if any(indicators[name] is None for name in ("time_effect_indicator", "space_effect_indicator"))
                 else indicators["time_effect_indicator"]+indicators["space_effect_indicator"])
    record = {"schema": CONFIRMATION_SCHEMA, "inputs": config, "status": "confirmation_completed" if len(groups) == 2 else "confirmation_incomplete",
              "confirmation_lock_json_sha256": episode.file_sha256(directory/(CONFIRMATION_STAGES[0]+".json")),
              "locked_forecast_sha256": locked["locked_forecast_sha256"],
              "locked_original_forecast": locked["locked_original_forecast"], "producer_identity": locked["producer_identity"],
              "arms": rows, "failures": failures, "resolution_effects": groups,
              "reference_measurement_json_sha256": episode.file_sha256(Path(locked["reference_directory"])/(CURVATURE_STAGES[2]+".json")),
              "reference_observed_effect": reference_effect, "numerical_indicators": indicators,
              "original_effect_forecast_gap": original_gap, "time_plus_space_movement_indicator": movements,
              "original_gap_exceeds_measured_numerical_movements": None if original_gap is None or movements is None else abs(original_gap)>movements,
              "gap_attribution": "compare raw gap with this one time and spatial indicator; no forced conclusion",
              "original_forecast_recalibrated": False, "new_response_coefficients_computed": False,
              "all_green_threshold": None, "aggregate_cpu_seconds": _spent(locked)+float(time.process_time()-start)}
    if (record["aggregate_cpu_seconds"] >= config["cpu_budget_seconds"]
        or any(row["status"] == "budget_stop" for row in rows)
        or any("TimeoutError" in item["error"] for item in failures)):
        record["status"] = "budget_stop"
    return seal(directory, CONFIRMATION_STAGES[1], record, arrays)


def confirm_curvature(output, reference, *, production=False, cpu_budget=CONFIRMATION_CPU_CAP):
    lock_curvature_confirmation(output, reference, production=production, cpu_budget=cpu_budget)
    return measure_confirmation(output)


def check_confirmation(output):
    directory = _confirmation_output(output); checked = []; locked = None
    for stage in CONFIRMATION_STAGES:
        if not (directory/(stage+".json")).exists():
            continue
        record, arrays = load(directory, stage)
        family.authenticate_producer({"producer_identity": record["producer_identity"]})
        if _digest(record["locked_original_forecast"]) != record["locked_forecast_sha256"]:
            raise ValueError("original locked forecast changed")
        if record.get("aggregate_cpu_seconds", 0.) > record["inputs"]["cpu_budget_seconds"] and record["status"] != "budget_stop":
            raise ValueError("completed confirmation exceeded its CPU allowance")
        if record["inputs"]["cases"] != confirmation_cases() or record["inputs"]["cpu_budget_seconds"] > CONFIRMATION_CPU_CAP:
            raise ValueError("numerical confirmation case or allowance contract changed")
        if stage == CONFIRMATION_STAGES[0]:
            locked = record
            reference, _source = load(record["reference_directory"], CURVATURE_STAGES[1])
            family.authenticate_producer({"producer_identity": reference["producer_identity"]})
            if reference["payload_sha256"] != record["reference_prediction_npz_sha256"] or episode.file_sha256(Path(record["reference_directory"])/(CURVATURE_STAGES[1]+".json")) != record["reference_prediction_json_sha256"]:
                raise ValueError("locked historical forecast identity changed")
        else:
            if locked is None or record["locked_original_forecast"] != locked["locked_original_forecast"]:
                raise ValueError("confirmation recalibrated original coefficients")
            if record["confirmation_lock_json_sha256"] != episode.file_sha256(directory/(CONFIRMATION_STAGES[0]+".json")):
                raise ValueError("confirmation lock predecessor changed")
            for arm in record["arms"]:
                prefix = f"a{arm['index']}_"
                if not np.array_equal(arrays[prefix+"W"], family.load_frozen_basis(arm["nf"])["W"]):
                    raise ValueError("control used another domain W")
                for component in ("phi0", "phi1"):
                    if not np.array_equal(arrays[prefix+"initial_"+component], arrays[prefix+"source_"+component]):
                        raise ValueError("control initial field differs from own source")
                if arm["target_tau"] != record["inputs"]["common_absolute_tau_target"]:
                    raise ValueError("control clock target changed")
        checked.append({"stage": stage, "status": record["status"]})
    if not checked:
        raise FileNotFoundError("no numerical confirmation stages")
    return {"schema": CONFIRMATION_SCHEMA, "ok": True, "checked": checked, "repaired": False, "evolved": False}
