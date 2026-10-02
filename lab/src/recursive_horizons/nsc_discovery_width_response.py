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


def frame_direction(pair, state, h):
    h = float(h)
    if h not in DERIVATIVE_H:
        raise ValueError("frame direction uses predeclared h=0.001 or 0.0005")
    plus0, plus1, plus_meta = source_columns(pair.grid, BASE_WIDTH+h)
    minus0, minus1, minus_meta = source_columns(pair.grid, BASE_WIDTH-h)
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
    report = {"h": h, "h_hex": h.hex(), "plus": plus_meta, "minus": minus_meta,
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
    record = dict(record, schema=SCHEMA, stage=stage, npz=npz_path.name,
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
    if payload_path != directory/(stage+".npz") or record["schema"] != SCHEMA or record["stage"] != stage:
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
