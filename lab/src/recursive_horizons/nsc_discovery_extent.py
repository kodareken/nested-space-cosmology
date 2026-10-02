"""Ambient-period intervention of the existing conformal Galerkin equations.

Only grid extent, complete frame and regional/checkpoint adapters are new.
The source is translated by two physical units, never dilated with period.
"""
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
import hashlib
import io
import json
from pathlib import Path
import time

import numpy as np

from . import nsc_discovery_backend as backend
from . import nsc_discovery_episode as episode
from . import nsc_discovery_family as family
from . import nsc_discovery_tidal as tidal
from . import nsc_discovery_translation as translation
from . import nsc_nested_parent_child as model
from . import nsc_regional_energy_exchange as regional
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

SCHEMA = "NSC-DISCOVERY-EXTENT-v1"
PERIODS = (8., 12.)
POINTS_PER_UNIT = 32
SHIFT = 2.
PARENT = (2., 6.)
CHILD = (3., 5.)
PROTECTED = (1., 7.)
CLOCKS = (3., 4., 5.)
CPU_CAP = 21600.
FORECAST_FACTOR = 1.5
CHUNK_LIMIT = 64*1024**2
CADENCE = 0.1
FIELDS = model.STATE_NAMES


def _json(value):
    return json.dumps(episode.jsonable(value), sort_keys=True, allow_nan=False, separators=(",", ":"))


def _hash(value):
    return hashlib.sha256(_json(value).encode()).hexdigest()


def settings(*, periods=PERIODS, points_per_unit=POINTS_PER_UNIT, prefix=0.3,
             stations=(8., 12.), step_cap=0.0005, cadence=CADENCE, cpu_budget=CPU_CAP,
             period16_admission=None):
    periods = tuple(float(x) for x in periods)
    if not periods or len(set(periods)) != len(periods) or any(x not in (8., 12., 16.) for x in periods):
        raise ValueError("periods must be distinct members of 8, 12 and admitted 16")
    if 16. in periods:
        if not period16_admission or period16_admission.get("extent_discriminates") is not True:
            raise PermissionError("period16 opens only after an explicit discriminating period8/12 assessment")
    density = int(points_per_unit)
    if density < 8 or density != points_per_unit:
        raise ValueError("integer fermion density of at least 8 is required")
    stations = tuple(float(x) for x in stations)
    if not stations or list(stations) != sorted(set(stations)) or any(not np.isfinite(x) or x <= prefix for x in stations):
        raise ValueError("stations must increase after the prefix")
    if not np.isfinite(prefix) or not 0 < prefix <= 0.3:
        raise ValueError("prefix must be positive and at most 0.3")
    if not np.isfinite(step_cap) or step_cap <= 0 or not 0 < cadence <= 0.1:
        raise ValueError("positive timestep and diagnostic cadence at most 0.1 are required")
    if not np.isfinite(cpu_budget) or not 0 < cpu_budget <= CPU_CAP:
        raise ValueError("aggregate CPU budget must lie in (0,21600]")
    result = {"periods": list(periods), "points_per_unit": density, "prefix": float(prefix),
              "stations": list(stations), "step_cap": float(step_cap), "cadence": float(cadence),
              "cpu_budget_seconds": float(cpu_budget), "forecast_factor": FORECAST_FACTOR,
              "translation": SHIFT, "parent": list(PARENT), "child": list(CHILD),
              "protected": list(PROTECTED), "clocks": list(CLOCKS),
              "source_dilated_with_period": False, "Q0_free": False,
              "carrier_k": float(coupling.CARRIER_K), "phase": float(coupling.CALIBRATION["phase"]),
              "locked_coefficients": coupling.locked_coefficients(), "kappa": int(coupling.KAPPA),
              "occupations": coupling.OCCUPATIONS.tolist(), "period16_admission": period16_admission,
              "nominal_dx_f": 1./density, "nominal_dx_q": 1./(4*density),
              "geometric_resolution_exactly_identical": False}
    result["inputs_sha256"] = _hash(result)
    return result


def _readonly(value):
    return family._readonly(value)


def build_pair(length, nf=None):
    """Native length-aware grid, source, observer and complete translated W."""
    length = float(length); nf = int(length*POINTS_PER_UNIT if nf is None else nf)
    if length not in (8., 12., 16.) or nf < 64 or nf % 2:
        raise ValueError("extent pair requires period8/12/16 and resolved even nf>=64")
    if abs(SHIFT*nf/length-round(SHIFT*nf/length)) > 1e-10:
        raise ValueError("physical translation must lie on the AP fermion lattice")
    grid = galerkin.build_grid(nf, quadrature=4*nf, length=length, gauge="conformal")
    phi0, phi1, metadata = model._separated_columns(grid)
    reference0, reference1, _ = coupling.prepare_rank6(grid.xi_f, grid.dx_f)
    W, coarse, child, parent, geometry_meta = model._geometry_frame(grid)
    reference = np.vstack((reference0, reference1)); source = np.vstack((phi0, phi1))
    raw = model.NestedPair(grid, _readonly(W), family._readonly(coarse, int), family._readonly(child, int),
                           family._readonly(parent, int), _readonly(reference0), _readonly(reference1),
                           _readonly(reference), _readonly(phi0), _readonly(phi1), _readonly(source),
                           _readonly(coupling.OCCUPATIONS), dict(metadata), dict(geometry_meta))
    dummy = model.encode_state(raw, galerkin.blank_state(grid, phi0, phi1))
    moved, _dummy = translation.translate_pair_state(raw, dummy, SHIFT)
    support = [[2., 3.], [3., 5.], [5., 6.]]
    columns = np.asarray(moved.source_columns)
    outside_g = (grid.xi_g < CHILD[0]) | (grid.xi_g > CHILD[1])
    outside_q = (grid.xi_q < CHILD[0]) | (grid.xi_q > CHILD[1])
    details = moved.geometry_map[:, child]
    fine_details = grid.A_g @ details
    geometry_meta = dict(moved.geometry_metadata,
        child_nodal_outside_power_fraction=(np.sum(details[outside_g]**2, axis=0)/np.sum(details**2, axis=0)).tolist(),
        child_detail_outside_power_fraction=(np.sum(fine_details[outside_q]**2, axis=0)/np.sum(fine_details**2, axis=0)).tolist(),
        period=length, basis="new complete frame of this physical domain", frozen_old_basis_reused=False)
    metadata = dict(moved.source_metadata, source_support_closures=support,
                    carrier_k=float(coupling.CARRIER_K), phase=float(coupling.CALIBRATION["phase"]),
                    local_carrier_origins=[2., 3., 5.], physical_widths=[1., 2., 1.],
                    source_dilated=False, gram_max=float(np.max(np.abs(columns.conj().T@columns-np.eye(6)))),
                    full_mean_current_retained=True, angular_inputs_changed=False, period=length)
    return replace(moved, source_metadata=metadata, geometry_metadata=geometry_meta)


def pair_record(pair):
    return {"length": float(pair.grid.length), "nf": int(pair.grid.nf), "nq": int(pair.grid.nq),
            "dx_g": float(pair.grid.dx_g), "dx_f": float(pair.grid.dx_f), "dx_q": float(pair.grid.dx_q),
            "coarse_indices": pair.geometry_coarse_indices.tolist(), "child_indices": pair.geometry_child_indices.tolist(),
            "parent_indices": pair.geometry_parent_indices.tolist(), "source_metadata": pair.source_metadata,
            "geometry_metadata": pair.geometry_metadata, "parent_interval": list(pair.parent_interval),
            "child_interval": list(pair.child_interval), "clock_locations": list(pair.clock_locations)}


def pair_from_arrays(arrays, record):
    """No nf-only or frozen period-eight reconstruction is permitted."""
    length = float(record["length"]); nf = int(record["nf"]); nq = int(record["nq"])
    if np.asarray(arrays["length"]).shape != (1,) or float(arrays["length"][0]) != length:
        raise ValueError("checkpoint physical period differs from its hashed length array")
    if int(arrays["quadrature"][0]) != nq or nq != 4*nf:
        raise ValueError("checkpoint quadrature domain differs")
    grid = galerkin.build_grid(nf, quadrature=nq, length=length, gauge="conformal")
    grid.fine = replace(grid.fine, occupations=np.array(arrays["source_weights"], copy=True))
    if arrays["W"].shape != (grid.ng, grid.ng):
        raise ValueError("checkpoint W does not span its own domain")
    reference = arrays["observer_columns"]; source = np.vstack((arrays["source_phi0"], arrays["source_phi1"]))
    return model.NestedPair(grid, _readonly(arrays["W"]), family._readonly(record["coarse_indices"], int),
        family._readonly(record["child_indices"], int), family._readonly(record["parent_indices"], int),
        _readonly(reference[:nf]), _readonly(reference[nf:]), _readonly(reference),
        _readonly(arrays["source_phi0"]), _readonly(arrays["source_phi1"]), _readonly(source),
        _readonly(arrays["source_weights"]), dict(record["source_metadata"]), dict(record["geometry_metadata"]),
        parent_interval=tuple(record["parent_interval"]), child_interval=tuple(record["child_interval"]),
        clock_locations=tuple(record["clock_locations"]))


def arrays_from(pair, state=None, clocks=None):
    arrays = {"length": np.array([pair.grid.length]), "quadrature": np.array([pair.grid.nq]),
              "W": pair.geometry_map, "source_phi0": pair.source_phi0, "source_phi1": pair.source_phi1,
              "observer_columns": pair.reference_columns, "source_weights": pair.weights}
    if state is not None:
        arrays.update({name: getattr(state, name) for name in FIELDS})
        arrays["normal_clocks"] = np.zeros(3) if clocks is None else np.asarray(clocks, dtype=float)
    return arrays


def state_from(arrays):
    return model.NestedState(*(np.array(arrays[name], copy=True) for name in FIELDS))


def source_pins(arrays):
    return {name: episode.array_sha256(arrays[name]) for name in
            ("length", "quadrature", "W", "source_phi0", "source_phi1", "observer_columns", "source_weights")}


def translation_audit(pair):
    """Independent integer AP roll and periodic force/covariance audit at L8."""
    if pair.grid.length != 8.:
        return {"applicable": False, "reason": "period8 reference audit"}
    old = model.build_pair(pair.grid.nf, quadrature=pair.grid.nq)
    cells = int(round(SHIFT/pair.grid.dx_f))
    expected0 = np.roll(old.source_phi0, cells, axis=0).copy(); expected0[:cells] *= -1
    expected1 = np.roll(old.source_phi1, cells, axis=0).copy(); expected1[:cells] *= -1
    expected = np.vstack((expected0, expected1)); actual = np.asarray(pair.source_columns)
    covariance_gap = float(np.max(np.abs((actual*pair.weights)@actual.conj().T
                                         -(expected*old.weights)@expected.conj().T)))
    old_state = model.encode_state(old, galerkin.blank_state(old.grid, old.source_phi0, old.source_phi1))
    new_state = model.encode_state(pair, galerkin.blank_state(pair.grid, pair.source_phi0, pair.source_phi1))
    _old_rate, old_bundle = model.rates(old, old_state, return_bundle=True)
    _new_rate, new_bundle = model.rates(pair, new_state, return_bundle=True)
    fine_cells = int(round(SHIFT/pair.grid.dx_q))
    forces = {name: float(np.max(np.abs(new_bundle["source"][name]-np.roll(old_bundle["source"][name], fine_cells))))
              for name in ("force_L", "force_Q", "force_beta")}
    scales = {name: max(1., float(np.max(np.abs(old_bundle["source"][name])))) for name in forces}
    passed = covariance_gap < 2e-9 and all(forces[name] <= 2e-9*scales[name] for name in forces)
    return {"applicable": True, "covariance_gap": covariance_gap, "forces_gap": forces,
            "forces_scale": scales, "passed": passed, "source_carrier_k": float(coupling.CARRIER_K),
            "old_geometry_reused_for_new_radius": False, "test_state": "constant-Q blank metric; source/operator audit"}


def _summary(value):
    value = np.asarray(value, dtype=float)
    if not np.isfinite(value).all():
        raise ValueError("nonfinite diagnostic values")
    scale = float(np.max(np.abs(value)))
    rms = 0. if scale == 0. else scale*float(np.sqrt(np.mean((value/scale)**2)))
    return {"min": float(np.min(value)), "max": float(np.max(value)), "max_abs": scale, "rms": rms}


def observe(pair, state, coordinate_time, clocks, mode="coupled", *, protected=False):
    """Underlying period-aware energy, tidal and integral owners, no old atlas."""
    jets = tidal.analytic_accelerations(pair, state, mode)
    fine, system, source = jets["fine"], jets["bundle"]["fine_system"], jets["bundle"]["source"]
    profile = tidal.profiles_on_grid(pair, jets)
    for name in ("R_h", "R4", "R_0101", "R_0202", "owned_W"):
        if not np.isfinite(profile["tides"][name]).all():
            raise ValueError("actual tidal diagnostic is nonfinite: " + name)
    ledger = regional.matter_ledger(system, fine)
    metric = model.metrics(pair, state)
    mass_columns = (np.abs(fine.phi0)**2+np.abs(fine.phi1)**2)*pair.weights[None, :]
    tagged = np.sum(mass_columns[:, [2, 3]], axis=1)
    density = np.sum(mass_columns, axis=1)/pair.grid.dx_q
    tagged_total = float(np.sum(tagged))
    tagged_child = model.interval_integral(pair.grid, tagged/pair.grid.dx_q, pair.child_interval)
    constraints = galerkin.constraint_diagnostics(pair.grid, jets["nodal"], source=source, fine_state=fine)
    columns = np.vstack((state.phi0, state.phi1))
    gram = columns.conj().T@columns
    car = galerkin.occupation_eigenvalues(state.phi0, state.phi1, pair.weights)
    window = lambda values, interval: model.interval_integral(pair.grid, np.asarray(values, dtype=float), interval)
    energies = {"total": float(coupling.gravity_energy(system, fine)+coupling.field_energy(system, fine)),
                "field": float(coupling.field_energy(system, fine)), "gravity": float(coupling.gravity_energy(system, fine))}
    worldline = int(np.argmin(np.abs(pair.grid.xi_q-pair.clock_locations[1])))
    tide_names = ("R_h", "R4", "R_0101", "R_0202", "owned_W")
    row = {"time": float(coordinate_time), "length": float(pair.grid.length), "mode": mode,
           "carrier_circuit_coordinate_time": float(pair.grid.length), "physical_principal_speed": 1.,
           "parent_interval": list(pair.parent_interval), "child_interval": list(pair.child_interval),
           "ambient_intervals": [[0., PARENT[0]], [PARENT[1], pair.grid.length]],
           "normal_clocks": np.asarray(clocks).tolist(), "clock_rates": metric["clock_rates"],
           "clock_locations": list(pair.clock_locations), "clock_quadrature": "endpoint_trapezoid_per_owned_RK4_step",
           "child_proper_length": metric["child_proper_length"], "parent_proper_length": metric["parent_proper_length"],
           "r": _summary(fine.r), "Q": _summary(fine.Q), "N": _summary(fine.r*fine.Q),
           "child_r_proper_mean": metric["child_r_proper_mean"],
           "child_Q_proper_mean": metric["child_Q_proper_mean"],
           "tagged_child_fraction": tagged_child/tagged_total,
           "tagged_source_columns": [2, 3], "tagged_total_probability": tagged_total,
           "total_probability": float(np.sum(mass_columns)),
           "all_source_child_probability": window(density, pair.child_interval),
           "gram_max": float(np.max(np.abs(gram-np.eye(6)))),
           "CAR_min": float(np.min(car)), "CAR_max": float(np.max(car)), "constraints": constraints,
           "source_rho": _summary(source["force_L"]/pair.grid.dx_q),
           "source_current": _summary(source["Pmom"]/pair.grid.dx_q),
           "shift_residual_mean": float(np.mean(coupling.shift_constraint(system, fine)+source["force_beta"]/pair.grid.dx_q)),
           "current_mean_deleted": False, "total_energy": energies,
           "child_matter_normal_energy": window(ledger["normal_energy_nodal"]/pair.grid.dx_q, pair.child_interval),
           "child_boundary_noether_flux": model._periodic_values(pair.grid, ledger["flux"]/pair.grid.dx_q, pair.child_interval).tolist(),
           "noether_flux": _summary(ledger["flux"]/pair.grid.dx_q),
           "actual_tides": {name: _summary(profile["tides"][name]) for name in tide_names},
           "worldline_actual_tides": {name: float(profile["tides"][name][worldline]) for name in tide_names},
           "actual_tides_finite": all(np.isfinite(profile["tides"][name]).all() for name in tide_names),
           "chi_substituted_for_curvature": False, "return_is_holding": False,
           "renewal_asserted": False, "continuum_certified": False,
           "frozen_mass_potential": "kappa times own handoff Q(x), independent of r directly"}
    samples = {}
    if protected:
        mask = (pair.grid.xi_q >= PROTECTED[0]) & (pair.grid.xi_q <= PROTECTED[1])
        dx = lambda value: tidal.metric.spectral_dx(value, pair.grid.length)
        values = {"r": fine.r, "Q": fine.Q, "N": fine.r*fine.Q, "q": fine.r*fine.Q,
                  "r_x": dx(fine.r), "Q_x": dx(fine.Q), "r_xx": dx(dx(fine.r)), "Q_xx": dx(dx(fine.Q)),
                  "rho": source["force_L"]/pair.grid.dx_q,
                  "hamilton_constraint": coupling.constraint_residuals(system, fine, source)["hamilton"],
                  "shift_constraint": coupling.constraint_residuals(system, fine, source)["momentum"],
                  **{name: profile["tides"][name] for name in tide_names}}
        samples = {"protected_x": pair.grid.xi_q[mask], **{"protected_"+name: np.asarray(value)[mask] for name, value in values.items()}}
        row["protected_interval"] = list(PROTECTED)
    return row, samples


def protected_matching(left, right):
    if not np.array_equal(left["protected_x"], right["protected_x"]):
        raise ValueError("protected comparison requires identical physical sample coordinates")
    gaps = {}
    for name in left:
        if name == "protected_x":
            continue
        a, b = left[name], right[name]
        absolute = float(np.max(np.abs(a-b))); scale = max(1., float(np.max(np.abs(a))))
        gaps[name.removeprefix("protected_")] = {"absolute_max": absolute, "relative_to_reference_scale": absolute/scale,
                                                "roundoff_equivalent": absolute <= 2e-9*scale}
    geometry_equal = all(gaps[name]["roundoff_equivalent"] for name in ("r", "N", "q", "Q", "r_x", "Q_x", "r_xx", "Q_xx"))
    return {"protected_interval": list(PROTECTED), "gaps": gaps,
            "exact_local_matching_guaranteed": False, "geometry_roundoff_equivalent": geometry_equal,
            "coupled_pure_causal_attribution": "unestablished" if geometry_equal else "ambiguous_due_to_changed_initial_local_data",
            "physical_smallness_bound": None, "numerical_equivalence_tolerance": 2e-9,
            "not_a_continuation_veto": True}


def _output(output):
    path = episode.assert_campaign_output(output)
    if path.is_relative_to(episode.REPO) and not any(part.startswith("nsc-discovery-extent-") for part in path.parts):
        raise PermissionError("extent evidence needs its own nsc-discovery-extent- successor directory")
    return path


def _write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    data = (_json(value)+"\n").encode()
    if len(data) > CHUNK_LIMIT:
        raise ValueError("JSON chunk exceeds 64 MiB")
    with path.open("xb") as handle:
        handle.write(data)
    path.chmod(0o444)


def commit(output, record, arrays):
    directory = _output(output); directory.mkdir(parents=True, exist_ok=True)
    case_id = record["case_id"]
    ordinal = episode._next_ordinal(directory, case_id)
    stem = f"{case_id}-{ordinal:06d}"
    buffer = io.BytesIO(); np.savez_compressed(buffer, **arrays); payload = buffer.getvalue()
    if any(not np.isfinite(value).all() for value in arrays.values()):
        raise ValueError("nonfinite checkpoint array")
    result = dict(record, schema=SCHEMA, ordinal=ordinal, npz=stem+".npz", json=stem+".json",
                  payload_sha256=hashlib.sha256(payload).hexdigest(), source_pins=source_pins(arrays),
                  array_sha256={name: episode.array_sha256(value) for name, value in arrays.items()},
                  immutable=True, momentum_representation=episode.CANONICAL_PI)
    if len(payload)+len(_json(result).encode()) > CHUNK_LIMIT:
        raise ValueError("checkpoint chunk exceeds 64 MiB")
    path = directory/result["npz"]
    with path.open("xb") as handle:
        handle.write(payload)
    path.chmod(0o444); _write(directory/result["json"], result)
    return result


def load(output, case_id, ordinal=None):
    directory = _output(output)
    ordinal = episode._next_ordinal(directory, case_id)-1 if ordinal is None else int(ordinal)
    path = directory/f"{case_id}-{ordinal:06d}.json"
    record = json.loads(path.read_text()); npz = directory/record["npz"]
    if record["schema"] != SCHEMA or record["case_id"] != case_id or record["ordinal"] != ordinal:
        raise ValueError("extent checkpoint identity mismatch")
    if episode.file_sha256(npz) != record["payload_sha256"]:
        raise ValueError("extent checkpoint payload mismatch")
    for item in (path, npz):
        if item.stat().st_mode & 0o222 or item.stat().st_size > CHUNK_LIMIT:
            raise ValueError("extent checkpoint is writable or oversized")
    with np.load(npz, allow_pickle=False) as saved:
        arrays = {name: np.array(saved[name], copy=True) for name in saved.files}
    if set(arrays) != set(record["array_sha256"]):
        raise ValueError("checkpoint array inventory mismatch")
    for name, value in arrays.items():
        if episode.array_sha256(value) != record["array_sha256"][name]:
            raise ValueError("checkpoint array mismatch: "+name)
    if source_pins(arrays) != record["source_pins"]:
        raise ValueError("checkpoint frozen source identity mismatch")
    return record, arrays


def _identity():
    identity = family.producer_identity()
    for path in (Path(__file__).resolve(), Path(translation.__file__).resolve(),
                 episode.LAB/"scripts"/"derive_nsc_discovery_extent.py", Path(regional.__file__).resolve()):
        relative = path.relative_to(episode.REPO).as_posix(); digest = episode.file_sha256(path)
        commit_hash = identity["working_tree_commit"]
        try:
            family.resolve_pinned_source_bytes(episode.REPO, relative, digest, commit=commit_hash)
        except (RuntimeError, ValueError):
            commit_hash = None
        if relative not in {item["path"] for item in identity["files"]}:
            identity["files"].append({"path": relative, "sha256": digest, "commit": commit_hash})
    identity["all_sources_commit_pinned"] = all(x["commit"] for x in identity["files"])
    return identity


def prepare(output, **kwargs):
    directory = _output(output)
    if (directory/"prepare.json").exists():
        raise FileExistsError("extent preparation already exists")
    config = settings(**kwargs); start = time.process_time(); cases = []; failures = []; protected = {}
    identity = _identity()
    for length in config["periods"]:
        case_id = f"L{int(length)}_initial"
        try:
            pair = build_pair(length, int(length*config["points_per_unit"]))
            audit = translation_audit(pair)
            if audit.get("applicable") and not audit["passed"]:
                raise ValueError("period8 translated source/operator audit failed")
            with backend.fft_thread_limit(1):
                state, solver = family.solve_prepared(pair)
                row, samples = observe(pair, state, 0., np.zeros(3), protected=True)
            arrays = arrays_from(pair, state); arrays.update(samples)
            record = dict(pair_record(pair), case_id=case_id, status="prepared", time=0., steps=0,
                          mode="coupled", solver=solver, observations=[row], translation_audit=audit,
                          input_sha256=config["inputs_sha256"], producer_identity=identity,
                          copied_radius=False, source_current_deleted=False)
            committed = commit(directory, record, arrays); cases.append(committed)
            protected[length] = samples
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
            failures.append({"length": length, "case_id": case_id, "status": "initial_preparation_failed",
                             "error": type(error).__name__+": "+str(error), "non_existence_claimed": False})
        if time.process_time()-start >= config["cpu_budget_seconds"]:
            break
    matching = None
    if 8. in protected and 12. in protected:
        matching = protected_matching(protected[8.], protected[12.])
    result = {"schema": SCHEMA, "stage": "prepare", "inputs": config, "cases": cases, "failures": failures,
              "initial_protected_matching": matching, "producer_identity": identity,
              "aggregate_cpu_seconds": float(time.process_time()-start), "executor_pools": 0,
              "initial_data_policy": "each extent solves its own radius; global periodic mean can change local data"}
    _write(directory/"prepare.json", result); return result


def advance(pair, state, *, current_time, target, step_cap, cadence, clocks, mode,
            cpu_allowance, forecast_factor=FORECAST_FACTOR, max_steps=None, production=False):
    if target > 0.01 and not production:
        raise PermissionError("extent trajectory above 0.01 belongs to the root production executor")
    pair, state, _info = episode.resolve_pair(pair, state, backend="fft")
    mark = float(current_time); clocks = np.array(clocks, copy=True); start = time.process_time()
    observations = []; steps = 0; cost = 0.; diag_cost = 0.; status = "target_reached"; stop = None
    stepper = episode.evolving_step if mode == "coupled" else episode.frozen_geometry_step
    rates = np.asarray(model.metrics(pair, state)["clock_rates"])
    next_diagnostic = min(target, episode.next_observation_time(mark, cadence))
    while mark < target-1e-12:
        if max_steps is not None and steps >= max_steps:
            status = "step_limit"; break
        with backend.fft_thread_limit(1):
            dt, restriction = episode.step_restriction(pair, state, step_cap)
        principal_dt = dt
        dt = min(dt, target-mark, next_diagnostic-mark)
        spent = time.process_time()-start
        remaining = int(np.ceil((target-mark)/max(principal_dt, 1e-15)))+int(np.ceil((target-mark)/cadence))
        forecast = spent+forecast_factor*(remaining*cost+int(np.ceil((target-mark)/cadence))*diag_cost)
        if spent >= cpu_allowance or forecast > cpu_allowance:
            status, stop = "budget_stop", {"forecast_cpu_seconds": forecast}; break
        before = time.process_time()
        try:
            with backend.fft_thread_limit(1):
                updated = stepper(pair, state, dt)
            if any(not np.isfinite(getattr(updated, name)).all() for name in FIELDS):
                raise FloatingPointError("nonfinite owned RK4 output")
            new_rates = np.asarray(model.metrics(pair, updated)["clock_rates"])
        except (coupling.PositiveChartExit, ValueError, FloatingPointError) as error:
            status, stop = "chart_or_numerical_stop", {"error": type(error).__name__+": "+str(error),
                "last_admissible_time": mark, "event_bracket": [mark, mark+dt], "physical_instability_claimed": False}; break
        clocks += 0.5*dt*(rates+new_rates); rates = new_rates; state = updated
        mark += dt; steps += 1; cost = max(cost, time.process_time()-before)
        if mark >= next_diagnostic-1e-12:
            before = time.process_time()
            try:
                with backend.fft_thread_limit(1):
                    row, _samples = observe(pair, state, mark, clocks, mode)
                observations.append(row)
            except (ValueError, FloatingPointError) as error:
                status, stop = "diagnostic_unavailable", {"error": type(error).__name__+": "+str(error)}; break
            diag_cost = max(diag_cost, time.process_time()-before)
            next_diagnostic = min(target, episode.next_observation_time(mark, cadence))
    return {"pair": pair, "state": state, "time": mark, "normal_clocks": clocks, "steps": steps,
            "status": status, "stop": stop, "observations": observations,
            "cpu_seconds": float(time.process_time()-start), "last_step_cpu": cost,
            "last_observation_cpu": diag_cost, "state_clamped": False,
            "clock_quadrature": "endpoint trapezoid per owned RK4 step"}


def materialize(output, *, production=False):
    directory = _output(output)
    if (directory/"prefix.json").exists():
        raise FileExistsError("extent prefix already exists")
    prepared = json.loads((directory/"prepare.json").read_text()); config = prepared["inputs"]
    family.authenticate_producer({"producer_identity": prepared["producer_identity"]})
    if config["prefix"] > 0.01 and not production:
        raise PermissionError("T=0.3 prefix requires root production execution")
    start = time.process_time(); cases = []; failures = []; protected = {}
    prior = prepared["aggregate_cpu_seconds"]
    for case in prepared["cases"]:
        record, arrays = load(directory, case["case_id"], 0)
        pair = pair_from_arrays(arrays, record); state = state_from(arrays)
        allowance = config["cpu_budget_seconds"]-prior-(time.process_time()-start)
        result = advance(pair, state, current_time=0., target=config["prefix"],
                         step_cap=config["step_cap"], cadence=config["cadence"], clocks=np.zeros(3),
                         mode="coupled", cpu_allowance=allowance, production=production)
        if result["status"] != "target_reached":
            failed = dict(record, case_id=f"L{int(record['length'])}_prefix_stop", status=result["status"],
                          time=result["time"], stop=result["stop"], observations=result["observations"])
            failures.append(commit(directory, failed, arrays_from(pair, result["state"], result["normal_clocks"]))); continue
        with backend.fft_thread_limit(1):
            row, samples = observe(result["pair"], result["state"], result["time"], result["normal_clocks"], protected=True)
        protected[record["length"]] = samples
        for mode in ("coupled", "frozen_geometry"):
            if mode == "coupled":
                branch_row = row
            else:
                with backend.fft_thread_limit(1):
                    branch_row, _unused = observe(result["pair"], result["state"], result["time"], result["normal_clocks"], mode)
            branch = dict(record, case_id=f"L{int(record['length'])}_{mode}", parent_case=record["case_id"],
                          status="handoff", mode=mode, time=result["time"], steps=0,
                          observations=[branch_row], prefix_steps=result["steps"],
                          prefix_state_sha256=episode.state_sha256(result["state"]),
                          initial_state_called_during_continuation=False, own_source_handoff=True, cumulative_child_cpu_seconds=0.)
            arrays = arrays_from(result["pair"], result["state"], result["normal_clocks"]); arrays.update(samples)
            cases.append(commit(directory, branch, arrays))
    matching = protected_matching(protected[8.], protected[12.]) if 8. in protected and 12. in protected else None
    result = {"schema": SCHEMA, "stage": "prefix", "inputs": config, "cases": cases, "failures": failures,
              "handoff_protected_matching": matching, "prepare_json_sha256": episode.file_sha256(directory/"prepare.json"),
              "producer_identity": prepared["producer_identity"], "executor_pools": 0,
              "aggregate_cpu_seconds": prior+float(time.process_time()-start)}
    _write(directory/"prefix.json", result); return result


def execute_case(payload):
    directory = payload["directory"]; case_id = payload["case_id"]
    start = time.process_time(); record, arrays = load(directory, case_id)
    pair = pair_from_arrays(arrays, record); state = state_from(arrays)
    prior_cpu = float(record.get("cumulative_child_cpu_seconds", 0.))
    config = payload["inputs"]; snapshots = []
    # Each station commits its own immutable observations and full last state.
    for target in config["stations"]:
        if target <= record["time"]+1e-12:
            continue
        remaining = payload["cpu_allowance"]-(time.process_time()-start)
        result = advance(pair, state, current_time=record["time"], target=target,
                         step_cap=config["step_cap"], cadence=config["cadence"],
                         clocks=arrays["normal_clocks"], mode=record["mode"], cpu_allowance=remaining,
                         max_steps=payload.get("max_steps"), production=payload["production"])
        new_arrays = arrays_from(result["pair"], result["state"], result["normal_clocks"])
        if source_pins(new_arrays) != record["source_pins"]:
            raise ValueError("continuation changed source, observer, W or period pins")
        updated = dict(record, time=result["time"], steps=record["steps"]+result["steps"], status=result["status"],
                       stop=result["stop"], observations=result["observations"], last_step_cpu=result["last_step_cpu"],
                       last_observation_cpu=result["last_observation_cpu"], station=target,
                       own_source_handoff=True, initial_state_called_during_continuation=False,
                       cumulative_child_cpu_seconds=prior_cpu+float(time.process_time()-start))
        committed = commit(directory, updated, new_arrays); snapshots.append(committed)
        record, arrays, pair, state = committed, new_arrays, result["pair"], result["state"]
        if result["status"] != "target_reached":
            break
    return {"case_id": case_id, "status": record["status"], "time": record["time"],
            "cpu_seconds": float(time.process_time()-start), "snapshots": snapshots}


def run(output, *, production=False, workers=6, max_steps=None, executor=None):
    directory = _output(output)
    if (directory/"run.json").exists():
        raise FileExistsError("extent run already exists; use a new successor or a separately reviewed continuation")
    prefix = json.loads((directory/"prefix.json").read_text()); config = prefix["inputs"]
    if max(config["stations"]) > 0.01 and not production:
        raise PermissionError("extent stations 8/12 are root production execution")
    family.authenticate_producer({"producer_identity": prefix["producer_identity"]})
    cases = prefix["cases"]
    if not cases:
        raise ValueError("no admissible own-source handoffs")
    prior_child_cpu = sum(float(load(directory, case["case_id"])[0].get("cumulative_child_cpu_seconds", 0.)) for case in cases)
    allowance = config["cpu_budget_seconds"]-prefix["aggregate_cpu_seconds"]-prior_child_cpu
    if allowance <= 0:
        raise ValueError("aggregate CPU budget exhausted before executor")
    count = min(int(workers), 6, len(cases))
    if count < 1:
        raise ValueError("workers must be positive")
    start = time.process_time(); results = []
    child_allowance = max(0., allowance-min(10., allowance*0.05))
    factory = ProcessPoolExecutor if executor is None else executor
    with factory(max_workers=count) as pool:
        pending = [pool.submit(execute_case, {"directory": str(directory), "case_id": case["case_id"],
                   "inputs": config, "cpu_allowance": child_allowance/len(cases), "production": production,
                   "max_steps": max_steps}) for case in cases]
        for future in as_completed(pending):
            results.append(future.result())
    result = {"schema": SCHEMA, "stage": "run", "inputs": config, "results": results,
              "prefix_json_sha256": episode.file_sha256(directory/"prefix.json"), "executor_pools": 1,
              "workers": count, "producer_identity": prefix["producer_identity"],
              "aggregate_cpu_seconds": prefix["aggregate_cpu_seconds"]+prior_child_cpu+sum(x["cpu_seconds"] for x in results)+(time.process_time()-start),
              "initial_matching": json.loads((directory/"prepare.json").read_text())["initial_protected_matching"],
              "handoff_matching": prefix["handoff_protected_matching"],
              "return_alone_is_holding": False, "echo_proof": False, "continuum_certified": False,
              "period16_automatically_authorized": False}
    result["assessment"] = assess(directory, results, result["initial_matching"], result["handoff_matching"])
    result["aggregate_cpu_seconds"] = (prefix["aggregate_cpu_seconds"]+prior_child_cpu+sum(x["cpu_seconds"] for x in results)
                                       +time.process_time()-start)
    _write(directory/"run.json", result); return result


def assess(directory, results, initial_matching=None, handoff_matching=None):
    """Measured timing/geometry comparisons, never an echo or holding proof."""
    profiles = {}; candidates = {}
    for result in results:
        identifier = result["case_id"]
        rows = []
        for ordinal in range(episode._next_ordinal(directory, identifier)):
            record, _arrays = load(directory, identifier, ordinal)
            rows.extend(record.get("observations", []))
        rows = list({float(row["time"]): row for row in rows}.values())
        rows.sort(key=lambda row: row["time"])
        profiles[identifier] = rows
        # Four units is the parent crossing time. Exclude the initial retained packet.
        late = [row for row in rows if row["time"] >= PARENT[1]-PARENT[0]]
        peaks = []
        for index in range(1, len(late)-1):
            before, row, after = late[index-1:index+2]
            if row["tagged_child_fraction"] > before["tagged_child_fraction"] and row["tagged_child_fraction"] >= after["tagged_child_fraction"]:
                peaks.append({"time": row["time"], "fraction": row["tagged_child_fraction"],
                              "bracket": [before["time"], after["time"]],
                              "child_proper_length": row["child_proper_length"], "actual_tides": row["actual_tides"]})
        maximum = max(late, key=lambda row: row["tagged_child_fraction"]) if late else None
        candidates[identifier] = {
            "sampled_post_parent_crossing_peaks": peaks,
            "maximum_late_sample": None if maximum is None else {
                "time": maximum["time"], "fraction": maximum["tagged_child_fraction"],
                "at_recorded_boundary": maximum is late[0] or maximum is late[-1],
                "child_proper_length": maximum["child_proper_length"], "actual_tides": maximum["actual_tides"]},
            "holding_inferred": False,
            "circuit_time_hypothesis": None if not rows else rows[0]["length"],
            "diagnostic_spacing_is_timing_indicator_not_continuum_bound": True}
    comparisons = []
    for mode in ("coupled", "frozen_geometry"):
        left, right = profiles.get("L8_"+mode, []), profiles.get("L12_"+mode, [])
        for target in (8., 12.):
            a = next((row for row in left if abs(row["time"]-target) < 1e-9), None)
            b = next((row for row in right if abs(row["time"]-target) < 1e-9), None)
            if a and b:
                comparisons.append({"time": target, "mode": mode,
                    "L12_minus_L8_tagged_fraction": b["tagged_child_fraction"]-a["tagged_child_fraction"],
                    "L12_minus_L8_child_proper_length": b["child_proper_length"]-a["child_proper_length"],
                    "actual_tides_L8": a["actual_tides"], "actual_tides_L12": b["actual_tides"]})
    feedback = []
    for length in (8, 12, 16):
        coupled, frozen = profiles.get(f"L{length}_coupled", []), profiles.get(f"L{length}_frozen_geometry", [])
        for a in coupled:
            b = next((row for row in frozen if abs(row["time"]-a["time"]) < 1e-9), None)
            if b:
                feedback.append({"length": length, "time": a["time"],
                    "coupled_minus_frozen_tagged_fraction": a["tagged_child_fraction"]-b["tagged_child_fraction"],
                    "coupled_minus_frozen_child_proper_length": a["child_proper_length"]-b["child_proper_length"]})
    return {"return_candidates": candidates, "same_time_extent_comparisons": comparisons,
            "same_extent_feedback_comparisons": feedback, "initial_matching": initial_matching,
            "handoff_matching": handoff_matching, "pure_causal_echo_proof": False,
            "extent_discriminates": None, "period16_requires_root_review": True,
            "physical_return_or_holding_conclusion_forced": False}


def audit_existing_period8(directory, case_id, candidate_pair, candidate_state, ordinal=None, *, candidate_time=None, candidate_mode=None):
    """Read-only comparator admission; never silently bind old geometry as new."""
    record, arrays = episode.load_checkpoint(directory, case_id, ordinal)
    if int(record["nf"]) != candidate_pair.grid.nf or candidate_pair.grid.length != 8.:
        return {"admissible": False, "reason": "different period8 resolution"}
    legacy = dict(record, length=8., nq=4*record["nf"], parent_interval=[0., 4.], child_interval=[1., 3.],
                  clock_locations=[1., 2., 3.])
    arrays = dict(arrays, length=np.array([8.]), quadrature=np.array([4*record["nf"]]))
    old = pair_from_arrays(arrays, legacy); state = episode.state_from_arrays(arrays, episode.CANONICAL_PI)
    moved, state = translation.translate_pair_state(old, state, SHIFT)
    before, _ = model._active(moved, state); after, _ = model._active(candidate_pair, candidate_state)
    gaps = {name: float(np.max(np.abs(getattr(before, name)-getattr(after, name)))) for name in FIELDS}
    scales = {name: max(1., float(np.max(np.abs(getattr(before, name))))) for name in FIELDS}
    legacy_mode = episode.normalize_control_mode(record.get("control_mode") or record.get("geometry"))
    candidate_mode = legacy_mode if candidate_mode is None else episode.normalize_control_mode(candidate_mode)
    _old_nodal, old_nodal_rate, _old_bundle, _ = episode.actual_control_rate(moved, state, legacy_mode)
    _new_nodal, candidate_nodal_rate, _new_bundle, _ = episode.actual_control_rate(candidate_pair, candidate_state, candidate_mode)
    old_fine_rate = galerkin.prolong_state(moved.grid, old_nodal_rate)
    candidate_fine_rate = galerkin.prolong_state(candidate_pair.grid, candidate_nodal_rate)
    rate_gaps = {name: float(np.max(np.abs(getattr(old_fine_rate, name)-getattr(candidate_fine_rate, name)))) for name in FIELDS}
    rate_scales = {name: max(1., float(np.max(np.abs(getattr(old_fine_rate, name))))) for name in FIELDS}
    old_means, candidate_means = model.metrics(moved, state), model.metrics(candidate_pair, candidate_state)
    means = {name: {"legacy": old_means[name], "candidate": candidate_means[name],
                   "absolute_gap": abs(float(old_means[name])-float(candidate_means[name]))}
             for name in ("child_proper_length", "parent_proper_length", "child_r_proper_mean", "child_Q_proper_mean")}
    covariance_gap = float(np.max(np.abs((moved.source_columns*moved.weights)@moved.source_columns.conj().T
                                         -(candidate_pair.source_columns*candidate_pair.weights)@candidate_pair.source_columns.conj().T)))
    numerical_match = (candidate_mode == legacy_mode and covariance_gap < 2e-9 and all(gaps[name] < 2e-9*scales[name] for name in gaps)
                       and all(rate_gaps[name] < 2e-9*rate_scales[name] for name in rate_gaps))
    time_matches = candidate_time is not None and abs(float(candidate_time)-float(record["coordinate_time"])) < 1e-12
    return {"admissible": numerical_match and time_matches, "numerical_source_state_rate_match": numerical_match,
            "candidate_time_matches": time_matches,
            "covariance_gap": covariance_gap, "physical_state_gaps": gaps,
            "physical_rate_gaps": rate_gaps, "physical_rate_scales": rate_scales,
            "regional_geometric_means": means, "legacy_mode": legacy_mode, "candidate_mode": candidate_mode,
            "legacy_payload_sha256": record["arrays_sha256"], "legacy_time": record["coordinate_time"],
            "candidate_time_requires_external_check": candidate_time is None, "old_record_rewritten": False}


def check(output):
    directory = _output(output); stages = []; config = None; cases = []
    for name in ("prepare", "prefix", "run"):
        path = directory/(name+".json")
        if not path.exists():
            continue
        if path.stat().st_mode & 0o222 or path.stat().st_size > CHUNK_LIMIT:
            raise ValueError("extent stage manifest is writable or oversized")
        record = json.loads(path.read_text())
        if record["schema"] != SCHEMA:
            raise ValueError("extent manifest schema mismatch")
        current = record["inputs"]
        if _hash({k:v for k,v in current.items() if k != "inputs_sha256"}) != current["inputs_sha256"]:
            raise ValueError("extent input digest mismatch")
        if config is not None and config != current:
            raise ValueError("extent stage changed numerical inputs")
        config = current
        family.authenticate_producer({"producer_identity": record["producer_identity"]})
        if name == "prefix" and record["prepare_json_sha256"] != episode.file_sha256(directory/"prepare.json"):
            raise ValueError("extent predecessor changed")
        if name == "run" and record["prefix_json_sha256"] != episode.file_sha256(directory/"prefix.json"):
            raise ValueError("extent prefix predecessor changed")
        stages.append(name)
    previous_pins = {}
    for path in sorted(directory.glob("L*-*.json")):
        stem = path.stem; case_id, ordinal = stem.rsplit("-", 1)
        record, arrays = load(directory, case_id, int(ordinal))
        if config is None or record["input_sha256"] != config["inputs_sha256"]:
            raise ValueError("checkpoint belongs to different numerical inputs")
        if case_id in previous_pins and record["source_pins"] != previous_pins[case_id]:
            raise ValueError("trajectory changed a frozen extent/source/frame pin")
        previous_pins[case_id] = record["source_pins"]
        pair = pair_from_arrays(arrays, record)
        if pair.grid.length not in config["periods"] or pair.grid.nf != int(pair.grid.length*config["points_per_unit"]):
            raise ValueError("checkpoint resolution/period mismatch")
        if pair.parent_interval != PARENT or pair.child_interval != CHILD or pair.clock_locations != CLOCKS:
            raise ValueError("checkpoint changed physical observation windows")
        if np.max(np.abs(arrays["W"].T@arrays["W"]-np.eye(pair.grid.ng))) > 1e-8:
            raise ValueError("checkpoint W is not complete orthonormal")
        if "phi0" in arrays:
            model._active(pair, state_from(arrays))
        cases.append({"case_id": case_id, "ordinal": int(ordinal), "length": pair.grid.length})
    if not stages:
        raise FileNotFoundError("no extent stages")
    return {"schema": SCHEMA, "ok": True, "stages": stages, "checked": cases,
            "repaired": False, "evolved": False}
