"""Independent full-carrier check of a sealed homogeneous-turn prediction.

Only the existing four-geometry leading RHS is advanced. The retained PREPv2
standing frame is folded with the adjoint AP isometry, never reselected. This
is a finite-CAR leading-EFT diagnostic, not a vacuum match or bounce claim.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import numpy as np
from threadpoolctl import threadpool_limits

from . import nsc_discovery_leading_einstein as leading
from . import nsc_discovery_leading_step_control as control
from . import provenance

SCHEMA = "NSC-DISCOVERY-HOMOGENEOUS-TURN-v1"
NF = 32
MAX_CPU = 20.
MAX_TIME = 4.
MAX_BYTES = 64*1024**2
ROOT = leading.ROOT


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def authenticate_lock(path):
    """Read the creation-only lock before source preparation or any time step."""
    path = Path(path).resolve()
    record = json.loads(path.read_text())
    if (record.get("schema") != SCHEMA or record.get("mode") != "locked_prediction"
            or record.get("locked_before_held_data") is not True):
        raise ValueError("full-field advancement requires a LOCKED prediction")
    prep = record["preparation"]
    for name, expected in prep["input_hashes"].items():
        if digest(name) != expected:
            raise ValueError("locked source input changed: "+name)
    # The lock may itself contain forecast output, but it is never used to
    # choose r0, the trajectory end, the crossing bracket, or a root seed.
    producers = record.get("producers", {})
    required = {str(Path(__file__).resolve().relative_to(ROOT)),
                str(Path(leading.__file__).resolve().relative_to(ROOT)),
                str(Path(control.__file__).resolve().relative_to(ROOT))}
    if not required.issubset(producers):
        raise ValueError("locked prediction lacks full-field/RHS/step producer pins")
    for name, expected in producers.items():
        if digest(ROOT/name) != expected:
            raise ValueError("locked live producer changed: "+name)
        if record.get("producing_commit"):
            provenance.resolve_pinned_source_bytes(ROOT, name, expected,
                                                   commit=record["producing_commit"])
    for key in ("input_record", "input_payload"):
        name = str(Path(prep[key]).resolve())
        if name not in prep["input_hashes"]:
            raise ValueError("locked preparation lacks input pin: "+key)
    if prep.get("source_case") != "uniform" or prep.get("initial_momenta_zero") is not True:
        raise ValueError("turn check requires the declared uniform zero-momentum preparation")
    return record, {"prediction_path": str(path), "prediction_sha256": digest(path),
                    "producers": producers, "input_hashes": prep["input_hashes"]}


def fold_ap_frame(phi0, phi1, *, nf=NF, period=8.):
    """Canonical AP restriction preserves column order, phase and density units."""
    phi0, phi1 = np.asarray(phi0, complex), np.asarray(phi1, complex)
    if phi0.shape != phi1.shape or phi0.ndim != 2 or phi0.shape[1] != 6:
        raise ValueError("actual source must have six paired spinor columns")
    source_nf = phi0.shape[0]
    if source_nf < nf:
        raise ValueError("fold requires a retained finer source")
    lift = leading.backend.SpinorCarrier(nf, source_nf, period, canonical=True)
    folded = tuple(lift.conj().T@value for value in (phi0, phi1))
    frame = np.vstack(folded)
    original = np.vstack((phi0, phi1))
    remainder = original-np.vstack(tuple(lift@value for value in folded))
    gram = frame.conj().T@frame
    metadata = {"source_nf": source_nf, "nf": nf, "period": float(period),
        "method": "adjoint_existing_canonical_AP_isometry",
        "discarded_tail_frobenius": float(np.linalg.norm(remainder)),
        "discarded_tail_max": float(np.max(abs(remainder))),
        "discarded_power": float(np.linalg.norm(remainder)**2),
        "Gram_max": float(np.max(abs(gram-np.eye(6)))),
        "density_units": "canonical_half_density; norm_one_per_column",
        "original_frame_sha256": array_digest(original),
        "folded_frame_sha256": array_digest(frame),
        "column_order_phase_preserved": True, "orthonormalization_applied": False,
        "new_eigenframe_selected": False}
    if metadata["Gram_max"] > 1e-10 or metadata["discarded_tail_frobenius"] > 1e-10:
        raise ValueError("source does not close in the declared NF32 AP band")
    return folded, metadata


def prepare_source(source_preparation, held_weights):
    """Prepare the own zero-momentum radius from actual folded field energy."""
    prep = source_preparation
    path = Path(prep["input_record"])
    saved = json.loads(path.read_text())
    payload = Path(prep["input_payload"])
    if digest(payload) != saved["payload_sha256"] or saved["nf"] != 128:
        raise ValueError("requires the actual authenticated PREPv2 NF128 source")
    if float(prep["period"]) != saved["period"] or float(prep["period"]) != 8.:
        raise ValueError("turn diagnostic requires the unchanged period-eight carrier")
    weights = np.asarray(held_weights, float)
    if weights.shape == (3,):
        weights = np.repeat(weights, 2)
    if (weights.shape != (6,) or not np.isfinite(weights).all()
            or np.min(weights) <= 0 or np.max(weights) > 1
            or not np.array_equal(weights[::2], weights[1::2])):
        raise ValueError("positive equal-pair CAR occupations required")
    with np.load(payload, allow_pickle=False) as arrays:
        name = prep["source_case"]+"_nodal_"
        Q = np.array(arrays[name+"Q"], copy=True)
        for key in ("chi", "p_chi", "p_Q", "p_r"):
            if np.any(arrays[name+key] != 0):
                raise ValueError("PREPv2 transfer requires zero "+key)
        folded, frame = fold_ap_frame(arrays[name+"phi0"], arrays[name+"phi1"])
    Q0 = float(prep["Q0"])
    if Q0 <= 0 or np.max(abs(Q-Q0)) > 1e-10:
        raise ValueError("declared Q0 disagrees with the actual uniform preparation")
    pair = leading.nested.build_pair(NF, coarse_modes=1, child_details=2,
        source_layout="override", columns_override=folded, occupations=weights)
    pair = leading.backend.make_fft_pair(pair)
    zeros = np.zeros(pair.grid.ng)
    nodal = leading.State(np.full(pair.grid.ng, Q0), np.ones(pair.grid.ng),
                          zeros.copy(), zeros.copy(), *folded)
    fine = leading.prolong(pair.grid, nodal)
    system = leading.active_system(pair.grid, fine)
    a, Z, mag = leading.coefficients(system)
    if system.multiplicity != prep["multiplicity"] or system.kappa != prep["kappa"]:
        raise ValueError("locked multiplicity/kappa differs from the full RHS")
    constants = {"a": a, "Z": Z, "mag": mag, "A": system.A,
                 "C_F": system.C_F, "flux": system.flux}
    if any(abs(constants[key]-float(prep["constants"][key])) > 1e-13
           for key in constants):
        raise ValueError("locked leading-action constants changed")
    source = leading.coupling.source_from_columns(system, fine)
    rho = source["force_L"]/pair.grid.dx_q
    field_energy = leading.coupling.field_energy(system, fine)
    rho_mean = float(np.mean(rho))
    if np.max(abs(rho-rho_mean)) > 1e-9:
        raise ValueError("folded actual source energy is not spatially uniform")
    if abs(field_energy-pair.grid.length*Q0*rho_mean) > 1e-10:
        raise ValueError("full source energy/density normalization mismatch")
    r_squared = -(mag+rho_mean/Q0)/a
    if not np.isfinite(r_squared) or r_squared <= 0:
        raise ValueError("actual source has no positive initial radius")
    nodal.r[:] = np.sqrt(r_squared)
    state = leading.encode(pair, nodal)
    constraints = leading.constraints(pair, state)
    if max(constraints["raw_C_max"], constraints["D_max"]) > 1e-8:
        raise ValueError("analytic initial radius fails full source C/D constraints")
    metadata = {"Q0": Q0, "r0": float(np.sqrt(r_squared)), "initial_momenta_zero": True,
        "radius_formula": "r0^2=-(mag+actual_mean_rho/Q0)/a",
        "actual_initial_field_energy": field_energy, "actual_mean_rho": rho_mean,
        "rho_nonuniform_max": float(np.max(abs(rho-rho_mean))),
        "initial_constraints": constraints, "frame": frame,
        "source_weights": weights.tolist(), "source_trace": float(np.sum(weights)),
        "source_weights_sha256": array_digest(weights),
        "initial_full_state_sha256": {name: array_digest(getattr(state, name))
                                      for name in leading.FIELDS},
        "source_eta": 1., "multiplicity": int(system.multiplicity),
        "multiplicity_applied_once_outside_H": source["multiplicity_applied_once"],
        "constants": constants, "old_radius_copied": False,
        "Q0_fixed_preparation_limitation": True}
    return pair, state, metadata


def momentum(pair, state):
    return float(np.mean(leading.decode(pair, state).p_Q))


def positive_turn_bracket(previous, following, armed):
    """The zero initial maximum cannot count as the later radius minimum."""
    return bool(armed and previous > 0 and following <= 0)


def rk4_with_ledgers(pair, state, dt):
    """Original full RHS/RK4 with passive clocks and field work at its stages."""
    states = [state]
    rates = [leading.rates(pair, state)]
    states.append(leading.combine(state, rates[0], dt/2))
    rates.append(leading.rates(pair, states[-1]))
    states.append(leading.combine(state, rates[1], dt/2))
    rates.append(leading.rates(pair, states[-1]))
    states.append(leading.combine(state, rates[2], dt))
    rates.append(leading.rates(pair, states[-1]))
    weights = (1, 2, 2, 1)
    result = leading.State(*(getattr(state, name)+dt*sum(
        weight*getattr(rate, name) for weight, rate in zip(weights, rates))/6
        for name in leading.FIELDS))
    leading.check_chart(pair, result)
    clocks = dt*sum(weight*leading.clock_rates(pair, stage)
                    for weight, stage in zip(weights, states))/6
    work = dt*sum(weight*rate.fieldwork_power for weight, rate in zip(weights, rates))/6
    return result, clocks, float(work)


def root_turn(pair, left, dt, *, deadline):
    """Partial original RK4 steps from the preceding strictly positive p_Q state."""
    if momentum(pair, left) <= 0:
        raise ValueError("turn root needs a preceding positive p_Q state")
    right, _clock, _work = rk4_with_ledgers(pair, left, dt)
    if momentum(pair, right) > 0:
        raise ValueError("turn root lacks a positive-to-negative bracket")
    lo, hi = 0., dt
    for _ in range(40):
        if time.process_time() >= deadline:
            raise TimeoutError("CPU cap during full-field event root")
        mid = (lo+hi)/2
        trial, _clock, _work = rk4_with_ledgers(pair, left, mid)
        if momentum(pair, trial) > 0:
            lo = mid
        else:
            hi = mid
        if hi-lo < 2e-12:
            break
    fraction = (lo+hi)/2
    event, clock, work = rk4_with_ledgers(pair, left, fraction)
    return event, fraction, clock, work, {"partial_RK4_root": True,
        "preceding_positive_p_Q": momentum(pair, left), "right_nonpositive_p_Q": momentum(pair, right),
        "root_dt_bracket": [lo, hi], "root_time_width": hi-lo,
        "preceding_positive_state_energy": leading.energy(pair, left),
        "right_nonpositive_state_energy": leading.energy(pair, right),
        "sampled_minimum_substituted": False}


def _advance(pair, initial, cap, *, deadline, maximum_time=MAX_TIME, before_advance=None):
    """Bounded single full-field arm; caller supplies the authenticated lock guard."""
    if before_advance is None:
        raise ValueError("full-field advancement requires a lock guard")
    before_advance()
    state = initial.copy()
    instant, steps, work = 0., 0, 0.
    clocks = np.zeros(3)
    armed = False
    start = time.process_time()
    status, root = "NO_POSITIVE_TO_NEGATIVE_EVENT_IN_DECLARED_INTERVAL", None
    last_restriction = None
    blocker = None
    while instant < maximum_time:
        if time.process_time() >= deadline:
            status = "CPU_BUDGET_STOP"
            break
        try:
            dt, last_restriction = control.step_restriction(pair, state, cap)
            dt = min(dt, maximum_time-instant)
            previous = momentum(pair, state)
            armed = armed or previous > 0
            following, dc, dw = rk4_with_ledgers(pair, state, dt)
        except leading.coupling.PositiveChartExit as exc:
            status, blocker = "POSITIVE_CHART_EXIT", str(exc)
            break
        if positive_turn_bracket(previous, momentum(pair, following), armed):
            try:
                state, fraction, dc, dw, root = root_turn(pair, state, dt, deadline=deadline)
            except TimeoutError:
                status = "CPU_BUDGET_STOP"
                break
            instant += fraction
            clocks += dc
            work += dw
            steps += 1
            status = "FIRST_POSITIVE_TO_NEGATIVE_P_Q_EVENT"
            break
        state = following
        clocks += dc
        work += dw
        instant += dt
        steps += 1
    rate, bundle = leading.rates(pair, state, return_bundle=True)
    observation = leading.observe(pair, state, instant)
    initial_energy = leading.energy(pair, initial)
    nodal = leading.decode(pair, state)
    record = {"source": "held", "status": status, "step_cap": float(cap), "nf": NF,
        "time": instant, "radius": float(np.mean(nodal.r)), "Q": float(np.mean(nodal.Q)),
        "p_Q": float(np.mean(nodal.p_Q)), "p_r": float(np.mean(nodal.p_r)),
        "p_Q_dot": float(np.mean(leading.decode(pair, rate).p_Q)),
        "p_Q_spatial_spread": float(np.max(abs(nodal.p_Q-np.mean(nodal.p_Q)))),
        "radius_spatial_spread": float(np.max(abs(nodal.r-np.mean(nodal.r)))),
        "normal_clocks": clocks.tolist(), "integrated_coordinate_fieldwork": work,
        "field_energy_change": observation["energy"]["field"]-initial_energy["field"],
        "fieldwork_balance_residual": observation["energy"]["field"]-initial_energy["field"]-work,
        "total_energy_drift": observation["energy"]["total"]-initial_energy["total"],
        "initial_energy": initial_energy, "blocker": blocker,
        "event_root": root, "initial_zero_event_ignored": True, "armed_after_positive_p_Q": armed,
        "steps": steps, "last_scaled_restriction": last_restriction,
        "CPU_seconds": time.process_time()-start, "diagnostics": observation,
        "original_full_RHS": leading.rates.__module__, "reduced_ODE_imported": False,
        "actual_metric_curvature_measured": True, "bounce_claim": False,
        "vacuum_matching_completed": False, "EFT_strong_curvature_controlled": False}
    arrays = {name: np.array(getattr(state, name), copy=True) for name in leading.FIELDS}
    arrays.update(normal_clocks=clocks, source_weights=pair.weights.copy(), W=pair.geometry_map.copy())
    return record, arrays


def measure_turn(source_preparation, held_weights, step_caps=(.001, .0005),
                 cpu_budget=MAX_CPU, deadline=None, *, locked_prediction=None):
    """Two full NF32 arms; only an authenticated immutable lock permits stepping."""
    if locked_prediction is None:
        raise ValueError("measure_turn requires a LOCKED prediction path")
    locked, binding = authenticate_lock(locked_prediction)
    if source_preparation != locked["preparation"]:
        raise ValueError("source preparation differs from the sealed lock")
    expected = np.asarray(source_preparation["held_pair_occupations"], float)
    actual = np.asarray(held_weights, float)
    if actual.shape == (3,):
        actual = np.repeat(actual, 2)
    if not np.array_equal(actual, expected):
        raise ValueError("held occupations differ from the sealed lock")
    if not 0 < cpu_budget <= MAX_CPU or tuple(step_caps) != (.001, .0005):
        raise ValueError("full-field check uses two declared caps and at most20 CPU seconds")
    start = time.process_time()
    # Leave a small bounded reserve for serializing the last valid state and
    # computing actual endpoint curvature/constraints if the arm stops.
    stop = min(start+max(0., cpu_budget-.15),
               float(deadline) if deadline is not None else start+cpu_budget)
    records, arrays = [], {}
    with threadpool_limits(limits=1), leading.backend.fft_thread_limit(1):
        pair, initial, metadata = prepare_source(source_preparation, actual)
        for name in leading.FIELDS:
            arrays["initial_"+name] = np.array(getattr(initial, name), copy=True)
        for index, cap in enumerate(step_caps):
            row, values = _advance(pair, initial, cap, deadline=stop,
                before_advance=lambda: authenticate_lock(locked_prediction))
            records.append(row)
            arrays.update({f"cap{index}_"+name: value for name, value in values.items()})
            if row["status"] == "CPU_BUDGET_STOP":
                break
    indicators = None
    if len(records) == 2 and all(row["event_root"] for row in records):
        indicators = {name+"_cap_movement": abs(records[0][name]-records[1][name])
                      for name in ("time", "radius", "Q", "p_Q_dot")}
    report = {"schema": SCHEMA, "backend": "independent_full_NF32_leading_RK4",
        "events": records, "numerical_indicators": indicators, "preparation": metadata,
        "binding": binding, "CPU_seconds": time.process_time()-start,
        "cpu_limit_seconds": float(cpu_budget), "deadline_is_process_CPU": True,
        "endpoint_diagnostic_CPU_reserve_seconds": .15,
        "array_geometry_representation": "W coefficients; p_Q/p_r are canonical pi=dx_g W.T p",
        "locked_before_advancement": True, "forecast_output_used_as_initial_or_event_data": False,
        "no_sampled_minimum": True, "no_universal_error_gate": True,
        "scope": "finite CAR leading EFT; strong curvature and vacuum matching remain open"}
    if sum(value.nbytes for value in arrays.values()) > MAX_BYTES:
        raise ValueError("full-field arrays exceed64 MiB")
    return report, arrays


def measure_locked_prediction(prediction_json_path, *, cpu_limit_seconds=MAX_CPU, output_prefix=None):
    """Producer-facing interface; external owner writes the returned immutable data."""
    locked, _binding = authenticate_lock(prediction_json_path)
    # The external producer owns naming/writing. Existing paths must never be
    # silently treated as permission to overwrite data.
    if output_prefix is not None:
        prefix = Path(output_prefix)
        if prefix.with_suffix(".json").exists() or prefix.with_suffix(".npz").exists():
            raise FileExistsError("full-field output is creation-only")
    return measure_turn(locked["preparation"], locked["preparation"]["held_pair_occupations"],
                        cpu_budget=cpu_limit_seconds, locked_prediction=prediction_json_path)
