#!/usr/bin/env python3
"""Cached same-source response, then two bounded space-confirmation cases.

The resident v1 nf256 witness remains primary. This successor authenticates
it without rewriting/reintegrating it. nf512 baseline/parent-only cases only
confirm the declared space movement on the same initial-subtracted curves.
Only the coordinating root agent launches this numerical campaign.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np

import derive_nsc_nested_parent_child as producer
from recursive_horizons import nsc_nested_parent_child as model
from recursive_horizons import nsc_nested_parent_child_response as response


LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-nested-parent-child-confirmation-v2.json"
PAYLOAD = OUT.with_suffix(".npz")
PROGRESS = LAB / "results/development/nsc-nested-parent-child-confirmation-progress.json"
PROGRESS_PAYLOAD = PROGRESS.with_suffix(".npz")
SCHEMA = "NSC-NESTED-PARENT-CHILD-CONFIRMATION-v2"
CPU_LIMIT = 1800.
PAYLOAD_LIMIT = 64 * 1024**2
LOCAL_CPU_LIMIT = 300.
FINE_NF = 512
STEP_CAP = .0005
TARGET = .3
FRAME_COUNT = 31
PRIMARY_CASE = "nf256_baseline_dt0.0005"
PRIMARY_PHYSICAL_EFFECTS = ("parent_to_child_state", "child_to_parent_state",
                          "parent_to_child_geometry", "child_to_parent_annulus_geometry")
DEPENDENCIES = tuple(dict.fromkeys((*producer.DEPENDENCIES, producer.OUT,
                                  producer.PAYLOAD, Path(__file__))))


def source_bindings():
    return {path.relative_to(LAB).as_posix(): producer.digest(path) for path in DEPENDENCIES}


def array_digest(values):
    values = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(values.dtype).encode())
    digest.update(str(values.shape).encode())
    digest.update(values.tobytes())
    return digest.hexdigest()


def interpolate_state(states, times, mark):
    """Same linear canonical state interpolation as the frozen v1 consumer."""
    mark = float(mark)
    if mark < times[0] - 1e-12 or mark > times[-1] + 1e-12:
        raise ValueError("operator query leaves the authenticated trajectory")
    upper = min(max(int(np.searchsorted(times, mark, side="right")), 1), len(times) - 1)
    lower = upper - 1
    fraction = np.clip((mark - times[lower]) / (times[upper] - times[lower]), 0., 1.)
    return model.NestedState(*(getattr(states[lower], name) * (1 - fraction)
                              + getattr(states[upper], name) * fraction for name in model.STATE_NAMES))


class CachedLinearHamiltonian:
    """Exact cached H for the owned conformal LINEAR-Q interpolation schedule.

    There is no additional state/geometry approximation: the existing H_G
    is affine in Q, with a constant sigma2 P kinetic term. Full matrices are
    a runtime cache only; they are not duplicated in the evidence payload.
    """
    def __init__(self, pair, states, times, *, check_budget=None):
        self.times = np.array(times, dtype=float, copy=True)
        if self.times.ndim != 1 or len(states) != len(self.times) or len(states) < 2 or np.any(np.diff(self.times) <= 0):
            raise ValueError("cached operator requires matching states and increasing times")
        self.pair, self.states = pair, states
        matrices = []
        started = time.process_time()
        first = model.hamiltonian(pair, states[0])
        self.one_matrix_cpu_seconds = time.process_time() - started
        self.forecast_cpu_seconds = 1.5 * len(states) * self.one_matrix_cpu_seconds
        if check_budget is not None:
            check_budget(self.forecast_cpu_seconds)
        matrices.append(first)
        for state in states[1:]:
            if check_budget is not None:
                check_budget(0.)
            matrices.append(model.hamiltonian(pair, state))
        self.matrices = np.stack(matrices)
        self.build_cpu_seconds = time.process_time() - started
        self.callback_calls = 0

    def __call__(self, mark):
        mark = float(mark)
        if mark < self.times[0] - 1e-12 or mark > self.times[-1] + 1e-12:
            raise ValueError("operator query leaves the authenticated trajectory")
        upper = min(max(int(np.searchsorted(self.times, mark, side="right")), 1), len(self.times) - 1)
        lower = upper - 1
        fraction = np.clip((mark - self.times[lower]) / (self.times[upper] - self.times[lower]), 0., 1.)
        self.callback_calls += 1
        return (1 - fraction) * self.matrices[lower] + fraction * self.matrices[upper]

    def validate_midpoints(self, *, relative_tolerance=1e-12, check_budget=None):
        maximum = 0.
        for left, right in zip(self.times[:-1], self.times[1:]):
            if check_budget is not None:
                check_budget(0.)
            midpoint = (left + right) / 2
            actual = model.hamiltonian(self.pair, interpolate_state(self.states, self.times, midpoint))
            relative = float(np.linalg.norm(self(midpoint) - actual) / max(np.linalg.norm(actual), 1e-30))
            maximum = max(maximum, relative)
        if maximum > relative_tolerance:
            raise ValueError("cached operator changes the frozen linear-Q schedule")
        return {"midpoints_checked": len(self.times) - 1, "relative_max": maximum,
                "relative_tolerance": relative_tolerance, "same_linear_Q_schedule": True,
                "frames": len(self.times), "runtime_cache_bytes": int(self.matrices.nbytes),
                "matrix_build_cpu_seconds": self.build_cpu_seconds,
                "one_matrix_cpu_seconds": self.one_matrix_cpu_seconds,
                "matrix_build_forecast_cpu_seconds": self.forecast_cpu_seconds,
                "cached_matrices_saved_in_payload": False}


def load_witness():
    """Authenticate all resident v1 data using its frozen producer verifier."""
    record = producer.verify()
    case = record["results"][PRIMARY_CASE]
    pair = model.build_pair(256, occupations=producer.WEIGHTS["baseline"])
    with np.load(producer.PAYLOAD, allow_pickle=False) as stored:
        times = np.array(stored[PRIMARY_CASE + "_times"], copy=True)
        states = [model.NestedState(*(np.array(stored[PRIMARY_CASE + "_" + field][i], copy=True)
                                    for field in model.STATE_NAMES)) for i in range(len(times))]
        if not np.array_equal(pair.reference_columns, stored[PRIMARY_CASE + "_reference_columns"]):
            raise ValueError("fixed T0 observer differs from authenticated v1")
    if times.shape != (FRAME_COUNT,) or abs(times[0]) > 1e-12 or abs(times[-1] - TARGET) > 1e-12:
        raise ValueError("primary v1 trajectory is not the declared complete witness")
    if producer.state_hash(states[0]) != case["summary"]["initial_state_sha256"] or producer.state_hash(states[-1]) != case["summary"]["final_state_sha256"]:
        raise ValueError("primary saved Cauchy state authentication failed")
    return record, pair, states, times


def _raw_curve(case, observable):
    rows = case["rows"]
    if observable == "child_probability":
        values = [row["windows"]["child"]["probability"] for row in rows]
    elif observable == "child_occupation":
        values = [row["child_occupation"] for row in rows]
    else:
        values = [row["metrics"][observable] for row in rows]
    return np.asarray(values, dtype=float)


def control_curve(case_baseline, case_parent, observable):
    baseline, parent = _raw_curve(case_baseline, observable), _raw_curve(case_parent, observable)
    if baseline.shape != parent.shape:
        raise ValueError("confirmation controls use different sample domains")
    times_baseline = np.array([row["time"] for row in case_baseline["rows"]])
    times_parent = np.array([row["time"] for row in case_parent["rows"]])
    if not np.allclose(times_baseline, times_parent, rtol=0, atol=1e-12):
        raise ValueError("confirmation controls use different coordinate clocks")
    return (parent - parent[0]) - (baseline - baseline[0])


def confirmation_comparisons(witness, results):
    """Fine512 space indicators; primary256 trajectory and observables unchanged."""
    keys = {"parent_to_child_state": "child_probability",
            "parent_to_child_modal_state": "child_occupation",
            "parent_to_child_geometry": "child_r_proper_mean"}
    outputs = {}
    for name, observable in keys.items():
        primary = control_curve(witness["results"][PRIMARY_CASE],
                                witness["results"]["nf256_parent_only_dt0.0005"], observable)
        fine = control_curve(results["nf512_baseline_dt0.0005"],
                             results["nf512_parent_only_dt0.0005"], observable)
        if primary.shape != fine.shape:
            raise ValueError("space confirmation does not match the primary sample domain")
        effect, primary_effect = float(np.max(abs(fine))), float(np.max(abs(primary)))
        gap = float(np.max(abs(fine - primary)))
        inherited_time = witness["two_way_effects"][name]["time_indicator"]
        outputs[name] = {"observable": observable, "primary_case": "nf256_parent_only_dt0.0005",
            "confirmation_case": "nf512_parent_only_dt0.0005", "primary_effect_max": primary_effect,
            "fine_confirmation_effect_max": effect, "effect_max": effect, "space_indicator": gap,
            "space_fraction": None if effect == 0 else gap / effect,
            "fraction_of_primary_v1_effect": None if primary_effect == 0 else gap / primary_effect,
            "primary_time_indicator": inherited_time,
            "primary_time_fraction": None if primary_effect == 0 else inherited_time / primary_effect,
            "time_confirmation_claimed": False, "initial_offset_removed": True,
            "matched_proper_time_claimed": False,
            "resolved_at_one_percent": bool(effect > 0 and gap < .01 * effect and
                                             primary_effect > 0 and inherited_time < .01 * primary_effect),
            "primary_curve": primary, "fine_confirmation_curve": fine}
    return outputs


def case_quality(case):
    summary = case["summary"]
    rows = case["rows"]
    energy_scale = max(1., max(abs(row["field_energy"]) for row in rows))
    return {"target_reached": bool(summary["target_reached"]),
        "positive_chart": bool(summary["minimum_r"] > 0 and summary["minimum_Q"] > 0),
        "CAR_admissible": bool(summary["maximum_gram_gap"] < 1e-8 and
                               summary["minimum_covariance_eigenvalue"] >= -1e-10 and
                               summary["maximum_covariance_eigenvalue"] <= 1 + 1e-8),
        "full_source_accounting": bool(summary["maximum_field_energy_closure"] < 1e-8 * energy_scale),
        "sampled_constraint_budget_positive": bool(summary["sampled_constraint_end_margin"] >= 0),
        "energy_trace_closure_scale": energy_scale}


def _fraction_ok(row):
    return row.get("effect", 0) > 0 and row.get("fraction") is not None and row["fraction"] < .01


def local_quality(result):
    child, parent = result["child"], result["parent"]
    tests = {"child_reduction": _fraction_ok(child["reduction_error"]),
             "child_full_midpoint": _fraction_ok(child["full_midpoint_indicator"]),
             "conditional_autonomous": _fraction_ok(child["conditional_versus_autonomous"]),
             "parent_reduction": _fraction_ok(parent["reduction_error"]),
             "parent_full_midpoint": _fraction_ok(parent["full_midpoint_indicator"])}
    # A zero coherence change or roundoff-only omission is not a physical effect.
    roundoff = 128 * np.finfo(float).eps * max(1., float(np.max(abs(child["occupation_full"]))))
    if child["coherence_error"]["effect"] > roundoff:
        tests["child_coherence"] = _fraction_ok(child["coherence_error"])
    active = {}
    for name, control in child["controls"].items():
        is_active = control["occupation_movement"] > roundoff
        active[name] = is_active
        if is_active:
            tests[name + "_own_error"] = _fraction_ok(control["own_reference_error"])
            tests[name + "_midpoint"] = _fraction_ok(control["reference_midpoint_indicator"])
            if control["coherence_movement"] > roundoff:
                tests[name + "_coherence_error"] = _fraction_ok(control["coherence_reference_error"])
    return {"checks": tests, "active_controls": active, "roundoff_activity_scale": roundoff,
            "all_within_one_percent": all(tests.values()), "observer_preserved": bool(result["observer_preserved"]),
            "all_source_components_retained": bool(result["all_source_components_retained"]),
            "no_dense_propagator_history": all(value["time_indexed_exterior_propagator_bytes"] == 0
                                               for value in result["allocation"].values())}


def goal_complete(record):
    """Exactly the four named physical effects, source accounting and local errors."""
    physical = record.get("primary_physical_effects", {})
    quality = record.get("case_quality", {})
    local = record.get("local_quality", {})
    return bool(record.get("source_hashes_unchanged") is True and
                record.get("cached_operator_identity", {}).get("same_linear_Q_schedule") is True and
                record.get("endpoint_nested_schur", {}).get("nested_direct_gap", 1) < 1e-10 and
                record.get("endpoint_nested_schur", {}).get("joined_direct_gap", 1) < 1e-10 and
                set(physical) == set(PRIMARY_PHYSICAL_EFFECTS) and
                all(value["resolved_at_one_percent"] for value in physical.values()) and
                set(record.get("results", {})) == {"nf512_baseline_dt0.0005", "nf512_parent_only_dt0.0005"} and
                set(quality) == set(record.get("required_case_names", ())) and PRIMARY_CASE in quality and
                len(quality) >= 3 and all(all(value[key] for key in ("target_reached", "positive_chart", "CAR_admissible", "full_source_accounting"))
                                        for value in quality.values()) and
                local.get("all_within_one_percent") is True and local.get("observer_preserved") is True and
                local.get("all_source_components_retained") is True and local.get("no_dense_propagator_history") is True and
                record.get("local_response_saved") is True and record.get("confirmation_states_saved") is True)


def pack_arrays(value, arrays, prefix):
    """Bind every numerical local curve without duplicating it as JSON lists."""
    if isinstance(value, np.ndarray):
        arrays[prefix] = np.array(value, copy=True)
        return {"array": prefix, "shape": list(value.shape), "dtype": str(value.dtype),
                "sha256": array_digest(value)}
    if isinstance(value, dict):
        return {str(key): pack_arrays(item, arrays, prefix + "_" + str(key)) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [pack_arrays(item, arrays, prefix + "_" + str(i)) for i, item in enumerate(value)]
    return producer.jsonable(value)


def unpack_arrays(value, arrays):
    if isinstance(value, dict) and set(value) == {"array", "shape", "dtype", "sha256"}:
        array = arrays[value["array"]]
        if list(array.shape) != value["shape"] or str(array.dtype) != value["dtype"] or array_digest(array) != value["sha256"]:
            raise ValueError("saved numerical curve binding changed: " + value["array"])
        return array
    if isinstance(value, dict):
        return {key: unpack_arrays(item, arrays) for key, item in value.items()}
    if isinstance(value, list):
        return [unpack_arrays(item, arrays) for item in value]
    return value


class Checkpoints:
    """Progress is separate; existing stopped and completed evidence is immutable."""
    def __init__(self, output=OUT, progress=PROGRESS, *, payload_limit=PAYLOAD_LIMIT):
        self.output, self.progress = Path(output), Path(progress)
        self.payload_limit = payload_limit
        for path in (self.output, self.output.with_suffix(".npz"), self.progress, self.progress.with_suffix(".npz")):
            if path.exists():
                raise FileExistsError("confirmation evidence already exists; use --check or a named successor")

    def save(self, record, arrays, *, final=False):
        destination = self.output if final else self.progress
        payload = destination.with_suffix(".npz")
        if final and (destination.exists() or payload.exists()):
            raise FileExistsError("refusing to overwrite final confirmation evidence")
        payload.parent.mkdir(parents=True, exist_ok=True)
        temporary = payload.with_suffix(".npz.tmp")
        with temporary.open("wb") as stream:
            np.savez_compressed(stream, **arrays)
        output = dict(record, payload_sha256=producer.digest(temporary),
                      payload_npz_bytes=temporary.stat().st_size,
                      final_artifact=bool(final), payload_bytes=0)
        for _ in range(4):
            text = json.dumps(producer.jsonable(output), indent=2, allow_nan=False) + "\n"
            total = temporary.stat().st_size + len(text.encode())
            if total == output["payload_bytes"]:
                break
            output["payload_bytes"] = total
        text = json.dumps(producer.jsonable(output), indent=2, allow_nan=False) + "\n"
        if temporary.stat().st_size + len(text.encode()) > self.payload_limit:
            temporary.unlink()
            raise RuntimeError("confirmation record plus payload exceeds 64 MiB")
        json_temporary = destination.with_suffix(".json.tmp")
        json_temporary.write_text(text)
        os.replace(temporary, payload)
        os.replace(json_temporary, destination)
        return output


def local_on_witness(pair, states, times, record, arrays, writer, check_budget, *, local_cpu_limit=LOCAL_CPU_LIMIT):
    print("begin cached same-source local response", PRIMARY_CASE, "T=.1->.3", flush=True)
    cache = CachedLinearHamiltonian(pair, states, times, check_budget=check_budget)
    record["cached_operator_identity"] = cache.validate_midpoints(check_budget=check_budget)
    hierarchy = response.nested_frame(pair.reference_columns)
    record["endpoint_nested_schur"] = pack_arrays(response.nested_schur(cache(times[-1]), hierarchy), arrays, "endpoint_schur")
    start = int(np.flatnonzero(np.isclose(times, .1, rtol=0, atol=1e-12))[0])
    actual_times, actual_states = times[start:], states[start:]
    actual = np.array([producer.state_columns(state) for state in actual_states])
    record["local_preparation"] = {"source": PRIMARY_CASE, "window": [.1, .3],
        "initial_state_sha256": producer.state_hash(actual_states[0]),
        "observer_sha256": array_digest(pair.reference_columns), "actual_Phi_t_start": True,
        "observer_reset": False, "geometry_rerun": False, "weights": pair.weights.tolist()}
    # Complete and checkpoint this named blocker before any new fine trajectory.
    for factor in (1, 2, 4):
        check_budget(0.)
        selected_times = np.linspace(.1, TARGET, factor * (len(actual_times) - 1) + 1)
        before = time.process_time()
        print("admitted cached local response", "output factor", factor, "CPU ceiling", local_cpu_limit, flush=True)
        kwargs = {"autonomous_columns": actual} if factor == 1 else {}
        result = response.retained_response(cache, selected_times, actual[0], pair.weights, hierarchy,
            substeps=2, reference_substeps=4, cpu_limit_s=local_cpu_limit, **kwargs)
        if factor > 1:
            child_actual = np.einsum("ij,tjk->tik", hierarchy.child.conj().T, actual)
            actual_covariance = np.einsum("tak,k,tbk->tab", child_actual, pair.weights, child_actual.conj())
            result["child"]["autonomous_covariance"] = actual_covariance
            full_at_frames = result["child"]["occupation_full"][::factor]
            effect = result["child"]["reduction_error"]["effect"]
            discrepancy = float(np.max(abs(full_at_frames - np.diagonal(actual_covariance, axis1=1, axis2=2).real)))
            result["child"]["conditional_versus_autonomous"] = {
                "movement": discrepancy, "effect": effect, "fraction": None if effect == 0 else discrepancy / effect,
                "within_one_percent": bool(effect > 0 and discrepancy < .01 * effect)}
        record.setdefault("local_trials", []).append({"output_factor": factor, "cpu_seconds": time.process_time() - before,
                                                     "quality": local_quality(result)})
        record["local_quality"] = local_quality(result)
        record["local_response"] = pack_arrays(result, arrays, "local_response")
        record["local_response_saved"] = True
        record["stage"] = "local_response_complete" if record["local_quality"]["all_within_one_percent"] else "local_response_unrefined"
        writer.save(record, arrays)
        print("cached local response finished", "output factor", factor,
              "all within one percent", record["local_quality"]["all_within_one_percent"], flush=True)
        if record["local_quality"]["all_within_one_percent"]:
            return
        forecast = 4 * (time.process_time() - before)
        check_budget(forecast)
    raise RuntimeError("local response remains unresolved after the declared bounded refinements")


def run(*, cpu_limit=CPU_LIMIT, local_cpu_limit=LOCAL_CPU_LIMIT):
    if not np.isfinite(cpu_limit) or not np.isfinite(local_cpu_limit) or cpu_limit <= 0 or cpu_limit > CPU_LIMIT or local_cpu_limit <= 0 or local_cpu_limit > LOCAL_CPU_LIMIT:
        raise ValueError("requested CPU ceilings exceed the admitted confirmation budget")
    writer = Checkpoints()
    started = time.process_time()
    before = source_bindings()
    witness, pair, states, times = load_witness()
    record = {"schema": SCHEMA, "status": "RUNNING", "stage": "authenticate_primary",
        "primary_witness": "nsc-nested-parent-child-v1", "predecessor_checkpoint": "ec0cd25",
        "primary_witness_case": PRIMARY_CASE, "primary_witness_reintegrated": False,
        "cpu_budget_seconds": cpu_limit, "local_cpu_budget_seconds": local_cpu_limit,
        "payload_cap_bytes": PAYLOAD_LIMIT, "source_bindings_before": before,
        "comparison_clock": "same global coordinate T; proper clocks remain separate recorded curves",
        "geometry_prescribed": False, "source_reset": False, "observer_reset": False,
        "results": {}, "forecasts": {},
        "case_quality": {name: case_quality(case) for name, case in witness["results"].items()},
        "required_case_names": list(witness["results"]) + ["nf512_baseline_dt0.0005", "nf512_parent_only_dt0.0005"],
        "local_response_saved": False, "confirmation_states_saved": False, "goal_complete": False}
    arrays = {}

    def check_budget(forecast=0.):
        record["cpu_seconds"] = time.process_time() - started
        if record["cpu_seconds"] + forecast + 10 >= cpu_limit:
            raise RuntimeError("confirmation CPU resource stop; separate progress checkpoint preserved")

    try:
        writer.save(record, arrays)
        local_on_witness(pair, states, times, record, arrays, writer, check_budget, local_cpu_limit=local_cpu_limit)
        for control in ("baseline", "parent_only"):
            name = f"nf512_{control}_dt0.0005"
            record["stage"] = "initial_state:" + name
            check_budget()
            fine_pair = model.build_pair(FINE_NF, occupations=producer.WEIGHTS[control])
            initial, preparation = model.initial_state(fine_pair)
            check_budget()
            hierarchy = response.nested_frame(fine_pair.reference_columns)
            dt, omega, _ = model.stable_timestep(fine_pair, initial, STEP_CAP)
            pilot_started = time.process_time()
            model.rk4_step(fine_pair, initial, dt)
            step_cpu = time.process_time() - pilot_started
            frame_started = time.process_time()
            first = producer.observe(fine_pair, initial, 0., hierarchy)
            frame_cpu = time.process_time() - frame_started
            forecast = float(1.5 * (np.ceil(TARGET / dt) * step_cpu + FRAME_COUNT * frame_cpu))
            record["forecasts"][name] = {"step_cpu_seconds": step_cpu, "frame_cpu_seconds": frame_cpu,
                "omega": omega, "actual_pilot_step": dt, "estimated_steps": int(np.ceil(TARGET / dt)),
                "forecast_cpu_seconds": forecast, "admitted": False}
            writer.save(record, arrays)
            check_budget(forecast)
            record["forecasts"][name]["admitted"] = True
            print("admitted", name, "CPU forecast", round(forecast, 3), flush=True)
            state, mark = initial.copy(), 0.
            frame_times = np.linspace(0., TARGET, FRAME_COUNT)
            case_states, rows = [initial.copy()], [first]
            record["stage"] = "evolve:" + name
            for frame_time in frame_times[1:]:
                while mark < frame_time - 1e-13:
                    check_budget()
                    step, _, _ = model.stable_timestep(fine_pair, state, STEP_CAP)
                    step = min(step, frame_time - mark)
                    state = model.rk4_step(fine_pair, state, step)
                    mark += step
                    for field in model.STATE_NAMES:
                        arrays["interrupted_checkpoint_" + field] = np.array(getattr(state, field), copy=True)
                    record["interrupted_case"], record["interrupted_time"] = name, mark
                case_states.append(state.copy())
                rows.append(producer.observe(fine_pair, state, mark, hierarchy))
                if rows[-1]["gram_gap"] > 1e-8 or rows[-1]["metrics"]["r_min"] <= 0 or rows[-1]["metrics"]["Q_min"] <= 0:
                    raise RuntimeError("fine confirmation lost its recorded chart/CAR admissibility")
                for field in model.STATE_NAMES:
                    arrays[name + "_" + field] = np.array([getattr(saved, field) for saved in case_states])
                arrays[name + "_times"] = frame_times[:len(case_states)]
                record["partial_case"] = {"name": name, "frames": len(case_states), "time": mark, "rows": rows}
                writer.save(record, arrays)
            case = {"nf": FINE_NF, "step_cap": STEP_CAP, "control": control, "preparation": preparation,
                    "source": fine_pair.source_metadata, "geometry_basis": fine_pair.geometry_metadata,
                    "rows": rows, "summary": producer.summary(fine_pair, initial, state, rows)}
            record["results"][name] = case
            record["case_quality"][name] = case_quality(case)
            arrays[name + "_reference_columns"] = fine_pair.reference_columns
            record.pop("partial_case", None)
            writer.save(record, arrays)
            print("completed", name, "time", mark, flush=True)
        record["confirmation_states_saved"] = True
        comparisons = confirmation_comparisons(witness, record["results"])
        record["space_confirmation"] = pack_arrays(comparisons, arrays, "space_confirmation")
        physical = {name: dict(witness["two_way_effects"][name]) for name in PRIMARY_PHYSICAL_EFFECTS}
        physical["parent_to_child_state"] = {key: value for key, value in comparisons["parent_to_child_state"].items()
                                             if not isinstance(value, np.ndarray)}
        record["primary_physical_effects"] = physical
        record["source_bindings_after"] = source_bindings()
        record["source_hashes_unchanged"] = before == record["source_bindings_after"]
        if not record["source_hashes_unchanged"]:
            raise RuntimeError("scientific source changed during confirmation")
        record["goal_complete"] = goal_complete(record)
        record["status"] = "MEASURED_CONFIRMED_NESTED_PAIR" if record["goal_complete"] else "MEASURED_UNRESOLVED_CONFIRMATION"
        record["stage"] = "complete_confirmation"
        record["cpu_seconds"] = time.process_time() - started
        writer.save(record, arrays, final=True)
    except (RuntimeError, ValueError, producer.coupling.PositiveChartExit) as error:
        record["status"], record["stop_reason"] = "STOPPED_CONFIRMATION", str(error)
        record["source_bindings_after"] = source_bindings()
        record["source_hashes_unchanged"] = before == record["source_bindings_after"]
        record["cpu_seconds"] = time.process_time() - started
        writer.save(record, arrays)
    return record


def verify(*, progress=False):
    path = PROGRESS if progress else OUT
    record = json.loads(path.read_text())
    if record.get("schema") != SCHEMA or producer.digest(path.with_suffix(".npz")) != record["payload_sha256"]:
        raise ValueError("invalid confirmation schema or payload")
    current = source_bindings()
    if record["source_bindings_before"] != current or record.get("source_bindings_after", current) != current:
        raise ValueError("confirmation or v1 source binding changed")
    witness = producer.verify()
    with np.load(path.with_suffix(".npz"), allow_pickle=False) as saved:
        arrays = {name: np.array(saved[name]) for name in saved.files}
    decoded = unpack_arrays(record, arrays)
    if decoded.get("local_response_saved"):
        if local_quality(decoded["local_response"]) != decoded["local_quality"]:
            raise ValueError("local response assessment does not replay")
        child = decoded["local_response"]["child"]
        for name, control in child["controls"].items():
            measured = float(np.max(abs(child[name + "_full"] - child["occupation_full"])))
            error = float(np.max(abs((child[name] - child["occupation"]) - (child[name + "_full"] - child["occupation_full"]))))
            if abs(measured - control["occupation_movement"]) > 1e-12 or abs(error - control["own_reference_error"]["movement"]) > 1e-12:
                raise ValueError("saved own omission effect/error does not replay: " + name)
    for name, case in decoded["results"].items():
        initial, final = [model.NestedState(*(arrays[name + "_" + field][i] for field in model.STATE_NAMES)) for i in (0, -1)]
        if producer.state_hash(initial) != case["summary"]["initial_state_sha256"] or producer.state_hash(final) != case["summary"]["final_state_sha256"]:
            raise ValueError("fine confirmation Cauchy binding changed: " + name)
        if case_quality(case) != decoded["case_quality"][name]:
            raise ValueError("fine confirmation chart/CAR/accounting assessment changed")
    if decoded.get("confirmation_states_saved"):
        replay = producer.jsonable(confirmation_comparisons(witness, decoded["results"]))
        if replay != producer.jsonable(decoded["space_confirmation"]):
            raise ValueError("space confirmation effect does not replay")
    if decoded.get("goal_complete", False) != goal_complete(decoded):
        raise ValueError("combined goal predicate does not replay")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--progress", action="store_true", help="verify the separate interrupted/running checkpoint")
    parser.add_argument("--cpu-limit", type=float, default=CPU_LIMIT)
    parser.add_argument("--local-cpu-limit", type=float, default=LOCAL_CPU_LIMIT)
    args = parser.parse_args()
    if args.progress and not args.check:
        parser.error("--progress is a read-only check option")
    result = verify(progress=args.progress) if args.check else run(cpu_limit=args.cpu_limit, local_cpu_limit=args.local_cpu_limit)
    print(result["status"], "CPU", result.get("cpu_seconds"), "goal complete", result.get("goal_complete"), flush=True)
