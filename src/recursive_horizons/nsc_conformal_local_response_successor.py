"""Same-realization local response on a finalized conformal continuation.

The original feasibility producer is unchanged. A 12-column linear solve
transports the retained/exterior parts of the same six source columns.
Their coherent sum is the actual state; their covariance sum omits only
the actual initial cross. No source or geometry is reintegrated.
"""
from __future__ import annotations

from collections import OrderedDict
import json
from pathlib import Path
import time

import numpy as np

from .nsc_conformal_local_response import (
    CASES, ConformalSchedule, PRIMARY, action_gate, full_split, initial_blocks,
    _covariance, _diagonal, _fraction, _max_abs,
)
from .nsc_coupled_local_response import (
    allocation_report, array_id, autonomous_series, evolve_full_columns,
    evolve_streamed, hermite_with_slopes, occupations_and_coherence,
    phase_error, sha256_file, signal_error,
)
from .nsc_spherical_coupling import KAPPA, OCCUPATIONS
from .nsc_spherical_galerkin_coupling import periodic_interpolation

LAB = Path(__file__).resolve().parents[2]
DEV = LAB / "results" / "development"
SOURCE_JSON = DEV / "nsc-spherical-conformal-episode-v2.json"
SOURCE_NPZ = SOURCE_JSON.with_suffix(".npz")
PREDECESSOR_JSON = DEV / "nsc-spherical-conformal-episode-v1.json"
PREDECESSOR_NPZ = PREDECESSOR_JSON.with_suffix(".npz")
TRANSPORT_JSON = DEV / "nsc-spherical-conformal-transport-v2.json"
V5_NPZ = DEV / "nsc-spherical-coupling-refinement-v5.npz"
FEASIBILITY_JSON = DEV / "nsc-conformal-local-response-feasibility-v1.json"
FEASIBILITY_NPZ = FEASIBILITY_JSON.with_suffix(".npz")
OUTPUT_JSON = DEV / "nsc-conformal-local-response-v2.json"
OUTPUT_NPZ = OUTPUT_JSON.with_suffix(".npz")
SCHEMA = "NSC-CONFORMAL-LOCAL-RESPONSE-v2"
WINDOW = (.2, .3)
PAYLOAD_LIMIT = 64 * 1024 * 1024
STATE_FIELDS = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")


def source_bindings():
    return {path.name: sha256_file(path) for path in (
        SOURCE_JSON, SOURCE_NPZ, PREDECESSOR_JSON, PREDECESSOR_NPZ,
        TRANSPORT_JSON, V5_NPZ, FEASIBILITY_JSON, FEASIBILITY_NPZ,
    )}


def readiness_errors(episode, transport, bindings):
    errors = []
    if episode.get("schema") != "NSC-SPHERICAL-CONFORMAL-EPISODE-v2" or episode.get("gauge") != "conformal":
        errors.append("source_schema_or_gauge")
    if episode.get("all_target_reached") is not True:
        errors.append("unfinished_episode")
    if episode.get("source_changed") is not False or episode.get("midtrajectory_reset") is not False:
        errors.append("source_or_state_reset")
    if episode.get("sealed_sources_unchanged_during_run") is not True:
        errors.append("source_binding")
    if episode.get("payload_sha256") != bindings[SOURCE_NPZ.name]:
        errors.append("source_payload_hash")
    predecessor = episode.get("predecessor_hashes", {})
    if predecessor.get("v1_npz") != bindings[PREDECESSOR_NPZ.name] or predecessor.get("v5_npz") != bindings[V5_NPZ.name]:
        errors.append("preparation_hash")
    if transport.get("schema") != "NSC-SPHERICAL-CONFORMAL-TRANSPORT-v2":
        errors.append("transport_schema")
    if transport.get("resolved_regional_surface_exchange") is not True or transport.get("maintained_localization") is not True:
        errors.append("regional_exchange_unresolved")
    physical_binding = transport.get("source_bindings", {})
    if physical_binding.get("episode_v2") != bindings[SOURCE_JSON.name] or physical_binding.get("payload_v2") != bindings[SOURCE_NPZ.name]:
        errors.append("transport_source_binding")
    for case in CASES:
        run = episode.get("results", {}).get(case, {})
        if run.get("target_reached") is not True or run.get("positive_chart") is not True:
            errors.append("chart_or_completion:" + case)
        if run.get("constraint_end_margin", -1) <= 0 or run.get("max_sampled_budget_violation", 1) > 1e-10:
            errors.append("sampled_constraint_budget:" + case)
        if run.get("gram_gap_max", 1) > 1e-6:
            errors.append("gram:" + case)
        if episode.get("handoffs", {}).get(case, {}).get("bitwise_equal") is not True:
            errors.append("handoff:" + case)
    return errors


def load_inputs(case=PRIMARY, window=WINDOW):
    if case not in CASES:
        raise ValueError("unrecorded conformal continuation case")
    bindings = source_bindings()
    episode = json.loads(SOURCE_JSON.read_text())
    transport = json.loads(TRANSPORT_JSON.read_text())
    errors = readiness_errors(episode, transport, bindings)
    if errors:
        raise RuntimeError("production local response held: " + ", ".join(errors))
    nf = int(case.split("_")[0][2:])
    prefix = "" if nf == 256 else "nf512_"
    start, end = map(float, window)
    with np.load(SOURCE_NPZ, allow_pickle=False) as source, np.load(PREDECESSOR_NPZ, allow_pickle=False) as predecessor, np.load(V5_NPZ, allow_pickle=False) as original:
        weights = np.array(original[prefix + "occupations"], copy=True)
        if not np.array_equal(weights, OCCUPATIONS):
            raise RuntimeError("original weights changed")
        for name in STATE_FIELDS:
            origin_key = ("columns_" if name.startswith("phi") else "geometry_") + name
            if not np.array_equal(source[f"nf{nf}_original_T0_{name}"], original[prefix + origin_key]):
                raise RuntimeError("original T0 preparation changed")
            if not np.array_equal(source[case + "_frames_" + name][0], predecessor[case + "_final_" + name]):
                raise RuntimeError("continuation is not the bitwise predecessor handoff")
        observer = np.vstack((source[f"nf{nf}_original_T0_phi0"], source[f"nf{nf}_original_T0_phi1"]))[:, :2].copy()
        times = np.array(source[case + "_frame_times"], copy=True)
        selected = (times >= start - 1e-12) & (times <= end + 1e-12)
        times = times[selected]
        if times.size < 3 or abs(times[0] - start) > 1e-12 or abs(times[-1] - end) > 1e-12:
            raise RuntimeError("requested comparison window is not stored")
        phi0 = np.array(source[case + "_frames_phi0"][selected], copy=True)
        phi1 = np.array(source[case + "_frames_phi1"][selected], copy=True)
        q = np.array(source[case + "_frames_Q"][selected], copy=True)
        slopes = np.array(source[case + "_frames_rate_Q"][selected], copy=True)
        lapse = np.array(source[case + "_frames_L"][selected], copy=True)
        lapse_rate = np.array(source[case + "_frames_L_dot"][selected], copy=True)
        if not np.array_equal(slopes, lapse_rate):
            raise RuntimeError("saved conformal lapse rate is not the actual lifted Q rate")
        proper_clock = np.array(source[case + "_proper_clock_frames"][selected], copy=True)
        inputs = {
            "case": case, "nf": nf, "nq": int(original[prefix + "nq"]), "length": float(original[prefix + "length"]),
            "weights": weights, "observer": observer, "columns": np.vstack((phi0[0], phi1[0])),
            "times": times, "Q": q, "slopes": slopes, "lapse": lapse,
            "phi0": phi0, "phi1": phi1, "proper_clock": proper_clock,
            "source_bindings": bindings, "episode": episode, "transport": transport,
        }
    if np.array_equal(observer, inputs["columns"][:, :2]):
        raise RuntimeError("observer reset to the transported segment basis")
    if source_bindings() != bindings:
        raise RuntimeError("source changed while the successor binder read it")
    return inputs


def shared_stream(hamiltonian, observer, times, columns, weights):
    """One linear solve for both pieces; the coherent sum keeps actual C_AE."""
    a0 = observer @ (observer.conj().T @ columns)
    e0 = columns - a0
    width = columns.shape[1]
    result, _, _ = evolve_streamed(hamiltonian, observer, times,
        np.concatenate((a0, e0), axis=1), np.tile(weights, 2), propagator_substeps=2)
    a, e = result["amplitudes"][:, :, :width], result["amplitudes"][:, :, width:]
    amplitudes = a + e
    occupation, coherence = occupations_and_coherence(amplitudes, weights)
    total = _covariance(amplitudes, weights)
    dropped = _covariance(a, weights) + _covariance(e, weights)
    drive_off, _ = occupations_and_coherence(a, weights)
    return {"result": result, "amplitudes": amplitudes, "occupation": occupation, "coherence": coherence,
            "covariance": total, "covariance_without_cross": dropped, "cross_covariance": total - dropped,
            "occupation_drive_off": drive_off,
            "original_column_history": result["history"][:, :, :width] + result["history"][:, :, width:]}


def fine_geometry_hamiltonian(schedule, values, slopes=None):
    """Same finite operator and initial state, changing only declared geometry."""
    values = np.asarray(values)
    masses = None if slopes is not None else np.stack([
        KAPPA * (schedule.columns_map.conj().T @ (row[:, None] * schedule.columns_map)) for row in values
    ])
    cache = OrderedDict()

    def hamiltonian(mark):
        key = float(mark).hex()
        if key in cache:
            cache.move_to_end(key)
            return cache[key]
        if slopes is not None:
            q = hermite_with_slopes(schedule.times, values, slopes, [mark])[0]
            if np.min(q) <= 0:
                raise RuntimeError("rate-Hermite geometry leaves the positive chart")
            mass = KAPPA * (schedule.columns_map.conj().T @ (q[:, None] * schedule.columns_map))
        else:
            index, fraction = schedule._interval(mark)
            mass = (1 - fraction) * masses[index] + fraction * masses[index + 1]
        mass = .5 * (mass + mass.conj().T)
        zero = np.zeros_like(mass)
        matrix = np.block([[zero, mass - 1j * schedule.momentum], [mass + 1j * schedule.momentum, zero]])
        cache[key] = matrix
        if len(cache) > 2:
            cache.popitem(last=False)
        return matrix

    return hamiltonian


def full_local(hamiltonian, inputs, times, substeps=2):
    propagated = evolve_full_columns(hamiltonian, times, inputs["columns"], substeps=substeps)
    projected = np.einsum("ij,tjk->tik", inputs["observer"].conj().T, propagated)
    return occupations_and_coherence(projected, inputs["weights"])


def physical_scope(inputs):
    """Use physical surface flux, including its already-owned nodal/dx map."""
    selected = {}
    start, end = inputs["times"][[0, -1]]
    for case in CASES:
        selected[case] = [row for row in inputs["episode"]["physical_series"][case]
                          if start - 1e-12 <= row["time"] <= end + 1e-12]
    surfaces = []
    for x, window, channel in ((0, 0, "normal_left_outward_flux"), (2, 0, "normal_right_outward_flux"), (4, 1, "normal_right_outward_flux")):
        measurements = {}
        for case, rows in selected.items():
            times = np.array([row["time"] for row in rows])
            values = np.abs([row["windows"][window][channel] for row in rows])
            dt = np.diff(times)
            trap = float(np.sum(.5 * dt * (values[:-1] + values[1:])))
            simpson = float(dt[0] / 3 * (values[0] + values[-1] + 4 * np.sum(values[1:-1:2]) + 2 * np.sum(values[2:-1:2])))
            measurements[case] = {"absolute_integral": trap, "simpson_indicator": abs(trap - simpson)}
        effect = measurements[PRIMARY]["absolute_integral"]
        indicators = {
            "time": abs(effect - measurements["nf512_dt_0.0010"]["absolute_integral"]),
            "space": abs(effect - measurements["nf256_dt_0.0005"]["absolute_integral"]),
            "frames": measurements[PRIMARY]["simpson_indicator"],
        }
        surfaces.append({"x": x, "quantity": "time-integrated absolute normal energy flux", "effect": effect,
                         "indicators": {name: _fraction(value, effect) for name, value in indicators.items()},
                         "per_case": measurements, "certified": False})
    rows = selected[PRIMARY]
    return {
        "surfaces": surfaces, "leader_region": 0,
        "leader_unchanged": all(row["shell_leader"] == 0 and row["probability_leader"] == 0 for row in rows),
        "shell_leader_share_min": min(row["shell_leader_share"] for row in rows),
        "probability_leader_share_min": min(row["probability_leader_share"] for row in rows),
        "proper_clock_start": float(inputs["proper_clock"][0]), "proper_clock_end": float(inputs["proper_clock"][-1]),
        "normal_clock_protocol": "continuous T0 normal worldline x=1; d tau=r Q dt; occupation is a separate fixed two-mode measurement",
        "flux_normalization": "source physical_series divides flux_nodal by the quadrature dx; energy flow is not occupation",
    }


def _save(record, arrays):
    if sum(np.asarray(value).nbytes for value in arrays.values()) > PAYLOAD_LIMIT:
        raise RuntimeError("consumer arrays exceed the separate 64 MiB payload cap")
    with OUTPUT_NPZ.with_suffix(".npz.tmp").open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    OUTPUT_NPZ.with_suffix(".npz.tmp").replace(OUTPUT_NPZ)
    record["payload_sha256"] = sha256_file(OUTPUT_NPZ)
    OUTPUT_JSON.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    for _ in range(3):
        record["payload_bytes"] = OUTPUT_JSON.stat().st_size + OUTPUT_NPZ.stat().st_size
        OUTPUT_JSON.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    if record["payload_bytes"] > PAYLOAD_LIMIT:
        raise RuntimeError("consumer payload exceeds its separate 64 MiB cap")


def execute(pilot_budget_s=300., confirmation_budget_s=1800.):
    if OUTPUT_JSON.exists() or OUTPUT_NPZ.exists():
        raise FileExistsError("new production record already exists; use readonly check or a separately named successor")
    started = time.process_time()
    inputs = load_inputs()
    schedule = ConformalSchedule(inputs)
    if _max_abs(schedule.fine_Q - inputs["lapse"]) > 1e-12:
        raise RuntimeError("prolonged geometry is not the saved conformal lapse")
    observer, columns, weights = inputs["observer"], inputs["columns"], inputs["weights"]
    blocks = initial_blocks(inputs)
    if not blocks["covariance_admissible"] or blocks["C_AE_frobenius"] <= 1e-8:
        raise RuntimeError("actual correlated preparation is inadmissible or inactive")
    if not 0 < pilot_budget_s <= 300 or not pilot_budget_s <= confirmation_budget_s <= 1800:
        raise ValueError("CPU caps must fit the declared 300/1800 ceilings")
    ceiling = float(pilot_budget_s)
    phase = "pilot"
    last_notice = time.process_time()

    def monitored_hamiltonian(mark):
        nonlocal last_notice
        elapsed = time.process_time() - started
        if elapsed > ceiling:
            raise RuntimeError("declared CPU ceiling exhausted during " + phase)
        if time.process_time() - last_notice >= 40:
            print(f"local response {phase}: CPU {elapsed:.1f}/{ceiling:.0f}", flush=True)
            last_notice = time.process_time()
        return schedule(mark)

    gates = [action_gate(schedule, columns, float(inputs["times"][0])),
             action_gate(schedule, np.vstack((inputs["phi0"][-1], inputs["phi1"][-1])), float(inputs["times"][-1]))]
    if max(gates) > 1e-10:
        raise RuntimeError("conformal operator fails independent Fourier action")
    pilot_times = inputs["times"][:3]
    pilot_start = time.process_time()
    pilot = shared_stream(monitored_hamiltonian, observer, pilot_times, columns, weights)
    pilot_cpu = time.process_time() - pilot_start
    direct, direct_occ, _ = evolve_streamed(monitored_hamiltonian, observer, pilot_times, columns, weights, propagator_substeps=2)
    superposition_error = _max_abs(pilot["amplitudes"] - direct["amplitudes"])
    if superposition_error > 1e-11:
        raise RuntimeError("shared column solve does not reproduce actual unsplit state")
    steps = inputs["times"].size - 1
    baseline_forecast = pilot_cpu * (steps / 2) ** 2
    full_plan_forecast = 1.15 * (baseline_forecast + pilot_cpu * (4 + 20))
    elapsed = time.process_time() - started
    if elapsed + full_plan_forecast > ceiling:
        if elapsed + full_plan_forecast > confirmation_budget_s:
            raise RuntimeError("measured production forecast exceeds confirmation ceiling")
        ceiling = float(confirmation_budget_s)
    print(f"pilot {pilot_cpu:.6f}s; shared full forecast {baseline_forecast:.3f}s; controls plan {full_plan_forecast:.3f}s; ceiling {ceiling:.0f}s", flush=True)
    occupation_auto, coherence_auto = autonomous_series(observer, inputs["phi0"], inputs["phi1"], weights)
    occupation_effect = _max_abs(occupation_auto - occupation_auto[0])
    coherence_effect = _max_abs(coherence_auto - coherence_auto[0])
    record = {
        "schema": SCHEMA, "status": "RUNNING", "case": PRIMARY, "gauge": "conformal", "time_window": list(WINDOW),
        "source_bindings": inputs["source_bindings"], "initial_blocks": blocks,
        "observer": "original T0 columns 0 and 1, no QR or phase reset", "weights": weights.tolist(),
        "initial_state": "actual Phi(0.2) on the bitwise continued conformal episode", "observer_id": array_id(observer),
        "initial_state_id": array_id(columns), "geometry_id": array_id(inputs["Q"]),
        "operator": "sigma2 P+kappa Q sigma1; original representative block, no isotropic copy factor",
        "geometry": "stored coarse Q prolonged to original quadrature; piecewise linear between stored frames",
        "operator_action_relative": gates, "pilot_superposition_max_abs": superposition_error,
        "pilot_cpu_seconds": pilot_cpu, "baseline_forecast_seconds": baseline_forecast,
        "full_plan_forecast_seconds": full_plan_forecast, "pilot_cpu_ceiling": pilot_budget_s,
        "confirmation_cpu_ceiling": confirmation_budget_s, "confirmation_admitted_by_forecast": ceiling > pilot_budget_s,
        "occupation_effect": occupation_effect, "coherence_effect": coherence_effect,
        "physical_scope": physical_scope(inputs), "production_geometry_rerun": False,
        "source_code_sha256": sha256_file(__file__),
        "frozen_feasibility_module_sha256": sha256_file(LAB / "src/recursive_horizons/nsc_conformal_local_response.py"),
        "reducer_sha256": sha256_file(LAB / "src/recursive_horizons/nsc_evolving_reduction.py"),
        "consumer_helpers_sha256": sha256_file(LAB / "src/recursive_horizons/nsc_coupled_local_response.py"),
        "continuum_constraint_certified": False, "state_error_certified": False, "stress_claimed": False, "renewal_claimed": False,
        "state_dependent_duhamel": "e<=e0+integral||delta H Phi W^1/2||F; fixed-observer covariance error<=2a e+e^2; unknown geometric defect not bounded here",
        "refinement_trials": [],
    }
    arrays = {"autonomous_times": inputs["times"], "occupation_autonomous": occupation_auto, "coherence_autonomous": coherence_auto}
    _save(record, arrays)
    for refinement in range(3):
        times = np.linspace(*WINDOW, steps * 2 ** refinement + 1)
        if time.process_time() - started + full_plan_forecast * 4 ** refinement > ceiling:
            raise RuntimeError("next complete refinement does not fit the CPU forecast")
        phase = f"shared {times.size - 1}-step solve"
        print(phase, flush=True)
        solve_started = time.process_time()
        shared = shared_stream(monitored_hamiltonian, observer, times, columns, weights)
        solve_cpu = time.process_time() - solve_started
        phase = "independent full split"
        independent = full_split(monitored_hamiltonian, observer, times, columns, weights, substeps=2)
        # Full coherence of the actual coherent sum, independently propagated.
        occupation_full, coherence_full = full_local(monitored_hamiltonian, inputs, times, substeps=2)
        phase = "memory omission and its temporal refinement"
        omitted, memory_off, _ = evolve_streamed(monitored_hamiltonian, observer, times, columns, weights, propagator_substeps=2, memory=False)
        fine_times = np.linspace(*WINDOW, 2 * (times.size - 1) + 1)
        _, memory_fine, _ = evolve_streamed(monitored_hamiltonian, observer, fine_times, columns, weights, propagator_substeps=2, memory=False)
        phase = "full midpoint refinement"
        full_refined, coherence_refined = full_local(monitored_hamiltonian, inputs, times, substeps=4)
        phase = "rate-Hermite geometry indicator"
        rate_hamiltonian = fine_geometry_hamiltonian(schedule, schedule.fine_Q, inputs["slopes"])
        rate_occ, rate_coh = full_local(rate_hamiltonian, inputs, times, substeps=2)
        drive_effect = _max_abs(shared["occupation_drive_off"] - occupation_full)
        cross_effect = _max_abs(_diagonal(shared["cross_covariance"]))
        cross_matrix_effect = float(np.max(np.linalg.norm(shared["cross_covariance"], axis=(1, 2))))
        memory_effect = _max_abs(memory_off - occupation_full)
        frame_indices = np.arange(inputs["times"].size) * 2 ** refinement
        metrics = {
            "reduction_occupation": _fraction(_max_abs(shared["occupation"] - occupation_full), occupation_effect),
            "reduction_coherence": _fraction(_max_abs(shared["coherence"] - coherence_full), coherence_effect),
            "conditional_versus_autonomous_occupation": _fraction(_max_abs(occupation_full[frame_indices] - occupation_auto), occupation_effect),
            "conditional_versus_autonomous_coherence": _fraction(_max_abs(coherence_full[frame_indices] - coherence_auto), coherence_effect),
            "drive_off_against_independent_full": _fraction(_max_abs(shared["occupation_drive_off"] - independent["occupation_drive_off"]), drive_effect),
            "cross_against_independent_full": _fraction(_max_abs(_diagonal(shared["cross_covariance"] - independent["cross_covariance"])), cross_effect),
            "cross_matrix_against_independent_full": _fraction(float(np.max(np.linalg.norm(shared["cross_covariance"] - independent["cross_covariance"], axis=(1, 2)))), cross_matrix_effect),
            "memory_off_temporal_indicator": _fraction(_max_abs(memory_off - memory_fine[::2]), memory_effect),
            "full_midpoint_refinement_occupation": _fraction(_max_abs(full_refined - occupation_full), occupation_effect),
            "full_midpoint_refinement_coherence": _fraction(_max_abs(coherence_refined - coherence_full), coherence_effect),
            "rate_Hermite_occupation": _fraction(_max_abs(rate_occ - occupation_full), occupation_effect),
            "rate_Hermite_coherence": _fraction(_max_abs(rate_coh - coherence_full), coherence_effect),
        }
        passed = all(row["within_one_percent"] for row in metrics.values())
        record["refinement_trials"].append({"output_steps": times.size - 1, "shared_solve_cpu_seconds": solve_cpu,
                                            "metrics": metrics, "all_within_one_percent": passed})
        record.update(occupation_error=signal_error(shared["occupation"], occupation_full),
            coherence_error=signal_error(shared["coherence"], coherence_full), phase_error=phase_error(shared["coherence"], coherence_full),
            conditional_versus_autonomous=signal_error(occupation_full[frame_indices], occupation_auto),
            numerical_controls=metrics,
            controls={"memory": {"occupation_movement": memory_effect, "history_norm": float(np.linalg.norm(omitted["history"]))},
                      "outside_drive": {"occupation_movement": drive_effect, "initial_exterior_norm": blocks["exterior_column_norm"]},
                      "initial_cross": {"occupation_movement": cross_effect, "matrix_frobenius_movement": cross_matrix_effect,
                                        "actual_initial_C_AE": blocks["C_AE_frobenius"], "synthetic_cross_substituted": False}},
            allocation=allocation_report(shared["result"]),
            original_six_column_history_shape=list(shared["original_column_history"].shape),
            original_six_column_history_norm=float(np.linalg.norm(shared["original_column_history"])),
            trapezoid_residual_max=float(shared["result"]["trapezoid_residual_max"]))
        arrays.update(times=times, occupation_retained=shared["occupation"], coherence_retained=shared["coherence"],
            occupation_full=occupation_full, coherence_full=coherence_full,
            occupation_memory_off=memory_off, memory_refined_times=fine_times, occupation_memory_off_refined=memory_fine,
            occupation_drive_off=shared["occupation_drive_off"], occupation_drive_off_full=independent["occupation_drive_off"],
            covariance_retained=shared["covariance"], covariance_without_cross=shared["covariance_without_cross"],
            cross_covariance_full=independent["cross_covariance"], occupation_full_midpoint_refined=full_refined,
            occupation_rate_Hermite=rate_occ, coherence_rate_Hermite=rate_coh)
        record["cpu_seconds"] = time.process_time() - started
        _save(record, arrays)
        if passed:
            break
    phase = "space and timestep geometry comparisons"
    geometry_comparisons = {}
    state_comparisons = {}
    for case in ("nf512_dt_0.0010", "nf256_dt_0.0005"):
        other = load_inputs(case)
        fine_other = other["Q"] @ periodic_interpolation(other["nf"] - 1, inputs["nq"], inputs["length"]).T
        alternate = fine_geometry_hamiltonian(schedule, fine_other)
        occupation, coherence = full_local(alternate, inputs, times, substeps=2)
        label = "time" if other["nf"] == inputs["nf"] else "space"
        geometry_comparisons[label] = {
            "case": case, "Q_gap_on_common_quadrature": _max_abs(fine_other - schedule.fine_Q),
            "occupation": _fraction(_max_abs(occupation - occupation_full), occupation_effect),
            "coherence": _fraction(_max_abs(coherence - coherence_full), coherence_effect),
            "initial_state_and_observer": "same primary Phi(0.2) and T0 observer; only geometry schedule differs"}
        other_occ, other_coh = autonomous_series(other["observer"], other["phi0"], other["phi1"], other["weights"])
        state_comparisons[label] = {"occupation": _fraction(_max_abs(other_occ - occupation_auto), occupation_effect),
                                   "coherence": _fraction(_max_abs(other_coh - coherence_auto), coherence_effect)}
        arrays["occupation_geometry_" + label] = occupation
        arrays["coherence_geometry_" + label] = coherence
    record["geometry_comparisons"] = geometry_comparisons
    record["autonomous_refinement_comparisons"] = state_comparisons
    record["all_numerical_movements_within_one_percent"] = (
        all(row["within_one_percent"] for row in record["numerical_controls"].values())
        and all(item[key]["within_one_percent"] for item in geometry_comparisons.values() for key in ("occupation", "coherence"))
        and all(item[key]["within_one_percent"] for item in state_comparisons.values() for key in ("occupation", "coherence")))
    record["status"] = "MEASURED_SAME_REALIZATION_LOCAL_RESPONSE" if record["all_numerical_movements_within_one_percent"] else "NUMERICAL_CONTROL_UNRESOLVED"
    record["scope"] = "maintained region-0 structure and resolved energy surfaces on the same saved realization; conditional fixed-observer occupation and causal omissions; no unknown-geometry state-error certificate"
    record["occupation_initial"] = occupation_auto[0].tolist()
    record["occupation_final"] = occupation_auto[-1].tolist()
    record["cache"] = {"dense_H_cache_limit": 2, "builds": schedule.builds, "hits": schedule.hits}
    record["cpu_seconds"] = time.process_time() - started
    if record["cpu_seconds"] > ceiling:
        raise RuntimeError("production exceeded its total admitted CPU ceiling")
    if source_bindings() != inputs["source_bindings"]:
        raise RuntimeError("source changed during production reduction")
    _save(record, arrays)
    return record


def consistency_errors(record, arrays, inputs):
    errors = []
    if record.get("schema") != SCHEMA or record.get("source_bindings") != inputs["source_bindings"]:
        errors.append("source_binding")
    for key, values in (("observer_id", inputs["observer"]), ("initial_state_id", inputs["columns"]), ("geometry_id", inputs["Q"])):
        if record.get(key) != array_id(values):
            errors.append(key)
    if record.get("production_geometry_rerun") is not False or record.get("state_error_certified") is not False:
        errors.append("domain_claim")
    if record.get("allocation", {}).get("time_indexed_exterior_propagator_bytes") != 0:
        errors.append("dense_W")
    measured = _max_abs(arrays["occupation_retained"] - arrays["occupation_full"])
    if abs(measured - record["occupation_error"]["max_abs"]) > 1e-12:
        errors.append("reduction_error")
    cross = arrays["covariance_retained"] - arrays["covariance_without_cross"]
    if abs(_max_abs(_diagonal(cross)) - record["controls"]["initial_cross"]["occupation_movement"]) > 1e-12:
        errors.append("cross_movement")
    for name, array_name in (("memory", "occupation_memory_off"), ("outside_drive", "occupation_drive_off")):
        if abs(_max_abs(arrays[array_name] - arrays["occupation_full"]) - record["controls"][name]["occupation_movement"]) > 1e-12:
            errors.append(name + "_movement")
    return errors


def verify():
    before = source_bindings()
    record = json.loads(OUTPUT_JSON.read_text())
    with np.load(OUTPUT_NPZ, allow_pickle=False) as stored:
        arrays = {name: np.array(stored[name]) for name in stored.files}
    errors = consistency_errors(record, arrays, load_inputs())
    if record.get("payload_sha256") != sha256_file(OUTPUT_NPZ):
        errors.append("payload_hash")
    if record.get("source_code_sha256") != sha256_file(__file__):
        errors.append("producer_hash")
    if source_bindings() != before:
        errors.append("readonly")
    if errors:
        raise RuntimeError("production check failed: " + ", ".join(errors))
    return record
