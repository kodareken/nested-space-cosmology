"""Conditional column reduction of a finalized same-action conformal episode.

The old prescribed-gauge consumer and sealed evidence are unchanged. Here
H=σ₂P+κQσ₁, the observer is the original T=0 pair, and the correlated
initial state is the transported Phi(0.025). No geometry is reintegrated.
"""
from __future__ import annotations

from collections import OrderedDict
import json
from pathlib import Path
import time

import numpy as np

from .nsc_coupled_local_response import (
    allocation_report, array_id, autonomous_series, column_cross_control,
    compare_retained, evolve_full_columns, evolve_streamed, occupations_and_coherence,
    phase_error, sha256_file, signal_error,
)
from .nsc_spherical_coupling import KAPPA, OCCUPATIONS, antiperiodic_momentum
from .nsc_spherical_galerkin_coupling import (
    antiperiodic_interpolation, canonical_column_map, periodic_interpolation,
)

LAB = Path(__file__).resolve().parents[2]
DEV = LAB / "results" / "development"
EPISODE_JSON = DEV / "nsc-spherical-conformal-episode-v1.json"
EPISODE_NPZ = EPISODE_JSON.with_suffix(".npz")
V5_NPZ = DEV / "nsc-spherical-coupling-refinement-v5.npz"
OUTPUT_JSON = DEV / "nsc-conformal-local-response-v1.json"
OUTPUT_NPZ = OUTPUT_JSON.with_suffix(".npz")
PROBE_JSON = DEV / "nsc-conformal-local-response-feasibility-v1.json"
PROBE_NPZ = PROBE_JSON.with_suffix(".npz")
CASES = ("nf256_dt_0.0010", "nf256_dt_0.0005", "nf512_dt_0.0010", "nf512_dt_0.0005")
PRIMARY = "nf512_dt_0.0005"
T0, T1 = .025, .05
PAYLOAD_LIMIT = 64 * 1024 * 1024


def _covariance(amplitudes, weights):
    return np.einsum("tak,k,tbk->tab", amplitudes, weights, amplitudes.conj())


def _diagonal(covariance):
    return np.diagonal(covariance, axis1=1, axis2=2).real


def _max_abs(values):
    return float(np.max(np.abs(values)))


def load_final_episode(case=PRIMARY):
    """Reject partial or failed assessments before constructing a reduction."""
    if case not in CASES:
        raise ValueError("case must name a recorded conformal run")
    before = {path.name: sha256_file(path) for path in (EPISODE_JSON, EPISODE_NPZ, V5_NPZ)}
    record = json.loads(EPISODE_JSON.read_text())
    readiness_errors = episode_readiness_errors(record)
    if readiness_errors:
        raise RuntimeError("conformal episode reduction held: " + ", ".join(readiness_errors))
    if before[EPISODE_NPZ.name] != record["payload_sha256"]:
        raise RuntimeError("episode JSON does not bind the complete NPZ")
    if before[V5_NPZ.name] != record["source_hashes_before"]["v5_npz"]:
        raise RuntimeError("original v5 preparation changed")
    nf = int(case.split("_")[0][2:])
    prefix = "" if nf == 256 else "nf512_"
    with np.load(V5_NPZ, allow_pickle=False) as original, np.load(EPISODE_NPZ, allow_pickle=False) as payload:
        weights = np.array(original[prefix + "occupations"], copy=True)
        length = float(original[prefix + "length"])
        nq = int(original[prefix + "nq"])
        for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
            original_key = ("columns_" if name.startswith("phi") else "geometry_") + name
            if not np.array_equal(payload[f"nf{nf}_initial_{name}"], original[prefix + original_key]):
                raise RuntimeError("initial conformal data are not the unchanged v5 source")
        original_columns = np.vstack((payload[f"nf{nf}_initial_phi0"], payload[f"nf{nf}_initial_phi1"]))
        observer = np.array(original_columns[:, :2], copy=True)
        times = np.array(payload[case + "_frame_times"], copy=True)
        selected = (times >= T0 - 1e-12) & (times <= T1 + 1e-12)
        times = times[selected]
        phi0 = np.array(payload[case + "_frames_phi0"][selected], copy=True)
        phi1 = np.array(payload[case + "_frames_phi1"][selected], copy=True)
        geometry = np.array(payload[case + "_frames_Q"][selected], copy=True)
        columns = np.vstack((phi0[0], phi1[0]))
    if not np.array_equal(weights, OCCUPATIONS):
        raise RuntimeError("the original Gaussian weights changed")
    if times.size != 6 or abs(times[0] - T0) > 1e-12 or abs(times[-1] - T1) > 1e-12:
        raise RuntimeError("the actual correlated segment is unavailable")
    if np.array_equal(observer, columns[:, :2]):
        raise RuntimeError("observer reset to the transported state")
    after = {path.name: sha256_file(path) for path in (EPISODE_JSON, EPISODE_NPZ, V5_NPZ)}
    if before != after:
        raise RuntimeError("episode changed while the consumer read it")
    return {
        "case": case, "nf": nf, "nq": nq, "length": length, "weights": weights,
        "observer": observer, "columns": columns, "times": times, "Q": geometry,
        "phi0": phi0, "phi1": phi1, "source_bindings": before,
        "episode_assessment": "measured finite constraint forcing; sampled integral, no continuum enclosure",
        "regional_exchange_resolved": False,
    }


def episode_readiness_errors(record):
    errors = []
    if record.get("schema") != "NSC-SPHERICAL-CONFORMAL-EPISODE-v1" or record.get("gauge") != "conformal":
        errors.append("schema_or_gauge")
    if record.get("verdict") != "MEASURED_SAME_ACTION_CONSTRAINT_FORCING" or record.get("full_target_reached") is not True:
        errors.append("unfinished_episode")
    if record.get("source_changed") is not False or record.get("initial_reset") is not False:
        errors.append("changed_preparation")
    if record.get("bound_sources_unchanged_during_run") is not True:
        errors.append("source_binding")
    for name in CASES:
        run = record.get("results", {}).get(name, {})
        if run.get("completed") is not True or run.get("positive_chart") is not True:
            errors.append("chart_or_completion:" + name)
        if run.get("constraint_end_margin", -1) <= 0 or run.get("max_sampled_budget_violation", 1) > 1e-10:
            errors.append("sampled_constraint_budget:" + name)
        if run.get("gram_gap_max", 1) > 1e-6 or run.get("sbp_identity_error_max", 1) > 1e-10:
            errors.append("numerical_admissibility:" + name)
    return errors


class ConformalSchedule:
    """Canonical Galerkin H(Q(t)); at most two dense H samples are cached."""
    def __init__(self, inputs):
        self.nf, self.nq, self.length = inputs["nf"], inputs["nq"], inputs["length"]
        self.times = inputs["times"]
        self.geometry_map = periodic_interpolation(self.nf - 1, self.nq, self.length)
        self.columns_map = canonical_column_map(
            antiperiodic_interpolation(self.nf, self.nq, self.length), self.nf, self.nq,
        )
        self.fine_Q = inputs["Q"] @ self.geometry_map.T
        if not np.isfinite(self.fine_Q).all() or np.min(self.fine_Q) <= 0:
            raise ValueError("stored conformal chart is not positive on quadrature")
        self.momentum, _metric = antiperiodic_momentum(self.nf, self.length)
        self.mass = np.stack([
            KAPPA * (self.columns_map.conj().T @ (q[:, None] * self.columns_map))
            for q in self.fine_Q
        ])
        self.mass = .5 * (self.mass + self.mass.conj().transpose(0, 2, 1))
        self.cache = OrderedDict()
        self.builds = self.hits = 0

    def radial(self, mark):
        index, fraction = self._interval(mark)
        return (1 - fraction) * self.fine_Q[index] + fraction * self.fine_Q[index + 1]

    def _interval(self, mark):
        mark = float(mark)
        if mark < self.times[0] - 1e-12 or mark > self.times[-1] + 1e-12:
            raise ValueError("query leaves the saved conformal segment")
        index = min(max(int(np.searchsorted(self.times, mark, side="right") - 1), 0), self.times.size - 2)
        fraction = np.clip((mark - self.times[index]) / (self.times[index + 1] - self.times[index]), 0, 1)
        return index, fraction

    def __call__(self, mark):
        key = float(mark).hex()
        if key in self.cache:
            self.hits += 1
            self.cache.move_to_end(key)
            return self.cache[key]
        index, fraction = self._interval(mark)
        mass = (1 - fraction) * self.mass[index] + fraction * self.mass[index + 1]
        zero = np.zeros_like(mass)
        matrix = np.block([[zero, mass - 1j * self.momentum], [mass + 1j * self.momentum, zero]])
        self.cache[key] = matrix
        if len(self.cache) > 2:
            self.cache.popitem(last=False)
        self.builds += 1
        return matrix


def action_gate(schedule, columns, mark):
    """Independent fine-grid AP Fourier action, then the canonical adjoint."""
    n = schedule.nq
    phase = np.exp(1j * np.pi * np.arange(n) / n)[:, None]
    wave = (np.fft.fftfreq(n) * n + .5) * 2 * np.pi / schedule.length

    def momentum(values):
        return phase * np.fft.ifft(wave[:, None] * np.fft.fft(values / phase, axis=0), axis=0)

    phi0 = schedule.columns_map @ columns[:schedule.nf]
    phi1 = schedule.columns_map @ columns[schedule.nf:]
    q = schedule.radial(mark)[:, None]
    image = np.vstack((
        schedule.columns_map.conj().T @ (-1j * momentum(phi1) + KAPPA * q * phi1),
        schedule.columns_map.conj().T @ (1j * momentum(phi0) + KAPPA * q * phi0),
    ))
    return float(np.linalg.norm(schedule(mark) @ columns - image) / np.linalg.norm(image))


def initial_blocks(inputs):
    observer, columns, weights = inputs["observer"], inputs["columns"], inputs["weights"]
    coefficients = observer.conj().T @ columns
    exterior = columns - observer @ coefficients
    parent = (coefficients * weights) @ coefficients.conj().T
    cross = (coefficients * weights) @ exterior.conj().T
    weighted = columns * np.sqrt(weights)
    spectrum = np.linalg.eigvalsh(weighted.conj().T @ weighted).real
    return {
        "C_AA_real": parent.real.tolist(), "C_AA_imag": parent.imag.tolist(),
        "C_AA_eigenvalues": np.linalg.eigvalsh(parent).real.tolist(),
        "C_AE_frobenius": float(np.linalg.norm(cross)),
        "exterior_column_norm": float(np.linalg.norm(exterior)),
        "covariance_eigenvalues": spectrum.tolist(),
        "covariance_admissible": bool(np.min(spectrum) >= -1e-12 and np.max(spectrum) <= 1 + 1e-12),
        "observer_gram_gap": _max_abs(observer.conj().T @ observer - np.eye(2)),
    }


def full_split(schedule, observer, times, columns, weights, substeps=2):
    a0 = observer @ (observer.conj().T @ columns)
    e0 = columns - a0
    propagated = evolve_full_columns(schedule, times, np.concatenate((a0, e0), axis=1), substeps=substeps)
    projected = np.einsum("ij,tjk->tik", observer.conj().T, propagated)
    a, e = projected[:, :, :columns.shape[1]], projected[:, :, columns.shape[1]:]
    total = _covariance(a + e, weights)
    dropped = _covariance(a, weights) + _covariance(e, weights)
    return {
        "cross_covariance": total - dropped,
        "occupation_drive_off": _diagonal(_covariance(a, weights)),
        "covariance_without_cross": dropped, "occupation_full": _diagonal(total),
    }


def _timed(callback):
    start = time.process_time()
    result = callback()
    return result, time.process_time() - start


def _fraction(error, effect):
    return {"movement": float(error), "effect": float(effect), "fraction": None if effect <= 0 else float(error / effect),
            "within_one_percent": bool(effect > 0 and error <= .01 * effect)}


def execute(budget_s=300.0, inputs=None):
    started = time.process_time()
    inputs = load_final_episode() if inputs is None else inputs
    if inputs.get("regional_exchange_resolved") is not True:
        raise RuntimeError("production reduction held: the initial conformal window is feasibility only; a finalized resolved-exchange successor is required")
    schedule = ConformalSchedule(inputs)
    observer, columns, weights = inputs["observer"], inputs["columns"], inputs["weights"]
    blocks = initial_blocks(inputs)
    if not blocks["covariance_admissible"] or blocks["C_AE_frobenius"] <= 1e-8:
        raise RuntimeError("actual correlated initial state is inadmissible or inactive")
    gate0 = action_gate(schedule, columns, T0)
    final_columns = np.vstack((inputs["phi0"][-1], inputs["phi1"][-1]))
    gate1 = action_gate(schedule, final_columns, T1)
    if max(gate0, gate1) > 1e-10:
        raise RuntimeError("conformal Hamiltonian does not reproduce the independent Fourier action")
    auto_occupation, auto_coherence = autonomous_series(observer, inputs["phi0"], inputs["phi1"], weights)
    effect = _max_abs(auto_occupation - auto_occupation[0])
    pilot, pilot_cpu = _timed(lambda: compare_retained(schedule, observer, inputs["times"][:3], columns, weights, substeps=2))
    steps = inputs["times"].size - 1
    forecast = pilot_cpu * (steps / 2) ** 2
    # One baseline, two omission runs, two cross-split runs, one refined
    # memoryless run and a full split. Attribute all cost quadratically.
    full_plan_forecast = 12.0 * forecast
    if time.process_time() - started + full_plan_forecast > budget_s:
        raise RuntimeError("measured full-control forecast exceeds the CPU ceiling")
    print(f"conformal response pilot {pilot_cpu:.6f}s; full plan forecast {full_plan_forecast:.6f}s", flush=True)
    trials = []
    final = None
    for refinement in range(3):
        times = np.linspace(T0, T1, steps * 2 ** refinement + 1)
        expected = full_plan_forecast * 4 ** refinement
        if time.process_time() - started + expected > budget_s:
            raise RuntimeError("refined control forecast exceeds the CPU ceiling")
        comparison, comparison_cpu = _timed(lambda: compare_retained(schedule, observer, times, columns, weights, substeps=2))
        omissions = {}
        controls = {}
        for name, flags in (("memory", {"memory": False}), ("outside_drive", {"outside_drive": False})):
            result, occupation, coherence = evolve_streamed(schedule, observer, times, columns, weights, propagator_substeps=2, **flags)
            controls[name] = occupation
            omissions[name] = {"movement": _max_abs(occupation - comparison["occupation_full"]),
                               "history_norm": float(np.linalg.norm(result["history"])),
                               "initial_exterior_norm": float(np.linalg.norm(result["initial_exterior"]))}
        split = column_cross_control(schedule, observer, times, columns, weights, substeps=2, baseline=comparison["result"])
        independent = full_split(schedule, observer, times, columns, weights)
        finer_times = np.linspace(T0, T1, 2 * (times.size - 1) + 1)
        _, memory_fine, _ = evolve_streamed(schedule, observer, finer_times, columns, weights, propagator_substeps=2, memory=False)
        full_fine = evolve_full_columns(schedule, times, columns, substeps=4)
        fine_projected = np.einsum("ij,tjk->tik", observer.conj().T, full_fine)
        full_refined_occupation, _ = occupations_and_coherence(fine_projected, weights)
        numerical = {
            "reduction": _fraction(comparison["occupation_error"]["max_abs"], effect),
            "full_midpoint_refinement": _fraction(_max_abs(full_refined_occupation - comparison["occupation_full"]), effect),
            "drive_off": _fraction(_max_abs(controls["outside_drive"] - independent["occupation_drive_off"]), omissions["outside_drive"]["movement"]),
            "memory_off_refinement_indicator": _fraction(_max_abs(memory_fine[::2] - controls["memory"]), omissions["memory"]["movement"]),
            "cross": _fraction(_max_abs(_diagonal(split["cross_covariance"] - independent["cross_covariance"])), split["occupation_diagonal_movement"]),
        }
        indices = np.arange(inputs["times"].size) * 2 ** refinement
        conditional_gap = signal_error(comparison["occupation_full"][indices], auto_occupation)
        numerical["conditional_versus_autonomous"] = _fraction(conditional_gap["max_abs"], effect)
        omissions["initial_cross"] = {"movement": split["occupation_diagonal_movement"],
            "separation_max_frobenius": split["separation_max_frobenius"],
            "algebraic_superposition_error": split["numerical_error"],
            "synthetic_cross_substituted": False}
        passed = all(row["within_one_percent"] for row in numerical.values())
        trials.append({"output_steps": times.size - 1, "comparison_cpu": comparison_cpu,
                       "numerical": numerical, "all_within_one_percent": passed})
        final = (times, comparison, controls, omissions, split, independent, numerical, indices, memory_fine, full_refined_occupation)
        if passed:
            break
    times, comparison, controls, omissions, split, independent, numerical, indices, memory_fine, full_refined_occupation = final
    passed = all(row["within_one_percent"] for row in numerical.values())
    partners = {}
    for case in CASES:
        other = load_final_episode(case)
        other_occ, other_coh = autonomous_series(other["observer"], other["phi0"], other["phi1"], other["weights"])
        partners[case] = {"occupation": _fraction(_max_abs(other_occ - auto_occupation), effect),
                          "coherence_max_gap": _max_abs(other_coh - auto_coherence)}
    arrays = {
        "times": times, "autonomous_times": inputs["times"],
        "occupation_retained": comparison["occupation_retained"], "occupation_full": comparison["occupation_full"],
        "occupation_autonomous": auto_occupation, "coherence_retained": comparison["coherence_retained"],
        "coherence_full": comparison["coherence_full"], "coherence_autonomous": auto_coherence,
        "occupation_memory_off": controls["memory"], "occupation_drive_off": controls["outside_drive"],
        "covariance_retained": split["covariance_total"], "covariance_without_cross": split["covariance_without_cross"],
        "cross_covariance_full": independent["cross_covariance"],
        "occupation_drive_off_full": independent["occupation_drive_off"],
        "occupation_memory_off_refined": memory_fine, "occupation_full_midpoint_refined": full_refined_occupation,
    }
    allocation = allocation_report(comparison["result"])
    if allocation["time_indexed_exterior_propagator_bytes"] != 0:
        raise RuntimeError("conformal reduction stored a dense propagator history")
    record = {
        "schema": "NSC-CONFORMAL-LOCAL-RESPONSE-v1",
        "status": "MEASURED_SAME_REALIZATION_LOCAL_RESPONSE" if passed else "NUMERICAL_CONTROL_NOT_RESOLVED",
        "case": PRIMARY, "gauge": "conformal", "time_window": [T0, T1], "substeps": 2,
        "observer": "original T=0 region-0 columns 0 and 1, no QR or phase reset",
        "initial_state": "actual transported Phi(0.025) on the generated conformal episode",
        "weights": weights.tolist(), "initial_blocks": blocks, "source_bindings": inputs["source_bindings"],
        "observer_id": array_id(observer), "initial_state_id": array_id(columns), "geometry_id": array_id(inputs["Q"]),
        "geometry": "piecewise-linear stored coarse Q, prolonged to the recorded quadrature; L=Q, beta=0",
        "hamiltonian": "sigma2 P + kappa Q sigma1; multiplicity only in the source energy/force",
        "episode_assessment": inputs["episode_assessment"], "action_gate": [gate0, gate1],
        "occupation_change": effect, "occupation_initial": auto_occupation[0].tolist(), "occupation_final": auto_occupation[-1].tolist(),
        "occupation_error": comparison["occupation_error"], "coherence_error": comparison["coherence_error"],
        "phase_error": phase_error(comparison["coherence_retained"], comparison["coherence_full"]),
        "conditional_versus_autonomous": signal_error(comparison["occupation_full"][indices], auto_occupation),
        "omissions": omissions, "own_numerical_controls": numerical, "trials": trials, "partners": partners,
        "method_error_within_one_percent": passed and all(row["occupation"]["within_one_percent"] for row in partners.values()),
        "history_norm": comparison["history_norm"], "initial_exterior_norm": comparison["initial_exterior_norm"],
        "trapezoid_residual_max": comparison["trapezoid_residual_max"], "allocation": allocation,
        "pilot_cpu_seconds": pilot_cpu, "full_plan_forecast_seconds": full_plan_forecast,
        "cpu_seconds": time.process_time() - started, "cpu_budget_seconds": budget_s,
        "cache": {"limit": 2, "builds": schedule.builds, "hits": schedule.hits, "dense_H_history_stored": False},
        "state_dependent_duhamel": "e<=e0+integral||delta H Phi W^1/2||_F; fixed-observer covariance error<=2a e+e^2",
        "geometry_state_error_bound": None, "conditional_geometry_only": True,
        "continuum_constraint_certified": False, "stress_claimed": False, "renewal_claimed": False,
        "production_geometry_rerun": False,
        "source_code_sha256": sha256_file(__file__),
        "reducer_sha256": sha256_file(LAB / "src/recursive_horizons/nsc_evolving_reduction.py"),
        "consumer_helpers_sha256": sha256_file(LAB / "src/recursive_horizons/nsc_coupled_local_response.py"),
    }
    if time.process_time() - started > budget_s:
        raise RuntimeError("conformal reduction exceeded the CPU ceiling")
    if sum(value.nbytes for value in arrays.values()) > PAYLOAD_LIMIT:
        raise RuntimeError("uncompressed local payload exceeds 64 MiB")
    np.savez_compressed(OUTPUT_NPZ, **arrays)
    record["payload_sha256"] = sha256_file(OUTPUT_NPZ)
    OUTPUT_JSON.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    for _ in range(3):
        record["payload_bytes"] = OUTPUT_JSON.stat().st_size + OUTPUT_NPZ.stat().st_size
        OUTPUT_JSON.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    if record["payload_bytes"] > PAYLOAD_LIMIT:
        raise RuntimeError("local payload exceeds 64 MiB")
    if {path.name: sha256_file(path) for path in (EPISODE_JSON, EPISODE_NPZ, V5_NPZ)} != inputs["source_bindings"]:
        raise RuntimeError("source episode changed during the local response")
    return record


def probe(budget_s=60.0):
    """Only two stored steps of the initial episode; no programme verdict."""
    started = time.process_time()
    inputs = load_final_episode()
    schedule = ConformalSchedule(inputs)
    observer, columns, weights = inputs["observer"], inputs["columns"], inputs["weights"]
    times = inputs["times"][:3]
    gate = action_gate(schedule, columns, float(times[0]))
    if gate > 1e-10:
        raise RuntimeError("conformal operator action failed the feasibility gate")
    comparison, pilot_cpu = _timed(lambda: compare_retained(schedule, observer, times, columns, weights, substeps=2))
    available_steps = inputs["times"].size - 1
    forecast = pilot_cpu * (available_steps / 2) ** 2
    if time.process_time() - started + 4 * pilot_cpu > budget_s:
        raise RuntimeError("probe controls forecast exceeds the feasibility budget")
    print(f"feasibility probe {pilot_cpu:.6f}s; five-step comparison forecast {forecast:.6f}s (not run)", flush=True)
    controls, arrays = {}, {"times": times, "occupation_full": comparison["occupation_full"],
        "occupation_retained": comparison["occupation_retained"], "coherence_full": comparison["coherence_full"],
        "coherence_retained": comparison["coherence_retained"]}
    for name, flags in (("memory", {"memory": False}), ("outside_drive", {"outside_drive": False})):
        result, occupation, _ = evolve_streamed(schedule, observer, times, columns, weights, propagator_substeps=2, **flags)
        controls[name] = {"movement": _max_abs(occupation - comparison["occupation_full"]),
                          "history_norm": float(np.linalg.norm(result["history"]))}
        arrays["occupation_" + name + "_off"] = occupation
    split = column_cross_control(schedule, observer, times, columns, weights, substeps=2, baseline=comparison["result"])
    independent = full_split(schedule, observer, times, columns, weights)
    arrays["covariance_retained"] = split["covariance_total"]
    arrays["covariance_without_cross"] = split["covariance_without_cross"]
    arrays["cross_covariance_full"] = independent["cross_covariance"]
    arrays["occupation_drive_off_full"] = independent["occupation_drive_off"]
    controls["initial_cross"] = {"movement": split["occupation_diagonal_movement"],
        "occupation_error_against_full_split": _max_abs(_diagonal(split["cross_covariance"] - independent["cross_covariance"])),
        "algebraic_superposition_error": split["numerical_error"], "synthetic_cross_substituted": False}
    controls["outside_drive"]["occupation_error_against_full_drive_off"] = _max_abs(
        arrays["occupation_outside_drive_off"] - independent["occupation_drive_off"])
    auto_occupation, _ = autonomous_series(observer, inputs["phi0"][:3], inputs["phi1"][:3], weights)
    arrays["occupation_autonomous"] = auto_occupation
    effect = _max_abs(auto_occupation - auto_occupation[0])
    record = {
        "schema": "NSC-CONFORMAL-LOCAL-RESPONSE-FEASIBILITY-v1", "status": "MEASURED_REDUCTION_FEASIBILITY",
        "case": PRIMARY, "gauge": "conformal", "time_window": times[[0, -1]].tolist(), "substeps": 2,
        "observer": "original T=0 region-0 pair; no phase or basis reset", "weights": weights.tolist(),
        "initial_state": "actual transported Phi(0.025)", "initial_blocks": initial_blocks(inputs),
        "occupation_change": effect, "occupation_error": comparison["occupation_error"],
        "conditional_versus_autonomous": signal_error(comparison["occupation_full"], auto_occupation),
        "controls": controls, "allocation": allocation_report(comparison["result"]),
        "operator_action_relative": gate, "pilot_cpu_seconds": pilot_cpu,
        "five_step_comparison_forecast_seconds": forecast, "full_segment_comparison_executed": False,
        "production_reduction_executed": False, "regional_throughflow_resolved": False,
        "programme_complete": False, "source_bindings": inputs["source_bindings"],
        "geometry_reintegrated": False, "observer_id": array_id(observer), "initial_state_id": array_id(columns),
        "cpu_seconds": time.process_time() - started, "cpu_budget_seconds": budget_s,
        "source_code_sha256": sha256_file(__file__), "stress_claimed": False, "renewal_claimed": False,
        "scope": "initial-episode reduction feasibility only; production waits for finalized resolved regional exchange",
    }
    if time.process_time() - started > budget_s:
        raise RuntimeError("feasibility controls exceeded their CPU ceiling")
    np.savez_compressed(PROBE_NPZ, **arrays)
    record["payload_sha256"] = sha256_file(PROBE_NPZ)
    PROBE_JSON.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    for _ in range(3):
        record["payload_bytes"] = PROBE_JSON.stat().st_size + PROBE_NPZ.stat().st_size
        PROBE_JSON.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    if record["payload_bytes"] > PAYLOAD_LIMIT:
        raise RuntimeError("feasibility payload exceeds 64 MiB")
    if {path.name: sha256_file(path) for path in (EPISODE_JSON, EPISODE_NPZ, V5_NPZ)} != inputs["source_bindings"]:
        raise RuntimeError("finalized source episode changed during the probe")
    return record
