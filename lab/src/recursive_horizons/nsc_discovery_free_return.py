"""Exact free-carrier transport compared with authenticated saved T3/T8 fields.

This applies a finite matrix exponential, not a new trajectory. Geometry and
canonical pi are decoded through W before prolongation; field columns remain
in their original retained antiperiodic band. All six weights and covariance
cross terms are retained.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import time

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import numpy as np
from threadpoolctl import threadpool_limits

from . import nsc_discovery_episode as episode
from . import nsc_discovery_observables as observer
from . import nsc_discovery_tidal as tidal
from . import nsc_nested_parent_child as nested
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

SCHEMA = "NSC-DISCOVERY-FREE-RETURN-v1"
REPO = episode.REPO
LAB = episode.LAB
PREDECESSOR = LAB / "results/development/nsc-discovery-episode-v1"
CROSSING = LAB / "results/development/nsc-discovery-crossing-v1"
OUTPUT = LAB / "results/development/nsc-discovery-free-return-v1.json"
PRIMARY = "nf256_coupled_dt0.0005"
CONTROLS = ("nf256_coupled_dt0.001", "nf128_coupled_dt0.0005",
            "nf128_coupled_dt0.001", "nf256_frozen_geometry_dt0.0005")
CPU_BUDGET = 10.0
CROSSING_COMMIT = "b7f0dab8de62d9a3472ac4e7b33351ab65d25847"


def relative(path):
    return str(Path(path).resolve().relative_to(REPO.resolve()))


def bound_path(name):
    if Path(name).is_absolute():
        raise ValueError("bound paths must be repository relative")
    path = (REPO / name).resolve()
    path.relative_to(REPO.resolve())
    return path


def hashes(paths):
    return {relative(path): tidal.sha256_file(path) for path in paths}


def source_binding():
    """Saved producing envelope plus current exact-operator/consumer source hashes."""
    path = CROSSING / "observed-run-binding.json"
    crossing = json.loads(path.read_text())
    if crossing.get("producing_commit") != CROSSING_COMMIT:
        raise ValueError("unexpected crossing producing commit")
    for group in ("producer_paths", "artifact_paths"):
        entries = crossing.get(group) or {}
        if not entries:
            raise ValueError("crossing lacks producing/source bindings")
        for name, value in entries.items():
            if tidal.sha256_file(bound_path(name)) != value["sha256"]:
                raise ValueError("crossing source/input changed: "+name)
    original = tidal.frozen_binding(PREDECESSOR, "episode")
    original["envelope"] = relative(original["envelope"])
    original["physics_hashes"] = {relative(name): value for name, value in original["physics_hashes"].items()}
    sources = list(crossing["producer_paths"])
    sources += [relative(path) for path in tidal.producer_hashes()]
    sources += [relative(observer.__file__), relative(__file__),
                "lab/scripts/derive_nsc_discovery_free_return.py",
                "lab/tests/test_nsc_discovery_free_return.py", "lab/docs/nsc-discovery-free-return.md"]
    return {
        "crossing_envelope": relative(path), "crossing_envelope_sha256": tidal.sha256_file(path),
        "producing_commit": crossing["producing_commit"], "binding_timing": crossing["binding_timing"],
        "artifact_hashes": {name: value["sha256"] for name, value in crossing["artifact_paths"].items()},
        "original_episode": original,
    }, hashes(bound_path(name) for name in set(sources))


def free_backend(grid):
    """Exact retained P=U_f^H P_fine U_f, with the original half-integer AP spectrum."""
    P = grid.U_f.conj().T @ grid.momentum @ grid.U_f
    spectrum, vectors = np.linalg.eigh(P)
    expected = 2*np.pi*galerkin.fermion_modes(grid.nf)/grid.length
    return {
        "P": P, "eigenvalues": spectrum, "eigenvectors": vectors,
        "hermitian_gap": float(np.max(abs(P-P.conj().T))),
        "AP_spectrum_gap": float(np.max(abs(spectrum-expected))),
        "period": float(grid.length),
    }


def propagate(backend, columns, delta_time):
    """U0(dt)=exp(-i sigma2 P dt), acting on all columns without covariance truncation."""
    phi = np.asarray(columns, dtype=complex)
    nf = backend["P"].shape[0]
    if phi.ndim != 2 or phi.shape[0] != 2*nf or not np.isfinite(phi).all():
        raise ValueError("columns must be finite with first dimension 2nf")
    V, wave = backend["eigenvectors"], backend["eigenvalues"]
    a, b = V.conj().T @ phi[:nf], V.conj().T @ phi[nf:]
    co, si = np.cos(wave*delta_time)[:, None], np.sin(wave*delta_time)[:, None]
    return np.vstack((V @ (co*a-si*b), V @ (si*a+co*b)))


def free_image(backend, columns):
    nf = backend["P"].shape[0]
    return np.vstack((-1j*backend["P"] @ columns[nf:], 1j*backend["P"] @ columns[:nf]))


def weighted_norm(columns, weights):
    return float(np.linalg.norm(np.asarray(columns)*np.sqrt(weights)[None, :]))


def covariance(columns, weights):
    return (columns*np.asarray(weights)[None, :]) @ columns.conj().T


def fixed_window_probability(grid, columns, weights, interval=(1.0, 3.0)):
    up0, up1 = grid.U_f @ columns[:grid.nf], grid.U_f @ columns[grid.nf:]
    mass = np.sum((abs(up0)**2+abs(up1)**2)*np.asarray(weights)[None, :], axis=1)
    return observer.reference_window_integral(grid, mass, interval)


def duhamel_bound(kappa, initial_norm, q_integral_upper_bound):
    """Conditional unitary bound; the caller must supply an actual integral upper bound."""
    if any(not math.isfinite(x) or x < 0 for x in (initial_norm, q_integral_upper_bound)):
        raise ValueError("norm and integral upper bound must be finite and nonnegative")
    return abs(float(kappa))*float(initial_norm)*float(q_integral_upper_bound)


def validate_join(original, handoff, endpoint):
    if (not original.get("array_sha256")
            or original["array_sha256"] != handoff.get("array_sha256")):
        raise ValueError("T3 state continuation array dictionaries differ")
    pointer = handoff.get("predecessor_chunk") or {}
    if (pointer.get("arrays_sha256") != original["arrays_sha256"]
            or pointer.get("ordinal") != original["ordinal"]
            or original["case_id"] != handoff["case_id"]
            or pointer.get("case_id") != original["case_id"]
            or endpoint["case_id"] != original["case_id"]
            or abs(original["coordinate_time"]-3) > 1e-8
            or abs(handoff["coordinate_time"]-3) > 1e-8
            or abs(endpoint["coordinate_time"]-8) > 1e-8
            or original["source_pins"]["W_sha256"] != endpoint["source_pins"]["W_sha256"]
            or original["source_pins"]["weights_sha256"] != endpoint["source_pins"]["weights_sha256"]):
        raise ValueError("saved state continuation pointer/time/source mismatch")
    return {"exact_T3_handoff": True, "all_array_hashes_equal": True,
            "no_state_geometry_clock_basis_reset": True}


def gram_report(columns, weights):
    gram = columns.conj().T @ columns
    eigen = np.linalg.eigvalsh(np.sqrt(weights[:, None]*weights[None, :])*gram)
    return {"gram_gap": float(np.max(abs(gram-np.eye(weights.size)))),
            "occupied_covariance_eigenvalues": eigen.tolist(), "other_covariance_eigenvalues": "zero"}


def measure_case(case_id, cache):
    original, initial_arrays = episode.load_checkpoint(PREDECESSOR, case_id, 2)
    handoff, _ = episode.load_checkpoint(CROSSING, case_id, 0)
    endpoint, final_arrays = episode.load_checkpoint(CROSSING, case_id, 1)
    joined = validate_join(original, handoff, endpoint)
    pair = episode.pair_from_arrays(initial_arrays, original)
    initial = episode.state_from_arrays(initial_arrays, original["momentum_representation"])
    final_pair = episode.pair_from_arrays(final_arrays, endpoint)
    final = episode.state_from_arrays(final_arrays, endpoint["momentum_representation"])
    # Q/r/chi are W coefficients; canonical pi additionally carries dx_g.
    # Never pass those coefficients directly into prolong_state.
    nodal_initial, nodal_final = nested.reconstruct_state(pair, initial), nested.reconstruct_state(final_pair, final)
    fine_initial = galerkin.prolong_state(pair.grid, nodal_initial)
    fine_final = galerkin.prolong_state(final_pair.grid, nodal_final)
    if min(float(fine_initial.Q.min()), float(fine_final.Q.min())) <= 0:
        raise ValueError("decoded saved metric left the positive conformal chart")
    weights = np.asarray(pair.weights)
    if not np.array_equal(weights, final_pair.weights):
        raise ValueError("saved source weights changed")
    backend = cache.setdefault(pair.grid.nf, None)
    if backend is None:
        backend = cache[pair.grid.nf] = free_backend(pair.grid)
    phi0, actual = np.vstack((initial.phi0, initial.phi1)), np.vstack((final.phi0, final.phi1))
    predicted = propagate(backend, phi0, 5.0)
    period = propagate(backend, phi0, backend["period"])
    norm0, actual_norm = weighted_norm(phi0, weights), weighted_norm(actual, weights)
    diff = weighted_norm(predicted-actual, weights)
    c_predicted, c_actual = covariance(predicted, weights), covariance(actual, weights)
    cov_diff = float(np.linalg.norm(c_predicted-c_actual))
    free_image0 = free_image(backend, phi0)
    actual_image0 = nested.apply_hamiltonian(pair, initial, phi0)
    up0, up1 = pair.grid.U_f @ phi0[:pair.grid.nf], pair.grid.U_f @ phi0[pair.grid.nf:]
    mass0 = pair.grid.fine.kappa*(pair.grid.U_f.conj().T @ (fine_initial.Q[:, None]*up1))
    mass1 = pair.grid.fine.kappa*(pair.grid.U_f.conj().T @ (fine_initial.Q[:, None]*up0))
    operator_gap = actual_image0-free_image0-np.vstack((mass0, mass1))
    # Purely algebraic zero-mass check on the saved positive metric; no grid,
    # coupling coefficient, state or trajectory is modified.
    zero0, zero1 = coupling.apply_dirac(up0, up1, fine_initial.Q, fine_initial.Q,
                                      np.zeros(pair.grid.nq), 0, pair.grid.momentum)
    zero_image = np.vstack((pair.grid.U_f.conj().T @ zero0, pair.grid.U_f.conj().T @ zero1))
    p_initial = fixed_window_probability(pair.grid, phi0, weights)
    p_predicted = fixed_window_probability(pair.grid, predicted, weights)
    p_actual = fixed_window_probability(pair.grid, actual, weights)
    latest = next(row for row in episode.read_observations(CROSSING, case_id) if abs(row["time"]-8) < 1e-8)
    return {
        "case_id": case_id, "nf": pair.grid.nf, "step_cap": original["step_cap"],
        "times": [3.0, 8.0], "delta_time": 5.0, "join": joined, "source_weights": weights.tolist(),
        "momentum_geometry_decoding": "nested.reconstruct_state before galerkin.prolong_state; field columns remain in retained AP band",
        "Q_sup_samples": {"3": float(np.max(abs(fine_initial.Q))), "8": float(np.max(abs(fine_final.Q)))},
        "weighted_field_norm_T3": norm0, "weighted_field_norm_T8": actual_norm,
        "weighted_field_error": diff, "weighted_field_relative_error": diff/max(actual_norm, 1e-300),
        "source_covariance_frobenius_error": cov_diff,
        "source_covariance_relative_error": cov_diff/max(float(np.linalg.norm(c_actual)), 1e-300),
        "source_covariance_trace_difference": float(np.trace(c_predicted-c_actual).real),
        "source_cross_terms_omitted": False,
        "fixed_child_window": {"interval": [1.0, 3.0], "probability_T3": p_initial,
                               "free_probability_T8": p_predicted, "saved_probability_T8": p_actual,
                               "absolute_free_saved_gap": abs(p_predicted-p_actual),
                               "saved_ledger_probability_T8": latest["child_probability"],
                               "saved_ledger_gap": abs(p_actual-latest["child_probability"])},
        "gram_CAR": {"initial": gram_report(phi0, weights), "free": gram_report(predicted, weights),
                     "saved": gram_report(actual, weights)},
        "backend": {"operator": "H0=sigma2 P; P=U_f^H P_fine U_f",
                    "method": "Hermitian eigensystem, exact finite matrix exponential",
                    "AP_modes": "half-integers", "period": backend["period"],
                    "hermitian_gap": backend["hermitian_gap"], "AP_spectrum_gap": backend["AP_spectrum_gap"],
                    "representative_multiplicity_applied": False},
        "controls": {"zero_duration_max_gap": float(np.max(abs(propagate(backend, phi0, 0)-phi0))),
                     "U0_period_plus_identity_weighted_gap": weighted_norm(period+phi0, weights),
                     "period_probability_gap": abs(fixed_window_probability(pair.grid, period, weights)-p_initial),
                     "actual_H_minus_free_minus_mass_max_gap": float(np.max(abs(operator_gap))),
                     "zero_mass_operator_max_gap": float(np.max(abs(zero_image-free_image0))),
                     "free_Gram_change_max": float(np.max(abs(predicted.conj().T@predicted-phi0.conj().T@phi0)))},
        "Duhamel": {"formula": "||Phi(T8)-U0(5)Phi(T3)||_c <= |kappa| ||Phi(T3)||_c integral_T3^T8 ||Q(t)||_infinity dt",
                    "coefficient": abs(pair.grid.fine.kappa)*norm0, "Gronwall_factor": False,
                    "integral_upper_bound": None, "certified_error_upper_bound": None,
                    "future_Q_history_known": False, "endpoint_samples_are_integral_bound": False,
                    "domain": "conditional exact semidiscrete field transport along the realized metric; saved RK4 drift is separate"},
    }


def measure(include_controls=True):
    started = time.process_time()
    binding, source_before = source_binding()
    cases = (PRIMARY,) + CONTROLS if include_controls else (PRIMARY,)
    paths = [CROSSING/"observed-run-binding.json", PREDECESSOR/"run-binding-envelope.json"]
    paths += [bound_path(name) for name in binding["artifact_hashes"]]
    for case_id in cases:
        paths += [PREDECESSOR/f"{case_id}-000002.{extension}" for extension in ("json", "npz")]
        paths += [CROSSING/f"{case_id}-{ordinal:06d}.{extension}" for ordinal in (0, 1) for extension in ("json", "npz")]
        paths += [CROSSING/f"{case_id}-observations.jsonl"]
    before = hashes(paths)
    with threadpool_limits(limits=1):
        cache = {}
        rows = []
        for case_id in cases:
            rows.append(measure_case(case_id, cache))
            if time.process_time()-started > CPU_BUDGET:
                raise RuntimeError("free-return consumer exceeded ten CPU seconds")
    source_after = hashes(bound_path(name) for name in source_before)
    after = hashes(bound_path(name) for name in before)
    if before != after or source_before != source_after:
        raise ValueError("free-return source or input changed during measurement")
    cpu = time.process_time()-started
    if cpu > CPU_BUDGET:
        raise RuntimeError("free-return consumer exceeded ten CPU seconds")
    return {
        "schema": SCHEMA, "primary_case": PRIMARY, "bindings": binding, "cases": rows,
        "interpretation": {"coupled_return": "closely reproduced by nearly free carrier transport",
                           "renewal_demonstrated": False, "free_transport_is_new_gravity_solution": False,
                           "frozen_return_alone_proves_free_limit": False,
                           "ambient_echo_proven": False, "new_trajectory": False,
                           "Duhamel_integral_certified": False},
        "source_hashes": source_after, "input_hashes": after, "input_sources_unchanged": True,
        "cpu_seconds": cpu, "cpu_budget_seconds": CPU_BUDGET, "numerical_threads": 1,
    }


def write_record(report, path=None):
    destination = OUTPUT if path is None else Path(path)
    if destination.resolve() != OUTPUT.resolve():
        raise ValueError("creation is restricted to the free-return-v1 owner")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False)+"\n")


def check_record(path=None):
    path = OUTPUT if path is None else Path(path)
    report = json.loads(path.read_text())
    if report.get("schema") != SCHEMA:
        raise ValueError("unexpected free-return schema")
    for group in ("source_hashes", "input_hashes"):
        if not report.get(group):
            raise ValueError("missing source/input bindings")
        for name, digest in report[group].items():
            if tidal.sha256_file(bound_path(name)) != digest:
                raise ValueError("bound hash changed: "+name)
    current = measure(include_controls=len(report["cases"]) > 1)
    if ([row["case_id"] for row in report["cases"]] != [row["case_id"] for row in current["cases"]]
            or report.get("interpretation") != current["interpretation"]):
        raise ValueError("free-return cases or interpretation changed")
    for saved, replayed in zip(report["cases"], current["cases"]):
        if saved["case_id"] != replayed["case_id"] or saved["join"] != replayed["join"]:
            raise ValueError("state continuation no longer matches")
        for key in ("weighted_field_relative_error", "source_covariance_relative_error"):
            if not math.isclose(saved[key], replayed[key], rel_tol=1e-8, abs_tol=1e-13):
                raise ValueError("free transport comparison does not reproduce")
        if not math.isclose(saved["fixed_child_window"]["free_probability_T8"],
                            replayed["fixed_child_window"]["free_probability_T8"], rel_tol=1e-10, abs_tol=1e-13):
            raise ValueError("free child probability does not reproduce")
    return {"ok": True, "recomputed_exact_free_operator": True, "evolved": False, "bytes_written": 0,
            "cpu_seconds": current["cpu_seconds"]}
