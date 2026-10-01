"""Independent triangular reference for the named memory-off control.

This causal negative control removes retained->exterior response while
keeping exterior->retained drive. Its generator is non-Hermitian and is
not proposed as a physical replacement for the coupled field model.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from .nsc_conformal_local_response import ConformalSchedule, _max_abs
from .nsc_conformal_local_response_successor import (
    DEV, OUTPUT_JSON as PRODUCTION_JSON, OUTPUT_NPZ as PRODUCTION_NPZ,
    load_inputs, source_bindings,
)
from .nsc_coupled_local_response import evolve_full_columns, occupations_and_coherence, sha256_file

OUTPUT_JSON = DEV / "nsc-conformal-memory-control-v1.json"
OUTPUT_NPZ = OUTPUT_JSON.with_suffix(".npz")


def triangular_generator(hamiltonian, observer):
    """G=H-P_E H P_A; no dense projector or propagator is allocated."""
    h = np.asarray(hamiltonian)
    v = np.asarray(observer)
    hv = h @ v
    a = v.conj().T @ hv
    return h - hv @ v.conj().T + v @ a @ v.conj().T


def execute(budget_s=30.):
    if OUTPUT_JSON.exists() or OUTPUT_NPZ.exists():
        raise FileExistsError("memory-control record already exists; use a separate successor")
    if not 0 < budget_s <= 30:
        raise ValueError("this supplementary reference has a 30 CPU-second cap")
    started = time.process_time()
    inputs = load_inputs()
    production_before = {path.name: sha256_file(path) for path in (PRODUCTION_JSON, PRODUCTION_NPZ)}
    production = json.loads(PRODUCTION_JSON.read_text())
    with np.load(PRODUCTION_NPZ, allow_pickle=False) as saved:
        times = np.array(saved["times"])
        memory_off = np.array(saved["occupation_memory_off"])
    schedule = ConformalSchedule(inputs)
    observer = inputs["observer"]

    def generator(mark):
        if time.process_time() - started > budget_s:
            raise RuntimeError("triangular reference exhausted its CPU budget")
        return triangular_generator(schedule(mark), observer)

    references = {}
    for substeps in (2, 4):
        columns = evolve_full_columns(generator, times, inputs["columns"], substeps=substeps)
        local = np.einsum("ij,tjk->tik", observer.conj().T, columns)
        occupation, coherence = occupations_and_coherence(local, inputs["weights"])
        references[substeps] = (occupation, coherence)
    error = _max_abs(memory_off - references[4][0])
    substep_indicator = _max_abs(references[2][0] - references[4][0])
    effect = float(production["controls"]["memory"]["occupation_movement"])
    g0 = generator(float(times[0]))
    h0 = schedule(float(times[0]))
    lower_block = g0 @ observer - observer @ (observer.conj().T @ g0 @ observer)
    upper_rows = observer.conj().T @ (g0 - h0)
    arrays = {"times": times, "occupation_reference_substeps_2": references[2][0],
              "occupation_reference_substeps_4": references[4][0],
              "coherence_reference_substeps_4": references[4][1]}
    record = {
        "schema": "NSC-CONFORMAL-MEMORY-CONTROL-v1", "status": "MEASURED_INDEPENDENT_MEMORY_REFERENCE",
        "case": inputs["case"], "time_window": [float(times[0]), float(times[-1])],
        "generator": "H-(I-V Vdagger) H V Vdagger", "negative_control": True,
        "physical_Hermitian_model_claimed": False, "geometry_rerun": False,
        "initial_state": "same actual Phi(0.2), all six original columns and Gaussian weights",
        "observer": "same fixed original T0 pair, no reset or phase change",
        "active_initial_C_AE": production["initial_blocks"]["C_AE_frobenius"],
        "memory_effect": effect, "streamed_memory_off_error_against_independent_reference": error,
        "error_fraction_of_memory_effect": error / effect,
        "reference_substeps_2_to_4_indicator": substep_indicator,
        "reference_indicator_fraction_of_memory_effect": substep_indicator / effect,
        "error_within_one_percent": bool(error < .01 * effect),
        "reference_indicator_within_one_percent": bool(substep_indicator < .01 * effect),
        "lower_exterior_from_retained_block_max_abs": _max_abs(lower_block),
        "retained_rows_unchanged_max_abs": _max_abs(upper_rows),
        "nonHermitian_frobenius": float(np.linalg.norm(g0 - g0.conj().T)),
        "dense_propagator_history_bytes": 0,
        "source_bindings": inputs["source_bindings"], "production_bindings": production_before,
        "cpu_seconds": time.process_time() - started, "cpu_budget_seconds": budget_s,
        "source_code_sha256": sha256_file(__file__), "time_path_enclosure_claimed": False,
    }
    if record["cpu_seconds"] > budget_s:
        raise RuntimeError("triangular reference exceeded its CPU budget")
    np.savez_compressed(OUTPUT_NPZ, **arrays)
    record["payload_sha256"] = sha256_file(OUTPUT_NPZ)
    OUTPUT_JSON.write_text(json.dumps(record, indent=2) + "\n")
    for _ in range(3):
        record["payload_bytes"] = OUTPUT_JSON.stat().st_size + OUTPUT_NPZ.stat().st_size
        OUTPUT_JSON.write_text(json.dumps(record, indent=2) + "\n")
    if source_bindings() != inputs["source_bindings"]:
        raise RuntimeError("source changed during the independent memory reference")
    if {path.name: sha256_file(path) for path in (PRODUCTION_JSON, PRODUCTION_NPZ)} != production_before:
        raise RuntimeError("memory reference changed frozen production bytes")
    return record
