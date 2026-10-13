#!/usr/bin/env python3
"""Whole-cone field core v1: hashed v4 payload, streamed cells, open source error.

This records the reusable accumulator and one restored v4 cell. It does not
run family 14_1 or write a field certificate.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import time
import tracemalloc

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from fractions import Fraction
from types import SimpleNamespace
from recursive_horizons.nsc_ks_evaluation_binding import write_bytes_atomic
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_current_field_cone import (
    SCHEMA, V4_PROFILE_IDENTITY, WholeConeAccumulator, WholeConeFieldConfig,
    accumulate_segment, finalize_family, packed_radius_tail,
    restore_v4_profile_payload, scientific_digest, scientific_record,
)
from recursive_horizons.nsc_ks_difference_residual_polynomial import monomial_keys
from recursive_horizons.nsc_ks_local_profile_fourier_bound import (
    enclose_profile_fourier_local, profile_payload_digest,
    serialize_local_profile_fourier,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


OUTPUT = ROOT / "results/development/nsc-ks-current-field-cone-v1.json"
CHECKPOINT = ROOT / "results/development/nsc-ks-current-field-cone-v1.checkpoint.json"
V4_PATH = ROOT / "results/development/nsc-ks-current-field-pilot-v4.json"
SOURCE_PATHS = (
    "scripts/validate_nsc_ks_current_field_cone_v1.py",
    "src/recursive_horizons/nsc_ks_current_field_cone.py",
    "src/recursive_horizons/nsc_ks_local_profile_fourier_bound.py",
    "tests/test_nsc_ks_current_field_cone.py",
    "tests/test_nsc_ks_local_profile_fourier_bound.py",
    "docs/nsc-ks-current-field-cone-v1.md",
)
INPUT_PATHS = (
    "results/development/nsc-ks-current-field-pilot-v4.json",
    "results/development/nsc-ks-current-trajectory-pilot-v2.json",
    "src/recursive_horizons/nsc_ks_difference_error.py",
    "src/recursive_horizons/nsc_ks_difference_residual_polynomial.py",
)


def digest(path):
    from recursive_horizons.nsc_frozen_sources import selected_source_sha256
    path = Path(path)
    relative = path.relative_to(ROOT).as_posix() if path.is_absolute() else path.as_posix()
    return selected_source_sha256(ROOT, relative)



# Capture constants belong to the measured 1024-node, 12-column cell 122.
_capture_TIME_DEGREE, _capture_RECIPROCAL_ORDER = 8, 4
_capture_CAPTURE = ROOT / 'results/development/nsc-ks-current-trajectory-pilot-v2.json'

_capture_EXPECTED_IDENTITY = '0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0'

_capture_EFFECTIVE_PERIOD = Fraction(1759218604441, 4398046511104)

def _capture__scalar(value):
    return float(np.asarray(value).reshape(()))

def _capture_load_capture():
    if not _capture_CAPTURE.exists():
        raise FileNotFoundError('parent current-trajectory capture is required')
    record = json.loads(_capture_CAPTURE.read_text())
    if record.get('schema') != 'NSC-KS-CURRENT-TRAJECTORY-PILOT-v2':
        raise ValueError('unexpected current-trajectory capture schema')
    for name, expected in {**record['source_hashes'], **record['input_hashes']}.items():
        if digest(name) != expected:
            raise ValueError('pilot capture dependency changed: ' + name)
    payload = ROOT / record['payload']['path']
    if digest(payload) != record['payload']['sha256']:
        raise ValueError('pilot capture payload changed')
    if record['cell_index'] != 122 or record['family'] != [14, 1] or record['grid_nodes'] != 1024 or (record['source_columns'] != 12) or (len(record['source_rows']) != 4) or (record['profile_identity'] != _capture_EXPECTED_IDENTITY):
        raise ValueError('capture is not the declared current Cheb32 family14_1 cell122 pilot')
    with np.load(payload, allow_pickle=False) as arrays:
        data = {key: arrays[key] for key in arrays.files}
    return (record, data)

def _capture_load_family(capture):
    history = json.loads((ROOT / capture['history_path']).read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    identity = profile_identity(family, include_normal_window=True)
    if identity != capture['profile_identity'] or identity != _capture_EXPECTED_IDENTITY:
        raise ValueError('current history identity differs from the captured cell')
    return family

def _capture_geometry_bounds(family, capture):
    from recursive_horizons.nsc_ks_current_history_bounds import radius_bounds
    from recursive_horizons.nsc_ks_radius_enclosure import axial_profile_bounds, reciprocal_radius_tail
    bounds = radius_bounds(family, capture['rho_up'])
    rejected = (axial_profile_bounds(family.functions[0])[0], axial_profile_bounds(family.functions[1])[0])
    tail = reciprocal_radius_tail(bounds['w']['profile_bounds'], bounds['U']['profile_bounds'], normal_support=bounds['normal_support'], reference_radius_lower=bounds['reference_radius_lower'], axial_lower=bounds['axial_lower'], absolute_angular=abs(capture['angular']), order=_capture_RECIPROCAL_ORDER)
    return (bounds, rejected, tail)

def _capture_segment_from_capture(capture, arrays):
    grid = np.asarray(arrays['computational_z'], float)
    if len(grid) != capture['grid_nodes']:
        raise ValueError('captured spatial grid changed')
    spacing = Fraction(grid[1]) - Fraction(grid[0])
    period_q = spacing * len(grid)
    period = float(period_q)
    origin = float(grid[0])
    if abs(period - float(_capture_EFFECTIVE_PERIOD)) > 1e-18 and period_q != _capture_EFFECTIVE_PERIOD:
        if abs(period - float(_capture_EFFECTIVE_PERIOD)) > 1e-12:
            raise ValueError('captured spatial period is not the measured source-Fourier period')
    segment = TrajectorySegment(_capture__scalar(arrays['rho_start']), _capture__scalar(arrays['rho_end']), np.asarray(arrays['start']), np.asarray(arrays['end']), np.asarray(arrays['dense_corrections']))
    expected = capture['cell_rho_interval']
    if [segment.rho_end, segment.rho_start] != expected:
        raise ValueError('saved cell covers another rho interval')
    weights = np.asarray(arrays['source_weights'], float)
    energies = np.asarray(arrays['source_energies'], float)
    if len(weights) != capture['source_columns'] or len(energies) != len(weights):
        raise ValueError('captured source columns changed')
    return (segment, grid, origin, period, period_q, weights, energies)

capture_owner = SimpleNamespace(
    CAPTURE=_capture_CAPTURE, OUTPUT=ROOT / "results/development/nsc-ks-current-field-pilot-v3.json",
    TIME_DEGREE=_capture_TIME_DEGREE, RECIPROCAL_ORDER=_capture_RECIPROCAL_ORDER,
    EXPECTED_IDENTITY=_capture_EXPECTED_IDENTITY, EFFECTIVE_PERIOD=_capture_EFFECTIVE_PERIOD,
    load_capture=_capture_load_capture, load_family=_capture_load_family,
    segment_from_capture=_capture_segment_from_capture, geometry_bounds=_capture_geometry_bounds)

def synthetic_segment(rho_start, rho_end, *, scale=1.0, nsrc=1, nz=8):
    size = 2 * nsrc * (1 + nz)
    start = np.zeros(size, complex)
    start[0] = 0.5 * scale
    start[1] = 0.25 * scale
    harmonic = np.array([1, 1j, -1, -1j] * (nz // 4), complex) * 0.01 * scale
    start[2:2 + nz] = harmonic
    start[2 + nz:2 + 2 * nz] = 2 * harmonic
    return TrajectorySegment(
        rho_start, rho_end, start, 2 * start, np.zeros((6, size), complex))


def small_slab():
    coefficients = np.zeros((2, 8))
    coefficients[0, :4] = (0.012, -0.003, 0.0015, -0.0004)
    coefficients[1, :4] = (0.021, 0.004, -0.001, 0.0003)
    family = LocalIncomingFamily(coefficients)
    origin = family.center - 0.2
    bits, order = 80, 1
    rows = enclose_profile_fourier_local(
        family, monomial_keys(order), origin, 0.4, retained_index=4,
        taylor_order=5, exterior_panels=16, transition_panels=16,
        interior_panels=16, max_derivative=2, flat_denominator=32,
        relative_tolerance_bits=20, absolute_tolerance_bits=30,
        max_depth=8, bits=bits)
    payload = serialize_local_profile_fourier(rows)
    digest_value = profile_payload_digest(payload)
    config = WholeConeFieldConfig(
        family=(14, 1),
        profile_identity=profile_identity(family, include_normal_window=True),
        bits=bits, time_degree=4, reciprocal_order=order, spatial_count=8,
        period_origin=origin, period_length_numerator=2,
        period_length_denominator=5, mass=0.4, angular=0.8,
        source_weights=(1.0,), source_energies=(0.7,),
        radius_tail=packed_radius_tail(family, 1.03, 0.8, order, bits=bits),
        profile_payload_sha256=digest_value,
        settings={
            "retained_index": 4, "taylor_order": 5, "exterior_panels": 16,
            "transition_panels": 16, "interior_panels": 16, "bits": bits,
        })
    step = 1.0 / 8192
    rho = [1.03]
    for _ in range(3):
        rho.append(rho[-1] - step)
    segments = (
        synthetic_segment(rho[0], rho[1]),
        synthetic_segment(rho[1], rho[2], scale=1.1),
        synthetic_segment(rho[2], rho[3], scale=0.9),
    )
    return family, payload, config, segments


def v4_cell():
    raw = V4_PATH.read_bytes()
    recorded = json.loads(raw)
    _profiles, payload_digest = restore_v4_profile_payload(
        recorded, expected_file_sha256=digest(V4_PATH), file_bytes=raw)
    capture, arrays = capture_owner.load_capture()
    family = capture_owner.load_family(capture)
    segment, _grid, origin, _period, period_q, weights, energies = capture_owner.segment_from_capture(
        capture, arrays)
    if (capture["cell_index"] != 122 or capture["family"] != [14, 1]
            or capture["profile_identity"] != V4_PROFILE_IDENTITY):
        raise ValueError("capture is not the declared family-14_1 cell 122")
    config = WholeConeFieldConfig(
        family=(14, 1), profile_identity=V4_PROFILE_IDENTITY, bits=90,
        time_degree=capture_owner.TIME_DEGREE, reciprocal_order=capture_owner.RECIPROCAL_ORDER,
        spatial_count=capture["grid_nodes"], period_origin=origin,
        period_length_numerator=period_q.numerator,
        period_length_denominator=period_q.denominator,
        mass=capture["mass"], angular=capture["angular"],
        source_weights=tuple(map(float, weights)),
        source_energies=tuple(map(float, energies)),
        radius_tail=packed_radius_tail(
            family, capture["rho_up"], capture["angular"], capture_owner.RECIPROCAL_ORDER,
            bits=90),
        profile_payload_sha256=payload_digest,
        settings=dict(recorded["settings"]))
    acc = WholeConeAccumulator(
        config, recorded["profile_coefficient_payload"], family)
    accumulate_segment(acc, segment, cell_index=122)
    result = finalize_family(
        acc,
        propagation={
            "reference_initial": 0, "reference_residual": 0,
            "reference_offdiagonal_integral": 0, "difference_initial": (0, 0),
            "offdiagonal_integral": 0, "Bz_integral": 0, "M_integral": 0,
            "Mz_integral": 0, "maximum_absolute_energy": float(np.max(np.abs(energies))),
        })
    recorded_total = [restored_upper(value) for value in
                      recorded["bounds"]["total_continuous_normalized_residual"]]
    computed = [restored_upper(value) for value in result.bounds["last_cell_total"]]
    if any(not (bound >= previous) for bound, previous in zip(computed, recorded_total)):
        raise ValueError("restored v4 cell is outside its directed enclosure")
    return {
        "cell_index": 122,
        "family": [14, 1],
        "profile_payload_sha256": payload_digest,
        "recorded_total_continuous_normalized_residual":
            recorded["bounds"]["total_continuous_normalized_residual"],
        "restored_last_cell_total": result.bounds["last_cell_total"],
        "restored_contains_recorded": True,
        "propagate_difference_error_applied": True,
        "difference_matter_error_applied": False,
        "source_accuracy_included": False,
        "physical_rho1_source_error": None,
        "coverage": result.coverage,
        "status": result.status,
        "config_digest": result.config_digest,
        "prefix_digest": result.prefix_digest,
    }


def compute():
    before = {path: digest(path) for path in (*SOURCE_PATHS, *INPUT_PATHS)}
    started_cpu, started_wall = time.process_time(), time.monotonic()
    tracemalloc.start()
    family, profile_payload, config, segments = small_slab()
    acc = WholeConeAccumulator(config, profile_payload, family)
    accumulate_segment(acc, segments[0], cell_index=0)
    checkpoint = acc.checkpoint()
    accumulate_segment(acc, segments[1], cell_index=1)
    accumulate_segment(acc, segments[2], cell_index=2)
    dense = WholeConeAccumulator(config, profile_payload, family)
    for index, segment in enumerate(segments):
        accumulate_segment(dense, segment, cell_index=index)
    streamed = finalize_family(acc)
    dense_result = finalize_family(dense)
    if streamed.to_mapping()["bounds"] != dense_result.to_mapping()["bounds"]:
        raise ValueError("dense and streamed small-slab totals differ")
    resumed = WholeConeAccumulator.from_checkpoint(
        checkpoint, config, profile_payload, family)
    resumed.replay_segment(segments[0], cell_index=0)
    accumulate_segment(resumed, segments[1], cell_index=1)
    accumulate_segment(resumed, segments[2], cell_index=2)
    resumed_result = finalize_family(resumed)
    if resumed_result.to_mapping()["bounds"] != dense_result.to_mapping()["bounds"]:
        raise ValueError("checkpoint resume totals differ from the one-shot slab")
    cell = v4_cell()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    after = {path: digest(path) for path in (*SOURCE_PATHS, *INPUT_PATHS)}
    if before != after:
        raise ValueError("whole-cone v1 source changed during evaluation")
    record = {
        "schema": SCHEMA,
        "status": (
            "OPEN: whole-cone field core and one restored v4 cell; family "
            "14_1 scientific run was not completed; physical rho=1 source "
            "error unresolved"),
        "method": (
            "hash-validated v4 local-Fourier payload, streamed "
            "history-minus-reference residual cells, directed totals, "
            "prefix-bound checkpoint resume"),
        "residual_integral_measure": (
            "S=X_xi-h*L*X; integral_abs_R_d_rho=integral_abs_S_d_xi; "
            "no_second_cell_width_factor"),
        "apis": [
            "WholeConeFieldConfig", "WholeConeAccumulator",
            "WholeConeCheckpoint", "WholeConeFieldResult",
            "accumulate_segment", "finalize_family",
            "restore_v4_profile_payload",
            "whole_cone_propagate_difference_error",
            "whole_cone_difference_matter_error",
        ],
        "settings": dict(config.settings),
        "v4_cell": cell,
        "small_slab": {
            "completed_cells": streamed.completed_cells,
            "bounds": streamed.bounds,
            "dense_streamed_agreement": True,
            "checkpoint_resume_agreement": True,
            "config_digest": streamed.config_digest,
            "prefix_digest": streamed.prefix_digest,
        },
        "coverage": {
            "whole_time_cell": True,
            "whole_spatial_period": True,
            "all_history_cells": False,
            "all_source_families": False,
            "family_14_1_scientific_run": False,
            "source_columns_v4_cell": 12,
        },
        "certificate_use": False,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "source_accuracy_included": False,
        "physical_rho1_source_error": None,
        "family_14_1_scientific_run": False,
        "checkpoint_path": str(CHECKPOINT.relative_to(ROOT)).replace("\\", "/"),
        "next_required_step": (
            "stream every accepted current-history cell for representative "
            "families after the 60 GiB free-space guard is met, then enclose "
            "physical rho=1 source error"),
        "source_hashes": {path: before[path] for path in SOURCE_PATHS},
        "input_hashes": {path: before[path] for path in INPUT_PATHS},
        "runtime": {
            "CPU_seconds": time.process_time() - started_cpu,
            "wall_seconds": time.monotonic() - started_wall,
            "python": platform.python_version(),
        },
        "memory": {
            "tracemalloc_current_bytes": current,
            "tracemalloc_peak_bytes": peak,
        },
    }
    checkpoint_mapping = checkpoint.to_mapping()
    checkpoint_bytes = (
        json.dumps(checkpoint_mapping, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode()
    record["checkpoint_sha256"] = sha256(checkpoint_bytes).hexdigest()
    record["scientific_digest"] = scientific_digest(record)
    return record, checkpoint_mapping, checkpoint_bytes


def display(value):
    print(json.dumps({
        "status": value["status"],
        "v4_cell": {
            "cell_index": value["v4_cell"]["cell_index"],
            "restored_contains_recorded":
                value["v4_cell"]["restored_contains_recorded"],
            "physical_rho1_source_error":
                value["v4_cell"]["physical_rho1_source_error"],
        },
        "coverage": value["coverage"],
        "family_14_1_scientific_run": value["family_14_1_scientific_run"],
        "runtime": value.get("runtime"),
        "memory": value.get("memory"),
    }, indent=2, sort_keys=True))


def check_recorded(recorded):
    if recorded.get("schema") != SCHEMA:
        raise ValueError("unexpected whole-cone v1 schema")
    for path, expected in {**recorded["source_hashes"],
                           **recorded["input_hashes"]}.items():
        if digest(path) != expected:
            raise ValueError("whole-cone v1 dependency changed: " + path)
    if recorded["certificate_use"] or recorded["family_14_1_scientific_run"]:
        raise ValueError("core v1 record cannot be family-14_1 certificate evidence")
    if recorded.get("residual_integral_measure") != (
            "S=X_xi-h*L*X; integral_abs_R_d_rho=integral_abs_S_d_xi; "
            "no_second_cell_width_factor"):
        raise ValueError("normalized residual integral convention changed")
    if (recorded["physical_EXISTENCE_certificate"]
            or recorded["physical_NONEXISTENCE_certificate"]
            or recorded["source_accuracy_included"]
            or recorded["physical_rho1_source_error"] is not None):
        raise ValueError("physical rho=1 source error is unresolved")
    if scientific_digest(recorded) != recorded["scientific_digest"]:
        raise ValueError("scientific digest changed")
    v4 = json.loads(V4_PATH.read_text())
    recorded_total = [restored_upper(value) for value in
                      v4["bounds"]["total_continuous_normalized_residual"]]
    computed = [restored_upper(value) for value in
                recorded["v4_cell"]["restored_last_cell_total"]]
    if any(not (bound >= previous) for bound, previous in zip(computed, recorded_total)):
        raise ValueError("stored v4 cell comparison left the directed enclosure")
    checkpoint = json.loads(CHECKPOINT.read_text())
    from recursive_horizons.nsc_ks_current_field_cone import WholeConeCheckpoint
    loaded = WholeConeCheckpoint.from_mapping(checkpoint)
    checkpoint_bytes = CHECKPOINT.read_bytes()
    if sha256(checkpoint_bytes).hexdigest() != recorded["checkpoint_sha256"]:
        raise ValueError("whole-cone checkpoint bytes changed")
    if loaded.physical_rho1_source_error is not None:
        raise ValueError("physical rho=1 source error is unresolved")
    display(recorded)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--run", action="store_true")
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    if args.check:
        if not OUTPUT.exists() or not CHECKPOINT.exists():
            raise FileNotFoundError("whole-cone v1 record is missing")
        check_recorded(json.loads(OUTPUT.read_text()))
        return
    result, checkpoint, checkpoint_bytes = compute()
    if args.record:
        if CHECKPOINT.exists():
            if CHECKPOINT.read_bytes() != checkpoint_bytes:
                raise ValueError("existing whole-cone checkpoint differs")
        else:
            write_bytes_atomic(CHECKPOINT, checkpoint_bytes, exclusive=True)
        if digest(CHECKPOINT) != result["checkpoint_sha256"]:
            raise ValueError("published checkpoint digest differs")
        record_bytes = (
            json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode()
        if OUTPUT.exists():
            if OUTPUT.read_bytes() != record_bytes:
                raise ValueError("existing whole-cone record differs")
        else:
            write_bytes_atomic(OUTPUT, record_bytes, exclusive=True)
    elif args.replay:
        recorded = json.loads(OUTPUT.read_text())
        if scientific_record(result) != scientific_record(recorded):
            raise ValueError("whole-cone v1 replay differs")
    display(result)


if __name__ == "__main__":
    main()
