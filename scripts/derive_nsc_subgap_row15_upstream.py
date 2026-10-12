#!/usr/bin/env python3
"""Replayable witness for original group14/low16_1 row 15 and its signed partner.

The selection and numerical settings are the mixed-transport pilot's. One
positive Jost solve and one Bloch capture build the witness. ``--check``
replays the saved phase prefix and Bloch trace without either solver. Source
covariance errors stay unweighted. No N/beta aggregate is formed. The local
gate stays OPEN.
"""
import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))

from pilot_nsc_massive_jost_mixed_transport import ENERGY, LOOSE_COVARIANCE, ROW
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import (
    file_digest, hex_float, implementation_hashes)
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_massive_jost_transport_bound import (
    original_dense_solution, phase_transport_segments)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_subgap_row_witness import (
    BLOCH_MAX_STEP, BLOCH_Y_START, DEGREE, DEFECT_SUBDIVISIONS, FRAME_ORDER,
    INNER_REQUEST, METRIC_TERMS, MISSING_STEP, PAYLOAD_NAMES, PHASE_COLUMNS,
    BITS, BLOCH_DEGREE, SOLVER, TUBE,
    bloch_record, finalize_bloch, load_witness_arrays, phase_prefix_from_dense,
    phase_record, prepare_phase, requested_inner_radius, require_exact_source,
    require_same_phase_segments, segments_of_prefix)
from recursive_horizons.nsc_subgap_source_covariance import capture_bloch


OUTPUT = 'results/development/nsc-subgap-row15-upstream-v1.json'
PAYLOAD = 'results/development/artifacts/nsc-subgap-row15-upstream-v1.npz'
PANEL = 'group14/low16_1'
PILOT = 'scripts/pilot_nsc_massive_jost_mixed_transport.py'
PHASE_SEPARATION = 1e-4
OWNERS = (
    'scripts/derive_nsc_subgap_row15_upstream.py',
    'src/recursive_horizons/nsc_subgap_row_witness.py',
    'tests/test_nsc_subgap_row_witness.py',
    'docs/nsc-subgap-row-witness.md',
)
LABEL_KEYS = (
    'schema', 'status', 'source_panel', 'source_row', 'group', 'batch_rows',
    'source_column_slice', 'positive_energy_hex', 'negative_energy_hex',
    'mass_hex', 'positive_angular_hex', 'negative_angular_hex', 'rho_up_hex',
    'source_kappa_hex', 'horizon_rho_hex', 'outer_radius_hex',
    'inner_request_offset_hex', 'source_signs', 'angular_signs',
    'positive_preparation_digest', 'negative_preparation_digest',
    'positive_source_digest', 'negative_source_digest',
    'original_quadrature_weight_hex', 'solver', 'transport', 'bloch_settings',
    'separation_scales', 'phase_columns', 'phase_prefix_rows', 'bloch_trace_rows',
    'phase_nfev', 'bloch_nfev', 'covariance_error_unweighted',
    'quadrature_weight_applied', 'negative_error_copied_from_positive',
    'negative_complement_applied', 'phase_component_only',
    'amplitude_component_retained', 'dop853_dense_output_reconstructed',
    'replay_uses_saved_phase_and_bloch_witness', 'preparation_jost_rerun_on_replay',
    'bloch_ode_rerun_on_replay', 'source_columns_replaced',
    'horizon_sewing_error', 'frame_error_included_in_initial_bloch',
    'separate_frame_error_addition_required', 'physical_rho1_source_error', 'n_beta_aggregate',
    'physical_upstream_budget_component', 'all_source_families',
    'method_too_loose', 'physical_local_gate', 'missing_step',
)


def canonical_json(record):
    return (json.dumps(record, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()


def load_binding():
    """Restore the pilot's one positive row and the archived signed partner."""
    archive = RetainedUpstreamArchive(ROOT)
    config = archive.meta['config']
    preparation = PairedHorizonSeedMap(
        config['horizon_rho'], config['surface_gravity'], config['omega'],
        config['horizon_offset'], config['scattering_tolerance'], config['outer_floor'])
    batches = [
        batch for batch, _channel in archive.family_entries((14, 1))
        if batch.original_panel == PANEL and tuple(batch.rows) == (0, 16)]
    if len(batches) != 2 or {batch.energy_sign for batch in batches} != {-1, 1}:
        raise ValueError('original signed upstream pair missing')
    positive = next(batch for batch in batches if batch.energy_sign > 0)
    negative = next(batch for batch in batches if batch.energy_sign < 0)
    columns = slice(3 * ROW, 3 * ROW + 3)
    energy = float(positive.source.energies[3 * ROW])
    if (energy != ENERGY or positive.mass != np.pi / 2
            or positive.angular != np.sqrt(5.0)):
        raise ValueError('original group14/low16_1 row 15 channel changed')
    if (positive.angular_sign != 1 or negative.angular_sign != -1
            or negative.mass != positive.mass or negative.angular != -positive.angular
            or positive.rho_up != negative.rho_up or not positive.rho_up >= 1.03
            or positive.group != 14 or negative.group != 14):
        raise ValueError('source label mismatch')
    if not np.all(positive.source.energies[columns] == energy):
        raise ValueError('source label mismatch')
    if not np.all(negative.source.energies[columns] == -energy):
        raise ValueError('source label mismatch')
    panel_energies = archive._arrays[PANEL + '/energies']
    panel_weights = archive._arrays[PANEL + '/weights']
    if float(panel_energies[ROW]) != energy or not float(panel_weights[ROW]) > 0:
        raise ValueError('source label mismatch')
    return {
        'archive': archive,
        'config': config,
        'preparation': preparation,
        'positive': positive,
        'negative': negative,
        'energy': energy,
        'columns': columns,
        'weight': float(panel_weights[ROW]),
        'positive_columns': np.array(positive.initial_columns[:, columns], dtype=np.complex128, copy=True),
        'negative_columns': np.array(negative.initial_columns[:, columns], dtype=np.complex128, copy=True),
        'positive_covariance': np.array(positive.source.covariance[columns, columns], dtype=np.complex128, copy=True),
        'negative_covariance': np.array(negative.source.covariance[columns, columns], dtype=np.complex128, copy=True),
    }


def witness_arrays(prefix, trace, binding, outer_radius, phase_nfev, bloch_nfev):
    return {
        'phase_prefix': np.ascontiguousarray(prefix, dtype=np.float64),
        'bloch_trace': np.ascontiguousarray(trace, dtype=np.float64),
        'positive_columns': np.ascontiguousarray(binding['positive_columns']),
        'negative_columns': np.ascontiguousarray(binding['negative_columns']),
        'positive_covariance': np.ascontiguousarray(binding['positive_covariance']),
        'negative_covariance': np.ascontiguousarray(binding['negative_covariance']),
        'outer_radius': np.array([float(outer_radius)], dtype=np.float64),
        'phase_nfev': np.array([int(phase_nfev)], dtype=np.int64),
        'bloch_nfev': np.array([int(bloch_nfev)], dtype=np.int64),
    }


def load_payload(raw, expected_sha256):
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('row witness payload hash changed')
    with np.load(BytesIO(raw), allow_pickle=False) as saved:
        arrays = {name: saved[name].copy() for name in saved.files}
    loaded = load_witness_arrays(arrays)
    if set(arrays) != set(PAYLOAD_NAMES):
        raise ValueError('complete row witness payload required')
    return loaded


def binding_fields(binding, outer_radius, phase_nfev, bloch_nfev, prefix_rows, bloch_rows):
    positive, negative = binding['positive'], binding['negative']
    horizon = float(binding['preparation'].horizon_rho)
    return {
        'schema': 'NSC-SUBGAP-ROW15-UPSTREAM-v1',
        'status': (
            'BOUNDED: original group14/low16_1 row 15 and its signed upstream '
            'partner; source budget and local gate OPEN'),
        'source_panel': PANEL,
        'source_row': ROW,
        'group': 14,
        'batch_rows': [0, 16],
        'source_column_slice': [3 * ROW, 3 * ROW + 3],
        'positive_energy_hex': hex_float(binding['energy']),
        'negative_energy_hex': hex_float(-binding['energy']),
        'mass_hex': hex_float(positive.mass),
        'positive_angular_hex': hex_float(positive.angular),
        'negative_angular_hex': hex_float(negative.angular),
        'rho_up_hex': hex_float(positive.rho_up),
        'source_kappa_hex': hex_float(binding['config']['surface_gravity']),
        'horizon_rho_hex': hex_float(horizon),
        'outer_radius_hex': hex_float(outer_radius),
        'inner_request_offset_hex': hex_float(INNER_REQUEST),
        'source_signs': [1, -1],
        'angular_signs': [1, -1],
        'positive_preparation_digest': positive.preparation_digest,
        'negative_preparation_digest': negative.preparation_digest,
        'positive_source_digest': positive.source.digest,
        'negative_source_digest': negative.source.digest,
        'original_quadrature_weight_hex': hex_float(binding['weight']),
        'solver': {
            'method': 'DOP853',
            'order': SOLVER['order'],
            'radial_collar_hex': hex_float(SOLVER['radial_collar']),
            'rtol_hex': hex_float(SOLVER['rtol']),
            'atol_hex': hex_float(SOLVER['atol']),
        },
        'transport': {
            'degree': DEGREE,
            'metric_terms': METRIC_TERMS,
            'defect_subdivisions': DEFECT_SUBDIVISIONS,
            'tube': TUBE,
            'bits': BITS,
            'frame_order': FRAME_ORDER,
        },
        'bloch_settings': {
            'y_start_hex': hex_float(BLOCH_Y_START),
            'max_step_hex': hex_float(BLOCH_MAX_STEP),
            'degree': BLOCH_DEGREE,
            'bits': BITS,
        },
        'separation_scales': {
            'phase_upper_hex': hex_float(PHASE_SEPARATION),
            'covariance_upper_hex': hex_float(LOOSE_COVARIANCE),
        },
        'phase_columns': list(PHASE_COLUMNS),
        'phase_prefix_rows': int(prefix_rows),
        'bloch_trace_rows': int(bloch_rows),
        'phase_nfev': int(phase_nfev),
        'bloch_nfev': int(bloch_nfev),
        'covariance_error_unweighted': True,
        'quadrature_weight_applied': False,
        'negative_error_copied_from_positive': False,
        'negative_complement_applied': True,
        'phase_component_only': True,
        'amplitude_component_retained': False,
        'dop853_dense_output_reconstructed': True,
        'replay_uses_saved_phase_and_bloch_witness': True,
        'preparation_jost_rerun_on_replay': False,
        'bloch_ode_rerun_on_replay': False,
        'source_columns_replaced': False,
        'horizon_sewing_error': None,  # no separate budget entry; frame balls are propagated
        'frame_error_included_in_initial_bloch': True,
        'separate_frame_error_addition_required': False,
        'physical_rho1_source_error': None,
        'n_beta_aggregate': None,
        'physical_upstream_budget_component': None,
        'all_source_families': False,
        'method_too_loose': False,
        'physical_local_gate': 'OPEN',
        'missing_step': MISSING_STEP,
    }


def require_saved_labels(saved, expected):
    if not isinstance(saved, dict):
        raise ValueError('source label mismatch')
    for key in LABEL_KEYS:
        if key not in saved or saved[key] != expected[key]:
            raise ValueError('source label mismatch')
    payload = saved.get('payload')
    if not isinstance(payload, dict) or payload.get('path') != PAYLOAD:
        raise ValueError('complete row witness payload required')
    if not isinstance(payload.get('sha256'), str) or not isinstance(payload.get('bytes'), int):
        raise ValueError('complete row witness payload required')


def require_separated(record):
    phase = float(restored_upper(record['phase_error_inner_upper']))
    positive = float(restored_upper(record['positive_covariance_error_upper']))
    negative = float(restored_upper(record['negative_covariance_error_upper']))
    if not (phase < PHASE_SEPARATION and positive < LOOSE_COVARIANCE
            and negative < LOOSE_COVARIANCE):
        raise ArithmeticError(
            'row 15 phase or covariance bound does not separate the archived source')
    return phase, positive, negative


def require_equal_record(saved, fresh):
    if saved == fresh:
        return
    changed = sorted(key for key in set(saved) | set(fresh) if saved.get(key) != fresh.get(key))
    preview = ', '.join(changed[:8])
    raise ValueError('row witness replay differs from the saved record: ' + preview)


def calculate(saved=None, audit=None):
    """Capture one positive row, or replay its saved phase and Bloch witnesses."""
    binding = load_binding()
    horizon = float(binding['preparation'].horizon_rho)
    inner = requested_inner_radius(horizon)
    background = binding['preparation'].background
    if saved is None:
        mode = solve_jost(
            background, binding['energy'], binding['positive'].mass,
            binding['positive'].angular, **SOLVER)
        dense = original_dense_solution(mode.run)
        prefix = phase_prefix_from_dense(dense, horizon, inner)
        require_same_phase_segments(
            phase_transport_segments(dense, horizon, inner),
            segments_of_prefix(prefix, horizon, inner))
        outer_radius = float(mode.outer_radius)
        phase_nfev = int(mode.run.nfev)
        captured = {
            'prefix': np.array(prefix, copy=True),
            'original_segments': phase_transport_segments(dense, horizon, inner),
            'horizon_rho': horizon,
            'inner_radius': inner,
            'outer_radius': outer_radius,
            'positive_columns': binding['positive_columns'],
            'negative_columns': binding['negative_columns'],
            'positive_covariance': binding['positive_covariance'],
            'negative_covariance': binding['negative_covariance'],
            'positive_angular': binding['positive'].angular,
            'negative_angular': binding['negative'].angular,
        }
        trace = None
        bloch_nfev = None
    else:
        saved_record, saved_raw = saved
        if not isinstance(saved_raw, (bytes, bytearray)):
            raise ValueError('complete row witness payload required')
        saved_raw = bytes(saved_raw)
        payload_meta = saved_record.get('payload') if isinstance(saved_record, dict) else None
        if not isinstance(payload_meta, dict) or not isinstance(payload_meta.get('sha256'), str):
            raise ValueError('complete row witness payload required')
        if payload_meta.get('bytes') != len(saved_raw):
            raise ValueError('row witness payload hash changed')
        loaded = load_payload(saved_raw, payload_meta['sha256'])
        require_exact_source(
            loaded, binding['positive_columns'], binding['negative_columns'],
            binding['positive_covariance'], binding['negative_covariance'])
        prefix = loaded['phase_prefix']
        trace = loaded['bloch_trace']
        outer_radius = loaded['outer_radius']
        phase_nfev = loaded['phase_nfev']
        bloch_nfev = loaded['bloch_nfev']
        labels = binding_fields(
            binding, outer_radius, phase_nfev, bloch_nfev, len(prefix), len(trace))
        require_saved_labels(saved_record, labels)
        captured = None
    prepared = prepare_phase(
        prefix, energy=binding['energy'], mass=binding['positive'].mass,
        angular=binding['positive'].angular, background=background,
        outer_radius=outer_radius, horizon_rho=horizon,
        kappa=binding['config']['surface_gravity'], rho_up=binding['positive'].rho_up)
    if trace is None:
        # First capture only. Replay passes the saved trace and never reaches this call.
        trace, bloch_nfev = capture_bloch(
            prepared.model, prepared.initial, BLOCH_Y_START,
            float(prepared.target.mid()), max_step=BLOCH_MAX_STEP)
        bloch_nfev = int(bloch_nfev)
    validated, positive_error, negative_error = finalize_bloch(
        prepared, trace, binding['positive_columns'], binding['positive_covariance'],
        binding['negative_columns'], binding['negative_covariance'])
    arrays = witness_arrays(
        prefix, trace, binding, outer_radius, phase_nfev, bloch_nfev)
    raw = deterministic_npz_bytes(arrays)
    record = binding_fields(
        binding, outer_radius, phase_nfev, bloch_nfev, len(prefix), len(trace))
    record.update(phase_record(prepared))
    record.update(bloch_record(validated, positive_error, negative_error))
    if record['phase_cells'] != record['phase_prefix_rows']:
        raise ArithmeticError('transported cell count does not match the saved prefix')
    if record['bloch_cells'] != record['bloch_trace_rows']:
        raise ArithmeticError('validated Bloch cells do not match the saved trace')
    require_separated(record)
    record['payload'] = {'path': PAYLOAD, 'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)}
    record['source_hashes'] = implementation_hashes(ROOT, owners=OWNERS)
    inputs = dict(binding['archive'].input_hashes)
    inputs[PILOT] = file_digest(ROOT / PILOT)
    record['input_hashes'] = inputs
    if saved is not None:
        require_equal_record(saved_record, record)
        if raw != saved_raw:
            raise ValueError('row witness replay differs from the saved record: payload')
    if audit is not None and captured is not None:
        audit.update(captured)
        audit['positive_error'] = positive_error
        audit['negative_error'] = negative_error
        audit['validated'] = validated
    return record, raw


def publish_witness(root, record, raw):
    root = Path(root)
    if (root / OUTPUT).exists() or (root / PAYLOAD).exists():
        raise FileExistsError('row witness already exists; use --check')
    publish_exclusive_file(root, PAYLOAD, raw)
    publish_exclusive_file(root, OUTPUT, canonical_json(record))


def verify_saved(root=ROOT):
    root = Path(root)
    saved_bytes = (root / OUTPUT).read_bytes()
    raw = (root / PAYLOAD).read_bytes()
    saved = json.loads(saved_bytes)
    record, new_raw = calculate(saved=(saved, raw))
    if saved_bytes != canonical_json(record) or raw != new_raw:
        raise ValueError('row witness proof, payload or dependencies changed')
    return record


def _summary(record, started, wall):
    phase, positive, negative = require_separated(record)
    return {
        'status': record['status'],
        'phase_cells': record['phase_cells'],
        'bloch_cells': record['bloch_cells'],
        'phase_error_upper': phase,
        'positive_covariance_error_upper': positive,
        'negative_covariance_error_upper': negative,
        'physical_local_gate': record['physical_local_gate'],
        'n_beta_aggregate': None,
        'cpu_seconds': time.process_time() - started,
        'wall_seconds': time.perf_counter() - wall,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    started, wall = time.process_time(), time.perf_counter()
    if args.record:
        if (ROOT / OUTPUT).exists() or (ROOT / PAYLOAD).exists():
            raise FileExistsError('row witness already exists; use --check')
        record, raw = calculate()
        publish_witness(ROOT, record, raw)
    else:
        record = verify_saved(ROOT)
    print(json.dumps(_summary(record, started, wall), indent=2))
