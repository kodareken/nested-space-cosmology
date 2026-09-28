#!/usr/bin/env python3
"""Bounded cell-150 successor using the existing source and cone method."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'scripts')]
import derive_nsc_ks_current_trajectory_pilot_v2 as capture_owner
import run_nsc_ks_whole_cone_field_v1 as campaign_owner
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_current_field_cone import WholeConeAccumulator, canonical_digest
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment

CAPTURE = 'results/development/nsc-ks-current-trajectory-ramp-cell150-v1.json'
PAYLOAD = 'results/development/artifacts/nsc-ks-current-trajectory-ramp-cell150-v1.npz'
OUTPUT = 'results/development/nsc-ks-field-ramp-cell150-v1.json'
SELF = 'scripts/validate_nsc_ks_field_ramp_pilot.py'
CELL = 150


def configure_capture():
    # Reuse the original producer as an explicit adapter. Its source, grid,
    # solver settings, lock and CPU cap remain unchanged. Historical outputs
    # retain their original paths and bytes. The successor binds this adapter.
    capture_owner.CELL = CELL
    capture_owner.OUTPUT = ROOT/CAPTURE
    capture_owner.PAYLOAD = ROOT/PAYLOAD
    capture_owner.OWNERS = (*capture_owner.OWNERS, SELF)


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def proof_digest(record):
    return canonical_digest({k: v for k, v in record.items()
                             if k not in ('runtime', 'scientific_digest')})


def validate_pilot_shape(record):
    if record.get('scientific_digest') != proof_digest(record):
        raise ValueError('ramp scientific digest changed')
    if (record['schema'] != 'NSC-KS-FIELD-RAMP-CELL-PILOT-v1'
            or record['all_history_cells'] or record['full_family_source_coverage']
            or record['physical_rho1_source_error'] is not None
            or record['physical_local_gate'] != 'OPEN'):
        raise ValueError('one-cell pilot cannot close the physical gate')
    if record['cell']['cell_index'] != CELL:
        raise ValueError('incorrect ramp cell')
    if record['cell']['integral'] != record['cell']['total']:
        raise ValueError('normalized residual integral must not acquire another cell width')


def enclose(cpu_cap):
    if (ROOT/OUTPUT).exists():
        raise FileExistsError('ramp proof exists; use --check')
    captured = capture_owner.check()
    if captured['cell_index'] != CELL:
        raise ValueError('incorrect ramp cell')
    bound = campaign_owner.bind_production(ROOT)
    config = bound['config']
    selected = bound['selected']
    if captured['profile_identity'] != bound['profile_identity']:
        raise ValueError('ramp history differs from campaign')
    with np.load(ROOT/PAYLOAD, allow_pickle=False) as saved:
        for key, expected in (
            ('source_weights', selected['source'].column_weights),
            ('source_energies', selected['source'].energies),
            ('source_covariance', selected['source'].covariance),
            ('initial_columns', selected['initial_columns']),
            ('computational_z', bound['grid']), ('target_z', bound['target']),
        ):
            if not np.array_equal(saved[key], expected):
                raise ValueError('ramp input binding mismatch: '+key)
        segment = TrajectorySegment(float(saved['rho_start'].item()), float(saved['rho_end'].item()),
            saved['start'].copy(), saved['end'].copy(), saved['dense_corrections'].copy())
    accumulator = WholeConeAccumulator(config, bound['payload'], bound['family'])
    started_cpu, started_wall = time.process_time(), time.perf_counter()
    def stop(*_):
        raise TimeoutError('ramp enclosure CPU cap; capture preserved')
    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, cpu_cap)
    try:
        cell = accumulator.accumulate_segment(segment, cell_index=CELL)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    record = {
        'schema': 'NSC-KS-FIELD-RAMP-CELL-PILOT-v1',
        'status': 'OPEN: one active-ramp cell enclosed; full trajectory/source budget incomplete',
        'cell': cell, 'config_digest': config.digest(),
        'profile_identity': bound['profile_identity'],
        'source_preparation_digest': captured['source_preparation_digest'],
        'family': [14, 1], 'selected_source_columns': len(config.source_weights),
        'all_history_cells': False, 'full_family_source_coverage': False,
        'physical_rho1_source_error': None, 'physical_local_gate': 'OPEN',
        'runtime': {'CPU_seconds': time.process_time()-started_cpu,
                    'wall_seconds': time.perf_counter()-started_wall, 'CPU_cap': cpu_cap},
        'source_hashes': {SELF: digest(SELF), **config.settings['implementation_hashes']},
        'input_hashes': {CAPTURE: digest(CAPTURE), PAYLOAD: digest(PAYLOAD),
                        campaign_owner.V4_RELATIVE: digest(campaign_owner.V4_RELATIVE),
                        **selected['archive_hashes']},
    }
    record['scientific_digest'] = proof_digest(record)
    publish_exclusive_file(ROOT, OUTPUT,
        (json.dumps(record, indent=2, sort_keys=True)+'\n').encode())
    return record


def check():
    capture_owner.check()
    record = json.loads((ROOT/OUTPUT).read_text())
    validate_pilot_shape(record)
    for path, expected in {**record['source_hashes'], **record['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('ramp dependency changed: '+path)
    if record['config_digest'] != campaign_owner.bind_production(ROOT)['config'].digest():
        raise ValueError('ramp numerical/source binding changed')
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--capture', action='store_true')
    modes.add_argument('--enclose', action='store_true')
    modes.add_argument('--check', action='store_true', help='authenticate saved result; no re-enclosure')
    parser.add_argument('--cpu-cap', type=float, default=900.)
    args = parser.parse_args()
    if not np.isfinite(args.cpu_cap) or not args.cpu_cap > 0:
        raise ValueError('positive finite CPU cap required')
    configure_capture()
    result = capture_owner.run() if args.capture else (enclose(args.cpu_cap) if args.enclose else check())
    print(json.dumps(result if args.capture else {
        'status': result['status'], 'cell': result['cell'], 'runtime': result['runtime'],
        'physical_local_gate': 'OPEN'}, indent=2))
