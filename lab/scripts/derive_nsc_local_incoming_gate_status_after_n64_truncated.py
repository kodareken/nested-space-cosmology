#!/usr/bin/env python3
"""Gate status after three evolved n=64 truncated steps. Still OPEN."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BEST = ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json'
STEPS = (
    ROOT/'results/development/nsc-ks-coupled-newton-n64-truncated-iterate7.json',
    ROOT/'results/development/nsc-ks-coupled-newton-n64-truncated-iterate8.json',
    ROOT/'results/development/nsc-ks-coupled-newton-n64-truncated-iterate9.json',
)
OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-after-n64-truncated.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_after_n64_truncated.py',
    'docs/nsc-local-incoming-gate-status-after-n64-truncated.md',
)
EXISTENCE_TOLERANCE = 3e-11
STALL = 'n64_truncated_N_offset_independent_of_step'


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    record = json.loads(path.read_text())
    if record['residual_improved_on_all_components']:
        raise ValueError('an n=64 step that improved both components is not this stall')
    if record['physical_EXISTENCE_certificate'] or record['physical_NONEXISTENCE_certificate']:
        raise ValueError('these n=64 steps are not a physical certificate')
    return record


def compute():
    best = json.loads(BEST.read_text())
    if not best['residual_improved_on_all_components']:
        raise ValueError('iterate6 remains the campaign best only if it improved both components')
    parent = np.asarray(best['action_gradient'], float)
    parent_max = np.asarray(best['constraint_maxima'], float)
    rows = []
    for path in STEPS:
        record = load(path)
        measured = np.asarray(record['constraint_maxima'], float)
        predicted = np.asarray(record['linear_predicted_all_node_residual_maxima'], float)
        gradient = np.asarray(record['action_gradient'], float)
        if gradient.shape != parent.shape:
            raise ValueError('n=64 residual must stay on the iterate6 nodes')
        shift = gradient[:, 0] - parent[:, 0]
        if np.any(measured <= EXISTENCE_TOLERANCE):
            raise ValueError('a component at 3e-11 would require an enclosed-epsilon certificate')
        rows.append({
            'path': str(path.relative_to(ROOT)),
            'profile_identity': record['profile_identity'],
            'kept_modes': record['kept_modes'],
            'step_scale': record['step_scale'],
            'step_norm': record['step_norm'],
            'predicted_maxima': predicted.tolist(),
            'measured_maxima': measured.tolist(),
            'prediction_minus_measured': record['prediction_minus_measured'],
            'N_shift_mean': float(np.mean(shift)),
            'N_shift_std': float(np.std(shift)),
            'residual_improved_on_all_components': False,
        })
    rank64 = [row for row in rows if row['kept_modes'] == 64]
    if len(rank64) != 2:
        raise ValueError('the stall uses the two evolved rank-64 steps')
    if abs(rank64[0]['N_shift_mean'] - rank64[1]['N_shift_mean']) > 1e-6:
        raise ValueError('rank-64 N offsets are not independent of the damping')
    if min(abs(row['prediction_minus_measured'][1]) for row in rank64) > 1e-9:
        raise ValueError('rank-64 beta must still track the linear prediction')
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-AFTER-N64-TRUNCATED-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: n=64 truncated steps leave an N offset and do not beat iterate6',
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; original fixed source',
        'history_class': 'local source-fixed pure-radius compatible histories',
        'interval': 'S(1)+[.12,.18]',
        'best_profile_identity': best['profile_identity'],
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'best_constraint_maxima': parent_max.tolist(),
        'evolved_steps': rows,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'enclosed_error': None,
        'named_search_stall': STALL,
        'next_named_owner': 'same declared (w,U) class on I; remove the N offset before another damped step',
        's5_retracted': True,
        'scope': {
            'open_is_not_done': True,
            'manuscript_and_public_release_blocked': True,
            'arxiv_blocked': True,
            'no_class_rewrite': True,
            'metric_timestep': False,
            'failed_optimization_is_not_nonexistence': True,
        },
        'source_hashes': {path: digest(ROOT/path) for path in OWNERS},
        'input_hashes': {str(path.relative_to(ROOT)): digest(path) for path in (BEST, *STEPS)},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('n=64 gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('n=64 gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'named_search_stall': result['named_search_stall'],
        'best_constraint_maxima': result['best_constraint_maxima'],
        'measured_maxima': [row['measured_maxima'] for row in result['evolved_steps']],
        'physical_local_gate': result['physical_local_gate'],
    }, indent=2))
