#!/usr/bin/env python3
"""Record the linear floor under the measured rank-29 Jacobian.

Well-conditioned truncations no longer beat both components by more than the
last prediction error. The tail that does is ill-conditioned and is not evolved.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n64_truncated_iterate as I64

PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-rank29.json'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-linear-floor.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_linear_floor.py',
    'docs/nsc-ks-n64-primal-norm-linear-floor.md',
)
PARENT_IDENTITY = 'c2404e366476e2771f2ca3778e32d6fd53edd062b0c8d3e876f0dfcacc491ea3'


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def scan():
    parent = json.loads(PARENT.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY:
        raise ValueError('linear-floor parent must be the measured rank-29 history')
    if not parent['newton_step_accepted']:
        raise ValueError('linear-floor parent must be an accepted step')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    measured = np.asarray(parent['constraint_maxima'], float)
    floor = np.abs(np.asarray(parent['prediction_minus_measured'], float))
    _delta, _residual, _matrix, spectrum = I64.rank_truncated_delta(gradient, jacobian, 1)
    rows = []
    for rank in range(1, int(spectrum.size) + 1):
        delta, residual, matrix, singular = I64.rank_truncated_delta(gradient, jacobian, rank)
        raw = float(np.linalg.norm(delta))
        ratio = float(singular[0] / singular[rank - 1])
        step = delta if raw <= 1.0 else delta / raw
        predicted = I64.predicted_maxima(residual, matrix, step)
        gain = measured - predicted
        rows.append({
            'kept_modes': rank,
            'raw_step_norm': raw,
            'used_step_norm': float(np.linalg.norm(step)),
            'singular_ratio': ratio,
            'predicted_all_node_residual_maxima': predicted.tolist(),
            'predicted_gain': gain.tolist(),
            'beats_parent_on_both_components': bool(np.all(gain > 0)),
            'gain_exceeds_last_prediction_error': bool(np.all(gain > floor)),
        })
    small = [row for row in rows if row['raw_step_norm'] <= 1e-6 and row['singular_ratio'] < 1e4]
    clearing = [row for row in small if row['gain_exceeds_last_prediction_error']]
    damped = None
    rank74 = next(row for row in rows if row['kept_modes'] == 74)
    delta, residual, matrix, _singular = I64.rank_truncated_delta(gradient, jacobian, 74)
    step = delta * (0.01 / rank74['raw_step_norm'])
    predicted = I64.predicted_maxima(residual, matrix, step)
    damped = {
        'kept_modes': 74,
        'step_norm': 0.01,
        'singular_ratio': rank74['singular_ratio'],
        'predicted_all_node_residual_maxima': predicted.tolist(),
        'predicted_gain': (measured - predicted).tolist(),
        'evolved': False,
        'reason': 'N gain stays inside the nonlinearity already measured at step norm 0.01',
    }
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-LINEAR-FLOOR-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: well-conditioned linear gain is inside the last prediction error',
        'stall': 'n64_primal_norm_linear_gain_inside_prediction_error',
        'parent_profile_identity': parent['profile_identity'],
        'parent_constraint_maxima': measured.tolist(),
        'last_prediction_error_abs': floor.tolist(),
        'best_small_clearing_rank': None if not clearing else min(row['kept_modes'] for row in clearing),
        'rank27_predicted_gain': next(row['predicted_gain'] for row in rows if row['kept_modes'] == 27),
        'rank74_raw_step_norm': rank74['raw_step_norm'],
        'rank74_singular_ratio': rank74['singular_ratio'],
        'rank80_raw_step_norm': next(row['raw_step_norm'] for row in rows if row['kept_modes'] == 80),
        'rank80_singular_ratio': next(row['singular_ratio'] for row in rows if row['kept_modes'] == 80),
        'damped_rank74_norm_0_01': damped,
        'rows': rows,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'enclosed_error': None,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {str(PARENT.relative_to(ROOT)): digest(PARENT)},
    }


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = scan()
    for key in (
        'stall', 'parent_constraint_maxima', 'last_prediction_error_abs',
        'rank27_predicted_gain', 'rank74_raw_step_norm', 'rank80_raw_step_norm',
        'damped_rank74_norm_0_01',
    ):
        if saved[key] != fresh[key]:
            raise ValueError('linear-floor register changed: '+key)
    if saved['physical_EXISTENCE_certificate'] or saved['physical_NONEXISTENCE_certificate']:
        raise ValueError('linear floor is not a gate certificate')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        record = scan()
        OUTPUT.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'stall': record['stall'],
            'best_small_clearing_rank': record['best_small_clearing_rank'],
            'rank27_predicted_gain': record['rank27_predicted_gain'],
            'rank74_raw_step_norm': record['rank74_raw_step_norm'],
            'rank80_raw_step_norm': record['rank80_raw_step_norm'],
            'damped_rank74_predicted_gain': record['damped_rank74_norm_0_01']['predicted_gain'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({'stall': saved['stall'], 'status': saved['status']}, indent=2))
