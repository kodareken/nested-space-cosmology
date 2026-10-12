#!/usr/bin/env python3
"""Record the sup-norm floor under the measured factor-200 history.

The well-conditioned block no longer clears the new prediction error by a
margin that survived the last measurement. Missing 1e-8 error terms cannot
reach the 3e-11 existence tolerance.
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
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n64_truncated_iterate as I64

PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-sup16.json'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-sup-floor.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_sup_floor.py',
    'docs/nsc-ks-n64-primal-norm-sup-floor.md',
)
PARENT_IDENTITY = '3ba9f4c2fddecddd55d2017431bff71768501b99084815e0d229967fd10e4ea3'
EXISTENCE = 3e-11
SCALE = 1e11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def largest_factor(matrix, residual, right, col, measured, floor, kept):
    nodes = measured.size and residual.size // 2
    basis = right[:kept].T / col[:, None]
    response = matrix @ basis
    low, high = 0.0, 30.0
    best = 0.0
    best_norm = None
    best_gain = None
    for _repeat in range(16):
        factor = 0.5 * (low + high)
        gain = factor * floor
        if np.any(measured <= gain):
            high = factor
            continue
        bounds = np.concatenate([
            np.full(nodes, (measured[0] - gain[0]) * SCALE),
            np.full(nodes, (measured[1] - gain[1]) * SCALE),
        ])
        scaled_residual = residual * SCALE
        scaled_response = response * SCALE
        rows = scaled_response.shape[0]
        objective = np.concatenate([np.zeros(kept), np.ones(kept)])
        inequalities = np.zeros((2 * rows + 2 * kept, 2 * kept))
        inequalities[:rows, :kept] = scaled_response
        inequalities[rows:2 * rows, :kept] = -scaled_response
        inequalities[2 * rows:2 * rows + kept, :kept] = np.eye(kept)
        inequalities[2 * rows:2 * rows + kept, kept:] = -np.eye(kept)
        inequalities[2 * rows + kept:, :kept] = -np.eye(kept)
        inequalities[2 * rows + kept:, kept:] = -np.eye(kept)
        limits = np.concatenate([
            bounds - scaled_residual,
            bounds + scaled_residual,
            np.zeros(2 * kept),
        ])
        solved = linprog(
            objective, A_ub=inequalities, b_ub=limits, bounds=(None, None),
            method='highs',
            options={'primal_feasibility_tolerance': 1e-10,
                     'dual_feasibility_tolerance': 1e-10})
        if not solved.success:
            high = factor
            continue
        step = basis @ solved.x[:kept]
        predicted = I64.predicted_maxima(residual, matrix, step)
        realized = measured - predicted
        norm = float(np.linalg.norm(step))
        if norm <= 1e-6 and np.all(realized >= gain * 0.999):
            low = factor
            best = factor
            best_norm = norm
            best_gain = realized
        else:
            high = factor
    return {
        'kept_modes': kept,
        'largest_gain_factor': best,
        'step_norm': best_norm,
        'predicted_gain': None if best_gain is None else best_gain.tolist(),
    }


def scan():
    parent = json.loads(PARENT.read_text())
    if parent['profile_identity'] != PARENT_IDENTITY:
        raise ValueError('sup-floor parent must be the measured factor-200 history')
    if not parent['newton_step_accepted']:
        raise ValueError('sup-floor parent must be an accepted step')
    gradient = np.asarray(parent['action_gradient'], float)
    jacobian = np.asarray(parent['history_jacobian'], float)
    _delta, residual, matrix, _spectrum = I64.rank_truncated_delta(gradient, jacobian, 1)
    row = np.linalg.norm(matrix, axis=1)
    col = np.linalg.norm(matrix, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    _left, _spectrum, right = np.linalg.svd(matrix / row[:, None] / col[None, :], full_matrices=False)
    measured = np.asarray(parent['constraint_maxima'], float)
    floor = np.abs(np.asarray(parent['prediction_minus_measured'], float))
    rows = [largest_factor(matrix, residual, right, col, measured, floor, kept)
            for kept in (16, 32, 48)]
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-SUP-FLOOR-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: well-conditioned sup-norm gain is inside the last N prediction error',
        'stall': 'n64_primal_norm_sup_gain_inside_prediction_error',
        'parent_profile_identity': parent['profile_identity'],
        'parent_constraint_maxima': measured.tolist(),
        'last_prediction_error_abs': floor.tolist(),
        'existence_tolerance': EXISTENCE,
        'error_term_at_1e-8_can_reach_existence': bool(np.all(measured + 1e-8 <= EXISTENCE)),
        'blocks': rows,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'enclosed_error': None,
        'next_owner': 'declared (w,U) class on I=S(1)+[0.12,0.18]',
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {str(PARENT.relative_to(ROOT)): digest(PARENT)},
    }


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = scan()
    for key in (
        'stall', 'parent_constraint_maxima', 'last_prediction_error_abs',
        'error_term_at_1e-8_can_reach_existence', 'blocks',
    ):
        if saved[key] != fresh[key]:
            raise ValueError('sup-floor register changed: '+key)
    if saved['physical_EXISTENCE_certificate'] or saved['physical_NONEXISTENCE_certificate']:
        raise ValueError('sup floor is not a gate certificate')
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
            'last_prediction_error_abs': record['last_prediction_error_abs'],
            'blocks': record['blocks'],
            'error_term_at_1e-8_can_reach_existence': record['error_term_at_1e-8_can_reach_existence'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({'stall': saved['stall'], 'status': saved['status']}, indent=2))
