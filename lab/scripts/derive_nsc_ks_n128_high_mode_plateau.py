#!/usr/bin/env python3
"""Stable-rank ceiling of the measured n=128 high-mode block.

The width-128 linear programs in the column owner return no feasible point.
This register keeps the numerical-rank screen, where the singular-value
ratio stays below 1e6. It does not evolve a history.
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
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
from derive_nsc_ks_n128_high_modes import max_n_gain, predicted_maxima, whiten

SOURCE = ROOT/'results/development/nsc-ks-n128-high-modes.json'
PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust03.json'
OUTPUT = ROOT/'results/development/nsc-ks-n128-high-mode-plateau.json'
BAR = 1e-10
RATIO = 1e6
HS = (1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 2e-4)
OWNERS = (
    'scripts/derive_nsc_ks_n128_high_mode_plateau.py',
    'docs/nsc-ks-n128-high-mode-plateau.md',
    'scripts/derive_nsc_ks_n128_high_modes.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    saved = json.loads(SOURCE.read_text())
    parent = json.loads(PARENT.read_text())
    if saved['stall'] != 'n128_high_mode_linear_gain_below_1e-10':
        raise ValueError('plateau parent stall changed')
    if saved['profile_identity'] != '3af275b0ee0d0d0237799290e7a76639513d47cc5f1bc6da1d4d33593392f888':
        raise ValueError('plateau lifted identity changed')
    if saved['probe_authorized'] or saved['newton_step_accepted']:
        raise ValueError('plateau source must not authorize a step')
    jacobian = np.asarray(saved['history_jacobian'], float)
    gradient = np.asarray(parent['action_gradient'], float)
    measured = np.asarray(saved['best_constraint_maxima'], float)
    if jacobian.shape != (256, gradient.shape[0], 2):
        raise ValueError('n=128 Jacobian shape changed')
    residual = np.concatenate([gradient[:, 0], gradient[:, 1]])
    matrix = np.vstack([jacobian[:, :, 0].T, jacobian[:, :, 1].T])
    high = np.concatenate([np.arange(64, 128), np.arange(192, 256)])
    block = matrix[:, high]
    row = np.linalg.norm(block, axis=1)
    col = np.linalg.norm(block, axis=0)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    _left, spectrum, right = np.linalg.svd(block / row[:, None] / col[None, :], full_matrices=False)
    kept = np.where(spectrum[0] / spectrum < RATIO)[0]
    if kept.size == 0:
        raise ValueError('n=128 high-mode block has no stable singular value')
    last = int(kept[-1]) + 1
    basis = whiten(right, col, 0, last)
    rows = []
    for h in HS:
        gain, direction = max_n_gain(block @ basis, residual, measured, h, 0.0)
        beta = 0.0
        if direction is not None and gain > 0.0:
            predicted = predicted_maxima(residual, block, h * (basis @ direction))
            beta = float(measured[1] - predicted[1])
        rows.append({
            'h': float(h),
            'n_gain': float(gain),
            'beta_gain': beta,
            'clears_bar_with_beta_moving': bool(gain >= BAR and beta > 1e-16),
        })
    plateau = rows[-1]['n_gain'] == rows[-2]['n_gain'] and rows[-1]['n_gain'] < BAR
    if any(row['clears_bar_with_beta_moving'] for row in rows) or not plateau:
        raise ValueError('stable-rank screen no longer supports the n=128 stall')
    return {
        'schema': 'NSC-KS-N128-HIGH-MODE-PLATEAU-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: stable-rank n=128 N gain plateaus at 9.899343e-11 with beta flat',
        'stall': 'n128_high_mode_linear_gain_below_1e-10',
        'profile_identity': saved['profile_identity'],
        'parent_profile_identity': saved['parent_profile_identity'],
        'best_constraint_maxima': saved['best_constraint_maxima'],
        'stable_rank': last,
        'singular_ratio_cut': RATIO,
        'prune_bar_on_n_gain': BAR,
        'rows': rows,
        'plateau_n_gain': rows[-1]['n_gain'],
        'plateau_beta_gain': rows[-1]['beta_gain'],
        'probe_authorized': False,
        'quad_coeff_measured': False,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'higher_n_authorized': False,
        'higher_n_note': 'n=256 is the same subclass pattern and is not authorized by this plateau',
        'uv_field_low_subgap_between_node': None,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(SOURCE.relative_to(ROOT)): digest(SOURCE),
            str(PARENT.relative_to(ROOT)): digest(PARENT),
        },
    }


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = compute()
    for key in ('stall', 'stable_rank', 'rows', 'plateau_n_gain', 'plateau_beta_gain', 'probe_authorized'):
        if saved[key] != fresh[key]:
            raise ValueError('n=128 plateau register changed: '+key)
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        if OUTPUT.exists():
            raise SystemExit('plateau register already exists')
        payload = compute()
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'stall': payload['stall'],
            'stable_rank': payload['stable_rank'],
            'plateau_n_gain': payload['plateau_n_gain'],
            'rows': payload['rows'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({'stall': saved['stall'], 'plateau_n_gain': saved['plateau_n_gain']}, indent=2))
