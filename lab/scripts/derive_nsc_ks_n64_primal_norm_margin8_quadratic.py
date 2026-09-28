#!/usr/bin/env python3
"""Record the rejected ball step whose N residual rose."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-beta48.json'
CHILD = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-margin8.json'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-margin8-quadratic.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_margin8_quadratic.py',
    'docs/nsc-ks-n64-primal-norm-margin8-quadratic.md',
)
EXISTENCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def record():
    parent = json.loads(PARENT.read_text())
    child = json.loads(CHILD.read_text())
    if child['newton_step_accepted'] or child['residual_improved_on_all_components']:
        raise ValueError('the ball step must stay rejected')
    if child['profile_identity'] == parent['profile_identity']:
        raise ValueError('the rejected history must differ from the best point')
    parent_max = parent['constraint_maxima']
    predicted = child['linear_predicted_all_node_residual_maxima']
    measured = child['constraint_maxima']
    if not (measured[0] > parent_max[0] and measured[1] < parent_max[1]):
        raise ValueError('the ball step must raise N and lower beta')
    step_norm = float(child['step_norm'])
    n_error = abs(float(child['prediction_minus_measured'][0]))
    n_gain = float(parent_max[0] - predicted[0])
    beta_gain = float(parent_max[1] - predicted[1])
    if n_error <= n_gain or step_norm <= 0.0 or n_gain <= 0.0 or beta_gain <= 0.0:
        raise ValueError('the ball step must have an N error above its predicted N gain')
    coefficient = n_error / step_norm**2
    parent_coefficient = abs(float(parent['prediction_minus_measured'][0])) / float(parent['step_norm'])**2
    n_slope = n_gain / step_norm
    beta_slope = beta_gain / step_norm
    break_even_norm = n_slope / coefficient
    break_even_n_gain = n_slope * break_even_norm
    break_even_beta_gain = beta_slope * break_even_norm
    if parent_max[0] - break_even_n_gain > EXISTENCE:
        certification_moved = False
    else:
        raise ValueError('break-even N gain would require a certification review')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-MARGIN8-QUADRATIC-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: beta followed the Jacobian and N rose inside the 6.5e-9 ball',
        'stall': 'n64_primal_norm_margin8_n_quadratic',
        'best_profile_identity': parent['profile_identity'],
        'best_constraint_maxima': parent_max,
        'rejected_profile_identity': child['profile_identity'],
        'rejected_constraint_maxima': measured,
        'linear_predicted_all_node_residual_maxima': predicted,
        'prediction_minus_measured': child['prediction_minus_measured'],
        'predicted_n_gain': n_gain,
        'predicted_beta_gain': beta_gain,
        'step_norm': step_norm,
        'singular_ratio': child['singular_ratio'],
        'parent_quadratic_n_coefficient': parent_coefficient,
        'direction_quadratic_n_coefficient': coefficient,
        'direction_coefficient_over_parent': coefficient / parent_coefficient,
        'break_even_step_norm': break_even_norm,
        'break_even_predicted_n_gain': break_even_n_gain,
        'break_even_predicted_beta_gain': break_even_beta_gain,
        'break_even_n_gain_can_reach_existence': certification_moved,
        'newton_step_accepted': False,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'enclosed_error': None,
        'rank_74_and_above': 'unevolved',
        'next_owner': 'declared (w,U) class on I=S(1)+[0.12,0.18]',
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(PARENT.relative_to(ROOT)): digest(PARENT),
            str(CHILD.relative_to(ROOT)): digest(CHILD),
        },
    }


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = record()
    for key in (
        'stall', 'best_constraint_maxima', 'rejected_constraint_maxima',
        'prediction_minus_measured', 'direction_quadratic_n_coefficient',
        'break_even_step_norm', 'break_even_predicted_n_gain',
        'break_even_n_gain_can_reach_existence', 'newton_step_accepted',
    ):
        if saved[key] != fresh[key]:
            raise ValueError('margin quadratic register changed: '+key)
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        if OUTPUT.exists():
            raise SystemExit('quadratic register already exists')
        payload = record()
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'stall': payload['stall'],
            'rejected_constraint_maxima': payload['rejected_constraint_maxima'],
            'direction_quadratic_n_coefficient': payload['direction_quadratic_n_coefficient'],
            'direction_coefficient_over_parent': payload['direction_coefficient_over_parent'],
            'break_even_step_norm': payload['break_even_step_norm'],
            'break_even_predicted_n_gain': payload['break_even_predicted_n_gain'],
            'break_even_n_gain_can_reach_existence': payload['break_even_n_gain_can_reach_existence'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({'stall': saved['stall'], 'status': saved['status']}, indent=2))
