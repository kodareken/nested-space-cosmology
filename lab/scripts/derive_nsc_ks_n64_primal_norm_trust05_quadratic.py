#!/usr/bin/env python3
"""Record the rejected five-times step whose N residual rose."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust03.json'
CHILD = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust05.json'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-trust05-quadratic.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_trust05_quadratic.py',
    'docs/nsc-ks-n64-primal-norm-trust05-quadratic.md',
)
EXISTENCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def record():
    parent = json.loads(PARENT.read_text())
    child = json.loads(CHILD.read_text())
    if not parent['newton_step_accepted']:
        raise ValueError('the shrink must stay the accepted parent')
    if child['newton_step_accepted'] or child['residual_improved_on_all_components']:
        raise ValueError('the five-times step must stay rejected')
    parent_max = parent['constraint_maxima']
    predicted = child['linear_predicted_all_node_residual_maxima']
    measured = child['constraint_maxima']
    if not (measured[0] > parent_max[0] and measured[1] < parent_max[1]):
        raise ValueError('the five-times step must raise N and lower beta')
    step_norm = float(child['step_norm'])
    n_error = abs(float(child['prediction_minus_measured'][0]))
    n_gain = float(parent_max[0] - predicted[0])
    beta_gain = float(parent_max[1] - predicted[1])
    if n_error <= n_gain or step_norm <= 0.0 or n_gain <= 0.0 or beta_gain <= 0.0:
        raise ValueError('the five-times step must have an N error above its predicted N gain')
    coefficient = n_error / step_norm**2
    slope = n_gain / step_norm
    break_even_norm = slope / coefficient
    break_even_n_gain = slope * break_even_norm
    if parent_max[0] - break_even_n_gain <= EXISTENCE:
        raise ValueError('break-even N gain would require a certification review')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-TRUST05-QUADRATIC-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: beta followed the Jacobian and N rose at five times the shrink norm',
        'stall': 'n64_primal_norm_trust05_n_rose',
        'best_profile_identity': parent['profile_identity'],
        'best_constraint_maxima': parent_max,
        'best_prediction_minus_measured': parent['prediction_minus_measured'],
        'rejected_profile_identity': child['profile_identity'],
        'rejected_constraint_maxima': measured,
        'linear_predicted_all_node_residual_maxima': predicted,
        'prediction_minus_measured': child['prediction_minus_measured'],
        'predicted_n_gain': n_gain,
        'predicted_beta_gain': beta_gain,
        'step_norm': step_norm,
        'singular_ratio': child['singular_ratio'],
        'n_error_power_used_for_authorization': child['n_error_power'],
        'direction_quadratic_n_coefficient': coefficient,
        'break_even_step_norm': break_even_norm,
        'break_even_predicted_n_gain': break_even_n_gain,
        'break_even_n_gain_can_reach_existence': False,
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
            raise ValueError('trust05 quadratic register changed: '+key)
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
            'best_constraint_maxima': payload['best_constraint_maxima'],
            'rejected_constraint_maxima': payload['rejected_constraint_maxima'],
            'direction_quadratic_n_coefficient': payload['direction_quadratic_n_coefficient'],
            'break_even_step_norm': payload['break_even_step_norm'],
            'break_even_predicted_n_gain': payload['break_even_predicted_n_gain'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({'stall': saved['stall'], 'status': saved['status']}, indent=2))
