#!/usr/bin/env python3
"""Record the rejected beta-primary step whose N residual rose."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-beta48.json'
CHILD = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-beta48b.json'
OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n64-primal-norm-beta-nonlinearity.json'
OWNERS = (
    'scripts/derive_nsc_ks_n64_primal_norm_beta_nonlinearity.py',
    'docs/nsc-ks-n64-primal-norm-beta-nonlinearity.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def record():
    parent = json.loads(PARENT.read_text())
    child = json.loads(CHILD.read_text())
    if child['newton_step_accepted']:
        raise ValueError('the nonlinear beta step must stay rejected')
    if child['residual_improved_on_all_components']:
        raise ValueError('the nonlinear beta step must not be an improvement on both components')
    return {
        'schema': 'NSC-KS-N64-PRIMAL-NORM-BETA-NONLINEARITY-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: beta followed the Jacobian and N rose',
        'stall': 'n64_primal_norm_beta_step_n_nonlinearity',
        'best_profile_identity': parent['profile_identity'],
        'best_constraint_maxima': parent['constraint_maxima'],
        'rejected_profile_identity': child['profile_identity'],
        'rejected_constraint_maxima': child['constraint_maxima'],
        'prediction_minus_measured': child['prediction_minus_measured'],
        'step_norm': child['step_norm'],
        'singular_ratio': child['singular_ratio'],
        'newton_step_accepted': False,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'enclosed_error': None,
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
        'prediction_minus_measured', 'newton_step_accepted',
    ):
        if saved[key] != fresh[key]:
            raise ValueError('beta-nonlinearity register changed: '+key)
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        payload = record()
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'stall': payload['stall'],
            'best_constraint_maxima': payload['best_constraint_maxima'],
            'rejected_constraint_maxima': payload['rejected_constraint_maxima'],
            'prediction_minus_measured': payload['prediction_minus_measured'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({'stall': saved['stall'], 'status': saved['status']}, indent=2))
