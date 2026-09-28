#!/usr/bin/env python3
"""Value-only merits of the saved campaign histories.

Loop numerics were not declared: the one-knob floor stays above 1e-10.
These rows use the base solver at 47 nodes. They are measured residuals,
not an EXISTENCE certificate.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_gate_numerical_floor as F

OUTPUT = ROOT/'results/development/nsc-ks-gate-reanchor.json'
POINTS = (
    ('trial', 'results/development/nsc-ks-coupled-newton-n32-truncated-trial.json',
     'results/development/nsc-ks-gate-value-reanchor-trial-n47.json'),
    ('iterate3', 'results/development/nsc-ks-coupled-newton-n32-truncated-iterate3.json',
     'results/development/nsc-ks-gate-value-reanchor-iterate3-n47.json'),
    ('iterate4', 'results/development/nsc-ks-coupled-newton-n32-truncated-iterate4.json',
     'results/development/nsc-ks-gate-value-reanchor-iterate4-n47.json'),
    ('iterate5', 'results/development/nsc-ks-coupled-newton-n32-truncated-iterate5.json',
     'results/development/nsc-ks-gate-value-reanchor-iterate5-n47.json'),
    ('iterate6', 'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json',
     'results/development/nsc-ks-gate-value-iterate6-47.json'),
    ('iterate7', 'results/development/nsc-ks-coupled-newton-n64-truncated-iterate7.json',
     'results/development/nsc-ks-gate-value-reanchor-iterate7-n47.json'),
    ('iterate8', 'results/development/nsc-ks-coupled-newton-n64-truncated-iterate8.json',
     'results/development/nsc-ks-gate-value-reanchor-iterate8-n47.json'),
)
OWNERS = (
    'scripts/derive_nsc_ks_gate_reanchor.py',
    'docs/nsc-ks-gate-reanchor.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    floor = json.loads(F.OUTPUT.read_text())
    if floor['loop_numerics'] is not None or floor['certificate_numerics'] is not None:
        raise ValueError('re-anchor expected undeclared loop numerics')
    rows = []
    pending = []
    for name, tangent_path, value_path in POINTS:
        tangent = json.loads((ROOT/tangent_path).read_text())
        value_file = ROOT/value_path
        if not value_file.exists():
            pending.append(name)
            continue
        value = json.loads(value_file.read_text())
        if value['profile_identity'] != tangent['profile_identity']:
            raise ValueError(name + ' value-only identity does not match the saved history')
        rows.append({
            'name': name,
            'profile_identity': value['profile_identity'],
            'tangent_mode_constraint_maxima': tangent['constraint_maxima'],
            'value_only_constraint_maxima': value['constraint_maxima'],
            'value_only_merit': value['merit'],
            'nodes': value['node_count'],
        })
    best = None if not rows else min(rows, key=lambda row: row['value_only_merit'])
    return {
        'schema': 'NSC-KS-GATE-REANCHOR-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: value-only merits at the base solver; loop numerics were not declared'
            if not pending else
            'OPEN: value-only re-anchor is incomplete'),
        'solver': 'base value-only: grid 64, degree 32, rtol 2e-11, 47 nodes',
        'loop_numerics': None,
        'rows': rows,
        'pending': pending,
        'campaign_best_name': None if best is None else best['name'],
        'campaign_best_merit': None if best is None else best['value_only_merit'],
        'campaign_best_profile_identity': None if best is None else best['profile_identity'],
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_gap': floor['named_gap'],
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {str(F.OUTPUT.relative_to(ROOT)): digest(F.OUTPUT)},
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
            raise FileExistsError('re-anchor register exists; use --check')
        if result['pending']:
            raise ValueError('re-anchor is missing ' + result['pending'][0])
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('re-anchor replay differs')
    print(json.dumps({
        'status': result['status'],
        'pending': result['pending'],
        'campaign_best_name': result['campaign_best_name'],
        'campaign_best_merit': result['campaign_best_merit'],
    }, indent=2))
