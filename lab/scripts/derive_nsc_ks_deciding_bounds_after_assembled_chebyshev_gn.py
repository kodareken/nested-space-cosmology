#!/usr/bin/env python3
"""Leave deciding bounds None after the assembled Chebyshev leftover floor."""
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
# Do not import iterate4/next5 or the deciding-budget module: they pull Dirac.
import derive_nsc_ks_assembled_chebyshev_gn_stall as S
BUDGET = ROOT/'results/development/nsc-ks-deciding-error-budget.json'
OUTPUT = ROOT/'results/development/nsc-ks-deciding-bounds-after-assembled-chebyshev-gn.json'
OWNERS = (
    'scripts/derive_nsc_ks_deciding_bounds_after_assembled_chebyshev_gn.py',
    'docs/nsc-ks-deciding-bounds-after-assembled-chebyshev-gn.md',
)
EXISTENCE_TOLERANCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def bind_register(register):
    for path, expected in {**register.get('source_hashes', {}),
                           **register.get('input_hashes', {})}.items():
        if digest(path) != expected:
            raise ValueError('Chebyshev-GN deciding-bounds input changed: '+path)


def compute():
    search = json.loads(S.OUTPUT.read_text())
    budget = json.loads(BUDGET.read_text())
    bind_register(search)
    bind_register(budget)
    leftover = np.asarray(search['assembled_unclipped_leftover'], float)
    if search['evolve_authorized'] or np.any(leftover <= 100 * EXISTENCE_TOLERANCE):
        raise ValueError('a leftover near 3e-11 would require spending the missing bounds')
    return {
        'schema': 'NSC-KS-DECIDING-BOUNDS-AFTER-ASSEMBLED-CHEBYSHEV-GN-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: assembled Chebyshev leftover stays residual-scale; '
            'UV, full between-node, field and low/subgap remain None'),
        'unclipped_all_node_leftover': leftover.tolist(),
        'geometry_between_node_remainder': budget['geometry_between_node_remainder'],
        'full_between_node_remainder': None,
        'changed_history_UV_tail': None,
        'field_error_bound': None,
        'baseline_low_subgap': None,
        'bounds_spent': False,
        'scope': {
            'physical_local_gate': 'OPEN',
            'missing_bounds_left_none': True,
            'missing_bounds_cannot_create_EXISTENCE_against_residual_scale_leftover': True,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(S.OUTPUT.relative_to(ROOT)): digest(S.OUTPUT),
            str(BUDGET.relative_to(ROOT)): digest(BUDGET),
        },
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
            raise FileExistsError('Chebyshev-GN deciding-bounds exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('Chebyshev-GN deciding-bounds replay differs')
    print(json.dumps({
        'status': result['status'],
        'bounds_spent': result['bounds_spent'],
        'changed_history_UV_tail': result['changed_history_UV_tail'],
        'field_error_bound': result['field_error_bound'],
        'baseline_low_subgap': result['baseline_low_subgap'],
    }, indent=2))
