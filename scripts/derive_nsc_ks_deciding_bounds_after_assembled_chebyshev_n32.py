#!/usr/bin/env python3
"""Deciding-budget status after n=32 assembled finite-image stall. No Dirac."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STALL = ROOT/'results/development/nsc-ks-assembled-chebyshev-n32-stall.json'
PREV_BOUNDS = ROOT/'results/development/nsc-ks-deciding-bounds-after-assembled-chebyshev-gn.json'
OUTPUT = ROOT/'results/development/nsc-ks-deciding-bounds-after-assembled-chebyshev-n32.json'
OWNERS = (
    'scripts/derive_nsc_ks_deciding_bounds_after_assembled_chebyshev_n32.py',
    'docs/nsc-ks-deciding-bounds-after-assembled-chebyshev-n32.md',
)
EXISTENCE_TOLERANCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    stall = json.loads(STALL.read_text())
    previous = json.loads(PREV_BOUNDS.read_text())
    for path, expected in previous['source_hashes'].items():
        if digest(path) != expected:
            raise ValueError('previous deciding-bounds input changed: '+path)
    leftover = stall['assembled_unclipped_leftover']
    gap = min(leftover) - EXISTENCE_TOLERANCE
    # Geometry between-node 1e-13 cannot move a gap of order 1e-3.
    # UV/field/low-subgap/full between-node stay None.
    return {
        'schema': 'NSC-KS-DECIDING-BOUNDS-AFTER-ASSEMBLED-CHEBYSHEV-N32-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: n=32 leftover gap stays far above 3e-11; '
            'no deciding epsilon spend'),
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'assembled_unclipped_leftover': leftover,
        'gap_above_existence_tolerance': gap,
        'geometry_between_node_bound': 1e-13,
        'geometry_between_node_can_move_gap': False,
        'changed_history_UV_bound': None,
        'full_between_node_bound': None,
        'field_error_bound': None,
        'low_subgap_error_bound': None,
        'epsilon_spend_authorized': False,
        'named_search_stall': stall['named_search_stall'],
        'physical_local_gate': 'OPEN',
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(STALL.relative_to(ROOT)): digest(STALL),
            str(PREV_BOUNDS.relative_to(ROOT)): digest(PREV_BOUNDS),
        },
        'reproducer': (
            'python3 scripts/derive_nsc_ks_deciding_bounds_after_assembled_chebyshev_n32.py '
            '--check'),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.record:
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    else:
        if json.loads(OUTPUT.read_text()) != result:
            raise ValueError('n=32 deciding-bounds replay differs')
    print(json.dumps({
        'status': result['status'],
        'gap_above_existence_tolerance': result['gap_above_existence_tolerance'],
        'epsilon_spend_authorized': result['epsilon_spend_authorized'],
        'changed_history_UV_bound': result['changed_history_UV_bound'],
        'full_between_node_bound': result['full_between_node_bound'],
        'field_error_bound': result['field_error_bound'],
        'low_subgap_error_bound': result['low_subgap_error_bound'],
    }, indent=2))
