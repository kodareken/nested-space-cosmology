#!/usr/bin/env python3
"""Record the assembled Chebyshev n=32 finite-image stall. No new Dirac."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
N32_OUTPUT = ROOT/'results/development/nsc-ks-coupled-newton-n32-assembled.json'
OUTPUT = ROOT/'results/development/nsc-ks-assembled-chebyshev-n32-stall.json'
OWNERS = (
    'scripts/derive_nsc_ks_assembled_chebyshev_n32_stall.py',
    'docs/nsc-ks-assembled-chebyshev-n32-stall.md',
)
EXISTENCE_TOLERANCE = 3e-11
PARENT_RESIDUAL = (0.004107193193038788, 0.005790788097514318)
EXPECTED_PADDED = '49223e613233a84ece27eb4efd544eb462b4343453d7dcbc8f2cf0cfec9dff68'
EXPECTED_PARENT = 'e15982e241bdc2244481491f38878b501b34ad3bbb208a54f40415d862482898'


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    trial = json.loads(N32_OUTPUT.read_text())
    if trial['profile_identity'] != EXPECTED_PADDED:
        raise ValueError('n=32 stall must sit on the zero-padded iterate4 history')
    if trial['parent_profile_identity'] != EXPECTED_PARENT:
        raise ValueError('n=32 stall parent must be the accepted iterate4 identity')
    if trial['completed_positive_families'] != 60 or trial['new_operator_solves'] != 60:
        raise ValueError('n=32 assembled must evolve all sixty high-mode families')
    leftover = trial['assembled_unclipped_leftover']
    measured = trial['constraint_maxima']
    if measured != list(PARENT_RESIDUAL):
        raise ValueError('n=32 residual must equal the iterate4 residual at the same g')
    if trial['evolve_authorized'] or trial['selected_proposal'] is not None:
        raise ValueError('finite-image stall must not authorize a trial')
    if trial['assembled_condition_number'] < 1e14:
        raise ValueError('expected an ill-conditioned n=32 assembled map')
    return {
        'schema': 'NSC-KS-ASSEMBLED-CHEBYSHEV-N32-STALL-v1',
        'accountable_author': 'Douglas Ek',
        'status': trial['status'],
        'profile_identity': trial['profile_identity'],
        'parent_profile_identity': trial['parent_profile_identity'],
        'constraint_maxima': measured,
        'assembled_unclipped_leftover': leftover,
        'assembled_clipped_leftover': trial['assembled_clipped_leftover'],
        'assembled_condition_number': trial['assembled_condition_number'],
        'assembled_unclipped_step_norm': trial['assembled_unclipped_step_norm'],
        'unclipped_leftover_reaches_existence_tolerance': bool(
            leftover[0] <= EXISTENCE_TOLERANCE and leftover[1] <= EXISTENCE_TOLERANCE),
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'evolve_authorized': False,
        'families_evolved': trial['families_evolved'],
        'high_mode_directions_evolved': trial['high_mode_directions_evolved'],
        'low_mode_directions_reused': trial['low_mode_directions_reused'],
        'named_search_stall': trial['named_search_stall'],
        'fas_c_cause': (
            'n32 assembled finite image; condition ~1e14; linear step unavailable; '
            'leftover equals residual and stays above 3e-11'),
        'slot_transfer': False,
        'geometry_only_map': False,
        'invented_matter_columns': False,
        's5_retracted': True,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(N32_OUTPUT.relative_to(ROOT)): digest(N32_OUTPUT)},
        'reproducer': 'python3 scripts/derive_nsc_ks_assembled_chebyshev_n32_stall.py --check',
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
        previous = json.loads(OUTPUT.read_text())
        if previous != result:
            raise ValueError('n=32 assembled stall replay differs')
    print(json.dumps({
        'status': result['status'],
        'assembled_unclipped_leftover': result['assembled_unclipped_leftover'],
        'assembled_condition_number': result['assembled_condition_number'],
        'evolve_authorized': result['evolve_authorized'],
        'families_evolved': result['families_evolved'],
        'named_search_stall': result['named_search_stall'],
    }, indent=2))
