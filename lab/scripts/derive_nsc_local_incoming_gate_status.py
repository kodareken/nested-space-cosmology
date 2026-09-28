#!/usr/bin/env python3
"""Record the local incoming-gate status: OPEN named stall, not a certificate."""
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
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
import derive_nsc_ks_coupled_newton_n16_damped_next5 as N5
import derive_nsc_ks_coupled_newton_n16_jacobian_refresh as R
import derive_nsc_ks_deciding_error_budget as B

OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status.py',
    'docs/nsc-local-incoming-gate-status.md',
)
STALL = 'n16_damped_frozen_jacobian_leftover_cannot_reach_existence'


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    trial = json.loads(I4.OUTPUT.read_text())
    next5 = json.loads(N5.OUTPUT.read_text())
    refresh = json.loads(R.OUTPUT.read_text())
    budget = json.loads(B.OUTPUT.read_text())
    for register in (trial, next5, refresh, budget):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('gate status input changed: '+path)
    if refresh['named_search_stall'] != STALL or next5['named_search_stall'] != STALL:
        raise ValueError('gate status requires the leftover stall')
    if refresh['scope']['physical_NONEXISTENCE_claimed'] or budget['scope']['physical_NONEXISTENCE_claimed']:
        raise ValueError('named stall must not be promoted to NON-EXISTENCE')
    leftover = refresh['unclipped_linear_predicted_all_node_residual_maxima']
    measured = trial['constraint_maxima']
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: named n=16 leftover stall; no physical gate certificate',
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; original fixed source',
        'history_class': 'local source-fixed pure-radius compatible histories',
        'interval': 'S(1)+[.12,.18]',
        'profile_identity': trial['profile_identity'],
        'measured_constraint_maxima': measured,
        'existence_tolerance': N5.EXISTENCE_TOLERANCE,
        'unclipped_linear_leftover': leftover,
        'named_search_stall': STALL,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'error_budget': {
            'geometry_between_node_remainder': budget['geometry_between_node_remainder'],
            'full_between_node_remainder': None,
            'changed_history_UV_tail': None,
            'field_error_bound': None,
            'baseline_low_subgap': None,
        },
        'next_named_owner': (
            'a necessary relation for the declared switched class, '
            'or a changed Jacobian/representation that can reach 3e-11'),
        'scope': {
            'solver_step_is_not_physical_EXISTENCE': True,
            'failed_optimization_is_not_NONEXISTENCE': True,
            'open_is_not_done': True,
            'manuscript_and_public_release_blocked': True,
            'arxiv_blocked': True,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
            str(N5.OUTPUT.relative_to(ROOT)): digest(N5.OUTPUT),
            str(R.OUTPUT.relative_to(ROOT)): digest(R.OUTPUT),
            str(B.OUTPUT.relative_to(ROOT)): digest(B.OUTPUT),
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
            raise FileExistsError('local incoming-gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('local incoming-gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'physical_EXISTENCE_certificate': result['physical_EXISTENCE_certificate'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'named_search_stall': result['named_search_stall'],
        'measured_constraint_maxima': result['measured_constraint_maxima'],
        'unclipped_linear_leftover': result['unclipped_linear_leftover'],
        'next_named_owner': result['next_named_owner'],
    }, indent=2))
