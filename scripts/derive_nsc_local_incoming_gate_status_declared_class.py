#!/usr/bin/env python3
"""Gate status for the predeclared class after leftover-null and cokernel."""
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
import derive_nsc_ks_deciding_bounds_after_cokernel as B
import derive_nsc_ks_declared_wu_cokernel as C
import derive_nsc_ks_fourier_leftover_null as N
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_local_incoming_gate_status_retract_s5 as R

OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-declared-class.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_declared_class.py',
    'docs/nsc-local-incoming-gate-status-declared-class.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    retract = json.loads(R.OUTPUT.read_text())
    leftover_null = json.loads(N.OUTPUT.read_text())
    cokernel = json.loads(C.OUTPUT.read_text())
    bounds = json.loads(B.OUTPUT.read_text())
    for register in (retract, leftover_null, cokernel, bounds):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('declared-class gate status input changed: '+path)
    existence = False
    nonexistence = bool(
        leftover_null['physical_NONEXISTENCE_certificate']
        or cokernel['physical_NONEXISTENCE_certificate'])
    closed = existence or nonexistence
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-DECLARED-CLASS-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'scoped NON-EXISTENCE of the predeclared (w,U) class'
            if nonexistence else
            'OPEN: finite leftover is not a class identity; no class rewrite'),
        'state_law': retract['state_law'],
        'history_class': retract['history_class'],
        'interval': retract['interval'],
        'profile_identity': retract['profile_identity'],
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'leftover_null_walk_sign_stable': leftover_null['walk_sign_stable'],
        'continuous_cokernel_trivial': cokernel['continuous_cokernel_trivial'],
        'is_L_or_flux_again': cokernel['is_L_or_flux_again'],
        's5_retracted': True,
        'physical_local_gate': 'NON-EXISTENCE' if nonexistence else 'OPEN',
        'physical_EXISTENCE_certificate': existence,
        'physical_NONEXISTENCE_certificate': nonexistence,
        'named_search_stall': (
            None if closed else 'finite_wu_geometry_image_cannot_reach_existence'),
        'next_named_owner': (
            None if closed else
            'same declared (w,U) class on I; no class rewrite'),
        'scope': {
            'open_is_not_done': not closed,
            'manuscript_and_public_release_blocked': not closed,
            'arxiv_blocked': not closed,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'no_class_rewrite': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(R.OUTPUT.relative_to(ROOT)): digest(R.OUTPUT),
            str(N.OUTPUT.relative_to(ROOT)): digest(N.OUTPUT),
            str(C.OUTPUT.relative_to(ROOT)): digest(C.OUTPUT),
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
            raise FileExistsError('declared-class gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('declared-class gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'named_search_stall': result['named_search_stall'],
        'next_named_owner': result['next_named_owner'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
    }, indent=2))
