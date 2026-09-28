#!/usr/bin/env python3
"""Record the post-L / n=32 local incoming-gate status: still OPEN."""
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
import derive_nsc_ks_deciding_bounds_after_n32 as B
import derive_nsc_ks_n16_geometry_killing_current as K
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_ks_n32_evolve_refusal as R

OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-after-n32.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_after_n32.py',
    'docs/nsc-local-incoming-gate-status-after-n32.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    killing = json.loads(K.OUTPUT.read_text())
    refusal = json.loads(R.OUTPUT.read_text())
    bounds = json.loads(B.OUTPUT.read_text())
    for register in (killing, refusal, bounds):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('post-n32 gate status input changed: '+path)
    if killing['physical_NONEXISTENCE_certificate'] or refusal['evolve_authorized']:
        raise ValueError('status owner is only for the refused certificate / refused evolve')
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-AFTER-N32-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: L(E) and n=32 geometry leftover classified; no physical gate certificate',
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; original fixed source',
        'history_class': 'local source-fixed pure-radius compatible histories',
        'interval': 'S(1)+[.12,.18]',
        'profile_identity': killing['profile_identity'],
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'L_walk_sign_stable': killing['walk_sign_stable'],
        'L_walk_min': killing['walk_L_min'],
        'L_walk_max': killing['walk_L_max'],
        'n32_unclipped_all_node_leftover': refusal['all_node_unclipped_leftover'],
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_search_stall': refusal['named_search_stall'],
        'next_named_owner': (
            'a representation or necessary relation outside Chebyshev n<=32 '
            'of the declared s w + s^3 U/6 class, or a locked current other '
            'than leftover v and L(E) on this walk'),
        'scope': {
            'open_is_not_done': True,
            'manuscript_and_public_release_blocked': True,
            'arxiv_blocked': True,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(K.OUTPUT.relative_to(ROOT)): digest(K.OUTPUT),
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
            raise FileExistsError('post-n32 gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('post-n32 gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'physical_EXISTENCE_certificate': result['physical_EXISTENCE_certificate'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'named_search_stall': result['named_search_stall'],
        'next_named_owner': result['next_named_owner'],
    }, indent=2))
