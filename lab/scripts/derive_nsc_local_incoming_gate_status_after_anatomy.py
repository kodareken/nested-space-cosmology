#!/usr/bin/env python3
"""Record the post-anatomy local incoming-gate status: still OPEN."""
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
import derive_nsc_ks_deciding_bounds_after_leftover as B
import derive_nsc_ks_n16_freed_jet_proposal as P
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_ks_n16_leftover_identification as ID

OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-after-anatomy.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_after_anatomy.py',
    'docs/nsc-local-incoming-gate-status-after-anatomy.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    anatomy = json.loads(A.OUTPUT.read_text())
    identification = json.loads(ID.OUTPUT.read_text())
    proposal = json.loads(P.OUTPUT.read_text())
    bounds = json.loads(B.OUTPUT.read_text())
    for register in (anatomy, identification, proposal, bounds):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('post-anatomy gate status input changed: '+path)
    if identification['physical_NONEXISTENCE_certificate'] or proposal['evolve_authorized']:
        raise ValueError('status owner is only for the refused certificate / refused evolve')
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-AFTER-ANATOMY-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: leftover classified; no physical gate certificate',
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; original fixed source',
        'history_class': 'local source-fixed pure-radius compatible histories',
        'interval': 'S(1)+[.12,.18]',
        'profile_identity': anatomy['profile_identity'],
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'solve_unclipped_leftover': anatomy['solve_unclipped_leftover'],
        'freed_jet_unclipped_all_node_leftover': proposal['unclipped_all_node_leftover'],
        'walk_projection_sign_stable': False,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_search_stall': proposal['named_search_stall'],
        'next_named_owner': (
            'a representation or necessary relation outside the n=16 '
            'jet-constrained and freed-jet linear leftovers'),
        'scope': {
            'open_is_not_done': True,
            'manuscript_and_public_release_blocked': True,
            'arxiv_blocked': True,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(A.OUTPUT.relative_to(ROOT)): digest(A.OUTPUT),
            str(ID.OUTPUT.relative_to(ROOT)): digest(ID.OUTPUT),
            str(P.OUTPUT.relative_to(ROOT)): digest(P.OUTPUT),
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
            raise FileExistsError('post-anatomy gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('post-anatomy gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'physical_EXISTENCE_certificate': result['physical_EXISTENCE_certificate'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'named_search_stall': result['named_search_stall'],
        'freed_jet_unclipped_all_node_leftover': result['freed_jet_unclipped_all_node_leftover'],
        'next_named_owner': result['next_named_owner'],
    }, indent=2))
