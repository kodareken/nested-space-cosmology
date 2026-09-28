#!/usr/bin/env python3
"""Retract s^5 and the class-wide stall name. Fourier bytes stay historical."""
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
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_local_incoming_gate_status_after_fourier as G

OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-retract-s5.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_retract_s5.py',
    'docs/nsc-local-incoming-gate-status-retract-s5.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    previous = json.loads(G.OUTPUT.read_text())
    for path, expected in {**previous['source_hashes'], **previous['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('s^5 retraction input changed: '+path)
    if previous['physical_EXISTENCE_certificate'] or previous['physical_NONEXISTENCE_certificate']:
        raise ValueError('retraction is only for an OPEN Fourier status')
    if 's^5' not in (previous.get('next_named_owner') or ''):
        raise ValueError('expected the overclaiming s^5 next owner to retract')
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-RETRACT-S5-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: finite Fourier leftover is representation evidence, '
            'not NON-EXISTENCE of the predeclared class; s^5 is retracted'),
        'state_law': previous['state_law'],
        'history_class': previous['history_class'],
        'interval': previous['interval'],
        'profile_identity': previous['profile_identity'],
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'fourier_unclipped_all_node_leftover': previous['fourier_unclipped_all_node_leftover'],
        'previous_named_search_stall': previous['named_search_stall'],
        'previous_next_named_owner': previous['next_named_owner'],
        'named_search_stall': 'finite_wu_geometry_image_cannot_reach_existence',
        'next_named_owner': (
            'Fourier leftover-null, then the continuous cokernel of the '
            'declared (w,U) geometry map; no class rewrite'),
        's5_retracted': True,
        'class_rewrite_is_not_next': True,
        'finite_leftover_is_not_class_identity': True,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'scope': {
            'open_is_not_done': True,
            'manuscript_and_public_release_blocked': True,
            'arxiv_blocked': True,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'historical_fourier_bytes_unchanged': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(G.OUTPUT.relative_to(ROOT)): digest(G.OUTPUT),
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
            raise FileExistsError('s^5 retraction exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('s^5 retraction replay differs')
    print(json.dumps({
        'status': result['status'],
        'named_search_stall': result['named_search_stall'],
        'next_named_owner': result['next_named_owner'],
        's5_retracted': result['s5_retracted'],
    }, indent=2))
