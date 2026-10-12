#!/usr/bin/env python3
"""Gate status after the assembled Fourier EXISTENCE search."""
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
import derive_nsc_ks_deciding_bounds_after_assembled_search as B
import derive_nsc_ks_fourier_assembled_search as S
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_local_incoming_gate_status_declared_class as G

OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-after-assembled-search.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_after_assembled_search.py',
    'docs/nsc-local-incoming-gate-status-after-assembled-search.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    previous = json.loads(G.OUTPUT.read_text())
    search = json.loads(S.OUTPUT.read_text())
    bounds = json.loads(B.OUTPUT.read_text())
    for register in (previous, search, bounds):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('assembled-search gate status input changed: '+path)
    existence = bool(search['evolve_authorized'] and search['unclipped_reaches_existence'])
    if existence:
        raise ValueError('a leftover at 3e-11 is not a measured EXISTENCE certificate')
    if search['families_evolved'] != 0:
        raise ValueError('this status owner records the Dirac-free search')
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-AFTER-ASSEMBLED-SEARCH-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: assembled Fourier leftover cannot reach EXISTENCE; '
            'finite image, not a class identity'),
        'state_law': previous['state_law'],
        'history_class': previous['history_class'],
        'interval': previous['interval'],
        'profile_identity': previous['profile_identity'],
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'assembled_unclipped_leftover': search['all_node_unclipped_leftover'],
        'assembled_transfer_faithful': search['assembled_transfer_faithful'],
        'families_evolved': search['families_evolved'],
        's5_retracted': True,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_search_stall': 'finite_wu_geometry_image_cannot_reach_existence',
        'next_named_owner': 'same declared (w,U) class on I; no class rewrite',
        'fas_c_cause': search['fas_c_cause'],
        'scope': {
            'open_is_not_done': True,
            'manuscript_and_public_release_blocked': True,
            'arxiv_blocked': True,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'no_class_rewrite': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(G.OUTPUT.relative_to(ROOT)): digest(G.OUTPUT),
            str(S.OUTPUT.relative_to(ROOT)): digest(S.OUTPUT),
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
            raise FileExistsError('assembled-search gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('assembled-search gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'named_search_stall': result['named_search_stall'],
        'next_named_owner': result['next_named_owner'],
        'assembled_unclipped_leftover': result['assembled_unclipped_leftover'],
        'families_evolved': result['families_evolved'],
    }, indent=2))
