#!/usr/bin/env python3
"""Gate status after the assembled Chebyshev n=16 leftover floor."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
# Do not import next5 or the assembled-search status module: they pull Dirac.
import derive_nsc_ks_assembled_chebyshev_gn_stall as S
import derive_nsc_ks_deciding_bounds_after_assembled_chebyshev_gn as B
PREVIOUS = ROOT/'results/development/nsc-local-incoming-gate-status-after-assembled-search.json'
OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-after-assembled-chebyshev-gn.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_after_assembled_chebyshev_gn.py',
    'docs/nsc-local-incoming-gate-status-after-assembled-chebyshev-gn.md',
)
EXISTENCE_TOLERANCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def bind_register(register):
    for path, expected in {**register.get('source_hashes', {}),
                           **register.get('input_hashes', {})}.items():
        if digest(path) != expected:
            raise ValueError('Chebyshev-GN gate status input changed: '+path)


def compute():
    previous = json.loads(PREVIOUS.read_text())
    search = json.loads(S.OUTPUT.read_text())
    bounds = json.loads(B.OUTPUT.read_text())
    bind_register(previous)
    bind_register(search)
    bind_register(bounds)
    if search['evolve_authorized'] or search['unclipped_leftover_reaches_existence_tolerance']:
        raise ValueError('a leftover at 3e-11 is not a measured EXISTENCE certificate')
    if search['families_evolved'] != 0:
        raise ValueError('this status owner records the Dirac-free leftover floor')
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-AFTER-ASSEMBLED-CHEBYSHEV-GN-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: assembled Chebyshev leftover cannot reach EXISTENCE; '
            'finite image, not a class identity'),
        'state_law': previous['state_law'],
        'history_class': previous['history_class'],
        'interval': previous['interval'],
        'profile_identity': previous['profile_identity'],
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'assembled_unclipped_leftover': search['assembled_unclipped_leftover'],
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
            str(PREVIOUS.relative_to(ROOT)): digest(PREVIOUS),
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
            raise FileExistsError('Chebyshev-GN gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('Chebyshev-GN gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'named_search_stall': result['named_search_stall'],
        'next_named_owner': result['next_named_owner'],
        'assembled_unclipped_leftover': result['assembled_unclipped_leftover'],
        'families_evolved': result['families_evolved'],
    }, indent=2))
