#!/usr/bin/env python3
"""Gate status after the assembled Chebyshev n=32 finite-image stall."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREVIOUS = ROOT/'results/development/nsc-local-incoming-gate-status-after-assembled-chebyshev-gn.json'
STALL = ROOT/'results/development/nsc-ks-assembled-chebyshev-n32-stall.json'
BOUNDS = ROOT/'results/development/nsc-ks-deciding-bounds-after-assembled-chebyshev-n32.json'
OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-after-assembled-chebyshev-n32.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_after_assembled_chebyshev_n32.py',
    'docs/nsc-local-incoming-gate-status-after-assembled-chebyshev-n32.md',
)
EXISTENCE_TOLERANCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def bind_register(register):
    for path, expected in {**register.get('source_hashes', {}),
                           **register.get('input_hashes', {})}.items():
        if digest(path) != expected:
            raise ValueError('n=32 gate status input changed: '+path)


def compute():
    previous = json.loads(PREVIOUS.read_text())
    stall = json.loads(STALL.read_text())
    bounds = json.loads(BOUNDS.read_text())
    bind_register(previous)
    bind_register(stall)
    bind_register(bounds)
    if stall['evolve_authorized'] or stall['unclipped_leftover_reaches_existence_tolerance']:
        raise ValueError('a finite-image leftover is not a measured EXISTENCE certificate')
    if bounds['epsilon_spend_authorized']:
        raise ValueError('deciding epsilon spend must stay off for this leftover gap')
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-AFTER-ASSEMBLED-CHEBYSHEV-N32-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: assembled Chebyshev n=32 leftover cannot reach EXISTENCE; '
            'finite image, not a class identity'),
        'state_law': previous['state_law'],
        'history_class': previous['history_class'],
        'interval': previous['interval'],
        'profile_identity': stall['profile_identity'],
        'parent_profile_identity': stall['parent_profile_identity'],
        'existence_tolerance': EXISTENCE_TOLERANCE,
        'assembled_unclipped_leftover': stall['assembled_unclipped_leftover'],
        'assembled_condition_number': stall['assembled_condition_number'],
        'families_evolved': stall['families_evolved'],
        'evolve_authorized': False,
        's5_retracted': True,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_search_stall': stall['named_search_stall'],
        'next_named_owner': 'same declared (w,U) class on I; no class rewrite',
        'fas_c_cause': stall['fas_c_cause'],
        'scope': {
            'open_is_not_done': True,
            'manuscript_and_public_release_blocked': True,
            'arxiv_blocked': True,
            'public_027_blocked': True,
            'new_field_or_source_evolutions': 60,
            'high_mode_only_evolutions': True,
            'metric_timestep': False,
            'no_class_rewrite': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(PREVIOUS.relative_to(ROOT)): digest(PREVIOUS),
            str(STALL.relative_to(ROOT)): digest(STALL),
            str(BOUNDS.relative_to(ROOT)): digest(BOUNDS),
        },
        'reproducer': (
            'python3 scripts/derive_nsc_local_incoming_gate_status_after_assembled_chebyshev_n32.py '
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
            raise ValueError('n=32 gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'assembled_unclipped_leftover': result['assembled_unclipped_leftover'],
        'evolve_authorized': result['evolve_authorized'],
        'families_evolved': result['families_evolved'],
        'named_search_stall': result['named_search_stall'],
        'physical_EXISTENCE_certificate': result['physical_EXISTENCE_certificate'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
    }, indent=2))
