#!/usr/bin/env python3
"""Record the local incoming-gate after Fourier leftover and flux currents."""
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
import derive_nsc_ks_deciding_bounds_after_fourier as B
import derive_nsc_ks_flux_first_integral as F
import derive_nsc_ks_fourier_wu_evolve_refusal as R
import derive_nsc_ks_n16_leftover_anatomy as A

OUTPUT = ROOT/'results/development/nsc-local-incoming-gate-status-after-fourier.json'
OWNERS = (
    'scripts/derive_nsc_local_incoming_gate_status_after_fourier.py',
    'docs/nsc-local-incoming-gate-status-after-fourier.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    flux = json.loads(F.OUTPUT.read_text())
    refusal = json.loads(R.OUTPUT.read_text())
    bounds = json.loads(B.OUTPUT.read_text())
    for register in (flux, refusal, bounds):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('post-Fourier gate status input changed: '+path)
    existence = False
    nonexistence = bool(flux['physical_NONEXISTENCE_certificate'])
    closed = existence or nonexistence
    stall = None if closed else 'declared_wu_image_cannot_reach_existence'
    return {
        'schema': 'NSC-LOCAL-INCOMING-GATE-STATUS-AFTER-FOURIER-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'scoped NON-EXISTENCE: flux current locks the declared (w,U) class'
            if nonexistence else
            'OPEN: Fourier leftover and flux currents classified; no physical gate certificate'),
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; original fixed source',
        'history_class': 'local source-fixed pure-radius compatible histories',
        'interval': 'S(1)+[.12,.18]',
        'profile_identity': flux['profile_identity'],
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'fourier_unclipped_all_node_leftover': refusal['all_node_unclipped_leftover'],
        'fourier_condition_number': refusal['all_node_condition_number'],
        'walk_integral_sign_stable': flux['walk_integral_sign_stable'],
        'walk_divergence_sign_stable': flux['walk_divergence_sign_stable'],
        'physical_local_gate': 'NON-EXISTENCE' if nonexistence else 'OPEN',
        'physical_EXISTENCE_certificate': existence,
        'physical_NONEXISTENCE_certificate': nonexistence,
        'named_search_stall': stall,
        'next_named_owner': (
            None if closed else
            'an extra normal power s^5 V(z) as a new undeclared class, not this (w,U) discretization'),
        'scope': {
            'open_is_not_done': not closed,
            'manuscript_and_public_release_blocked': not closed,
            'arxiv_blocked': not closed,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(F.OUTPUT.relative_to(ROOT)): digest(F.OUTPUT),
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
            raise FileExistsError('post-Fourier gate status exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('post-Fourier gate status replay differs')
    print(json.dumps({
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'physical_EXISTENCE_certificate': result['physical_EXISTENCE_certificate'],
        'physical_NONEXISTENCE_certificate': result['physical_NONEXISTENCE_certificate'],
        'named_search_stall': result['named_search_stall'],
        'next_named_owner': result['next_named_owner'],
    }, indent=2))
