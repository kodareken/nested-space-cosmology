#!/usr/bin/env python3
"""Record the freed-jet n=16 proposal and refuse evolution."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped as D
import derive_nsc_ks_n16_leftover_identification as ID
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-n16-freed-jet-proposal.json'
OWNERS = (
    'scripts/derive_nsc_ks_n16_freed_jet_proposal.py',
    'docs/nsc-ks-n16-freed-jet-proposal.md',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    identification = json.loads(ID.OUTPUT.read_text())
    for path, expected in {**identification['source_hashes'],
                           **identification['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('freed-jet proposal input changed: '+path)
    if identification['physical_NONEXISTENCE_certificate']:
        raise ValueError('freed-jet proposal is only for a failed NON-EXISTENCE branch')
    leftover = np.asarray(identification['dropped_jet_unclipped_all_node_leftover'], float)
    clipped = np.asarray(identification['dropped_jet_clipped_all_node_leftover'], float)
    if identification['evolve_authorized']:
        raise ValueError('this owner records a refused evolve; do not hide an authorized trial')
    if np.any(leftover <= ID.A.EXISTENCE_TOLERANCE):
        raise ValueError('a leftover that reaches 3e-11 would require evolution')
    trial = json.loads(ID.I4.OUTPUT.read_text())
    family = LocalIncomingFamily(np.asarray(trial['history']['coefficients'], float))
    return {
        'schema': 'NSC-KS-N16-FREED-JET-PROPOSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: freed-jet linear leftover stays about 1e-4; '
            'no family evolved'),
        'current_profile_identity': profile_identity(family),
        'unclipped_all_node_leftover': leftover.tolist(),
        'clipped_all_node_leftover': clipped.tolist(),
        'unclipped_step_norm': identification['dropped_jet_step_norm'],
        'clipped_step_norm': identification['dropped_jet_clipped_step_norm'],
        'existence_tolerance': ID.A.EXISTENCE_TOLERANCE,
        'evolve_authorized': False,
        'named_search_stall': 'n16_freed_jet_leftover_cannot_reach_existence',
        'forbidden_identities': list(D.FORBIDDEN_IDENTITIES),
        'scope': {
            'physical_local_gate': 'OPEN',
            'prediction_is_not_a_nonlinear_residual': True,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
            'metric_timestep': False,
            'full_scale_n16_not_evolved': True,
            'n8_next_not_evolved': True,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(ID.OUTPUT.relative_to(ROOT)): digest(ID.OUTPUT)},
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
            raise FileExistsError('freed-jet proposal exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('freed-jet proposal replay differs')
    print(json.dumps({
        'status': result['status'],
        'evolve_authorized': result['evolve_authorized'],
        'unclipped_all_node_leftover': result['unclipped_all_node_leftover'],
        'clipped_all_node_leftover': result['clipped_all_node_leftover'],
        'named_search_stall': result['named_search_stall'],
    }, indent=2))
