#!/usr/bin/env python3
"""Leave deciding error bounds None after the Fourier leftover floor."""
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
import derive_nsc_ks_deciding_error_budget as B
import derive_nsc_ks_flux_first_integral as F
import derive_nsc_ks_fourier_wu_evolve_refusal as R
import derive_nsc_ks_n16_leftover_anatomy as A

OUTPUT = ROOT/'results/development/nsc-ks-deciding-bounds-after-fourier.json'
OWNERS = (
    'scripts/derive_nsc_ks_deciding_bounds_after_fourier.py',
    'docs/nsc-ks-deciding-bounds-after-fourier.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    refusal = json.loads(R.OUTPUT.read_text())
    flux = json.loads(F.OUTPUT.read_text())
    budget = json.loads(B.OUTPUT.read_text())
    for register in (refusal, flux, budget):
        for path, expected in {**register.get('source_hashes', {}),
                               **register.get('input_hashes', {})}.items():
            if digest(path) != expected:
                raise ValueError('Fourier deciding-bounds input changed: '+path)
    leftover = np.asarray(refusal['all_node_unclipped_leftover'], float)
    near_existence = bool(np.any(leftover <= 100 * A.EXISTENCE_TOLERANCE))
    claimed_gap = bool(flux['physical_NONEXISTENCE_certificate'])
    if (near_existence or claimed_gap) and not flux['physical_NONEXISTENCE_certificate']:
        raise ValueError('a leftover near 3e-11 would require spending the missing bounds')
    spent = False
    return {
        'schema': 'NSC-KS-DECIDING-BOUNDS-AFTER-FOURIER-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: leftover stays on the residual scale; '
            'UV, full between-node, field and low/subgap remain None'),
        'unclipped_all_node_leftover': leftover.tolist(),
        'geometry_between_node_remainder': budget['geometry_between_node_remainder'],
        'full_between_node_remainder': None,
        'changed_history_UV_tail': None,
        'field_error_bound': None,
        'baseline_low_subgap': None,
        'bounds_spent': spent,
        'scope': {
            'physical_local_gate': 'OPEN',
            'missing_bounds_left_none': True,
            'missing_bounds_cannot_create_EXISTENCE_against_residual_scale_leftover': True,
            'new_field_or_source_evolutions': 0,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': flux['physical_NONEXISTENCE_certificate'],
            'metric_timestep': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(R.OUTPUT.relative_to(ROOT)): digest(R.OUTPUT),
            str(F.OUTPUT.relative_to(ROOT)): digest(F.OUTPUT),
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
            raise FileExistsError('Fourier deciding-bounds exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('Fourier deciding-bounds replay differs')
    print(json.dumps({
        'status': result['status'],
        'bounds_spent': result['bounds_spent'],
        'geometry_between_node_remainder': result['geometry_between_node_remainder'],
        'changed_history_UV_tail': result['changed_history_UV_tail'],
        'field_error_bound': result['field_error_bound'],
        'baseline_low_subgap': result['baseline_low_subgap'],
    }, indent=2))
