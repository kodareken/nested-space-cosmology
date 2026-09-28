#!/usr/bin/env python3
"""Refuse the n=32 retarded evolve after a 1e-4 geometry leftover."""
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
import derive_nsc_ks_n16_leftover_anatomy as A
import derive_nsc_ks_n32_geometry_leftover as G

OUTPUT = ROOT/'results/development/nsc-ks-n32-evolve-refusal.json'
OWNERS = (
    'scripts/derive_nsc_ks_n32_evolve_refusal.py',
    'docs/nsc-ks-n32-evolve-refusal.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    leftover = json.loads(G.OUTPUT.read_text())
    for path, expected in {**leftover['source_hashes'], **leftover['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('n=32 evolve refusal input changed: '+path)
    if leftover['evolve_authorized']:
        raise ValueError('this owner records only a refused n=32 evolve')
    unclipped = np.asarray(leftover['all_node_unclipped_leftover'], float)
    if np.any(unclipped <= A.EXISTENCE_TOLERANCE):
        raise ValueError('a leftover that reaches 3e-11 would require evolution')
    return {
        'schema': 'NSC-KS-N32-EVOLVE-REFUSAL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: n=32 retarded evolve refused; leftover cannot reach 3e-11',
        'evolve_authorized': False,
        'families_evolved': 0,
        'all_node_unclipped_leftover': leftover['all_node_unclipped_leftover'],
        'all_node_condition_number': leftover['all_node_condition_number'],
        'named_search_stall': leftover['named_search_stall'],
        'existence_tolerance': A.EXISTENCE_TOLERANCE,
        'scope': leftover['scope'],
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(G.OUTPUT.relative_to(ROOT)): digest(G.OUTPUT)},
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
            raise FileExistsError('n=32 evolve refusal exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('n=32 evolve refusal replay differs')
    print(json.dumps({
        'status': result['status'],
        'evolve_authorized': result['evolve_authorized'],
        'families_evolved': result['families_evolved'],
        'named_search_stall': result['named_search_stall'],
    }, indent=2))
