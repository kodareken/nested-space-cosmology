#!/usr/bin/env python3
"""Record the N drift between two saved evolutions of one history.

Both sums are the iterate6 history b4b53f0b…. One carried 64 tangent
directions. The other froze that history and carried 64 further modes.
No family is evolved. The difference is an indicator, not an error bound
and not a physical residual.
"""
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
from recursive_horizons.nsc_frozen_sources import (
    selected_source_sha256,selected_source_context,frozen_source_context)
from recursive_horizons.nsc_ks_campaign_inputs import require_campaign_raw_inputs
OUTPUT = ROOT/'results/development/nsc-ks-gate-ruler-drift.json'
ITERATE6 = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n32-truncated-iterate6'
HIGH = ROOT/'results/development/artifacts/nsc-ks-coupled-newton-n64-high-modes'
IDENTITY = 'b4b53f0bba7290cfde14af741c82fda9af1b49a13f7a3b207a7230f39ea106e9'
HISTORICAL_OWNERS = (
    'scripts/derive_nsc_ks_gate_ruler_drift.py',
    'docs/nsc-ks-gate-ruler-drift.md',
)
OWNERS=tuple('docs/active-code-map.md' if name.startswith('docs/') else name
             for name in HISTORICAL_OWNERS)+('src/recursive_horizons/nsc_ks_campaign_inputs.py',)


def owned_sources():
    return HISTORICAL_OWNERS if selected_source_context() is not None else OWNERS



def digest(path):
    path = Path(path)
    relative=path.relative_to(ROOT).as_posix() if path.is_absolute() else path.as_posix()
    return selected_source_sha256(ROOT,relative)


def load_sum(directory):
    require_campaign_raw_inputs(ROOT,directory)
    files = sorted(Path(directory).glob('family-*.json'))
    if len(files) != 60:
        raise ValueError(f'{directory.name} does not have 60 family records')
    total = None
    relatives = []
    changes = []
    for path in files:
        record = json.loads(path.read_text())
        if record['profile_identity'] != IDENTITY:
            raise ValueError(f'{path.name} is not the iterate6 history')
        payload = ROOT/record['payload']['path']
        if digest(payload) != record['payload']['sha256']:
            raise ValueError(f'{path.name} payload hash changed')
        with np.load(payload, allow_pickle=False) as handle:
            change = np.asarray(handle['paired_matter_change'], float)
        if change.shape[1] != 2:
            raise ValueError(f'{path.name} matter change is not (nodes, 2)')
        changes.append(change)
        total = change.copy() if total is None else total + change
    return total, changes


def component_relative(reference, other):
    rows = []
    for left, right in zip(reference, other):
        scale = np.maximum(np.max(np.abs(left), axis=0), 1e-30)
        rows.append(np.max(np.abs(right - left), axis=0) / scale)
    rows = np.asarray(rows, float)
    return {
        'median': np.median(rows, axis=0).tolist(),
        'maximum': np.max(rows, axis=0).tolist(),
    }


def compute(*,source_context=None):
    if source_context is not None:
        with frozen_source_context(source_context):
            return compute()
    iterate6, changes6 = load_sum(ITERATE6)
    high, changes_high = load_sum(HIGH)
    if iterate6.shape != high.shape:
        raise ValueError('the two matter sums do not share a node grid')
    difference = high - iterate6
    maxima = np.max(np.abs(difference), axis=0)
    relative = component_relative(changes6, changes_high)
    return {
        'schema': 'NSC-KS-GATE-RULER-DRIFT-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: same-history matter sums differ by 5.78e-5 in N; '
            'N below that scale is not reproducible on these numerics'),
        'profile_identity': IDENTITY,
        'components': ['N', 'beta'],
        'node_count': int(iterate6.shape[0]),
        'positive_families': 60,
        'sum_difference_max_abs': maxima.tolist(),
        'per_family_relative': relative,
        'indicator_not_a_bound': True,
        'drift_subtracted_as_physics': False,
        'families_evolved': 0,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_gap': 'value_only_evaluator_then_numerical_floor',
        'source_hashes': {p: digest(p) for p in owned_sources()},
        'input_directories': [
            str(ITERATE6.relative_to(ROOT)),
            str(HIGH.relative_to(ROOT)),
        ],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    parser.add_argument('--source-context',help='full original source commit for read-only --check')
    args = parser.parse_args()
    if args.source_context and not args.check:parser.error('--source-context requires --check')
    result = compute(source_context=args.source_context)
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('ruler-drift register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('ruler-drift replay differs')
    print(json.dumps({
        'status': result['status'],
        'sum_difference_max_abs': result['sum_difference_max_abs'],
        'per_family_relative': result['per_family_relative'],
        'families_evolved': result['families_evolved'],
    }, indent=2))
