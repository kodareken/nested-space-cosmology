#!/usr/bin/env python3
"""Fas C cause for the value-only ruler. No new search step is evolved."""
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

OUTPUT = ROOT/'results/development/nsc-ks-gate-search-stall.json'
FLOOR = ROOT/'results/development/nsc-ks-gate-numerical-floor.json'
REANCHOR = ROOT/'results/development/nsc-ks-gate-reanchor.json'
OWNERS = (
    'scripts/derive_nsc_ks_gate_search_stall.py',
    'docs/nsc-ks-gate-search-stall.md',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    floor = json.loads(FLOOR.read_text())
    reanchor = json.loads(REANCHOR.read_text())
    if floor['loop_numerics'] is not None:
        raise ValueError('a declared loop setting would authorize a search')
    movement = float(floor['largest_one_knob_movement'])
    if movement <= 1e-10:
        raise ValueError('a stable floor would not be classified as numerical error')
    return {
        'schema': 'NSC-KS-GATE-SEARCH-STALL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'OPEN: value-only floor tracks the solver; search is not authorized',
        'stall_name': 'value_only_indicator_floor_above_1e-10',
        'fas_c_cause': 'numerical error',
        'next_owner': 'the same declared (w, U) class on I, after the ruler is fixed',
        'families_evolved': 0,
        'campaign_best_name': reanchor['campaign_best_name'],
        'campaign_best_merit': reanchor['campaign_best_merit'],
        'largest_one_knob_movement': movement,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'named_gap': 'value_only_indicator_floor_above_1e-10',
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(FLOOR.relative_to(ROOT)): digest(FLOOR),
            str(REANCHOR.relative_to(ROOT)): digest(REANCHOR),
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
            raise FileExistsError('search-stall register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('search-stall replay differs')
    print(json.dumps({
        'status': result['status'],
        'stall_name': result['stall_name'],
        'fas_c_cause': result['fas_c_cause'],
        'families_evolved': result['families_evolved'],
    }, indent=2))
