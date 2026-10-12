#!/usr/bin/env python3
"""Successor low/subgap decision after the existing group12/32 refinements."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_source_inventory import RetainedSourceInventory
from recursive_horizons.nsc_ks_upstream_error import classify_low_subgap_windows

OUTPUT = ROOT / 'results/development/nsc-ks-low-subgap-decision-v2.json'
INPUTS = (
    'results/development/nsc-ks-low-subgap-decision.json',
    'results/development/nsc-ks-source-inventory.json',
    'results/development/nsc-incoming-source-update-v5.json',
    'results/development/nsc-incoming-retained-order24-bound.json',
    'results/development/nsc-incoming-group32-order24-bound.json',
    'results/development/nsc-incoming-subgap-source.json',
    'results/development/nsc-incoming-low-energy-vacuum-enclosure.json',
)
OWNERS = (
    'scripts/derive_nsc_ks_low_subgap_decision_v2.py',
    'tests/test_nsc_ks_low_subgap_decision_v2.py',
    'docs/nsc-ks-low-subgap-decision-v2.md',
)


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def calculate():
    inventory_record = json.loads((ROOT / INPUTS[1]).read_text())
    v5 = json.loads((ROOT / INPUTS[2]).read_text())
    group12 = json.loads((ROOT / INPUTS[3]).read_text())
    group32 = json.loads((ROOT / INPUTS[4]).read_text())
    with RetainedSourceInventory(ROOT) as inventory:
        classified = classify_low_subgap_windows(
            inventory_record, inventory.positive_panels(), v5['coverage'])
        provenance = {item['path']: item['sha256'] for item in inventory.provenance}
    remaining = [row for row in classified['windows']
                 if row.get('source_reconstruction') == 'unbounded'
                 and row.get('baseline_action_region') == 'not_covered']
    return {
        'schema': 'NSC-KS-LOW-SUBGAP-DECISION-v2',
        'status': (
            'OPEN: refined group12/group32 action windows reused; physical low/subgap '
            'rho=1 column reconstruction and transport remain unbounded'),
        'supersedes_decision_only': INPUTS[0],
        'v5_covered_action_bound_N_beta':
            v5['partial_error_budget']['covered_spectral_regions_action_error_upper'],
        'reused_successors': {
            'group12_order24_lapse_error_upper':
                group12['energy']['lapse_action_error_upper'],
            'group32_combined_action_error_upper':
                group32['windows']['combined']['conditional_constraint_action_error_upper'],
            'group32_validated_lower_window': v5['coverage']['32']['validated_lower_window'],
        },
        'historical_coarse_4p998e_minus9_is_current_remaining_bound': False,
        'coverage': classified,
        'remaining_uncovered_reconstruction_windows': remaining,
        'remaining_uncovered_reconstruction_window_count': len(remaining),
        'remaining_uncovered_reconstruction_rows': sum(row['n_rows'] for row in remaining),
        'remaining_low_subgap_source_error_bound': None,
        'next_method': (
            'directed rho=1 physical-column reconstruction errors on panels below each '
            'V5 joined threshold, followed by correlated transport through the final history'),
        'redo_group12_or_group32_order24_middle_bounds': False,
        'physical_EXISTENCE_certificate': False,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {path: digest(path) for path in INPUTS} | provenance,
    }


def check(value):
    if value.get('schema') != 'NSC-KS-LOW-SUBGAP-DECISION-v2':
        raise ValueError('unexpected low/subgap successor schema')
    for path, expected in {**value['source_hashes'], **value['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('low/subgap successor dependency changed: ' + path)
    if value['remaining_low_subgap_source_error_bound'] is not None:
        raise ValueError('successor may not invent a missing physical-column bound')
    if value['historical_coarse_4p998e_minus9_is_current_remaining_bound']:
        raise ValueError('superseded coarse middle estimate cannot be current low/subgap error')
    if value['redo_group12_or_group32_order24_middle_bounds']:
        raise ValueError('already-certified middle windows must be reused')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = calculate() if args.record else json.loads(OUTPUT.read_text())
    if args.record:
        publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)),
            (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode())
    check(value)
    print(json.dumps({key: value[key] for key in (
        'status', 'reused_successors', 'remaining_uncovered_reconstruction_window_count',
        'remaining_uncovered_reconstruction_rows',
        'remaining_low_subgap_source_error_bound', 'next_method')}, indent=2))
