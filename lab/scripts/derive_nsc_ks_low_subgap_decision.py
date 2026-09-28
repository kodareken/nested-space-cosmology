#!/usr/bin/env python3
"""Record the reusable and still-unbounded low/subgap source windows."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_source_inventory import RetainedSourceInventory
from recursive_horizons.nsc_ks_upstream_error import (
    REMAINING_N_AFTER_KNOWN, SCHEMA_DECISION, classify_low_subgap_windows)

OUTPUT = ROOT / 'results/development/nsc-ks-low-subgap-decision.json'
INVENTORY = 'results/development/nsc-ks-source-inventory.json'
V5 = 'results/development/nsc-incoming-source-update-v5.json'
COARSE = 'results/development/nsc-incoming-energy-error-budget.json'
OWNERS = (
    'scripts/derive_nsc_ks_low_subgap_decision.py',
    'src/recursive_horizons/nsc_ks_upstream_error.py',
    'tests/test_nsc_ks_upstream_error.py',
    'docs/nsc-ks-low-subgap-decision.md',
)


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def calculate():
    inventory_record = json.loads((ROOT / INVENTORY).read_text())
    v5 = json.loads((ROOT / V5).read_text())
    coarse = json.loads((ROOT / COARSE).read_text())
    with RetainedSourceInventory(ROOT) as inventory:
        classified = classify_low_subgap_windows(
            inventory_record, inventory.positive_panels(), v5['coverage'])
        source_inputs = {item['path']: item['sha256'] for item in inventory.provenance}
    coarse_n = coarse['aggregate'][
        'conditional_order24_group14_plus_other31_and_thermal_lapse_bound']
    return {
        'schema': SCHEMA_DECISION,
        'status': (
            'OPEN: covered baseline regions and ell=0 zero response reusable; '
            'remaining low/subgap source reconstruction has no directed bound'),
        'coverage': classified,
        'covered_region_action_bound_N_beta':
            v5['partial_error_budget']['covered_spectral_regions_action_error_upper'],
        'existing_conditional_coarse_lapse_bound': coarse_n,
        'remaining_N_after_known_components': REMAINING_N_AFTER_KNOWN,
        'coarse_bound_fits_remaining_N': coarse_n <= REMAINING_N_AFTER_KNOWN,
        'dominant_existing_coarse_groups': [12, 32, 14],
        'decision': (
            'derive window-specific source reconstruction bounds for groups12/32 '
            'and complete group14; do not scale the unchanged coarse enclosure'),
        'remaining_low_subgap_source_error_bound': None,
        'physical_EXISTENCE_certificate': False,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {INVENTORY: digest(INVENTORY), V5: digest(V5),
                         COARSE: digest(COARSE), **source_inputs},
    }


def check(value):
    if value.get('schema') != SCHEMA_DECISION:
        raise ValueError('unexpected low/subgap decision schema')
    for path, expected in {**value['source_hashes'], **value['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('low/subgap dependency changed: ' + path)
    if value['remaining_low_subgap_source_error_bound'] is not None:
        raise ValueError('decision record may not invent a missing source bound')
    if value['coarse_bound_fits_remaining_N']:
        raise ValueError('existing coarse enclosure must not be promoted into the gate')


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
        'status', 'existing_conditional_coarse_lapse_bound',
        'remaining_N_after_known_components', 'coarse_bound_fits_remaining_N',
        'remaining_low_subgap_source_error_bound')}, indent=2))
