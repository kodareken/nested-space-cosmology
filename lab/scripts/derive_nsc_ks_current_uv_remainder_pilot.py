#!/usr/bin/env python3
"""Record/replay the current finite-history UV order and conditional M=4 gap."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_finite_history_uv_remainder import (
    CAUCHY_RECORD, HISTORY_RECORD, RESTART_RECORD,
    V2_RECORD, V3_RECORD, V4_RECORD, current_history_remainder_pilot,
    finite_history_uv_remainder_report, frozen_coefficient_records)

OUTPUT = ROOT / 'results/development/nsc-ks-current-uv-remainder-pilot.json'
INVENTORY_RECORD = 'results/development/nsc-ks-source-inventory.json'
OWNERS = (
    'scripts/derive_nsc_ks_current_uv_remainder_pilot.py',
    'src/recursive_horizons/nsc_ks_finite_history_uv_remainder.py',
    'tests/test_nsc_ks_finite_history_uv_remainder.py',
    'docs/nsc-ks-current-uv-remainder.md',
)
INPUTS = (V2_RECORD, V3_RECORD, V4_RECORD, HISTORY_RECORD,
          CAUCHY_RECORD, RESTART_RECORD, INVENTORY_RECORD)


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def plain(value):
    if isinstance(value, dict) or hasattr(value, 'items'):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def calculate():
    report = plain(finite_history_uv_remainder_report())
    pilot = plain(current_history_remainder_pilot(ROOT))
    frozen = plain(frozen_coefficient_records(ROOT))
    return {
        'schema': 'NSC-KS-CURRENT-UV-REMAINDER-PILOT-v1',
        'status': (
            'OPEN: first noncancelling paired order e^-3 and M=4 defect verified; '
            'thermal tail bounded; numerical C4/C_M and vacuum tail remain missing'),
        'coefficient_records': frozen,
        'remainder_proof': report,
        'current_history_pilot': pilot,
        'first_noncancelling_paired_order': 3,
        'thermal_tail_below_provisional_allocation': all(
            row['thermal_below_allocation'] for row in pilot['provisional_comparison']),
        'numerical_C4': None,
        'numerical_C_M': None,
        'vacuum_integrated_tail_N_beta': None,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {path: digest(path) for path in INPUTS},
        'next_required_primitive': {
            'upstream_affine_H2_remainder_orders_gt4': None,
            'current_history_L0_A4_H2_integrals': None,
            'transported_major_A2_A3_and_n4': None,
        },
    }


def check(value):
    if value.get('schema') != 'NSC-KS-CURRENT-UV-REMAINDER-PILOT-v1':
        raise ValueError('unexpected UV remainder pilot schema')
    for path, expected in {**value['source_hashes'], **value['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('UV remainder dependency changed: ' + path)
    if value['numerical_C4'] is not None or value['numerical_C_M'] is not None:
        raise ValueError('pilot may not invent C4 or C_M')
    if value['vacuum_integrated_tail_N_beta'] is not None:
        raise ValueError('leading cancellation is not a vacuum-tail bound')


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
        'status', 'first_noncancelling_paired_order',
        'thermal_tail_below_provisional_allocation', 'numerical_C4',
        'numerical_C_M', 'vacuum_integrated_tail_N_beta',
        'next_required_primitive')}, indent=2))
