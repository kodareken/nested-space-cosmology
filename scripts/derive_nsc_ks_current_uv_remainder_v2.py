#!/usr/bin/env python3
"""Record/replay the changed-history UV remainder successor.

The v1 remainder pilot JSON and schema stay immutable. This successor extends
the M=4 owner: it represents A2, A3, n4, C4, C_M, current-history L0 A4 H2
integrals, the upstream higher-order H2 remainder, transported majors and the
vacuum N/beta tail, and binds those names to the evolved state law, current
profile, horizon preparation, source inventory and KS subtraction. Only
quantities available in a small authenticated control are computed.
Unavailable values remain null and OPEN. The E^{-2} cancellation and E^{-3}
first allowed paired order are preserved and are not tail bounds.
"""
import argparse
from collections.abc import Mapping
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.nsc_ks_evaluation_binding import write_bytes_atomic
from recursive_horizons.nsc_ks_finite_history_uv_remainder import (
    CAUCHY_RECORD, CUTOFF_BRIDGE_CONTROL, HISTORY_RECORD, INVENTORY_RECORD,
    PILOT_RECORD, PILOT_SCHEMA, PREPARATION_OWNER, RESTART_RECORD,
    STATE_LAW_OWNER, SUBTRACTION_OWNER, SUCCESSOR_SCHEMA, V2_RECORD,
    V3_RECORD, V4_RECORD, changed_history_uv_successor_report,
    frozen_coefficient_records, validate_changed_history_uv_successor)

OUTPUT = ROOT / 'results/development/nsc-ks-current-uv-remainder-v2.json'
INVENTORY_PAYLOAD = 'results/development/artifacts/nsc-ks-source-inventory.npz'
OWNERS = (
    'scripts/derive_nsc_ks_current_uv_remainder_v2.py',
    'src/recursive_horizons/nsc_ks_finite_history_uv_remainder.py',
    'tests/test_nsc_ks_finite_history_uv_remainder.py',
    'tests/test_nsc_ks_current_uv_remainder_v2.py',
    'docs/nsc-ks-current-uv-remainder.md',
)
INPUTS = (
    PILOT_RECORD,
    'scripts/derive_nsc_ks_current_uv_remainder_pilot.py',
    V2_RECORD, V3_RECORD, V4_RECORD, HISTORY_RECORD, CAUCHY_RECORD,
    RESTART_RECORD, INVENTORY_RECORD, INVENTORY_PAYLOAD, CUTOFF_BRIDGE_CONTROL,
    STATE_LAW_OWNER, PREPARATION_OWNER, SUBTRACTION_OWNER,
)


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def plain(value):
    if isinstance(value, Mapping):
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
    report = plain(changed_history_uv_successor_report(ROOT))
    frozen = plain(frozen_coefficient_records(ROOT))
    return {
        'schema': SUCCESSOR_SCHEMA,
        'status': report['status'],
        'supersedes_record_only': PILOT_RECORD,
        'supersedes_schema_only': PILOT_SCHEMA,
        'coefficient_records': frozen,
        'remainder': report,
        'first_noncancelling_paired_order': report['first_noncancelling_paired_order'],
        'leading_paired_e_minus2_cancels_for_equal_mu': report[
            'leading_paired_e_minus2_cancels_for_equal_mu'],
        'leading_e_minus2_cancellation_is_not_a_tail_bound': True,
        'e_minus3_first_allowed_order_is_not_a_tail_bound': True,
        'numerical_C4': None,
        'numerical_C_M': None,
        'vacuum_integrated_tail_N_beta': None,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'v1_pilot_rewritten': False,
        'v1_pilot_sha256': report['v1_pilot_sha256'],
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {path: digest(path) for path in INPUTS},
    }


def check(value):
    if value.get('schema') != SUCCESSOR_SCHEMA:
        raise ValueError('unexpected UV remainder successor schema')
    for path, expected in {**value['source_hashes'], **value['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('UV remainder successor dependency changed: ' + path)
    if value['v1_pilot_sha256'] != digest(PILOT_RECORD):
        raise ValueError('historical v1 remainder pilot bytes changed')
    if value['v1_pilot_rewritten'] or json.loads(
            (ROOT / PILOT_RECORD).read_text()).get('schema') != PILOT_SCHEMA:
        raise ValueError('historical v1 remainder pilot schema is immutable')
    if value['numerical_C4'] is not None or value['numerical_C_M'] is not None:
        raise ValueError('successor may not invent C4 or C_M')
    if value['vacuum_integrated_tail_N_beta'] is not None:
        raise ValueError('leading cancellation is not a vacuum-tail bound')
    if value['physical_EXISTENCE_certificate'] or value['physical_NONEXISTENCE_certificate']:
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if not value['leading_e_minus2_cancellation_is_not_a_tail_bound']:
        raise ValueError('E^{-2} cancellation is not a tail bound')
    if not value['e_minus3_first_allowed_order_is_not_a_tail_bound']:
        raise ValueError('E^{-3} first allowed order is not a tail bound')
    validate_changed_history_uv_successor(value['remainder'], ROOT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    computed = calculate()
    if args.check and not OUTPUT.is_file():
        raise FileNotFoundError(
            'UV remainder successor record is missing; --check cannot create evidence')
    if args.record:
        payload = (
            json.dumps(computed, indent=2, sort_keys=True, allow_nan=False) + '\n'
        ).encode()
        if OUTPUT.is_file():
            if OUTPUT.read_bytes() != payload:
                raise ValueError('existing UV remainder successor record differs')
        else:
            write_bytes_atomic(OUTPUT, payload, exclusive=True)
        value = computed
    else:
        value = json.loads(OUTPUT.read_text())
        if value != computed:
            raise ValueError('UV remainder successor replay differs')
    check(value)
    print(json.dumps({
        'status': value['status'],
        'schema': value['schema'],
        'first_noncancelling_paired_order': value['first_noncancelling_paired_order'],
        'numerical_C4': value['numerical_C4'],
        'numerical_C_M': value['numerical_C_M'],
        'vacuum_integrated_tail_N_beta': value['vacuum_integrated_tail_N_beta'],
        'v1_pilot_rewritten': value['v1_pilot_rewritten'],
        'quantity_names': sorted(value['remainder']['quantities']),
    }, indent=2))
