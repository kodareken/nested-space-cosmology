#!/usr/bin/env python3
"""Record/replay the changed-history UV remainder v3 successor.

The v1 remainder pilot and v2 successor JSON stay immutable. This successor
owns current-history characteristic transports of major A2, A3 and n4 under
the switched state law, forms the production history-minus-reference e^{-3}
C4 definition, and records newly finite values only where proved. Numerical
C4, C_M and the vacuum tail remain OPEN. Runtime is outside the scientific
digest.
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
from recursive_horizons.nsc_ks_current_uv_report import (
    CAUCHY_RECORD, CUTOFF_BRIDGE_CONTROL, HISTORY_RECORD, INVENTORY_RECORD,
    PILOT_RECORD, PILOT_SCHEMA, PREPARATION_OWNER, RESTART_RECORD,
    STATE_LAW_OWNER, SUBTRACTION_OWNER, SUCCESSOR_SCHEMA, V2_RECORD,
    V2_REMAINDER_RECORD, V3_RECORD, V3_SUCCESSOR_SCHEMA, V4_RECORD,
    changed_history_uv_v3_report, frozen_coefficient_records,
    validate_changed_history_uv_v3)
from recursive_horizons.nsc_ks_current_uv_transport import transport_arrays
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = ROOT / 'results/development/nsc-ks-current-uv-remainder-v3-rebound.json'
PAYLOAD = ROOT / 'results/development/artifacts/nsc-ks-current-uv-remainder-v3.npz'
INVENTORY_PAYLOAD = 'results/development/artifacts/nsc-ks-source-inventory.npz'
OWNERS = (
    'src/recursive_horizons/nsc_ks_current_uv_report.py',
    'scripts/derive_nsc_ks_current_uv_remainder_v3.py',
    'src/recursive_horizons/nsc_ks_current_uv_transport.py',
    'src/recursive_horizons/nsc_ks_finite_history_uv_remainder.py',
    'tests/test_nsc_ks_current_uv_remainder_v3.py',
    'docs/nsc-ks-current-uv-remainder-v3.md',
)
INPUTS = (
    PILOT_RECORD,
    V2_REMAINDER_RECORD,
    'scripts/derive_nsc_ks_current_uv_remainder_pilot.py',
    'scripts/derive_nsc_ks_current_uv_remainder_v2.py',
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
    report = plain(changed_history_uv_v3_report(ROOT))
    frozen = plain(frozen_coefficient_records(ROOT))
    arrays = transport_arrays(report['transport'])
    payload = deterministic_npz_bytes(arrays)
    return {
        'schema': V3_SUCCESSOR_SCHEMA,
        'status': report['status'],
        'supersedes_record_only': V2_REMAINDER_RECORD,
        'supersedes_schema_only': SUCCESSOR_SCHEMA,
        'historical_v1_schema': PILOT_SCHEMA,
        'coefficient_records': frozen,
        'remainder': report,
        'proved': report['proved'],
        'next_primitive': report['next_primitive'],
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
        'v2_record_rewritten': False,
        'v1_pilot_sha256': report['v1_pilot_sha256'],
        'v2_record_sha256': report['v2_record_sha256'],
        'payload': {
            'path': 'results/development/artifacts/nsc-ks-current-uv-remainder-v3.npz',
            'sha256': sha256(payload).hexdigest(),
            'bytes': len(payload),
        },
        'budget_implication': report['budget_implication'],
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {path: digest(path) for path in INPUTS},
    }, payload


def check(value, payload=None):
    if value.get('schema') != V3_SUCCESSOR_SCHEMA:
        raise ValueError('unexpected UV remainder v3 schema')
    for path, expected in {**value['source_hashes'], **value['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('UV remainder v3 dependency changed: ' + path)
    if value['v1_pilot_sha256'] != digest(PILOT_RECORD):
        raise ValueError('historical v1 remainder pilot bytes changed')
    if value['v2_record_sha256'] != digest(V2_REMAINDER_RECORD):
        raise ValueError('historical v2 remainder successor bytes changed')
    if value['v1_pilot_rewritten'] or json.loads(
            (ROOT / PILOT_RECORD).read_text()).get('schema') != PILOT_SCHEMA:
        raise ValueError('historical v1 remainder pilot schema is immutable')
    if value['v2_record_rewritten'] or json.loads(
            (ROOT / V2_REMAINDER_RECORD).read_text()).get('schema') != SUCCESSOR_SCHEMA:
        raise ValueError('historical v2 remainder successor schema is immutable')
    if value['numerical_C4'] is not None or value['numerical_C_M'] is not None:
        raise ValueError('v3 successor may not invent C4 or C_M')
    if value['vacuum_integrated_tail_N_beta'] is not None:
        raise ValueError('leading cancellation is not a vacuum-tail bound')
    if value['physical_EXISTENCE_certificate'] or value['physical_NONEXISTENCE_certificate']:
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if not value['leading_e_minus2_cancellation_is_not_a_tail_bound']:
        raise ValueError('E^{-2} cancellation is not a tail bound')
    if not value['e_minus3_first_allowed_order_is_not_a_tail_bound']:
        raise ValueError('E^{-3} first allowed order is not a tail bound')
    if payload is not None:
        if sha256(payload).hexdigest() != value['payload']['sha256']:
            raise ValueError('v3 UV remainder payload digest differs')
        if len(payload) != value['payload']['bytes']:
            raise ValueError('v3 UV remainder payload length differs')
    validate_changed_history_uv_v3(value['remainder'], ROOT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    computed, payload = calculate()
    if args.check and (not OUTPUT.is_file() or not PAYLOAD.is_file()):
        raise FileNotFoundError(
            'UV remainder v3 record or payload is missing; --check cannot create evidence')
    if args.record:
        encoded = (
            json.dumps(computed, indent=2, sort_keys=True, allow_nan=False) + '\n'
        ).encode()
        # Rebound record deliberately reuses the immutable original payload.
        if PAYLOAD.is_file():
            if PAYLOAD.read_bytes() != payload:
                raise ValueError('existing UV remainder v3 payload differs')
        else:
            write_bytes_atomic(PAYLOAD, payload, exclusive=True)
        if OUTPUT.is_file():
            if OUTPUT.read_bytes() != encoded:
                raise ValueError('existing UV remainder v3 record differs')
        else:
            write_bytes_atomic(OUTPUT, encoded, exclusive=True)
        value = computed
    else:
        value = json.loads(OUTPUT.read_text())
        recorded_payload = PAYLOAD.read_bytes()
        if value != computed:
            raise ValueError('UV remainder v3 replay differs')
        if recorded_payload != payload:
            raise ValueError('UV remainder v3 payload replay differs')
        payload = recorded_payload
    check(value, payload)
    print(json.dumps({
        'status': value['status'],
        'schema': value['schema'],
        'first_noncancelling_paired_order': value['first_noncancelling_paired_order'],
        'proved': value['proved'],
        'numerical_C4': value['numerical_C4'],
        'numerical_C_M': value['numerical_C_M'],
        'vacuum_integrated_tail_N_beta': value['vacuum_integrated_tail_N_beta'],
        'next_primitive': value['next_primitive'],
        'v1_pilot_rewritten': value['v1_pilot_rewritten'],
        'v2_record_rewritten': value['v2_record_rewritten'],
        'payload_sha256': value['payload']['sha256'],
    }, indent=2))
