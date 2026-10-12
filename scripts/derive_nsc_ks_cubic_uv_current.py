#!/usr/bin/env python3
"""Record/replay the generated massless e^{-3} N,beta coefficient.

Historical UV remainder and invariance JSON stay immutable. This successor
records the actual Pauli contraction of the generated A3/n4 transport on the
live 0b0e4ced history for original group 1. Diagnostic quadrature is not a
remainder-validated C4. C_M and the local gate remain OPEN.
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
from recursive_horizons.nsc_ks_current_uv_transport import (
    CAUCHY_RECORD,
    HISTORY_RECORD,
)
from recursive_horizons.nsc_ks_evaluation_binding import (
    implementation_hashes,
    write_bytes_atomic,
)
from recursive_horizons.nsc_ks_finite_history_uv_coefficients import (
    CUTOFF_BRIDGE_CONTROL,
    INVENTORY_RECORD,
)
from recursive_horizons.nsc_ks_cubic_uv_current import (
    SCHEMA,
    current_history_cubic_uv_report,
    validate_cubic_uv_report,
)

OUTPUT = ROOT / 'results/development/nsc-ks-cubic-uv-current-v1.json'
OWNERS = (
    'src/recursive_horizons/nsc_ks_cubic_uv_current.py',
    'scripts/derive_nsc_ks_cubic_uv_current.py',
    'tests/test_nsc_ks_cubic_uv_current.py',
    'docs/nsc-ks-cubic-uv-current.md',
)
INPUTS = (
    HISTORY_RECORD,
    CAUCHY_RECORD,
    INVENTORY_RECORD,
    'results/development/artifacts/nsc-ks-source-inventory.npz',
    CUTOFF_BRIDGE_CONTROL,
    'src/recursive_horizons/nsc_ks_current_uv_transport.py',
    'src/recursive_horizons/nsc_ks_uv_minor_a3.py',
    'src/recursive_horizons/nsc_ks_uv_upstream_invariance.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'results/development/nsc-ks-uv-upstream-invariance-v1.json',
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
    if isinstance(value, np.ndarray):
        return [plain(item) for item in value.tolist()]
    return value


def calculate():
    report = plain(current_history_cubic_uv_report(ROOT))
    return {
        **report,
        'status': (
            'OPEN: generated massless e^{-3} N,beta difference on I is '
            'constructed from actual Pauli vertices and consistent A3/n4 '
            'transport; diagnostic quadrature is unvalidated without a '
            'remainder proof; C4, C_M and the local gate remain missing'),
        'source_hashes': implementation_hashes(ROOT, owners=OWNERS),
        'input_hashes': {path: digest(path) for path in INPUTS},
    }


def check(value):
    if value.get('schema') != SCHEMA:
        raise ValueError('unexpected cubic UV schema')
    for path, expected in {**value['source_hashes'], **value['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('cubic UV dependency changed: ' + path)
    if (value['v1_pilot_rewritten'] or value['v2_record_rewritten']
            or value['v3_record_rewritten'] or value['invariance_record_rewritten']):
        raise ValueError('historical UV remainder/invariance bytes are immutable')
    if value['numerical_C4'] is not None or value['numerical_C_M'] is not None:
        raise ValueError('cubic successor may not invent a validated C4 or C_M')
    if value['validated_cubic_coefficient']:
        raise ValueError('diagnostic quadrature is not a remainder-validated C4')
    if value['vacuum_integrated_tail_N_beta'] is not None:
        raise ValueError('leading cancellation is not a vacuum-tail bound')
    validate_cubic_uv_report(value, ROOT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    computed = calculate()
    encoded = (
        json.dumps(computed, indent=2, sort_keys=True, allow_nan=False) + '\n'
    ).encode()
    if args.check and not OUTPUT.is_file():
        raise FileNotFoundError(
            'cubic UV record is missing; --check cannot create evidence')
    if args.record:
        if OUTPUT.is_file():
            if OUTPUT.read_bytes() != encoded:
                raise ValueError('existing cubic UV record differs')
        else:
            write_bytes_atomic(OUTPUT, encoded, exclusive=True)
        value = computed
    else:
        value = json.loads(OUTPUT.read_text())
        if value != computed:
            raise ValueError('cubic UV replay differs')
    check(value)
    diagnostic = value['diagnostic_quadrature']
    print(json.dumps({
        'status': value['status'],
        'schema': value['schema'],
        'profile_identity': value['profile_identity'][:16],
        'I_max_abs_N': diagnostic['I_max_abs_N'],
        'I_max_abs_beta': diagnostic['I_max_abs_beta'],
        'I_integral_N': diagnostic['I_integral_N'],
        'I_integral_beta': diagnostic['I_integral_beta'],
        'refinement_discrepancy_N': diagnostic['refinement_discrepancy_N'],
        'refinement_discrepancy_beta': diagnostic['refinement_discrepancy_beta'],
        'validated_cubic_coefficient': value['validated_cubic_coefficient'],
        'numerical_C4': value['numerical_C4'],
        'numerical_C_M': value['numerical_C_M'],
        'physical_local_gate': value['physical_local_gate'],
    }, indent=2))
