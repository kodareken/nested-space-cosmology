#!/usr/bin/env python3
"""Record/replay the minor-A3 jet bound on the live current history.

Historical UV remainder v1/v2/v3 JSON stay immutable. This successor records
second geometry/phase jets, minor A2 and the jet piece of minor A3. Complete
minor A3, C4, C_M and the vacuum tail remain OPEN. Runtime is outside the
scientific digest.
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
from recursive_horizons.nsc_ks_evaluation_binding import write_bytes_atomic, implementation_hashes
from recursive_horizons.nsc_ks_finite_history_uv_coefficients import (
    CUTOFF_BRIDGE_CONTROL,
    INVENTORY_RECORD,
)
from recursive_horizons.nsc_ks_uv_minor_a3 import (
    SCHEMA,
    current_history_minor_a3_report,
    validate_minor_a3_report,
)

OUTPUT = ROOT / 'results/development/nsc-ks-uv-minor-a3-v1.json'
OWNERS = (
    'src/recursive_horizons/nsc_ks_uv_minor_a3.py',
    'scripts/derive_nsc_ks_uv_minor_a3.py',
    'tests/test_nsc_ks_uv_minor_a3.py',
    'docs/nsc-ks-uv-minor-a3.md',
)
INPUTS = (
    HISTORY_RECORD,
    CAUCHY_RECORD,
    INVENTORY_RECORD,
    'results/development/artifacts/nsc-ks-source-inventory.npz',
    CUTOFF_BRIDGE_CONTROL,
    'src/recursive_horizons/nsc_ks_current_uv_transport.py',
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
    report = plain(current_history_minor_a3_report(ROOT))
    return {
        **report,
        'status': (
            'OPEN: minor A3 jet piece and second geometry/phase jets enclosed '
            'on 0b0e4ced; complete history-minus-reference minor A3 retains '
            'unowned upstream major A2 times delta(1/r); C4, C_M and the UV '
            'tail remain missing'),
        'source_hashes': implementation_hashes(ROOT, owners=OWNERS),
        'input_hashes': {path: digest(path) for path in INPUTS},
    }


def check(value):
    if value.get('schema') != SCHEMA:
        raise ValueError('unexpected minor A3 schema')
    for path, expected in {**value['source_hashes'], **value['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('minor A3 dependency changed: ' + path)
    if value['v1_pilot_rewritten'] or value['v2_record_rewritten'] or value[
            'v3_record_rewritten']:
        raise ValueError('historical UV remainder bytes are immutable')
    if value['numerical_C4'] is not None or value['numerical_C_M'] is not None:
        raise ValueError('minor A3 successor may not invent C4 or C_M')
    if value['vacuum_integrated_tail_N_beta'] is not None:
        raise ValueError('leading cancellation is not a vacuum-tail bound')
    if not value['leading_e_minus2_cancellation_is_not_this_value']:
        raise ValueError('E^{-2} cancellation is not this bound')
    validate_minor_a3_report(value, ROOT)


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
            'minor A3 record is missing; --check cannot create evidence')
    if args.record:
        if OUTPUT.is_file():
            if OUTPUT.read_bytes() != encoded:
                raise ValueError('existing minor A3 record differs')
        else:
            write_bytes_atomic(OUTPUT, encoded, exclusive=True)
        value = computed
    else:
        value = json.loads(OUTPUT.read_text())
        if value != computed:
            raise ValueError('minor A3 replay differs')
    check(value)
    print(json.dumps({
        'status': value['status'],
        'schema': value['schema'],
        'profile_identity': value['profile_identity'][:16],
        'proved_jet': value['majorants'][
            'history_minus_reference_minor_A3_jet_abs_upper']['binary64_upper'],
        'missing_kernel': value['majorants'][
            'upstream_major_A2_difference_kernel_abs_upper']['binary64_upper'],
        'numerical_C4': value['numerical_C4'],
        'numerical_C_M': value['numerical_C_M'],
        'physical_local_gate': value['physical_local_gate'],
    }, indent=2))
