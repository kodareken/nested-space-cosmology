#!/usr/bin/env python3
"""Record/replay the surface-scoped common-upstream UV identity."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_uv_upstream_invariance import upstream_normalization_identities

OUTPUT = 'results/development/nsc-ks-uv-upstream-invariance-v1.json'
OWNERS = (
    'src/recursive_horizons/nsc_ks_uv_upstream_invariance.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'scripts/derive_nsc_ks_uv_upstream_invariance.py',
    'tests/test_nsc_ks_uv_upstream_invariance.py',
    'docs/nsc-ks-uv-upstream-invariance.md',
)
INPUTS = (
    'results/development/nsc-ks-gate-history-lm-broyden.json',
    'results/development/nsc-ks-uv-minor-a3-v1.json',
)


def calculate():
    return {
        'schema': 'NSC-KS-UV-UPSTREAM-INVARIANCE-v1',
        'status': 'EXACT FORMAL IDENTITY on Sigma; numerical C4 and physical gate OPEN',
        **upstream_normalization_identities(),
        'class': 'delta r=chi(s_normal)*(s_normal*w+s_normal^3*U/6); fixed a; Sigma rho=1',
        'assumptions': {
            'same_source_and_upstream_neighborhood': True,
            'z_homogeneous_upstream_scalar_constants': True,
            'fixed_major_A1_upstream': True,
            'matched_incoming_a_and_r': True,
            'all_higher_coefficients_transformed_consistently': True,
            'occupied_vacuum_power_series_only': True,
        },
        'previous_minor_A3_record_rewritten': False,
        'source_hashes': implementation_hashes(ROOT, owners=OWNERS),
        'input_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.record:
        publish_exclusive_file(ROOT, OUTPUT,
            (json.dumps(result, sort_keys=True, indent=2)+'\n').encode())
    elif json.loads((ROOT/OUTPUT).read_text()) != result:
        raise ValueError('upstream invariance record or dependencies changed')
    print(json.dumps({'status': result['status'], 'exact_identities': len(result['identities']),
        'n2_up_assumed_zero': result['invariance_uses_n2_up_equals_zero'],
        'numerical_C4': None, 'numerical_C_M': None, 'physical_local_gate': 'OPEN'}, indent=2))
