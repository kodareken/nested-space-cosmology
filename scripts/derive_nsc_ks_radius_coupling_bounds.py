#!/usr/bin/env python3
"""Record two continuous field-method bounds for the current saved history."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_current_history_bounds import rational_record
from recursive_horizons.nsc_ks_radius_coupling_bounds import history_coupling_bounds
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

HISTORY = 'results/development/nsc-ks-gate-history-lm-broyden.json'
INVENTORY = 'results/development/nsc-ks-source-inventory.json'
OUTPUT = 'results/development/nsc-ks-radius-coupling-bounds-v1.json'
SOURCES = (
    'src/recursive_horizons/nsc_ks_radius_coupling_bounds.py',
    'src/recursive_horizons/nsc_ks_current_history_bounds.py',
    'src/recursive_horizons/nsc_ks_residual_error.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
    'scripts/derive_nsc_ks_radius_coupling_bounds.py',
    'tests/test_nsc_ks_radius_coupling_bounds.py',
    'docs/nsc-ks-radius-coupling-bounds.md',
)


def compute():
    history = json.loads((ROOT / HISTORY).read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    identity = profile_identity(family, include_normal_window=True)
    if identity != history['profile_identity']:
        raise ValueError('history identity mismatch')
    # Whole supported slab: safely contains recorded preparation rho_up.
    # Unit angular bounds scale exactly with the absolute channel label.
    from fractions import Fraction
    bounds = history_coupling_bounds(family, Fraction(33, 32), 1)
    return {
        'schema': 'NSC-KS-RADIUS-COUPLING-BOUNDS-v1',
        'status': 'ENCLOSED: continuous M and M_z integrals per unit absolute angular label',
        'profile_identity': identity,
        'scope': 'upstream half slab; fixed pure-radius history; per |ell|',
        'bounds_per_absolute_angular': rational_record(bounds),
        'physical_source_error_included': False,
        'field_error_component_closed': False,
        'physical_local_gate': 'OPEN',
        'source_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
        'input_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in (HISTORY, INVENTORY)},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = compute()
    if args.record:
        publish_exclusive_file(ROOT, OUTPUT, (json.dumps(value, indent=2, sort_keys=True)+'\n').encode())
    elif value != json.loads((ROOT/OUTPUT).read_text()):
        raise ValueError('radius coupling record differs')
    print(json.dumps({k: value[k] for k in ('status', 'bounds_per_absolute_angular', 'physical_local_gate')}, indent=2))
