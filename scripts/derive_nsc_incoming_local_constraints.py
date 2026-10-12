#!/usr/bin/env python3
"""Record only the new local Euler-force binding and its focused controls."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'tests'))
from test_nsc_incoming_local_constraints import (
    make_owners, inherited_light_check, compact_neck_check,
    weak_variation_check, radius_check,
)


OUTPUT = ROOT/'results/development/nsc-incoming-local-constraints.json'
SOURCES = (
    'src/recursive_horizons/nsc_incoming_local_constraints.py',
    'tests/test_nsc_incoming_local_constraints.py',
    'scripts/derive_nsc_incoming_local_constraints.py',
    'docs/nsc-incoming-local-constraints.md',
)
INPUTS = (
    'results/development/nsc-subgap-history-response.json',
    'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
    'src/recursive_horizons/nsc_light_restoration_action.py',
    'src/recursive_horizons/nsc_magnetic_light_reference.py',
    'src/recursive_horizons/nsc_spherical_local_history.py',
    'src/recursive_horizons/nsc_spatial_reference_symbol.py',
    'src/recursive_horizons/nsc_angular_stress.py',
    'src/recursive_horizons/nsc_horizon_source.py',
    'src/recursive_horizons/nsc_charged_ctp_neck.py',
    'scripts/derive_nsc_transmitting_boundary_binding.py',
)


def make_record():
    owners = make_owners()
    checks = {
        'incoming_light_tensor': inherited_light_check(owners),
        'committed_compact_neck_tensor': compact_neck_check(owners),
        'independent_nonhomogeneous_weak_variation': weak_variation_check(owners),
        'independent_numerical_radii': radius_check(owners),
    }
    baseline = checks['incoming_light_tensor'].pop('local_result')
    control = checks['independent_nonhomogeneous_weak_variation'].pop('local_result')
    maxima = {k: v['maximum_error'] for k, v in checks.items()}
    maxima.update({
        'Cauchy_resolution': max(baseline['maximum_numerical_indicator'], control['maximum_numerical_indicator']),
        'Euler_bulk_identity': max(baseline['Euler_bulk_identity_residual'], control['Euler_bulk_identity_residual']),
    })
    if max(maxima.values()) > 3e-11:
        raise ArithmeticError('new local-force verification exceeded 3e-11')
    sha = lambda path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    return {
        'schema': 'NSC-INCOMING-LOCAL-CONSTRAINTS-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'PASS: local lapse/shift Euler forces of the existing action',
        'source_hashes': {p: sha(p) for p in SOURCES},
        'input_hashes': {p: sha(p) for p in INPUTS},
        'locked_inputs': owners[1].ledger,
        'baseline': baseline, 'mixed_normal_jet_control': control,
        'independent_checks': checks, 'maxima': maxima,
        'numerical_verification_tolerance': 3e-11,
        'constraint_stationarity_tolerance': 3e-11,
        'full_constraint_residual_or_pass': None,
        'scope': {
            'equations_evaluated': ['local action Euler N', 'local action Euler beta'],
            'physical_state_reference_band_terms_included': False,
            'spatial_metric_equations_checked_only_as_homogeneous_normalization_controls': True,
            'physical_history_selected': False, 'normal_jets_selected': False,
            'constraint_solution_or_stationarity_claimed': False,
            'endpoint_variations_completed': False,
            'canonical_coordinates': 'same KS spatial foliation and incoming intrinsic frame',
            'differentiation': 'finite Taylor germ plus complex Cauchy coefficients; no physical patch evolution',
            'force_sign': 'local_force = -delta S_local/d(N,beta), per dT dz',
            'allocation': 'light geometry once; locked compact action once; no occupation or band insertion',
        },
        'reproducer': 'python3 scripts/derive_nsc_incoming_local_constraints.py --check',
        'comparison': {'fields': 'all', 'float_atol': 3e-13, 'float_rtol': 3e-13,
                       'exact': 'types, hashes, labels, structure and scope'},
    }


def json_value(value):
    # On some NumPy builds longdouble.tolist() remains a NumPy scalar.
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError('unsupported record value '+type(value).__name__)


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.prepare and OUTPUT.exists():
        raise FileExistsError('existing record is not overwritten: '+str(OUTPUT))
    data = json.loads(json.dumps(make_record(), default=json_value))
    if args.prepare:
        OUTPUT.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(json.loads(OUTPUT.read_text()), data)
    print(json.dumps({k: data[k] for k in ('status', 'maxima', 'full_constraint_residual_or_pass')}, indent=2))


if __name__ == '__main__':
    main()
