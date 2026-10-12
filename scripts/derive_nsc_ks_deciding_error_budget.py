#!/usr/bin/env python3
"""Enclose iterate4 geometry remainder; keep deciding leftovers None."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import BarycentricInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_between_node_geometry as G
import derive_nsc_ks_coupled_newton_n16_damped_iterate4 as I4
import derive_nsc_ks_coupled_newton_n16_jacobian_refresh as R
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget

OUTPUT = ROOT/'results/development/nsc-ks-deciding-error-budget.json'
OWNERS = (
    'scripts/derive_nsc_ks_deciding_error_budget.py',
    'docs/nsc-ks-deciding-error-budget.md',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def compute():
    trial = json.loads(I4.OUTPUT.read_text())
    refresh = json.loads(R.OUTPUT.read_text())
    for path, expected in {**trial['source_hashes'], **refresh['source_hashes'],
                           **refresh['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('deciding budget input changed: '+path)
    if trial['profile_identity'] != I4.EXPECTED_IDENTITY:
        raise ValueError('deciding budget must use the fifth clipped history')
    ctx = I4.context()
    family = ctx['family']
    if profile_identity(family, include_normal_window=True) != trial['profile_identity']:
        raise ValueError('deciding budget history identity changed')
    z = np.asarray(trial['z'], float)
    if len(z) != 47 or not np.array_equal(z, ctx['target']):
        raise ValueError('owned 47-node union required')
    residual = np.asarray(trial['action_gradient'], float)
    geo_nodes = G.geometry_at(family, z, ctx['coeff'])
    dense = family.collocation_nodes(G.DENSE_COUNT)
    if dense[0] < family.interval[0] - 1e-15 or dense[-1] > family.interval[1] + 1e-15:
        raise ValueError('dense nodes must stay inside the declared interval I')
    geo_dense = G.geometry_at(family, dense, ctx['coeff'])
    interpolants = [BarycentricInterpolator(z, geo_nodes[:, k]) for k in range(2)]
    geo_interp = np.column_stack([fn(dense) for fn in interpolants])
    geometry_remainder_raw = np.max(abs(geo_dense - geo_interp), axis=0)
    if np.any(~np.isfinite(geometry_remainder_raw)) or np.any(geometry_remainder_raw < 0):
        raise ValueError('finite nonnegative geometry interpolation remainder required')
    geometry_remainder = np.full(2, 1e-13) if np.all(geometry_remainder_raw < 1e-13) else np.nextafter(
        geometry_remainder_raw + 1e-14, np.inf)
    leftover = np.asarray(refresh['unclipped_linear_predicted_all_node_residual_maxima'], float)
    budget = local_error_budget(between_node_remainder=None)
    return {
        'schema': 'NSC-KS-DECIDING-ERROR-BUDGET-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: iterate4 geometry remainder enclosed; '
            'UV, full between-node, field and low/subgap remain None'),
        'profile_identity': trial['profile_identity'],
        'interval': list(family.interval),
        'owned_nodes': 47,
        'dense_nodes': G.DENSE_COUNT,
        'measured_constraint_maxima': np.asarray(trial['constraint_maxima'], float).tolist(),
        'unclipped_linear_leftover': leftover.tolist(),
        'geometry_between_node_remainder': geometry_remainder.tolist(),
        'full_between_node_remainder': None,
        'changed_history_UV_tail': None,
        'field_error_bound': None,
        'baseline_low_subgap': None,
        'residual_chebyshev_tail_max_last4': refresh['residual_chebyshev_tail_max_last4'],
        'error_budget': {name: None if value is None else value.tolist()
                         for name, value in budget.items()},
        'scope': {
            'physical_local_gate': 'OPEN',
            'geometry_remainder_reused_pattern': True,
            'missing_bounds_left_none': True,
            'identity_is_not_a_source_error_bound': True,
            'missing_bounds_cannot_create_EXISTENCE_against_1e-4_leftover': True,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {
            str(I4.OUTPUT.relative_to(ROOT)): digest(I4.OUTPUT),
            str(R.OUTPUT.relative_to(ROOT)): digest(R.OUTPUT),
        },
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('deciding error-budget register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('deciding error-budget replay differs')
    print(json.dumps({
        'status': result['status'],
        'geometry_between_node_remainder': result['geometry_between_node_remainder'],
        'full_between_node_remainder': result['full_between_node_remainder'],
        'changed_history_UV_tail': result['changed_history_UV_tail'],
        'field_error_bound': result['field_error_bound'],
        'baseline_low_subgap': result['baseline_low_subgap'],
    }, indent=2))
