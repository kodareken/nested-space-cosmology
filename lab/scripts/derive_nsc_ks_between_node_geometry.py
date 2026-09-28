#!/usr/bin/env python3
"""Enclose the geometry interpolation remainder on I; keep the full remainder None.

The 47-node residual is not sup_I. This owner evaluates the local/reference
geometry term densely and compares it to the barycentric interpolant of the
same term on the owned 47-node union. That difference is the geometry
between-node remainder. The matter/source remainder stays None, so the full
gate between-node bound stays None.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.chebyshev import chebfit
from scipy.interpolate import BarycentricInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n16_damped_iterate2 as I2
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-between-node-geometry.json'
DENSE_COUNT = 129
OWNERS = (
    'scripts/derive_nsc_ks_between_node_geometry.py',
    'docs/nsc-ks-between-node-geometry.md',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def geometry_at(family, z, coeff):
    slots, tangents = compatible_history_slots(family.metric(), z, 32)
    return surface_geometry_response(slots, tangents, coeff)['action_gradient_change']


def chebyshev_tail(z, values, center, halfwidth):
    x = (np.asarray(z, float) - center) / halfwidth
    coeff = chebfit(x, np.asarray(values, float), deg=min(len(z) - 1, 46))
    return coeff, float(np.max(abs(coeff[-4:])))


def compute():
    trial = json.loads(I2.OUTPUT.read_text())
    if trial['completed_positive_families'] != trial['required_positive_families']:
        raise ValueError('between-node geometry requires a complete clipped residual')
    ctx = I2.context()
    family = ctx['family']
    if profile_identity(family) != trial['profile_identity']:
        raise ValueError('between-node geometry must use the latest complete clipped history')
    z = np.asarray(trial['z'], float)
    if len(z) != 47 or not np.array_equal(z, ctx['target']):
        raise ValueError('owned 47-node union required')
    residual = np.asarray(trial['action_gradient'], float)
    geo_nodes = geometry_at(family, z, ctx['coeff'])
    dense = family.collocation_nodes(DENSE_COUNT)
    if dense[0] < family.interval[0] - 1e-15 or dense[-1] > family.interval[1] + 1e-15:
        raise ValueError('dense nodes must stay inside the declared interval I')
    geo_dense = geometry_at(family, dense, ctx['coeff'])
    interpolants = [BarycentricInterpolator(z, geo_nodes[:, k]) for k in range(2)]
    geo_interp = np.column_stack([fn(dense) for fn in interpolants])
    geometry_remainder_raw = np.max(abs(geo_dense - geo_interp), axis=0)
    if np.any(~np.isfinite(geometry_remainder_raw)) or np.any(geometry_remainder_raw < 0):
        raise ValueError('finite nonnegative geometry interpolation remainder required')
    geometry_remainder = np.full(2, 1e-13) if np.all(geometry_remainder_raw < 1e-13) else np.nextafter(
        geometry_remainder_raw + 1e-14, np.inf)
    residual_interp = [BarycentricInterpolator(z, residual[:, k]) for k in range(2)]
    residual_dense_interp = np.column_stack([fn(dense) for fn in residual_interp])
    residual_interp_maxima = np.round(np.max(abs(residual_dense_interp), axis=0), 12)
    center = family.center
    halfwidth = family.axial_inner
    tails = [round(chebyshev_tail(z, residual[:, k], center, halfwidth)[1], 12) for k in range(2)]
    budget = local_error_budget(
        between_node_remainder=None,
    )
    return {
        'schema': 'NSC-KS-BETWEEN-NODE-GEOMETRY-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: geometry between-node remainder enclosed; '
            'full residual between-node bound remains None'),
        'profile_identity': trial['profile_identity'],
        'interval': list(family.interval),
        'owned_nodes': 47,
        'dense_nodes': DENSE_COUNT,
        'measured_nodal_residual_maxima': np.max(abs(residual), axis=0).tolist(),
        'interpolated_residual_maxima_on_dense_nodes': residual_interp_maxima.tolist(),
        'residual_chebyshev_tail_max_last4': tails,
        'geometry_between_node_remainder': geometry_remainder.tolist(),
        'full_between_node_remainder': None,
        'changed_history_tail': None,
        'field_error_bound': None,
        'baseline_low_subgap': None,
        'error_budget': {name: None if value is None else value.tolist()
                         for name, value in budget.items()},
        'scope': {
            'physical_local_gate': 'OPEN',
            'geometry_remainder_is_not_the_full_residual': True,
            'dense_nodes_do_not_replace_a_derivative_majorant_for_matter': True,
            'missing_bounds_left_none': True,
            'new_field_or_source_evolutions': 0,
            'metric_timestep': False,
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_claimed': False,
        },
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {str(I2.OUTPUT.relative_to(ROOT)): digest(I2.OUTPUT)},
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
            raise FileExistsError('between-node geometry register exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('between-node geometry replay differs')
    print(json.dumps({
        'status': result['status'],
        'geometry_between_node_remainder': result['geometry_between_node_remainder'],
        'full_between_node_remainder': result['full_between_node_remainder'],
        'measured_nodal_residual_maxima': result['measured_nodal_residual_maxima'],
    }, indent=2))
