#!/usr/bin/env python3
"""Exact analytic geometry bounds for the captured nonzero history; no solve."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_local_prepared_response as P
from recursive_horizons.nsc_ks_radius_enclosure import axial_profile_bounds, reciprocal_radius_tail
from recursive_horizons.nsc_ks_residual_error import pure_radius_commutator_integrals

OUTPUT = ROOT / 'results/development/nsc-ks-radius-bound.json'
OWNED = ('scripts/derive_nsc_ks_radius_bound.py', 'src/recursive_horizons/nsc_ks_radius_enclosure.py',
         'docs/nsc-ks-radius-enclosure.md', 'src/recursive_horizons/nsc_ks_residual_error.py',
         'docs/nsc-ks-residual-error.md')
INPUTS = ('results/development/nsc-ks-trajectory-control.json',
          'scripts/derive_nsc_local_prepared_response.py', 'src/recursive_horizons/nsc_local_incoming_family.py',
          'src/recursive_horizons/nsc_compatible_history_geometry.py', 'src/recursive_horizons/nsc_lorentzian.py')


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def pack(q):
    q = Q(q)
    number = float(q)
    if Q(number) < q:
        number = math.nextafter(number, math.inf)
    return {'exact_rational': str(q), 'binary64_upper': number}


def record():
    control = json.loads((ROOT / INPUTS[0]).read_text())
    payload = control['payload']
    if digest(payload['path']) != payload['sha256']:
        raise ValueError('captured source/history artifact changed')
    with np.load(ROOT / payload['path'], allow_pickle=False) as arrays:
        meta = json.loads(arrays['metadata_json'].tobytes())
    upper_rho = Q(33, 32)
    if not 1 < Q(meta['rho_up']) <= upper_rho:
        raise ValueError('slab exceeds the explicit geometry-bound domain')
    # beta^2 is decreasing for rho>0. atan(33/32)=pi/4+atan(1/65).
    # pi>157/50 and atan(1/65)<1/65 bound a^2=beta^2-1 below.
    a2_lower = 3*((1+upper_rho**2)*(Q(157, 200)-Q(1, 65))-upper_rho)-1
    a_lower, rref_lower = Q(4, 5), Q(7, 5)
    if a2_lower <= a_lower**2 or rref_lower**2 >= 2:
        raise ArithmeticError('elementary background lower bounds failed')
    family = P.family(P.ALPHA)
    W, U = [Q(0)]*3, [Q(0)]*3
    for amplitude, direction in zip(family.amplitudes, family.directions):
        for total, profile in ((W, direction.w), (U, direction.U)):
            bound = axial_profile_bounds(profile)
            for j in range(3):
                total[j] += abs(Q(amplitude))*bound[j]
    sigma = max(Q(d.outer_radius) for d in family.directions)
    ell = abs(Q(meta['parameters']['angular']))
    tail = reciprocal_radius_tail(W, U, normal_support=sigma,
        reference_radius_lower=rref_lower, axial_lower=a_lower, absolute_angular=ell, order=4)
    K1, K2 = pure_radius_commutator_integrals(W, U, sigma, tail['radius_lower'], ell)
    return {'schema': 'NSC-KS-RADIUS-BOUND-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: exact geometry remainder bounds for the control; physical local gate OPEN',
        'control': {'history_amplitude': P.ALPHA, 'local_interval': 'S(1)+[.12,.18]',
                    'prepared_source_digest': meta['preparation'], 'series_order': 4},
        'bounds': {'rho_upper': pack(upper_rho), 'axial_squared_lower': pack(a2_lower),
                   'axial_lower': pack(a_lower), 'reference_radius_lower': pack(rref_lower),
                   'w_sup_orders_0_1_2': [pack(v) for v in W], 'U_sup_orders_0_1_2': [pack(v) for v in U],
                   'radius_lower': pack(tail['radius_lower']), 'radius_ratio_upper': pack(tail['radius_ratio_upper']),
                   'potential_remainder_sup_orders_0_1_2': [pack(v) for v in tail['potential_derivative_tail_bounds']],
                   'B_z_integral_upper': pack(K1), 'B_zz_integral_upper': pack(K2)},
        'scope': {'finite_history_geometry_bound': True, 'profile_evaluation_roundoff_included': False,
                  'continuous_field_residual_bound': None, 'source_error_bound': None,
                  'physical_local_gate': 'OPEN', 'new_field_runs': 0, 'metric_timestep': False},
        'source_hashes': {p: digest(p) for p in OWNED}, 'input_hashes': {p: digest(p) for p in INPUTS},
        'reused_artifact': payload, 'reproducer': 'python3 scripts/derive_nsc_ks_radius_bound.py --check'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--record', action='store_true')
    group.add_argument('--check', action='store_true')
    result = record()
    if parser.parse_args().record:
        if OUTPUT.exists():
            raise FileExistsError('radius bound exists; use --check')
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('radius bound or its inputs no longer replay')
    print(json.dumps({'status': result['status'], 'potential_tail_orders_0_1_2':
        [v['binary64_upper'] for v in result['bounds']['potential_remainder_sup_orders_0_1_2']],
        'B_z_integral_upper': result['bounds']['B_z_integral_upper']['binary64_upper'],
        'B_zz_integral_upper': result['bounds']['B_zz_integral_upper']['binary64_upper']}, indent=2))
