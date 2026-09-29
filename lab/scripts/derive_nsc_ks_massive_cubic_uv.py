#!/usr/bin/env python3
"""Record one directed group14 massive cubic coefficient and replay it."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_current_uv_transport import (
    ARCHIVE_RHO_UP, CAUCHY_RECORD, CURRENT_HISTORY_IDENTITY, HISTORY_RECORD,
)
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_massive_cubic_uv import (
    CUTOFF_BRIDGE_RECORD, INVENTORY_RECORD, SCHEMA, MassiveCubicGeometry,
    enclose_massive_characteristic, massive_cubic_identities,
    original_group14_labels, require_original_history,
)

OUTPUT = 'results/development/nsc-ks-massive-cubic-uv-v1.json'
OWNERS = (
    'scripts/derive_nsc_ks_massive_cubic_uv.py',
    'src/recursive_horizons/nsc_ks_massive_cubic_uv.py',
    'tests/test_nsc_ks_massive_cubic_uv.py',
    'docs/nsc-ks-massive-cubic-uv.md',
)
INPUTS = (HISTORY_RECORD, CAUCHY_RECORD, CUTOFF_BRIDGE_RECORD, INVENTORY_RECORD)
PILOT_CELLS = 1024
PILOT_ORDER = 8


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def interval(value):
    return {'lower': exact_upper(value.lower()), 'upper': exact_upper(value.upper())}


def calculate():
    proven = massive_cubic_identities()
    if any(value != '0' for value in proven['residuals'].values()):
        raise ValueError('massive cubic residual is not an exact zero')
    for name in ('higher_uv_remainder_C_M', 'uniform_C4_on_I', 'complete_UV_tail'):
        if proven[name] is not None:
            raise ValueError(name + ' was presented as constructed')
    if proven['physical_local_gate'] != 'OPEN':
        raise ValueError('massive cubic record must leave the local gate open')
    if proven['massless_h_equation_imposed']:
        raise ValueError('massless h transport was imposed on the massive coefficient')
    labels = original_group14_labels(ROOT)
    if (labels['group'] != 14 or labels['multiplicity_per_signed_family'] != 12
            or labels['cutoff'] != 160 or not labels['positive_mass']):
        raise ValueError('group14 multiplicity, cutoff or mass label changed')
    family, identity = require_original_history(ROOT)
    if identity != CURRENT_HISTORY_IDENTITY:
        raise ValueError('original current history required')
    rows = []
    started = time.process_time()
    with ctx.workprec(160):
        mass, angular = labels['mass_ball'], labels['angular_ball']
        model = MassiveCubicGeometry(family, angular, mass, bits=160)
        incoming = arb(family.center)
        weight = 2 * arb(labels['multiplicity_per_signed_family']) / (2 * arb.pi())
        paired_n, paired_beta = arb(0), arb(0)
        for sign in (1, -1):
            bound = enclose_massive_characteristic(
                model, incoming, sign, ARCHIVE_RHO_UP,
                cells=PILOT_CELLS, order=PILOT_ORDER)
            if (bound['C_M'] is not None or bound['complete_UV_tail'] is not None
                    or bound['entire_incoming_interval'] or bound['physical_local_gate'] != 'OPEN'):
                raise ValueError('point pilot was labeled as tail or interval closure')
            rows.append({
                'cells': PILOT_CELLS, 'order': PILOT_ORDER, 'sign': sign,
                'accepted_cells': bound['accepted_cells'],
                'values': {key: interval(bound[key]) for key in
                           ('q_z', 'q_zz', 'J3', 'N_bracket3', 'surface_mass_correction')},
            })
            paired_n -= weight * bound['N_bracket3']
            paired_beta += weight * bound['J3']
        cutoff = arb(labels['cutoff'])
        tail_divisor = 2 * cutoff**2
        result = {
            'schema': SCHEMA,
            'status': (
                'ENCLOSED: formal massive cubic coefficient at one incoming '
                'point of original group 14; not a UV tail; physical gate OPEN'),
            'profile_identity': identity,
            'channel': {
                'group': 14,
                'positive_mass': True,
                'mass_hex': float(labels['mass']).hex(),
                'absolute_angular_hex': float(labels['absolute_angular']).hex(),
                'multiplicity_per_signed_family': 12,
                'degeneracy': int(labels['degeneracy']),
                'quadrature_cutoff': 160.0,
                'quadrature_measure_by_sign': {'1': 160.0, '-1': 160.0},
                'archived_thermal_splits_not_a_tail': [160.0, 320.0],
            },
            'bits': 160,
            'rho_up_hex': float(ARCHIVE_RHO_UP).hex(),
            'incoming_point_hex': float(family.center).hex(),
            'mass_enclosure': interval(mass),
            'angular_enclosure': interval(angular),
            'cases': rows,
            'paired_action_cubic': {'N': interval(paired_n), 'beta': interval(paired_beta)},
            'formal_leading_term_at_cutoff_160': {
                'N': interval(paired_n / tail_divisor),
                'beta': interval(paired_beta / tail_divisor),
                'is_complete_uv_tail': False,
            },
            'leading_term_is_complete_uv_tail': False,
            'angular_evenness_used': True,
            'source_signs': [1, -1],
            'measure': (
                'mu=12 twice for the angular pair, divided by 2 pi; '
                'source signs separate; N minus and beta plus'),
            'initial_difference_values': [0, 0, 0],
            'initial_physical_normalization_claimed_zero': False,
            'massless_h_equation_imposed': False,
            'identities': {
                'residuals': dict(proven['residuals']),
                'Ts_qz': proven['Ts_qz'],
                'Ts_qzz': proven['Ts_qzz'],
                'massless_Ts_J3': proven['massless_Ts_J3'],
                'Ts_J3_massive': proven['Ts_J3_massive'],
                'Sigma_delta_N_bracket': proven['Sigma_delta_N_bracket'],
                'L_h_parent': proven['L_h_parent'],
                'factorization': dict(proven['factorization']),
                'action_N': proven['action_N'],
                'action_beta': proven['action_beta'],
            },
            'coverage': {
                'group': 14, 'incoming_points': 1, 'source_signs': [1, -1],
                'cells': PILOT_CELLS, 'order': PILOT_ORDER,
                'entire_incoming_interval': False, 'all_source_families': False,
                'factorization_integrated_for_all_groups': False,
            },
            'uniform_C4_on_I': None,
            'higher_uv_remainder_C_M': None,
            'physical_source_error': None,
            'physical_UV_budget_component': None,
            'complete_UV_tail': None,
            'physical_local_gate': 'OPEN',
            'source_inventory_unchanged': True,
            'source_hashes': implementation_hashes(ROOT, owners=OWNERS),
            'input_hashes': {path: digest(path) for path in INPUTS},
        }
    return result, time.process_time() - started


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result, cpu = calculate()
    encoded = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()
    if args.record:
        publish_exclusive_file(ROOT, OUTPUT, encoded)
    elif (ROOT / OUTPUT).read_bytes() != encoded:
        raise ValueError('massive cubic coefficient or dependencies changed')
    def endpoint(row):
        return {end: float(restored_upper(value)) for end, value in row.items()}
    print(json.dumps({
        'status': result['status'],
        'pilot_cpu_seconds': cpu,
        'accepted_cells': [row['accepted_cells'] for row in result['cases']],
        'paired_action_cubic': {key: endpoint(row)
                                for key, row in result['paired_action_cubic'].items()},
        'formal_leading_term_at_cutoff_160': {
            key: endpoint(row) for key, row in result['formal_leading_term_at_cutoff_160'].items()
            if key != 'is_complete_uv_tail'},
        'leading_term_is_complete_uv_tail': False,
        'higher_uv_remainder_C_M': result['higher_uv_remainder_C_M'],
        'uniform_C4_on_I': result['uniform_C4_on_I'],
        'complete_UV_tail': result['complete_UV_tail'],
        'physical_UV_budget_component': result['physical_UV_budget_component'],
        'physical_local_gate': result['physical_local_gate'],
    }, indent=2))
