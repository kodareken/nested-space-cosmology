#!/usr/bin/env python3
"""Join certified high intervals and infinite tails without changing sources."""
import argparse
import json
from pathlib import Path
import sys

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from derive_nsc_incoming_source_update_v3 import authenticate, sha
from recursive_horizons.nsc_incoming_source_quadrature_bound import unpack, pack
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _hi, _up_float

OUTPUT = ROOT/'results/development/nsc-incoming-source-update-v4.json'
INPUTS = ('results/development/nsc-incoming-source-update-v3.json',
          'results/development/nsc-incoming-tail-quadrature-batch.json',
          'scripts/derive_nsc_incoming_source_update_v3.py',
          'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py')
SOURCES = ('scripts/derive_nsc_incoming_source_update_v4.py',
           'tests/test_nsc_incoming_source_update_v4.py',
           'docs/nsc-incoming-source-update-v4.md')


def join_coverage(previous, tail):
    rows = tail['groups']
    if len(rows) != 32 or {r['group'] for r in rows} != set(range(1, 33)):
        raise ValueError('exactly32 distinct tail groups required')
    if set(previous['physical_coverage']) != {str(g) for g in range(1, 33)}:
        raise ValueError('exactly32 finite high groups required')
    coverage = {}
    for row in rows:
        group = row['group']; pieces = previous['physical_coverage'][str(group)]
        end = max(p['interval'][1] for p in pieces)
        start = min(p['interval'][0] for p in pieces)
        if row['lower_energy'] != end or not row['archived_source_unchanged']:
            raise ValueError('tail must meet its unchanged finite high region exactly')
        coverage[group] = {'finite_high': [start, end], 'tail': [end, 'infinity'],
                           'joined': [start, 'infinity'], 'actual_signs': row['actual_signs']}
    if tail['aggregate']['archived_total_order13_source'] != previous['matter']['parts']['paired_approximate_tail']:
        raise ValueError('tail certificate must bound the exact source values already assembled')
    return coverage


def make_record():
    previous, tail = [json.loads((ROOT/p).read_text()) for p in INPUTS[:2]]
    for record in (previous, tail):
        authenticate(record)
    coverage = join_coverage(previous, tail)
    if not tail['batch_complete'] or tail['failures']:
        raise ValueError('completed authenticated tail batch required')
    with _precision(80):
        high = [mp.iv.mpf(v) for v in previous['partial_error_budget']['finite_high_plus_source_sum_rounding_action_error_upper']]
        tail_bound = [unpack(v) for v in tail['aggregate']['combined_action_upper_binary']]
        if any(_hi(v) <= 0 for v in tail_bound):
            raise ValueError('retain the positive physical thermal remainder, including shift')
        total = [a+b for a, b in zip(high, tail_bound)]
        budget = {'constraint_order': ['N', 'beta'],
                  'high_finite_plus_source_sum_rounding_action_error_upper': [_up_float(v) for v in high],
                  'tail_physical_numerical_thermal_action_error_upper': [_up_float(v) for v in tail_bound],
                  'tail_action_error_upper_binary': [pack(v) for v in tail_bound],
                  'covered_high_to_infinity_action_error_upper': [_up_float(v) for v in total],
                  'covered_high_to_infinity_action_error_upper_binary': [pack(v) for v in total],
                  'stationarity_tolerance': 3e-11,
                  'covered_domain_status': 'PASS' if all(_hi(v) < mp.mpf('3e-11') for v in total) else 'OPEN',
                  'remaining_low_subgap_source_error_bound': None,
                  'full_source_error_bound': None, 'full_source_status': 'OPEN',
                  'finite_high_physical_pressure_error_certified': False,
                  'joint_local_reference_action_arithmetic_certified': False}
    matter = {**previous['matter'], 'unresolved': [item for item in previous['matter']['unresolved']
              if item != 'tail numerical integration, distinct from physical mode error']}
    return {'schema': 'NSC-INCOMING-SOURCE-UPDATE-v4', 'accountable_author': 'Douglas Ek',
            'status': budget['covered_domain_status']+': high constraint-source budget through infinity certified; lower source and constraints OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
            'coverage': coverage, 'matter': matter, 'baseline': previous['baseline'],
            'mixed_normal_probe': previous['mixed_normal_probe'], 'partial_error_budget': budget,
            'source_value_change': [0., 0., 0., 0.],
            'scope': {'original_tail_source_preserved': True, 'tail_physical_error_counted_once': True,
                      'finite_high_source_preserved': True, 'thermal_declared_zero': False,
                      'scientific_producers_run': False, 'physical_IV_selected': False,
                      'metric_evolution': False, 'Gamma_rest_assigned': False,
                      'extended_stationarity': 'OPEN', 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_source_update_v4.py --check'}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write and OUTPUT.exists():
        raise FileExistsError('existing source-budget receipt is not overwritten')
    data = json.loads(json.dumps(make_record()))
    if args.write:
        OUTPUT.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    elif json.loads(OUTPUT.read_text()) != data:
        raise ValueError('source-budget replay differs')
    print(json.dumps({'status': data['status'], 'budget': data['partial_error_budget']}, indent=2))


if __name__ == '__main__':
    main()
