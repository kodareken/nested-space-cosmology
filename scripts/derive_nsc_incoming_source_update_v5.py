#!/usr/bin/env python3
"""Compose the certified group32 lower window; no field producer is rerun."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from derive_nsc_incoming_source_update_v3 import authenticate, sha
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _hi, _up_float

OUTPUT = ROOT/'results/development/nsc-incoming-source-update-v5.json'
INPUTS = ('results/development/nsc-incoming-source-update-v4.json',
          'results/development/nsc-incoming-validated-energy-window.json',
          'scripts/derive_nsc_incoming_source_update_v3.py',
          'src/recursive_horizons/nsc_incoming_joint_constraints.py',
          'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py')
SOURCES = ('scripts/derive_nsc_incoming_source_update_v5.py',
           'tests/test_nsc_incoming_source_update_v5.py',
           'docs/nsc-incoming-source-update-v5.md')


def validate_window(previous, window):
    if window['group'] != 32 or window['energy_interval'] != [24, 32]:
        raise ValueError('owned group32 lower cell required')
    if previous['coverage']['32']['joined'] != [32, 'infinity']:
        raise ValueError('new lower cell must meet, not overlap, prior coverage')
    if not window['complete_window_budget_pass']:
        raise ValueError('complete field, quadrature and thermal certificate required')
    if window['kernel_order'] != ['rho', 'p_parallel', 'T01', 'p_perp']:
        raise ValueError('same source component convention required')
    if window['source']['endpoint_column_normalized'] or window['source']['exact_vacuum_current'] != 0:
        raise ValueError('unnormalized polynomial and explicit vacuum-current identity required')
    if window['positive_thermal_stress_error_upper'][2] <= 0:
        raise ValueError('positive physical thermal current remainder required')
    coverage = deepcopy(previous['coverage'])
    coverage['32']['validated_lower_window'] = [24, 32]
    coverage['32']['joined'] = [24, 'infinity']
    return coverage


def make_record():
    previous, window = [json.loads((ROOT/p).read_text()) for p in INPUTS[:2]]
    for record in (previous, window):
        authenticate(record)
    coverage = validate_window(previous, window)
    old = np.asarray(window['original_source'])
    new = np.asarray(window['validated_vacuum_source'])
    delta = new-old
    prior_stress = np.asarray(previous['matter']['stress_approximant'])
    stress = prior_stress+delta
    action = source_action_gradient(incoming_cauchy_jets(), delta)
    residuals = {'explicit_delta': float(np.max(abs(delta-window['explicit_source_delta'])))}
    matter = deepcopy(previous['matter'])
    if 'group32_low_24_32_replacement' in matter['parts']:
        raise ValueError('lower-window replacement was already applied')
    matter['parts']['group32_low_24_32_replacement'] = delta.tolist()
    matter['previous_stress_approximant'] = prior_stress.tolist()
    matter['stress_approximant'] = stress.tolist()
    cases = {}
    for name in ('baseline', 'mixed_normal_probe'):
        before = previous[name]
        values = np.asarray(before['action_gradient_approximant'])+action[:2]
        cases[name] = {**before, 'action_gradient_approximant': values.tolist(),
                       'force_approximant': (-values).tolist(),
                       'matter_baseline_action_gradient': (np.asarray(before['matter_baseline_action_gradient'])+action[:2]).tolist(),
                       'maximum_absolute_approximant': float(max(abs(values)))}
        residuals[name+'_increment_identity'] = float(max(abs(values-before['action_gradient_approximant']-action[:2])))
    if max(residuals.values()) > 3e-13:
        raise ArithmeticError('lower-window source composition failed')
    with _precision(80):
        # Incremental rounding only: v4 already carries its own source-sum
        # rounding. Fixed old LOW values are removed once, not re-certified.
        rounding = [abs(mp.iv.mpf(prior_stress[i])+mp.iv.mpf(new[i])-mp.iv.mpf(old[i])-mp.iv.mpf(stress[i])) for i in range(4)]
        a = mp.iv.sqrt(3*mp.iv.pi/2-4)
        extra = [8*mp.iv.pi*a*rounding[0], 8*mp.iv.pi*a*a*rounding[2]]
        old_budget = [mp.iv.mpf(v) for v in previous['partial_error_budget']['covered_high_to_infinity_action_error_upper']]
        window_budget = [mp.iv.mpf(v) for v in window['constraint_action_error_upper']]
        total = [old_budget[i]+window_budget[i]+extra[i] for i in range(2)]
        budget = {'constraint_order': ['N', 'beta'],
                  'previous_covered_action_error_upper': [_up_float(v) for v in old_budget],
                  'new_window_complete_action_error_upper': [_up_float(v) for v in window_budget],
                  'incremental_source_sum_rounding_action_error_upper': [_up_float(v) for v in extra],
                  'covered_spectral_regions_action_error_upper': [_up_float(v) for v in total],
                  'stationarity_tolerance': 3e-11,
                  'covered_domain_status': 'PASS' if all(_hi(v) < mp.mpf('3e-11') for v in total) else 'OPEN',
                  'full_source_error_bound': None, 'full_source_status': 'OPEN',
                  'remaining_low_subgap_source_error_bound': None,
                  'whole_action_arithmetic_certified': False,
                  'all_covered_pressure_errors_certified': False}
    return {'schema': 'NSC-INCOMING-SOURCE-UPDATE-v5', 'accountable_author': 'Douglas Ek',
            'status': budget['covered_domain_status']+': certified group32 lower cell composed; full numerical source and stationarity OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
            'coverage': coverage, 'matter': matter, **cases,
            'new_window': {'group': 32, 'interval': [24, 32], 'original_source': old.tolist(),
                           'validated_vacuum_source': new.tolist(), 'replacement_delta': delta.tolist(),
                           'positive_thermal_stress_error_upper': window['positive_thermal_stress_error_upper']},
            'action_gradient_increment': action.tolist(), 'partial_error_budget': budget,
            'residuals': residuals, 'composition_tolerance': 3e-13,
            'scope': {'new_cell_replaced_once': True, 'prior_high_and_tail_values_preserved': True,
                      'LLL_and_group13_subgap_values_preserved': True, 'thermal_declared_zero': False,
                      'old_endpoint_columns_normalized': False, 'scientific_producers_run': False,
                      'physical_initial_data_selected': False, 'metric_evolution': False,
                      'Gamma_rest_assigned': False, 'extended_stationarity': 'OPEN', 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_source_update_v5.py --check'}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write and OUTPUT.exists():
        raise FileExistsError('existing source receipt is not overwritten')
    data = json.loads(json.dumps(make_record()))
    if args.write:
        OUTPUT.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    elif json.loads(OUTPUT.read_text()) != data:
        raise ValueError('lower-window source composition replay differs')
    print(json.dumps({k: data[k] for k in ('status', 'action_gradient_increment', 'partial_error_budget', 'residuals')}, indent=2))


if __name__ == '__main__':
    main()
