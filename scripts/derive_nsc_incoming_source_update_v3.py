#!/usr/bin/env python3
"""Join completed finite-high source and error receipts; no producers run."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_joint_constraints import baseline_matter_source, source_action_gradient
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _up_float, _lo, _hi

OUTPUT = ROOT/'results/development/nsc-incoming-source-update-v3.json'
NAMES = ('spectral-source', 'subgap-source', 'source-tail', 'joint-constraints',
         'high-band-quadrature-batch', 'high-band-batch-bound', 'centered-order24',
         'retained-order24-bound', 'middle-order-correction', 'retained-order-correction',
         'low-high-source', 'window-reuse', 'group32-order24-bound', 'vacuum-tail-bound')
SOURCES = ('scripts/derive_nsc_incoming_source_update_v3.py',
           'tests/test_nsc_incoming_source_update_v3.py', 'docs/nsc-incoming-source-update-v3.md')
OWNERS = ('src/recursive_horizons/nsc_incoming_joint_constraints.py',
          'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py')


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def authenticate(record):
    """Check declared producer/input hashes and every nested payload receipt."""
    for key in ('source_hashes', 'input_hashes', 'numerical_owner_and_input_hashes', 'prepared_dependency_hashes'):
        for path, digest in record.get(key, {}).items():
            if sha(path) != digest:
                raise ValueError('composition dependency changed: '+path)
    def visit(value):
        if isinstance(value, dict):
            if 'path' in value and 'sha256' in value and sha(value['path']) != value['sha256']:
                raise ValueError('composition payload changed: '+value['path'])
            for item in value.values():
                visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)
    visit(record)


def expected_interval(group):
    return [16 if group in (13, 14) else 32, 320 if group in (10, 11, 12, 31, 32) else 160]


def validate_coverage(direct, physical):
    if set(direct) != {str(g) for g in range(1, 33)} or set(physical) != set(range(1, 33)):
        raise ValueError('exactly32 high source and physical groups required')
    for group in range(1, 33):
        row = direct[str(group)]
        if row['group'] != group or row['source_order'] != 24 or row['old_baseline_quadrature_carried']:
            raise ValueError('matching direct order24 source required')
        if row['interval'] != expected_interval(group):
            raise ValueError('source interval differs from owned high region')
        pieces = sorted(physical[group], key=lambda p: p['interval'])
        if not pieces or pieces[0]['interval'][0] != row['interval'][0] or pieces[-1]['interval'][1] != row['interval'][1]:
            raise ValueError('physical certificate union differs from direct source')
        for piece in pieces:
            if piece['order'] != 24 or piece['interval'][0] >= piece['interval'][1]:
                raise ValueError('matching physical order24 certificate required')
        if any(a['interval'][1] != b['interval'][0] for a, b in zip(pieces, pieces[1:])):
            raise ValueError('physical certificate overlap or gap')


def physical_pieces(records):
    r = records
    pieces = {}
    def add(group, interval, energy, thermal, owner):
        if energy['physical_Riccati_order'] != 24 or energy['vacuum_current_error'] != 0:
            raise ValueError('owned energy/current certificate required')
        pieces.setdefault(group, []).append({'interval': interval, 'order': 24,
            'vacuum_density_error_upper': energy['density_error_upper'],
            'thermal_stress_error_upper': thermal, 'owner': owner})
    for row in r['high-band-batch-bound']['groups']:
        group = row['group']
        if group in pieces or group in (12, 14, 22, 32):
            raise ValueError('duplicate or overlapping physical group')
        window = row['windows']['middle' if group == 13 else 'combined']
        add(group, window['interval'], window['vacuum_energy'], window['thermal_stress_upper'], 'high-band-batch-bound')
    add(14, [16, 160], r['centered-order24']['energy'],
        r['middle-order-correction']['thermal_scattering_stress_error_upper'], 'centered-order24 + middle-order-correction')
    add(12, r['retained-order24-bound']['energy_interval'], r['retained-order24-bound']['energy'],
        r['retained-order-correction']['thermal_scattering_bound'], 'retained-order24-bound + retained-order-correction')
    add(22, r['low-high-source']['interval'], r['low-high-source']['physical_mode_energy_bound'],
        r['low-high-source']['omitted_thermal_stress_upper'], 'low-high-source')
    for name in ('group12_low', 'group22_middle'):
        row = r['window-reuse']['windows'][name]
        add(row['group'], row['interval'], row['interval_binding']['energy'],
            row['thermal_scattering_stress_upper'], 'window-reuse/'+name)
    row = r['group32-order24-bound']['windows']['combined']
    add(32, row['energy_interval'], row['vacuum_energy'], row['thermal_stress_upper'], 'group32-order24-bound/combined')
    return pieces


def make_record():
    records = {n: json.loads((ROOT/f'results/development/nsc-incoming-{n}.json').read_text()) for n in NAMES}
    for record in records.values():
        authenticate(record)
    direct = records['high-band-quadrature-batch']
    physical = physical_pieces(records)
    validate_coverage(direct['groups'], physical)
    old_matter = baseline_matter_source(*(records[n] for n in NAMES[:3]))
    aggregate = direct['aggregate']
    old_high = np.asarray(aggregate['original_high_source'])
    new_high = np.asarray(aggregate['direct_high_vacuum_source'])
    delta = new_high-old_high
    stress = np.asarray(old_matter['stress_approximant'])+delta
    action = source_action_gradient(incoming_cauchy_jets(), delta)
    residuals = {'replacement_delta': float(np.max(abs(delta-aggregate['replacement_delta_from_original_high_source']))),
                 'group_original_sum': float(np.max(abs(sum(np.asarray(v['original_high_source']) for v in direct['groups'].values())-old_high))),
                 'group_direct_sum': float(np.max(abs(sum(np.asarray(v['direct_high_vacuum_source']) for v in direct['groups'].values())-new_high)))}
    cases = {}
    for name in ('baseline', 'mixed_normal_probe'):
        old = records['joint-constraints'][name]
        value = np.asarray(old['action_gradient_approximant'])+action[:2]
        cases[name] = {**old, 'action_gradient_approximant': value.tolist(), 'force_approximant': (-value).tolist(),
            'matter_baseline_action_gradient': (np.asarray(old['matter_baseline_action_gradient'])+action[:2]).tolist(),
            'maximum_absolute_approximant': float(max(abs(value)))}
        residuals[name+'_replacement_identity'] = float(max(abs(value-old['action_gradient_approximant']-action[:2])))
    if max(residuals.values()) > 3e-13:
        raise ArithmeticError('original high replacement composition failed')
    with _precision(60):
        a = mp.iv.sqrt(3*mp.iv.pi/2-4)
        lapse_factor, shift_factor = 8*mp.iv.pi*a, 8*mp.iv.pi*a*a
        density = mp.iv.mpf(0); current = mp.iv.mpf(0)
        for rows in physical.values():
            for row in rows:
                density += mp.iv.mpf(row['vacuum_density_error_upper'])+mp.iv.mpf(row['thermal_stress_error_upper'][0])
                current += mp.iv.mpf(row['thermal_stress_error_upper'][2])
        mode = [lapse_factor*density, shift_factor*current]
        numeric = [lapse_factor*mp.iv.mpf(aggregate['stress_error_upper'][0]),
                   shift_factor*mp.iv.mpf(aggregate['stress_error_upper'][2])]
        finite = [mode[i]+numeric[i] for i in range(2)]
        # Bound embedding the high replacement in the recorded source sum.
        # The other pieces are fixed binary approximants here, not certified
        # physical sources; their physical errors remain independently OPEN.
        exact_assembly = [sum((mp.iv.mpf(part[i]) for part in old_matter['parts'].values()), mp.iv.mpf(0))
                          -mp.iv.mpf(old_high[i])+mp.iv.mpf(new_high[i]) for i in range(4)]
        assembly_error = [abs(value-mp.iv.mpf(stress[i])) for i, value in enumerate(exact_assembly)]
        assembly_action = [lapse_factor*assembly_error[0], shift_factor*assembly_error[2]]
        tail = records['vacuum-tail-bound']
        tail_physical = [lapse_factor*(mp.iv.mpf(tail['aggregate']['vacuum_tail_error_upper'][0])+mp.iv.mpf(tail['thermal_source_bound']['decimal_upper_bounds'][0])),
                         shift_factor*mp.iv.mpf(tail['thermal_source_bound']['decimal_upper_bounds'][2])]
        budget = {'finite_high_physical_mode_and_thermal_action_error_upper': [_up_float(x) for x in mode],
                  'finite_high_quadrature_and_arithmetic_action_error_upper': [_up_float(x) for x in numeric],
                  'finite_high_total_action_error_upper': [_up_float(x) for x in finite],
                  'source_sum_rounding_action_error_upper': [_up_float(x) for x in assembly_action],
                  'finite_high_plus_source_sum_rounding_action_error_upper': [_up_float(finite[i]+assembly_action[i]) for i in range(2)],
                  'finite_high_status': 'PASS' if all(_hi(x) < mp.mpf('3e-11') for x in finite) else 'OPEN',
                  'tail_physical_mode_and_thermal_action_error_upper': [_up_float(x) for x in tail_physical],
                  'tail_numerical_integration_error_bound': None,
                  'remaining_low_subgap_source_error_bound': None,
                  'full_source_error_bound': None, 'stationarity_tolerance': 3e-11,
                  'full_source_status': 'OPEN', 'physical_pressure_error_certified': False,
                  'constraint_order': ['N', 'beta']}
    matter = {**old_matter, 'parts': {**{k: v.tolist() for k, v in old_matter['parts'].items()},
                                    'direct_high_replacement': delta.tolist()},
              'original_stress_approximant': old_matter['stress_approximant'].tolist(),
              'stress_approximant': stress.tolist(), 'original_high_source_removed': old_high.tolist(),
              'direct_high_vacuum_source_added': new_high.tolist(), 'explicit_high_replacement_delta': delta.tolist(),
              'full_source_error_bound': None,
              'unresolved': ['remaining low/subgap physical and numerical errors',
                             'tail numerical integration, distinct from physical mode error',
                             'retained angular/compact scope versus complete physical inventory']}
    paths = [f'results/development/nsc-incoming-{n}.json' for n in NAMES]
    return {'schema': 'NSC-INCOMING-SOURCE-UPDATE-v3', 'accountable_author': 'Douglas Ek',
            'status': budget['finite_high_status']+': finite high constraint-source budget joined; full source and stationarity OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in (*paths, *OWNERS)},
            'physical_coverage': physical, 'direct_group_count': 32, 'matter': matter, **cases,
            'action_gradient_increment_from_original': action.tolist(), 'partial_error_budget': budget,
            'residuals': residuals, 'composition_tolerance': 3e-13,
            'scope': {'original_finite_high_replaced_once': True, 'prior_high_corrections_added': False,
                      'LLL_state_count': 1, 'local_light_geometry_count': 1, 'thermal_declared_zero': False,
                      'new_source_law_or_coupling': False, 'new_scientific_producers_run': False,
                      'physical_initial_data_selected': False, 'metric_evolution': False,
                      'Gamma_rest_assigned': False, 'extended_stationarity': 'OPEN', 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_source_update_v3.py --check'}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write and OUTPUT.exists():
        raise FileExistsError('existing source update is not overwritten')
    data = json.loads(json.dumps(make_record()))
    if args.write:
        OUTPUT.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    elif json.loads(OUTPUT.read_text()) != data:
        raise ValueError('saved high-source composition differs')
    print(json.dumps({k: data[k] for k in ('status', 'baseline', 'partial_error_budget', 'residuals')}, indent=2))


if __name__ == '__main__':
    main()
