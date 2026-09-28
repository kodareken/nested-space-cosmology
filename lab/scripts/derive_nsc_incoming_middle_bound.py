#!/usr/bin/env python3
"""Evaluate/replay finite middle bounds; no radial or field generator."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_middle_bound import (
    INPUT_RECORDS, STATIONARITY_TOLERANCE, authenticated_middle_inputs,
    vacuum_middle_bound, thermal_middle_difference_bound, positive_sum,
    constraint_action_bound,
)

OUTPUT = 'results/development/nsc-incoming-middle-bound.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_middle_bound.py',
           'tests/test_nsc_incoming_middle_bound.py',
           'scripts/derive_nsc_incoming_middle_bound.py',
           'docs/nsc-incoming-middle-bound.md')
OWNERS = ('src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_joint_constraints.py',
          'scripts/derive_nsc_incoming_source_tail.py',
          'src/recursive_horizons/nsc_pg_archived_high_energy_modes.py',
          'src/recursive_horizons/nsc_incoming_high_energy_source.py')


def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()


def make_record():
    inputs = authenticated_middle_inputs(ROOT); groups = []
    for channel, radial, interval in zip(inputs['channels'][1:], inputs['radial_groups'], inputs['middle_intervals']):
        group = channel['index']; left, right = interval['left'], interval['right']
        if interval['group'] != group: raise ValueError('middle interval and actual retained group differ')
        vacuum = vacuum_middle_bound(channel, radial, left, right)
        thermal = thermal_middle_difference_bound(channel, inputs['config'], left, right)
        total = positive_sum([vacuum, thermal]); action = constraint_action_bound(total)
        groups.append({'group': group, 'middle_interval': interval,
                       'vacuum_stress_error_upper': vacuum,
                       'thermal_scattering_stress_error_upper': thermal,
                       'combined_stress_error_upper': total,
                       'constraint_action_error_upper': action,
                       'fits_stationarity_tolerance': max(action) <= STATIONARITY_TOLERANCE})
    if [g['group'] for g in groups] != list(range(1, 33)):
        raise ValueError('all32 retained middle bands required')
    vacuum = positive_sum([g['vacuum_stress_error_upper'] for g in groups])
    thermal = positive_sum([g['thermal_scattering_stress_error_upper'] for g in groups])
    combined = positive_sum([vacuum, thermal]); action = constraint_action_bound(combined)
    fits = max(action) <= STATIONARITY_TOLERANCE
    return {'schema': 'NSC-INCOMING-MIDDLE-BOUND-v1', 'accountable_author': 'Douglas Ek',
        'status': ('PASS: middle-band enclosure fits existing constraint tolerance; full source OPEN' if fits else
                   'OPEN: valid middle-band enclosure is too broad for existing constraint tolerance'),
        'valid_enclosure': True, 'physical_NON_EXISTENCE_claimed': False,
        'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'], 'constraint_order': ['N', 'beta'],
        'groups': groups, 'group_count': 32,
        'aggregate': {'vacuum_stress_error_upper': vacuum,
                      'thermal_scattering_stress_error_upper': thermal,
                      'combined_stress_error_upper': combined,
                      'constraint_action_error_upper': action,
                      'stationarity_tolerance': STATIONARITY_TOLERANCE,
                      'fits_stationarity_tolerance': fits,
                      'ratio_to_stationarity_tolerance': [v/STATIONARITY_TOLERANCE for v in action]},
        'density_bound_ranking': [{'group': g['group'], 'density_error_upper': g['combined_stress_error_upper'][0],
                                  'lapse_action_error_upper': g['constraint_action_error_upper'][0]}
                                 for g in sorted(groups, key=lambda g: g['combined_stress_error_upper'][0], reverse=True)],
        'residuals': {'middle_endpoint_measure': inputs['endpoint_measure_residual'],
                      'selected_panel_endpoint_or_label_mismatches': 0},
        'verification_tolerances': {'middle_endpoint_measure': 3e-11,
                                    'selected_panel_endpoint_or_label_mismatches': 0},
        'thermal_bound': {'trace_norm_of_each_departure': '2*exp(-pi*E/kappa)+exp(-2*pi*E/(Omega*kappa))',
                          'difference_rule': 'triangle bound: actual and archived-recipe departures each satisfy b(E), hence difference <=2*b(E)',
                          'source_horizon_coherence_retained_in_bound': True,
                          'finite_exponential_primitives': True},
        'scope': {'radial_coefficient_preparation_rerun': False,
                  'field_or_scattering_solve': False, 'new_energy_partition': False,
                  'same_affine_horizon_vacuum_definition': True,
                  'middle_recipe_uses_finite_offset_initialization': False,
                  'finite_start_error_assigned_to_infinite_tail': False,
                  'low_subgap_modal_accuracy_covered': False, 'quadrature_error_covered': False,
                  'angular_compact_complement_covered': False, 'changed_normal_jets_covered': False,
                  'full_source_convergence_claimed': False, 'constraint_root_selected': False,
                  'metric_evolution': False, 'Gamma_rest_assigned': False,
                  'scales_seeds_state_stress_or_archived_outputs_changed': False},
        'next_gap': ('Low-panel quadrature and modal accuracy, remaining subgap refinements, and complete source matching remain.' if fits else
                     'The interval/Riccati defect enclosure on the ranked finite middle bands is too broad. A sharper physical approximation bound is required; the enclosure does not establish the actual error. Low-panel quadrature and modal accuracy remain independent gaps.'),
        'input_payloads': inputs['input_payloads'],
        'input_hashes': {p: digest(p) for p in (*INPUT_RECORDS, *OWNERS)},
        'source_hashes': {p: digest(p) for p in SOURCES},
        'reproducer': 'python3 scripts/derive_nsc_incoming_middle_bound.py --check'}


def main():
    parser = argparse.ArgumentParser(description=__doc__); modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write', action='store_true'); modes.add_argument('--check', action='store_true')
    args = parser.parse_args(); result = make_record()
    if args.write: (ROOT/OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    elif result != json.loads((ROOT/OUTPUT).read_text()):
        raise ValueError('middle bound differs from authenticated finite-primitive replay')
    print(json.dumps({k: result[k] for k in ('status', 'aggregate', 'residuals')}, indent=2))
    # OPEN is an honest physical error-budget result, not failed replay.
    return 0


if __name__ == '__main__': raise SystemExit(main())
