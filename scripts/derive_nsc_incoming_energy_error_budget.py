#!/usr/bin/env python3
"""Freeze/replay the targeted incoming energy-error experiments."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_projector_energy_bound import projected_energy_bound, rank_one_energy_identity
from recursive_horizons.nsc_incoming_middle_bound import vacuum_middle_bound
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _up_float

OUTPUT = ROOT/'results/development/nsc-incoming-energy-error-budget.json'
SOURCES = (
    'src/recursive_horizons/nsc_incoming_defect_taylor_bound.py',
    'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
    'src/recursive_horizons/nsc_incoming_higher_order_defect.py',
    'tests/test_nsc_incoming_defect_taylor_bound.py',
    'tests/test_nsc_incoming_projector_energy_bound.py',
    'tests/test_nsc_incoming_higher_order_defect.py',
    'scripts/derive_nsc_incoming_energy_error_budget.py',
    'docs/nsc-incoming-energy-error-budget.md',
)
INPUTS = (
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-compact-matched-restart.json',
    'results/development/nsc-incoming-vacuum-tail-bound.json',
    'results/development/nsc-incoming-middle-bound.json',
    'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
    'src/recursive_horizons/nsc_incoming_middle_bound.py',
    'src/recursive_horizons/nsc_incoming_source_tail.py',
    'src/recursive_horizons/nsc_pg_high_energy.py',
)
CONTROL_FILES = {
    'natural16': 'nsc-group14-interval16-pilot.json',
    'centered4': 'nsc-group14-centered4-pilot.json',
    'point16': 'nsc-group14-point-diagnostic.json',
    'point24': 'nsc-group14-order24-point-diagnostic.json',
    'refined128_order16': 'nsc-group14-energy-128.json',
    'refined128_order24': 'nsc-group14-energy-order24-128.json',
    'all_group_baseline': 'nsc-all-group-projector-energy-baseline.json',
}


def sha(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def read(path): return json.loads((ROOT/path).read_text())
def signature(): return {p: sha(p) for p in (*SOURCES, *INPUTS)}


def authenticate_inputs():
    for path in INPUTS[:4]:
        data = read(path)
        for key in ('source_hashes', 'input_hashes'):
            for p, digest in data.get(key, {}).items():
                if sha(p) != digest: raise ValueError('upstream owner changed: '+p)
        if 'payload' in data and sha(data['payload']['path']) != data['payload']['sha256']:
            raise ValueError('upstream coefficient artifact changed')


def prepare(directory):
    authenticate_inputs()
    data = {name: json.loads((directory/path).read_text()) for name, path in CONTROL_FILES.items()}
    for p, digest in data['refined128_order24']['producer_hashes'].items():
        if sha(p) != digest: raise ValueError('order24 producer changed during calculation: '+p)
    payload = {'replay_signature': signature(), 'controls': data,
               'control_file_hashes': {k: hashlib.sha256((directory/p).read_bytes()).hexdigest() for k, p in CONTROL_FILES.items()},
               'point_diagnostics_are_bounds': False,
               'coefficient_generation': 'targeted group14 interval and numerical-order experiments; no mode or scattering solve'}
    raw = (json.dumps(payload, sort_keys=True, indent=2)+'\n').encode()
    digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-energy-error-budget.{digest}.json'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-address collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def replay(path):
    authenticate_inputs()
    payload = json.loads(path.read_text())
    if payload['replay_signature'] != signature(): raise ValueError('energy replay dependency changed')
    controls = payload['controls']
    channels = read(INPUTS[0])['channels']
    old = read(INPUTS[2])['groups']; middle = read(INPUTS[3])
    baseline = controls['all_group_baseline']
    if [row['group'] for row in baseline] != list(range(1, 33)):
        raise ValueError('all32 retained baseline groups required')
    residual = 0.; groups = []
    for row, radial in zip(baseline, old):
        group = row['group']; left = 16. if group in (13, 14) else 40.
        current = projected_energy_bound(channels[group], radial, row['local_cross'], left, radial['lower'])
        residual = max(residual, abs(current['lapse_action_error_upper']-row['energy']['lapse_action_error_upper']))
        groups.append({'group': group, **current})
    cross16 = baseline[13]['local_cross']; cases = {}
    source_cases = (
        ('natural8_order16', old[13], cross16),
        ('natural16_order16', controls['natural16']['group14_intervals16'], cross16),
        ('centered4_order16', controls['centered4']['bound'], cross16),
        ('natural128_order16', controls['refined128_order16']['radial'], cross16),
        ('natural128_order24', controls['refined128_order24']['radial'], controls['refined128_order24']['local_cross']),
    )
    for name, radial, cross in source_cases:
        energy = projected_energy_bound(channels[14], radial, cross, 16., 160.)
        cases[name] = {'numerical_order': cross['physical_Riccati_order'], 'radial_intervals': radial['intervals'],
                       'energy': energy, 'within_single_component_tolerance': energy['lapse_action_error_upper'] <= 3e-11}
        if cross['physical_Riccati_order'] == 16:
            cases[name]['generic_stress_bound'] = vacuum_middle_bound(channels[14], radial, 16., 160.)
    comparison = controls['refined128_order24']['energy']
    residual = max(residual, abs(cases['natural128_order24']['energy']['lapse_action_error_upper']-comparison['lapse_action_error_upper']))
    excess = 0.
    for row in controls['centered4']['bound']['per_sign']:
        for centered, natural in zip(row['I16_to_I32_upper'], row['natural_I16_to_I32_upper']):
            if natural: excess = max(excess, centered/natural-1)
            elif centered: raise ValueError('nonzero intersection bound against an exact zero')
    with _precision(40):
        old_sum = sum((mp.iv.mpf(row['lapse_action_error_upper']) for row in groups), mp.iv.mpf(0))
        other = sum((mp.iv.mpf(row['lapse_action_error_upper']) for row in groups if row['group'] != 14), mp.iv.mpf(0))
        new = mp.iv.mpf(cases['natural128_order24']['energy']['lapse_action_error_upper'])
        thermal = mp.iv.mpf(middle['aggregate']['thermal_scattering_stress_error_upper'][0])
        a = mp.iv.sqrt(3*mp.iv.pi/2-4)
        combined = _up_float(other+new+8*mp.iv.pi*a*thermal)
        aggregate = {'unchanged_order16_lapse_bound': _up_float(old_sum),
                     'other31_groups_lapse_bound': _up_float(other),
                     'conditional_order24_group14_plus_other31_and_thermal_lapse_bound': combined,
                     'fits_stationarity_tolerance': combined <= 3e-11,
                     'stationarity_tolerance': 3e-11,
                     'group14_source_correction_must_be_composed_separately': True}
    residuals = {'stored_bound_replay': residual, 'centered_bound_excess_over_natural': max(0., excess),
                 'rank_one_identity': 0.}
    if residual > 1e-24 or excess > 3e-15:
        raise ArithmeticError('energy bound replay or intersection failed')
    passes14 = cases['natural128_order24']['within_single_component_tolerance']
    return {'schema': 'NSC-INCOMING-ENERGY-ERROR-BUDGET-v1', 'accountable_author': 'Douglas Ek',
            'status': ('PASS' if passes14 else 'OPEN')+': group14 order24 vacuum energy component; complete source OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
            'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()},
            'identity': rank_one_energy_identity(), 'group14_cases': cases,
            'all_group_order16_energy_bounds': groups, 'aggregate': aggregate, 'residuals': residuals,
            'verification_tolerances': {'stored_bound_replay': 1e-24, 'centered_bound_excess': 3e-15, 'exact_identity': 0.},
            'point_diagnostic_scope': 'exploratory samples only; not used in any certified bound',
            'scope': {'physical_Riccati_order16_preserved': True, 'order24_new_approximation_explicit': True,
                      'physical_state_or_action_changed': False, 'source_correction_applied': False,
                      'full_pressure_accuracy_certified': False, 'low_subgap_source_accuracy_certified': False,
                      'constraints_solved': False, 'metric_evolution': False, 'Gamma_rest_assigned': False,
                      'old_mode_generators_rerun': False, 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_energy_error_budget.py --check'}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true')
    parser.add_argument('--controls', type=Path, default=Path('/tmp'))
    args = parser.parse_args()
    if args.prepare:
        if OUTPUT.exists(): raise FileExistsError('existing energy-error record is not overwritten')
        path = prepare(args.controls)
    else:
        old = json.loads(OUTPUT.read_text()); path = ROOT/old['payload']['path']
        if sha(old['payload']['path']) != old['payload']['sha256']: raise ValueError('energy artifact changed')
    data = replay(path)
    if args.prepare: OUTPUT.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(old, data)
    print(json.dumps({k: data[k] for k in ('status', 'group14_cases', 'aggregate', 'residuals')}, indent=2))


if __name__ == '__main__': main()
