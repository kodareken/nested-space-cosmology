#!/usr/bin/env python3
"""Freeze/replay one group12 centered coefficient pilot without old solvers."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_projector_energy_bound import projected_energy_bound
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _up_float

OUTPUT = ROOT/'results/development/nsc-incoming-centered-middle.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_centered_middle.py',
           'scripts/derive_nsc_incoming_centered_middle.py',
           'tests/test_nsc_incoming_centered_middle.py', 'docs/nsc-incoming-centered-middle.md')
INPUTS = ('results/development/nsc-mode-resolved-cauchy-state.json',
          'results/development/nsc-incoming-energy-error-budget.json',
          'results/development/nsc-incoming-middle-bound.json',
          'results/development/nsc-compact-matched-restart.json',
          'src/recursive_horizons/nsc_incoming_centered_order24.py',
          'src/recursive_horizons/nsc_incoming_defect_taylor_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py')


def sha(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def read(path): return json.loads((ROOT/path).read_text())
def signature(): return {p: sha(p) for p in (*SOURCES, *INPUTS)}


def authenticate():
    for p in INPUTS[:4]:
        record = read(p)
        for key in ('source_hashes', 'input_hashes'):
            for name, digest in record.get(key, {}).items():
                if sha(name) != digest: raise ValueError('upstream owner changed: '+name)
        if 'payload' in record and sha(record['payload']['path']) != record['payload']['sha256']:
            raise ValueError('upstream coefficient artifact changed')


def prepare(directory):
    authenticate()
    controls = {name: json.loads((directory/path).read_text()) for name, path in (
        ('natural', 'nsc-group12-energy-128.json'), ('centered', 'nsc-group12-centered-middle.json'))}
    for control in controls.values():
        for p, digest in control['producer_hashes'].items():
            if sha(p) != digest: raise ValueError('pilot producer changed: '+p)
    payload = {'signature': signature(), 'controls': controls,
               'group': 12, 'same_physical_projector_order': 16}
    raw = (json.dumps(payload, indent=2, sort_keys=True)+'\n').encode()
    digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-centered-middle.{digest}.json'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-address collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def replay(path):
    authenticate()
    payload = json.loads(path.read_text())
    if payload['signature'] != signature(): raise ValueError('centered middle dependency changed')
    channel = read(INPUTS[0])['channels'][12]
    original = read(INPUTS[1])['all_group_order16_energy_bounds'][11]
    source_mid = read(INPUTS[2])['groups'][11]
    if source_mid['middle_interval']['left'] != 40. or source_mid['middle_interval']['right'] != 320.:
        raise ValueError('authenticated group12 middle interval changed')
    results = {}; maximum = 0.; monotonic = 0.; natural_coefficient_change = 0.
    for name, control in payload['controls'].items():
        radial = control['radial']
        if radial['group'] != 12 or radial['lower'] != 320. or radial['intervals'] != 128:
            raise ValueError('only fixed group12/128cell middle pilot is owned')
        energy = projected_energy_bound(channel, radial, control['local_cross'], 40., 320.)
        maximum = max(maximum, abs(energy['lapse_action_error_upper']-control['energy']['lapse_action_error_upper']))
        results[name] = energy
    natural = payload['controls']['natural']['radial']['per_sign']
    centered = payload['controls']['centered']['radial']['per_sign']
    for old, new in zip(natural, centered):
        if old['sign'] != new['sign']: raise ValueError('same angular sign order required')
        for a, b, refined in zip(old['I16_to_I32_upper'], new['natural_I16_to_I32_upper'], new['I16_to_I32_upper']):
            if a: natural_coefficient_change = max(natural_coefficient_change, abs(a-b)/a)
            elif b: raise ValueError('unexpected nonzero natural bound')
            if b: monotonic = max(monotonic, refined/b-1)
            elif refined: raise ValueError('intersection increased an exact zero')
    with _precision(40):
        thermal = source_mid['thermal_scattering_stress_error_upper']
        a = mp.iv.sqrt(3*mp.iv.pi/2-4)
        combined = _up_float(mp.iv.mpf(results['centered']['lapse_action_error_upper'])+8*mp.iv.pi*a*thermal[0])
    residuals = {'energy_contraction_replay': maximum,
                 'natural_coefficient_relative_change': natural_coefficient_change,
                 'centered_excess_over_natural': max(0., monotonic)}
    if maximum > 1e-24 or max(natural_coefficient_change, monotonic) > 3e-14:
        raise ArithmeticError('centered coefficient replay checks failed')
    return {'schema': 'NSC-INCOMING-CENTERED-MIDDLE-v1', 'accountable_author': 'Douglas Ek',
            'status': ('PASS' if combined <= 3e-11 else 'OPEN')+': group12 middle energy component; full source OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
            'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()},
            'group': 12, 'energy_interval': [40., 320.], 'same_numerical_order': 16,
            'old8_cell_energy_bound': original, 'new128_cell_energy_bounds': results,
            'thermal_scattering_stress_error_upper': thermal,
            'combined_group12_lapse_error_upper': combined, 'stationarity_tolerance': 3e-11,
            'fits_component_budget': combined <= 3e-11, 'residuals': residuals,
            'verification_tolerances': {'replay': 1e-24, 'natural_relative': 3e-14, 'intersection_relative': 3e-14},
            'scope': {'state_action_scales_changed': False, 'source_values_changed': False,
                      'mode_scattering_solves': 0, 'constraints_solved': False,
                      'low_subgap_or_quadrature_error_covered': False, 'full_source_accuracy': 'OPEN',
                      'metric_evolution': False, 'physical_NON_EXISTENCE_claimed': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_centered_middle.py --check'}


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true')
    parser.add_argument('--controls', type=Path, default=Path('/tmp')); args = parser.parse_args()
    if args.prepare:
        if OUTPUT.exists(): raise FileExistsError('existing centered-middle record is not overwritten')
        path = prepare(args.controls)
    else:
        old = read(OUTPUT.relative_to(ROOT)); path = ROOT/old['payload']['path']
        if sha(old['payload']['path']) != old['payload']['sha256']: raise ValueError('pilot artifact changed')
    result = replay(path)
    if args.prepare: OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(old, result)
    print(json.dumps({k: result[k] for k in ('status', 'combined_group12_lapse_error_upper', 'residuals')}, indent=2))


if __name__ == '__main__': main()
