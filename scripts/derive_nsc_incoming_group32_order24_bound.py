#!/usr/bin/env python3
"""One radial order24 pass, three explicit saved-coefficient contractions."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_retained_order24_bound import retained_order24_bound
from recursive_horizons.nsc_incoming_projector_energy_bound import local_cross_product_coefficients, projected_energy_bound
from recursive_horizons.nsc_incoming_middle_bound import (
    INPUT_RECORDS, authenticated_middle_inputs, thermal_middle_difference_bound, constraint_action_bound,
)
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _up_float

OUTPUT = ROOT/'results/development/nsc-incoming-group32-order24-bound.json'
WINDOWS = {'low': (32., 40.), 'middle': (40., 320.), 'combined': (32., 320.)}
SOURCES = ('scripts/derive_nsc_incoming_group32_order24_bound.py',
           'tests/test_nsc_incoming_group32_order24_bound.py',
           'docs/nsc-incoming-group32-order24-bound.md')
OWNERS = ('src/recursive_horizons/nsc_incoming_retained_order24_bound.py',
          'src/recursive_horizons/nsc_incoming_centered_order24.py',
          'src/recursive_horizons/nsc_incoming_defect_taylor_bound.py',
          'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py')


def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT/p).read_text())
def signature(): return {p: sha(p) for p in (*SOURCES, *INPUT_RECORDS, *OWNERS)}


def authenticated_windows():
    inherited = authenticated_middle_inputs(ROOT); channel = inherited['channels'][32]
    middle = inherited['middle_intervals'][31]
    if middle['group'] != 32 or (middle['left'], middle['right']) != WINDOWS['middle']:
        raise ValueError('actual group32 middle interval must remain[40,320]')
    finite = read('results/development/nsc-incoming-spectral-source.json')
    retained = read('results/development/nsc-pg-retained-covariance.json')
    provenance = []; maximum = 0.
    with np.load(ROOT/finite['payload']['path'], allow_pickle=False) as source, np.load(ROOT/retained['payload']['path'], allow_pickle=False) as original:
        for sign in (1, -1):
            prefix = f'low_ref/32_{sign}'
            names = [p['panel'] for p in finite['groups']['32']['signed_inventory'][str(sign)]['pieces']]
            if prefix not in names: raise ValueError('the actual selected refined LOW panel is required')
            meta = json.loads(original[prefix+'/metadata_json'].tobytes())
            if meta['channel'] != channel or meta['sign'] != sign or meta['upper'] != 40.:
                raise ValueError('owned signed LOW label or upper endpoint changed')
            if WINDOWS['low'] not in list(zip(meta['edges'], meta['edges'][1:])):
                raise ValueError('LOW[32,40] must be one actual metadata cell')
            # LOW metadata is inside the retained content-addressed payload;
            # its producer hashes are authenticated by the inherited owner.
            E, w = source[prefix+'/energies'], source[prefix+'/weights']
            if not np.array_equal(E, original[prefix+'/energy']) or not np.array_equal(w, original[prefix+'/weight']):
                raise ValueError('selected LOW nodes differ from their authenticated producer')
            mask = (E > 32)&(E < 40)
            if mask.sum() != 32 or np.any(w[mask] <= 0): raise ValueError('actual positive32-node LOW cell required')
            residual = abs(float(w[mask].sum())-8.); maximum = max(maximum, residual)
            if residual > 3e-13: raise ValueError('LOW cell measure mismatch')
            provenance.append({'sign': sign, 'panel': prefix, 'interval': list(WINDOWS['low']),
                               'rows': int(mask.sum()), 'measure_residual': residual})
    return {'channel': channel, 'config': inherited['config'], 'middle': middle,
            'low': provenance, 'low_measure_residual': maximum,
            'input_payloads': inherited['input_payloads']}


def bind_window(radial, name):
    """Explicit endpoint-guard adapter, preserving every radial coefficient."""
    if name not in WINDOWS: raise ValueError('only authenticated LOW, middle or union window allowed')
    if (radial.get('group'), radial.get('lower'), radial.get('physical_Riccati_order'),
        radial.get('intervals'), radial.get('centered_remainder_depth')) != (32, 320., 24, 128, 4):
        raise ValueError('original group32 order24 radial certificate required')
    left, right = WINDOWS[name]
    return {**radial, 'lower': right, 'original_radial_endpoint': 320.,
            'energy_contraction_window': [left, right],
            'energy_window_rebinding': 'only endpoint metadata changes; radial coefficients are energy-independent'}


def contract(owned, radial, cross):
    windows = {}
    for name, (left, right) in WINDOWS.items():
        bound = projected_energy_bound(owned['channel'], bind_window(radial, name), cross, left, right)
        thermal = thermal_middle_difference_bound(owned['channel'], owned['config'], left, right)
        thermal_action = constraint_action_bound(thermal)
        with _precision(40):
            combined = [_up_float(mp.iv.mpf(bound['lapse_action_error_upper'])+mp.iv.mpf(thermal_action[0])), thermal_action[1]]
        windows[name] = {'energy_interval': [left, right], 'vacuum_energy': bound,
                         'thermal_stress_upper': thermal, 'conditional_constraint_action_error_upper': combined,
                         'within_component_tolerance': max(combined) <= 3e-11,
                         'explicit_matching_order24_source_correction_required': True}
    return windows


def prepare():
    owned = authenticated_windows(); before = signature()
    cross = local_cross_product_coefficients(owned['channel'], order=24)
    radial = retained_order24_bound(owned['config'], owned['channel'], progress=lambda n: print(f'group32 order24: {n}/128 cells', flush=True))
    windows = contract(owned, radial, cross)
    if signature() != before: raise ValueError('producer changed during single radial pass')
    data = {'signature': before, 'provenance': owned, 'radial': radial, 'local_cross': cross, 'windows': windows,
            'radial_preparation_passes': 1}
    raw = (json.dumps(data, indent=2, sort_keys=True)+'\n').encode(); digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-group32-order24-bound.{digest}.json'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-address collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def replay(path):
    owned = authenticated_windows(); payload = json.loads(path.read_text())
    if payload['signature'] != signature() or payload['provenance'] != owned or payload['radial_preparation_passes'] != 1:
        raise ValueError('producer or single-pass window provenance changed')
    radial = payload['radial']; windows = contract(owned, radial, payload['local_cross'])
    replay_error = max(abs(windows[name]['vacuum_energy'][key]-payload['windows'][name]['vacuum_energy'][key])
                       for name in WINDOWS for key in ('density_error_upper', 'lapse_action_error_upper'))
    excess = 0.
    for row in radial['per_sign']:
        for refined, natural in zip(row['coefficient_integral_upper'], row['natural_coefficient_integral_upper']):
            if natural: excess = max(excess, refined/natural-1)
            elif refined: raise ArithmeticError('nonzero intersection against exact zero')
    additivity = max(abs(windows['combined']['vacuum_energy'][key]-windows['low']['vacuum_energy'][key]-windows['middle']['vacuum_energy'][key])
                     for key in ('density_error_upper', 'lapse_action_error_upper'))
    if replay_error > 1e-24 or additivity > 1e-24 or excess > 3e-15:
        raise ArithmeticError('coefficient replay, interval additivity or intersection failed')
    fits = windows['combined']['within_component_tolerance']
    return {'schema': 'NSC-INCOMING-GROUP32-ORDER24-BOUND-v1', 'accountable_author': 'Douglas Ek',
            'status': ('PASS' if fits else 'OPEN')+': group32 LOW and middle energy bounds conditional on matching order24 source corrections; full source OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES},
            'input_hashes': {p: sha(p) for p in (*INPUT_RECORDS, *OWNERS)},
            'input_payloads': owned['input_payloads'],
            'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()},
            'group': 32, 'radial_preparation_passes': 1, 'radial': radial,
            'local_cross': payload['local_cross'], 'windows': windows,
            'window_provenance': {'LOW': owned['low'], 'middle': owned['middle']},
            'rebindings': {name: {'original_radial_endpoint': 320., 'contraction_interval': list(window),
                                'interface_endpoint': window[1], 'radial_coefficients_changed': False}
                           for name, window in WINDOWS.items()},
            'stationarity_tolerance': 3e-11,
            'residuals': {'saved_coefficient_replay': replay_error, 'interval_additivity': additivity,
                          'centered_excess_over_natural': max(0., excess), 'LOW_cell_measure': owned['low_measure_residual']},
            'verification_tolerances': {'saved_coefficient_replay': 1e-24, 'interval_additivity': 1e-24,
                                        'centered_excess_over_natural': 3e-15, 'LOW_cell_measure': 3e-13},
            'scope': {'matching_order24_source_corrections_required': True, 'source_corrections_applied': False,
                      'archived_order16_accuracy_certified': False, 'combined_is_alternative_to_sum_of_windows': True,
                      'source_quadrature_accuracy_certified': False, 'other_groups_calculated': False,
                      'constraints_solved': False, 'physical_state_or_action_changed': False,
                      'metric_evolution': False, 'Gamma_rest_assigned': False, 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_group32_order24_bound.py --check'}


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true'); args = parser.parse_args()
    if args.prepare:
        if OUTPUT.exists(): raise FileExistsError('existing certificate is not overwritten')
        path = prepare()
    else:
        old = json.loads(OUTPUT.read_text()); path = ROOT/old['payload']['path']
        if sha(old['payload']['path']) != old['payload']['sha256']: raise ValueError('artifact changed')
    result = replay(path)
    if args.prepare: OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(old, result)
    print(json.dumps({k: result[k] for k in ('status', 'windows', 'residuals')}, indent=2))


if __name__ == '__main__': main()
