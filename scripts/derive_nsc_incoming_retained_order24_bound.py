#!/usr/bin/env python3
"""One group12 order24 certificate; saved-coefficient replay never prepares."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_retained_order24_bound import retained_order24_bound
from recursive_horizons.nsc_incoming_projector_energy_bound import local_cross_product_coefficients, projected_energy_bound
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _up_float

OUTPUT = ROOT/'results/development/nsc-incoming-retained-order24-bound.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_retained_order24_bound.py',
           'tests/test_nsc_incoming_retained_order24_bound.py',
           'scripts/derive_nsc_incoming_retained_order24_bound.py',
           'docs/nsc-incoming-retained-order24-bound.md')
INPUTS = ('results/development/nsc-incoming-centered-middle.json',
          'results/development/nsc-incoming-middle-bound.json',
          'results/development/nsc-mode-resolved-cauchy-state.json',
          'results/development/nsc-compact-matched-restart.json',
          'src/recursive_horizons/nsc_incoming_centered_order24.py',
          'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_defect_taylor_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py')


def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT/p).read_text())
def signature(): return {p: sha(p) for p in (*SOURCES, *INPUTS)}


def inputs():
    for path in INPUTS[:4]:
        data = read(path)
        for kind in ('source_hashes', 'input_hashes'):
            for p, digest in data.get(kind, {}).items():
                if sha(p) != digest: raise ValueError('upstream owner changed: '+p)
        if 'payload' in data and sha(data['payload']['path']) != data['payload']['sha256']:
            raise ValueError('upstream payload changed')
    baseline = read(INPUTS[0]); middle = read(INPUTS[1])['groups'][11]
    if baseline['group'] != 12 or baseline['energy_interval'] != [40., 320.]:
        raise ValueError('same group12 interval required')
    if middle['group'] != 12 or middle['middle_interval']['left'] != 40. or middle['middle_interval']['right'] != 320.:
        raise ValueError('authenticated middle metadata interval changed')
    channel = read(INPUTS[2])['channels'][12]
    if channel['index'] != 12: raise ValueError('only group12 is calculated by this record')
    return channel, read(INPUTS[3])['scattering_provenance']['config'], baseline, middle


def prepare():
    channel, config, _, _ = inputs(); before = signature()
    cross = local_cross_product_coefficients(channel, order=24)
    radial = retained_order24_bound(config, channel, progress=lambda n: print(f'group12 order24: {n}/128 cells', flush=True))
    energy = projected_energy_bound(channel, radial, cross, 40., 320.)
    if before != signature(): raise ValueError('producer changed during preparation')
    payload = {'signature': before, 'radial': radial, 'local_cross': cross, 'energy': energy,
               'channel': channel, 'energy_interval': [40., 320.]}
    raw = (json.dumps(payload, indent=2, sort_keys=True)+'\n').encode()
    digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-retained-order24-bound.{digest}.json'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-address collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def replay(path):
    channel, _, previous, middle = inputs(); payload = json.loads(path.read_text())
    if payload['signature'] != signature() or payload['channel'] != channel or payload['energy_interval'] != [40., 320.]:
        raise ValueError('prepared channel, interval or producer changed')
    radial = payload['radial']; cross = payload['local_cross']
    if radial['group'] != 12 or radial['intervals'] != 128 or radial['centered_remainder_depth'] != 4:
        raise ValueError('single group12 fixed128/depth4 certificate required')
    energy = projected_energy_bound(channel, radial, cross, 40., 320.)
    residual = max(abs(energy[k]-payload['energy'][k]) for k in energy if type(energy[k]) is float)
    excess = 0.
    for row in radial['per_sign']:
        for refined, natural in zip(row['coefficient_integral_upper'], row['natural_coefficient_integral_upper']):
            if natural: excess = max(excess, refined/natural-1)
            elif refined: raise ArithmeticError('nonzero centered bound against exact zero')
    with _precision(40):
        a = mp.iv.sqrt(3*mp.iv.pi/2-4)
        thermal = middle['thermal_scattering_stress_error_upper'][0]
        combined = _up_float(mp.iv.mpf(energy['lapse_action_error_upper'])+8*mp.iv.pi*a*thermal)
    if residual > 1e-24 or excess > 3e-15: raise ArithmeticError('coefficient replay or intersection failed')
    passes = combined <= 3e-11
    return {'schema': 'NSC-INCOMING-RETAINED-ORDER24-BOUND-v1', 'accountable_author': 'Douglas Ek',
            'status': ('PASS' if passes else 'OPEN')+': group12 energy bound conditional on explicit order24 source correction; full source OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
            'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()},
            'group': 12, 'energy_interval': [40., 320.], 'radial': radial, 'energy': energy,
            'local_cross': cross, 'thermal_scattering_density_upper': thermal,
            'conditional_combined_lapse_error_upper': combined,
            'old_order16_lapse_error_upper': previous['combined_group12_lapse_error_upper'],
            'stationarity_tolerance': 3e-11, 'fits_conditional_component_tolerance': passes,
            'residuals': {'energy_contraction_replay': residual, 'centered_excess_over_natural': max(0., excess)},
            'verification_tolerances': {'energy_contraction_replay': 1e-24, 'intersection_relative': 3e-15},
            'scope': {'explicit_order24_source_correction_required': True, 'source_correction_applied': False,
                      'archived_order16_accuracy_certified': False, 'source_quadrature_accuracy_certified': False,
                      'other_groups_calculated': False, 'low_subgap_accuracy_certified': False,
                      'constraints_solved': False, 'physical_state_or_action_changed': False,
                      'metric_evolution': False, 'Gamma_rest_assigned': False, 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_retained_order24_bound.py --check'}


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
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
    print(json.dumps({k: result[k] for k in ('status', 'energy', 'conditional_combined_lapse_error_upper', 'residuals')}, indent=2))


if __name__ == '__main__': main()
