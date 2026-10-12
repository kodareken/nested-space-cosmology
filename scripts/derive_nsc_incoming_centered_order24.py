#!/usr/bin/env python3
"""Prepare once; replay centered group14 order24 coefficients without recurrence."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_centered_order24 import CenteredOrder24Geometry
from recursive_horizons.nsc_incoming_projector_energy_bound import projected_energy_bound

OUTPUT = ROOT/'results/development/nsc-incoming-centered-order24.json'
SOURCES = (
    'src/recursive_horizons/nsc_incoming_centered_order24.py',
    'tests/test_nsc_incoming_centered_order24.py',
    'scripts/derive_nsc_incoming_centered_order24.py',
    'docs/nsc-incoming-centered-order24.md',
)
INPUTS = (
    'results/development/nsc-incoming-energy-error-budget.json',
    'results/development/nsc-incoming-middle-order-correction.json',
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-compact-matched-restart.json',
    'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
    'src/recursive_horizons/nsc_incoming_defect_taylor_bound.py',
    'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
    'src/recursive_horizons/nsc_incoming_middle_bound.py',
)


def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT/p).read_text())
def signature(): return {p: sha(p) for p in (*SOURCES, *INPUTS)}


def inputs():
    for path in INPUTS[:4]:
        record = read(path)
        for kind in ('source_hashes', 'input_hashes'):
            for p, digest in record.get(kind, {}).items():
                if sha(p) != digest: raise ValueError('upstream owner changed: '+p)
        if 'payload' in record and sha(record['payload']['path']) != record['payload']['sha256']:
            raise ValueError('upstream payload changed')
    energy = read(INPUTS[0]); payload = read(energy['payload']['path'])
    previous = payload['controls']['refined128_order24']
    channel = read(INPUTS[2])['channels'][14]
    config = read(INPUTS[3])['scattering_provenance']['config']
    return channel, config, previous, energy['payload']


def prepare():
    channel, config, previous, inherited = inputs(); before = signature()
    owner = CenteredOrder24Geometry(config)
    radial = owner.bound_channel(channel, progress=lambda n: print(f'group14 centered order24: {n}/128 cells', flush=True))
    if before != signature(): raise ValueError('dependency changed during certificate')
    data = {'signature': before, 'radial': radial, 'local_cross': previous['local_cross'],
            'previous_natural_energy': previous['energy'], 'inherited_energy_payload': inherited}
    raw = (json.dumps(data, sort_keys=True, indent=2)+'\n').encode()
    digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-centered-order24.{digest}.json'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-address collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def replay(path):
    channel, _, previous, inherited = inputs()
    payload = json.loads(path.read_text())
    if payload['signature'] != signature(): raise ValueError('centered certificate dependency changed')
    if payload['local_cross'] != previous['local_cross'] or payload['inherited_energy_payload'] != inherited:
        raise ValueError('unchanged local cross coefficients required')
    radial = payload['radial']
    if radial['physical_Riccati_order'] != 24 or radial['centered_remainder_depth'] != 4 or radial['intervals'] != 128:
        raise ValueError('fixed order24 depth4 partition128 required')
    energy = projected_energy_bound(channel, radial, payload['local_cross'], 16., 160.)
    natural = {**radial, 'per_sign': [{'sign': r['sign'], 'coefficient_integral_upper': r['natural_coefficient_integral_upper']}
                                    for r in radial['per_sign']]}
    same_natural = projected_energy_bound(channel, natural, payload['local_cross'], 16., 160.)
    excess = max(0., max(c/n-1 for row in radial['per_sign'] for c, n in zip(
        row['coefficient_integral_upper'], row['natural_coefficient_integral_upper']) if n))
    replay_residual = abs(same_natural['lapse_action_error_upper']-previous['energy']['lapse_action_error_upper'])
    if excess > 3e-15 or replay_residual > 1e-24:
        raise ArithmeticError('natural-owner agreement or intersection failed')
    passes = energy['lapse_action_error_upper'] <= 3e-11
    return {'schema': 'NSC-INCOMING-CENTERED-ORDER24-v1', 'accountable_author': 'Douglas Ek',
            'status': ('PASS' if passes else 'OPEN')+': group14 order24 vacuum energy component; complete source OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
            'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()},
            'inherited_energy_payload': inherited, 'radial': radial, 'energy': energy,
            'previous_natural_energy': previous['energy'],
            'lapse_bound_improvement_factor': previous['energy']['lapse_action_error_upper']/energy['lapse_action_error_upper'],
            'stationarity_tolerance': 3e-11, 'within_single_component_tolerance': passes,
            'residuals': {'natural_owner_lapse_agreement': replay_residual, 'centered_excess_over_natural': excess},
            'verification_tolerances': {'natural_owner_lapse_agreement': 1e-24, 'centered_excess_over_natural': 3e-15},
            'scope': {'unchanged_explicit_order24_projector': True, 'group14_source_correction_applied': False,
                      'thermal_or_quadrature_bound_included': False, 'other31_groups_included': False,
                      'full_source_convergence': False, 'constraints_solved': False,
                      'metric_evolution': False, 'Gamma_rest_assigned': False, 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_centered_order24.py --check'}


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.prepare:
        if OUTPUT.exists(): raise FileExistsError('existing certificate is not overwritten')
        path = prepare()
    else:
        old = json.loads(OUTPUT.read_text()); path = ROOT/old['payload']['path']
        if sha(old['payload']['path']) != old['payload']['sha256']: raise ValueError('coefficient artifact changed')
    result = replay(path)
    if args.prepare: OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(old, result)
    print(json.dumps({k: result[k] for k in ('status', 'energy', 'lapse_bound_improvement_factor', 'residuals')}, indent=2))


if __name__ == '__main__': main()
