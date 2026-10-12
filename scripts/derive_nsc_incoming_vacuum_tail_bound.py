#!/usr/bin/env python3
"""Prepare/replay the directed retained vacuum tail without field evolution."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_vacuum_tail_bound import (
    IncomingVacuumTailBoundGeometry, aggregate_vacuum_tail_bounds,
    vacuum_tail_from_coefficients, recurrence_cancellation_identity,
)

OUTPUT = 'results/development/nsc-incoming-vacuum-tail-bound.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
           'tests/test_nsc_incoming_vacuum_tail_bound.py',
           'scripts/derive_nsc_incoming_vacuum_tail_bound.py',
           'docs/nsc-incoming-vacuum-tail-bound.md')
INPUTS = ('results/development/nsc-mode-resolved-cauchy-state.json',
          'results/development/nsc-compact-matched-restart.json',
          'results/development/nsc-incoming-source-tail.json',
          'src/recursive_horizons/nsc_pg_high_energy.py',
          'src/recursive_horizons/nsc_lorentzian.py',
          'src/recursive_horizons/nsc_incoming_state_moments.py')


def read(path): return json.loads((ROOT/path).read_text())
def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def signature(): return {p: digest(p) for p in (*SOURCES, *INPUTS)}
def encoded(value): return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()


def inputs():
    tail = read(INPUTS[2])
    for key in ('source_hashes', 'input_hashes'):
        for path, expected in tail[key].items():
            if digest(path) != expected: raise ValueError('paired-tail dependency changed: '+path)
    payload = tail['payload']
    if digest(payload['path']) != payload['sha256']:
        raise ValueError('paired-tail artifact changed')
    channels = read(INPUTS[0])['channels']
    if len(channels) != 33 or [c['index'] for c in channels] != list(range(33)):
        raise ValueError('owned 33-group inventory including the separately treated LLL required')
    config = read(INPUTS[1])['scattering_provenance']['config']
    return channels, config, tail


def prepare(*, budget_seconds=300.):
    if not 0 < budget_seconds <= 600:
        raise ValueError('explicit preparation budget in (0,600] seconds required')
    channels, config, _ = inputs(); before = signature(); start = perf_counter()
    owner = IncomingVacuumTailBoundGeometry(config, intervals=8, precision=40)
    groups = []
    for channel in channels[1:]:
        if perf_counter()-start > budget_seconds:
            raise TimeoutError('bounded interval preparation exceeded its budget before next group')
        group = channel['index']; lower = 320. if group in (10, 11, 12, 31, 32) else 160.
        groups.append(owner.bound_channel(channel, lower))
        if group % 8 == 0: print(f'directed vacuum tail {group}/32; interval algebra only', flush=True)
    if signature() != before: raise ValueError('owner or input changed during preparation')
    payload = {'schema': 'NSC-INCOMING-VACUUM-TAIL-BOUND-ARTIFACT-v1',
               'signature': before, 'groups': groups, 'intervals': 8, 'precision': 40,
               'mpmath_version': mp.__version__,
               'finite_offset_initialization': {'horizon_offset': config['horizon_offset'],
                    'additional_projector_term_bound': None,
                    'archive_accuracy_certified': False}}
    raw = encoded(payload); h = sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-vacuum-tail-bound.{h}.json'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-addressed collision')
    if not path.exists(): path.write_bytes(raw)
    print(f'preparation seconds: {perf_counter()-start:.3f}', flush=True)
    return path


def make_record(path):
    path = Path(path); channels, config, tail = inputs()
    payload = json.loads(path.read_text())
    if payload['signature'] != signature(): raise ValueError('vacuum-tail owner or input changed')
    groups = payload['groups']; maximum_replay = 0.
    if payload['intervals'] != 8 or payload['precision'] != 40:
        raise ValueError('declared directed cell count and precision required')
    for group in groups:
        replay = vacuum_tail_from_coefficients(channels[group['group']], group['lower'], group['per_sign'],
                                              precision=payload['precision'])
        residual = max(abs(a-b) for a, b in zip(replay, group['vacuum_tail_error_upper']))
        maximum_replay = max(maximum_replay, residual)
        if replay != group['vacuum_tail_error_upper']:
            raise ValueError('coefficient-integral replay differs from recorded group bound')
        if group['lower_order_cancellation'] != recurrence_cancellation_identity():
            raise ValueError('owned lower-order recurrence cancellation changed')
    aggregate = aggregate_vacuum_tail_bounds(groups)
    failed_groups = [g['group'] for g in groups if max(g['vacuum_tail_error_upper']) > 3e-12]
    within = aggregate['within_component_budget'] and not failed_groups
    return {'schema': 'NSC-INCOMING-VACUUM-TAIL-BOUND-v1', 'accountable_author': 'Douglas Ek',
        'status': ('PASS: affine-horizon retained vacuum-tail error bounded; finite-offset archive link OPEN'
                   if within else 'OPEN: directed vacuum-tail bound exceeds the declared component budget'),
        'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'],
        'aggregate': aggregate, 'groups': groups, 'failures': failed_groups,
        'residuals': {'coefficient_bound_replay': maximum_replay,
                      'lower_order_recurrence': 0.},
        'tolerances': {'coefficient_bound_replay': 0., 'lower_order_recurrence': 0.,
                       'per_group_component_upper': 3e-12, 'all_groups_component_upper': 32*3e-12},
        'scope': {'same_fixed_incoming_surface': 'rho=1; a^2=3*pi/2-4; r^2=2',
            'vacuum_definition': 'declared affine-horizon occupied projector of the same stationary Dirac operator',
            'horizon_initial_projector_difference': 0.,
            'finite_offset_archive_horizon_offset': config['horizon_offset'],
            'finite_offset_initial_projector_term_bound': None,
            'finite_offset_archive_accuracy_certified': False,
            'stored_order16_tail_approximation_has_this_limiting_vacuum_error_bound': within,
            'finite_middle_band_accuracy_certified': False, 'subgap_accuracy_certified': False,
            'angular_compact_complement_certified': False, 'changed_normal_jets_covered': False,
            'full_source_convergence_claimed': False, 'constraints_solved': False,
            'Gamma_rest_assigned': False, 'metric_evolution': False, 'new_field_or_scattering_solves': 0,
            'seeds_stress_scales_and_archived_outputs_changed': False},
        'remaining_initial_term': {
            'definition': '2*sum_group,sign factor * integral_Emax^infinity ||V_A(E)|| * ||P_delta(E)-P_affine(E)||op dE',
            'value_or_bound': None,
            'minimal_missing_connection': 'Bound the actual finite-offset initialization projector difference with integrable energy-weighted decay, or keep the archived-mode accuracy claim separate.'},
        'paired_approximation_tail': tail['total_tail_by_order']['13'],
        'thermal_source_bound': tail['thermal_source_bound'],
        'input_hashes': {p: digest(p) for p in INPUTS},
        'source_hashes': {p: digest(p) for p in SOURCES},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': sha256(path.read_bytes()).hexdigest(),
                    'bytes': path.stat().st_size},
        'reproducer': 'python3 scripts/derive_nsc_incoming_vacuum_tail_bound.py --check'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true')
    parser.add_argument('--budget-seconds', type=float, default=300.)
    args = parser.parse_args()
    if args.prepare:
        record = make_record(prepare(budget_seconds=args.budget_seconds)); (ROOT/OUTPUT).write_bytes(encoded(record))
    else:
        previous = read(OUTPUT); item = previous['payload']; path = ROOT/item['path']
        if digest(item['path']) != item['sha256']: raise ValueError('vacuum-tail bound artifact changed')
        record = make_record(path)
        if record != previous: raise ValueError('record differs from authenticated coefficient-bound replay')
    print(json.dumps({'status': record['status'], 'aggregate': record['aggregate'],
                      'residuals': record['residuals'], 'failures': record['failures']}, indent=2))
    return 0 if not record['failures'] and record['aggregate']['within_component_budget'] else 1


if __name__ == '__main__': raise SystemExit(main())
