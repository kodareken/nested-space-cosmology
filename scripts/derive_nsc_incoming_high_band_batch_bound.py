#!/usr/bin/env python3
"""Resumable one-pass high-band certificates, first group6 then <=3 workers."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import json
import multiprocessing
from pathlib import Path
import sys

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_high_band_batch_bound import (
    GROUPS, CachedGeometry, authenticated_batch, contract_channel, sign_symmetry_identity, interval_sign_control,
)
from recursive_horizons.nsc_incoming_projector_energy_bound import local_cross_product_coefficients
from recursive_horizons.nsc_incoming_middle_bound import INPUT_RECORDS
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _up_float

OUTPUT = ROOT/'results/development/nsc-incoming-high-band-batch-bound.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_high_band_batch_bound.py',
           'tests/test_nsc_incoming_high_band_batch_bound.py',
           'scripts/derive_nsc_incoming_high_band_batch_bound.py',
           'docs/nsc-incoming-high-band-batch-bound.md')
ARCHIVED = ('results/development/nsc-incoming-centered-order24.json',
            'results/development/nsc-incoming-retained-order24-bound.json',
            'results/development/nsc-incoming-low-high-source.json',
            'results/development/nsc-incoming-group32-order24-bound.json')
OWNERS = ('src/recursive_horizons/nsc_incoming_centered_order24.py',
          'src/recursive_horizons/nsc_incoming_defect_taylor_bound.py',
          'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py')


def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def signature(): return {p: digest(p) for p in (*SOURCES, *INPUT_RECORDS, *OWNERS, *ARCHIVED)}


def archived_sign_check():
    groups = []
    for path in ARCHIVED:
        record = json.loads((ROOT/path).read_text())
        for field in ('source_hashes', 'input_hashes'):
            for p, expected in record[field].items():
                if digest(p) != expected: raise ValueError('archived symmetry witness changed: '+p)
        radial = record.get('radial', record.get('radial_certificate'))
        first, second = radial['per_sign']
        if (first['sign'], second['sign']) != (1, -1): raise ValueError('both actual archived signs required')
        if first['coefficient_integral_upper'] != second['coefficient_integral_upper']:
            raise ArithmeticError('archived signed norm bounds differ')
        groups.append(radial['group'])
    return {'groups': groups, 'norm_array_difference': 0., 'existing_bounds_rerun': False}


def save_receipt(group, data):
    raw = (json.dumps(data, sort_keys=True, indent=2)+'\n').encode(); digest_value = sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-high-band-group{group:02d}.{digest_value}.json'
    if path.exists():
        if path.read_bytes() != raw: raise ValueError('content-address collision')
    else:
        temporary = path.with_suffix('.tmp'); temporary.write_bytes(raw); temporary.replace(path)
    return {'path': str(path.relative_to(ROOT)), 'sha256': digest_value}


def prior_receipt(group, expected):
    found = []
    for path in (ROOT/'results/development/artifacts').glob(f'nsc-incoming-high-band-group{group:02d}.*.json'):
        raw = path.read_bytes(); data = json.loads(raw)
        if data['signature'] != expected: continue
        if data['group'] != group or sha256(raw).hexdigest() != path.name.split('.')[-2]:
            raise ValueError('existing group receipt provenance mismatch')
        found.append({'path': str(path.relative_to(ROOT)), 'sha256': sha256(raw).hexdigest()})
    if len(found) > 1: raise ValueError('multiple compatible completed receipts; no arbitrary selection')
    return found[0] if found else None


_WORK = None


def worker_init(owned, expected):
    global _WORK
    _WORK = (owned, expected, CachedGeometry(owned['config']))


def calculate_group(group):
    owned, expected, geometry = _WORK
    if signature() != expected: raise ValueError('producer changed before group calculation')
    channel = owned['channels'][group]; windows = owned['windows'][str(group)]
    radial = geometry.bound(channel, progress=lambda n: print(f'group{group}: {n}/128 radial cells', flush=True))
    cross = local_cross_product_coefficients(channel, order=24)
    result = contract_channel(channel, owned['config'], windows, radial, cross)
    if signature() != expected: raise ValueError('producer changed during group calculation')
    receipt = {'schema': 'NSC-HIGH-BAND-GROUP-RECEIPT-v1', 'signature': expected,
               'group': group, 'channel': channel, 'windows_provenance': windows,
               'radial': radial, 'local_cross': cross, 'result': result,
               'radial_passes': 1, 'angular_norm_recurrences': 1,
               'source_corrections_applied': False}
    spec = save_receipt(group, receipt)
    print(json.dumps({'completed_group': group, 'status': result['status'],
                      'constraint_action_error_upper': result['aggregate_constraint_action_error_upper'],
                      'receipt': spec}), flush=True)
    return spec


def prepare(workers):
    owned = authenticated_batch(ROOT); expected = signature()
    proof = {'exact': sign_symmetry_identity(), 'interval_control': interval_sign_control(),
             'archived_norm_control': archived_sign_check()}
    receipts = {g: r for g in GROUPS if (r := prior_receipt(g, expected)) is not None}
    if 6 not in receipts:
        worker_init(owned, expected); receipts[6] = calculate_group(6)
    first = json.loads((ROOT/receipts[6]['path']).read_text())
    if max(first['result']['aggregate_constraint_action_error_upper']) <= 3e-11:
        pending = [g for g in GROUPS if g not in receipts]
        if pending:
            with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn'),
                                     initializer=worker_init, initargs=(owned, expected)) as pool:
                futures = {pool.submit(calculate_group, g): g for g in pending}
                for future in as_completed(futures): receipts[futures[future]] = future.result()
    if signature() != expected: raise ValueError('batch producer changed')
    return {'signature': expected, 'proof': proof, 'receipts': {str(g): receipts[g] for g in sorted(receipts)},
            'input_payloads': owned['input_payloads'], 'workers_max': workers,
            'group6_gate_closed': max(first['result']['aggregate_constraint_action_error_upper']) <= 3e-11}


def replay(manifest):
    expected = signature()
    if manifest['signature'] != expected: raise ValueError('batch signature changed')
    owned = authenticated_batch(ROOT)
    if manifest['input_payloads'] != owned['input_payloads']: raise ValueError('actual inputs changed')
    if manifest['proof']['exact'] != sign_symmetry_identity() or manifest['proof']['archived_norm_control'] != archived_sign_check():
        raise ValueError('sign proof differs from owned evidence')
    if manifest['proof']['interval_control']['maximum_endpoint_difference'] != 0.:
        raise ValueError('independent interval control did not pass')
    groups = []; replay_error = 0.; count = 0
    for key, spec in manifest['receipts'].items():
        if digest(spec['path']) != spec['sha256']: raise ValueError('completed receipt changed')
        receipt = json.loads((ROOT/spec['path']).read_text()); group = int(key)
        if group not in GROUPS or receipt['group'] != group or receipt['signature'] != expected:
            raise ValueError('receipt group/signature mismatch')
        channel = owned['channels'][group]; windows = owned['windows'][str(group)]
        if receipt['channel'] != channel or receipt['windows_provenance'] != windows or receipt['radial_passes'] != 1:
            raise ValueError('receipt physical inventory or single-pass claim mismatch')
        current = contract_channel(channel, owned['config'], windows, receipt['radial'], receipt['local_cross'])
        if current != receipt['result']: raise ValueError('saved contraction changed')
        groups.append({**current, 'receipt': spec}); count += len(windows)
    groups.sort(key=lambda r: r['group'])
    with _precision(40):
        total = [_up_float(sum((mp.iv.mpf(g['aggregate_constraint_action_error_upper'][i]) for g in groups), mp.iv.mpf(0))) for i in range(2)]
    complete = [g['group'] for g in groups] == list(GROUPS)
    fits = complete and max(total) <= 3e-11
    return {'schema': 'NSC-INCOMING-HIGH-BAND-BATCH-BOUND-v1', 'accountable_author': 'Douglas Ek',
            'status': ('PASS' if fits else 'OPEN')+': retained high-band energy budget conditional on explicit order24 source corrections; full source OPEN',
            'source_hashes': {p: digest(p) for p in SOURCES},
            'input_hashes': {p: digest(p) for p in (*INPUT_RECORDS, *OWNERS, *ARCHIVED)},
            'manifest': manifest, 'groups': groups, 'group_count': len(groups), 'source_window_count': count,
            'requested_groups': list(GROUPS), 'batch_complete': complete,
            'aggregate_constraint_action_error_upper': total, 'stationarity_tolerance': 3e-11,
            'within_aggregate_tolerance': fits, 'combined_intervals_added_once_per_group': True,
            'residuals': {'saved_coefficient_replay': replay_error,
                          'maximum_interval_additivity': max(g['residuals']['interval_additivity'] for g in groups),
                          'maximum_centered_excess_over_natural': max(g['residuals']['centered_excess_over_natural'] for g in groups),
                          'sign_interval_endpoint_control': manifest['proof']['interval_control']['maximum_endpoint_difference']},
            'verification_tolerances': {'saved_coefficient_replay': 0., 'interval_additivity': 1e-24,
                                        'centered_excess_over_natural': 3e-15, 'sign_interval_endpoint_control': 0.},
            'scope': {'physical_source_signs_remain_distinct': True, 'occupation_symmetry_assumed': False,
                      'completed_groups_recomputed': False, 'source_corrections_applied': False,
                      'source_quadrature_certified': False, 'modes_or_scattering_solved': 0,
                      'metric_evolution': False, 'Gamma_rest_assigned': False, 'full_source_convergence': False,
                      'constraints_solved': False, 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_high_band_batch_bound.py --check'}


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=3); args = parser.parse_args()
    if not 1 <= args.workers <= 3: raise ValueError('one to three CPU workers only')
    if args.prepare:
        if OUTPUT.exists(): raise FileExistsError('completed batch record is not overwritten')
        result = replay(prepare(args.workers)); OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    else:
        old = json.loads(OUTPUT.read_text()); result = replay(old['manifest'])
        if result != old: raise ValueError('batch replay differs')
    print(json.dumps({k: result[k] for k in ('status','group_count','source_window_count','aggregate_constraint_action_error_upper','residuals')}, indent=2))


if __name__ == '__main__': main()
