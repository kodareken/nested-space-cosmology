#!/usr/bin/env python3
"""Resumable full direct high-band quadrature; fixed resolution, at most3 CPUs."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
import sys
import tempfile

from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_high_band_quadrature_batch import (
    GROUPS, prepare_direct_group, replay_direct_group, aggregate_groups,
)

OUTPUT = ROOT/'results/development/nsc-incoming-high-band-quadrature-batch.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_high_band_quadrature_batch.py',
           'scripts/derive_nsc_incoming_high_band_quadrature_batch.py',
           'tests/test_nsc_incoming_high_band_quadrature_batch.py',
           'docs/nsc-incoming-high-band-quadrature-batch.md')
HELPER = 'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py'
INPUTS = (HELPER, 'results/development/nsc-incoming-source-quadrature-bound.json',
          'results/development/nsc-incoming-vacuum-tail-bound.json',
          'results/development/nsc-incoming-spectral-source.json',
          'results/development/nsc-pg-retained-covariance.json',
          'results/development/nsc-pg-group13-covariance.json',
          'results/development/nsc-mode-resolved-cauchy-state.json',
          'results/development/nsc-compact-matched-restart.json',
          'src/recursive_horizons/nsc_incoming_window_refinement.py',
          'src/recursive_horizons/nsc_incoming_source_tail.py',
          'src/recursive_horizons/nsc_incoming_source_assembly.py',
          'src/recursive_horizons/nsc_incoming_middle_order_correction.py',
          'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
          'src/recursive_horizons/nsc_incoming_middle_bound.py',
          'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
          'src/recursive_horizons/nsc_incoming_state_moments.py',
          'src/recursive_horizons/nsc_compact_ctp_neck.py')


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def signature():
    pilot = json.loads((ROOT/'results/development/nsc-incoming-source-quadrature-bound.json').read_text())
    if digest(HELPER) != pilot['source_hashes'][HELPER]:
        raise ValueError('frozen directed source primitive changed')
    return {p: digest(p) for p in (*SOURCES, *INPUTS)}


def cache_directory():
    tag = sha256(json.dumps(signature(), sort_keys=True).encode()).hexdigest()[:20]
    path = Path(tempfile.gettempdir())/('nsc-high-band-direct-quadrature-'+tag)
    path.mkdir(exist_ok=True)
    return path


def read_group(spec):
    if spec['group'] not in GROUPS or digest(spec['path']) != spec['sha256']:
        raise ValueError('direct quadrature artifact changed')
    payload = json.loads((ROOT/spec['path']).read_text())
    if payload['signature'] != signature() or payload['prepared']['owned_union']['group'] != spec['group']:
        raise ValueError('direct quadrature producer, input or label changed')
    replayed = replay_direct_group(ROOT, payload['prepared'])
    if replayed != payload['checked']:
        raise ValueError('saved direct quadrature result differs from interval replay')
    return replayed


def prepare_one(group):
    directory = cache_directory(); pointer = directory/f'group{group:02d}.json'
    if pointer.exists():
        spec = json.loads(pointer.read_text()); checked = read_group(spec)
        print(json.dumps({'group': group, 'resumed': True, 'lapse_error_upper': checked['lapse_error_upper']}), flush=True)
        return spec, checked
    before = signature()
    def progress(row):
        print(json.dumps({'quadrature_group': group, **row}), flush=True)
    with threadpool_limits(limits=1):
        prepared = prepare_direct_group(ROOT, group, progress=progress)
        checked = replay_direct_group(ROOT, prepared)
    if signature() != before:
        raise ValueError('direct quadrature producer changed while preparing channel')
    raw = json.dumps({'schema': 'NSC-DIRECT-HIGH-BAND-QUADRATURE-GROUP-v1',
                      'signature': before, 'prepared': prepared, 'checked': checked},
                     sort_keys=True, allow_nan=False).encode()
    h = sha256(raw).hexdigest()
    relative = f'results/development/artifacts/nsc-incoming-high-band-quadrature-group{group:02d}.{h}.json'
    target = ROOT/relative
    if target.exists() and target.read_bytes() != raw:
        raise ValueError('content-addressed collision')
    if not target.exists():
        target.write_bytes(raw)
    spec = {'group': group, 'path': relative, 'sha256': h, 'bytes': len(raw)}
    pointer.write_text(json.dumps(spec, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'direct_quadrature_group_complete': group,
                      'lapse_error_upper': checked['lapse_error_upper'],
                      'direct_source': checked['direct_high_vacuum_source'],
                      'status': 'PASS' if checked['numerical_accuracy_pass'] else 'OPEN',
                      'artifact': relative}), flush=True)
    return spec, checked


def make_record(specs):
    if [s['group'] for s in specs] != list(GROUPS):
        raise ValueError('exact32-group inventory, group6 first, required')
    groups = {str(spec['group']): read_group(spec) for spec in specs}
    aggregate = aggregate_groups(groups)
    return {'schema': 'NSC-DIRECT-HIGH-BAND-QUADRATURE-BATCH-v1', 'accountable_author': 'Douglas Ek',
            'status': ('PASS: all32 direct high-band numerical source integrals certified; physical source OPEN'
                       if aggregate['numerical_budget_pass'] else 'OPEN: aggregate numerical error exceeds3e-11'),
            'source_hashes': {p: digest(p) for p in SOURCES}, 'input_hashes': {p: digest(p) for p in INPUTS},
            'group_payloads': specs, 'groups': groups, 'aggregate': aggregate,
            'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'],
            'inventory': {'groups': list(GROUPS), 'disjoint_windows': 32,
                          'early_middle_groups': [13, 14], 'early_interval': [16, 160],
                          'ordinary_lower': 32, 'upper320_groups': [10, 11, 12, 31, 32],
                          'remaining_upper': 160},
            'method': {'purpose': 'direct order24-minus-ad4, not order24-minus16',
                       'frozen_primitive_kind': 'low_vacuum', 'Gauss_nodes': 48,
                       'interval_digits': 80, 'circle_arcs': 32, 'disk_radius': 8,
                       'cell_halfwidth': 4, 'maximum_CPU_workers': 3,
                       'one_fixed_resolution_per_channel': True},
            'scope': {'radial_or_physical_mode_solves': 0, 'direct_high_source_complete': True,
                      'old16_baseline_quadrature_carried': False,
                      'old_source_arrays_or_corrections_modified': False,
                      'thermal_or_physical_projector_error_included': False,
                      'low_energy_regions_below_declared_union_included': False,
                      'infinite_tails_included': False, 'new_coupling_or_source_law': False,
                      'physical_IV_or_root_selected': False, 'metric_evolution': False,
                      'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_high_band_quadrature_batch.py --check'}


def prepare(workers):
    if workers not in (1, 2, 3):
        raise ValueError('at most3 CPU workers required')
    print('resumable direct quadrature cache: '+str(cache_directory()), flush=True)
    pilot_spec, pilot = prepare_one(6)
    if not pilot['numerical_accuracy_pass']:
        raise ArithmeticError('group6 full direct numerical pilot OPEN; other31 groups not dispatched')
    specs = {6: pilot_spec}
    print('GROUP6 FULL DIRECT NUMERICAL GATE PASS; remaining31 groups now authorized', flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(prepare_one, group): group for group in GROUPS[1:]}
        for future in as_completed(pending):
            spec, _ = future.result(); specs[spec['group']] = spec
            print(f'direct quadrature completed receipts {len(specs)}/32', flush=True)
    return make_record([specs[g] for g in GROUPS])


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--pilot', action='store_true')
    modes.add_argument('--prepare', action='store_true')
    modes.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=3)
    args = parser.parse_args()
    if args.pilot:
        prepare_one(6)
        return
    if args.prepare:
        if OUTPUT.exists():
            raise FileExistsError('existing direct quadrature batch record is not overwritten')
        result = prepare(args.workers)
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    else:
        old = json.loads(OUTPUT.read_text())
        result = make_record(old['group_payloads'])
        if old != result:
            raise ValueError('direct high-band batch replay differs')
    print(json.dumps({'status': result['status'], 'aggregate': result['aggregate']}, indent=2))


if __name__ == '__main__':
    main()
