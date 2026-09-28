#!/usr/bin/env python3
"""Prepare remaining source groups once; resume authenticated per-group receipts."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_high_band_batch_source import GROUPS, EXCLUDED, TOLERANCES, prepare_group_source, replay_group_source
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets

OUTPUT = 'results/development/nsc-incoming-high-band-batch-source.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_high_band_batch_source.py',
           'tests/test_nsc_incoming_high_band_batch_source.py',
           'scripts/derive_nsc_incoming_high_band_batch_source.py',
           'docs/nsc-incoming-high-band-batch-source.md')
HELPER = 'src/recursive_horizons/nsc_incoming_window_refinement.py'
HELPER_SHA = '391b201fc27b1f9a1fe00d09dc03e815d70a8e1cf00bef8cda832e05d4314c9b'
INPUTS = (HELPER,
    'results/development/nsc-incoming-spectral-source.json',
    'results/development/nsc-pg-retained-covariance.json',
    'results/development/nsc-pg-group13-covariance.json',
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-compact-matched-restart.json',
    'src/recursive_horizons/nsc_incoming_source_assembly.py',
    'src/recursive_horizons/nsc_incoming_source_tail.py',
    'src/recursive_horizons/nsc_incoming_middle_order_correction.py',
    'src/recursive_horizons/nsc_incoming_state_moments.py',
    'src/recursive_horizons/nsc_incoming_middle_bound.py',
    'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
    'src/recursive_horizons/nsc_incoming_joint_constraints.py',
    'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
    'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
)


def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def signature():
    if digest(HELPER) != HELPER_SHA: raise ValueError('frozen source-window helper changed')
    return {p:digest(p) for p in (*SOURCES,*INPUTS)}


def cache_directory():
    tag = sha256(json.dumps(signature(), sort_keys=True).encode()).hexdigest()[:20]
    path = Path(tempfile.gettempdir())/('nsc-high-band-batch-source-'+tag); path.mkdir(exist_ok=True)
    return path


def read_group(spec):
    path = ROOT/spec['path']
    if digest(spec['path']) != spec['sha256']: raise ValueError('per-group source artifact changed')
    with np.load(path, allow_pickle=False) as a: arrays = {k:a[k].copy() for k in a.files}
    metadata = json.loads(arrays['metadata_json'].tobytes())
    if metadata['signature'] != signature() or metadata['prepared']['group'] != spec['group']:
        raise ValueError('per-group producer, input or label changed')
    checked = replay_group_source(ROOT, arrays, metadata['prepared'])
    if checked != metadata['checked']: raise ValueError('per-group numerical replay differs from prepared receipt')
    return checked


def prepare_one(group):
    if group not in GROUPS: raise ValueError('excluded/completed group must not be prepared')
    directory = cache_directory(); pointer = directory/f'group{group:02d}.json'
    if pointer.exists():
        spec = json.loads(pointer.read_text()); checked = read_group(spec)
        print(f'group{group} resumed authenticated source artifact; no preparation', flush=True)
        return spec, checked
    before = signature()
    with threadpool_limits(limits=1):
        arrays, prepared = prepare_group_source(ROOT, group)
        checked = replay_group_source(ROOT, arrays, prepared)
    if signature() != before: raise ValueError('source producer/input changed while preparing group')
    arrays['metadata_json'] = np.frombuffer(json.dumps({'schema':'NSC-HIGH-BAND-SOURCE-GROUP-v1',
        'signature':before, 'prepared':prepared, 'checked':checked}, sort_keys=True, allow_nan=False).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); h = sha256(raw).hexdigest()
    relative = f'results/development/artifacts/nsc-incoming-high-band-source-group{group:02d}.{h}.npz'; path = ROOT/relative
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-addressed collision')
    if not path.exists(): path.write_bytes(raw)
    spec = {'group':group,'path':relative,'sha256':h,'bytes':len(raw)}
    pointer.write_text(json.dumps(spec, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'source_group_complete':group,'status':checked['status'],
                      'explicit_source_delta':checked['explicit_source_delta'],'artifact':relative,'failures':checked['failures']}), flush=True)
    return spec, checked


def prepare(workers):
    if workers not in (1,2): raise ValueError('one or two CPU workers only')
    print('resumable source cache: '+str(cache_directory()), flush=True)
    first_spec, first = prepare_one(6)
    if first['failures']: raise ArithmeticError('group6 source gate failed; remaining27 groups were not dispatched')
    print('GROUP6 SOURCE GATE PASS; dispatching only the remaining27 groups', flush=True)
    specs = {6:first_spec}; groups = {6:first}
    if workers == 1:
        for g in GROUPS[1:]:
            spec, result = prepare_one(g); specs[g]=spec; groups[g]=result
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            pending = {pool.submit(prepare_one,g):g for g in GROUPS[1:]}
            for future in as_completed(pending):
                spec, result = future.result(); g=spec['group']; specs[g]=spec; groups[g]=result
                print(f'source batch completed receipts {len(specs)}/28', flush=True)
    return make_record([specs[g] for g in GROUPS])


def make_record(specs):
    if [s['group'] for s in specs] != list(GROUPS): raise ValueError('exactly28 remaining groups, group6 first, required')
    groups = {str(spec['group']):read_group(spec) for spec in specs}
    if groups['6']['failures']: raise ValueError('initial group6 source gate must pass')
    windows = {}; maxima = dict.fromkeys(TOLERANCES,0.); failures = {}
    for g in GROUPS:
        result = groups[str(g)]
        if result['failures']: failures[str(g)] = result['failures']
        if set(result['windows']) != ({'middle'} if g==13 else {'low','middle'}): raise ValueError('overlap or missing source window')
        for name, row in result['windows'].items():
            windows[f'group{g}_{name}'] = row
            for key,value in row['residuals'].items(): maxima[key] = max(maxima[key],value)
    if len(windows)!=55: raise ValueError('exactly55 disjoint source windows required')
    total = sum((np.array(groups[str(g)]['explicit_source_delta']) for g in GROUPS), np.zeros(4))
    return {'schema':'NSC-INCOMING-HIGH-BAND-BATCH-SOURCE-v1','accountable_author':'Douglas Ek',
        'status':'PASS: remaining28 source groups and55 windows prepared; physical certificates separate' if not failures else 'OPEN: source numerical controls failed',
        'groups':groups,'windows':windows,'group_order':list(GROUPS),'excluded_completed_groups':list(EXCLUDED),
        'kernel_order':['rho','p_parallel','T01','p_perp'],'raw_metric_order':['N','beta','a','r'],
        'explicit_source_delta':total.tolist(),
        'additive_action_gradient':source_action_gradient(incoming_cauchy_jets(),total).tolist(),
        'maxima':maxima,'numerical_tolerances':dict(TOLERANCES),'failures':failures,
        'physical_order24_error_bound':None,'physical_error_status':'OPEN: matching independent radial batch is composed separately',
        'scope':{'groups_calculated':28,'disjoint_windows':55,'group13_has_no_overlapping_LOW_window':True,
            'group6_pass_required_before_remaining27':True,'per_group_resume_artifacts':True,'maximum_worker_processes':2,
            'new_radial_mode_or_scattering_runs':0,'old_source_arrays_or_completed_groups_changed':False,
            'physical_state_action_scales_changed':False,'thermal_physically_zero':False,'pressure_mode_error_certified':False,
            'rigorous_source_quadrature_error_certified':False,'full_source_convergence_claimed':False,
            'metric_step':False,'Gamma_rest_assigned':False},
        'group6_gate_payload':specs[0],'payloads':specs,
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
        'reproducer':'python3 scripts/derive_nsc_incoming_high_band_batch_source.py --check'}


def main():
    parser = argparse.ArgumentParser(); mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true'); mode.add_argument('--check',action='store_true')
    parser.add_argument('--workers',type=int,choices=(1,2),default=2); args=parser.parse_args()
    if args.prepare:
        result=prepare(args.workers); (ROOT/OUTPUT).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:
        prior=json.loads((ROOT/OUTPUT).read_text()); result=make_record(prior['payloads'])
        if result!=prior: raise ValueError('batch record differs from per-group authenticated replay')
    print(json.dumps({k:result[k] for k in ('status','explicit_source_delta','additive_action_gradient','maxima','failures')},indent=2))
    return 1 if result['failures'] else 0


if __name__=='__main__':raise SystemExit(main())
