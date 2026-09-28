#!/usr/bin/env python3
"""Full authorized group32 E24–32 field panel and complete source certificate."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import sys
import tempfile

import mpmath as mp
import numpy as np
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_validated_energy_window import (
    WindowObstruction,propagate_panel_sign,replay_panel_sign,reference_integral,source_from_endpoints,
)
from recursive_horizons.nsc_incoming_window_refinement import authenticated_window
from recursive_horizons.nsc_incoming_middle_bound import thermal_middle_difference_bound,constraint_action_bound
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision
from recursive_horizons.nsc_incoming_source_quadrature_bound import float_enclosure

OUTPUT=ROOT/'results/development/nsc-incoming-validated-energy-window.json'
SOURCES=('src/recursive_horizons/nsc_incoming_validated_energy_window.py',
         'scripts/derive_nsc_incoming_validated_energy_window.py',
         'tests/test_nsc_incoming_validated_energy_window.py',
         'docs/nsc-incoming-validated-energy-window.md')
INPUTS=('results/development/nsc-incoming-energy-panel-preflight.json',
        'results/development/nsc-incoming-matched-horizon-initializer.json',
        'results/development/nsc-incoming-source-quadrature-bound.json',
        'results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-compact-matched-restart.json',
        'results/development/nsc-incoming-spectral-source.json',
        'results/development/nsc-pg-retained-covariance.json',
        'src/recursive_horizons/nsc_incoming_energy_panel.py',
        'src/recursive_horizons/nsc_incoming_validated_vacuum_propagation.py',
        'src/recursive_horizons/nsc_incoming_matched_horizon_initializer.py',
        'src/recursive_horizons/nsc_incoming_window_refinement.py',
        'src/recursive_horizons/nsc_incoming_source_assembly.py',
        'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
        'src/recursive_horizons/nsc_incoming_source_tail.py',
        'src/recursive_horizons/nsc_incoming_middle_bound.py',
        'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
        'src/recursive_horizons/nsc_incoming_state_moments.py',
        'src/recursive_horizons/nsc_compact_ctp_neck.py')


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def signature():return{p:digest(p) for p in(*SOURCES,*INPUTS)}


def cache():
    key=sha256(json.dumps(signature(),sort_keys=True).encode()).hexdigest()
    path=Path(tempfile.gettempdir())/('nsc-full-E24-32-'+key[:20]);path.mkdir(exist_ok=True)
    return path,key


def authenticated_record(path):
    record=json.loads((ROOT/path).read_text())
    for field in('source_hashes','input_hashes'):
        for p,h in record.get(field,{}).items():
            if digest(p)!=h:raise ValueError('reused owner/input changed: '+p)
    if 'payload'in record and digest(record['payload']['path'])!=record['payload']['sha256']:
        raise ValueError('reused artifact bytes changed')
    return record


def inputs():
    record=authenticated_record(INPUTS[0])
    preflight=json.loads((ROOT/record['payload']['path']).read_text())['preflight']
    if not all(c['step']['local_gate_pass'] for c in preflight['cases'].values()):raise ValueError('successful frozen preflight required')
    c=json.loads((ROOT/INPUTS[3]).read_text())['channels'][32]
    cfg=json.loads((ROOT/INPUTS[4]).read_text())['scattering_provenance']['config']
    return c,cfg,preflight


def write_artifact(prefix,data):
    raw=json.dumps(data,sort_keys=True,allow_nan=False).encode();h=sha256(raw).hexdigest()
    relative=f'results/development/artifacts/{prefix}.{h}.json';path=ROOT/relative
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content-addressed collision')
    if not path.exists():path.write_bytes(raw)
    return{'path':relative,'sha256':h,'bytes':len(raw)}


def read_sign(spec):
    if digest(spec['path'])!=spec['sha256']:raise ValueError('energy-window sign artifact changed')
    data=json.loads((ROOT/spec['path']).read_text())
    if data['signature']!=signature() or data['panel']['binding']['sign']!=spec['sign']:
        raise ValueError('full panel producer/input/sign changed')
    checked=replay_panel_sign(data['panel'])
    if checked!=data['checked']:raise ValueError('uniform field residual replay differs')
    return data['panel'],checked


def prepare_sign(sign):
    directory,key=cache();pointer=directory/f'sign{sign}.receipt.json'
    if pointer.exists():
        spec=json.loads(pointer.read_text());_,checked=read_sign(spec)
        print(json.dumps({'resumed_complete_sign':sign,'field_lapse_error_upper':checked['field_lapse_error_upper']}),flush=True)
        return spec
    c,cfg,preflight=inputs();before=signature()
    budget=preflight['future_work']['planning_allowance_seconds']/2
    print(json.dumps({'started_sign':sign,'cpu_budget_seconds':budget,'checkpoint':str(directory/f'sign{sign}.checkpoint.json')}),flush=True)
    with threadpool_limits(limits=1):
        panel=propagate_panel_sign(c,cfg,sign,preflight,directory=directory,key=key,cpu_budget=budget,
                                   progress=lambda row:print(json.dumps(row),flush=True))
        checked=replay_panel_sign(panel)
    if signature()!=before:raise ValueError('window producer/input changed during field run')
    spec=write_artifact(f'nsc-incoming-validated-energy-sign{sign}',{'signature':before,'panel':panel,'checked':checked})
    spec['sign']=sign;pointer.write_text(json.dumps(spec,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'completed_sign':sign,'field_lapse_error_upper':checked['field_lapse_error_upper'],
                      'cpu_seconds':checked['cpu_seconds'],'artifact':spec['path']}),flush=True)
    return spec


def prepare_observable(specs):
    c,cfg,_=inputs();directory,_=cache();pointer=directory/'observable.receipt.json'
    if pointer.exists():
        spec=json.loads(pointer.read_text())
        if digest(spec['path'])!=spec['sha256']:raise ValueError('cached observable artifact changed')
        data=json.loads((ROOT/spec['path']).read_text())
        if data['signature']!=signature() or data['sign_specs']!=specs:raise ValueError('cached observable belongs to changed field data')
        return spec
    before=signature();panels={str(s['sign']):read_sign(s)[0] for s in specs}
    prior=authenticated_record(INPUTS[2]);rules=json.loads((ROOT/prior['payload']['path']).read_text())
    rule=rules['windows']['group6_low']['rule']
    with _precision(80):
        references={str(sign):reference_integral(c,sign,rule) for sign in(1,-1)}
        source=source_from_endpoints(c,panels,references)
    owned=authenticated_window(ROOT,32,'low_vacuum',24.,32.)
    old=sum((owned['factor']*np.einsum('n,nv->v',v['weights'],v['kernels']) for v in owned['selected'].values()),np.zeros(4))
    if signature()!=before:raise ValueError('observable producer/input changed')
    data={'signature':before,'sign_specs':specs,'references':references,'source':source,
          'selected_panels':owned['descriptors'],'input_payloads':owned['input_payloads'],
          'old_source':old.tolist(),'thermal_bound':thermal_middle_difference_bound(c,cfg,24.,32.)}
    spec=write_artifact('nsc-incoming-validated-energy-observable',data)
    pointer.write_text(json.dumps(spec,sort_keys=True,indent=2)+'\n')
    return spec


def make_record(specs,observable_spec):
    if [s['sign'] for s in specs]!=[1,-1]:raise ValueError('both actual signs exactly once required')
    c,cfg,_=inputs();signed={str(s['sign']):read_sign(s) for s in specs}
    if digest(observable_spec['path'])!=observable_spec['sha256']:raise ValueError('window observable bytes changed')
    data=json.loads((ROOT/observable_spec['path']).read_text())
    if data['signature']!=signature() or data['sign_specs']!=specs:raise ValueError('window observable binding changed')
    owned=authenticated_window(ROOT,32,'low_vacuum',24.,32.)
    old=sum((owned['factor']*np.einsum('n,nv->v',v['weights'],v['kernels']) for v in owned['selected'].values()),np.zeros(4))
    if data['selected_panels']!=owned['descriptors'] or data['input_payloads']!=owned['input_payloads'] or data['old_source']!=old.tolist():
        raise ValueError('actual original LOW24–32 source changed')
    thermal=thermal_middle_difference_bound(c,cfg,24.,32.)
    if thermal!=data['thermal_bound']:raise ValueError('frozen positive thermal law changed')
    with _precision(80):
        source=source_from_endpoints(c,{s:p[0] for s,p in signed.items()},data['references'])
        if source!=data['source']:raise ValueError('complete source contraction differs')
        field=[sum((mp.iv.mpf(p[1]['field_stress_error_upper'][i]) for p in signed.values()),mp.iv.mpf(0)) for i in range(4)]
        raw_current_bound=float_enclosure(field[2])[1]
        field[2]=mp.iv.mpf(0)  # Explicit exact pure-vacuum current identity.
        total=[field[i]+mp.iv.mpf(source['source_arithmetic_error_upper'][i])+mp.iv.mpf(source['reference_analytic_quadrature_error_upper'][i])+mp.iv.mpf(thermal[i]) for i in range(4)]
        stress=[float_enclosure(v)[1] for v in total]
        field_bounds=[float_enclosure(v)[1] if v else 0. for v in field]
    actions=constraint_action_bound(stress);passed=max(actions)<=3e-12
    new=np.array(source['direct_vacuum_source']);delta=new-old
    return {'schema':'NSC-INCOMING-VALIDATED-ENERGY-WINDOW-v1','accountable_author':'Douglas Ek',
            'status':('PASS: complete validated group32 E24–32 vacuum source plus thermal remainder; other LOW regions OPEN' if passed else
                      'OPEN: fixed group32 energy-window enclosure exceeds3e-12'),
            'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
            'sign_payloads':specs,'observable_payload':observable_spec,
            'group':32,'energy_interval':[24,32],'selected_original_panels':owned['descriptors'],
            'kernel_order':['rho','p_parallel','T01','p_perp'],'original_source':old.tolist(),
            'validated_vacuum_source':new.tolist(),'explicit_source_delta':delta.tolist(),
            'source':source,'field_certificates':{s:p[1] for s,p in signed.items()},
            'field_stress_error_upper':field_bounds,'raw_polynomial_current_field_error_upper':raw_current_bound,
            'positive_thermal_stress_error_upper':thermal,'complete_stress_error_upper':stress,
            'constraint_action_error_upper':actions,'window_tolerance':3e-12,'complete_window_budget_pass':passed,
            'CPU_seconds':sum(p[1]['cpu_seconds'] for p in signed.values()),
            'scope':{'complete_real_energy_window_field_and_integration_proof':True,
                     'preflight_first_step_reused_each_sign':True,'pointE24_propagation_rerun':False,
                     'raw_endpoint_polynomial_normalized':False,'same_frozen_binary_mass_angular':True,
                     'same_ad4_reference_and_thermal_law':True,'old_source_or_columns_modified':False,
                     'region_disjoint_from_previous_high_windows':True,
                     'whole_source_or_joint_constraint_solution':False,'physical_IV_or_metric_timestep':False,
                     'push_or_PDF':False},
            'next_connection':'compose this explicit disjoint LOW-cell replacement once; other low and subgap source errors remain open',
            'reproducer':'python3 scripts/derive_nsc_incoming_validated_energy_window.py --check'}


def main():
    p=argparse.ArgumentParser();m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--prepare',action='store_true');m.add_argument('--check',action='store_true')
    p.add_argument('--workers',type=int,default=2);args=p.parse_args()
    if args.workers not in(1,2):raise ValueError('at most2 sign workers permitted')
    if args.prepare:
        if OUTPUT.exists():raise FileExistsError('completed window record is not overwritten')
        print('resumable full-window cache: '+str(cache()[0]),flush=True)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:specs=list(pool.map(prepare_sign,(1,-1)))
        observable=prepare_observable(specs);result=make_record(specs,observable)
        OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:
        old=json.loads(OUTPUT.read_text());result=make_record(old['sign_payloads'],old['observable_payload'])
        if result!=old:raise ValueError('complete window record replay differs')
    print(json.dumps({'status':result['status'],'source':result['validated_vacuum_source'],
                      'delta':result['explicit_source_delta'],'constraint_action_error_upper':result['constraint_action_error_upper']},indent=2))


if __name__=='__main__':main()
