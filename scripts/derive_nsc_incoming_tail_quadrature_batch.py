#!/usr/bin/env python3
"""Fixed all-retained tail certificate batch, reusing the completed group22."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from hashlib import sha256
import json
import multiprocessing
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_tail_quadrature_batch import endpoint,prepare_tail,replay_tail,thermal_tail,summarize_tail
from recursive_horizons.nsc_incoming_source_quadrature_bound import pack,unpack
from recursive_horizons.nsc_incoming_middle_bound import _parameters
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision,_hi,_up_float

OUTPUT=ROOT/'results/development/nsc-incoming-tail-quadrature-batch.json'
SOURCES=('src/recursive_horizons/nsc_incoming_tail_quadrature_batch.py',
         'scripts/derive_nsc_incoming_tail_quadrature_batch.py',
         'tests/test_nsc_incoming_tail_quadrature_batch.py','docs/nsc-incoming-tail-quadrature-batch.md')
INPUTS=('results/development/nsc-incoming-tail-quadrature-bound.json',
        'results/development/nsc-incoming-source-tail.json',
        'results/development/nsc-incoming-vacuum-tail-bound.json',
        'results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-compact-matched-restart.json',
        'src/recursive_horizons/nsc_incoming_tail_quadrature_bound.py',
        'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
        'src/recursive_horizons/nsc_incoming_projector_energy_bound.py',
        'src/recursive_horizons/nsc_incoming_middle_bound.py',
        'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
        'src/recursive_horizons/nsc_compact_ctp_neck.py')


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def signature():return {p:digest(p) for p in (*SOURCES,*INPUTS)}


def inputs():
    records=[read(p) for p in INPUTS[:5]]
    for record in records:
        for field in ('source_hashes','input_hashes'):
            for p,expected in record.get(field,{}).items():
                if digest(p)!=expected:raise ValueError('tail-batch input changed: '+p)
        if 'payload' in record and digest(record['payload']['path'])!=record['payload']['sha256']:
            raise ValueError('inherited artifact changed')
    pilot,tail,vacuum,inventory,restart=records
    payload=read(pilot['payload']['path'])
    if payload['versions']!={'mpmath':mp.__version__,'sympy':sp.__version__}:raise ValueError('pilot numerical backend versions changed')
    if [c['index'] for c in inventory['channels']]!=list(range(33)):raise ValueError('actual retained inventory required')
    for group in range(1,33):
        lower=endpoint(inventory['channels'][group])
        if tail['groups'][str(group)]['lower']!=lower or vacuum['groups'][group-1]['lower']!=lower:
            raise ValueError('owned tail endpoints differ')
    return {'pilot_record':pilot,'pilot_payload':payload,'tail':tail,'vacuum':vacuum,
            'channels':inventory['channels'],'config':restart['scattering_provenance']['config']}


def save_group(group,data):
    raw=(json.dumps(data,sort_keys=True,indent=2)+'\n').encode();h=sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-incoming-tail-quadrature-group{group:02d}.{h}.json'
    if path.exists():
        if path.read_bytes()!=raw:raise ValueError('content-address collision')
    else:
        temporary=path.with_suffix('.tmp');temporary.write_bytes(raw);temporary.replace(path)
    return {'group':group,'kind':'new_group','path':str(path.relative_to(ROOT)),'sha256':h}


def saved_group(group,expected):
    found=[]
    for path in (ROOT/'results/development/artifacts').glob(f'nsc-incoming-tail-quadrature-group{group:02d}.*.json'):
        raw=path.read_bytes();data=json.loads(raw)
        if data['signature']!=expected:continue
        if data['group']!=group or sha256(raw).hexdigest()!=path.name.split('.')[-2]:raise ValueError('saved group receipt invalid')
        found.append({'group':group,'kind':'new_group','path':str(path.relative_to(ROOT)),'sha256':sha256(raw).hexdigest()})
    if len(found)>1:raise ValueError('multiple compatible group receipts')
    return found[0] if found else None


_WORK=None


def worker_init(owned,expected):
    global _WORK
    _WORK=(owned,expected)


def calculate(group):
    owned,expected=_WORK
    if group==22:raise ValueError('pilot22 must not be regenerated')
    if signature()!=expected:raise ValueError('producer changed before group')
    channel=owned['channels'][group];pilot=owned['pilot_payload']['integral']
    integral=prepare_tail(channel,pilot['rule'])
    if integral['analytic_adapter']!=pilot['analytic_adapter']:raise ValueError('analytic continuation differs from certified pilot')
    thermal=thermal_tail(channel,owned['config']);contracted=replay_tail(integral,channel)
    summary=summarize_tail(contracted,thermal,owned['tail']['groups'][str(group)]['tail_by_order']['13'],
                           owned['vacuum']['groups'][group-1]['vacuum_tail_error_upper'],channel)
    if signature()!=expected:raise ValueError('producer changed during group')
    spec=save_group(group,{'signature':expected,'group':group,'channel':channel,
                           'integral':integral,'thermal':thermal,'summary':summary})
    print(json.dumps({'completed_group':group,'status':summary['status'],
                      'combined_action_error_upper':summary['combined_action_error_upper'],'receipt':spec}),flush=True)
    return spec


def prepare(workers):
    owned=inputs();expected=signature();pilot=owned['pilot_record']['payload']
    receipts={22:{'group':22,'kind':'reused_pilot',**pilot}}
    for group in range(1,33):
        if group!=22 and (old:=saved_group(group,expected)) is not None:receipts[group]=old
    pending=[g for g in range(1,33) if g not in receipts];failures={}
    if pending:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
                                 initializer=worker_init,initargs=(owned,expected)) as pool:
            futures={pool.submit(calculate,g):g for g in pending}
            for future in as_completed(futures):
                group=futures[future]
                try:receipts[group]=future.result()
                except Exception as error:
                    failures[str(group)]=str(error)
                    for other in futures:other.cancel()
                    break
    return {'signature':expected,'receipts':{str(g):receipts[g] for g in sorted(receipts)},
            'workers_max':workers,'reused_groups':[22],'new_group_count':sum(g!=22 for g in receipts),
            'failures':failures}


def replay(manifest):
    owned=inputs()
    if manifest['signature']!=signature():raise ValueError('batch signature changed')
    groups=[];intervals={}
    for key,spec in manifest['receipts'].items():
        group=int(key)
        if digest(spec['path'])!=spec['sha256']:raise ValueError('group payload changed')
        data=read(spec['path']);channel=owned['channels'][group]
        if spec['kind']=='reused_pilot':
            if group!=22 or spec['path']!=owned['pilot_record']['payload']['path']:raise ValueError('only authenticated22 pilot may be reused')
        elif data['signature']!=signature() or data['group']!=group or data['channel']!=channel:
            raise ValueError('new group receipt inventory/signature changed')
        if data['integral']['analytic_adapter']!=owned['pilot_payload']['integral']['analytic_adapter']:
            raise ValueError('generated analytic continuation changed')
        integral=replay_tail(data['integral'],channel);thermal=thermal_tail(channel,owned['config'])
        if thermal!=data['thermal']:raise ValueError('outward thermal replay changed')
        summary=summarize_tail(integral,thermal,owned['tail']['groups'][key]['tail_by_order']['13'],
                               owned['vacuum']['groups'][group-1]['vacuum_tail_error_upper'],channel)
        if spec['kind']=='new_group' and summary!=data['summary']:raise ValueError('saved group contraction differs')
        if spec['kind']=='reused_pilot' and integral!=owned['pilot_record']['certified_integral']:
            raise ValueError('generic replay differs from frozen pilot')
        groups.append({**summary,'receipt':spec});intervals[key]=integral
    groups.sort(key=lambda g:g['group'])
    complete=[g['group'] for g in groups]==list(range(1,33))
    with _precision(80):
        numeric=[sum((unpack(g['numerical_stress_error_upper_binary'][i]) for g in groups),mp.iv.mpf(0)) for i in range(4)]
        actual_sum=[sum((mp.iv.mpf(g['archived_order13_source'][i]) for g in groups),mp.iv.mpf(0)) for i in range(4)]
        archived_total=owned['tail']['total_tail_by_order']['13']
        rounding=[mp.iv.mpf(_hi(abs(mp.iv.mpf(archived_total[i])-actual_sum[i]))) if complete else mp.iv.mpf(0) for i in range(4)]
        physical=[sum((mp.iv.mpf(g['physical_stress_error_upper'][i]) for g in groups),mp.iv.mpf(0)) for i in range(4)]
        thermal=[sum((unpack(g['thermal_tail']['intervals'][i]) for g in groups),mp.iv.mpf(0)) for i in range(4)]
        a,r,_,_,_=_parameters(owned['channels'][1])
        num_action=[4*mp.iv.pi*a*r*r*(numeric[0]+rounding[0]),4*mp.iv.pi*a*a*r*r*(numeric[2]+rounding[2])]
        full_action=[num_action[0]+4*mp.iv.pi*a*r*r*(physical[0]+thermal[0]),
                     num_action[1]+4*mp.iv.pi*a*a*r*r*(physical[2]+thermal[2])]
        up=lambda v:0. if _hi(v)==0 else _up_float(v)
        aggregate={'archived_total_order13_source':archived_total,'archived_source_unchanged':True,
                   'group_sum_rounding_stress_upper':[up(v) for v in rounding],
                   'numerical_action_error_upper':[up(v) for v in num_action],
                   'combined_action_error_upper':[up(v) for v in full_action],
                   'combined_action_upper_binary':[pack(mp.iv.mpf(_hi(v))) for v in full_action],
                   'combined_action_upper_display':[mp.nstr(_hi(v),25) for v in full_action],
                   'thermal_stress_upper_binary':[pack(mp.iv.mpf(_hi(v))) for v in thermal],
                   'thermal_stress_upper_display':[mp.nstr(_hi(v),25) for v in thermal],
                   'stationarity_tolerance':3e-11,'within_budget':complete and max(_hi(v) for v in full_action)<=mp.mpf('3e-11')}
    return {'schema':'NSC-INCOMING-TAIL-QUADRATURE-BATCH-v1','accountable_author':'Douglas Ek',
            'status':('PASS' if aggregate['within_budget'] and not manifest['failures'] else 'OPEN')+': archived retained infinite-tail numerical plus physical budget; full source OPEN',
            'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
            'manifest':manifest,'groups':groups,'group_count':len(groups),'batch_complete':complete,
            'aggregate':aggregate,'failures':manifest['failures'],
            'residuals':{'saved_integral_replay':0.,'thermal_replay':0.,'pilot_replay':0.,'exact_endpoint_certificate':0.},
            'verification_tolerances':{'saved_integral_replay':0.,'thermal_replay':0.,'pilot_replay':0.,'exact_endpoint_certificate':0.},
            'scope':{'riccati_order':16,'reference_order':4,'old_Taylor_or_radial_or_mode_generators_rerun':False,
                     'pilot22_recomputed':False,'finite_offset_low_modes_certified':False,
                     'low_subgap_or_inventory_complement_certified':False,'source_values_changed':False,
                     'constraints_solved':False,'metric_evolution':False,'Gamma_rest_assigned':False,'push_or_PDF':False},
            'reproducer':'python3 scripts/derive_nsc_incoming_tail_quadrature_batch.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--check',action='store_true')
    parser.add_argument('--workers',type=int,default=3);args=parser.parse_args()
    if not 1<=args.workers<=3:raise ValueError('at most3 CPU workers required')
    if args.prepare:
        if OUTPUT.exists():raise FileExistsError('existing batch record is not overwritten')
        result=replay(prepare(args.workers));OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:
        old=json.loads(OUTPUT.read_text());result=replay(old['manifest'])
        if old!=result:raise ValueError('batch replay differs')
    print(json.dumps({k:result[k] for k in ('status','group_count','aggregate','failures','residuals')},indent=2))


if __name__=='__main__':main()
