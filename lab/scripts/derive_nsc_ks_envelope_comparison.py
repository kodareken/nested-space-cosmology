#!/usr/bin/env python3
"""Bounded independent KS/PG comparison of the same actual prepared history.

Reuse all PG arrays, source fibers and upstream columns. At most5 envelope
solves and120 CPU seconds: reference64; alpha=.001 at64/128;256 only if
the space indicator fails; final temporal-tolerance comparison. No new PG
or spectral-source run. --check replays only stored restrictions.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_local_prepared_response as P
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_ks_source_envelope import (
    evolve_ks_source_envelope,computational_z_grid,usual_axial_support,
)
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_retarded_radial_response import archived_amplitudes

OUTPUT='results/development/nsc-ks-envelope-comparison.json'
ARTIFACT='results/development/artifacts/nsc-ks-envelope-comparison.npz'
OWNED=('scripts/derive_nsc_ks_envelope_comparison.py','docs/nsc-ks-envelope-comparison.md')
INPUTS=(*P.INPUTS,*P.OWNED,'results/development/nsc-local-cf4-response.json',
        'src/recursive_horizons/nsc_ks_source_envelope.py')
CAP=120.
TARGET=1e-11
PG_AGREEMENT=5e-9
BASE_OPTIONS={'rtol':2e-12,'atol':2e-14,'max_step':.001}
TIGHT_OPTIONS={'rtol':2e-13,'atol':2e-15,'max_step':.0005}


class BudgetReached(RuntimeError):pass
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}
def maximum(x):return float(np.max(abs(x),initial=0.))


def matter(arrays,name,parameters):
    F,Fz,dF,dFz=(arrays[name+'/'+k] for k in ('F','Fz','dF','dFz'))
    w=arrays['source/column_weights'];C=arrays['source/covariance']
    return source_column_matter(F*w,Fz*w,C,dF*w,dFz*w,**parameters)


def analyze(arrays,meta):
    rows={};values={};checks={}
    for name in meta['completed']:
        value=matter(arrays,name,meta['parameters']);values[name]=value
        rows[name]={'matter':value['action_gradient'].tolist(),
                    'matter_tangent':value['action_gradient_tangent'].tolist(),
                    'CPU_seconds':meta['cases'][name]['CPU_seconds'],
                    'function_evaluations':meta['cases'][name]['function_evaluations'],
                    'accepted_steps':meta['cases'][name]['accepted_steps'],
                    'computational_norm_residual':meta['cases'][name]['computational_norm_residual']}
    ref=matter(arrays,'analytic_reference',meta['parameters'])['action_gradient']
    if 'reference64' in values:
        drift=np.max(abs(values['reference64']['action_gradient']-ref),axis=0)
        checks['reference_matter_drift']=drift.tolist()
        checks['reference_passed']=bool(np.max(drift)<=3e-11)
    for coarse,fine in (('n64','n128'),('n128','n256')):
        if coarse in values and fine in values:
            change=np.max(abs(values[fine]['action_gradient']-values[coarse]['action_gradient']),axis=0)
            checks[coarse+'_'+fine]={'matter_change':change.tolist(),'passes_indicator':bool(np.max(change)<=TARGET)}
    if 'tight' in values:
        base=meta['tight_base']
        change=np.max(abs(values['tight']['action_gradient']-values[base]['action_gradient']),axis=0)
        checks['time_tolerance_change']=change.tolist()
        checks['time_indicator_passed']=bool(np.max(change)<=TARGET)
        pg=matter(arrays,'pg_cf4_240',meta['parameters'])
        delta=np.max(abs(values['tight']['action_gradient']-pg['action_gradient']),axis=0)
        checks['PG_CF4_matter_difference']=delta.tolist()
        checks['PG_agreement_indicator_passed']=bool(np.max(delta)<=PG_AGREEMENT)
        checks['PG_field_difference']=maximum(arrays['tight/F']-arrays['pg_cf4_240/F'])
        checks['PG_axial_difference']=maximum(arrays['tight/Fz']-arrays['pg_cf4_240/Fz'])
        checks['PG_tangent_difference']=maximum(arrays['tight/dF']-arrays['pg_cf4_240/dF'])
        checks['PG_axial_tangent_difference']=maximum(arrays['tight/dFz']-arrays['pg_cf4_240/dFz'])
        checks['local_matter_correction_max']=np.max(abs(values['tight']['action_gradient']-ref),axis=0).tolist()
    checks['same_fixed_preparation']=len({v['preparation'] for v in meta['cases'].values()})<=1
    return {'per_case':rows,'checks':checks}


def record(arrays,meta,payload):
    return {'schema':'NSC-KS-ENVELOPE-COMPARISON-v1','accountable_author':'Douglas Ek',
            'status':'OPEN: independent prepared-state method comparison; physical local gate uncertified',
            'control':{'history_amplitude':P.ALPHA,'source_family':'14_1',
                       'rho_up':meta['rho_up'],'rho_sigma':1.,'local_interval':'S(1)+[.12,.18]',
                       'computational_length':.4,'space_indicator_target':TARGET,'time_indicator_target':TARGET,
                       'reference_target':3e-11,'PG_agreement_indicator':PG_AGREEMENT,
                       'base_options':BASE_OPTIONS,'tight_options':TIGHT_OPTIONS},
            'measurements':analyze(arrays,meta),'runtime':meta['runtime'],
            'scope':{'physical_local_gate':'OPEN','full_source_error':None,'continuum_error':None,
                     'between_node_error':None,'indicators_are_error_bounds':False,
                     'new_PG_runs':0,'new_source_preparation':False,'source_changed':False,
                     'periodic_Pi_boundary_condition':False,'source_energy_used_as_outgoing_momentum':False,
                     'stress_drift_subtracted':False,'metric_timestep':False},
            'source_hashes':{p:meta['signature'][p] for p in OWNED},
            'input_hashes':{p:meta['signature'][p] for p in INPUTS},
            'reused_artifacts':meta['reused_artifacts'],'payload':payload,
            'reproducer':'python3 scripts/derive_nsc_ks_envelope_comparison.py --check'}


def publish(arrays,meta):
    raw=deterministic_npz_bytes({**arrays,'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)})
    (ROOT/ARTIFACT).write_bytes(raw)
    result=record(arrays,meta,{'path':ARTIFACT,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    (ROOT/OUTPUT).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return result


def inputs():
    records,a,channel=P.load_inputs()
    pg=json.loads((ROOT/'results/development/nsc-local-cf4-response.json').read_text())
    if digest(pg['payload']['path'])!=pg['payload']['sha256']:raise ValueError('stored PG payload changed')
    with np.load(ROOT/pg['payload']['path'],allow_pickle=False) as f:old={k:f[k] for k in f.files}
    return records,a,channel,pg,old


def run():
    if (ROOT/OUTPUT).exists() or (ROOT/ARTIFACT).exists():raise FileExistsError('existing comparison; use --check')
    sig=signatures();records,a,channel,pg,old=inputs()
    index=int(np.flatnonzero(a['x']>=1.03)[0]);rho_up=float(a['x'][index])
    initial=archived_amplitudes(a['x'],a['initial_fields'],a['source/energies'],index).transpose(1,0,2).reshape(2,12)
    source=FixedSourcePreparation(a['source/covariance'],a['source/column_weights'],np.repeat(a['source/energies'],3))
    target=a['times'][40:]+chart_coordinates(1.)[1]
    arrays={'target_z':target,'initial_canonical_columns':initial,
            'source/covariance':source.covariance,'source/column_weights':source.column_weights,
            'source/energies':source.energies}
    for key,oldkey in (('F','reference_columns'),('Fz','reference_axial_columns')):
        arrays['analytic_reference/'+key]=a[oldkey][40:]
    arrays['analytic_reference/dF']=arrays['analytic_reference/dFz']=np.zeros((0,len(target),2,12),complex)
    for key in ('F','Fz'):arrays['pg_cf4_240/'+key]=old['cf4_240/'+key][160::4]
    for key in ('dF','dFz'):arrays['pg_cf4_240/'+key]=old['cf4_240/'+key][:,160::4]
    coeff=P.coefficients_from_records(records[2],records[3])
    meta={'signature':sig,'completed':[],'attempted':[],'cases':{},'rho_up':rho_up,
          'parameters':{'mass':channel['compact_mass'],'angular':channel['angular_eigenvalue'],
                        'axial_scale':coeff['a'],'radius':coeff['r'],'multiplicity':channel['copy_count']*channel['degeneracy']/2},
          'reused_artifacts':[records[0]['payload'],pg['payload']]}
    cpu=time.process_time();wall=time.monotonic();cap=False
    def stop(*_):raise BudgetReached('KS comparison120 CPU-second cap')
    def checkpoint():
        meta['runtime']={'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,
                         'CPU_cap':CAP,'attempted':list(meta['attempted']),'completed':list(meta['completed']),'cap_reached':cap}
        if signatures()!=sig:raise ValueError('comparison owner/input changed')
        return publish(arrays,meta)
    def case(name,n,alpha,options,tangents):
        begin=time.process_time();meta['attempted'].append(name)
        result=evolve_ks_source_envelope(source,initial,P.family(alpha),computational_z_grid(n),target,
            channel['compact_mass'],channel['angular_eigenvalue'],rho_up,
            axial_support=usual_axial_support(),tangents=tangents,**options)
        for key,value in (('F',result.columns),('Fz',result.axial_columns),
                          ('dF',result.column_tangents),('dFz',result.axial_tangents)):
            arrays[name+'/'+key]=value
        meta['cases'][name]={'CPU_seconds':time.process_time()-begin,'preparation':result.fixed_preparation_digest,
                            'binding':result.binding.fingerprint,**{k:result.diagnostics[k] for k in
                            ('function_evaluations','accepted_steps','computational_norm_residual',
                             'continuum_speed_distance','padding','global_radius_bound_supplied')}}
        meta['completed'].append(name);checkpoint();print(name+' complete CPU '+str(meta['cases'][name]['CPU_seconds']),flush=True)
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    try:
        case('reference64',64,0.,BASE_OPTIONS,'zero')
        case('n64',64,P.ALPHA,BASE_OPTIONS,'all')
        case('n128',128,P.ALPHA,BASE_OPTIONS,'all')
        if not analyze(arrays,meta)['checks']['n64_n128']['passes_indicator']:
            case('n256',256,P.ALPHA,BASE_OPTIONS,'all');base='n256';n=256
        else:base='n128';n=128
        meta['tight_base']=base
        case('tight',n,P.ALPHA,TIGHT_OPTIONS,'all')
    except BudgetReached:cap=True
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    return checkpoint()


def check():
    saved=json.loads((ROOT/OUTPUT).read_text());inputs()
    if digest(ARTIFACT)!=saved['payload']['sha256']:raise ValueError('saved envelope artifact changed')
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if meta['signature']!=signatures():raise ValueError('envelope comparison owner/input changed')
    if record(arrays,meta,saved['payload'])!=saved:raise ValueError('envelope comparison replay differs')
    return saved


if __name__=='__main__':
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    result=run() if parser.parse_args().run else check()
    print(json.dumps({'status':result['status'],'checks':result['measurements']['checks'],'runtime':result['runtime']},indent=2))
