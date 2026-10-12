#!/usr/bin/env python3
"""Two512-node checks of the remaining local envelope spatial indicator."""
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
import derive_nsc_ks_envelope_comparison as C
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_source_envelope import evolve_ks_source_envelope,computational_z_grid,usual_axial_support
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-ks-envelope-resolution.json'
ARTIFACT='results/development/artifacts/nsc-ks-envelope-resolution.npz'
OWNED=('scripts/derive_nsc_ks_envelope_resolution.py','docs/nsc-ks-envelope-resolution.md')
INPUTS=(*C.OWNED,*C.INPUTS,C.OUTPUT)
TIGHTER={'rtol':5e-14,'atol':5e-16,'max_step':.00025}
CAP=30.


class BudgetReached(RuntimeError):pass
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}


def analyze(arrays,meta):
    values={name:C.matter(arrays,name,meta['parameters']) for name in ['previous256',*meta['completed']]}
    checks={}
    for first,second,label in [('previous256','n512','space'),('n512','tighter512','time')]:
        if second in values and first in values:
            delta=np.max(abs(values[second]['action_gradient']-values[first]['action_gradient']),axis=0)
            tangent=np.max(abs(values[second]['action_gradient_tangent']-values[first]['action_gradient_tangent']),axis=(0,1))
            checks[label]={'matter_change':delta.tolist(),'tangent_change':tangent.tolist(),
                           'target':1e-11,'passes_indicator':bool(np.max(delta)<=1e-11)}
    return checks


def record(arrays,meta,payload):
    checks=analyze(arrays,meta)
    passed=all(k in checks and checks[k]['passes_indicator'] for k in ('space','time'))
    return {'schema':'NSC-KS-ENVELOPE-RESOLUTION-v1','accountable_author':'Douglas Ek',
            'status':('PASS' if passed else 'OPEN')+': local envelope resolution indicators; physical gate OPEN',
            'control':{'points':512,'max_runs':2,'history_amplitude':C.P.ALPHA,
                       'local_interval':'S(1)+[.12,.18]','options':[C.TIGHT_OPTIONS,TIGHTER]},
            'measurements':checks,'runtime':meta['runtime'],
            'scope':{'physical_local_gate':'OPEN','full_source_error':None,'continuum_error':None,
                     'between_node_error':None,'indicators_are_error_bounds':False,
                     'source_changed':False,'new_PG_runs':0,'new_source_preparation':False,
                     'stress_drift_subtracted':False,'metric_timestep':False},
            'source_hashes':{p:meta['signature'][p] for p in OWNED},
            'input_hashes':{p:meta['signature'][p] for p in INPUTS},
            'reused_artifacts':meta['reused_artifacts'],'payload':payload,
            'reproducer':'python3 scripts/derive_nsc_ks_envelope_resolution.py --check'}


def publish(arrays,meta):
    raw=deterministic_npz_bytes({**arrays,'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)})
    (ROOT/ARTIFACT).write_bytes(raw)
    result=record(arrays,meta,{'path':ARTIFACT,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    (ROOT/OUTPUT).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return result


def run():
    if (ROOT/OUTPUT).exists() or (ROOT/ARTIFACT).exists():raise FileExistsError('existing resolution record; use --check')
    sig=signatures();old=C.check()
    with np.load(ROOT/old['payload']['path'],allow_pickle=False) as f:previous={k:f[k] for k in f.files}
    oldmeta=json.loads(previous['metadata_json'].tobytes())
    arrays={k:v for k,v in previous.items() if k.startswith('source/') or k in ('target_z','initial_canonical_columns')}
    for key in ('F','Fz','dF','dFz'):arrays['previous256/'+key]=previous['tight/'+key]
    source=FixedSourcePreparation(arrays['source/covariance'],arrays['source/column_weights'],arrays['source/energies'])
    meta={'signature':sig,'parameters':oldmeta['parameters'],'rho_up':oldmeta['rho_up'],
          'completed':[],'attempted':[],'cases':{},'reused_artifacts':[old['payload']]}
    cpu=time.process_time();wall=time.monotonic();cap=False
    def stop(*_):raise BudgetReached('two512-node checks30 CPU-second cap')
    def checkpoint():
        meta['runtime']={'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,
                         'CPU_cap':CAP,'completed':list(meta['completed']),'attempted':list(meta['attempted']),'cap_reached':cap}
        if signatures()!=sig:raise ValueError('resolution owner/input changed')
        return publish(arrays,meta)
    prior=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    try:
        for name,options in (('n512',C.TIGHT_OPTIONS),('tighter512',TIGHTER)):
            begin=time.process_time();meta['attempted'].append(name)
            result=evolve_ks_source_envelope(source,arrays['initial_canonical_columns'],C.P.family(C.P.ALPHA),
                computational_z_grid(512),arrays['target_z'],meta['parameters']['mass'],meta['parameters']['angular'],meta['rho_up'],
                axial_support=usual_axial_support(),tangents='all',**options)
            for key,value in (('F',result.columns),('Fz',result.axial_columns),('dF',result.column_tangents),('dFz',result.axial_tangents)):
                arrays[name+'/'+key]=value
            if result.fixed_preparation_digest!=oldmeta['cases']['tight']['preparation']:raise ValueError('fixed upstream preparation changed')
            meta['cases'][name]={'CPU_seconds':time.process_time()-begin,'preparation':result.fixed_preparation_digest,
                                'binding':result.binding.fingerprint,'function_evaluations':result.diagnostics['function_evaluations']}
            meta['completed'].append(name);checkpoint()
    except BudgetReached:cap=True
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,prior)
    return checkpoint()


def check():
    saved=json.loads((ROOT/OUTPUT).read_text());C.check()
    if digest(ARTIFACT)!=saved['payload']['sha256']:raise ValueError('resolution payload changed')
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as f:arrays={k:f[k] for k in f.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if meta['signature']!=signatures():raise ValueError('resolution source/input changed')
    if record(arrays,meta,saved['payload'])!=saved:raise ValueError('resolution array replay differs')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    result=run() if p.parse_args().run else check()
    print(json.dumps({k:result[k] for k in ('status','measurements','runtime')},indent=2))
