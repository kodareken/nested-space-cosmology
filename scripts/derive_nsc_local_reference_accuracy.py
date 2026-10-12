#!/usr/bin/env python3
"""One h=.0005 reference check on the approved I; replay never propagates.

Reuse the completed951/3801 buffer arrays. Their I drift is above3e-11 and
scales approximately with h^4. Sample missing nodes with the owned radial
equation, preserve every3801 input column, and evolve once only to tau=.18.
CPU cap90seconds. This is a partial reference indicator, never a full gate.
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
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_boundary_buffer import (
    characteristic_from_ks_path, continue_canonical_amplitudes,
    exact_harmonic_reference, maxabs, weighted_matter,
)
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator
from recursive_horizons.nsc_transmitting_history_modes import reference_fields_on_grid,characteristic_columns

OUTPUT='results/development/nsc-local-reference-accuracy.json'
ARTIFACT='results/development/artifacts/nsc-local-reference-accuracy.npz'
OWNED=('scripts/derive_nsc_local_reference_accuracy.py','docs/nsc-local-reference-accuracy.md')
INPUTS=('results/development/nsc-incoming-boundary-buffer.json',
        'results/development/nsc-retarded-response-resolution.json',
        'src/recursive_horizons/nsc_incoming_boundary_buffer.py',
        'src/recursive_horizons/nsc_transmitting_history_modes.py',
        'src/recursive_horizons/nsc_transmitting_history_jets.py',
        'src/recursive_horizons/nsc_retarded_radial_response.py',
        'src/recursive_horizons/nsc_evolved_incoming_constraints.py')
CAP=90.
TOL=3e-11


class BudgetReached(RuntimeError):pass
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}


def load_inputs():
    receipts=[json.loads((ROOT/p).read_text()) for p in INPUTS[:2]]
    data=[]
    for r in receipts:
        p=r['payload']
        if digest(p['path'])!=p['sha256']:raise ValueError('input payload changed')
        with np.load(ROOT/p['path'],allow_pickle=False) as a:data.append({k:a[k] for k in a.files})
    return receipts,data


def analysis(arrays,meta):
    if 'columns' not in arrays:return {'reference_indicator':'OPEN','sampled_local_drift':None}
    parameters=meta['parameters']
    C,weights=arrays['source/covariance'],arrays['source/column_weights']
    value=weighted_matter(arrays['columns'],arrays['pde_axial_columns'],C,weights,parameters)
    ref=weighted_matter(arrays['reference_columns'],arrays['reference_axial_columns'],C,weights,parameters)
    error=value-ref
    # Exactly the last21 nodes cover [.12,.18], with both endpoints included.
    mask=np.arange(len(arrays['times']))>=40
    local=np.max(abs(error[mask]),axis=0)
    phase=maxabs(arrays['phase']-arrays['expected_phase'])
    passed=bool(np.max(local)<=TOL and phase<=TOL and meta['prefix_bytes_equal'])
    return {'reference_indicator':'PASS' if passed else 'OPEN',
            'sampled_local_drift':local.tolist(),'local_node_count':int(mask.sum()),
            'local_time_nodes':arrays['times'][mask].tolist(),
            'phase_residual':phase,'flux_algebra':meta['flux_algebra'],
            'full_computed_drift':np.max(abs(error),axis=0).tolist(),
            'prefix_bytes_equal':meta['prefix_bytes_equal'],
            'new_sampling_common_node_residual':meta['new_sampling_common_node_residual'],
            'maximum_incoming_field_drift':maxabs(arrays['columns']-arrays['reference_columns']),
            'maximum_incoming_PDE_derivative_drift':maxabs(arrays['pde_axial_columns']-arrays['reference_axial_columns'])}


def record(arrays,meta,payload):
    return {'schema':'NSC-LOCAL-REFERENCE-ACCURACY-v1','accountable_author':'Douglas Ek',
            'status':analysis(arrays,meta)['reference_indicator']+': sampled four-energy reference control; physical local gate OPEN',
            'domain':{'rho':[-2.,1.8],'points':7601,'spacing':.0005,'PG_samples':61,
                      'local_interval':'S(1)+[.12,.18]','physical_duration':None},
            'measurements':analysis(arrays,meta),'tolerance':TOL,'runtime':meta['runtime'],
            'scope':{'full_source_error_bound':None,'continuum_field_error_bound':None,
                     'between_sample_error_bound':None,'nonzero_history_tested':False,
                     'source_law_changed':False,'columns_normalized':False,'stress_drift_subtracted':False,
                     'frozen_C0_imposed':False,'metric_timestep':False,'physical_gate':'OPEN'},
            'source_hashes':{p:meta['signature'][p] for p in OWNED},
            'input_hashes':{p:meta['signature'][p] for p in INPUTS},
            'reused_artifacts':meta['reused_artifacts'],'payload':payload,
            'reproducer':'python3 scripts/derive_nsc_local_reference_accuracy.py --check'}


def run():
    if (ROOT/OUTPUT).exists() or (ROOT/ARTIFACT).exists():raise FileExistsError('existing result; use --check')
    sig=signatures();receipts,(buffer,resolution)=load_inputs()
    meta={'signature':sig,'reused_artifacts':[r['payload'] for r in receipts],
          'parameters':json.loads(buffer['metadata_json'].tobytes())['parameters']}
    arrays={};cpu=time.process_time();wall=time.monotonic();attempts=completed=0;timed_out=False
    def stop(*_):raise BudgetReached('local reference90 CPU-second cap')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    try:
        E=buffer['source/energies'];m=float(buffer['source/mass'].item());ell=float(buffer['source/angular'].item())
        x=np.linspace(-2.,1.8,7601);x[::2]=buffer['fine/x']
        x_old=np.linspace(-2.,1.2,6401);x_old[::2]=buffer['original/x'];x[:6401]=x_old
        sampled=characteristic_columns(reference_fields_on_grid(E,m,ell,resolution['source/at_zero'],x_old))
        common=sampled.reshape(2,6401,12)[:,::2]
        res=maxabs(common-buffer['original/phi'].reshape(2,3201,12))
        meta['new_sampling_common_node_residual']=res
        if res>TOL:raise ArithmeticError('same-source new sampling disagrees with existing prefix')
        new=x[6401:]
        A,solver=continue_canonical_amplitudes(E,m,ell,buffer['extension/start_amplitudes'],new)
        ext=characteristic_from_ks_path(new,E,A)
        phi=np.concatenate((sampled.reshape(2,6401,12),ext.reshape(2,len(new),12)),axis=1)
        phi[:,::2]=buffer['fine/phi'].reshape(2,3801,12)
        meta['prefix_bytes_equal']=bool(np.array_equal(phi[:,::2],buffer['fine/phi'].reshape(2,3801,12))
                                       and np.array_equal(x[::2],buffer['fine/x']))
        phi=phi.reshape(15202,12)
        arrays.update({k:v for k,v in buffer.items() if k.startswith('source/')})
        arrays.update(x=x,initial_fields=phi,times=np.linspace(0.,.18,61))
        owner=FourthOrderModePropagator(x,m,ell);attempts=1
        result=exact_harmonic_reference(owner,phi,E,arrays['times']);completed=1
        meta['flux_algebra']=result['flux_algebra']
        arrays.update({k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    except BudgetReached:timed_out=True
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    meta['runtime']={'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,
                     'CPU_cap':CAP,'field_runs_attempted':attempts,'field_runs_completed':completed,'cap_reached':timed_out}
    if signatures()!=sig:raise ValueError('owner/input changed during computation')
    arrays['metadata_json']=np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);(ROOT/ARTIFACT).write_bytes(raw)
    out=record(arrays,meta,{'path':ARTIFACT,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    (ROOT/OUTPUT).write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return out


def check():
    saved=json.loads((ROOT/OUTPUT).read_text());receipts,(buffer,_)=load_inputs()
    if digest(ARTIFACT)!=saved['payload']['sha256']:raise ValueError('saved payload changed')
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if meta['signature']!=signatures():raise ValueError('owner/input changed')
    if 'initial_fields' in arrays:
        if not np.array_equal(arrays['initial_fields'].reshape(2,7601,12)[:,::2],buffer['fine/phi'].reshape(2,3801,12)):
            raise ValueError('old source-column bytes changed')
    if record(arrays,meta,saved['payload'])!=saved:raise ValueError('array-only replay differs')
    return saved


if __name__=='__main__':
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    out=run() if parser.parse_args().run else check()
    print(json.dumps({k:out[k] for k in ('status','measurements','runtime')},indent=2))
