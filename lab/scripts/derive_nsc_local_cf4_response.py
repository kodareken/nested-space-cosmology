#!/usr/bin/env python3
"""Bounded CF4 time-accuracy comparison on the already tested physical control.

Use60/120 steps, and240 only if their local matter difference exceeds1e-11.
At most three solves and360 CPU seconds. No source or spatial resampling.
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
from recursive_horizons.nsc_cf4_prepared_state import prepare_cf4_incoming
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator

OUTPUT='results/development/nsc-local-cf4-response.json'
ARTIFACT='results/development/artifacts/nsc-local-cf4-response.npz'
OWNED=('scripts/derive_nsc_local_cf4_response.py','docs/nsc-local-cf4-response.md',
        'src/recursive_horizons/nsc_cf4_prepared_state.py','tests/test_nsc_cf4_prepared_state.py')
INPUTS=(*P.OWNED,*P.INPUTS,'results/development/nsc-local-prepared-response.json')
CAP=360.
TARGET=1e-11


class BudgetReached(RuntimeError):pass
def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}


def analyze(arrays,meta):
    outputs={};rows={};changes=[]
    for steps in meta['completed']:
        name='cf4_'+str(steps)
        out=P.evaluate(arrays,name,P.ALPHA,meta['coefficients'],np.asarray(meta['baseline']),meta['parameters'])
        outputs[steps]=out;mask=np.arange(steps+1)>=2*steps//3
        rows[name]={'sampled_local_constraint_approximant_max':np.max(abs(out['residual'][mask]),axis=0).tolist(),
                    'local_matter_correction_max':np.max(abs(out['matter_correction'][mask]),axis=0).tolist(),
                    'phase_residual':meta['cases'][name]['phase_residual'],
                    'flux_residual':meta['cases'][name]['flux_residual']}
        if steps//2 in outputs:
            coarse=outputs[steps//2];start=steps//3
            value=np.max(abs(out['matter'][::2][start:]-coarse['matter'][start:]),axis=0)
            tangent=np.max(abs(out['matter_tangent'][0,::2][start:]-coarse['matter_tangent'][0,start:]),axis=0)
            changes.append({'steps':[steps//2,steps],'local_matter_change':value.tolist(),
                            'local_tangent_change':tangent.tolist(),'passes_indicator':bool(np.max(value)<=TARGET)})
    return {'per_case':rows,'time_comparisons':changes,'indicator_target':TARGET,
            'last_time_indicator_passed':bool(changes and changes[-1]['passes_indicator']),
            'same_fixed_preparation':len({c['preparation'] for c in meta['cases'].values()})<=1}


def record(arrays,meta,payload):
    info=analyze(arrays,meta)
    return {'schema':'NSC-LOCAL-CF4-RESPONSE-v1','accountable_author':'Douglas Ek',
            'status':('PASS' if info['last_time_indicator_passed'] else 'OPEN')+': CF4 time indicator; physical local gate OPEN',
            'control':{'profile':'same alpha=.001 w=T1*plateau,U=0 as midpoint control',
                       'source_family':'14_1','points':7601,'time_interval':[0.,.18],
                       'local_interval':'S(1)+[.12,.18]','steps_max':[60,120,240],
                       'method':'Blanes-Moan equation43 on the same augmented field/tangent/source matrix'},
            'measurements':info,'runtime':meta['runtime'],
            'scope':{'physical_local_gate':'OPEN','full_source_error':None,'continuum_field_error':None,
                     'between_node_error':None,'time_difference_is_error_bound':False,
                     'source_changed':False,'new_spatial_sampling':False,'stress_drift_subtracted':False,
                     'frozen_C0_imposed':False,'metric_timestep':False},
            'source_hashes':{p:meta['signature'][p] for p in OWNED},
            'input_hashes':{p:meta['signature'][p] for p in INPUTS},
            'reused_artifacts':meta['reused_artifacts'],'payload':payload,
            'reproducer':'python3 scripts/derive_nsc_local_cf4_response.py --check'}


def write(arrays,meta):
    raw=deterministic_npz_bytes({**arrays,'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)})
    (ROOT/ARTIFACT).write_bytes(raw)
    out=record(arrays,meta,{'path':ARTIFACT,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    (ROOT/OUTPUT).write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return out


def run():
    if (ROOT/OUTPUT).exists() or (ROOT/ARTIFACT).exists():raise FileExistsError('existing control; use --check')
    sig=signatures();records,a,channel=P.load_inputs();old=P.check()
    coeff=P.coefficients_from_records(records[2],records[3])
    source=FixedSourcePreparation(a['source/covariance'],a['source/column_weights'],np.repeat(a['source/energies'],3))
    arrays={k:v for k,v in a.items() if k.startswith('source/')}
    arrays['source/energies_repeated']=source.energies;arrays['reference_F0']=a['reference_columns'][0]
    meta={'signature':sig,'completed':[],'attempted':[],'cases':{},
          'coefficients':{k:np.asarray(v).tolist() for k,v in coeff.items()},
          'baseline':records[1]['baseline']['action_gradient_approximant'],
          'parameters':{'mass':channel['compact_mass'],'angular':channel['angular_eigenvalue'],
                        'axial_scale':coeff['a'],'radius':coeff['r'],'multiplicity':channel['copy_count']*channel['degeneracy']/2},
          'reused_artifacts':[records[0]['payload'],old['payload']]}
    cpu=time.process_time();wall=time.monotonic();cap_reached=False
    def stop(*_):raise BudgetReached('CF4 comparison360 CPU-second cap')
    def checkpoint():
        meta['runtime']={'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,
                         'CPU_cap':CAP,'attempted':list(meta['attempted']),'completed':list(meta['completed']),'cap_reached':cap_reached}
        if signatures()!=sig:raise ValueError('owner/input changed during CF4 comparison')
        return write(arrays,meta)
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    try:
        owner=FourthOrderModePropagator(a['x'],channel['compact_mass'],channel['angular_eigenvalue'])
        for steps in (60,120,240):
            if steps==240 and analyze(arrays,meta)['last_time_indicator_passed']:break
            name='cf4_'+str(steps);begin=time.process_time();meta['attempted'].append(steps)
            times=np.linspace(0.,.18,steps+1)
            result=prepare_cf4_incoming(owner,times,P.family(P.ALPHA),a['initial_fields'],source)
            arrays.update({name+'/times':times,name+'/F':result.state.columns,name+'/Fz':result.axial_columns,
                           name+'/dF':result.state.column_tangents,name+'/dFz':result.axial_tangents})
            meta['cases'][name]={'CPU_seconds':time.process_time()-begin,
                                'preparation':result.fixed_preparation_digest,'history':result.sampled_history_fingerprint,
                                'Gauss_history':result.gauss_fingerprint,
                                'phase_residual':result.diagnostics['exact_harmonic_phase_residual'],
                                'flux_residual':result.diagnostics['norm_flux_residual']}
            meta['completed'].append(steps);del result
            checkpoint();print(name+' complete CPU '+str(meta['cases'][name]['CPU_seconds']),flush=True)
    except BudgetReached:cap_reached=True
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    return checkpoint()


def check():
    saved=json.loads((ROOT/OUTPUT).read_text());P.load_inputs()
    if digest(ARTIFACT)!=saved['payload']['sha256']:raise ValueError('CF4 artifact changed')
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if meta['signature']!=signatures():raise ValueError('CF4 owner/input changed')
    if record(arrays,meta,saved['payload'])!=saved:raise ValueError('CF4 array-only replay differs')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    out=run() if p.parse_args().run else check()
    print(json.dumps({k:out[k] for k in ('status','measurements','runtime')},indent=2))
