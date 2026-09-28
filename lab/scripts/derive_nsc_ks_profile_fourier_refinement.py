#!/usr/bin/env python3
"""Sharper geometric Fourier remainder; no source or field evolution."""
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
import flint
from flint import arb,acb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_ks_profile_fourier_bound as C
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily
from recursive_horizons.nsc_ks_profile_fourier_bound import enclose_profile_fourier,alias_and_tail_bounds
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-ks-profile-fourier-refinement.json'
PAYLOAD='results/development/artifacts/nsc-ks-profile-fourier-refinement.npz'
OWNED=('scripts/derive_nsc_ks_profile_fourier_refinement.py','docs/nsc-ks-profile-fourier-refinement.md')
INPUTS=tuple(dict.fromkeys((*C.OWNED,*C.INPUTS,C.OUTPUT)))
CAP,BITS=120.,120
SETTINGS={'derivative_order':16,'transition_panels':1024,'quadrature_points':65536,
          'retained_index':8192,'max_derivative':2,'bits':BITS}


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}


def measure(arrays,meta):
    rows={}
    with ctx.workprec(BITS):
        for key,data in meta['profiles'].items():
            B=restored_upper(data['derivative_l1_upper'])
            alias,tails=alias_and_tail_bounds(B,arb(meta['period_length']),SETTINGS['derivative_order'],
                SETTINGS['quadrature_points'],SETTINGS['retained_index'],SETTINGS['max_derivative'])
            if exact_upper(alias)!=data['alias_error_upper'] or [exact_upper(v) for v in tails]!=data['tail_upper']:
                raise ValueError('refined alias/tail replay differs')
            coefficients=[C.unpack_ball(row) for row in arrays[key+'/coefficients']]
            K=SETTINGS['retained_index']
            if len(coefficients)!=2*K+1 or any(not v.is_finite() for v in coefficients):
                raise ValueError('refined coefficient inventory is incomplete')
            if not coefficients[K].imag.contains(0) or any(not coefficients[K+k].overlaps(coefficients[K-k].conjugate()) for k in range(1,K+1)):
                raise ValueError('refined real-profile conjugacy failed')
            old=[C.unpack_ball(row) for row in arrays[key+'/previous_coefficients']]
            oldK=(len(old)-1)//2
            if any(not coefficients[K+k].overlaps(old[oldK+k]) for k in range(-oldK,oldK+1)):
                raise ValueError('the two certified profile intervals disagree')
            rows[key]={'retained_coefficients':len(coefficients),'derivative_cells':data['derivative_cells'],
                'tail_sup_orders_0_1_2':[float(v) for v in tails],
                'alias_error_upper':float(alias),'previous_coefficients_compatible':True}
    return rows


def record(arrays,meta,payload):
    return {'schema':'NSC-KS-PROFILE-FOURIER-REFINEMENT-v1','accountable_author':'Douglas Ek',
        'status':'PASS: tighter certified axial-profile enclosure; physical local gate OPEN',
        'settings':SETTINGS,'measurements':measure(arrays,meta),'runtime':meta['runtime'],
        'scope':{'source_energy_UV_certificate':False,'continuous_state_error_bound':None,
                 'physical_local_gate':'OPEN','field_runs':0,'source_runs':0,
                 'metric_timestep':False,'geometric_profile_only':True},
        'source_hashes':{p:meta['signature'][p] for p in OWNED},
        'input_hashes':{p:meta['signature'][p] for p in INPUTS},'reused_artifacts':meta['reused_artifacts'],
        'payload':payload,'reproducer':'python scripts/derive_nsc_ks_profile_fourier_refinement.py --check'}


def run():
    if (ROOT/OUTPUT).exists() or (ROOT/PAYLOAD).exists():raise FileExistsError('refinement exists; use --check')
    sig=signatures();previous=C.check()
    with np.load(ROOT/previous['payload']['path'],allow_pickle=False) as f:old={k:f[k] for k in f.files}
    oldmeta=json.loads(old['metadata_json'].tobytes())
    trajectory=oldmeta['reused_artifact']
    if digest(trajectory['path'])!=trajectory['sha256']:raise ValueError('trajectory input changed')
    with np.load(ROOT/trajectory['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
    prepared,family,_=C.C.unpack(a,json.loads(a['metadata_json'].tobytes()))
    cpu,wall=time.process_time(),time.monotonic()
    def stop(*_):raise RuntimeError('geometric Fourier refinement exceeded120 CPU seconds')
    prior=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    try:
        result=enclose_profile_fourier(AnalyticRadiusFamily(family),[(n,0) for n in range(1,5)],
            oldmeta['period_origin'],oldmeta['period_length'],**SETTINGS)
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,prior)
    arrays={};profiles={}
    with ctx.workprec(BITS):
        for (p,q),row in result.items():
            key=f'profile_{p}_{q}'
            arrays[key+'/coefficients']=np.asarray([C.pack_ball(v) for v in row.coefficients])
            arrays[key+'/previous_coefficients']=old[key+'/coefficients']
            profiles[key]={'powers':[p,q],'derivative_l1_upper':exact_upper(row.derivative_l1_upper),
                'alias_error_upper':exact_upper(row.alias_error_upper),
                'tail_upper':[exact_upper(v) for v in row.uniform_tail_bounds],
                'derivative_cells':row.derivative_cells}
    meta={'signature':sig,'profiles':profiles,'period_origin':oldmeta['period_origin'],'period_length':oldmeta['period_length'],
          'reused_artifacts':[previous['payload'],trajectory],'prepared_source_digest':prepared.fixed_preparation_digest,
          'python_flint_version':flint.__version__,
          'runtime':{'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,'CPU_cap':CAP}}
    if signatures()!=sig:raise ValueError('profile refinement source/input changed')
    arrays['metadata_json']=np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);payload={'path':PAYLOAD,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)}
    output=record(arrays,meta,payload)
    (ROOT/PAYLOAD).write_bytes(raw);(ROOT/OUTPUT).write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return output


def check():
    saved=json.loads((ROOT/OUTPUT).read_text());C.check()
    if digest(PAYLOAD)!=saved['payload']['sha256']:raise ValueError('refined Fourier payload changed')
    with np.load(ROOT/PAYLOAD,allow_pickle=False) as f:a={k:f[k] for k in f.files}
    meta=json.loads(a['metadata_json'].tobytes())
    if meta['signature']!=signatures():raise ValueError('refined Fourier source/input changed')
    if record(a,meta,saved['payload'])!=saved:raise ValueError('refined Fourier replay differs')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    result=run() if p.parse_args().run else check()
    print(json.dumps({k:result[k] for k in ('status','measurements','runtime')},indent=2))
