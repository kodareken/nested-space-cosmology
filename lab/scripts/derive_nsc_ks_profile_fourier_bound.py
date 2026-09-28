#!/usr/bin/env python3
"""One certified axial-profile preparation; no source or history evolution."""
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
import derive_nsc_ks_trajectory_control as C
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily
from recursive_horizons.nsc_ks_profile_fourier_bound import enclose_profile_fourier,alias_and_tail_bounds
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-ks-profile-fourier-bound.json'
PAYLOAD='results/development/artifacts/nsc-ks-profile-fourier-bound.npz'
OWNED=('scripts/derive_nsc_ks_profile_fourier_bound.py','src/recursive_horizons/nsc_ks_profile_fourier_bound.py',
       'src/recursive_horizons/nsc_ks_ball_geometry.py','src/recursive_horizons/nsc_ks_ball_operator.py',
       'src/recursive_horizons/nsc_ks_ball_trajectory.py','requirements-validation.txt',
       'docs/nsc-ks-profile-fourier-bound.md','docs/nsc-ks-ball-geometry.md')
INPUTS=tuple(dict.fromkeys((*C.OWNED,*C.INPUTS,C.OUTPUT)))
CAP,BITS=120.,120
SETTINGS={'derivative_order':8,'transition_panels':256,'quadrature_points':16384,
          'retained_index':2048,'max_derivative':2,'bits':BITS}


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}


def pack_ball(value):
    result=[]
    for part in (value.real,value.imag):
        mid,rad=part.mid().man_exp(),part.rad().man_exp()
        result.append([str(mid[0]),str(mid[1]),str(rad[0]),str(rad[1])])
    return result


def unpack_ball(row):
    return acb(*(arb((int(c[0]),int(c[1])),(int(c[2]),int(c[3]))) for c in row))


def measure(arrays,meta):
    rows={}
    with ctx.workprec(BITS):
        for key,data in meta['profiles'].items():
            B=restored_upper(data['derivative_l1_upper'])
            alias,tail=alias_and_tail_bounds(B,arb(meta['period_length']),SETTINGS['derivative_order'],
                SETTINGS['quadrature_points'],SETTINGS['retained_index'],SETTINGS['max_derivative'])
            if exact_upper(alias)!=data['alias_error_upper'] or [exact_upper(v) for v in tail]!=data['tail_upper']:
                raise ValueError('Fourier alias/tail inequality replay differs')
            values=[unpack_ball(row) for row in arrays[key+'/coefficients']]
            if len(values)!=2*SETTINGS['retained_index']+1 or any(not v.is_finite() for v in values):
                raise ValueError('finite retained Fourier coefficient inventory required')
            K=SETTINGS['retained_index']
            if not values[K].imag.contains(0) or any(not values[K+k].overlaps(values[K-k].conjugate()) for k in range(1,K+1)):
                raise ValueError('real profile Fourier conjugacy lost')
            rows[key]={'derivative_l1_upper':float(B),'retained_coefficients':len(values),
                       'alias_error_upper':float(alias),
                       'tail_sup_orders_0_1_2':[float(v) for v in tail],
                       'derivative_cells':data['derivative_cells']}
    return rows


def record(arrays,meta,payload):
    return {'schema':'NSC-KS-PROFILE-FOURIER-BOUND-v1','accountable_author':'Douglas Ek',
        'status':'PASS: certified axial-profile Fourier enclosure; physical local gate OPEN',
        'settings':SETTINGS,'measurements':measure(arrays,meta),'runtime':meta['runtime'],
        'scope':{'history_amplitude':C.C.R.C.P.ALPHA,'local_interval':'S(1)+[.12,.18]',
                  'source_energy_UV_certificate':False,'continuous_state_error_bound':None,
                  'physical_local_gate':'OPEN','field_runs':0,'source_runs':0,'metric_timestep':False,
                  'geometric_profile_only':True,'empirical_Fourier_decay_fit':False},
        'source_hashes':{p:meta['signature'][p] for p in OWNED},
        'input_hashes':{p:meta['signature'][p] for p in INPUTS},'reused_artifact':meta['reused_artifact'],
        'payload':payload,'reproducer':'python scripts/derive_nsc_ks_profile_fourier_bound.py --check'}


def run():
    if (ROOT/OUTPUT).exists() or (ROOT/PAYLOAD).exists():raise FileExistsError('profile bound exists; use --check')
    sig=signatures();old=json.loads((ROOT/C.OUTPUT).read_text())
    if digest(old['payload']['path'])!=old['payload']['sha256']:raise ValueError('trajectory input changed')
    with np.load(ROOT/old['payload']['path'],allow_pickle=False) as loaded:a={k:loaded[k] for k in loaded.files}
    oldmeta=json.loads(a['metadata_json'].tobytes());prepared,family,_=C.unpack(a,oldmeta)
    model=AnalyticRadiusFamily(family);mesh=prepared.binding.computational_z
    origin,length=float(mesh[0]),float((mesh[1]-mesh[0])*len(mesh))
    cpu,wall=time.process_time(),time.monotonic()
    def stop(*_):raise RuntimeError('axial profile enclosure exceeded120 CPU seconds')
    prior=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    try:
        result=enclose_profile_fourier(model,[(n,0) for n in range(1,5)],origin,length,**SETTINGS)
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,prior)
    arrays={};profiles={}
    with ctx.workprec(BITS):
        for (p,q),row in result.items():
            key=f'profile_{p}_{q}'
            arrays[key+'/coefficients']=np.asarray([pack_ball(v) for v in row.coefficients])
            profiles[key]={'powers':[p,q],'derivative_l1_upper':exact_upper(row.derivative_l1_upper),
                'alias_error_upper':exact_upper(row.alias_error_upper),
                'tail_upper':[exact_upper(v) for v in row.uniform_tail_bounds],
                'derivative_cells':row.derivative_cells}
    meta={'signature':sig,'profiles':profiles,'period_origin':origin,'period_length':length,
          'reused_artifact':old['payload'],'prepared_source_digest':prepared.fixed_preparation_digest,
          'python_flint_version':flint.__version__,
          'runtime':{'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,'CPU_cap':CAP}}
    if signatures()!=sig:raise ValueError('profile enclosure source/input changed')
    arrays['metadata_json']=np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);payload={'path':PAYLOAD,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)}
    output=record(arrays,meta,payload)
    (ROOT/PAYLOAD).write_bytes(raw);(ROOT/OUTPUT).write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return output


def check():
    saved=json.loads((ROOT/OUTPUT).read_text())
    if digest(PAYLOAD)!=saved['payload']['sha256']:raise ValueError('profile Fourier payload changed')
    with np.load(ROOT/PAYLOAD,allow_pickle=False) as loaded:a={k:loaded[k] for k in loaded.files}
    meta=json.loads(a['metadata_json'].tobytes())
    if meta['signature']!=signatures():raise ValueError('profile Fourier source/input changed')
    if record(a,meta,saved['payload'])!=saved:raise ValueError('profile Fourier record replay differs')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    result=run() if p.parse_args().run else check()
    print(json.dumps({k:result[k] for k in ('status','measurements','runtime')},indent=2))
