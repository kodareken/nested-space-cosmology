#!/usr/bin/env python3
"""Certified axial profiles on the fine trajectory's exact numerical period."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_ks_fine_trajectory as F
import derive_nsc_ks_profile_fourier_bound as P
import derive_nsc_ks_profile_fourier_refinement as Q
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily
from recursive_horizons.nsc_ks_profile_identity import profile_identity,profile_description
from recursive_horizons.nsc_ks_profile_fourier_bound import enclose_profile_fourier,alias_and_tail_bounds
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT=ROOT/'results/development/nsc-ks-fine-profile-table.json'
PAYLOAD=ROOT/'results/development/artifacts/nsc-ks-fine-profile-table.npz'
OWNED=('scripts/derive_nsc_ks_fine_profile_table.py','docs/nsc-ks-fine-profile-table.md')
INPUTS=tuple(dict.fromkeys((*Q.OWNED,*Q.INPUTS,*F.OWNED,str(F.OUTPUT.relative_to(ROOT)),Q.OUTPUT)))
SETTINGS=Q.SETTINGS


def digest(path):
    path=Path(path);return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}


def measures(arrays,meta):
    rows={}
    with ctx.workprec(SETTINGS['bits']):
        for key,proof in meta['profiles'].items():
            alias,tails=alias_and_tail_bounds(restored_upper(proof['derivative_l1_upper']),arb(meta['period_length']),
                SETTINGS['derivative_order'],SETTINGS['quadrature_points'],SETTINGS['retained_index'],2)
            if exact_upper(alias)!=proof['alias_error_upper'] or [exact_upper(v) for v in tails]!=proof['tail_upper']:
                raise ValueError('fine-period Fourier bounds no longer replay')
            coef=[P.unpack_ball(v) for v in arrays[key+'/coefficients']];K=SETTINGS['retained_index']
            if len(coef)!=2*K+1 or any(not v.is_finite() for v in coef):raise ValueError('incomplete Fourier table')
            if not coef[K].imag.contains(0) or any(not coef[K+k].overlaps(coef[K-k].conjugate()) for k in range(1,K+1)):
                raise ValueError('real-profile conjugacy failed')
            rows[key]={'tail_sup_orders_0_1_2':[float(v) for v in tails],'alias_upper':float(alias)}
    return rows


def record(arrays,meta,payload):
    return {'schema':'NSC-KS-FINE-PROFILE-TABLE-v1','accountable_author':'Douglas Ek',
        'status':'PASS: fine-period geometric profile enclosure; field and physical gates OPEN',
        'settings':SETTINGS,'measurements':measures(arrays,meta),'runtime':meta['runtime'],
        'period_length_exact':'205/512','axial_profile_identity':meta['axial_profile_identity'],
        'source_preparation_digest':meta['source_preparation_digest'],
        'scope':{'physical_profiles_changed':False,'physical_local_gate':'OPEN','field_runs':0,'source_runs':0,
                 'source_energy_UV_certificate':False,'metric_timestep':False},
        'source_hashes':{p:meta['signature'][p] for p in OWNED},'input_hashes':{p:meta['signature'][p] for p in INPUTS},
        'reused_artifacts':meta['reused_artifacts'],'payload':payload,
        'reproducer':'python scripts/derive_nsc_ks_fine_profile_table.py --check'}


def run():
    if OUTPUT.exists() or PAYLOAD.exists():raise FileExistsError('fine profile table exists; use --check')
    sig=signatures();fine=json.loads(F.OUTPUT.read_text());old=Q.check()
    if not fine['complete_capture']:raise ValueError('fine trajectory incomplete')
    if digest(fine['payload']['path'])!=fine['payload']['sha256']:raise ValueError('fine prepared data changed')
    with np.load(ROOT/fine['payload']['path'],allow_pickle=False) as f:
        fm=json.loads(f['metadata_json'].tobytes());grid=f['computational_z']
    with np.load(ROOT/old['payload']['path'],allow_pickle=False) as f:oldmeta=json.loads(f['metadata_json'].tobytes())
    family=F.T.C.R.C.P.family(F.T.C.R.C.P.ALPHA)
    identity=profile_identity(family,include_normal_window=False)
    if identity!=fine['axial_profile_identity'] or profile_description(family)!=fm['history_description']:
        raise ValueError('analytic profile identity differs from the fine trajectory')
    length=float((grid[1]-grid[0])*len(grid));origin=float(grid[0])
    if length!=205/512:raise ValueError('fine grid period changed')
    cpu,wall=time.process_time(),time.monotonic()
    result=enclose_profile_fourier(AnalyticRadiusFamily(family),[(n,0) for n in range(1,5)],origin,length,**SETTINGS)
    arrays={};proofs={}
    with ctx.workprec(SETTINGS['bits']):
        for (p,q),item in result.items():
            key=f'profile_{p}_{q}'
            if exact_upper(item.derivative_l1_upper)!=oldmeta['profiles'][key]['derivative_l1_upper']:
                raise ValueError('physical derivative bound changed with the numerical period')
            arrays[key+'/coefficients']=np.asarray([P.pack_ball(v) for v in item.coefficients])
            proofs[key]={'powers':[p,q],'derivative_l1_upper':exact_upper(item.derivative_l1_upper),
                'alias_error_upper':exact_upper(item.alias_error_upper),
                'tail_upper':[exact_upper(v) for v in item.uniform_tail_bounds]}
    meta={'signature':sig,'profiles':proofs,'period_origin':origin,'period_length':length,
          'axial_profile_identity':identity,'profile_description':profile_description(family,include_normal_window=False),
          'source_preparation_digest':fine['source_preparation_digest'],
          'reused_artifacts':[fine['payload'],old['payload']],
          'runtime':{'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall}}
    if sig!=signatures():raise ValueError('fine profile owner/input changed')
    arrays['metadata_json']=np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);payload={'path':str(PAYLOAD.relative_to(ROOT)),
        'sha256':sha256(raw).hexdigest(),'bytes':len(raw)}
    output=record(arrays,meta,payload);PAYLOAD.write_bytes(raw)
    OUTPUT.write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return output


def check():
    saved=json.loads(OUTPUT.read_text())
    if digest(PAYLOAD)!=saved['payload']['sha256']:raise ValueError('fine profile payload changed')
    with np.load(PAYLOAD,allow_pickle=False) as f:a={k:f[k] for k in f.files}
    meta=json.loads(a['metadata_json'].tobytes())
    if meta['signature']!=signatures():raise ValueError('fine profile owner/input changed')
    if record(a,meta,saved['payload'])!=saved:raise ValueError('fine profile replay differs')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    result=run() if p.parse_args().run else check()
    print(json.dumps({k:result[k] for k in ('status','measurements','runtime')},indent=2))
