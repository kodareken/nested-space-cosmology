#!/usr/bin/env python3
"""One continuous residual-cell certificate from saved fields; no evolution."""
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
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_ks_trajectory_control as T
import derive_nsc_ks_profile_fourier_bound as P
import derive_nsc_ks_profile_fourier_refinement as Q
from recursive_horizons.nsc_ks_ball_trajectory import BallFourierSegment,exact_upper,restored_upper
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily,time_operator_enclosure
from recursive_horizons.nsc_ks_residual_polynomial import ResidualPolynomial
from recursive_horizons.nsc_ks_fourier_residual_bound import (
    polynomial_residual_bounds,operator_remainder_bounds,profile_sup_bounds,
)

OUTPUT=ROOT/'results/development/nsc-ks-continuous-residual-cell.json'
OWNED=('scripts/derive_nsc_ks_continuous_residual_cell.py','docs/nsc-ks-continuous-residual-cell.md',
       'src/recursive_horizons/nsc_ks_fourier_residual_bound.py',
       'src/recursive_horizons/nsc_ks_residual_polynomial.py')
INPUTS=tuple(dict.fromkeys((*Q.OWNED,*Q.INPUTS,Q.OUTPUT,'results/development/nsc-ks-radius-bound.json')))
SEGMENT,BITS,CAP=30,120,90.


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}


def compute():
    sig=signatures()
    trajectory=json.loads((ROOT/T.OUTPUT).read_text());fourier=Q.check()
    for source in (trajectory,fourier):
        if digest(source['payload']['path'])!=source['payload']['sha256']:raise ValueError('proof input payload changed')
    with np.load(ROOT/trajectory['payload']['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
    m=json.loads(a['metadata_json'].tobytes());prepared,family,segments=T.unpack(a,m)
    with np.load(ROOT/fourier['payload']['path'],allow_pickle=False) as f:b={k:f[k] for k in f.files}
    bm=json.loads(b['metadata_json'].tobytes())
    if bm['prepared_source_digest']!=prepared.fixed_preparation_digest:raise ValueError('fixed source binding changed')
    mesh=prepared.binding.computational_z
    if bm['period_origin']!=float(mesh[0]) or bm['period_length']!=float((mesh[1]-mesh[0])*len(mesh)):
        raise ValueError('Fourier origins or numerical periods differ')
    prepared.require_history(family)
    cpu,wall=time.process_time(),time.monotonic()
    def stop(*_):raise RuntimeError('one continuous residual certificate exceeded90 CPU seconds')
    prior=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    try:
        with ctx.workprec(BITS):
            profiles={tuple(row['powers']):{'coefficients':[P.unpack_ball(v) for v in b[key+'/coefficients']],
                'tail':[restored_upper(v) for v in row['tail_upper']]} for key,row in bm['profiles'].items()}
            field=BallFourierSegment(segments[SEGMENT],prepared.column_weights,len(mesh),bm['period_length'],bits=BITS)
            model=AnalyticRadiusFamily(family)
            polys,errors=time_operator_enclosure(model,field.rho_start,field.rho_end,prepared.angular,degree=8)
            terms={k:v for k,v in polys.items() if isinstance(k,tuple)}
            residual=ResidualPolynomial(field,polys['inv_a2'],polys['inv_a'],polys['inv_ar'],terms,
                prepared.mass,prepared.angular,prepared.source_energies)
            polynomial=polynomial_residual_bounds(residual,profiles)
            radius=json.loads((ROOT/'results/development/nsc-ks-radius-bound.json').read_text())
            if radius['control']['prepared_source_digest']!=prepared.fixed_preparation_digest:
                raise ValueError('radius remainder uses another preparation')
            tails=[arb(v['exact_rational']) for v in radius['bounds']['potential_remainder_sup_orders_0_1_2']]
            remainder=operator_remainder_bounds(field,errors,profile_sup_bounds(profiles,field.length),tails,
                prepared.mass,prepared.angular,prepared.source_energies)
            total=[(x+y).upper() for x,y in zip(polynomial,remainder)]
            bounds={'polynomial_and_profile_tail':[exact_upper(v) for v in polynomial],
                    'time_coefficient_and_radius_remainder':[exact_upper(v) for v in remainder],
                    'total_continuous_normalized_residual':[exact_upper(v) for v in total]}
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,prior)
    if sig!=signatures():raise ValueError('residual certificate source/input changed')
    return {'schema':'NSC-KS-CONTINUOUS-RESIDUAL-CELL-v1','accountable_author':'Douglas Ek',
        'status':'PASS: one continuous residual cell; full-history field bound and physical local gate OPEN',
        'bounds':bounds,'coverage':{'segment_index':SEGMENT,'total_segments':len(segments),
             'rho_interval':[segments[SEGMENT].rho_end,segments[SEGMENT].rho_start],
             'spatial_period_origin':bm['period_origin'],'spatial_period_length':bm['period_length'],
             'whole_spatial_period':True,'whole_time_cell':True,'all_history_segments':False},
        'bindings':{'prepared_source_digest':prepared.fixed_preparation_digest,'history_fingerprint':prepared.binding.fingerprint,
                    'history_amplitude':T.C.R.C.P.ALPHA,'source_family':'14_1'},
        'scope':{'physical_local_gate':'OPEN','continuous_full_history_field_error':None,
                 'source_error_bound':None,'source_energy_UV_bound':None,'field_runs':0,
                 'source_runs':0,'metric_timestep':False,'sample_max_used_as_bound':False},
        'runtime':{'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,'CPU_cap':CAP},
        'source_hashes':{p:sig[p] for p in OWNED},'input_hashes':{p:sig[p] for p in INPUTS},
        'reused_artifacts':[trajectory['payload'],fourier['payload']],
        'reproducer':'python scripts/derive_nsc_ks_continuous_residual_cell.py --check'}


def check():
    saved=json.loads(OUTPUT.read_text())
    sig=signatures()
    if any(sig[p]!=v for p,v in {**saved['source_hashes'],**saved['input_hashes']}.items()):
        raise ValueError('continuous residual owner/input binding changed')
    for item in saved['reused_artifacts']:
        if digest(item['path'])!=item['sha256']:raise ValueError('continuous residual data changed')
    with ctx.workprec(BITS):
        for a,b,t in zip(saved['bounds']['polynomial_and_profile_tail'],
                         saved['bounds']['time_coefficient_and_radius_remainder'],
                         saved['bounds']['total_continuous_normalized_residual']):
            if not restored_upper(t)>=(restored_upper(a)+restored_upper(b)).upper():
                raise ValueError('continuous residual contributions are not enclosed')
    if saved['coverage']['all_history_segments'] or saved['scope']['continuous_full_history_field_error'] is not None:
        raise ValueError('a single residual cell cannot certify the entire field evolution')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true');g.add_argument('--recompute',action='store_true')
    args=p.parse_args()
    if args.run:
        if OUTPUT.exists():raise FileExistsError('continuous residual cell exists; use --check or --recompute')
        result=compute();OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif args.recompute:
        saved=check();result=compute()
        if result['bounds']!=saved['bounds']:raise ValueError('continuous residual recomputation differs')
    else:result=check()
    print(json.dumps({'status':result['status'],'coverage':result['coverage'],
        'continuous_residual_upper_display':[float(restored_upper(v)) for v in result['bounds']['total_continuous_normalized_residual']],
        'runtime':result['runtime']},indent=2))
