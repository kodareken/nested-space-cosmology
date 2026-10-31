"""Bounded passive onset roots from three strictly authenticated GH histories."""
from concurrent.futures import ProcessPoolExecutor
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_trapped_parent import GHSurfaceGeometry

ROOT=Path(__file__).resolve().parents[1]
BUDGET=180.


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculation(arguments):
    degree,deadline=arguments
    started=time.perf_counter()
    def guard():
        if time.time()>deadline:raise RuntimeError('180-second onset-analysis wall budget exhausted')
    guard()
    recordpath=ROOT/f'.build/valve-e1-L{degree}-extended-v1.json'
    historypath=ROOT/f'.build/valve-e1-L{degree}-extended-v1-history.npz'
    record=json.loads(recordpath.read_text(encoding='utf-8'))
    with threadpool_limits(limits=2):
        initial=GHSurfaceGeometry.from_history(historypath,record,coordinate_t=0.,retain_history=True)
        history=initial.retained_history;system=initial.system;dof=system.dof
        t0=float(history.bounds[0][0]);proper0=float(initial.state[dof+8])
        cache={};queries=0
        def reading(t):
            nonlocal queries
            guard();t=float(t)
            if t in cache:return cache[t]
            state=history.evaluate(t)
            provider=GHSurfaceGeometry(system,state,coordinate_t=t,provenance=initial.provenance)
            raw=state[dof:dof+4];norm=float(np.linalg.norm(raw));X=raw/norm
            data=provider(X[None,:])
            Pi=data['Pi'][0];eigen=data['K_eigenvalues'][0];W=float(data['W'][0])
            F=float((Pi@Pi+eigen@eigen/system.parent.kappa)/W)
            theta=-float(np.sum(eigen))
            s=float(system.sstar+data['q'][0,0]/math.sqrt(system.K))
            ns=float(Pi[0]/math.sqrt(system.K))
            angle=system.observer_angle(state)
            if angle is None or not -1.e-8<=angle<math.pi:
                raise ValueError('event analysis left the first positive wrapped phase')
            N=float(data['lapse'][0]);beta=data['shift'][0];h=data['h'][0]
            g=np.zeros((4,4));g[1:,1:]=h;g[0,1:]=g[1:,0]=h@beta
            g[0,0]=-N*N+beta@h@beta
            u=state[dof+4:dof+8].copy();unorm=float(u@g@u);u/=math.sqrt(-unorm)
            p=N*Pi+beta@data['dq'][0]
            properrates=u[0]*p+u[1:]@data['dq'][0]
            mu2=400*math.exp(-12*s);chi=float(data['q'][0,1])
            row=dict(coordinate_t=t,coordinate_duration=t-t0,
                observer_proper_duration=float(state[dof+8]-proper0),F=F,
                slice_normal_theta_over_H=theta/system.r,
                external_Jordan_normal_theta_over_H=math.exp(6*s)*(theta-18*ns)/system.r,
                s=s,normal_s_rate=ns,normal_s_rate_over_H=ns/system.r,
                chi=chi,normal_chi_rate=float(Pi[1]),geodesic_chi_rate=float(properrates[1]),
                geodesic_s_rate=float(properrates[0]/math.sqrt(system.K)),
                canonical_normal_momenta=Pi.tolist(),physical_K_eigenvalues=eigen.tolist(),
                potential_W=W,oscillator_action_loading=float(.5*(Pi[1]**2+mu2*chi**2)*math.exp(6*s)),
                observer_angle=float(angle),protocol_pass_count=float(angle/(2*math.pi*system.parent.offset)),
                observer_raw_point_norm_error=abs(norm-1),geodesic_input_unit_norm_error=abs(unorm+1))
            cache[t]=row;queries+=1
            return row
        physical_end=float(history.metadata['physical_coordinate_endpoint'])
        tend=brentq(lambda t:float(history.evaluate(t)[dof+8]-proper0)-.075,t0,physical_end,
                     xtol=1.e-8,rtol=1.e-12)
        sample_times=np.linspace(t0,tend,24)
        samples=[reading(t) for t in sample_times]
        def first_event(key,threshold,direction):
            for a,b in zip(samples[:-1],samples[1:]):
                va,vb=a[key]-threshold,b[key]-threshold
                crossing=(va<0<=vb) if direction>0 else (va>0>=vb)
                if crossing:
                    t=brentq(lambda value:reading(value)[key]-threshold,a['coordinate_t'],b['coordinate_t'],
                              xtol=1.e-8,rtol=1.e-12)
                    return dict(found=True,bracket_coordinate_t=[a['coordinate_t'],b['coordinate_t']],
                                bracket_values=[a[key],b[key]],reading=reading(t))
            return dict(found=False,reason='No first-pulse bracket in the declared proper-duration interval.')
        focusing=first_event('F',1.,1)
        contraction=first_event('slice_normal_theta_over_H',0.,-1)
        guard()
    return dict(degree=degree,producer_commit=record['producer_commit'],
        producer_record=recordpath.relative_to(ROOT).as_posix(),producer_record_sha256=sha(recordpath),
        history=historypath.relative_to(ROOT).as_posix(),history_sha256=initial.provenance['cache_sha256'],
        strict_source_binding=initial.provenance['source_binding'],source_pin_count=len(record['source_sha256']),
        physical_history_endpoint=physical_end,search_proper_interval=[0.,.075],
        initial_sample_count=len(samples),point_geometry_queries=queries,history_load_count=1,
        H_reference=system.r,alpha=system.parent.offset,
        first_F_equals_one=focusing,first_normal_expansion_positive_to_negative=contraction,
        initial_reading=samples[0],final_search_reading=samples[-1],elapsed_seconds=time.perf_counter()-started)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degrees',type=int,nargs='+',default=[16,20,24])
    parser.add_argument('--budget-seconds',type=float,default=BUDGET)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    if not 0<args.budget_seconds<=600:raise ValueError('bounded positive wall allowance required')
    started=time.perf_counter();deadline=time.time()+args.budget_seconds
    with ProcessPoolExecutor(max_workers=min(3,len(args.degrees))) as executor:
        rows=list(executor.map(calculation,[(L,deadline) for L in args.degrees]))
    result=dict(schema='nsc-atlas-valve-onset-events-v1',script_sha256=sha(__file__),
        elapsed_seconds=time.perf_counter()-started,budget_seconds=args.budget_seconds,
        numerical_time_root_tolerance=1.e-8,cohorts=rows,
        scope=dict(F='Geodesic-relaunch focusing proxy at the tracked observer position; not a trapped-surface threshold. F=(sum Pi^2+sum K_eigenvalues^2/kappa)/W.',
                   theta='Actual slice-normal expansion -trK at the tracked observer position; not the full geodesic-bundle expansion.',
                   pass_count='Early positive observer_angle/(2*pi*alpha); a synchronization protocol count, not nesting depth.',
                   state='Passive evaluation of authenticated actual DOP853 dense polynomials; no PDE restart or sparse interpolation.',
                   precision='Resolution spreads are empirical indicators; no continuum enclosure or formation theorem.'))
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    for row in rows:
        print(json.dumps(dict(degree=row['degree'],F_event=row['first_F_equals_one'],theta_event=row['first_normal_expansion_positive_to_negative']),sort_keys=True),flush=True)
    print('Output SHA256',sha(output),'elapsed',result['elapsed_seconds'],flush=True)


if __name__=='__main__':main()
