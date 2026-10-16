#!/usr/bin/env python3
"""Authenticate original source interpolation before a common pulse update."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate,producer,plain
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_trace import SourceSnapshot
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('scripts/check_nsc_atlas_quiet_source_front.py',
         'scripts/diagnose_nsc_atlas_boundary_rows.py',
         'src/recursive_horizons/nsc_atlas_characteristic_trace.py',
         'src/recursive_horizons/nsc_atlas_reference_wall_gauge.py')


def polynomial(x0,x1,y0,d0,y1,d1):
    """Exact binary-input Hermite coefficients in t=(x-x0)/(x1-x0)."""
    w=Q(x1)-Q(x0);a,b,c,d=map(Q,(y0,d0,y1,d1))
    return (a,w*b,3*(c-a)-w*(2*b+d),-2*(c-a)+w*(b+d))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    budget=Budget(maximum_evaluations=37577,cpu_seconds=20.)
    _,_,_,identity=authenticate(ROOT,budget)
    path=producer.BASE_INPUT;raw=(ROOT/path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=identity['source_bindings'][path]:
        raise ValueError('original source payload differs from its authenticated binding')
    base=json.loads(raw);geometry=base['geometry'];data=base['initial_data']['finite']
    rows={k:np.asarray(v,dtype=float)[::2] for k,v in data['exact_rows'].items()}
    x=rows['x'];front=geometry['support_x'][0]
    j=int(np.searchsorted(x,front,side='right')-1)
    if len(x)!=129 or not x[j]<front<x[j+1]:
        raise ValueError('expected original 129 knots and an interior support front')
    polynomials=[]
    for i in range(j+1):
        polynomials.append(polynomial(float(x[i]),float(x[i+1]),
            float(rows['h'][i]),float(rows['h_x'][i]),float(rows['h'][i+1]),float(rows['h_x'][i+1])))
        budget.check(256)
    nonzero=[i for i,p in enumerate(polynomials) if any(p)]
    snapshot=SourceSnapshot(geometry,alpha=data['alpha'],x_nodes=x,source_rows=rows,budget=budget)
    query=np.array([x[j],x[j]+(front-x[j])/2,front])
    values=snapshot.sample(query)
    if np.any(values['c']!=0.) or np.any(values['c_x']!=0.):
        raise ValueError('selected scalar quiet queries are not exactly zero')
    p=polynomials[-1];w=Q(float(x[j+1]))-Q(float(x[j]));defects=[]
    # Compare cancellation against the coefficient scale, not a nearly zero
    # result. This is a floating evaluation control, not an interval enclosure.
    scales=[sum(abs(v) for v in p),sum(k*abs(p[k]) for k in range(1,4))/w,
            sum(k*(k-1)*abs(p[k]) for k in range(2,4))/(w*w)]
    tolerances={key:128*np.finfo(float).eps*float(scale) for key,scale in
                zip(('h','h_x','h_xx'),scales)}
    for n,q in enumerate(query):
        t=(Q(float(q))-Q(float(x[j])))/w
        exact=[sum(p[k]*t**k for k in range(4)),
               sum(k*p[k]*t**(k-1) for k in range(1,4))/w,
               sum(k*(k-1)*p[k]*t**(k-2) for k in range(2,4))/(w*w)]
        defects.append({key:float(values[key][n])-float(v)
                        for key,v in zip(('h','h_x','h_xx'),exact)})
        budget.check(1024)
        for key in tolerances:
            if abs(defects[-1][key])>tolerances[key]:
                raise ValueError('Hermite evaluation differs from exact binary polynomial: '+str(
                    dict(query=n,field=key,defect=defects[-1][key],tolerance=tolerances[key],
                         charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started)))
    gauge=ReferenceWallGauge(geometry,order=64,budget=budget)
    wall=gauge.sample(query[[0,2]])
    fields=('x','h','h_x','h_xx','r','r_x','r_xx','c','c_x',
            'source_constraint','source_h_xx_defect')
    result=dict(schema='nsc-atlas-quiet-source-front-v1',alpha=data['alpha'],
        knots=129,front=front,crossing_interval=dict(index=j,left=float(x[j]),right=float(x[j+1])),
        exact_radial_polynomial_audit=dict(intervals=j+1,nonzero_interval_indices=nonzero,
            crossing_coefficients=[str(v) for v in p],
            zero_on_entire_pre_front_interval=not nonzero),
        source_samples={k:values[k] for k in fields},
        exact_polynomial_evaluation_defects=defects,
        exact_polynomial_evaluation_tolerances=tolerances,
        reference_wall={k:wall[k] for k in ('x','r','T','s','W','W_x','T_x','r_x_total')},
        input=dict(path=path,sha256=hashlib.sha256(raw).hexdigest(),
            authenticated_input_files=identity['authenticated_input_files'],
            operator_sha256=identity['sha256'],manifest_sha256=identity['manifest_sha256']),
        source_bindings={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
        scope='same retained source; exact scalar quietness and binary radial polynomial check; no source-error enclosure or trajectory',
        execution=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
            source_queries=3,gauge_queries=2,source_integrations=0,time_steps=0),
        prior_attempt=dict(outcome='stopped at exact-polynomial comparison using result-relative tolerance',
            charged_values_upper_bound=12423,cpu_seconds_upper_bound=20.,
            actual_cpu_unrecorded=True,source_integrations=0,time_steps=0),
        limits=dict(charged_values=37577,cpu_seconds=20.,
            cumulative_charged_values_upper_bound=50000,cumulative_cpu_seconds_upper_bound=40.))
    budget.check()
    output=json.dumps(plain(result),indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.output:
        args.output.write_text(output,encoding='utf-8',newline='\n')
    else: print(output,end='')


if __name__=='__main__': main()
