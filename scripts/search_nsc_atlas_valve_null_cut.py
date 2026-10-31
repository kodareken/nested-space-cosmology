"""Search genuine time-tilted closed cuts of one inherited null sheet.

The optimized parameters describe the tested surface, never the parent
action, source or evolution. A sampled candidate is not a child admission.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_null_parent import GHNullHistorySampler
from recursive_horizons.nsc_atlas_null_cuts import evaluate_spacetime_cut
from recursive_horizons.nsc_atlas_null_transport import DenseNullBundle
from recursive_horizons.nsc_atlas_trapped_surfaces import _angular_basis

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sheet',type=Path,required=True)
    parser.add_argument('--parent-record',type=Path)
    parser.add_argument('--time-degree',type=int,default=2)
    parser.add_argument('--time-parameterization',choices=('harmonic','logit'),default='harmonic')
    parser.add_argument('--cut-degree',type=int,default=8)
    parser.add_argument('--initial-coordinate-time',type=float,default=.3018)
    parser.add_argument('--warm-start',type=Path)
    parser.add_argument('--budget-seconds',type=float,default=300.)
    parser.add_argument('--maxiter',type=int,default=30)
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    if not 0<=args.time_degree<=6 or args.cut_degree<args.time_degree:
        raise ValueError('bounded time degree and adequate cut reconstruction required')
    source_names=['scripts/search_nsc_atlas_valve_null_cut.py','src/recursive_horizons/nsc_atlas_null_cuts.py',
        'src/recursive_horizons/nsc_atlas_null_parent.py','src/recursive_horizons/nsc_atlas_null_transport.py']
    before={name:sha(ROOT/name) for name in source_names}
    commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    source_bound=all(hashlib.sha256(subprocess.check_output(['git','show',commit+':'+name])).hexdigest()==value
        for name,value in before.items())
    started=time.perf_counter()
    sheet_path=(ROOT/args.sheet).resolve();sheet=json.loads(sheet_path.read_bytes())
    ray=DenseNullBundle.from_sheet_record(sheet_path)
    parent_path=Path(sheet['parent_history']['record_path']) if args.parent_record is None else ROOT/args.parent_record
    if sha(parent_path)!=sheet['parent_history']['record_sha256']:raise ValueError('parent record identity differs from the sheet')
    sampler=GHNullHistorySampler.from_files(parent_path)
    order=tuple(ray.order)
    if args.cut_degree>=order[0] or 2*args.cut_degree+1>order[1]:
        raise ValueError('stored ray directions do not resolve requested cut degree')
    angular=_angular_basis(args.time_degree,order);Y=angular['values'];r=sampler.r
    bound=ray.stored_polynomial_noncaustic_bound()
    if not bound['all_sampled_rays_certified']:
        raise ValueError('stored ray polynomial noncaustic bound is unresolved; diagnose before searching cuts')
    buffer=1e-6*(ray.end-ray.start)
    tlo=ray.start+buffer;thi=ray.end-buffer
    if not tlo<args.initial_coordinate_time<thi:raise ValueError('initial numerical cut must lie within the recorded sheet')
    fraction=(args.initial_coordinate_time-tlo)/(thi-tlo)
    coeff=np.zeros((args.time_degree+1)**2)
    coeff[0]=math.sqrt(4*math.pi)*(args.initial_coordinate_time if args.time_parameterization=='harmonic'
        else math.log(fraction/(1-fraction)))
    warm_binding=None
    if args.warm_start:
        warm_path=(ROOT/args.warm_start).resolve();warm=json.loads(warm_path.read_bytes())
        if warm['sheet_sha256']!=sha(sheet_path) or warm['numerical_time_degree']!=args.time_degree or warm['time_parameterization']!=args.time_parameterization:
            raise ValueError('warm cut belongs to a different sheet or time family')
        coeff=np.asarray(warm['best_cut']['time_coefficients'],dtype=float)
        if args.time_parameterization=='harmonic':
            nodes=Y@coeff
            shift=max(0.,tlo+buffer-float(nodes.min()))+min(0.,thi-buffer-float(nodes.max()))
            coeff=coeff.copy();coeff[0]+=math.sqrt(4*math.pi)*shift
        warm_binding=dict(path=warm_path.relative_to(ROOT).as_posix(),sha256=sha(warm_path))
    cache={};best=None;rejected=[];evaluations=0

    def evaluate(coefficients):
        nonlocal best,evaluations
        if time.perf_counter()-started>args.budget_seconds:raise TimeoutError('time-tilted cut search budget exhausted')
        key=np.asarray(coefficients,dtype=float).tobytes()
        if key in cache:return cache[key]
        times=(Y@coefficients if args.time_parameterization=='harmonic' else tlo+(thi-tlo)*expit(Y@coefficients))
        # Search constraints use an interior buffer. Finite-difference
        # trials may leave that search buffer while still lying in the
        # ACTUAL recorded domain; those are valid geometry evaluations.
        if np.min(times)<ray.start or np.max(times)>ray.end:
            theta=np.full(2*len(times),1e3)
            row=dict(rejected_reason='time profile leaves the recorded ray domain',time_coefficients=list(map(float,coefficients)))
            if len(rejected)<12:rejected.append(row)
            evaluations+=1;cache[key]=(theta,row)
            return theta,row
        state=ray.events(times)
        try:
            cut=evaluate_spacetime_cut(times,state['X'],
                lambda t,X:sampler.event_slice(t,X,chunk=256),frame_scale=r,
                angular_degree=args.cut_degree,order=order,outgoing_k=state['k'],
                ray_area_determinant=np.linalg.det(state['D']),expansion_tolerance=1e-6*r)
            theta=np.r_[cut.theta_plus/r,cut.theta_minus/r]
            row=cut.summary();row.update(time_coefficients=list(map(float,coefficients)),
                maximum_expansion_over_H_reference=float(theta.max()),
                actual_time_node_range=[float(times.min()),float(times.max())],
                physical_child_admitted=False)
            affine=np.trace(np.linalg.solve(state['D'].swapaxes(1,2),state['P'].swapaxes(1,2)).swapaxes(1,2),axis1=1,axis2=2)
            gap=cut.theta_plus-affine/cut.normal_energy
            row['Sachs_geometric_theta_RMS_gap_over_H']=float(np.sqrt(np.mean(gap**2))/r)
            row['Sachs_geometric_theta_max_gap_over_H']=float(np.max(abs(gap))/r)
            row['integrated_Ricci_source_gap_range']=[float(np.min(state['integrated_Rkk']-state['integrated_source_Rkk'])),
                float(np.max(state['integrated_Rkk']-state['integrated_source_Rkk']))]
            if best is None or row['maximum_expansion_over_H_reference']<best['maximum_expansion_over_H_reference']:
                best=row
        except ValueError as error:
            theta=np.full(2*len(times),1e3)
            row=dict(rejected_reason=str(error),time_coefficients=list(map(float,coefficients)))
            if len(rejected)<12:rejected.append(row)
        evaluations+=1;cache[key]=(theta,row)
        return theta,row

    status='completed_finite_search';reason=None;optimum=None
    with threadpool_limits(limits=args.threads):
        first,initial=evaluate(coeff)
        variables=np.r_[coeff,float(first.max())+.01]
        limits=([(tlo*math.sqrt(4*math.pi),thi*math.sqrt(4*math.pi))]+[(-2*(thi-tlo),2*(thi-tlo))]*(len(coeff)-1)+[(-1000.,1000.)]
            if args.time_parameterization=='harmonic' else [(-30.,30.)]+[(-12.,12.)]*(len(coeff)-1)+[(-1000.,1000.)])
        def constraints(x):return x[-1]-evaluate(x[:-1])[0]
        conditions=[{'type':'ineq','fun':constraints}]
        if args.time_parameterization=='harmonic':
            conditions.extend([{'type':'ineq','fun':lambda x:Y@x[:-1]-tlo},
                               {'type':'ineq','fun':lambda x:thi-Y@x[:-1]}])
        try:
            optimum=minimize(lambda x:x[-1],variables,jac=lambda x:np.r_[np.zeros(len(coeff)),1.],
                method='SLSQP',bounds=limits,constraints=conditions,
                options={'maxiter':args.maxiter,'ftol':1e-8,'disp':False})
            evaluate(optimum.x[:-1])
        except TimeoutError as error:status='budget_exhausted';reason=str(error)
        validations=[]
        if best is not None:
            for degree in sorted(set([max(args.time_degree,args.cut_degree-2),args.cut_degree])):
                times=(Y@np.asarray(best['time_coefficients']) if args.time_parameterization=='harmonic'
                    else tlo+(thi-tlo)*expit(Y@np.asarray(best['time_coefficients'])))
                state=ray.events(times)
                try:
                    cut=evaluate_spacetime_cut(times,state['X'],lambda t,X:sampler.event_slice(t,X,chunk=256),
                        frame_scale=r,angular_degree=degree,order=order,outgoing_k=state['k'],
                        ray_area_determinant=np.linalg.det(state['D']),expansion_tolerance=1e-6*r)
                    validations.append(cut.summary())
                except ValueError as error:validations.append(dict(angular_degree=degree,rejected_reason=str(error)))
    result=dict(schema='nsc-valve-time-tilted-cut-v1',status=status,reason=reason,
        sheet_path=sheet_path.relative_to(ROOT).as_posix(),sheet_sha256=sha(sheet_path),
        parent_history=sampler.provenance,ray_dense_cache=sheet['ray_dense_cache'],
        numerical_time_degree=args.time_degree,cut_degree=args.cut_degree,angular_order=list(order),
        time_parameterization=args.time_parameterization,
        warm_start=warm_binding,search_time_buffer=buffer,
        source_sha256=before,source_bound_to_commit=source_bound,
        source_unchanged_during_run=before=={name:sha(ROOT/name) for name in source_names},
        producer_commit=commit,
        noncaustic_bound=bound,initial_cut=initial,best_cut=best,rejected_examples=rejected,
        alternate_reconstruction=validations,evaluations=evaluations,elapsed_seconds=time.perf_counter()-started,
        optimizer=None if optimum is None else dict(success=bool(optimum.success),message=str(optimum.message),
            iterations=int(optimum.nit),reported_maximum_expansion=float(optimum.fun)),
        physical_child_admitted=False,
        scope='Numerical surface/time parameters vary on one fixed inherited null sheet. No action/state retuning; finite search misses prove no universal nonexistence. Sampled negative expansions would still require independent metric/angular validation and receiving-region construction.')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print('Time-tilted cut:',status,'best max expansion/H=',None if best is None else best['maximum_expansion_over_H_reference'],
          'evaluations=',evaluations,'seconds=',result['elapsed_seconds'],flush=True)


if __name__=='__main__':main()
