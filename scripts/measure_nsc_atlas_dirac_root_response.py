"""Bounded full-radial, finite-split vacuum comparison on the actual atlas.

The state is the declared correlated cylinder preparation. The first signed
angular sector is complete radially but is not the angular sum, a coincident
renormalized stress, or the source-induced gravitational potential ratio.
"""
from __future__ import annotations

import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy import constants
from scipy.integrate import quad_vec
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_preparation import CylinderPreparation
from recursive_horizons.nsc_atlas_dirac_retarded import (
    retarded_null_responses,retarded_response,observation_moment_residual,
    volterra_source_observation_jets,differentiated_retarded_response,endpoint_hilbert_moment)
from recursive_horizons.nsc_atlas_dirac_root_pairing import (
    NullResponse,pair_cylinder_responses,logarithmic_coefficients)

ROOT=Path(__file__).resolve().parents[1]
CPU_LIMIT=240.
SOURCES=tuple('src/recursive_horizons/'+name+'.py' for name in (
    'nsc_atlas_dirac_retarded','nsc_atlas_dirac_root_pairing','nsc_atlas_dirac_preparation',
    'nsc_atlas_dirac_mean','nsc_atlas_dirac_scattering','nsc_atlas_dirac_radar',
    'nsc_atlas_cap_wall_response','nsc_atlas_causal_reference','nsc_atlas_dirac_evolution'))+(
    'tests/test_nsc_atlas_dirac_retarded.py','tests/test_nsc_atlas_dirac_root_pairing.py',
    'tests/test_nsc_atlas_dirac_preparation.py',
    'tests/test_nsc_atlas_dirac_mean.py',
    'scripts/measure_nsc_atlas_dirac_root_response.py')


def difference(a,b):
    if a.lengths!=b.lengths:
        raise ValueError('relative root pairing requires identical null endpoints')
    return NullResponse(a.lengths,tuple(x-y for x,y in zip(a.green,b.green)),
        a.direct-b.direct,a.direct_derivative-b.direct_derivative)


def relative_pair(a,b,c,d,mu):
    # A C_root B^dagger-S C_root T^dagger
    #       = (A-S) C_root B^dagger + S C_root (B-T)^dagger.
    # Cancel common direct distributions before doing numerical integration.
    first=pair_cylinder_responses(difference(a,c),b,mu,epsabs=2e-12)
    second=pair_cylinder_responses(c,difference(b,d),mu,epsabs=2e-12)
    return (first['kernel']+second['kernel'],
        first['pv_quadrature_estimate']+second['pv_quadrature_estimate'],
        first['bessel_series_tail_bound']+second['bessel_series_tail_bound'])


def endpoint(preparation,U,V,n,integral_cache,*,order=24,workers=8,observer_jets=False,
             observation_parameter_jet=False):
    du,dv=U/(n-.5),V/(n-.5)
    u=(np.arange(n)+.5)*du;v=(np.arange(n)+.5)*dv
    fields=preparation.fields(*np.meshgrid(u,v,indexing='ij'))
    if float(preparation.fields(U,V)['preparation_weight'])!=1:
        raise ValueError('physical endpoints must lie in the matched operator region')
    cache_key=(id(preparation),U,V,observer_jets)
    integrals=[];errors=[]
    for leg,length in enumerate((V,U)) if cache_key not in integral_cache else ():
        def integrand(s):
            if observer_jets:
                f=preparation.potential_null_jets(U,s) if leg==0 else preparation.potential_null_jets(s,V)
                W=float(f['W']);derivative=float(f['Wu' if leg==0 else 'Wv'])
                second=float(f['Wuu' if leg==0 else 'Wvv'])
                return np.array([W*W,W*derivative,derivative**2+W*second])
            f=preparation.fields(U,s) if leg==0 else preparation.fields(s,V)
            W=float(f['potential'])
            derivative=float((f['potential_t']-f['potential_x'])/2 if leg==0 else
                             (f['potential_t']+f['potential_x'])/2)
            return np.array([W*W,W*derivative])
        points=[t for t in (.1*preparation.duration,.2*preparation.duration,
            .9*preparation.duration,U if leg==0 else V) if 0<t<length]
        value,error=quad_vec(integrand,0,length,points=sorted(set(points)),
            epsabs=2e-12,epsrel=2e-11,limit=160)
        integrals.append(value);errors.append(float(error))
    if cache_key not in integral_cache:
        integral_cache[cache_key]=(integrals,errors)
    integrals,errors=integral_cache[cache_key]
    result=retarded_null_responses(fields['potential'],du,dv,
        potential_u=(fields['potential_t']-fields['potential_x'])/2,
        potential_v=(fields['potential_t']+fields['potential_x'])/2,
        boundary_square_integrals=[v[0] for v in integrals],
        boundary_mixed_integrals=[v[1] for v in integrals],order=order,workers=workers,
        observation_parameter_jet=observation_parameter_jet)
    result['endpoint_quadrature_estimate']=max(errors)
    if observer_jets:
        jets=preparation.potential_null_jets(U,V)
        result['source_endpoint_observation_derivative']=volterra_source_observation_jets(
            *(float(jets[key]) for key in ('W','Wu','Wv','Wuu','Wuv','Wvv')),
            [v[0] for v in integrals],[v[1] for v in integrals],[v[2] for v in integrals])
    return result


def packed(matrix):
    return dict(real=matrix.real.tolist(),imag=matrix.imag.tolist())


def statistic_row(n,order,fraction,a,r,separation,separation_error,K,DL,DR):
    kernel,ke,kb=K;dl,de,db=DL;dr,re,rb=DR
    derivative=(dl-dr)/2
    q=float((-1j*np.trace(np.diag([1.,-1.])@derivative)).real)
    s=float(np.trace(np.array([[0.,1.],[1.,0.]])@kernel).real)
    radial=2/(4*np.pi*a*a*r*r);angular=2/(8*np.pi*a*r**3)
    pressure=radial*q-angular*s;ell=10000.
    return dict(cells=n,polynomial_order=order,separation_over_S=fraction,
        proper_separation_m=float(separation)*ell,
        proper_separation_quadrature_estimate_m=float(separation_error)*ell,
        nested_minus_cap_kernel=packed(kernel),nested_minus_cap_split_derivative=packed(derivative),
        finite_split_anisotropy_difference=pressure,
        nominal_finite_split_pressure_Pa=pressure*constants.hbar*constants.c/ell**4,
        pv_pressure_quadrature_estimate=2*radial*(de+re)/2+2*angular*ke,
        bessel_pressure_tail_bound=2*radial*(db+rb)/2+2*angular*kb)


def pair_worker(job):
    key,channel,arguments=job
    with threadpool_limits(limits=1):
        try:
            value=relative_pair(*arguments);error=None
        except Exception as failure:
            value=None;error=type(failure).__name__+': '+str(failure)
    return dict(key=key,channel=channel,value=value,error=error,pid=os.getpid(),
        process_CPU_total=time.process_time())


def polynomial_calculation(chart,preparations,a,r,rows,evidence,worker_totals,budget,*,source_comparison=False):
    S=chart.duration;U,V=.8*S,.7*S;T,x=(U+V)/2,(V-U)/2
    jobs=[];separations={};integral_cache={};identity_errors={}
    # The field phase finishes before the process phase starts. Only small
    # response polynomials cross the process boundary; field arrays do not.
    for fraction in (.02,.01):
        eta=fraction*S;maps={}
        for g,preparation in preparations.items():
            for side,point in enumerate(((U-eta/2,V+eta/2),(U+eta/2,V-eta/2))):
                budget();sample=endpoint(preparation,*point,128,integral_cache,order=48);budget()
                maps[g,side]=sample
                evidence.append(dict(cells=128,separation_over_S=fraction,geometry=g,side=side,
                    quadrature_estimate=sample['endpoint_quadrature_estimate'],
                    reconstruction_projection_residual=sample['cell_projection_residual'].tolist(),
                    response_sha256=hashlib.sha256(b''.join(c.tobytes() for key in ('cell_field','cell_spatial')
                        for c in sample[key].green)).hexdigest()))
                print(json.dumps(dict(phase='field',cells=128,separation_over_S=fraction,geometry=g,side=side)),flush=True)
        separations[fraction]=quad_vec(lambda y:float(chart.fields(T-y,T+y,region='cap')['lapse']),
            x-eta/2,x+eta/2,epsabs=1e-13,epsrel=1e-12)
        for order in (24,48):
            for source_jet in ((False,True) if source_comparison else (False,)):
                fields=[retarded_response(maps[g,side],'field',order,source_jet=source_jet)
                    for g,side in (('nested',0),('nested',1),('cap',0),('cap',1))]
                spatial=[retarded_response(maps[g,side],'spatial',order)
                    for g,side in (('nested',0),('nested',1),('cap',0),('cap',1))]
                errors=[]
                for f,d in zip(fields,spatial):
                    W=float(preparations['cap'].fields(*f.lengths)['potential'])
                    coefficient=logarithmic_coefficients(f,d)
                    errors.extend((abs(coefficient['q_alpha']-W*W/(2*np.pi)),
                        abs(coefficient['angular']+W/np.pi)))
                key=(128,order,fraction,source_jet);identity_errors[key]=max(errors)
                if source_jet and identity_errors[key]>1e-9:
                    raise RuntimeError('Source-compatible logarithmic identity failed: '+str((key,identity_errors[key])))
                lf,rf,ls,rs=fields;ld,rd,sd,td=spatial
                arguments=((lf,rf,ls,rs,1/.45),(ld,rf,sd,rs,1/.45),(lf,rd,ls,td,1/.45))
                jobs.extend((key,channel,arg) for channel,arg in enumerate(arguments))
    values={};pending={}
    try:
        with ProcessPoolExecutor(max_workers=8,mp_context=multiprocessing.get_context('spawn')) as pool:
            pending={pool.submit(pair_worker,job):job for job in jobs}
            try:
                for future in as_completed(pending):
                    result=future.result()
                    worker_totals[result['pid']]=max(worker_totals.get(result['pid'],0.),result['process_CPU_total'])
                    budget()
                    if result['error'] is not None:
                        raise RuntimeError(result['error'])
                    key=tuple(result['key']);channels=values.setdefault(key,{})
                    channels[result['channel']]=result['value']
                    if len(channels)==3:
                        n,order,fraction,source_jet=key;distance,error=separations[fraction]
                        row=statistic_row(n,order,fraction,a,r,distance,error,channels[0],channels[1],channels[2])
                        row.update(source_endpoint_derivative=source_jet,endpoint_log_identity_max_error=identity_errors[key])
                        rows.append(row)
                        print(json.dumps(dict(phase='covariance',cells=n,polynomial_order=order,
                            separation_over_S=fraction,source_endpoint_derivative=source_jet)),flush=True)
            except Exception:
                for future in pending:
                    future.cancel()
                raise
    finally:
        # The pool has stopped before these receipts are collected, even
        # when another case failed. A future result is cached, not rerun.
        for future in pending:
            if future.done() and not future.cancelled() and future.exception() is None:
                result=future.result()
                worker_totals[result['pid']]=max(worker_totals.get(result['pid'],0.),result['process_CPU_total'])


def diagnose_source_logs(commit,hashes):
    """One fixed-grid coefficient diagnosis; no new bilocal/finite stress."""
    started=time.process_time();wall=time.perf_counter()
    chart=WallRadarChart();S=chart.duration;U,V=.8*S,.7*S
    metric=chart.fields(U,V,region='cap');a=float(metric['lapse']);r=float(metric['radius'])
    preparations={g:CylinderPreparation(chart,duration=S,root_radius=.45,region=g,angular_k=1)
        for g in ('nested','cap')}
    samples={};cache={}
    for g,preparation in preparations.items():
        samples[g]=endpoint(preparation,U,V,128,cache,order=48)
        if time.process_time()-started>30:
            raise RuntimeError('Source-log diagnosis exhausted its 30 process-CPU seconds')
    path=ROOT/'results/summaries/nsc-atlas-dirac-root-response-v2.json';raw=path.read_bytes()
    if raw!=subprocess.check_output(['git','show',f'{commit}:{path.relative_to(ROOT).as_posix()}'],cwd=ROOT):
        raise ValueError('freeze the diagnostic comparison record first')
    previous=json.loads(raw)
    observed={(v['polynomial_order'],v['separation_over_S']):v['finite_split_anisotropy_difference']
        for v in previous['cases']}
    rows=[];radial=2/(4*np.pi*a*a*r*r);angular=2/(8*np.pi*a*r**3)
    for order in (24,48):
        coefficients={};corrected={};errors={}
        for g,sample in samples.items():
            f=retarded_response(sample,'field',order);d=retarded_response(sample,'spatial',order)
            coefficients[g]=logarithmic_coefficients(f,d)
            corrected[g]=logarithmic_coefficients(f,d,source_slopes=sample['source_endpoint_values'])
            errors[g]=[float(abs(f.value(leg,L,derivative=1)[leg]-sample['source_endpoint_values'][leg,leg]))
                for leg,L in enumerate(f.lengths)]
        def pressure(values):
            return radial*(values['nested']['q_alpha']-values['cap']['q_alpha'])\
                -angular*(values['nested']['angular']-values['cap']['angular'])
        residue=pressure(coefficients)
        rows.append(dict(polynomial_order=order,represented_log_coefficients=coefficients,
            exact_source_log_coefficients=corrected,diagonal_source_jet_error=errors,
            relative_pressure_log_coefficient=residue,
            coefficient_with_exact_source_slopes=pressure(corrected),
            predicted_change_for_halved_separation=residue*np.log(2),
            observed_v2_change=observed[order,.01]-observed[order,.02]))
    result=dict(diagnostic='source_endpoint_log_residue_only',method_commit=commit,
        source_fingerprint=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest(),
        comparison_sha256=hashlib.sha256(raw).hexdigest(),cells=128,angular_k=1,
        detector_U_over_S=.8,detector_V_over_S=.7,rows=rows,
        compute=dict(process_CPU_seconds=time.process_time()-started,CPU_limit_seconds=30,
            wall_seconds=time.perf_counter()-wall),
        scope='Coefficient substitution diagnoses a missing identity; no covariance or finite stress is changed.')
    print(json.dumps(result,sort_keys=True,allow_nan=False),flush=True)


def diagnose_observer_moments(commit,hashes):
    """Fixed-probe compatibility on the same maps; no vacuum contraction."""
    started=time.process_time();wall=time.perf_counter()
    chart=WallRadarChart();S=chart.duration;U,V=.8*S,.7*S;degrees=(0,1,2,4,8)
    preparations={g:CylinderPreparation(chart,duration=S,root_radius=.45,region=g,angular_k=1)
        for g in ('nested','cap')}
    cache={};residuals={}
    for fraction in (.02,.01):
        eta=fraction*S;maps={}
        for g,preparation in preparations.items():
            for side,point in enumerate(((U-eta/2,V+eta/2),(U+eta/2,V-eta/2))):
                maps[g,side]=endpoint(preparation,*point,128,cache,order=48)
                if time.process_time()-started>60:
                    raise RuntimeError('Observer-moment diagnosis exhausted its 60 process-CPU seconds')
        for order in (24,48):
            for source_jet in (False,True):
                for g in preparations:
                    fields=[retarded_response(maps[g,side],'field',order,source_jet=source_jet) for side in (0,1)]
                    spatial=[retarded_response(maps[g,side],'spatial',order) for side in (0,1)]
                    residuals[g,order,source_jet,fraction]=observation_moment_residual(
                        fields[0],spatial[0],fields[1],spatial[1],separation=eta,
                        reference_duration=S,degrees=degrees)
                residuals['nested_minus_cap',order,source_jet,fraction]=(
                    residuals['nested',order,source_jet,fraction]-residuals['cap',order,source_jet,fraction])
    def norm(value):
        return np.max(abs(value),axis=(0,2))
    rows=[];ratios=[];floor=1e-11
    for g in ('nested','cap','nested_minus_cap'):
        for order in (24,48):
            correction={f:residuals[g,order,True,f]-residuals[g,order,False,f] for f in (.02,.01)}
            for source_jet in (False,True):
                for fraction in (.02,.01):
                    rows.append(dict(geometry=g,polynomial_order=order,
                        source_endpoint_derivative=source_jet,separation_over_S=fraction,
                        max_abs_residual_by_probe=norm(residuals[g,order,source_jet,fraction]).tolist(),
                        correction_induced_residual_by_probe=norm(correction[fraction]).tolist() if source_jet else None))
                large=norm(residuals[g,order,source_jet,.02]);small=norm(residuals[g,order,source_jet,.01])
                ratios.append(dict(geometry=g,polynomial_order=order,source_endpoint_derivative=source_jet,
                    large_over_small_by_probe=[float(a/b) if max(a,b)>floor and b>0 else None
                        for a,b in zip(large,small)]))
    result=dict(diagnostic='ordinary_kernel_observer_moment_compatibility',method_commit=commit,
        source_fingerprint=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest(),
        cells=128,angular_k=1,detector_U_over_S=.8,detector_V_over_S=.7,
        probe_degrees=degrees,fixed_source_reference_duration=S,reporting_floor=floor,
        rows=rows,split_ratios=ratios,
        compute=dict(process_CPU_seconds=time.process_time()-started,CPU_limit_seconds=60,
            wall_seconds=time.perf_counter()-wall,mode_workers=8,retained_full_arrays=False),
        scope='Necessary ordinary-kernel derivative check; central-difference, field and projection errors coexist. No stress bound or new finite vacuum amplitude.')
    print(json.dumps(result,sort_keys=True,allow_nan=False),flush=True)


def diagnose_projection_derivative(commit,hashes,*,discrete=False):
    """Conditional projection derivative versus independently evolved row."""
    started=time.process_time();wall=time.perf_counter()
    chart=WallRadarChart();S=chart.duration;U,V=.8*S,.7*S;degrees=(0,1,2,4,8)
    preparations={g:CylinderPreparation(chart,duration=S,root_radius=.45,region=g,angular_k=1)
        for g in ('nested','cap')}
    cache={};rows=[];residuals={};endpoint_errors=[];quadrature_errors=[]
    for fraction in (.02,.01):
        eta=fraction*S;maps={}
        for g,preparation in preparations.items():
            for side,point in enumerate(((U-eta/2,V+eta/2),(U+eta/2,V-eta/2))):
                maps[g,side]=endpoint(preparation,*point,128,cache,order=48,observer_jets=True,
                    observation_parameter_jet=discrete)
                quadrature_errors.append(maps[g,side]['endpoint_quadrature_estimate'])
                if time.process_time()-started>90:
                    raise RuntimeError('Projection derivative diagnosis exhausted its 90 process-CPU seconds')
        for order in (24,48):
            commutators={}
            for g in preparations:
                fields=[retarded_response(maps[g,side],'field',order,source_jet=True) for side in (0,1)]
                old=[retarded_response(maps[g,side],'spatial',order) for side in (0,1)]
                new=[differentiated_retarded_response(maps[g,side],order,discrete=discrete) for side in (0,1)]
                commutators[g]=np.array([[new[side].green[leg]-old[side].green[leg]
                    for leg in (0,1)] for side in (0,1)])
                for side in (0,1):
                    for leg in (0,1):
                        endpoint_errors.append(float(np.max(abs(commutators[g][side,leg].sum(axis=0)))))
                for method,spatial in (('independent_spatial',old),('projection_derivative',new)):
                    residuals[g,order,fraction,method]=observation_moment_residual(
                        fields[0],spatial[0],fields[1],spatial[1],
                        separation=eta,reference_duration=S,degrees=degrees)
            commutators['nested_minus_cap']=commutators['nested']-commutators['cap']
            for method in ('independent_spatial','projection_derivative'):
                residuals['nested_minus_cap',order,fraction,method]=(
                    residuals['nested',order,fraction,method]-residuals['cap',order,fraction,method])
            for g in ('nested','cap','nested_minus_cap'):
                moment=np.array([[endpoint_hilbert_moment(commutators[g][side,leg])
                    for leg in (0,1)] for side in (0,1)])
                rows.append(dict(geometry=g,polynomial_order=order,separation_over_S=fraction,
                    residuals_by_probe={method:np.max(abs(residuals[g,order,fraction,method]),axis=(0,2)).tolist()
                        for method in ('independent_spatial','projection_derivative')},
                    commutator_endpoint_subtracted_Hilbert_norm=float(np.max(abs(moment))),
                    commutator_mean_norm=float(np.max(abs(commutators[g][:,:,0])))))
    ratios=[]
    for g in ('nested','cap','nested_minus_cap'):
        for order in (24,48):
            for method in ('independent_spatial','projection_derivative'):
                a=np.max(abs(residuals[g,order,.02,method]),axis=(0,2))
                b=np.max(abs(residuals[g,order,.01,method]),axis=(0,2))
                ratios.append(dict(geometry=g,polynomial_order=order,method=method,
                    large_over_small_by_probe=[float(x/y) if max(x,y)>1e-11 and y else None for x,y in zip(a,b)]))
    if max(endpoint_errors)>1e-9:
        raise RuntimeError('Differentiated projection failed its exact endpoint trace')
    result=dict(diagnostic='finite_map_projection_derivative' if discrete else
        'conditional_constrained_projection_derivative',method_commit=commit,
        source_fingerprint=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest(),
        cells=128,angular_k=1,detector_U_over_S=.8,detector_V_over_S=.7,
        probe_degrees=degrees,fixed_source_reference_duration=S,rows=rows,split_ratios=ratios,
        maximum_endpoint_trace_error=max(endpoint_errors),maximum_ray_quadrature_estimate=max(quadrature_errors),
        compute=dict(process_CPU_seconds=time.process_time()-started,CPU_limit_seconds=90,
            wall_seconds=time.perf_counter()-wall,mode_workers=8,retained_full_arrays=False),
        scope=('Analytic derivative of the declared finite reconstruction; continuum error remains uncontrolled. '
            if discrete else 'Exact continuum projection identity evaluated on cell maps; not a proof of their derivative. ')
            +'No covariance, finite stress, angular completion or source-induced slip.')
    print(json.dumps(result,sort_keys=True,allow_nan=False),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--stage',choices=('initial','polynomial','source-jets'),default='initial')
    parser.add_argument('--diagnose-source-logs',action='store_true')
    parser.add_argument('--diagnose-observer-moments',action='store_true')
    parser.add_argument('--diagnose-projection-derivative',action='store_true')
    parser.add_argument('--diagnose-discrete-derivative',action='store_true')
    args=parser.parse_args()
    diagnostic=any((args.diagnose_source_logs,args.diagnose_observer_moments,
        args.diagnose_projection_derivative,args.diagnose_discrete_derivative))
    if sum((args.diagnose_source_logs,args.diagnose_observer_moments,
            args.diagnose_projection_derivative,args.diagnose_discrete_derivative))>1:
        parser.error('choose one coefficient, moment or projection diagnosis')
    if diagnostic and (args.output is not None or args.stage!='initial'):
        parser.error('diagnosis prints its small table; do not combine it with output or stage')
    if not diagnostic and args.output is None:
        parser.error('output required for a physical comparison')
    if args.output is not None and (args.output.exists() or not args.output.parent.is_dir()):
        raise ValueError('unused output path and existing parent required')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    hashes={}
    for name in SOURCES:
        raw=(ROOT/name).read_bytes()
        if raw!=subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT):
            raise ValueError('freeze the response method before execution: '+name)
        hashes[name]=hashlib.sha256(raw).hexdigest()
    if args.diagnose_source_logs:
        diagnose_source_logs(commit,hashes)
        return
    if args.diagnose_observer_moments:
        diagnose_observer_moments(commit,hashes)
        return
    if args.diagnose_projection_derivative:
        diagnose_projection_derivative(commit,hashes)
        return
    if args.diagnose_discrete_derivative:
        diagnose_projection_derivative(commit,hashes,discrete=True)
        return
    cpu_start=time.process_time();wall_start=time.perf_counter();worker_totals={}
    def budget():
        if time.process_time()-cpu_start+sum(worker_totals.values())>CPU_LIMIT:
            raise RuntimeError('Full-root response exhausted its 240 process-CPU seconds; stop before another calculation')
    chart=WallRadarChart();S=chart.duration;U,V=.8*S,.7*S
    T,x=(U+V)/2,(V-U)/2;ell=10000.;mu=1/.45
    metric=chart.fields(U,V,region='cap');a=float(metric['lapse']);r=float(metric['radius'])
    preparations={g:CylinderPreparation(chart,duration=S,root_radius=.45,region=g,angular_k=1)
                  for g in ('nested','cap')}
    rows=[];endpoint_evidence=[];stopped=None;integral_cache={}
    baseline=None;baseline_gaps={}
    if args.stage in ('polynomial','source-jets'):
        path=ROOT/('results/summaries/nsc-atlas-dirac-root-response-v1.json' if args.stage=='polynomial' else
            'results/summaries/nsc-atlas-dirac-root-response-v2.json')
        raw=path.read_bytes()
        if raw!=subprocess.check_output(['git','show',f'{commit}:{path.relative_to(ROOT).as_posix()}'],cwd=ROOT):
            raise ValueError('freeze the numerical baseline before comparison')
        baseline=dict(path=path.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(raw).hexdigest(),
            absolute_tolerance=1e-9)
        baseline_values=json.loads(raw)
    try:
        with threadpool_limits(limits=1):
            if args.stage in ('polynomial','source-jets'):
                polynomial_calculation(chart,preparations,a,r,rows,endpoint_evidence,worker_totals,budget,
                    source_comparison=args.stage=='source-jets')
                expected={(v['polynomial_order'],v['separation_over_S']):v['finite_split_anisotropy_difference']
                    for v in baseline_values['cases'] if v['cells']==128}
                baseline_orders=(24,) if args.stage=='polynomial' else (24,48)
                baseline_gaps={f"{v['polynomial_order']}:{v['separation_over_S']}":abs(v['finite_split_anisotropy_difference']
                    -expected[v['polynomial_order'],v['separation_over_S']]) for v in rows
                    if v['polynomial_order'] in baseline_orders and not v.get('source_endpoint_derivative',False)}
                if len(baseline_gaps)!=2*len(baseline_orders) or max(baseline_gaps.values())>baseline['absolute_tolerance']:
                    raise RuntimeError('Parallel/order-24 baseline mismatch: '+str(baseline_gaps))
                budget()
            for n in (64,128) if args.stage=='initial' else ():
                for fraction in (.02,.01):
                    eta=fraction*S
                    points=((U-eta/2,V+eta/2),(U+eta/2,V-eta/2))
                    maps={}
                    for g,preparation in preparations.items():
                        for side,point in enumerate(points):
                            budget();sample=endpoint(preparation,*point,n,integral_cache);budget()
                            maps[g,side]=sample
                            endpoint_evidence.append(dict(cells=n,separation_over_S=fraction,
                                geometry=g,side=side,quadrature_estimate=sample['endpoint_quadrature_estimate'],
                                reconstruction_projection_residual=sample['cell_projection_residual'].tolist(),
                                response_sha256=hashlib.sha256(b''.join(
                                    c.tobytes() for key in ('cell_field','cell_spatial')
                                    for c in sample[key].green)).hexdigest()))
                            print(json.dumps(dict(cells=n,separation_over_S=fraction,geometry=g,side=side)),flush=True)
                    separation,separation_error=quad_vec(lambda y:float(chart.fields(T-y,T+y,region='cap')['lapse']),
                        x-eta/2,x+eta/2,epsabs=1e-13,epsrel=1e-12)
                    for order in (12,24):
                        budget()
                        lf,rf,ls,rs=[retarded_response(maps[g,side],'field',order)
                            for g,side in (('nested',0),('nested',1),('cap',0),('cap',1))]
                        ld,rd,sd,td=[retarded_response(maps[g,side],'spatial',order)
                            for g,side in (('nested',0),('nested',1),('cap',0),('cap',1))]
                        K=relative_pair(lf,rf,ls,rs,mu)
                        DL=relative_pair(ld,rf,sd,rs,mu)
                        DR=relative_pair(lf,rd,ls,td,mu)
                        rows.append(statistic_row(n,order,fraction,a,r,separation,separation_error,K,DL,DR))
                        budget()
    except Exception as error:
        stopped=type(error).__name__+': '+str(error)
    comparisons={}
    grids=(64,128) if args.stage=='initial' else (128,)
    orders=(12,24) if args.stage=='initial' else (24,48)
    if stopped is None:
        lookup={(v['cells'],v['polynomial_order'],v['separation_over_S'],v.get('source_endpoint_derivative',False)):
            v['finite_split_anisotropy_difference'] for v in rows}
        if args.stage=='source-jets':
            comparisons=dict(grid_spread=None,
                polynomial_spread={str(s):{str(f):abs(lookup[128,48,f,s]-lookup[128,24,f,s])
                    for f in (.02,.01)} for s in (False,True)},
                separation_spread={str(s):{str(p):abs(lookup[128,p,.01,s]-lookup[128,p,.02,s])
                    for p in orders} for s in (False,True)},
                source_constraint_change={f'{p}:{f}':lookup[128,p,f,True]-lookup[128,p,f,False]
                    for p in orders for f in (.02,.01)},baseline_absolute_gaps=baseline_gaps)
        else:
            comparisons=dict(grid_spread={f'{p}:{f}':abs(lookup[128,p,f,False]-lookup[64,p,f,False])
                for p in orders for f in (.02,.01)} if len(grids)==2 else None,
                polynomial_spread={f'{n}:{f}':abs(lookup[n,orders[-1],f,False]-lookup[n,orders[0],f,False])
                    for n in grids for f in (.02,.01)},
                separation_spread={f'{n}:{p}':abs(lookup[n,p,.01,False]-lookup[n,p,.02,False]) for n in grids for p in orders},
                baseline_absolute_gaps=baseline_gaps)
    version={'initial':'v1','polynomial':'v2','source-jets':'v3'}[args.stage]
    result=dict(schema='nsc-atlas-dirac-full-root-finite-split-'+version,method_commit=commit,stage=args.stage,
        sources_sha256=hashes,status='finite_split_reference_only' if stopped is None else 'stopped_before_complete_comparison',
        stopped_reason=stopped,domain=dict(wall_clock_duration=S,root_radius=.45,
            detector_U_over_S=.8,detector_V_over_S=.7,physical_reference_length_m=ell,
            detector_radius_over_ell=r,physical_radar_lapse=a,central_mass=0.,angular_sectors=[1],
            magnetic_degeneracy=2,root_radial_spectrum='complete analytic covariance, no momentum cutoff',
            state='same correlated product-cylinder parent and fixed auxiliary switch in both geometries'),
        cases=sorted(rows,key=lambda v:(v['cells'],v['polynomial_order'],v['separation_over_S'],
            v.get('source_endpoint_derivative',False))),
        endpoint_evidence=endpoint_evidence,numerical_spreads=comparisons,baseline_reference=baseline,
        limits=['One signed angular sector only; the angular tail is not computed.',
            'Finite nonzero split in the declared canonical radar frame; no coincident renormalized stress.',
            'The pressure-shaped bilocal statistic is not full Psi/Phi or its compact-source derivative.',
            'Grid, polynomial and separation spreads are diagnostics, not certified continuum error bounds.',
            'Bessel tail and PV estimates exclude field representation and floating arithmetic.',
            'A corrected leading logarithm does not establish the finite observation-derivative compatibility or coincidence part.',
            'The physical formation of the supplied preparation remains open; no P001 grid ran.'],
        compute=dict(mode_workers=8,covariance_processes=8 if args.stage!='initial' else 1,
            parallel_phases='separate field and covariance phases',BLAS_threads=1,
            process_CPU_seconds=time.process_time()-cpu_start+sum(worker_totals.values()),
            child_process_CPU_seconds=sum(worker_totals.values()),observed_child_processes=len(worker_totals),
            CPU_limit_seconds=CPU_LIMIT,wall_seconds=time.perf_counter()-wall_start,
            live_field_array_ceiling_bytes=256*1024**2,retained_full_arrays=False))
    raw=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>24000:
        raise RuntimeError('Compact full-root summary exceeds its 24000-byte ceiling')
    with args.output.open('x',encoding='utf-8',newline='\n') as stream:
        stream.write(raw)
    print(json.dumps(dict(output=str(args.output),bytes=len(raw.encode()),status=result['status'],
        stopped_reason=stopped,compute=result['compute'],spreads=comparisons)),flush=True)
    if stopped is not None:
        raise SystemExit(1)


if __name__=='__main__':
    main()
