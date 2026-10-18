"""Unexecuted preparation for a clocked causal geometry-resolution application.

Runner failure-gate controls are unfinished. The missing control file in
EXTRA_SOURCES deliberately prevents authentication and a physical application.
The current priority is the cap--wall slip and P001 mass-gap diagnosis.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from io import BytesIO
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_info,threadpool_limits

import measure_nsc_conformal_causal_response as old
from recursive_horizons.nsc_conformal_causal_resolution import (
    CANONICAL_FIELDS,combine_variations,clocked_canonical_jvp,
    fixed_receiver_clock_redshift,initial_packet,observer_clock_direction,
    observer_values,proper_time_reference,reconstructed_tangent,
    resolution_rk4_step,retained_projection,to_physical_tangent,
)
from recursive_horizons.nsc_spherical_coupling import CauchyState
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    build_grid,compose_fine_hamiltonian,constraint_diagnostics,
)


ROOT=Path(__file__).resolve().parents[1]
RESOLUTION='results/summaries/nsc-conformal-resolution-v1.json'
PHOTON='results/summaries/nsc-conformal-causal-response-v2.json'
SETUP_CPU,ANALYSIS_CPU,TOTAL_CPU=140.0,160.0,6000.0
BYTE_LIMIT,RAW_LIMIT=24000,2000000
GATE,CANCELLATION_LIMIT=1e-8,1e8
CASES=({'name':'nf512_dt_0.0010','dt':.001,'cpu_seconds':1900.0},
       {'name':'nf512_dt_0.0005','dt':.0005,'cpu_seconds':3800.0})
EXTRA_SOURCES=(
    'scripts/measure_nsc_conformal_causal_resolution.py',
    'tests/test_measure_nsc_conformal_causal_resolution.py',
    'src/recursive_horizons/nsc_conformal_causal_resolution.py',
    'tests/test_nsc_conformal_causal_resolution.py',
    'src/recursive_horizons/nsc_conformal_resolution.py',RESOLUTION,PHOTON,
)


def authenticate(geometry,reference_raw):
    commit,bindings,episode=old.authenticate(geometry)
    for name in EXTRA_SOURCES:
        raw=(ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT)!=raw:
            raise ValueError('freeze the causal resolution method before application: '+name)
        bindings[name]=hashlib.sha256(raw).hexdigest()
    reference=json.loads((ROOT/PHOTON).read_text(encoding='utf-8'))
    reduction=json.loads((ROOT/RESOLUTION).read_text(encoding='utf-8'))
    if (reference['status']!=old.COMPLETED or old.sha256(reference_raw)!=reference['raw']['sha256']
            or reference['input']['payload_sha256']!=old.EXPECTED_PAYLOAD
            or reduction['input']['payload_sha256']!=old.EXPECTED_PAYLOAD):
        raise ValueError('the original geometry/tangent cache is not authenticated')
    for name,digest in reference['source_sha256'].items():
        if bindings.get(name)!=digest:raise ValueError('the old causal source differs: '+name)
    domain,projector=reduction['domain'],reduction['projector']
    if (domain['target_observation_used'] is not False or domain['reference_T']!=.05
            or (domain['nf'],domain['ng'],domain['nq'])!=(512,511,2048)
            or projector['retained_modes']!=1 or projector['unresolved_modes']!=510
            or projector['constant_mode_residual']>GATE
            or projector['smallest_excluded_eigenvalue']*domain['reference_length']**2<=1):
        raise ValueError('the fixed physical reference no longer selects only the constant geometry mode')
    return commit,bindings,reference,reduction


def block_gap(actual,reference):
    return old.block_relative(actual,reference)


def closure_report(packet,reference):
    combined=reconstructed_tangent(packet)
    if any(not np.isfinite(getattr(v,n)).all()
           for v in (combined,reference,packet.initial_drive,packet.feedback) for n in CANONICAL_FIELDS):
        raise ValueError('a causal-resolution variation is nonfinite')
    gaps=[block_gap(getattr(combined,n),getattr(reference,n)) for n in CANONICAL_FIELDS]
    amplification=[]
    for n in CANONICAL_FIELDS:
        h,f=getattr(packet.initial_drive,n),getattr(packet.feedback,n)
        scale=max(1.0,float(np.max(abs(getattr(reference,n)))))
        amplification.append(float((np.max(abs(h))+np.max(abs(f)))/scale))
    return {'max_abs':[g['max_abs'] for g in gaps],'relative':[g['relative'] for g in gaps],
            'amplification':amplification,
            'accepted':bool(max(g['relative'] for g in gaps)<=GATE
                            and max(amplification)<=CANCELLATION_LIMIT)}


def event_sample(grid,packet,rate,index):
    row=grid.A_g[index]
    lapse=observer_values(packet.state,row)[2]
    parts={'retained':packet.retained,'initial_drive':packet.initial_drive,'feedback':packet.feedback}
    return {'N':lapse,
        'coordinate_dlogN':observer_clock_direction(packet.state,packet.coordinate_tangent,row)/lapse,
        'dotlogN':observer_clock_direction(packet.state,rate,row)/lapse,
        'part_dlogN':{n:observer_clock_direction(packet.state,to_physical_tangent(packet.state,v),row)/lapse
                     for n,v in parts.items()}}


def photon_rows(receiver,emitters):
    if receiver is None or len(emitters)!=2:return []
    rows=[]
    for direction,e in sorted(emitters.items()):
        logvalue=fixed_receiver_clock_redshift(receiver['coordinate_dlogN'],e['coordinate_dlogN'],
            receiver['dotlogN'],e['dotlogN'],receiver['clock_tangent'],receiver['N'])
        parts={n:receiver['part_dlogN'][n]-e['part_dlogN'][n] for n in receiver['part_dlogN']}
        parts['clock']=e['dotlogN']*(receiver['anchor_shift']-e['anchor_shift'])
        assembled=sum(parts.values());gap=abs(assembled-logvalue)
        factor=receiver['N']/e['N']
        rows.append({'direction':direction,'one_plus_z':factor,
            'coordinate_log_derivative':receiver['coordinate_dlogN']-e['coordinate_dlogN'],
            'proper_log_derivative':logvalue,'proper_redshift_derivative':factor*logvalue,
            'proper_log_parts':parts,'assembly_gap':gap,
            'accepted':bool(gap<=GATE*max(1,abs(logvalue),abs(assembled)))})
    return rows


def diagnostic(grid,packet,row,mark,steps):
    rate,bundle=compose_fine_hamiltonian(grid,packet.state)
    reference=proper_time_reference(packet,rate,row)
    closure=closure_report(packet,reference)
    constraints=constraint_diagnostics(grid,packet.state,bundle['source'],bundle['fine_state'])
    def apply(v):return clocked_canonical_jvp(grid,packet.state,v,row,rate=rate,bundle=bundle)
    initial=apply(packet.initial_drive)
    coupled=apply(combine_variations((1,packet.feedback),(1,packet.retained)))
    instantaneous=apply(packet.retained)
    lapse=observer_values(packet.state,row)[2]
    def small(v):
        p=retained_projection(v)
        return [float(np.mean(getattr(p,n))/lapse) for n in CANONICAL_FIELDS[:6]]+[
            float(np.sqrt(np.mean(abs(getattr(p,n))**2))/lapse) for n in CANONICAL_FIELDS[6:]]
    rates={'instantaneous':small(instantaneous),'initial_drive':small(initial),
           'feedback':small(combine_variations((1,coupled),(-1,instantaneous)))}
    physical=to_physical_tangent(packet.state,reconstructed_tangent(packet))
    info={'T':mark,'steps':steps,'observer_tau':packet.clock,'delta_tau':packet.clock_tangent,
        'observer_N':lapse,'closure':closure,'coarse_rates_per_tau':rates,
        'full_hamilton_max':float(constraints['full_hamilton_max']),
        'full_momentum_max':float(constraints['full_momentum_max']),
        'gram_max':float(np.max(abs(old.column_gram(packet.state.phi0,packet.state.phi1)-np.eye(6)))),
        'coordinate_tangent_gram_max':float(np.max(abs(old.tangent_gram(packet.state,packet.coordinate_tangent)))),
        'proper_tangent_gram_max':float(np.max(abs(old.tangent_gram(packet.state,physical))))}
    return info,rate


def evolve_case(grid,state,tangent,spec,started):
    clock=old.CpuClock(spec['cpu_seconds'],started)
    dt=float(spec['dt']);required=int(round(old.TARGET_TIME/dt))
    marks={int(round(t/dt)):t for t in old.DIAGNOSTIC_TIMES}
    nodes=old.production_nodes(grid.xi_q);row=grid.A_g[nodes['receiver']]
    packet=initial_packet(state,tangent)
    diagnostics=[];checkpoints={};receiver=None;emitters={};forecast=None;reason=None;steps=0
    def record():
        nonlocal receiver,reason
        mark=marks[steps]
        info,rate=diagnostic(grid,packet,row,mark,steps);diagnostics.append(info)
        if not info['closure']['accepted']:reason='closure_or_cancellation_gate_failed'
        if mark in (.05,.3):checkpoints[f'{mark:.3f}']=old.pack_state(packet.state)
        shift=packet.clock_tangent/info['observer_N']
        if mark==old.EMISSION_TIME:
            for d,index in nodes['emitters'].items():
                emitters[d]={**event_sample(grid,packet,rate,index),'anchor_shift':shift}
        if mark==old.RECEPTION_TIME:
            receiver={**event_sample(grid,packet,rate,nodes['receiver']),
                'clock_tangent':packet.clock_tangent,'anchor_shift':shift}
        print(f"{spec['name']}: T={mark:g}, steps={steps}, CPU={clock.used:.3f}s, "
              f"closure={max(info['closure']['relative']):.3g}",flush=True)
    try:
        clock.check('initial diagnostic');record();clock.check('after initial diagnostic')
        while steps<required and reason is None:
            clock.check('before resolution step');before=time.process_time()
            packet=resolution_rk4_step(grid,packet,dt,row,checkpoint=lambda:clock.check('resolution stage'))
            cost=time.process_time()-before;steps+=1
            if steps in marks:record()
            if steps==1:
                accepted,projected=old.forecast_accepts(clock.used,cost,1,required,clock.seconds)
                forecast={'first_step_cpu':cost,'projected_cpu':projected,'accepted':bool(accepted)}
                print(f"{spec['name']}: first-step forecast {projected:.1f}/{clock.seconds:g} CPU s",flush=True)
                if not accepted:reason='first_step_forecast_exceeds_budget'
        clock.check('after final accepted diagnostic')
    except Exception as error:
        reason=f'{type(error).__name__}: {error}'
    reached=reason is None and steps==required
    rays=photon_rows(receiver,emitters)
    report={'name':spec['name'],'dt':dt,'T':steps*dt,'steps':steps,'target_reached':reached,
        'reason':reason,'CPU_seconds':clock.used,'forecast':forecast,'diagnostics':diagnostics,'rays':rays,
        'closure_accepted':bool(diagnostics and all(r['closure']['accepted'] for r in diagnostics)),
        'photon_accepted':bool(reached and len(rays)==2 and all(r['accepted'] for r in rays))}
    cache={'name':spec['name'],'time':steps*dt,'target_reached':reached,
        'state':old.pack_state(packet.state),'tangent':old.pack_state(packet.coordinate_tangent),
        'checkpoints':checkpoints,
        'split':{n:{f:np.array(getattr(getattr(packet,n),f),copy=True) for f in CANONICAL_FIELDS}
                 for n in ('retained','initial_drive','feedback')}}
    return {'report':old.plain(report),'cache':cache}


def run_case(payload):
    for name in old.THREAD_VARIABLES:os.environ[name]='4'
    started=time.process_time()
    try:
        with threadpool_limits(limits=4):
            grid=build_grid(512,quadrature=2048,gauge='conformal')
            state=CauchyState(**{n:np.array(payload['state'][n],copy=True) for n in old.STATE_NAMES})
            tangent=CauchyState(**{n:np.array(payload['tangent'][n],copy=True) for n in old.STATE_NAMES})
            result=evolve_case(grid,state,tangent,payload['spec'],started)
            result['report']['threads']=[p['num_threads'] for p in threadpool_info()]
            return result
    except Exception as error:
        return {'report':{'name':payload['spec']['name'],'target_reached':False,
            'CPU_seconds':time.process_time()-started,'reason':f'{type(error).__name__}: {error}'},'cache':None}


def case_succeeded(report):
    spec=next((s for s in CASES if s['name']==report.get('name')),None)
    cost=report.get('CPU_seconds')
    if not isinstance(cost,(float,int)) or not np.isfinite(cost) or cost<0:return False
    if (spec is None or report.get('dt')!=spec['dt']
            or report.get('steps')!=int(round(old.TARGET_TIME/spec['dt']))
            or abs(report.get('T',-1)-old.TARGET_TIME)>1e-12
            or cost>spec['cpu_seconds']
            or report.get('reason') is not None):return False
    diagnostics=report.get('diagnostics',[])
    if [r['T'] for r in diagnostics]!=list(old.DIAGNOSTIC_TIMES):return False
    for row in diagnostics:
        for name,limit in [('relative',GATE),('amplification',CANCELLATION_LIMIT)]:
            values=np.asarray(row['closure'][name])
            if values.shape!=(len(CANONICAL_FIELDS),) or not np.isfinite(values).all() or max(values)>limit:return False
    rays=report.get('rays',[])
    if len(rays)!=2 or {r['direction'] for r in rays}!={-1,1}:return False
    for ray in rays:
        if (not np.isfinite([ray['proper_redshift_derivative'],ray['proper_log_derivative'],
                             ray['one_plus_z'],ray['assembly_gap']]).all()
                or ray['one_plus_z']<=0 or not ray['accepted']):return False
    return bool(report.get('target_reached') and report.get('closure_accepted')
        and report.get('photon_accepted') and report.get('baseline',{}).get('within_tolerance')
        and report.get('old_tangent_replay',{}).get('accepted')
        and report.get('coordinate_photon_replay',{}).get('accepted'))


def summary_complete(reports,setup_cpu,analysis_cpu):
    costs=[r.get('CPU_seconds') for r in reports]
    known=all(isinstance(c,(float,int)) and np.isfinite(c) and c>=0 for c in costs)
    total=setup_cpu+analysis_cpu+sum(costs) if known else None
    passed=(len(reports)==2 and {r.get('name') for r in reports}=={s['name'] for s in CASES}
        and all(case_succeeded(r) for r in reports) and setup_cpu<=SETUP_CPU
        and analysis_cpu<=ANALYSIS_CPU and total is not None and total<=TOTAL_CPU)
    return bool(passed),total


def prepare_inputs(geometry,clock):
    """Regenerate the owned initial direction; reuse its frozen JVP certificate."""
    clock.check('before grid')
    grid=build_grid(512,quadrature=2048,gauge='conformal')
    old.production_nodes(grid.xi_q)
    with np.load(geometry,allow_pickle=False) as data:state=old._archive_state(data,grid)
    prepared=old.preparation_tangent(grid,state,checkpoint=lambda:clock.check('initial direction'))
    tangent=old.initial_tangent(prepared)
    couplings={name:float(getattr(grid.fine,name)) for name in ('A','C_W','C_F','flux')}
    couplings.update(kappa=int(grid.fine.kappa),multiplicity=int(grid.fine.multiplicity),
        occupations=[float(v) for v in grid.fine.occupations],changed=False)
    clock.check('after initial direction')
    return state,tangent,prepared.report,couplings


def replay_tangent(reference_raw,result):
    report,cache=result['report'],result.get('cache')
    if not cache or not report.get('target_reached'):return {'accepted':False}
    with np.load(reference_raw,allow_pickle=False) as data:
        gaps=[block_gap(cache['tangent'][n],data[report['name']+'_u_'+n]) for n in old.STATE_NAMES]
    return {'accepted':bool(max(g['relative'] for g in gaps)<=GATE),
            'max_abs':[g['max_abs'] for g in gaps],'relative':[g['relative'] for g in gaps]}


def write_raw(path,results):
    values={}
    for result in results:
        cache=result.get('cache')
        if not cache:continue
        prefix=cache['name'];values[prefix+'_T']=np.float64(cache['time'])
        for group in ('state','tangent'):
            for name,value in cache[group].items():values[prefix+'_'+group+'_'+name]=value
        for group,fields in cache['split'].items():
            for name,value in fields.items():values[prefix+'_'+group+'_'+name]=value
    if not values:return None
    stream=BytesIO();np.savez_compressed(stream,**values)
    if stream.tell()>RAW_LIMIT:raise ValueError('final raw cache exceeds its allocation')
    path.write_bytes(stream.getvalue())
    return {'path':str(path),'bytes':stream.tell(),'sha256':old.sha256(path),'full_history':False}


def temporal_comparison(reports):
    if len(reports)!=2 or not all(case_succeeded(r) for r in reports):return None
    coarse,fine=sorted(reports,key=lambda r:r['dt'],reverse=True)
    rows=[]
    for a,b in zip(coarse['rays'],fine['rays']):
        rows.append({'direction':a['direction'],'coarse':a['proper_redshift_derivative'],
            'fine':b['proper_redshift_derivative'],
            'absolute_difference':abs(a['proper_redshift_derivative']-b['proper_redshift_derivative'])})
    return {'proper_redshift_derivatives':rows,
        'receiver_tau_spread':abs(coarse['diagnostics'][-1]['observer_tau']-fine['diagnostics'][-1]['observer_tau']),
        'empirical_only':True,'continuum_or_physical_bound':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry',type=Path,default=old.DEFAULT_GEOMETRY)
    parser.add_argument('--reference-raw',type=Path,default=old.DEFAULT_RAW)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--raw',type=Path,default=ROOT/'.build/results/nsc-conformal-causal-resolution-v1.npz')
    args=parser.parse_args();old.require_fresh(args.output,'summary');old.require_fresh(args.raw,'raw')
    started=time.process_time();setup=old.CpuClock(SETUP_CPU,started)
    commit,bindings,old_result,resolution=authenticate(args.geometry,args.reference_raw)
    with threadpool_limits(limits=8):
        state,tangent,preparation,couplings=prepare_inputs(args.geometry,setup)
    if couplings!=old_result['couplings']:raise ValueError('the original coupled parameters differ')
    preflight=old_result['preflight']
    if not preflight['accepted']:raise ValueError('the frozen original JVP certificate is not accepted')
    setup_cpu=setup.used;setup.check('after setup')
    result={'schema':'nsc-conformal-causal-resolution-v1','status':'running',
        'source_commit':commit,'source_sha256':bindings,'couplings':couplings,
        'input':{'payload_sha256':old.EXPECTED_PAYLOAD,'old_tangent_raw_sha256':old_result['raw']['sha256']},
        'reference':{'L_ref':resolution['domain']['reference_length'],'reference_T':.05,
            'retained_geometry_modes':1,'unresolved_geometry_modes':510,
            'analytic_constant_projector':True,'old_constant_mode_residual':resolution['projector']['constant_mode_residual'],
            'all_fermions_retained':True,'observer_x':1.0,'clock':'observer proper time',
            'topology':'periodic S1 x S2','source_or_state_changed':False},
        'initial_constraint':preparation['background']['full_hamilton_max_abs'],
        'initial_clock_fraction':preparation['observer']['delta_N_over_N_per_epsilon'],
        'inherited_preflight':{'accepted':preflight['accepted'],
            'central_difference_accepted':preflight['central_difference_accepted'],
            'affine_radius_relative':preflight['initial_affine_radius']['relative'],'rerun':False},
        'setup_CPU_seconds':setup_cpu,
        'canonical_fields':list(CANONICAL_FIELDS),'physical_fields':list(old.STATE_NAMES),
        'coarse_rate_fields':list(CANONICAL_FIELDS),
        'coarse_rate_measure':'geometry means and separate fermion RMS, per epsilon and observer proper time',
        'limits':{'total_CPU_seconds':TOTAL_CPU,'setup_CPU_seconds':SETUP_CPU,
            'coarse_CPU_seconds':1900,'fine_CPU_seconds':3800,'analysis_CPU_seconds':ANALYSIS_CPU,
            'workers':2,'BLAS_threads_each':4,'summary_bytes':BYTE_LIMIT,'raw_bytes':RAW_LIMIT,
            'relative_closure_and_replay_gate':GATE,'normalized_split_amplification_limit':CANCELLATION_LIMIT},
        'cases':[]}
    def save():
        raw=json.dumps(old.plain(result),sort_keys=True,separators=(',',':'),allow_nan=False)+'\n'
        if len(raw.encode('utf-8'))>BYTE_LIMIT:raise ValueError('the compact summary exceeds its allocation')
        args.output.write_text(raw,encoding='utf-8',newline='\n')
    save()
    payloads=[{'state':old.pack_state(state),'tangent':old.pack_state(tangent),'spec':c} for c in CASES]
    results=[]
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(run_case,p) for p in payloads]
        for future,payload in zip(futures,payloads):
            try:results.append(future.result())
            except Exception as error:results.append({'report':{'name':payload['spec']['name'],
                'target_reached':False,'CPU_seconds':None,'reason':f'{type(error).__name__}: {error}'},'cache':None})
    analysis_started=time.process_time()
    baselines=old.baseline_reports(args.geometry,results)
    for item,baseline in zip(results,baselines):
        item['report']['baseline']={'within_tolerance':baseline['within_tolerance'],
            'max_relative_by_time':{t:max(r['fields'][n]['max_abs']/r['fields'][n]['scale'] for n in old.STATE_NAMES)
                for t,r in baseline['times'].items() if 'fields' in r}}
        item['report']['old_tangent_replay']=replay_tangent(args.reference_raw,item)
        prior=next(r for r in old_result['cases'] if r['name']==item['report']['name'])
        ray_gaps=[]
        for ray in item['report'].get('rays',[]):
            before=next(r for r in prior['rays'] if r['direction']==ray['direction'])
            ray_gaps.append(abs(ray['coordinate_log_derivative']-before['log_redshift_derivative']))
        item['report']['coordinate_photon_replay']={'accepted':len(ray_gaps)==2 and max(ray_gaps)<=GATE,
                                                   'absolute_gaps':ray_gaps}
    reports=[r['report'] for r in results]
    result.update(cases=reports,temporal_comparison=temporal_comparison(reports),raw=write_raw(args.raw,results))
    analysis_cpu=time.process_time()-analysis_started
    complete,total=summary_complete(reports,setup_cpu,analysis_cpu)
    result.update(status='completed_causal_geometric_resolution' if complete else 'partial_causal_geometric_resolution',
        analysis_CPU_seconds=analysis_cpu,total_CPU_seconds=total,
        scope=('Finite original local canonical conformal model, same state and unreset constraints. '
            'No nonlinear filtered trajectory. Retarded memory is the exact finite variational Schur split. '
            'The fixed-receiver-clock photon response is conditional on this supplied preparation. '
            'No full spectral-complement causal state, continuum error, physical-scale calibration, '
            'astrophysical comparison or absolute cosmological Lambda is inferred.'),
        uncertainty='Two fixed time steps give empirical temporal indicators including the baseline clock target; no total error bound.',
        counts={'new_spatial_grids':0,'eigensolves':0,'root_solves':0,'source_fits':0,'constraint_resets':0,
            'old_producer_reruns':0,'old_JVP_preflight_reruns':0,'physical_time_steps':[300,600]},
        reproduce=('python scripts/run_python.py scripts/measure_nsc_conformal_causal_resolution.py '
            '--output .build/causal-resolution-summary.json --raw .build/results/causal-resolution-replay.npz'))
    save()
    print(json.dumps({'status':result['status'],'total_CPU_seconds':total,
                      'temporal':result['temporal_comparison']}),flush=True)
    return 0 if complete else 1


if __name__=='__main__':raise SystemExit(main())
