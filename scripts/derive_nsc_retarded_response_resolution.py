#!/usr/bin/env python3
"""Bounded new-node sampling and targeted resolution; replay is array-only."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
import importlib.util
import json
from pathlib import Path
import signal
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('frozen_retarded_pilot',ROOT/'scripts/derive_nsc_retarded_compatible_response.py')
P=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(P)
from recursive_horizons.nsc_transmitting_history_modes import reference_fields_on_grid,characteristic_columns

OUTPUT='results/development/nsc-retarded-response-resolution.json'
ARTIFACT='results/development/artifacts/nsc-retarded-response-resolution.npz'
OWNED=('scripts/derive_nsc_retarded_response_resolution.py','tests/test_nsc_retarded_response_resolution.py','docs/nsc-retarded-response-resolution.md')
EXTRA=(P.OUTPUT,'src/recursive_horizons/nsc_common_time_bulk_split.py')
CPU_CAP=120.;RESAMPLING_TOL=3e-11
# Analytic finite-speed support formula supplied by the independent diagnosis;
# these decimal endpoints are diagnostic estimates, not directed bounds.
CAUSAL_TRACE_INTERVAL=(.05404795482456391,.24595204517543606)


def signature():
    return {**P.signature(),**{p:P.digest(p) for p in (*OWNED,*EXTRA)}}


def subgrid(x,phi,points):
    if points<2 or (len(x)-1)%(points-1) or phi.shape!=(2*len(x),12):
        raise ValueError('nested source grid and twelve spin-major columns required')
    stride=(len(x)-1)//(points-1)
    return x[::stride],phi.reshape(2,len(x),12)[:,::stride].reshape(2*points,12)


def resampling_error(x,phi,oldx,oldphi):
    xx,ff=subgrid(x,phi,801)
    if not np.array_equal(xx,oldx):raise ValueError('new reference grid does not exactly nest archived801')
    return P.MAXABS(ff-oldphi)


def one_case(points,steps,xmax,phimax,data,guard):
    x,phi=subgrid(xmax,phimax,points)
    owner=P.FourthOrderModePropagator(x,data['mass'].item(),data['angular'].item())
    direction=P.CompatibleRadiusDirection(P.axial_bump(guard['axial_center']),lambda z,n:0.,.007,.03)
    provider=P.CachedZeroRadiusProvider(x,direction)
    times=np.linspace(0.,.3,steps+1);energy=np.repeat(data['energies'],3)
    initial_tangents=np.zeros((1,*phi.shape),complex)
    right=np.array([points-1,2*points-1]);N,beta,q,r=owner.reference
    sat=-np.array([N[-1]/q[-1]-beta[-1],-N[-1]/q[-1]-beta[-1]])/owner.weights[-1]
    def incoming(t):
        forcing=np.zeros_like(phi);forcing[right]=sat[:,None]*phi[right]*np.exp(-1j*energy*t)[None,:]
        return forcing,np.zeros_like(initial_tangents)
    ix=np.flatnonzero(x==1.)
    if len(ix)!=1:raise ValueError('exact rho1 node required')
    rows=np.array([ix[0],points+ix[0]])
    run=P.evolve_prepared_field_jets(owner,times,provider,phi,initial_tangents,incoming,sample_rows=rows)
    domain=P.TransmittingDiracSeamDomain(lapse=N[ix[0]],radial_scale=q[ix[0]],shift=beta[ix[0]],radius=r[ix[0]],surface_id='fixed-rho1-retarded-resolution')
    _,inverse,_=domain.trace_map(np.ones(1));restriction=inverse[0]@P.MODE_TO_CURRENT/(r[ix[0]]*np.sqrt(q[ix[0]]))
    F=np.einsum('ab,tbc->tac',restriction,run['sampled_fields'])
    dF=np.einsum('ab,tbc->tac',restriction,run['sampled_tangent_fields'][0])
    rootw=np.repeat(np.sqrt(data['weights']/(2*np.pi)),3)
    weightedF=(F*rootw).reshape(-1,12);weighteddF=(dF*rootw).reshape(-1,12)
    kernels={}
    for name,key in (('covariance','source'),('projector','projector')):
        C=P.block_diag(*data[key])
        kernels[name+'_tangent']=P.restriction_covariance_tangent(weightedF,weighteddF,C,np.zeros_like(C))
    diagnostics={k:run[k] for k in ('norm_flux_residual','metric_flux_tangent_residual','auxiliary_identity_residual','tangent_boundary_trace')}
    diagnostics.update(initial_tangent_norm=P.MAXABS(initial_tangents),incident_tangent_norm=0.,source_covariance_tangent_norm=0.,
        initial_metric_change=P.MAXABS(provider.metric-owner.reference),source_closed_column_norm=P.MAXABS(phi[:,2::3]),
        source_covariance_closed_column_norm=P.MAXABS(data['source'][:,:,2]),
        source_covariance_hermiticity=P.MAXABS(data['source']-data['source'].swapaxes(-1,-2).conj()),
        frame_map_parameter_derivative=0.,stationary_phase_reapplied=False)
    return {'x':x,'times':times,'sample_rows':rows,'field':run['evolved_field'],'field_tangent':run['tangent_fields'][0],
        'trace':F,'trace_tangent':dF,'characteristic_trace_tangent':run['sampled_tangent_fields'][0],
        'restriction':restriction,**kernels,'diagnostics_json':np.frombuffer(json.dumps(diagnostics,sort_keys=True).encode(),np.uint8)}


def compare(c,f,kind,limit):
    rows={}
    for name,key in (('field','field_tangent'),('trace','trace_tangent'),('covariance','covariance_tangent')):
        coarse=c[key];fine=f[key]
        if kind=='spatial' and name=='field':
            fine=subgrid(f['x'],fine,len(c['x']))[1]
        if kind=='temporal':
            stride=(len(f['times'])-1)//(len(c['times'])-1)
            if not np.array_equal(c['times'],f['times'][::stride]):raise ValueError('time nodes are not nested')
            if name=='trace':fine=fine[::stride]
            if name=='covariance':fine=fine.reshape(len(f['times']),2,len(f['times']),2)[::stride,:,::stride,:].reshape(coarse.shape)
        error=P.MAXABS(coarse-fine);norm=P.MAXABS(fine);ratio=error/norm if norm else None
        rows[name]={'max_absolute_difference':error,'fine_compared_response_norm':norm,'difference_over_fine_norm':ratio,
                    'relative_indicator_limit':limit,'passes':ratio is not None and ratio<limit}
    return rows


def case_arrays(arrays,name):return {k[len(name)+1:]:v for k,v in arrays.items() if k.startswith(name+'/')}


def peak(values,times):
    magnitudes=abs(values)
    index=np.unravel_index(np.argmax(magnitudes),magnitudes.shape)
    return {'magnitude':float(magnitudes[index]),'PG_time':float(times[index[0]]),'source_column':int(index[-1])}


def outside_causal_support(values,times):
    outside=(times<CAUSAL_TRACE_INTERVAL[0])|(times>CAUSAL_TRACE_INTERVAL[1])
    weights=np.ones(len(times));weights[[0,-1]]=.5
    # Uniform-time trapezoidal sampled energy; common dt cancels in the ratio.
    energy=np.sum(abs(values)**2,axis=tuple(range(1,values.ndim)))
    total=float(weights@energy);outer=float(weights[outside]@energy[outside])
    return {'maximum_entry':P.MAXABS(values[outside]),'sampled_energy_fraction':outer/total if total else 0.,
            'interval_estimate':list(CAUSAL_TRACE_INTERVAL),'directed_bound':False,'acceptance_gate':False}


def analyze(arrays,meta,oldarrays):
    summaries={};comparisons={};historical={}
    for case in meta['completed']:
        a=case_arrays(arrays,case);x,phi=subgrid(arrays['reference/x'],arrays['reference/phi'],len(a['x']))
        E=np.repeat(arrays['source/energies'],3);t=a['times']
        stationary=phi*np.exp(-1j*E*t[-1])[None,:]
        trace0=a['restriction']@phi[a['sample_rows']]
        stationary_trace=trace0[None,:,:]*np.exp(-1j*t[:,None,None]*E[None,None,:])
        late=t>meta['support']['PG_time_support'][1]
        slow=a['characteristic_trace_tangent'][:,0,:]
        summaries[case]={'response_max_norms':{n:P.MAXABS(a[k]) for n,k in (('field','field_tangent'),('trace','trace_tangent'),('covariance','covariance_tangent'))},
            'covariance_hermiticity_residual':P.MAXABS(a['covariance_tangent']-a['covariance_tangent'].conj().T),
            'sampled_projector_kernel_tangent_norm':P.MAXABS(a['projector_tangent']),
            'sampled_projector_kernel_hermiticity_residual':P.MAXABS(a['projector_tangent']-a['projector_tangent'].conj().T),
            'stationary_reference_field_drift':P.MAXABS(a['field']-stationary),
            'stationary_reference_trace_drift':P.MAXABS(a['trace']-stationary_trace),
            'slow_characteristic_trace_peak':peak(slow,t),
            'slow_characteristic_trace_after_pulse':peak(slow[late],t[late]),
            'late_window_start':meta['support']['PG_time_support'][1],
            'characteristic_trace_outside_causal_support':outside_causal_support(a['characteristic_trace_tangent'],t),
            'slow_trace_outside_causal_support':outside_causal_support(slow,t),
            'diagnostics':json.loads(a['diagnostics_json'].tobytes())}
    if '3201_128' in meta['completed']:
        comparisons['spatial']=compare(case_arrays(arrays,'1601_128'),case_arrays(arrays,'3201_128'),'spatial',.25)
    if '3201_256' in meta['completed']:
        comparisons['temporal']=compare(case_arrays(arrays,'3201_128'),case_arrays(arrays,'3201_256'),'temporal',.1)
    for case in ('1601_128','3201_128'):
        if case in meta['completed']:
            historical['801_128_to_'+case]=compare(case_arrays(oldarrays,'801_128'),case_arrays(arrays,case),'spatial',.25)
    algebra=max((r['diagnostics'][k] for r in summaries.values() for k in ('norm_flux_residual','metric_flux_tangent_residual','auxiliary_identity_residual')),default=0.)
    herm=max((max(r['covariance_hermiticity_residual'],r['sampled_projector_kernel_hermiticity_residual']) for r in summaries.values()),default=0.)
    passed=(len(comparisons)==2 and all(v['passes'] for comp in comparisons.values() for v in comp.values())
            and algebra<3e-11 and herm<3e-11 and meta.get('resampling_residual',float('inf'))<RESAMPLING_TOL and not meta['cpu_budget_exceeded'])
    return summaries,comparisons,historical,passed,algebra,herm


def record_from_arrays(arrays,meta,digest,size,oldarrays):
    summaries,comparisons,historical,passed,algebra,herm=analyze(arrays,meta,oldarrays)
    return {'schema':'NSC-RETARDED-RESPONSE-RESOLUTION-v1','accountable_author':'Douglas Ek',
        'status':'PASS: sampled resolution indicators only' if passed else 'OPEN: targeted sampled resolution unresolved',
        'control':{'family':'14_1','amplitude':0.,'cases_completed':meta['completed'],'support':meta.get('support'),
                   'resampling_tolerance':RESAMPLING_TOL,'spatial_indicator_limit':.25,'temporal_indicator_limit':.1,'algebra_tolerance':3e-11,
                   'reference_sampling_calls':meta['reference_sampling_calls'],'source_law':'unchanged original signed frequencies, weights and coherent source fibers'},
        'residuals':{'archived801_resampling':meta.get('resampling_residual'),'maximum_flux_and_auxiliary_algebra':algebra,'maximum_kernel_hermiticity':herm},
        'per_case':summaries,'comparison':comparisons,'historical_context_only':historical,
        'runtime':{k:meta[k] for k in ('CPU_seconds','wall_seconds','cpu_budget_seconds','cpu_budget_exceeded','case_CPU_seconds','resampling_CPU_seconds','stop_reason')},
        'scope':{'continuum_error_bound':None,'full_continuous_C0_match':'OPEN','full_CAR_verified':False,'sampled_projector_tangent_reported':True,
                 'stationary_drift_subtracted':False,'new_horizon_scattering_source_state_solve':False,'metric_evolution':False,'physical_constraint_solution':False,
                 'same_pulse_and_source_as_pilot':True,'differences_are_indicators':True},
        'source_hashes':{p:meta['signature'][p] for p in OWNED},'input_hashes':{p:h for p,h in meta['signature'].items() if p not in OWNED},
        'reused_artifacts':meta['reused_artifacts'],'payload':{'path':ARTIFACT,'sha256':digest,'bytes':size},
        'reproducer':'python3 scripts/derive_nsc_retarded_response_resolution.py --check'}


def publish(arrays,meta,oldarrays):
    encoded={**arrays,'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)}
    raw=P.deterministic_npz_bytes(encoded);(ROOT/ARTIFACT).write_bytes(raw)
    record=record_from_arrays(encoded,meta,P.sha256(raw).hexdigest(),len(raw),oldarrays)
    (ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return record


def run():
    start=time.process_time();wall=time.monotonic();before=signature();arrays={};oldarrays={}
    meta={'signature':before,'completed':[],'case_CPU_seconds':{},'cpu_budget_seconds':CPU_CAP,'cpu_budget_exceeded':False,
          'reference_sampling_calls':0,'resampling_CPU_seconds':None,'reused_artifacts':[],'stop_reason':'checkpoint'}
    def stop(signum,frame):raise P.CPUExceeded('targeted resolution CPU budget exhausted')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CPU_CAP)
    def checkpoint(reason):
        meta.update(CPU_seconds=time.process_time()-start,wall_seconds=time.monotonic()-wall,stop_reason=reason)
        if signature()!=before:raise ValueError('resolution owner/input changed during bounded run')
        return publish(arrays,meta,oldarrays)
    try:
        pilot=P.check();oldx,oldphi,data,records=P.inputs();meta['support']=P.support_guard()
        with np.load(ROOT/pilot['payload']['path'],allow_pickle=False) as a:oldarrays={k:a[k] for k in a.files if k.startswith('801_128/')}
        with np.load(ROOT/records[1]['payload']['path'],allow_pickle=False) as a:at_zero=a['14_1/at_zero']
        meta['reused_artifacts']=[pilot['payload'],*[r['payload'] for r in records]]
        arrays.update({'source/'+k:v for k,v in data.items()});arrays['source/at_zero']=at_zero
        x=np.linspace(-2.,1.2,3201);begin=time.process_time();meta['reference_sampling_calls']=1
        phi=characteristic_columns(reference_fields_on_grid(data['energies'],data['mass'].item(),data['angular'].item(),at_zero,x))
        meta['resampling_CPU_seconds']=time.process_time()-begin
        arrays.update({'reference/x':x,'reference/phi':phi})
        meta['resampling_residual']=resampling_error(x,phi,oldx,oldphi)
        if meta['resampling_residual']>=RESAMPLING_TOL:return checkpoint('OPEN: archived reference mismatch; no propagation')
        print('reference resampling matched: '+str(meta['resampling_residual']),flush=True)
        for points,steps in ((1601,128),(3201,128),(3201,256)):
            if steps==256:
                spatial=compare(case_arrays(arrays,'1601_128'),case_arrays(arrays,'3201_128'),'spatial',.25)
                if not all(v['passes'] for v in spatial.values()):return checkpoint('OPEN: spatial indicator failed; temporal solve not authorized')
            case=f'{points}_{steps}';begin=time.process_time()
            result=one_case(points,steps,x,phi,data,meta['support'])
            arrays.update({case+'/'+k:v for k,v in result.items()});meta['completed'].append(case)
            meta['case_CPU_seconds'][case]=time.process_time()-begin
            checkpoint('completed '+case)
            print(case+' checkpoint CPU '+str(meta['case_CPU_seconds'][case]),flush=True)
        meta['stop_reason']='completed all three prescribed cases'
    except P.CPUExceeded:
        meta['cpu_budget_exceeded']=True;meta['stop_reason']='OPEN: CPU budget exhausted; no automatic extension'
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    return checkpoint(meta['stop_reason'])


def check():
    record=json.loads((ROOT/OUTPUT).read_text());P.check()
    if P.digest(ARTIFACT)!=record['payload']['sha256']:raise ValueError('resolution artifact changed')
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if signature()!=meta['signature']:raise ValueError('resolution producers/inputs changed')
    oldx,oldphi,data,records=P.inputs()
    if 'reference/phi' in arrays:
        residual=resampling_error(arrays['reference/x'],arrays['reference/phi'],oldx,oldphi)
        if residual!=meta['resampling_residual']:raise ValueError('reference resampling residual changed')
    for k,v in data.items():
        if not np.array_equal(arrays['source/'+k],v):raise ValueError('original source data changed')
    with np.load(ROOT/records[1]['payload']['path'],allow_pickle=False) as a:
        if not np.array_equal(arrays['source/at_zero'],a['14_1/at_zero']):raise ValueError('reference initial columns changed')
    pilot=json.loads((ROOT/P.OUTPUT).read_text())
    with np.load(ROOT/pilot['payload']['path'],allow_pickle=False) as a:oldarrays={k:a[k] for k in a.files if k.startswith('801_128/')}
    replay=record_from_arrays(arrays,meta,P.digest(ARTIFACT),(ROOT/ARTIFACT).stat().st_size,oldarrays)
    if replay!=record:raise ValueError('resolution receipt differs from saved arrays')
    return record


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    record=run() if args.run else check()
    print(json.dumps({k:record[k] for k in ('status','residuals','comparison','runtime')},indent=2))


if __name__=='__main__':main()
