#!/usr/bin/env python3
"""One bounded retarded prepared response; replay never propagates fields."""
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
from scipy.linalg import block_diag
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_compatible_history_geometry import CompatibleRadiusDirection,CompatibleIncomingMetric
from recursive_horizons.nsc_prepared_history_jets import evolve_prepared_field_jets,restriction_covariance_tangent
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator
from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain,MODE_TO_CURRENT
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_transmitting_resolvent import profile
from recursive_horizons.nsc_lorentzian import geometry
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-retarded-compatible-response.json'
ARTIFACT='results/development/artifacts/nsc-retarded-compatible-response.npz'
OWNED=('scripts/derive_nsc_retarded_compatible_response.py','tests/test_nsc_retarded_compatible_response.py','docs/nsc-retarded-compatible-response.md')
RECORDS=('results/development/nsc-transmitting-history-modes.json','results/development/nsc-pg-ctp-mode-jets.json','results/development/nsc-compatible-prepared-history.json')
OWNERS=tuple('src/recursive_horizons/'+n+'.py' for n in ('nsc_compatible_history_geometry','nsc_prepared_history_jets',
    'nsc_transmitting_history_jets','nsc_transmitting_history_modes','nsc_transmitting_dirac_domain','nsc_ks_spacetime_variation',
    'nsc_pg_ks_metric_pullback','nsc_transmitting_resolvent','nsc_lorentzian','nsc_mode_resolved_cauchy_state'))+('docs/nsc-retarded-compatible-preparation.md',)
CASES=((401,64),(801,64),(801,128));CPU_CAP=120.
MAXABS=lambda x:float(np.max(abs(x)))


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signature():return {p:digest(p) for p in (*OWNED,*RECORDS,*OWNERS)}

def inputs():
    records=[json.loads((ROOT/p).read_text()) for p in RECORDS]
    for record in records[:2]:
        for section in ('source_hashes','input_hashes'):
            for p,h in record[section].items():
                if digest(p)!=h:raise ValueError('retarded input owner changed: '+p)
        if digest(record['payload']['path'])!=record['payload']['sha256']:raise ValueError('retarded input artifact changed')
    for p,h in records[2]['source_sha256'].items():
        if digest(p)!=h:raise ValueError('prepared adapter receipt changed: '+p)
    with np.load(ROOT/records[0]['payload']['path'],allow_pickle=False) as a:
        x=a['14_1/x'];phi=a['14_1/reference_field'];meta=json.loads(a['metadata_json'].tobytes())
    if meta['mode_payload']!=records[1]['payload']:raise ValueError('reference/source artifact binding differs')
    with np.load(ROOT/records[1]['payload']['path'],allow_pickle=False) as a:
        data={k:a['14_1/'+k] for k in ('energies','weights','source','projector','mass','angular')}
    if x.shape!=(801,) or phi.shape!=(1602,12) or not np.array_equal(np.flatnonzero(x==1.),[750]):
        raise ValueError('authenticated group14 grid/columns changed')
    if data['energies'].shape!=(4,) or np.any(data['weights']<=0) or np.any(abs(data['energies'])>=data['mass'].item()):
        raise ValueError('original four real signed below-threshold controls required')
    if not np.array_equal(data['projector'],np.broadcast_to(np.diag([1.,1.,0.]),(4,3,3))):
        raise ValueError('closed infinity source projector changed')
    return x,phi,data,records[:2]


def axial_bump(center,halfwidth=.06):
    """Existing normalized profile and coherent derivatives through order3."""
    def value(z,order):
        if order not in (0,1,2,3):raise ValueError('owned incoming derivative order0..3 required')
        u=(z-center)/halfwidth
        if abs(u)>=1:return 0.
        b,bp=profile(u)
        if order==0:return b
        if order==1:return bp/halfwidth
        d=1-u*u;g1=-2*u/d**2;g2=-2/d**2-8*u*u/d**3
        if order==2:return b*(g1*g1+g2)/halfwidth**2
        g3=-24*u/d**3-48*u**3/d**4
        return b*(g1**3+3*g1*g2+g3)/halfwidth**3
    return value


def support_guard():
    T1,S1=chart_coordinates(1.)
    lower=brentq(lambda rho:chart_coordinates(rho)[0]-T1-.03,.8,1.,xtol=5e-15)
    upper=brentq(lambda rho:chart_coordinates(rho)[0]-T1+.03,1.,1.15,xtol=5e-15)
    start=S1+.09-chart_coordinates(upper)[1];stop=S1+.21-chart_coordinates(lower)[1]
    if not (-2<lower<1<upper<1.2 and 0<start<stop<.3):raise ValueError('retarded support does not fit supplied domain/time')
    # Pure radius directions retain these reference principal velocities.
    speeds=np.array([1.,-1.])[:,None]-geometry(np.array([-2.,lower,1.,upper,1.2]))[0]
    if np.max(speeds)>=0:raise ValueError('strict left-going reference characteristics required')
    return {'radial_support':[float(lower),float(upper)],'PG_time_support':[float(start),float(stop)],
        'normal_inner':.007,'normal_outer':.03,'axial_center':float(S1+.15),'axial_halfwidth':.06,
        'T1':float(T1),'S1':float(S1),'initial_time':0.,'final_time':.3,
        'initial_support_margin':float(start),'final_support_margin':float(.3-stop),
        'right_support_margin':float(1.2-upper),'maximum_checked_characteristic_speed':float(np.max(speeds)),
        'guard_kind':'numerical support bounds with strict margins; analytic monotonicity and fixed-radius-only principal speeds'}


class CachedZeroRadiusProvider:
    """Cache fixed grid/chart only; use the owned radius basis at every time."""
    def __init__(self,x,direction):
        self.x=np.asarray(x);self.direction=direction
        self.original=CompatibleIncomingMetric((0.,),(direction,))
        self.metric=self.original.values(0.,self.x)
        coordinates=np.array([chart_coordinates(float(rho)) for rho in self.x])
        self.normal=coordinates[:,0]-chart_coordinates(1.)[0];self.shift=coordinates[:,1]
    def values(self,time,x):
        if not np.array_equal(x,self.x):raise ValueError('cached grid changed')
        return self.metric.copy()
    def log_directions(self,time,x,metric):
        if not np.array_equal(x,self.x) or not np.array_equal(metric,self.metric):raise ValueError('cached zero metric/grid changed')
        d=np.zeros((1,4,len(x)));d[0,3]=self.direction.value(self.normal,time+self.shift)/metric[3]
        return d


def one_case(points,steps,x801,phi801,data,guard):
    stride=800//(points-1);x=x801[::stride]
    phi=phi801.reshape(2,801,12)[:,::stride].reshape(2*points,12)
    owner=FourthOrderModePropagator(x,data['mass'].item(),data['angular'].item())
    direction=CompatibleRadiusDirection(axial_bump(guard['axial_center']),lambda z,n:0.,.007,.03)
    provider=CachedZeroRadiusProvider(x,direction)
    times=np.linspace(0.,.3,steps+1);E=np.repeat(data['energies'],3)
    initial_tangents=np.zeros((1,*phi.shape),complex)
    right=np.array([points-1,2*points-1]);N,beta,q,r=owner.reference
    speed=np.array([N[-1]/q[-1]-beta[-1],-N[-1]/q[-1]-beta[-1]])
    sat=-speed/owner.weights[-1]
    def incoming(t):
        forcing=np.zeros_like(phi);forcing[right]=sat[:,None]*phi[right]*np.exp(-1j*E*t)[None,:]
        return forcing,np.zeros_like(initial_tangents)
    ix=np.flatnonzero(x==1.)
    if len(ix)!=1:raise ValueError('exact rho1 node required')
    sample_rows=np.array([ix[0],points+ix[0]])
    run=evolve_prepared_field_jets(owner,times,provider,phi,initial_tangents,incoming,sample_rows=sample_rows)
    domain=TransmittingDiracSeamDomain(N[ix[0]],q[ix[0]],beta[ix[0]],r[ix[0]],surface_id='fixed-rho1-retarded-control')
    _,inverse,_=domain.trace_map(np.ones(1));restriction=inverse[0]@MODE_TO_CURRENT/(r[ix[0]]*np.sqrt(q[ix[0]]))
    F=np.einsum('ab,tbc->tac',restriction,run['sampled_fields'])
    dF=np.einsum('ab,tbc->tac',restriction,run['sampled_tangent_fields'][0])
    analytic=phi[sample_rows][None,:,:]*np.exp(-1j*times[:,None,None]*E[None,None,:])
    stationary=np.einsum('ab,tbc->tac',restriction,analytic)
    rootw=np.repeat(np.sqrt(data['weights']/(2*np.pi)),3);C=block_diag(*data['source'])
    weightedF=(F*rootw).reshape(-1,12);weighteddF=(dF*rootw).reshape(-1,12)
    dkernel=restriction_covariance_tangent(weightedF,weighteddF,C,np.zeros_like(C))
    P=block_diag(*data['projector'])
    dprojector=restriction_covariance_tangent(weightedF,weighteddF,P,np.zeros_like(P))
    diagnostics={k:run[k] for k in ('norm_flux_residual','metric_flux_tangent_residual','auxiliary_identity_residual','tangent_boundary_trace')}
    diagnostics.update(initial_tangent_norm=MAXABS(initial_tangents),incident_tangent_norm=0.,source_covariance_tangent_norm=0.,
        initial_metric_change=MAXABS(provider.metric-owner.reference),
        source_closed_column_norm=MAXABS(phi[:,2::3]),source_covariance_closed_column_norm=MAXABS(data['source'][:,:,2]),
        source_covariance_hermiticity=MAXABS(C-C.conj().T),
        maximum_characteristic_speed=float(np.max(1.-beta)),
        frame_map_parameter_derivative=0.,stationary_phase_reapplied=False)
    return {'x':x,'times':times,'sample_rows':sample_rows,'field':run['evolved_field'],'field_tangent':run['tangent_fields'][0],
        'trace':F,'trace_tangent':dF,'stationary_trace':stationary,
        'stationary_field':phi*np.exp(-1j*E*.3)[None,:],
        'covariance_tangent':dkernel,'projector_tangent':dprojector,'diagnostics_json':np.frombuffer(json.dumps(diagnostics,sort_keys=True).encode(),np.uint8)}


def frame_constructor_correction_identity():
    """Fixed-rho1 map equality, not equality of the misordered orbit metrics."""
    N,q,r=1.,1.,float(np.sqrt(2.));beta=float(geometry(1.)[0])
    legacy=TransmittingDiracSeamDomain(N,beta,q,r)
    corrected=TransmittingDiracSeamDomain(N,q,beta,r)
    old=legacy.trace_map(np.ones(1));new=corrected.trace_map(np.ones(1))
    equality=all(np.array_equal(a,b) for a,b in zip(old,new))
    oldR=old[1][0]@MODE_TO_CURRENT/(r*np.sqrt(q))
    newR=new[1][0]@MODE_TO_CURRENT/(r*np.sqrt(q))
    if not equality or not np.array_equal(oldR,newR):
        raise ArithmeticError('constructor correction changes the executed restriction; cannot reuse output')
    return {'trace_map_triples_bitwise_equal':True,'restriction_map_bitwise_equal':True,
        'restriction_residual':MAXABS(oldR-newR),
        'reason':'the used trace/inverse/Gram depend on N,r,q*beta, unchanged by this swap; the external r*sqrt(q) denominator was already correct',
        'coordinate_normal_or_metric_equality_claimed':False}


def analyze(arrays,meta):
    summaries={}
    for case in meta['completed']:
        get=lambda key:arrays[case+'/'+key]
        norms={name:MAXABS(get(key)) for name,key in (('field','field_tangent'),('trace','trace_tangent'),('covariance','covariance_tangent'))}
        summaries[case]={'response_max_norms':norms,
            'covariance_hermiticity_residual':MAXABS(get('covariance_tangent')-get('covariance_tangent').conj().T),
            'sampled_projector_kernel_tangent_norm':MAXABS(get('projector_tangent')),
            'sampled_projector_kernel_hermiticity_residual':MAXABS(get('projector_tangent')-get('projector_tangent').conj().T),
            'stationary_reference_field_drift':MAXABS(get('field')-get('stationary_field')),
            'stationary_reference_trace_drift':MAXABS(get('trace')-get('stationary_trace')),
            'diagnostics':json.loads(get('diagnostics_json').tobytes())}
    comparisons={};passed=len(meta['completed'])==3
    if passed:
        for label,coarse,fine,limit in [('spatial','401_64','801_64',.25),('temporal','801_64','801_128',.10)]:
            rows={}
            for name,key in [('field','field_tangent'),('trace','trace_tangent'),('covariance','covariance_tangent')]:
                c=arrays[coarse+'/'+key];f=arrays[fine+'/'+key]
                if label=='spatial' and name=='field':f=f.reshape(2,801,12)[:,::2].reshape(802,12)
                if label=='temporal' and name=='trace':f=f[::2]
                if label=='temporal' and name=='covariance':f=f.reshape(129,2,129,2)[::2,:,::2,:].reshape(130,130)
                difference=MAXABS(c-f);norm=MAXABS(f);ratio=difference/norm if norm else None
                ok=ratio is not None and ratio<limit
                rows[name]={'max_absolute_difference':difference,'fine_compared_response_norm':norm,
                    'difference_over_fine_norm':ratio,'relative_indicator_limit':limit,'passes':ok}
                passed=passed and ok
            comparisons[label]=rows
    algebra=max((r['diagnostics'][k] for r in summaries.values() for k in ('norm_flux_residual','metric_flux_tangent_residual','auxiliary_identity_residual')),default=0.)
    herm=max((r['covariance_hermiticity_residual'] for r in summaries.values()),default=0.)
    passed=passed and algebra<3e-11 and herm<3e-11 and not meta['cpu_budget_exceeded']
    return summaries,comparisons,passed,algebra,herm


def publish(arrays,meta):
    current=signature()
    meta.setdefault('execution_signature',dict(current))
    meta['verification_signature']=dict(current)
    meta['signature']=current
    arrays={**arrays,'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)}
    raw=deterministic_npz_bytes(arrays);path=ROOT/ARTIFACT;path.write_bytes(raw)
    record=record_from_arrays(arrays,meta,sha256(raw).hexdigest(),len(raw))
    (ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return record


def record_from_arrays(arrays,meta,artifact_hash,artifact_bytes):
    summaries,comparisons,passed,algebra,herm=analyze(arrays,meta)
    return {'schema':'NSC-RETARDED-COMPATIBLE-RESPONSE-v1','accountable_author':'Douglas Ek',
        'status':'PASS: sampled retarded response indicators only' if passed else 'OPEN: bounded sampled response unresolved or checkpointed',
        'control':{'family':'14_1','amplitude':0.,'cases_requested':[list(x) for x in CASES],'cases_completed':meta['completed'],
            'support':meta['support'],'weights':'sqrt(w/(2*pi)) on each source column; no multiplicity factor added',
            'source_law':'unchanged original coherent source blocks; all four infinity columns closed',
            'initial_and_incident_tangents':'explicit0 by unchanged entire past/upstream support',
            'response':'actual rho1 two-spin time trace, z=tau+S1; no stationary phase or interpolation',
            'spatial_indicator_limit':.25,'temporal_indicator_limit':.10,'algebra_tolerance':3e-11},
        'per_case':summaries,'comparison':comparisons,
        'residuals':{'maximum_flux_and_auxiliary_algebra':algebra,'maximum_covariance_hermiticity':herm},
        'runtime':{k:meta[k] for k in ('CPU_seconds','wall_seconds','cpu_budget_seconds','cpu_budget_exceeded','case_CPU_seconds')},
        'scope':{'causal_prepared_off_shell_family':True,'full_continuous_C0_match':'OPEN',
            'four_energy_contribution_decides_C0':False,'continuum_error_bound':None,'differences_are_indicators':True,
            'complete_canonical_CAR_kernel_tangent':'zero on fixed intrinsic Sigma after full source resolution',
            'sampled_projector_tangent':'reported without enforcing zero or subtracting a correction',
            'parent_or_constraint_solution_selected':False,'stress_or_metric_equations_solved':False,
            'old_radial_source_or_field_producers_rerun':False,'new_source_law_or_Gamma_term':False,
            'metric_timestep':False,'publication':False},
        'source_hashes':{p:meta['verification_signature'][p] for p in OWNED},
        'verification_source_hashes':{p:meta['verification_signature'][p] for p in OWNED},
        'execution_source_hashes':{p:meta['execution_signature'][p] for p in OWNED},
        'postprocessing_correction':meta.get('postprocessing_correction'),
        'input_hashes':{p:meta['signature'][p] for p in (*RECORDS,*OWNERS)},
        'reused_artifacts':meta['reused_artifacts'],
        'payload':{'path':ARTIFACT,'sha256':artifact_hash,'bytes':artifact_bytes},
        'reproducer':'python3 scripts/derive_nsc_retarded_compatible_response.py --check'}


class CPUExceeded(RuntimeError):pass

def run():
    start_cpu=time.process_time();start_wall=time.monotonic();before=signature()
    arrays={};meta={'completed':[],'case_CPU_seconds':{},'cpu_budget_seconds':CPU_CAP,'cpu_budget_exceeded':False}
    def stop(signum,frame):raise CPUExceeded('120 second producer CPU budget exhausted')
    old_handler=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CPU_CAP)
    try:
        x,phi,data,records=inputs();meta['support']=support_guard();meta['reused_artifacts']=[r['payload'] for r in records]
        arrays.update({'source/'+key:value for key,value in data.items()})
        print('support checked; three prescribed solves only',flush=True)
        for points,steps in CASES:
            case=f'{points}_{steps}';begin=time.process_time()
            result=one_case(points,steps,x,phi,data,meta['support'])
            arrays.update({case+'/'+k:v for k,v in result.items()});meta['completed'].append(case)
            meta['case_CPU_seconds'][case]=time.process_time()-begin
            meta.update(CPU_seconds=time.process_time()-start_cpu,wall_seconds=time.monotonic()-start_wall)
            if signature()!=before:raise ValueError('input/producer changed during bounded pilot')
            publish(arrays,meta)
            print(case+' checkpoint, CPU '+str(meta['case_CPU_seconds'][case]),flush=True)
    except CPUExceeded:
        meta['cpu_budget_exceeded']=True
        if 'support' not in meta:raise
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,old_handler)
    meta.update(CPU_seconds=time.process_time()-start_cpu,wall_seconds=time.monotonic()-start_wall)
    if signature()!=before:raise ValueError('input/producer changed during bounded pilot')
    return publish(arrays,meta)


def check():
    old=json.loads((ROOT/OUTPUT).read_text());path=ROOT/old['payload']['path']
    if digest(old['payload']['path'])!=old['payload']['sha256']:raise ValueError('response artifact changed')
    with np.load(path,allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if meta['verification_signature']!=signature() or meta['signature']!=meta['verification_signature']:
        raise ValueError('response verification producer/input changed')
    correction=meta.get('postprocessing_correction')
    if correction:
        if correction['map_identity']!=frame_constructor_correction_identity():raise ValueError('map correction identity changed')
        for key,h in correction['unchanged_array_sha256'].items():
            if sha256(arrays[key].tobytes()).hexdigest()!=h:raise ValueError('executed numerical array changed: '+key)
        for p in OWNED:
            if sha256(arrays['executed_sources/'+p].tobytes()).hexdigest()!=meta['execution_signature'][p]:
                raise ValueError('executed source provenance changed: '+p)
    record=record_from_arrays(arrays,meta,digest(old['payload']['path']),path.stat().st_size)
    if record!=old:raise ValueError('response receipt differs from replay')
    return record


def main():
    p=argparse.ArgumentParser();mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run',action='store_true');mode.add_argument('--check',action='store_true');args=p.parse_args()
    record=run() if args.run else check()
    print(json.dumps({k:record[k] for k in ('status','per_case','comparison','runtime')},indent=2))

if __name__=='__main__':main()
