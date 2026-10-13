#!/usr/bin/env python3
"""Two bounded short-shell solves; authenticated replay never solves an ODE."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np
from recursive_horizons.nsc_frozen_sources import selected_source_context

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_retarded_radial_response import resolution_control as R
import derive_nsc_retarded_fourier_response as F
from recursive_horizons.nsc_retarded_radial_response import archived_amplitudes,pulse_fourier,local_frame_checks,short_response

OUTPUT='results/development/nsc-retarded-radial-response.json'
ARTIFACT='results/development/artifacts/nsc-retarded-radial-response.npz'
HISTORICAL_OWNED=('src/recursive_horizons/nsc_retarded_radial_response.py','scripts/derive_nsc_retarded_radial_response.py',
       'tests/test_nsc_retarded_radial_response.py','docs/nsc-retarded-radial-response.md')
OWNED=('src/recursive_horizons/nsc_retarded_radial_response.py','scripts/derive_nsc_retarded_radial_response.py',
       'tests/test_nsc_retarded_radial_response.py','docs/nsc-retarded-fixed-c0-certificate.md')


def owned_sources():
    return HISTORICAL_OWNED if selected_source_context() is not None else OWNED

FOURIER='results/development/nsc-retarded-fourier-response.json'
TOL=3e-11;CPU_CAP=30.;SETTINGS=((2e-10,2e-12),(2e-13,2e-15))


def signature():return {**R.signature(),**F.signature(),**{p:R.P.digest(p) for p in (*owned_sources(),R.OUTPUT,FOURIER)}}
def unpack(value):return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def inputs():
    record=R.check();fourier=json.loads((ROOT/FOURIER).read_text())
    if not record['status'].startswith('PASS:') or not fourier['status'].startswith('PASS:'):
        raise ValueError('completed source resolution and Fourier indicators required')
    for p,h in fourier['source_sha256'].items():
        if R.P.digest(p)!=h:raise ValueError('Fourier owner/input changed: '+p)
    if fourier['input_payload']!=record['payload']:raise ValueError('Fourier/reference payload binding differs')
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as a:
        data={k:a[k] for k in a.files if k.startswith(('source/','reference/'))}
    return record,fourier,data


def pair_comparison(actual,reference):
    rows=[]
    for o in range(4):
        for i in range(4):
            error=F.maxabs(actual[o,i]-reference[o,i]);norm=F.maxabs(reference[o,i]);ratio=error/norm if norm else None
            rows.append({'output_index':o,'input_index':i,'difference':error,'reference_norm':norm,
                         'relative_difference':ratio,'passes':ratio is not None and ratio<.01})
    return {'maximum_difference':F.maxabs(actual-reference),'relative_max_difference':F.maxabs(actual-reference)/F.maxabs(reference),
            'relative_indicator_limit':.01,'every_pair_passes':all(v['passes'] for v in rows),'pairs':rows}


def analyze(arrays,meta,fourier):
    rows={};f=arrays.get('reference/incoming');C=arrays.get('source/source');P=arrays.get('source/projector')
    Bfine=None;Dfine=None;checks=dict(meta.get('local_checks',{}))
    for label in meta['completed']:
        A=arrays[label+'/homogeneous'];B=arrays[label+'/response'];D=F.covariance_pairs(B,f,C);CAR=F.covariance_pairs(B,f,P)
        K=np.array([B[i,i]@f[i].conj().T for i in range(4)])
        Ksolve=np.array([np.linalg.solve(f[i,:,:2].T,B[i,i,:,:2].T).T for i in range(4)])
        row={'homogeneous_endpoint_residual':F.maxabs(A-f),'mode_response_max':F.maxabs(B),'covariance_pair_max':F.maxabs(D),
             'CAR_pair_residual':F.maxabs(CAR),'pair_hermiticity_residual':F.maxabs(D-D.transpose(1,0,3,2).conj()),
             'diagonal_K_antiHermitian':F.maxabs(K+K.swapaxes(-1,-2).conj()),
             'diagonal_K_inverse_antiHermitian':F.maxabs(Ksolve+Ksolve.swapaxes(-1,-2).conj()),
             'diagonal_K_definition_difference':F.maxabs(K-Ksolve),
             'source_coisometry_residual':F.maxabs(f@f.swapaxes(-1,-2).conj()-np.eye(2)),
             'closed_response_column':F.maxabs(B[:,:,:,2]),'solver':meta['solver'][label]}
        rows[label]=row
        if label=='fine':
            Bfine,Dfine=B,D
            checks.update({k:row[k] for k in ('homogeneous_endpoint_residual','CAR_pair_residual','pair_hermiticity_residual','closed_response_column')})
    comparison={};positive={};refinement={}
    if Bfine is not None:
        archived=fourier['per_case']['3201_256']
        comparison={'mode':pair_comparison(Bfine,unpack(archived['mode_pair_response'])),
                    'covariance':pair_comparison(Dfine,unpack(archived['covariance_pairs']))}
        if 'coarse' in meta['completed']:
            Bcoarse=arrays['coarse/response'];Dcoarse=F.covariance_pairs(Bcoarse,f,C)
            refinement={'mode':F.maxabs(Bfine-Bcoarse),'covariance':F.maxabs(Dfine-Dcoarse),
                        'homogeneous':F.maxabs(arrays['fine/homogeneous']-arrays['coarse/homogeneous'])}
            checks.update({'refinement_'+k:v for k,v in refinement.items()})
        indices=np.flatnonzero(arrays['source/energies']>0);Pv=np.broadcast_to(np.diag([0.,1.,0.]),C.shape)
        vacuum=F.covariance_pairs(Bfine,f,Pv);remainder=F.covariance_pairs(Bfine,f,C-Pv)
        select=lambda value:value[indices[:,None],indices[None,:]]
        positive={'energy_indices':indices.tolist(),'full_pairs':F.pack(select(Dfine)),'Pv_pairs':F.pack(select(vacuum)),
                  'retained_remainder_pairs':F.pack(select(remainder)),
                  'decomposition_residual':F.maxabs(Dfine-vacuum-remainder),
                  'full_max':F.maxabs(select(Dfine)),'Pv_max':F.maxabs(select(vacuum)),
                  'retained_remainder_max':F.maxabs(select(remainder)),
                  'phase_independent_physical_bound':None}
    passed=(meta['completed']==['coarse','fine'] and max(checks.values(),default=float('inf'))<TOL
            and all(v['every_pair_passes'] for v in comparison.values()) and meta.get('Fourier_rule_difference',float('inf'))<TOL
            and not meta['cpu_budget_exceeded'])
    return rows,checks,comparison,refinement,positive,passed


def record_from_arrays(arrays,meta,digest,size,fourier):
    rows,checks,comparison,refinement,positive,passed=analyze(arrays,meta,fourier)
    return {'schema':'NSC-RETARDED-RADIAL-RESPONSE-v1','accountable_author':'Douglas Ek',
        'status':'PASS: conditional short-shell numerical indicators; physical preparation OPEN' if passed else 'OPEN: short-shell numerical indicators unresolved',
        'control':{'family':'14_1','energies':arrays.get('source/energies',np.array([])).tolist(),'upstream':meta.get('upstream'),
                   'upstream_index':meta.get('upstream_index'),'support':meta.get('support'),'completed':meta['completed'],
                   'algebra_and_refinement_tolerance':TOL,'pair_relative_indicator_limit':.01,
                   'Fourier_orders':[96,192],'Fourier_rule_difference':meta.get('Fourier_rule_difference'),
                   'source_weights_used':False,'angular_multiplicity_used':False},
        'per_solve':rows,'residuals':checks,'comparison_to_saved_PDE_fourier':comparison,'radial_refinement':refinement,
        'positive_energy_decomposition':positive,'runtime':{k:meta[k] for k in ('CPU_seconds','wall_seconds','cpu_budget_seconds','cpu_budget_exceeded','stop_reason')},
        'scope':{'physical_low_energy_preparation':'OPEN','rigorous_response_error_bound':None,'full_C0_exclusion':'OPEN',
                 'global_PDE_horizon_scattering_solve':False,'source_state_replaced':False,'metric_evolution':False,
                 'stationarity_or_constraint_solution':False,'tolerance_changes_are_indicators':True,'CAR_residual_subtracted':False},
        'source_hashes':{p:meta['signature'][p] for p in owned_sources()},'input_hashes':{p:h for p,h in meta['signature'].items() if p not in owned_sources()},
        'input_payload':meta.get('input_payload'),'payload':{'path':ARTIFACT,'sha256':digest,'bytes':size},
        'reproducer':'python3 scripts/derive_nsc_retarded_radial_response.py --check'}


def publish(arrays,meta,fourier):
    content={**arrays,'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)}
    raw=R.P.deterministic_npz_bytes(content);(ROOT/ARTIFACT).write_bytes(raw)
    record=record_from_arrays(content,meta,R.P.sha256(raw).hexdigest(),len(raw),fourier)
    (ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return record


def run():
    start=time.process_time();wall=time.monotonic();before=signature();arrays={};fourier={}
    meta={'signature':before,'completed':[],'solver':{},'cpu_budget_seconds':CPU_CAP,'cpu_budget_exceeded':False,'stop_reason':'checkpoint'}
    def stop(signum,frame):raise R.P.CPUExceeded('radial response CPU budget exhausted')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CPU_CAP)
    def checkpoint(reason):
        meta.update(CPU_seconds=time.process_time()-start,wall_seconds=time.monotonic()-wall,stop_reason=reason)
        if signature()!=before:raise ValueError('radial producer/input changed during execution')
        return publish(arrays,meta,fourier)
    try:
        parent,fourier,data=inputs();x=data['reference/x'];phi=data['reference/phi'];E=data['source/energies']
        index=int(np.flatnonzero(x>=1.03)[0]);rho=float(x[index]);end=int(np.flatnonzero(x==1.)[0])
        guard=parent['control']['support']
        if not rho>guard['radial_support'][1]:raise ValueError('upstream node is not above the entire pulse')
        meta.update(upstream=rho,upstream_index=index,support=guard,input_payload=parent['payload'])
        arrays.update({k:v for k,v in data.items() if k.startswith('source/')})
        A0=archived_amplitudes(x,phi,E,index);incoming=archived_amplitudes(x,phi,E,end)
        arrays.update({'reference/upstream':A0,'reference/incoming':incoming})
        mass=data['source/mass'].item();angular=data['source/angular'].item()
        meta['local_checks']=local_frame_checks((1.,(1+rho)/2,rho),E,mass,angular)
        for order in (96,192):
            arrays['Fourier/'+str(order)]=pulse_fourier(E[:,None]-E[None,:],guard['axial_center'],guard['axial_halfwidth'],order)
        meta['Fourier_rule_difference']=F.maxabs(arrays['Fourier/96']-arrays['Fourier/192'])
        if max(meta['local_checks'].values())>=TOL:return checkpoint('OPEN: local frame algebra failed; no solves')
        for label,(rtol,atol) in zip(('coarse','fine'),SETTINGS,strict=True):
            begin=time.process_time()
            A,B,info=short_response(rho,E,mass,angular,A0,arrays['Fourier/192'],rtol=rtol,atol=atol)
            info['CPU_seconds']=time.process_time()-begin;meta['solver'][label]=info
            arrays.update({label+'/homogeneous':A,label+'/response':B});meta['completed'].append(label)
            checkpoint('completed '+label)
            print(label+' short solve CPU '+str(info['CPU_seconds']),flush=True)
        meta['stop_reason']='completed exactly two prescribed short-shell solves'
    except R.P.CPUExceeded:
        meta['cpu_budget_exceeded']=True;meta['stop_reason']='OPEN: CPU cap reached; no automatic extension'
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    return checkpoint(meta['stop_reason'])


def archived_map_replay_envelope(x, phi, energies, index):
    """Propagate the existing chart tolerance and binary64 map arithmetic.

    chart_coordinates uses epsabs=epsrel=2e-13 for S.  The phase is Lipschitz:
    |exp(i E S1)-exp(i E S2)| <= |E| |S1-S2|.  Two replays therefore charge
    twice that fixed requested quadrature accuracy.  A conservative gamma_128
    covers the small coframe/dot chain, with its geometric conditioning and
    uncancelled absolute-product scale.  This is replay accuracy, not a
    certified quadrature or physical source-error bound.
    """
    from recursive_horizons.nsc_retarded_radial_response import canonical_map
    from recursive_horizons.nsc_transmitting_dirac_domain import MODE_TO_CURRENT
    from recursive_horizons.nsc_lorentzian import geometry
    from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
    rho = float(x[index]); beta = float(geometry(rho)[0])
    scale = np.max(abs(canonical_map(rho)) @ abs(MODE_TO_CURRENT) @ abs(phi[[index, len(x)+index]]))
    eps = np.finfo(np.float64).eps
    gamma = 128 * eps / (1 - 128 * eps)
    geometric_condition = 1 + beta**2/(beta**2-1) + 1/(beta*(1-1/beta**2))
    S = chart_coordinates(rho)[1]
    coordinate_error = max(2e-13, 2e-13 * abs(S))
    phase_error = float(np.max(abs(energies))) * coordinate_error
    bound = 2 * float(scale) * (1+gamma) * (phase_error + gamma*geometric_condition)
    return float(np.nextafter(bound, np.inf)) if bound else 0.


def inverse_diagnostic_replay_envelopes(arrays, label):
    """2x2 partial-pivot LU forward error from its backward-error model.

    Partial pivoting has growth <=2 for a 2x2 matrix. gamma_(8*(3n+1))
    conservatively counts complex LU/triangular arithmetic. Propagate that
    through cond_inf(A), provided the resulting relative error is <1, and
    bound two implementations. Anti-Hermitian sums and the direct-dot
    comparison add their own operand-scale arithmetic terms.
    """
    f=arrays['reference/incoming']; B=arrays[label+'/response']
    eps=np.finfo(np.float64).eps; n=2
    gamma=8*(3*n+1)*eps/(1-8*(3*n+1)*eps)
    solve_error=direct_error=scale=0.
    for i in range(4):
        A=f[i,:,:2].T; rhs=B[i,i,:,:2].T
        solution=np.linalg.solve(A,rhs)
        condition=float(np.real(np.linalg.cond(A,np.inf)))
        relative=2*gamma*condition
        if not np.isfinite(relative) or relative>=1:
            raise ValueError('ill-conditioned canonical input replay')
        norm=float(np.linalg.norm(solution,np.inf))
        solve_error=max(solve_error,relative/(1-relative)*norm)
        magnitude=float(np.max(abs(B[i,i]) @ abs(f[i].conj().T)))
        direct_error=max(direct_error,gamma*magnitude)
        scale=max(scale,norm,magnitude)
    return {
        'diagonal_K_inverse_antiHermitian':float(np.nextafter(4*solve_error+16*eps*scale,np.inf)),
        'diagonal_K_definition_difference':float(np.nextafter(2*(solve_error+direct_error)+16*eps*scale,np.inf)),
    }


def validate_replay_record(old, record, arrays):
    from copy import deepcopy
    from recursive_horizons.nsc_weighted_receipt_replay import exact_record_equal
    strict=deepcopy(record)
    for label in old['per_solve']:
        for field,bound in inverse_diagnostic_replay_envelopes(arrays,label).items():
            stored=old['per_solve'][label][field]; computed=record['per_solve'][label][field]
            if not all(type(v) is float and np.isfinite(v) and v>=0 for v in (stored,computed)) or abs(stored-computed)>bound:
                raise ValueError('radial inverse diagnostic exceeds replay arithmetic envelope')
            strict['per_solve'][label][field]=stored
    if not exact_record_equal(strict,old):raise ValueError('radial saved-array receipt differs')


def validate_archived_map(expected, saved, x, phi, energies, index):
    bound=archived_map_replay_envelope(x,phi,energies,index)
    if (not isinstance(expected,np.ndarray) or not isinstance(saved,np.ndarray)
            or expected.dtype!=np.dtype(complex) or saved.dtype!=expected.dtype
            or expected.shape!=saved.shape or not np.isfinite(expected).all() or not np.isfinite(saved).all() or F.maxabs(expected-saved)>bound):
        raise ValueError('canonical archived input map exceeds replay arithmetic envelope')


def check():
    old=json.loads((ROOT/OUTPUT).read_text());parent,fourier,data=inputs()
    if R.P.digest(ARTIFACT)!=old['payload']['sha256']:raise ValueError('radial artifact changed')
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if meta['signature']!=signature() or meta['input_payload']!=parent['payload']:raise ValueError('radial producer/input binding changed')
    for k,v in data.items():
        if k.startswith('source/') and not np.array_equal(arrays[k],v):raise ValueError('radial source input changed')
    for key,index in (('upstream',meta['upstream_index']),('incoming',int(np.flatnonzero(data['reference/x']==1.)[0]))):
        expected=archived_amplitudes(data['reference/x'],data['reference/phi'],data['source/energies'],index)
        validate_archived_map(expected,arrays['reference/'+key],data['reference/x'],data['reference/phi'],data['source/energies'],index)
    record=record_from_arrays(arrays,meta,R.P.digest(ARTIFACT),(ROOT/ARTIFACT).stat().st_size,fourier)
    validate_replay_record(old,record,arrays)
    return old


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    record=run() if args.run else check()
    print(json.dumps({k:record[k] for k in ('status','residuals','radial_refinement','runtime')},indent=2))


if __name__=='__main__':main()
