#!/usr/bin/env python3
"""One bounded source-phase accuracy diagnostic; check only replays arrays."""
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
from scipy.interpolate import make_interp_spline
from scipy.sparse import bmat,csr_matrix,diags
from scipy.sparse.linalg import expm_multiply

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_evolved_incoming_state as STATE
from derive_nsc_incoming_source_update_v3 import authenticate
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation,rho1_restriction_map
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter,column_derivatives
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-exact-source-phase.json'
ARTIFACT='results/development/artifacts/nsc-exact-source-phase.npz'
SOURCES=('scripts/derive_nsc_exact_source_phase.py','docs/nsc-exact-source-phase.md')
INPUTS=('results/development/nsc-evolved-incoming-state.json',
        'results/development/nsc-evolved-incoming-constraints.json',
        'results/development/nsc-incoming-fixed-transfer.json',
        'src/recursive_horizons/nsc_transmitting_history_jets.py',
        'src/recursive_horizons/nsc_transmitting_history_modes.py',
        'src/recursive_horizons/nsc_evolved_incoming_state.py',
        'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
        'scripts/derive_nsc_evolved_incoming_state.py',
        'scripts/derive_nsc_incoming_source_update_v3.py')
MAX=lambda x:float(np.max(abs(x),initial=0.))


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*SOURCES,*INPUTS)}


def inputs():
    records=[json.loads((ROOT/p).read_text()) for p in INPUTS[:3]]
    for record in records:authenticate(record)
    return records


def analyze(a,meta):
    F,analytic=a['columns'],a['reference_columns']
    C,w=a['source_covariance'],a['column_weights']
    zero=np.zeros((0,*F.shape),complex)
    params=meta['parameters']
    Dz=make_interp_spline(a['z'],np.eye(len(a['z'])),k=3).derivative()(a['z'])
    interp_z,_=column_derivatives(F,zero,Dz)
    value=source_column_matter(F*w,a['pde_axial_columns']*w,C,zero,zero,**params)['action_gradient']
    interp=source_column_matter(F*w,interp_z*w,C,zero,zero,**params)['action_gradient']
    reference=source_column_matter(analytic*w,a['reference_axial_columns']*w,C,zero,zero,**params)['action_gradient']
    pde_error=np.max(abs(value-reference),axis=0)
    interp_error=np.max(abs(interp-reference),axis=0)
    old=np.asarray(meta['old_midpoint_forcing_drift'])
    ratios=pde_error/old
    return {'exact_harmonic_phase_residual':MAX(a['phase']-a['expected_phase']),
        'incoming_column_drift':MAX(F-analytic),'pde_axial_column_drift':MAX(a['pde_axial_columns']-a['reference_axial_columns']),
        'PDE_matter_drift':pde_error.tolist(),'interpolated_trace_matter_drift':interp_error.tolist(),
        'PDE_over_old_midpoint_drift':ratios.tolist(),
        'reference_values':reference.tolist(),'PDE_values':value.tolist(),
        'numerical_effect_resolved':bool(np.max(ratios)<.01)}


def record(a,meta,path,hash_value):
    checks=analyze(a,meta)
    passed=(checks['numerical_effect_resolved'] and checks['exact_harmonic_phase_residual']<3e-11
            and meta['flux_algebra']<3e-11)
    return {'schema':'NSC-EXACT-SOURCE-PHASE-v1','accountable_author':'Douglas Ek',
        'status':'PASS: exact-source-phase numerical effect resolved; physical accuracy OPEN' if passed else
                 'OPEN: exact-source-phase diagnostic did not close the targeted numerical effect',
        'control':{'family':'14_1','grid_points':801,'sample_times':65,'history_amplitude':0.,
            'field_runs':1,'CPU_budget_seconds':30.,'old_midpoint_forcing_drift':meta['old_midpoint_forcing_drift'],
            'minimum_improvement_factor':100.,'phase_and_flux_tolerance':3e-11,
            'method':'same SBP L0 and SAT amplitude, exact harmonic phase block; one constant augmented exponential',
            'axial_derivative':'actual semi-discrete PDE at fixed rho1, equivalent to partial_z on Sigma',
            'time_coordinates':'inherited numerical [0,.3], not a selected physical duration'},
        'measurements':checks,'runtime':{'CPU_seconds':meta['CPU_seconds'],'wall_seconds':meta['wall_seconds']},
        'scope':{'source_law_changed':False,'source_phase_fitted':False,'interior_counterforce':False,
            'stress_drift_subtracted':False,'incoming_C0_imposed':False,'metric_evolution':False,
            'physical_constraint_solution':False,'full_source_error_bound':None,
            'continuum_field_error_bound':None,'new_horizon_or_source_grid':False,
            'nonzero_history_accuracy_verified':False,'publication':False},
        'source_hashes':{p:meta['signature'][p] for p in SOURCES},
        'input_hashes':{p:meta['signature'][p] for p in INPUTS},
        'reused_artifacts':meta['reused_artifacts'],
        'payload':{'path':path,'sha256':hash_value},
        'reproducer':'python3 scripts/derive_nsc_exact_source_phase.py --check'}


def run():
    before=signatures();receipts=inputs()
    start=time.process_time();wall=time.monotonic()
    def stop(*_):raise TimeoutError('30 CPU-second exact-phase diagnostic budget exhausted')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,30.)
    try:
        x,phi,data,reused=STATE.RETARDED.inputs()
        owner=STATE.FourthOrderModePropagator(x,data['mass'].item(),data['angular'].item())
        _,_,incoming=STATE.incident(owner,phi,data,0)
        forcing,_=incoming(0.);L,algebra=owner.generator(owner.reference)
        energy=np.repeat(data['energies'],3);rows,ns=phi.shape
        generator=bmat([[L,csr_matrix(forcing)],[None,diags(-1j*energy)]],format='csr')
        initial=np.vstack((phi,np.eye(ns,dtype=complex)))
        times=np.linspace(0.,.3,65)
        solution=expm_multiply(generator,initial,start=0.,stop=.3,num=65,endpoint=True,
                              traceA=generator.diagonal().sum())
        ix=int(np.flatnonzero(x==1.)[0]);take=np.array([ix,len(x)+ix])
        N,beta,q,r=owner.reference[:,ix]
        R=rho1_restriction_map(N,q,beta,r)
        F=np.einsum('ab,tbs->tas',R,solution[:,take,:])
        Fz=np.array([R@(L[take]@v[:rows]+forcing[take]@v[rows:]) for v in solution])
        expected=np.eye(ns)[None]*np.exp(-1j*times[:,None,None]*energy[None,None,:])
        reference=(R@phi[take])[None]*np.exp(-1j*times[:,None,None]*energy)
        source=FixedSourcePreparation.from_signed_blocks(data['energies'],data['weights'],data['source'])
        channel=receipts[2]['actual_source_bindings']['channel']
        if channel['compact_mass']!=owner.mass or channel['angular_eigenvalue']!=owner.angular:
            raise ValueError('same inherited signed family required')
        meta={'signature':before,'CPU_seconds':time.process_time()-start,'wall_seconds':time.monotonic()-wall,
            'flux_algebra':algebra['norm_flux_algebra'],
            'parameters':{'mass':owner.mass,'angular':owner.angular,'axial_scale':np.sqrt(3*np.pi/2-4),
                'radius':r,'multiplicity':channel['copy_count']*channel['degeneracy']/2},
            'old_midpoint_forcing_drift':receipts[1]['reference_drift_separation']['reference_field_and_derivative_drift_indicator'],
            'reused_artifacts':[r['payload'] for r in reused]}
        a={'columns':F,'pde_axial_columns':Fz,'reference_columns':reference,
           'reference_axial_columns':-1j*energy*reference,'phase':solution[:,rows:],'expected_phase':expected,
           'source_covariance':source.covariance,'column_weights':source.column_weights,
           'z':times+receipts[0]['control']['support']['S1'],
           'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)}
        if signatures()!=before:raise ValueError('diagnostic owners changed during execution')
        raw=deterministic_npz_bytes(a);(ROOT/ARTIFACT).write_bytes(raw)
        result=record(a,meta,ARTIFACT,sha256(raw).hexdigest())
        (ROOT/OUTPUT).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
        return result
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)


def check():
    saved=json.loads((ROOT/OUTPUT).read_text());authenticate(saved)
    with np.load(ROOT/saved['payload']['path'],allow_pickle=False) as archive:
        a={k:archive[k] for k in archive.files}
    meta=json.loads(a['metadata_json'].tobytes())
    actual=record(a,meta,saved['payload']['path'],saved['payload']['sha256'])
    if saved!=actual:raise ValueError('exact-source-phase replay differs')
    return actual


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    result=run() if args.run else check()
    print(json.dumps({'status':result['status'],'runtime':result['runtime'],
        'measurements':{k:v for k,v in result['measurements'].items() if not k.endswith('values')}},indent=2))


if __name__=='__main__':main()
