#!/usr/bin/env python3
"""One bounded operator solve, compared to saved actual-source column evolution."""
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
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_ks_difference_control as C
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_energy_propagator import evolve_energy_propagator, KSEnergyPropagator
from recursive_horizons.nsc_ks_source_envelope import (
    computational_z_grid, usual_axial_support, KSEnvelopeBinding, _sample_axial_profiles)
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT=ROOT/'results/development/nsc-ks-energy-propagator-control.json'
PAYLOAD=ROOT/'results/development/artifacts/nsc-ks-energy-propagator-control.npz'
INTERVAL=(-320.,320.)
DEGREE=48
CPU_CAP=60.
OWNED=('scripts/derive_nsc_ks_energy_propagator_control.py',
       'src/recursive_horizons/nsc_ks_energy_propagator.py','docs/nsc-ks-energy-propagator.md',
       'src/recursive_horizons/nsc_ks_difference_envelope.py',
       'src/recursive_horizons/nsc_ks_source_envelope.py',
       'src/recursive_horizons/nsc_evolved_incoming_constraints.py')


def digest(path):
    path=Path(path);return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def inputs():
    record=json.loads((ROOT/C.OUTPUT).read_text())
    # Check this input's complete recorded dependencies, without rerunning its solves.
    for mapping in ('source_hashes','input_hashes'):
        for path,expected in record[mapping].items():
            if digest(path)!=expected:raise ValueError(f'input owner changed: {path}')
    payload=record['payload']
    if digest(payload['path'])!=payload['sha256']:raise ValueError('input payload changed')
    with np.load(ROOT/payload['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
    meta=json.loads(a['metadata_json'].tobytes())
    return a,meta,{C.OUTPUT:digest(C.OUTPUT),payload['path']:payload['sha256']}


def source(a):
    return FixedSourcePreparation(a['source/covariance'],a['source/column_weights'],a['source/energies'])


def analyze(operator,a,meta):
    state=operator.apply(source(a),a['initial_canonical_columns'])
    family=C.R.C.P.family(C.R.C.P.ALPHA)
    state.require_history(family)
    if state.fixed_preparation_digest!=meta['cases']['tighter512']['preparation']:
        raise ValueError('original source preparation changed')
    keys={'F':'columns','Fz':'axial_columns','dF':'column_tangents','dFz':'axial_tangents'}
    residuals={key:float(np.max(abs(getattr(state,name)-a['tighter512/'+key]),initial=0)) for key,name in keys.items()}
    weights=a['source/column_weights'];cov=a['source/covariance']
    old=source_column_matter(*(a['tighter512/'+key]*weights for key in ('F','Fz')),cov,
                            *(a['tighter512/'+key]*weights for key in ('dF','dFz')),**meta['parameters'])
    new=source_column_matter(state.weighted_columns,state.weighted_axial_columns,cov,
                            state.weighted_column_tangents,state.weighted_axial_tangents,**meta['parameters'])
    matter=np.max(abs(old['action_gradient']-new['action_gradient']),axis=0)
    tangent=np.max(abs(old['action_gradient_tangent']-new['action_gradient_tangent']),axis=(0,1))
    passed=bool(np.max(matter)<=1e-11 and np.max(tangent)<=3e-8)
    return {'field_comparison':residuals,'matter_N_beta_comparison':matter.tolist(),
            'tangent_N_beta_comparison':tangent.tolist(),'matter_indicator_tolerance':1e-11,
            'tangent_indicator_tolerance':3e-8,'numerical_indicator_pass':passed}


def restore(a,meta):
    family=C.R.C.P.family(C.R.C.P.ALPHA);grid=computational_z_grid(512)
    w,U=_sample_axial_profiles(family.directions,grid)
    binding=KSEnvelopeBinding(grid,w,U,tuple(family.amplitudes),
        tuple(d.inner_radius for d in family.directions),tuple(d.outer_radius for d in family.directions),
        usual_axial_support(),meta['rho_up'],1.,**C.R.TIGHTER)
    return KSEnergyPropagator(INTERVAL,a['energies'],a['z'],*(a[key] for key in
        ('reference','difference','difference_z','tangent','tangent_z')),
        meta['parameters']['mass'],meta['parameters']['angular'],meta['rho_up'],binding,{})


def run():
    if OUTPUT.exists() or PAYLOAD.exists():raise FileExistsError('control exists; use --check')
    a,meta,input_hashes=inputs();sig={p:digest(p) for p in OWNED}
    def stop(*_):raise TimeoutError('one operator solve exceeded 60 CPU seconds')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CPU_CAP)
    cpu=time.process_time()
    try:
        op=evolve_energy_propagator(INTERVAL,DEGREE,C.R.C.P.family(C.R.C.P.ALPHA),
            computational_z_grid(512),a['target_z'],meta['parameters']['mass'],meta['parameters']['angular'],
            meta['rho_up'],axial_support=usual_axial_support(),**C.R.TIGHTER)
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    cpu=time.process_time()-cpu
    measured=analyze(op,a,meta)
    if sig!={p:digest(p) for p in OWNED}:raise ValueError('owner changed during run')
    raw=deterministic_npz_bytes({k:np.asarray(getattr(op,k)) for k in
        ('energies','z','reference','difference','difference_z','tangent','tangent_z')})
    result={'schema':'NSC-KS-ENERGY-PROPAGATOR-CONTROL-v1','accountable_author':'Douglas Ek',
        'status':('PASS' if measured['numerical_indicator_pass'] else 'OPEN')+': actual-source operator comparison; physical gate OPEN',
        'control':{'source_family':'14_1','operator_energy_interval':INTERVAL,'degree':DEGREE,
            'original_source_energy_labels':np.unique(a['source/energies']).tolist(),
            'local_interval':'S(1)+[.12,.18]','history_amplitude':C.R.C.P.ALPHA,
            'operator_digest':op.digest,'fixed_preparation_digest':meta['cases']['tighter512']['preparation']},
        'measurements':measured,'runtime':{'CPU_seconds':cpu,'CPU_cap':CPU_CAP,'new_operator_solves':1,
            'function_evaluations':op.diagnostics['function_evaluations']},
        'scope':{'physical_local_gate':'OPEN','indicators_are_error_bounds':False,
            'full_source_accuracy':None,'source_interpolated':False,'identity_is_operator_basis':True,
            'changed_history_UV_bound':None,'metric_timestep':False,'stress_drift_subtracted':False},
        'source_hashes':sig,'input_hashes':input_hashes,
        'payload':{'path':str(PAYLOAD.relative_to(ROOT)),'sha256':sha256(raw).hexdigest(),'bytes':len(raw)},
        'reproducer':'python3 scripts/derive_nsc_ks_energy_propagator_control.py --check'}
    PAYLOAD.write_bytes(raw);OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    return result


def check():
    r=json.loads(OUTPUT.read_text());a,meta,current=inputs()
    for path,expected in {**r['source_hashes'],**r['input_hashes']}.items():
        if digest(path)!=expected:raise ValueError(f'changed dependency: {path}')
    if current!=r['input_hashes'] or digest(PAYLOAD)!=r['payload']['sha256']:raise ValueError('changed input/payload')
    with np.load(PAYLOAD,allow_pickle=False) as f:op=restore({k:f[k] for k in f.files},meta)
    if op.digest!=r['control']['operator_digest']:raise ValueError('operator binding changed')
    if analyze(op,a,meta)!=r['measurements']:raise ValueError('comparison replay differs')
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    args=parser.parse_args();r=run() if args.run else check()
    print(json.dumps({k:r[k] for k in ('status','measurements','runtime')},indent=2))
