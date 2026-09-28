#!/usr/bin/env python3
"""Reuse the saved 14_1 operator on its entire retained finite source inventory."""
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
import derive_nsc_ks_energy_propagator_control as O
from derive_nsc_evolved_incoming_constraints import coefficients_from_records
from recursive_horizons.nsc_ks_source_inventory import ReferenceSourcePanel
from recursive_horizons.nsc_ks_batched_constraints import (
    KSConstraintAccumulator, KSUpstreamBatch, bind_upstream_batch, signed_incoming_pair)
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT=ROOT/'results/development/nsc-ks-retained-response-control.json'
PAYLOAD=ROOT/'results/development/artifacts/nsc-ks-retained-response-control.npz'
INVENTORY='results/development/nsc-ks-source-inventory.json'
COEFFICIENT='results/development/nsc-incoming-surface-coefficients.json'
BRANCH='results/development/nsc-incoming-surface-regular-branch.json'
BASELINE='results/development/nsc-incoming-source-update-v5.json'
OWNED=('scripts/derive_nsc_ks_retained_response_control.py','docs/nsc-ks-retained-response-control.md',
       'src/recursive_horizons/nsc_ks_batched_constraints.py','src/recursive_horizons/nsc_ks_signed_state.py',
       'src/recursive_horizons/nsc_ks_local_constraints.py','src/recursive_horizons/nsc_ks_matter_difference.py')
CPU_CAP=120.
BATCH_ROWS=16


def digest(path):return O.digest(path)


def inputs():
    O.check();old,oldmeta,_=O.inputs()
    with np.load(O.PAYLOAD,allow_pickle=False) as f:op=O.restore({k:f[k] for k in f.files},oldmeta)
    record=json.loads((ROOT/INVENTORY).read_text())
    for path,expected in record['source_hashes'].items():
        if digest(path)!=expected:raise ValueError(f'source inventory owner changed: {path}')
    if digest(record['payload']['path'])!=record['payload']['sha256']:raise ValueError('source inventory changed')
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
    meta=json.loads(a['metadata_json'].tobytes())
    panels=[]
    for name,p in meta['panels'].items():
        if p['group']==14 and p['angular_sign']==1:
            panels.append(ReferenceSourcePanel(name,14,1,p['mass'],p['angular_magnitude'],
                *(a[name+'/'+key] for key in ('energies','weights','covariance','amplitudes_at_one')),p['provenance']))
    if not panels:raise ValueError('retained group14 positive-angular panels required')
    coeff=coefficients_from_records(json.loads((ROOT/COEFFICIENT).read_text()),json.loads((ROOT/BRANCH).read_text()))
    baseline=json.loads((ROOT/BASELINE).read_text())
    hashes={p:digest(p) for p in (INVENTORY,record['payload']['path'],COEFFICIENT,BRANCH,BASELINE,
        str(O.OUTPUT.relative_to(ROOT)),str(O.PAYLOAD.relative_to(ROOT)))}
    return op,meta,panels,coeff,baseline,hashes


def evaluate(op,meta,panels,coeff,baseline,saved=None):
    family=O.C.R.C.P.family(O.C.R.C.P.ALPHA)
    budget=local_error_budget(baseline_covered_regions=baseline['partial_error_budget']['covered_spectral_regions_action_error_upper'])
    acc=KSConstraintAccumulator(family,baseline['baseline']['action_gradient_approximant'],coeff,
        (float(op.z[0]),float(op.z[-1])),budget)
    for group in (0,13,23):acc.declare_radius_invariant_channel({**meta['channels'][group],'angular_sign':1})
    output={};receipts=[];rows=0
    for panel in panels:
        for part in panel.batches(BATCH_ROWS):
            prefix='upstream/'+part.name
            if saved is None:
                batch=bind_upstream_batch(part,rho_up=op.rho_up,**O.C.R.TIGHTER)
                receipt={k:getattr(batch,k) for k in ('panel_name','original_panel','rows','group',
                    'angular_sign','energy_sign','mass','angular','rho_up','preparation_digest')}
                receipt.update(continuation=dict(batch.continuation),provenance=dict(batch.provenance))
            else:
                receipt=saved[1][len(receipts)]
                src=FixedSourcePreparation.from_signed_blocks(part.energies,part.weights,part.covariance)
                batch=KSUpstreamBatch(**receipt,source=src,initial_columns=saved[0][prefix])
            receipts.append(receipt);output[prefix]=batch.initial_columns
            state=op.apply(batch.source,batch.initial_columns)
            negative=part.negative_partner(meta['config'])
            negsource=FixedSourcePreparation.from_signed_blocks(negative.energies,negative.weights,negative.covariance)
            pos,pbatch,neg,nbatch=signed_incoming_pair(state,batch,negsource)
            acc.add(pos,pbatch,{**meta['channels'][14],'angular_sign':1})
            acc.add(neg,nbatch,{**meta['channels'][14],'angular_sign':-1})
            rows+=len(part.energies)
    result=acc.finalize()
    for key in ('action_gradient','history_jacobian','local_reference_gradient_change','local_reference_gradient_tangent'):
        output[key]=result[key]
    for key,values in result['family_corrections'].items():output['correction/'+key]=values
    for key,values in result['family_matter_tangents'].items():output['tangent/'+key]=values
    output['z']=op.z
    measured={'positive_rows':rows,'signed_rows':2*rows,'positive_panels':len(panels),
        'batches':len(result['batch_records']),'sampled_signed_families':sorted(result['family_records']),
        'family_matter_change_maxima':{k:np.max(abs(v),axis=0).tolist() for k,v in result['family_corrections'].items()},
        'paired_matter_change_maxima':np.max(abs(sum(result['family_corrections'].values())),axis=0).tolist(),
        'partial_constraint_maxima':np.max(abs(result['action_gradient']),axis=0).tolist(),
        'retained_quadrature_measure':float(sum(p.weights.sum() for p in panels)),
        'full_source_error_bound':None,'certified_constraint_residual':None}
    return output,receipts,measured


def run():
    if OUTPUT.exists() or PAYLOAD.exists():raise FileExistsError('response control exists; use --check')
    op,meta,panels,coeff,baseline,hashes=inputs();signature={p:digest(p) for p in OWNED}
    def stop(*_):raise TimeoutError('retained response control exceeded120 CPU seconds')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CPU_CAP)
    cpu=time.process_time()
    try:arrays,receipts,measured=evaluate(op,meta,panels,coeff,baseline)
    finally:signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    cpu=time.process_time()-cpu
    if signature!={p:digest(p) for p in OWNED}:raise ValueError('owner changed during run')
    arrays['metadata_json']=np.frombuffer(json.dumps({'upstream_batches':receipts},sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays)
    r={'schema':'NSC-KS-RETAINED-RESPONSE-CONTROL-v1','accountable_author':'Douglas Ek',
       'status':'OPEN: two signed finite-source family responses; incomplete physical gate',
       'measurements':measured,'state_law':'C_Sigma[g]=U_g C_up U_gdagger; original fixed source',
       'control':{'history_amplitude':O.C.R.C.P.ALPHA,'local_interval':'S(1)+[.12,.18]',
           'operator_digest':op.digest,'upstream_fixed_before_history':True,
           'energy_folding_factor':1,'baseline_and_geometry_counted_once':True},
       'scope':{'physical_local_gate':'OPEN','constraint_root_claimed':False,
           'complete_changed_history_source_sum':False,'source_preparation_error_bound':None,
           'low_subgap_error_bound':None,'changed_history_UV_bound':None,
           'node_evolution_and_interpolation_error_bound':None,'between_node_error_bound':None,
           'zero_radius_response_groups':[0,13,23],'ell0_baseline_retained':True,
           'metric_timestep':False,'stress_drift_subtracted':False},
       'runtime':{'CPU_seconds':cpu,'CPU_cap':CPU_CAP,'new_field_operator_solves':0,
           'source_upstream_continuations':len(receipts)},
       'source_hashes':signature,'input_hashes':hashes,
       'payload':{'path':str(PAYLOAD.relative_to(ROOT)),'sha256':sha256(raw).hexdigest(),'bytes':len(raw)},
       'reproducer':'python3 scripts/derive_nsc_ks_retained_response_control.py --check'}
    PAYLOAD.write_bytes(raw);OUTPUT.write_text(json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return r


def check():
    r=json.loads(OUTPUT.read_text());op,meta,panels,coeff,baseline,hashes=inputs()
    if hashes!=r['input_hashes']:raise ValueError('response inputs changed')
    for p,h in r['source_hashes'].items():
        if digest(p)!=h:raise ValueError('response owner changed: '+p)
    if digest(PAYLOAD)!=r['payload']['sha256']:raise ValueError('response payload changed')
    with np.load(PAYLOAD,allow_pickle=False) as f:a={k:f[k] for k in f.files}
    saved=json.loads(a['metadata_json'].tobytes())
    current,_,measured=evaluate(op,meta,panels,coeff,baseline,(a,saved['upstream_batches']))
    for k,v in current.items():
        if not np.array_equal(v,a[k]):raise ValueError('response replay changed: '+k)
    if measured!=r['measurements']:raise ValueError('response measurements changed')
    return r


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();r=run() if args.run else check()
    print(json.dumps({k:r[k] for k in ('status','measurements','runtime')},indent=2))
